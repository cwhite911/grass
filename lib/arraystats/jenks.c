/* Jenks natural breaks classification (optimal univariate k-means)
 *
 * Computes the partition of sorted values into classes minimizing the sum
 * of within-class variances, using exact dynamic programming in
 * O(classes * n log n) time. The speedup over the classic O(classes * n^2)
 * Jenks-Caspall matrix algorithm comes from the concave Monge property of
 * the within-class sum-of-squares cost: the optimal split position is
 * nondecreasing along each DP layer, so each layer can be filled by
 * divide and conquer (Wang & Song 2011, Ckmeans.1d.dp). Independent
 * subranges of a layer are computed in parallel with OpenMP tasks.
 */

#include <float.h>

#include <grass/arraystats.h>
#include <grass/glocale.h>

/* Below this range size the divide and conquer recursion stays serial to
 * avoid OpenMP task overhead. */
#define JENKS_TASK_MIN 2048

struct jenks_ctx {
    int m;           /* number of distinct values */
    const double *W; /* prefix sums of weights, length m + 1 */
    const double *S; /* prefix sums of weighted shifted values */
    const double *Q; /* prefix sums of weighted squared shifted values */
    double *Dprev;   /* previous DP layer: optimal cost with c - 1 classes */
    double *Dcur;    /* current DP layer being filled */
    int *J;          /* argmin split index, J[(c - 2) * m + i] for layer c */
};

/* Weighted sum of squared deviations of distinct values j..i (inclusive)
 * from their weighted mean, in O(1) via prefix sums. Weights are positive,
 * so the denominator never vanishes. Clamped to zero because cancellation
 * can produce small negative results. */
static double ssd(const struct jenks_ctx *ctx, int j, int i)
{
    double ws = ctx->W[i + 1] - ctx->W[j];
    double s = ctx->S[i + 1] - ctx->S[j];
    double cost = (ctx->Q[i + 1] - ctx->Q[j]) - s * s / ws;

    return cost > 0.0 ? cost : 0.0;
}

/* Fill Dcur[lo..hi] of layer c, knowing the optimal split index of every
 * i in the range lies in [jlo, jhi]. Computes the middle entry by a linear
 * scan, then recurses on both halves with the search range narrowed by the
 * monotonicity of the argmin. The halves are independent, so they run as
 * OpenMP tasks when large enough. */
static void fill_row(struct jenks_ctx *ctx, int c, int lo, int hi, int jlo,
                     int jhi)
{
    int mid, jmax, bestj, j;
    double best;

    if (lo > hi)
        return;

    mid = lo + (hi - lo) / 2;
    /* Class c covering j..mid must not be empty. */
    jmax = (mid < jhi) ? mid : jhi;
    best = DBL_MAX;
    bestj = jlo;
    for (j = jlo; j <= jmax; j++) {
        /* Strict comparison keeps the smallest argmin: required for the
         * monotonicity bound and for results independent of thread count. */
        double cost = ctx->Dprev[j - 1] + ssd(ctx, j, mid);

        if (cost < best) {
            best = cost;
            bestj = j;
        }
    }
    ctx->Dcur[mid] = best;
    ctx->J[(size_t)(c - 2) * ctx->m + mid] = bestj;

    if (hi - lo > JENKS_TASK_MIN) {
#pragma omp task default(none) firstprivate(ctx, c, lo, mid, jlo, bestj)
        fill_row(ctx, c, lo, mid - 1, jlo, bestj);
#pragma omp task default(none) firstprivate(ctx, c, mid, hi, bestj, jhi)
        fill_row(ctx, c, mid + 1, hi, bestj, jhi);
#pragma omp taskwait
    }
    else {
        fill_row(ctx, c, lo, mid - 1, jlo, bestj);
        fill_row(ctx, c, mid + 1, hi, bestj, jhi);
    }
}

/* Optimal breaks for m >= 2 distinct values x[] with positive weights w[]
 * into k classes, 2 <= k <= m. Returns the goodness of variance fit. */
static double jenks_core(const double x[], const double w[], int m, int k,
                         double classbreaks[])
{
    struct jenks_ctx ctx;
    double *W, *S, *Q;
    double mu, wsum, tss, wcss, gvf;
    int i, c;

    /* Shift values by their weighted mean before accumulating the prefix
     * sums: computing Q - S^2 / W on raw large-magnitude values (e.g.
     * projected coordinates) loses all significant digits. */
    wsum = 0.0;
    mu = 0.0;
    for (i = 0; i < m; i++) {
        wsum += w[i];
        mu += w[i] * x[i];
    }
    mu /= wsum;

    W = G_malloc((m + 1) * sizeof(double));
    S = G_malloc((m + 1) * sizeof(double));
    Q = G_malloc((m + 1) * sizeof(double));
    W[0] = S[0] = Q[0] = 0.0;
    for (i = 0; i < m; i++) {
        double y = x[i] - mu;

        W[i + 1] = W[i] + w[i];
        S[i + 1] = S[i] + w[i] * y;
        Q[i + 1] = Q[i] + w[i] * y * y;
    }

    ctx.m = m;
    ctx.W = W;
    ctx.S = S;
    ctx.Q = Q;
    ctx.Dcur = G_malloc(m * sizeof(double));
    ctx.Dprev = G_malloc(m * sizeof(double));
    ctx.J = G_malloc((size_t)(k - 1) * m * sizeof(int));

    /* Layer 1: a single class covering values 0..i. */
    for (i = 0; i < m; i++)
        ctx.Dcur[i] = ssd(&ctx, 0, i);

    /* Layers are sequential (each reads the previous one); parallelism
     * lives in the divide and conquer inside fill_row. One thread runs the
     * layer loop, the others execute the spawned tasks. */
#pragma omp parallel default(none) shared(ctx, k, m)
    {
#pragma omp single
        {
            int layer;

            for (layer = 2; layer <= k; layer++) {
                double *tmp = ctx.Dprev;

                ctx.Dprev = ctx.Dcur;
                ctx.Dcur = tmp;
                /* Entries below layer - 1 stay undefined; later layers
                 * never read them. */
                fill_row(&ctx, layer, layer - 1, m - 1, layer - 1, m - 1);
            }
        }
    }

    tss = ssd(&ctx, 0, m - 1);
    wcss = ctx.Dcur[m - 1];

    /* Each break is the last value of its class: values are compared
     * against breaks with <= (see AS_class_frequencies). */
    i = m - 1;
    for (c = k; c >= 2; c--) {
        int j = ctx.J[(size_t)(c - 2) * m + i];

        classbreaks[c - 2] = x[j - 1];
        i = j - 1;
    }

    G_free(W);
    G_free(S);
    G_free(Q);
    G_free(ctx.Dcur);
    G_free(ctx.Dprev);
    G_free(ctx.J);

    gvf = tss > 0.0 ? 1.0 - wcss / tss : 1.0;
    /* The caller treats 0 as failure; a nonpositive fit can only arise
     * from rounding and is indistinguishable from "no fit". */
    if (gvf < DBL_EPSILON)
        gvf = DBL_EPSILON;

    return gvf;
}

/* Handles the degenerate class counts, then runs the DP core. */
static double jenks_breaks(const double x[], const double w[], int m,
                           int *nbreaks, double classbreaks[])
{
    int i;

    if (m < 1)
        return 0.0;

    if (m == 1) {
        G_warning(_("All values are identical, no class breaks possible"));
        *nbreaks = 0;
        return 1.0;
    }

    if (*nbreaks >= m) {
        G_warning(_("Number of classes reduced from %d to the number of "
                    "distinct values (%d)"),
                  *nbreaks + 1, m);
        *nbreaks = m - 1;
        for (i = 0; i < *nbreaks; i++)
            classbreaks[i] = x[i];
        return 1.0;
    }

    return jenks_core(x, w, m, *nbreaks + 1, classbreaks);
}

/*!
 * \brief Jenks natural breaks classification
 *
 * Computes the class breaks minimizing the sum of within-class variances
 * (Fisher-Jenks optimal partition). Breaks are values from the data array;
 * a value belongs to a class when it is less than or equal to the class
 * break. When there are fewer distinct values than requested classes, the
 * number of breaks is reduced accordingly.
 *
 * \param data array of values, sorted in ascending order
 * \param count number of values
 * \param[in,out] nbreaks number of requested breaks (classes - 1); reduced
 *                        when fewer distinct values exist
 * \param[out] classbreaks array of *nbreaks class break values
 *
 * \return goodness of variance fit in (0, 1], or 0.0 on failure
 */
double AS_class_jenks(const double data[], int count, int *nbreaks,
                      double classbreaks[])
{
    double *x, *w;
    double finfo;
    int i, m;

    if (count < 1)
        return 0.0;
    if (*nbreaks < 1)
        return 1.0;

    /* Collapse runs of equal values into (value, count) pairs; the DP cost
     * depends only on distinct values and their multiplicities. */
    x = G_malloc(count * sizeof(double));
    w = G_malloc(count * sizeof(double));
    m = 0;
    for (i = 0; i < count; i++) {
        if (m > 0 && data[i] == x[m - 1])
            w[m - 1] += 1.0;
        else {
            x[m] = data[i];
            w[m] = 1.0;
            m++;
        }
    }

    finfo = jenks_breaks(x, w, m, nbreaks, classbreaks);

    G_free(x);
    G_free(w);

    return finfo;
}

/*!
 * \brief Jenks natural breaks classification of weighted values
 *
 * Weighted variant of AS_class_jenks() for data given as distinct values
 * with multiplicities, e.g. a histogram of raster cell values. A value
 * with weight q contributes like q repetitions of that value.
 *
 * \param values array of values, sorted in ascending order
 * \param weights weight of each value; entries with a weight of zero or
 *                less are ignored
 * \param count number of values
 * \param[in,out] nbreaks number of requested breaks (classes - 1); reduced
 *                        when fewer distinct values exist
 * \param[out] classbreaks array of *nbreaks class break values
 *
 * \return goodness of variance fit in (0, 1], or 0.0 on failure
 */
double AS_class_jenks_weighted(const double values[], const double weights[],
                               int count, int *nbreaks, double classbreaks[])
{
    double *x, *w;
    double finfo;
    int i, m;

    if (count < 1)
        return 0.0;
    if (*nbreaks < 1)
        return 1.0;

    x = G_malloc(count * sizeof(double));
    w = G_malloc(count * sizeof(double));
    m = 0;
    for (i = 0; i < count; i++) {
        if (weights[i] <= 0.0)
            continue;
        if (m > 0 && values[i] == x[m - 1])
            w[m - 1] += weights[i];
        else {
            x[m] = values[i];
            w[m] = weights[i];
            m++;
        }
    }

    finfo = jenks_breaks(x, w, m, nbreaks, classbreaks);

    G_free(x);
    G_free(w);

    return finfo;
}
