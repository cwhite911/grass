#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <grass/gis.h>
#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/bitmap.h>
#include <grass/linkm.h>
#include <grass/vector.h>
#include <grass/glocale.h>
#include <grass/parson.h>

double compute_continuity_flux_term(double **h, double **u, double **v,
                                    double dx, double dy, int i, int j)
{
    // Compute the continuity flux term using the FDM
    double dh_dx = (h[i + 1][j] - h[i - 1][j]) / (2.0 * dx);
    double dh_dy = (h[i][j + 1] - h[i][j - 1]) / (2.0 * dy);

    double dh_dt = -((u[i + 1][j] - u[i - 1][j]) / (2.0 * dx) +
                     (v[i][j + 1] - v[i][j - 1]) / (2.0 * dy)) -
                   h[i][j] * (dh_dx + dh_dy);

    return dh_dt;
}

// Finite Difference Method (FDM) using x-momentum equation
// ∂(hu)/∂t + ∂(huv)/∂x + ∂(huv)/∂y + gh(∂h/∂x) = -ghf ∂z/∂x + τx
double compute_x_momentum_flux_term(double **h, double **u, double **v,
                                    double dx, double dy, double g, int i,
                                    int j)
{
    // Compute the x-momentum flux term using the FDM
    double du_dx = (u[i + 1][j] - u[i - 1][j]) / (2.0 * dx);
    double du_dy = (u[i][j + 1] - u[i][j - 1]) / (2.0 * dy);
    double dh_dx = (h[i + 1][j] - h[i - 1][j]) / (2.0 * dx);
    double dh_dy = (h[i][j + 1] - h[i][j - 1]) / (2.0 * dy);
    double dz_dx = (z[i + 1][j] - z[i - 1][j]) / (2.0 * dx);

    // Compute the bed shear stress in the x-direction using the FDM
    double tau_x =
        -g * h[i][j] * (du_dx + dh_dy) +
        2.0 * g * h[i][j] * f[i][j] * (u[i][j] * dz_dx + v[i][j] * dz_dy) +
        g * h[i][j] * h[i][j] *
            (du_dx * du_dx + 2.0 * du_dy * dh_dx + dh_dy * dh_dy);

    double du_dt =
        -((h[i][j] * u[i + 1][j] - h[i][j] * u[i - 1][j]) / (2.0 * dx) +
          (h[i][j] * u[i][j + 1] - h[i][j] * u[i][j - 1]) / (2.0 * dy)) +
        g * h[i][j] * dh_dx - g * h[i][j] * f[i][j] * dz_dx +
        tau_x[i][j]; // tau_x is the bed shear stress in the x-direction

    return du_dt;
}

// Finite Difference Method (FDM) using y-momentum equation
// ∂(hv)/∂t + ∂(huv)/∂x + ∂(hvv)/∂y + gh(∂h/∂y) = -ghf ∂z/∂y + τy
double compute_y_momentum_flux_term(double **h, double **u, double **v,
                                    double dx, double dy, double g, int i,
                                    int j)
{
    // Compute the y-momentum flux term using the FDM
    double dv_dx = (v[i + 1][j] - v[i - 1][j]) / (2.0 * dx);
    double dv_dy = (v[i][j + 1] - v[i][j - 1]) / (2.0 * dy);
    double dh_dy = (h[i][j + 1] - h[i][j - 1]) / (2.0 * dy);

    double dv_dt =
        -((h[i][j] * v[i + 1][j] - h[i][j] * v[i - 1][j]) / (2.0 * dx) +
          (h[i][j] * v[i][j + 1] - h[i][j] * v[i][j - 1]) / (2.0 * dy)) +
        g * h[i][j] * dh_dy - g * h[i][j] * f[i][j] * dz_dy +
        tau_y[i][j]; // tau_y is the bed shear stress in the y-direction

    return dv_dt;
}

double **compute_green_function(int N, int M)
{
    double **green_function = G_allocate_2d_d_raster(N, M);
    double decay_rate =
        0.2; /* specify the decay rate for the Green's function */

    // Compute the Green's function contributions for each grid cell
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            // Compute the distance from the source cell (0, 0) to the current
            // cell (i, j)
            double distance = sqrt(i * i + j * j);
            // Compute the Green's function value based on the exponential decay
            green_function[i][j] = exp(-decay_rate * distance);
        }
    }

    return green_function;
}

void apply_green_function_contributions(double **h, double **u, double **v,
                                        double **green_function, int N, int M,
                                        double **h_temp, double **u_temp,
                                        double **v_temp)
{
    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; j < M - 1; j++) {
            // Apply Green's function contribution to update water depth and
            // velocity
            h[i][j] += green_function[i][j] * h_temp[i][j];
            u[i][j] += green_function[i][j] * u_temp[i][j];
            v[i][j] += green_function[i][j] * v_temp[i][j];
        }
    }
}

// Function to apply rainfall source term for a given time step
void apply_rainfall_source(double **h, double rainfall_rate, double dt, int i,
                           int j)
{
    // Compute the increase in water depth due to rainfall during the time step
    double dh_rainfall = rainfall_rate * dt;
    // Apply the rainfall source to the water depth at grid cell (i, j)
    h[i][j] += dh_rainfall;
}

// Function to apply inflow source term for a given time step
void apply_inflow_source(double **h, double inflow_rate, double dt, int i,
                         int j)
{
    // Compute the increase in water depth due to inflow during the time step
    double dh_inflow = inflow_rate * dt;
    // Apply the inflow source to the water depth at grid cell (i, j)
    h[i][j] += dh_inflow;
}

// Function to validate model results with observed data and established models
void validate_model_results(double **model_output, double **observed_data,
                            double **established_model_output, int N, int M)
{
    // Calculate validation metrics
    double sum_squared_diff_model_observed = 0.0;
    double sum_squared_diff_model_established = 0.0;
    double sum_squared_diff_observed_established = 0.0;
    double sum_squared_diff_observed_mean = 0.0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            double diff_model_observed =
                model_output[i][j] - observed_data[i][j];
            double diff_model_established =
                model_output[i][j] - established_model_output[i][j];
            double diff_observed_established =
                observed_data[i][j] - established_model_output[i][j];
            double diff_observed_mean =
                observed_data[i][j] - mean_observed_data;

            sum_squared_diff_model_observed +=
                diff_model_observed * diff_model_observed;
            sum_squared_diff_model_established +=
                diff_model_established * diff_model_established;
            sum_squared_diff_observed_established +=
                diff_observed_established * diff_observed_established;
            sum_squared_diff_observed_mean +=
                diff_observed_mean * diff_observed_mean;
        }
    }

    // Calculate validation metrics such as Root Mean Squared Error (RMSE),
    // Nash-Sutcliffe Efficiency (NSE), etc.
    double rmse_model_observed =
        sqrt(sum_squared_diff_model_observed / (N * M));
    double rmse_model_established =
        sqrt(sum_squared_diff_model_established / (N * M));
    double rmse_observed_established =
        sqrt(sum_squared_diff_observed_established / (N * M));
    double nse = 1.0 - (sum_squared_diff_model_observed /
                        sum_squared_diff_observed_mean);

    // Print or store the validation metrics
    printf("Validation Metrics:\n");
    printf("RMSE (Model vs. Observed): %.6f\n", rmse_model_observed);
    printf("RMSE (Model vs. Established Model): %.6f\n",
           rmse_model_established);
    printf("RMSE (Observed vs. Established Model): %.6f\n",
           rmse_observed_established);
    printf("Nash-Sutcliffe Efficiency (NSE): %.6f\n", nse);
}

// Function to output the final water depth and velocity profiles as JSON
void output_simulation_results_as_rasters(double **h, int N, int M,
                                          int num_time_steps)
{
    // Open a new GRASS GIS raster map layer for each time step
    for (int t = 0; t < num_time_steps; t++) {
        char raster_name[100];
        sprintf(raster_name, "water_depth_timestep_%03d", t);

        Rast_put_d_cell_title(raster_name, "Water Depth at Time Step", 1);

        // Write the water depth values to the raster map
        RASTER_MAP_TYPE data_type = CELL_TYPE;
        RASTER_MAP_TYPE data_null_type = CELL_TYPE;
        struct Cell_head cellhd;
        DCELL *buffer;

        Rast_init_colors();
        Rast_open_new(&cellhd, raster_name, data_type, data_null_type);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                buffer[j] = h[t][i * M + j];
            }
            Rast_put_row(cellhd, buffer, data_type);
        }

        Rast_close(cellhd);
    }
}

// Function to output the final water depth and velocity profiles
void output_simulation_results(double **h, double **u, double **v, int N, int M)
{
    // Output water depth profile to a file
    FILE *h_file = fopen("water_depth_profile.txt", "w");
    if (h_file != NULL) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                fprintf(h_file, "%.6f ", h[i][j]);
            }
            fprintf(h_file, "\n");
        }
        fclose(h_file);
        printf("Water depth profile has been saved to "
               "'water_depth_profile.txt'.\n");
    }
    else {
        printf("Error: Unable to save water depth profile.\n");
    }

    // Output x-velocity profile to a file
    FILE *u_file = fopen("x_velocity_profile.txt", "w");
    if (u_file != NULL) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                fprintf(u_file, "%.6f ", u[i][j]);
            }
            fprintf(u_file, "\n");
        }
        fclose(u_file);
        printf(
            "X-velocity profile has been saved to 'x_velocity_profile.txt'.\n");
    }
    else {
        printf("Error: Unable to save x-velocity profile.\n");
    }

    // Output y-velocity profile to a file
    FILE *v_file = fopen("y_velocity_profile.txt", "w");
    if (v_file != NULL) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                fprintf(v_file, "%.6f ", v[i][j]);
            }
            fprintf(v_file, "\n");
        }
        fclose(v_file);
        printf(
            "Y-velocity profile has been saved to 'y_velocity_profile.txt'.\n");
    }
    else {
        printf("Error: Unable to save y-velocity profile.\n");
    }
}

// Function to create a JSON array from a 2D double array
JSON_Value *create_json_array(double **data, int N, int M)
{
    JSON_Value *array = json_value_init_array();
    JSON_Array *json_array = json_value_get_array(array);

    for (int i = 0; i < N; i++) {
        JSON_Value *row = json_value_init_array();
        JSON_Array *json_row = json_value_get_array(row);

        for (int j = 0; j < M; j++) {
            json_array_append_number(json_row, data[i][j]);
        }
        json_array_append_value(json_array, row);
    }

    return array;
}

// Function to output the final water depth and velocity profiles as JSON
void output_simulation_results_as_json(double **h, double **u, double **v,
                                       int N, int M)
{
    JSON_Value *root_value = json_value_init_object();
    JSON_Object *root_object = json_value_get_object(root_value);

    // Create JSON arrays for water depth, x-velocity, and y-velocity profiles
    JSON_Value *h_array = create_json_array(h, N, M);
    JSON_Value *u_array = create_json_array(u, N, M);
    JSON_Value *v_array = create_json_array(v, N, M);

    // Add JSON arrays to the root object
    json_object_set_value(root_object, "water_depth_profile", h_array);
    json_object_set_value(root_object, "x_velocity_profile", u_array);
    json_object_set_value(root_object, "y_velocity_profile", v_array);

    // Serialize the JSON to a string
    char *json_string = json_serialize_to_string_pretty(root_value);

    // Output JSON string to a file
    FILE *json_file = fopen("simulation_results.json", "w");
    if (json_file != NULL) {
        fprintf(json_file, "%s", json_string);
        fclose(json_file);
        printf("Simulation results have been saved to "
               "'simulation_results.json'.\n");
    }
    else {
        printf("Error: Unable to save simulation results.\n");
    }

    // Free the JSON objects
    json_value_free(root_value);
    json_free_serialized_string(json_string);
}

int main(int argc, char *argv[])
{
    // Initialize GRASS GIS environment
    G_gisinit(argv[0]);

    // Set up parameters and grid
    int N = 1000; /* Number of rows */
    int M = 1000; /* Number of columns */

    double dx = 1.0; /* Spatial resolution in x direction */
    double dy = 1.0; /* Spatial resolution in y direction */
    double dt = 2.0; /* Time step size */
    double T = 60.0; /* Total simulation time */
    double g = 9.81; /* Gravitational acceleration */
    int num_time_steps = T / dt;

    // Initialize water depth, velocity components, and Green's function
    double *h = G_allocate_2d_d_raster(N, M);
    double *u = G_allocate_2d_d_raster(N, M);
    double *v = G_allocate_2d_d_raster(N, M);
    double *green_function = G_allocate_2d_d_raster(N, M);

    // Define initial and boundary conditions
    /* Implement functions to set initial and boundary conditions */
    double rainfall_rate = 50.0; /* Rainfall rate */
    double inflow_rate = 0.0;    /* Inflow rate */

    // Time stepping loop
    for (int t = 0; t < num_time_steps; t++) {
        // Initialize temporary arrays for updating h, u, v at the current time
        // step
        double **h_temp = G_allocate_2d_d_raster(N, M);
        double **u_temp = G_allocate_2d_d_raster(N, M);
        double **v_temp = G_allocate_2d_d_raster(N, M);

        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < M - 1; j++) {
                // Computevalidate_model_results
                // Compute flux terms using FDM for x-momentum and y-momentum
                // equations
                double du_dt =
                    compute_x_momentum_flux_term(h, u, v, dx, dy, g, i, j);
                double dv_dt =
                    compute_y_momentum_flux_term(h, u, v, dx, dy, g, i, j);
                double dh_dt =
                    compute_continuity_flux_term(h, u, v, dx, dy, i, j);

                // Update water depth and velocity using FDM
                h_temp[i][j] = h[i][j] + dt * dh_dt;
                u_temp[i][j] = u[i][j] + dt * du_dt;
                v_temp[i][j] = v[i][j] + dt * dv_dt;
            }
        }

        // Apply Green's function contributions
        /* Implement functions to compute and apply Green's function
         * contributions */
        apply_green_function_contributions(h, u, v, green_function, N, M,
                                           h_temp, u_temp, v_temp);

        // Free temporary arrays
        G_free_2d_d_raster(h_temp);
        G_free_2d_d_raster(u_temp);
        G_free_2d_d_raster(v_temp);

        // Optional: Incorporate source terms (e.g., rainfall, inflow) for this
        // time step
        /* Implement functions to apply source terms for this time step */
        // Apply rainfall source for this time step
        // apply_rainfall_source(h, rainfall_rate, dt, i_rainfall, j_rainfall);

        // Apply inflow source for this time step
        // apply_inflow_source(h, inflow_rate, dt, i_inflow, j_inflow);
    }

    // Model Validation: Compare simulation results with observed data or other
    // established models
    /* Implement validation functions */
    // validate_model_results(h, observed_data, established_model_output, N, M);

    // Output the final water depth and velocity profiles as the simulation
    // results
    /* Implement functions to output the simulation results */

    output_simulation_results_as_rasters(h, N, M, num_time_steps);
    output_simulation_results(h, u, v, N, M);
    output_simulation_results_as_json(h, u, v, N, M);

    // Free memory and clean up
    G_free_2d_d_raster(h);
    G_free_2d_d_raster(u);
    G_free_2d_d_raster(v);
    G_free_2d_d_raster(green_function);

    // Close GRASS GIS environment
    G_gisinit(argv[0]);

    return 0;
}
