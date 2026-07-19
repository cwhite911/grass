import pytest

from grass.exceptions import CalledModuleError


def test_jenks_breaks_exact(tools):
    """Integer map: breaks fall exactly on the last value of each cluster."""
    result = tools.r_class_breaks(input="clusters", nclasses=3, format="json")
    assert result["algorithm"] == "jen"
    assert result["classes"] == 3
    assert result["breaks"] == [14, 104]
    assert [i["cell_count"] for i in result["intervals"]] == [3900, 4000, 2100]
    assert result["goodness_of_variance_fit"] == pytest.approx(1.0, abs=1e-3)


def test_jenks_breaks_float_map(tools):
    """Floating-point map: breaks approximate the cluster gaps to bin width."""
    result = tools.r_class_breaks(input="clusters_float", nclasses=3, format="json")
    binwidth = (1004 - 10) / 10000
    breaks = result["breaks"]
    assert len(breaks) == 2
    assert breaks[0] == pytest.approx(14, abs=binwidth)
    assert breaks[1] == pytest.approx(104, abs=binwidth)
    assert [i["cell_count"] for i in result["intervals"]] == [3900, 4000, 2100]


def test_reference_values(tools):
    """Match reference breaks computed with jenkspy 0.4.1 on the same values.

    The raster holds each value of vals once (10x10 cells); expected breaks
    come from jenkspy.jenks_breaks(vals, n_classes=4)[1:-1].
    """
    vals = [
        2,
        3,
        5,
        8,
        9,
        12,
        13,
        14,
        18,
        19,
        21,
        22,
        24,
        25,
        28,
        30,
        31,
        34,
        36,
        38,
        41,
        43,
        44,
        47,
        49,
        51,
        53,
        56,
        57,
        59,
        101,
        103,
        104,
        107,
        109,
        111,
        113,
        116,
        117,
        119,
        121,
        123,
        125,
        127,
        129,
        131,
        133,
        135,
        137,
        139,
        201,
        203,
        205,
        207,
        209,
        211,
        213,
        215,
        217,
        219,
        221,
        223,
        225,
        227,
        229,
        231,
        233,
        235,
        237,
        239,
        241,
        243,
        245,
        247,
        249,
        251,
        253,
        255,
        257,
        259,
        301,
        303,
        305,
        307,
        309,
        311,
        313,
        315,
        317,
        319,
        321,
        323,
        325,
        327,
        329,
        331,
        333,
        335,
        337,
        339,
    ]
    tools.g_region(s=0, n=10, w=0, e=10, res=1)
    rows = "\n".join(
        " ".join(str(v) for v in vals[r * 10 : (r + 1) * 10]) for r in range(10)
    )
    from io import StringIO

    tools.r_in_ascii(
        input=StringIO(
            "north: 10\nsouth: 0\neast: 10\nwest: 0\nrows: 10\ncols: 10\n" + rows
        ),
        output="refmap",
    )
    result = tools.r_class_breaks(input="refmap", nclasses=4, format="json")
    tools.g_region(s=0, n=100, w=0, e=100, res=1)
    assert result["breaks"] == [59, 139, 259]


def test_thread_count_independence(tools):
    """Results must be identical regardless of the number of threads."""
    runs = [
        tools.r_class_breaks(input="clusters", nclasses=5, format="json", nprocs=n)
        for n in (1, 4)
    ]
    assert runs[0].text == runs[1].text


def test_breaks_only_flag(tools):
    """The -b flag prints just the comma-separated breaks."""
    result = tools.r_class_breaks(
        input="clusters", nclasses=3, flags="b", format="plain"
    )
    breaks = [float(v) for v in result.text.split(",")]
    assert breaks == [14, 104]


def test_breaks_only_json(tools):
    result = tools.r_class_breaks(
        input="clusters", nclasses=3, flags="b", format="json"
    )
    assert list(result) == [14, 104]


def test_recode_rules_roundtrip(tools):
    """The -r rules feed r.recode and produce the expected classes."""
    rules = tools.r_class_breaks(input="clusters", nclasses=3, flags="r").text
    lines = rules.strip().splitlines()
    assert len(lines) == 3
    assert lines[0].endswith(":3")
    from io import StringIO

    tools.r_recode(input="clusters", output="classes", rules=StringIO(rules))
    stats = tools.r_stats(input="classes", flags="cn", format="json")
    counts = {int(r["categories"][0]["category"]): r["count"] for r in stats}
    assert counts == {1: 3900, 2: 4000, 3: 2100}


def test_interval_algorithm(tools):
    """Equal intervals split the range 10..1004 into equal parts."""
    result = tools.r_class_breaks(
        input="clusters", nclasses=2, algorithm="int", format="json"
    )
    assert result["breaks"] == [pytest.approx((10 + 1004) / 2)]


def test_quantile_algorithm(tools):
    """Quantile break is the value whose cumulative count reaches half.

    Cumulative counts around the median: 4700 cells hold values up to 100,
    5500 up to 101, so the median cell value is 101.
    """
    result = tools.r_class_breaks(
        input="clusters", nclasses=2, algorithm="qua", format="json"
    )
    assert result["breaks"] == [101]
    counts = [i["cell_count"] for i in result["intervals"]]
    assert counts == [5500, 4500]


def test_more_classes_than_values(tools):
    """Class count is reduced to the number of distinct values."""
    result = tools.r_class_breaks(input="constant", nclasses=5, format="json")
    assert result["classes"] == 1
    assert result["breaks"] == []
    assert result["intervals"][0]["cell_count"] == 10000


def test_all_null_map_fails(tools):
    tools.r_mapcalc(expression="nullmap = null()")
    with pytest.raises(CalledModuleError):
        tools.r_class_breaks(input="nullmap", nclasses=3, format="json")


def test_csv_output(tools):
    result = tools.r_class_breaks(input="clusters", nclasses=3, format="csv")
    lines = result.text.strip().splitlines()
    assert lines[0] == "from,to,cell_count"
    assert len(lines) == 4


def test_plain_output(tools):
    result = tools.r_class_breaks(input="clusters", nclasses=3, format="plain")
    assert "Goodness of variance fit" in result.text
    assert "Cell count" in result.text
