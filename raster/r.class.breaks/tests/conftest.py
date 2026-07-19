import os

import pytest

import grass.script as gs
from grass.tools import Tools


@pytest.fixture(scope="module")
def session(tmp_path_factory):
    """Active session in an XY project with test rasters.

    The clusters raster has three well separated groups of values, so the
    natural breaks for three classes are known exactly: cells in columns
    1-39 hold values 10-14, columns 40-79 hold 100-104, and columns 80-100
    hold 1000-1004.
    """
    project = tmp_path_factory.mktemp("r_class_breaks") / "xy_test"
    gs.create_project(project)
    with (
        gs.setup.init(project, env=os.environ.copy()) as session,
        Tools(session=session) as tools,
    ):
        tools.g_region(s=0, n=100, w=0, e=100, res=1)
        tools.r_mapcalc(
            expression="clusters = if(col() < 40, 10, if(col() < 80, 100, "
            "1000)) + row() % 5"
        )
        tools.r_mapcalc(
            expression="clusters_float = float(if(col() < 40, 10, "
            "if(col() < 80, 100, 1000)) + row() % 5)"
        )
        tools.r_mapcalc(expression="constant = 42")
        yield session


@pytest.fixture
def tools(session):
    return Tools(session=session)
