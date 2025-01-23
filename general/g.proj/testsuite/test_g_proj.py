"""g.proj tests

(C) 2023 by the GRASS Development Team
This program is free software under the GNU General Public
License (>=v2). Read the file COPYING that comes with GRASS
for details.

:author: Anna Petrasova & Corey White
"""

import json
from grass.gunittest.case import TestCase
from grass.gunittest.main import test
from grass.gunittest.gmodules import SimpleModule


def validate_json_output(data, expected):
    """Validate the returned JSON output with the expected json output."""
    assert data
    assert (
        expected.keys() == data.keys()
    ), f"test failed: expected {expected} but got {data.keys()}"


class GProjWKTTestCase(TestCase):
    """Test g.proj with WKT output"""

    def test_wkt_output(self):
        """Test if g.proj returns WKT"""
        module = SimpleModule("g.proj", flags="w")
        self.assertModule(module)
        self.assertIn("PROJCRS", module.outputs.stdout)

    def test_verify_datum_output(self):
        """Test if g.proj returns verified datum and transformation parameters"""
        module = SimpleModule("g.proj", flags="d")
        self.assertModule(module)
        self.assertIn("GRASS datum code: nad83", module.outputs.stdout)

    def test_json_no_output_format(self):
        """Test if g.proj returns JSON"""
        module = SimpleModule("g.proj", format="json", flags="")
        self.assertModuleFail(module)

    def test_wtk_json_output(self):
        """Test if g.proj returns wtk JSON"""
        module = SimpleModule("g.proj", format="json", flags="w")
        self.runModule(module)
        results = json.loads(module.outputs.stdout)
        expected = {
            "wkt": 'PROJCRS["unknown",\n    BASEGEOGCRS["grs80",\n        DATUM["North American Datum 1983",\n            ELLIPSOID["Geodetic_Reference_System_1980",6378137,298.257222101,\n                LENGTHUNIT["metre",1]],\n            ID["EPSG",6269]],\n        PRIMEM["Greenwich",0,\n            ANGLEUNIT["degree",0.0174532925199433],\n            ID["EPSG",8901]]],\n    CONVERSION["unnamed",\n        METHOD["Lambert Conic Conformal (2SP)",\n            ID["EPSG",9802]],\n        PARAMETER["Latitude of false origin",33.75,\n            ANGLEUNIT["degree",0.0174532925199433],\n            ID["EPSG",8821]],\n        PARAMETER["Longitude of false origin",-79,\n            ANGLEUNIT["degree",0.0174532925199433],\n            ID["EPSG",8822]],\n        PARAMETER["Latitude of 1st standard parallel",36.1666666666667,\n            ANGLEUNIT["degree",0.0174532925199433],\n            ID["EPSG",8823]],\n        PARAMETER["Latitude of 2nd standard parallel",34.3333333333333,\n            ANGLEUNIT["degree",0.0174532925199433],\n            ID["EPSG",8824]],\n        PARAMETER["Easting at false origin",609601.22,\n            LENGTHUNIT["metre",1],\n            ID["EPSG",8826]],\n        PARAMETER["Northing at false origin",0,\n            LENGTHUNIT["metre",1],\n            ID["EPSG",8827]]],\n    CS[Cartesian,2],\n        AXIS["easting",east,\n            ORDER[1],\n            LENGTHUNIT["metre",1,\n                ID["EPSG",9001]]],\n        AXIS["northing",north,\n            ORDER[2],\n            LENGTHUNIT["metre",1,\n                ID["EPSG",9001]]]]'
        }
        validate_json_output(results, expected)

    def test_wtk_json_flat_output(self):
        """Test if g.proj returns wtk with flat output JSON"""
        module = SimpleModule("g.proj", format="json", flags="wf")
        self.runModule(module)
        results = json.loads(module.outputs.stdout)
        expected = {
            "wkt": 'PROJCRS["unknown",BASEGEOGCRS["grs80",DATUM["North American Datum 1983",ELLIPSOID["Geodetic_Reference_System_1980",6378137,298.257222101,LENGTHUNIT["metre",1]],ID["EPSG",6269]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]]],CONVERSION["unnamed",METHOD["Lambert Conic Conformal (2SP)",ID["EPSG",9802]],PARAMETER["Latitude of false origin",33.75,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8821]],PARAMETER["Longitude of false origin",-79,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8822]],PARAMETER["Latitude of 1st standard parallel",36.1666666666667,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8823]],PARAMETER["Latitude of 2nd standard parallel",34.3333333333333,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8824]],PARAMETER["Easting at false origin",609601.22,LENGTHUNIT["metre",1],ID["EPSG",8826]],PARAMETER["Northing at false origin",0,LENGTHUNIT["metre",1],ID["EPSG",8827]]],CS[Cartesian,2],AXIS["easting",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["northing",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]'
        }
        validate_json_output(results, expected)

    def test_wtk_esri_json_output(self):
        """Test if g.proj returns wtk with ESRI-style format JSON"""
        module = SimpleModule("g.proj", format="json", flags="we")
        self.runModule(module)
        results = json.loads(module.outputs.stdout)
        expected = {
            "wkt": 'PROJCS["unknown",\n    GEOGCS["GCS_grs80",\n        DATUM["D_North_American_1983",\n            SPHEROID["Geodetic_Reference_System_1980",6378137.0,298.257222101]],\n        PRIMEM["Greenwich",0.0],\n        UNIT["Degree",0.0174532925199433]],\n    PROJECTION["Lambert_Conformal_Conic"],\n    PARAMETER["False_Easting",609601.22],\n    PARAMETER["False_Northing",0.0],\n    PARAMETER["Central_Meridian",-79.0],\n    PARAMETER["Standard_Parallel_1",36.1666666666667],\n    PARAMETER["Standard_Parallel_2",34.3333333333333],\n    PARAMETER["Latitude_Of_Origin",33.75],\n    UNIT["Meter",1.0]]'
        }
        validate_json_output(results, expected)

    def test_grass_format_json_output(self):
        """Test if g.proj returns GRASS format JSON"""
        module = SimpleModule("g.proj", format="json", flags="p")
        self.runModule(module)
        results = json.loads(module.outputs.stdout)
        expected = {
            "name": "Lambert Conformal Conic",
            "proj": "lcc",
            "datum": "nad83",
            "a": "6378137.0",
            "es": "0.006694380022900787",
            "lat_1": "36.16666666666666",
            "lat_2": "34.33333333333334",
            "lat_0": "33.75",
            "lon_0": "-79",
            "x_0": "609601.22",
            "y_0": "0",
            "no_defs": "defined",
            "unit": "Meter",
            "units": "Meters",
            "meters": "1",
        }
        validate_json_output(results, expected)

    def test_proj4_format_json_output(self):
        """Test if g.proj returns PROJ.4 format JSON"""
        module = SimpleModule("g.proj", format="json", flags="j")
        self.runModule(module)
        results = json.loads(module.outputs.stdout)
        expected = {
            "proj4": "+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +no_defs +a=6378137 +rf=298.257222101 +towgs84=0.000,0.000,0.000 +type=crs  +to_meter=1"
        }
        validate_json_output(results, expected)

    def test_proj4_format_json_flat_output(self):
        """Test if g.proj returns PROJ.4 format with flat output JSON"""
        module = SimpleModule("g.proj", format="json", flags="jf")
        self.runModule(module)
        results = json.loads(module.outputs.stdout)
        expected = {
            "proj4": "+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +no_defs +a=6378137 +rf=298.257222101 +towgs84=0.000,0.000,0.000 +type=crs  +to_meter=1"
        }
        validate_json_output(results, expected)


if __name__ == "__main__":
    test()
