"""
Name:       r.what test
Purpose:    Tests r.what module and its options.

Author:     Shubham Sharma, Google Code-in 2018
Copyright:  (C) 2018 by Shubham Sharma and the GRASS Development Team
Licence:    This program is free software under the GNU General Public
            License (>=v2). Read the file COPYING that comes with GRASS
            for details.
"""

from grass.gunittest.case import TestCase
from grass.gunittest.gmodules import SimpleModule
import os


class TestRasterWhat(TestCase):
    map1 = 'boundary_county_500m'
    map2 = 'landuse96_28m,aspect'
    coordinates = (633614.08, 224125.12, 632972.36, 225382.87)
    points = 'comm_colleges'
    refrence_points="""145096.8591495|154534.264883875||39
616341.4371495|146049.750883875||51
410595.7191495|174301.828883875||71
734153.6871495|169168.437883875||107
706338.2501495|54889.417883875||129
758009.7501495|112019.898883875||133
754002.7501495|200902.234883875||147
704771.7501495|183364.484883875||191
399187.0631495|220018.859883875||35
685098.9371495|33282.089883875||19
577750.8131495|257153.109883875||1
794095.5621495|199742.671883875||13
634688.2501495|100629.616883875||17
287638.7811495|207582.624883875||21
366218.5321495|222940.625883875||23
385212.4371495|236593.109883875||27
628137.4371495|63995.550883875||47
782600.5631495|152698.890883875||49
502813.9381495|235232.577883875||57
705922.6251495|136589.359883875||61
620397.8131495|246847.640883875||63
738465.3751495|237233.983883875||65
708944.7501495|247632.296883875||127
526666.6871495|249780.312883875||81
733439.3741495|298005.311883875||83
253886.0321495|204702.171883875||87
298337.5621495|178131.233883875||89
788989.3121495|284544.562883875||91
438378.3431495|227157.890883875||97
227406.0621495|190462.640883875||99
668782.8121495|195518.718883875||101
596325.2501495|190688.374883875||105
777034.8131495|229280.030883875||117
337714.9681495|219991.202883875||111
442913.7181495|164364.608883875||119
528779.3741495|179793.890883875||123
571926.0621495|162885.202883875||125
801852.0621495|154744.077883875||137
611371.9381495|295204.467883875||145
534602.9371495|213605.812883875||151
608960.3121495|102612.546883875||155
541667.1251495|293002.124883875||157
475819.8131495|210685.359883875||159
343685.3441495|182222.358883875||161
670160.3121495|140070.249883875||163
497214.4691495|177053.639883875||167
455564.7811495|293568.937883875||171
636209.2501495|210858.014883875||183
413565.1871495|267534.562883875||193
709485.0001495|220209.202883875||195
815182.1881495|110335.429883875||31
495151.3131495|258002.374883875||67
859129.4371495|288849.968883875||139
545200.0621495|127391.304883875||153
655627.6871495|285597.342883875||181
382174.7811495|171687.280883875||45
511235.1561495|135372.358883875||7
332533.5941495|242831.139883875||121
"""
    refrence_coordinates="""633614.08|224125.12||2|209.5939
632972.36|225382.87||15|140.7571
"""
    refrence_cats="""1|145096.8591495|154534.264883875||39
2|616341.4371495|146049.750883875||51
3|410595.7191495|174301.828883875||71
4|734153.6871495|169168.437883875||107
5|706338.2501495|54889.417883875||129
6|758009.7501495|112019.898883875||133
7|754002.7501495|200902.234883875||147
8|704771.7501495|183364.484883875||191
9|399187.0631495|220018.859883875||35
10|685098.9371495|33282.089883875||19
11|577750.8131495|257153.109883875||1
12|794095.5621495|199742.671883875||13
13|634688.2501495|100629.616883875||17
14|287638.7811495|207582.624883875||21
15|366218.5321495|222940.625883875||23
16|385212.4371495|236593.109883875||27
17|628137.4371495|63995.550883875||47
18|782600.5631495|152698.890883875||49
19|502813.9381495|235232.577883875||57
20|705922.6251495|136589.359883875||61
21|620397.8131495|246847.640883875||63
22|738465.3751495|237233.983883875||65
23|708944.7501495|247632.296883875||127
24|526666.6871495|249780.312883875||81
25|733439.3741495|298005.311883875||83
26|253886.0321495|204702.171883875||87
27|298337.5621495|178131.233883875||89
28|788989.3121495|284544.562883875||91
29|438378.3431495|227157.890883875||97
30|227406.0621495|190462.640883875||99
31|668782.8121495|195518.718883875||101
32|596325.2501495|190688.374883875||105
33|777034.8131495|229280.030883875||117
34|337714.9681495|219991.202883875||111
35|442913.7181495|164364.608883875||119
36|528779.3741495|179793.890883875||123
37|571926.0621495|162885.202883875||125
38|801852.0621495|154744.077883875||137
39|611371.9381495|295204.467883875||145
40|534602.9371495|213605.812883875||151
41|608960.3121495|102612.546883875||155
42|541667.1251495|293002.124883875||157
43|475819.8131495|210685.359883875||159
44|343685.3441495|182222.358883875||161
45|670160.3121495|140070.249883875||163
46|497214.4691495|177053.639883875||167
47|455564.7811495|293568.937883875||171
48|636209.2501495|210858.014883875||183
49|413565.1871495|267534.562883875||193
50|709485.0001495|220209.202883875||195
51|815182.1881495|110335.429883875||31
52|495151.3131495|258002.374883875||67
53|859129.4371495|288849.968883875||139
54|545200.0621495|127391.304883875||153
55|655627.6871495|285597.342883875||181
56|382174.7811495|171687.280883875||45
57|511235.1561495|135372.358883875||7
58|332533.5941495|242831.139883875||121
"""
    refrence_csv="""easting,northing,site_name,boundary_county_500m
145096.8591495,154534.264883875,,39
616341.4371495,146049.750883875,,51
410595.7191495,174301.828883875,,71
734153.6871495,169168.437883875,,107
706338.2501495,54889.417883875,,129
758009.7501495,112019.898883875,,133
754002.7501495,200902.234883875,,147
704771.7501495,183364.484883875,,191
399187.0631495,220018.859883875,,35
685098.9371495,33282.089883875,,19
577750.8131495,257153.109883875,,1
794095.5621495,199742.671883875,,13
634688.2501495,100629.616883875,,17
287638.7811495,207582.624883875,,21
366218.5321495,222940.625883875,,23
385212.4371495,236593.109883875,,27
628137.4371495,63995.550883875,,47
782600.5631495,152698.890883875,,49
502813.9381495,235232.577883875,,57
705922.6251495,136589.359883875,,61
620397.8131495,246847.640883875,,63
738465.3751495,237233.983883875,,65
708944.7501495,247632.296883875,,127
526666.6871495,249780.312883875,,81
733439.3741495,298005.311883875,,83
253886.0321495,204702.171883875,,87
298337.5621495,178131.233883875,,89
788989.3121495,284544.562883875,,91
438378.3431495,227157.890883875,,97
227406.0621495,190462.640883875,,99
668782.8121495,195518.718883875,,101
596325.2501495,190688.374883875,,105
777034.8131495,229280.030883875,,117
337714.9681495,219991.202883875,,111
442913.7181495,164364.608883875,,119
528779.3741495,179793.890883875,,123
571926.0621495,162885.202883875,,125
801852.0621495,154744.077883875,,137
611371.9381495,295204.467883875,,145
534602.9371495,213605.812883875,,151
608960.3121495,102612.546883875,,155
541667.1251495,293002.124883875,,157
475819.8131495,210685.359883875,,159
343685.3441495,182222.358883875,,161
670160.3121495,140070.249883875,,163
497214.4691495,177053.639883875,,167
455564.7811495,293568.937883875,,171
636209.2501495,210858.014883875,,183
413565.1871495,267534.562883875,,193
709485.0001495,220209.202883875,,195
815182.1881495,110335.429883875,,31
495151.3131495,258002.374883875,,67
859129.4371495,288849.968883875,,139
545200.0621495,127391.304883875,,153
655627.6871495,285597.342883875,,181
382174.7811495,171687.280883875,,45
511235.1561495,135372.358883875,,7
332533.5941495,242831.139883875,,121"""

    refrence_flag_i="""145096.8591495|154534.264883875||39
616341.4371495|146049.750883875||51
410595.7191495|174301.828883875||71
734153.6871495|169168.437883875||107
706338.2501495|54889.417883875||129
758009.7501495|112019.898883875||133
754002.7501495|200902.234883875||147
704771.7501495|183364.484883875||191
399187.0631495|220018.859883875||35
685098.9371495|33282.089883875||19
577750.8131495|257153.109883875||1
794095.5621495|199742.671883875||13
634688.2501495|100629.616883875||17
287638.7811495|207582.624883875||21
366218.5321495|222940.625883875||23
385212.4371495|236593.109883875||27
628137.4371495|63995.550883875||47
782600.5631495|152698.890883875||49
502813.9381495|235232.577883875||57
705922.6251495|136589.359883875||61
620397.8131495|246847.640883875||63
738465.3751495|237233.983883875||65
708944.7501495|247632.296883875||127
526666.6871495|249780.312883875||81
733439.3741495|298005.311883875||83
253886.0321495|204702.171883875||87
298337.5621495|178131.233883875||89
788989.3121495|284544.562883875||91
438378.3431495|227157.890883875||97
227406.0621495|190462.640883875||99
668782.8121495|195518.718883875||101
596325.2501495|190688.374883875||105
777034.8131495|229280.030883875||117
337714.9681495|219991.202883875||111
442913.7181495|164364.608883875||119
528779.3741495|179793.890883875||123
571926.0621495|162885.202883875||125
801852.0621495|154744.077883875||137
611371.9381495|295204.467883875||145
534602.9371495|213605.812883875||151
608960.3121495|102612.546883875||155
541667.1251495|293002.124883875||157
475819.8131495|210685.359883875||159
343685.3441495|182222.358883875||161
670160.3121495|140070.249883875||163
497214.4691495|177053.639883875||167
455564.7811495|293568.937883875||171
636209.2501495|210858.014883875||183
413565.1871495|267534.562883875||193
709485.0001495|220209.202883875||195
815182.1881495|110335.429883875||31
495151.3131495|258002.374883875||67
859129.4371495|288849.968883875||139
545200.0621495|127391.304883875||153
655627.6871495|285597.342883875||181
382174.7811495|171687.280883875||45
511235.1561495|135372.358883875||7
332533.5941495|242831.139883875||121
"""
    refrence_flag_f="""145096.8591495|154534.264883875||39|
616341.4371495|146049.750883875||51|
410595.7191495|174301.828883875||71|
734153.6871495|169168.437883875||107|
706338.2501495|54889.417883875||129|
758009.7501495|112019.898883875||133|
754002.7501495|200902.234883875||147|
704771.7501495|183364.484883875||191|
399187.0631495|220018.859883875||35|
685098.9371495|33282.089883875||19|
577750.8131495|257153.109883875||1|
794095.5621495|199742.671883875||13|
634688.2501495|100629.616883875||17|
287638.7811495|207582.624883875||21|
366218.5321495|222940.625883875||23|
385212.4371495|236593.109883875||27|
628137.4371495|63995.550883875||47|
782600.5631495|152698.890883875||49|
502813.9381495|235232.577883875||57|
705922.6251495|136589.359883875||61|
620397.8131495|246847.640883875||63|
738465.3751495|237233.983883875||65|
708944.7501495|247632.296883875||127|
526666.6871495|249780.312883875||81|
733439.3741495|298005.311883875||83|
253886.0321495|204702.171883875||87|
298337.5621495|178131.233883875||89|
788989.3121495|284544.562883875||91|
438378.3431495|227157.890883875||97|
227406.0621495|190462.640883875||99|
668782.8121495|195518.718883875||101|
596325.2501495|190688.374883875||105|
777034.8131495|229280.030883875||117|
337714.9681495|219991.202883875||111|
442913.7181495|164364.608883875||119|
528779.3741495|179793.890883875||123|
571926.0621495|162885.202883875||125|
801852.0621495|154744.077883875||137|
611371.9381495|295204.467883875||145|
534602.9371495|213605.812883875||151|
608960.3121495|102612.546883875||155|
541667.1251495|293002.124883875||157|
475819.8131495|210685.359883875||159|
343685.3441495|182222.358883875||161|
670160.3121495|140070.249883875||163|
497214.4691495|177053.639883875||167|
455564.7811495|293568.937883875||171|
636209.2501495|210858.014883875||183|
413565.1871495|267534.562883875||193|
709485.0001495|220209.202883875||195|
815182.1881495|110335.429883875||31|
495151.3131495|258002.374883875||67|
859129.4371495|288849.968883875||139|
545200.0621495|127391.304883875||153|
655627.6871495|285597.342883875||181|
382174.7811495|171687.280883875||45|
511235.1561495|135372.358883875||7|
332533.5941495|242831.139883875||121|
"""
    refrence_flag_r="""145096.8591495|154534.264883875||39|006:255:000
616341.4371495|146049.750883875||51|000:255:071
410595.7191495|174301.828883875||71|000:255:199
734153.6871495|169168.437883875||107|000:080:255
706338.2501495|54889.417883875||129|061:000:255
758009.7501495|112019.898883875||133|087:000:255
754002.7501495|200902.234883875||147|176:000:255
704771.7501495|183364.484883875||191|255:000:052
399187.0631495|220018.859883875||35|031:255:000
685098.9371495|33282.089883875||19|134:255:000
577750.8131495|257153.109883875||1|249:255:000
794095.5621495|199742.671883875||13|172:255:000
634688.2501495|100629.616883875||17|147:255:000
287638.7811495|207582.624883875||21|121:255:000
366218.5321495|222940.625883875||23|108:255:000
385212.4371495|236593.109883875||27|083:255:000
628137.4371495|63995.550883875||47|000:255:046
782600.5631495|152698.890883875||49|000:255:058
502813.9381495|235232.577883875||57|000:255:110
705922.6251495|136589.359883875||61|000:255:135
620397.8131495|246847.640883875||63|000:255:148
738465.3751495|237233.983883875||65|000:255:161
708944.7501495|247632.296883875||127|048:000:255
526666.6871495|249780.312883875||81|000:247:255
733439.3741495|298005.311883875||83|000:234:255
253886.0321495|204702.171883875||87|000:208:255
298337.5621495|178131.233883875||89|000:195:255
788989.3121495|284544.562883875||91|000:182:255
438378.3431495|227157.890883875||97|000:144:255
227406.0621495|190462.640883875||99|000:131:255
668782.8121495|195518.718883875||101|000:118:255
596325.2501495|190688.374883875||105|000:093:255
777034.8131495|229280.030883875||117|000:016:255
337714.9681495|219991.202883875||111|000:054:255
442913.7181495|164364.608883875||119|000:003:255
528779.3741495|179793.890883875||123|023:000:255
571926.0621495|162885.202883875||125|035:000:255
801852.0621495|154744.077883875||137|112:000:255
611371.9381495|295204.467883875||145|164:000:255
534602.9371495|213605.812883875||151|202:000:255
608960.3121495|102612.546883875||155|228:000:255
541667.1251495|293002.124883875||157|240:000:255
475819.8131495|210685.359883875||159|253:000:255
343685.3441495|182222.358883875||161|255:000:244
670160.3121495|140070.249883875||163|255:000:231
497214.4691495|177053.639883875||167|255:000:206
455564.7811495|293568.937883875||171|255:000:180
636209.2501495|210858.014883875||183|255:000:103
413565.1871495|267534.562883875||193|255:000:039
709485.0001495|220209.202883875||195|255:000:026
815182.1881495|110335.429883875||31|057:255:000
495151.3131495|258002.374883875||67|000:255:174
859129.4371495|288849.968883875||139|125:000:255
545200.0621495|127391.304883875||153|215:000:255
655627.6871495|285597.342883875||181|255:000:116
382174.7811495|171687.280883875||45|000:255:033
511235.1561495|135372.358883875||7|211:255:000
332533.5941495|242831.139883875||121|010:000:255
"""
    refrence_cache="""145096.8591495|154534.264883875||39
616341.4371495|146049.750883875||51
410595.7191495|174301.828883875||71
734153.6871495|169168.437883875||107
706338.2501495|54889.417883875||129
758009.7501495|112019.898883875||133
754002.7501495|200902.234883875||147
704771.7501495|183364.484883875||191
399187.0631495|220018.859883875||35
685098.9371495|33282.089883875||19
577750.8131495|257153.109883875||1
794095.5621495|199742.671883875||13
634688.2501495|100629.616883875||17
287638.7811495|207582.624883875||21
366218.5321495|222940.625883875||23
385212.4371495|236593.109883875||27
628137.4371495|63995.550883875||47
782600.5631495|152698.890883875||49
502813.9381495|235232.577883875||57
705922.6251495|136589.359883875||61
620397.8131495|246847.640883875||63
738465.3751495|237233.983883875||65
708944.7501495|247632.296883875||127
526666.6871495|249780.312883875||81
733439.3741495|298005.311883875||83
253886.0321495|204702.171883875||87
298337.5621495|178131.233883875||89
788989.3121495|284544.562883875||91
438378.3431495|227157.890883875||97
227406.0621495|190462.640883875||99
668782.8121495|195518.718883875||101
596325.2501495|190688.374883875||105
777034.8131495|229280.030883875||117
337714.9681495|219991.202883875||111
442913.7181495|164364.608883875||119
528779.3741495|179793.890883875||123
571926.0621495|162885.202883875||125
801852.0621495|154744.077883875||137
611371.9381495|295204.467883875||145
534602.9371495|213605.812883875||151
608960.3121495|102612.546883875||155
541667.1251495|293002.124883875||157
475819.8131495|210685.359883875||159
343685.3441495|182222.358883875||161
670160.3121495|140070.249883875||163
497214.4691495|177053.639883875||167
455564.7811495|293568.937883875||171
636209.2501495|210858.014883875||183
413565.1871495|267534.562883875||193
709485.0001495|220209.202883875||195
815182.1881495|110335.429883875||31
495151.3131495|258002.374883875||67
859129.4371495|288849.968883875||139
545200.0621495|127391.304883875||153
655627.6871495|285597.342883875||181
382174.7811495|171687.280883875||45
511235.1561495|135372.358883875||7
332533.5941495|242831.139883875||121
"""

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule('g.region', raster=cls.map1, flags='p')

    @classmethod
    def tearDownClass(cls):
        cls.del_temp_region()
        if os.path.isfile('result.csv'):
            os.remove('result.csv')

    def test_raster_what_points(self):
        """Testing r.what runs successfully with input coordinates given as a vector points map"""
        module = SimpleModule('r.what', map=self.map1, points=self.points)
        module.run()
        self.assertLooksLike(actual=str(module.outputs.stdout), reference=self.refrence_points,
                             msg="test_raster_what_points did't run successfully")

    def test_raster_what_coordinates(self):
        """Testing r.what runs successfully with input coordinates given as an option and with multiple maps"""
        module = SimpleModule('r.what', map=self.map2, coordinates=self.coordinates)
        module.run()
        self.assertLooksLike(actual=str(module.outputs.stdout), reference=self.refrence_coordinates,
                            msg="test_raster_what_coordinates did't run successfully")


    def test_raster_what_cats(self):
        """Testing r.what runs successfully with input coordinates given as a vector points map with cats and flag v"""
        module = SimpleModule('r.what', map=self.map1, points=self.points, flags='v')
        module.run()
        self.assertLooksLike(actual=str(module.outputs.stdout), reference=self.refrence_cats,
                             msg="test_raster_what_cats did't run successfully")

    def test_raster_what_csv(self):
        """Testing r.what runs successfully with input coordinates given as a vector points map, output into CSV file and flag n"""
        self.assertModule('r.what', map=self.map1, points=self.points, separator='comma', output='result.csv', flags='n')
        self.assertFileExists(filename='result.csv', msg="CSV file was not created")
        if os.path.isfile('result.csv'):
            file = open("result.csv", "r")
            fileData = file.read()
            self.assertLooksLike(actual=fileData, reference=self.refrence_csv,
                                 msg="test_raster_what_csv did't run successfully")
            file.close()

    def test_raster_what_points_flag_i(self):
        """Testing r.what runs successfully with flag i"""
        module = SimpleModule('r.what', map=self.map1, points=self.points, flags='i')
        module.run()
        self.assertLooksLike(actual=str(module.outputs.stdout), reference=self.refrence_flag_i,
                             msg="test_raster_what_cats did't run successfully")


    def test_raster_what_points_flag_f(self):
        """Testing r.what runs successfully with flag f"""
        module = SimpleModule('r.what', map=self.map1, points=self.points, flags='f')
        module.run()
        self.assertLooksLike(actual=str(module.outputs.stdout), reference=self.refrence_flag_f,
                             msg="test_raster_what_cats did't run successfully")

    def test_raster_what_points_flag_r(self):
        """Testing r.what runs successfully with flag r"""
        module = SimpleModule('r.what', map=self.map1, points=self.points, flags='r')
        module.run()
        self.assertLooksLike(actual=str(module.outputs.stdout), reference=self.refrence_flag_r,
                             msg="test_raster_what_cats did't run successfully")

    def test_raster_what_cache(self):
        """Testing r.what runs successfully with cache"""
        module = SimpleModule('r.what', map=self.map1, points=self.points, flags='c', cache=500)
        module.run()
        self.assertLooksLike(actual=str(module.outputs.stdout), reference=self.refrence_cache,
                             msg="test_raster_what_cats did't run successfully")


if __name__ == '__main__':
    from grass.gunittest.main import test
    test()