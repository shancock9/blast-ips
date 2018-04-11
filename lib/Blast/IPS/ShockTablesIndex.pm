package Blast::IPS::ShockTablesIndex;

# MIT License
# Copyright (c) 2018 Steven Hancock
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

use strict;
use warnings;
use 5.006;

use Carp;
our $rtables_info;

##################################################################
# Index to Shock Front Tables
##################################################################

BEGIN {

    # Created with: ./table_data_to_hash.pl
    # Wed Apr 11 06:11:56 2018   Inspiron-3668

    # symmetry = 0,1,2 for Plane, Cylindrical, Spherical
    # gamma = ideal gas gamma
    # Max-Error = the sum of max absolute FD + MOC + interpolation errors;
    # N=number of points from center to shock in Finite Difference run
    # FD-Error = maximum error in FD run by comparing with N/2 run
    # Interp-Error = maximum cubic interpolation error in overpressure
    # MOC-Error = maximum error in MOC run
    # Energy-Error = energy error in the FD run at 0.1 shock overpressure ratio
    # rs2 = radius of formation of second shock
    # zs2 = value of z for formation of second shock

    # table_name => [symmetry, gamma, Max-Error, N, Energy-Error, R-FD-MOC,
    # FD-Error, MOC-Error, Interp-Error, rs2, zs2],

    $rtables_info = {
        'S1.1' => [
            2, 1.1, 1.154994779502965E-06, 8000, 2.03E-07, 10, 3E-07,
            3.54994779502965E-07, 5E-07, 1590.90000697035, -0.430996019022666,
        ],
        'S1.15' => [
            2,        1.15,    1.036580500502762E-06, 8000,
            2.01E-07, 10,      4.74280500502762E-07,  6.23E-08,
            5E-07,    1817.49, -0.486164267299174,
        ],
        'S1.17' => [
            2,                      1.17,
            1.2842848916807499E-06, 8000,
            2E-07,                  10,
            3.30684504964059E-07,   4.53600386716691E-07,
            5E-07,                  1885.23504981054,
            -0.503986001046336,
        ],
        'S1.2' => [
            2, 1.2, 1.0799999999999998E-06, 8000, 2E-07, 10, 2.9E-07, 2.9E-07,
            5E-07, 1966.9304550491, -0.527578104660463,
        ],
        'S1.23' => [
            2,                     1.23,
            1.147112976355249E-06, 8000,
            2E-07,                 10,
            3.16886598146956E-07,  3.30226378208293E-07,
            5E-07,                 2027.77994386824,
            -0.548197622950195,
        ],
        'S1.25' => [
            2, 1.25, 1.1E-06, 8000, 2E-07, 10, 3.00445937758316E-07, 3.19E-07,
            5E-07, 2058.37246962197, -0.560613036407896,
        ],
        'S1.3' => [
            2, 1.3, 8.9E-07, 8000, 1.99E-07, 10, 2.2E-07, 1.7E-07, 5E-07,
            2103.19830388092, -0.587909502251664,
        ],
        'S1.35' => [
            2, 1.35, 1.014458473371616E-06, 8000, 2E-07, 10,
            1.91458473371616E-07, 3.23E-07, 5E-07, 2111.32990152372,
            -0.611012207905581,
        ],
        'S1.4' => [
            2, 1.4, 7.21E-07, 32000, 1.27E-08, 10, 5.4E-08, 1.67E-07, 5E-07,
            2090.536164577, -0.630899471186397,
        ],
        'S1.45' => [
            2, 1.45, 9.03179743603947E-07, 8000, 2.06E-07, 10,
            1.50179743603947E-07, 2.53E-07, 5E-07, 2047.39655039366,
            -0.648234950678994,
        ],
        'S1.5' => [
            2, 1.5, 9.714139578770349E-07, 8000, 2.08E-07, 10,
            1.41413957877035E-07, 3.3E-07, 5E-07, 1987.48025731162,
            -0.663499435733238,
        ],
        'S1.55' => [
            2,                     1.55,
            8.926474931760189E-07, 8000,
            2.1E-07,               10,
            1.29761083063803E-07,  2.62886410112216E-07,
            5E-07,                 1915.39166708926,
            -0.677046343668842,
        ],
        'S1.6' => [
            2, 1.6, 1.0884682653590609E-06, 8000, 2.1E-07, 10,
            1.28468265359061E-07, 4.6E-07, 5E-07, 1834.82586890454,
            -0.689144044780536,
        ],
        'S1.65' => [
            2, 1.65, 7.52645137387977E-07, 8000, 2.13E-07, 10,
            1.23645137387977E-07, 1.29E-07, 5E-07, 1749.07382559716,
            -0.700011158984465,
        ],
        'S1.667' => [
            2, 1.667, 8.2E-07, 8000, 2.14E-07, 10, 1.2E-07, 2E-07, 5E-07,
            1718.96873720766, -0.70344977418314,
        ],
        'S1.7' => [
            2,                    1.7,
            7.95180307735647E-07, 8000,
            2.15E-07,             10,
            1.18475757204806E-07, 1.76704550530841E-07,
            5E-07,                1660.14437065796,
            -0.709802810940815,
        ],
        'S1.8' => [
            2, 1.8, 7.544054782890239E-07, 8000, 2.17E-07, 10,
            1.07405478289024E-07, 1.47E-07, 5E-07, 1481.26775852019,
            -0.726723649148212,
        ],
        'S1.9' => [
            2,                     1.9,
            9.882225872339029E-07, 4000,
            8.8E-07,               10,
            3.71226335005304E-07,  1.16996252228599E-07,
            5E-07,                 1308.7083462168,
            -0.740734256395012,
        ],
        'S2' => [
            2, 2, 8.27E-07, 8000, 2.2E-07, 10, 9.7E-08, 2.3E-07, 5E-07,
            1148.4754981706, -0.75244193378212,
        ],
        'S2.2' => [
            2,                    2.2,
            9.94869976636274E-07, 4000,
            8.96E-07,             10,
            3.16415492918054E-07, 1.7845448371822E-07,
            5E-07,                874.4323272327,
            -0.770582210410918,
        ],
        'S2.5' => [
            2, 2.5, 9.82047301331223E-07, 4000, 9.1E-07, 10,
            2.78047301331223E-07, 2.04E-07, 5E-07, 574.687329593191,
            -0.788451647385764,
        ],
        'S3' => [
            2, 3, 1.4435908072865733E-06, 8000, 2.3E-07, 10,
            6.35908072865732E-08, 8.8E-07, 5E-07, 289.43005770563,
            -0.803012949195847,
        ],
        'S4' => [
            2, 4, 9.1E-07, 4000, 9.36E-07, 10, 1.8E-07, 2.3E-07, 5E-07,
            197.699960698537, -0.904330236782954,
        ],
        'S5' => [
            2, 5, 1.1250487205292939E-06, 4000, 9.2E-07, 10,
            1.55048720529294E-07, 4.7E-07, 5E-07, 27.6123237136561,
            -0.781573364765684,
        ],
        'S6' => [
            2, 6, 9.97E-07, 4000, 9.53E-07, 10, 1.4E-07, 3.57E-07, 5E-07,
            8.90793490986525, -0.734405840270103,
        ],
        'S6.5' => [
            2,                     6.5,
            8.267675831824938E-07, 8000,
            3.78E-07,              4,
            3.64124647470968E-08,  2.90355118435397E-07,
            5E-07,                 0.59027593350696,
            -0.838579307321098,
        ],
        'C1.1' => [
            1, 1.1, 1.756333057406828E-06, 8000, 2.1E-07, 20,
            7.01333057406828E-07, 5.55E-07, 5E-07, 31.3641162051672,
            -0.581240059301524,
        ],
        'C1.15' => [
            1, 1.15, 1.061E-06, 8000, 1.92E-07, 30, 3.23E-07, 2.38E-07, 5E-07,
            39.6343014919069, -0.700087571960808,
        ],
        'C1.17' => [
            1, 1.17, 1.3551196244667319E-06, 4000, 7.48E-07, 30,
            7.05119624466732E-07, 1.5E-07, 5E-07, 42.8447595772785,
            -0.74161219858533,
        ],
        'C1.2' => [
            1, 1.2, 1.032E-06, 8000, 1.8E-07, 30, 1.82E-07, 3.5E-07, 5E-07,
            46.613364650253, -0.794904822764063,
        ],
        'C1.23' => [
            1, 1.23, 1.619339279001318E-06, 4000, 7.2E-07, 30, 6.7E-07,
            4.49339279001318E-07, 5E-07, 50.6639953859758, -0.845290182174436,
        ],
        'C1.25' => [
            1, 1.25, 2.0100000000000002E-06, 8000, 1.82E-07, 30, 1.1E-06,
            4.1E-07, 5E-07, 52.9761489451367, -0.87534673773403,
        ],
        'C1.3' => [
            1, 1.3, 7.66E-07, 8000, 1.87E-07, 50, 1.6E-07, 1.06E-07, 5E-07,
            59.8246075603305, -0.949926864449608,
        ],
        'C1.35' => [
            1, 1.35, 8.71E-07, 8000, 1.92E-07, 30, 1.53E-07, 2.18E-07, 5E-07,
            64.0707364752969, -1.0069747620868,
        ],
        'C1.4' => [
            1, 1.4, 7.31E-07, 16000, 4.91E-08, 50, 5.2E-08, 1.79E-07, 5E-07,
            69.0122514583704, -1.06269039232466,
        ],
        'C1.45' => [
            1, 1.45, 7E-07, 8000, 2.01E-07, 30, 1.44E-07, 5.2E-08, 5E-07,
            73.3886943631669, -1.112455227271,
        ],
        'C1.49' => [
            1, 1.49, 1.383785845268209E-06, 4000, 8.21E-07, 30,
            5.63785845268209E-07, 3.2E-07, 5E-07, 76.7481691015437,
            -1.14958050865576,
        ],
        'C1.51' => [
            1,                      1.51,
            1.1793805256347179E-06, 4000,
            8.29E-07,               30,
            5.57999273749961E-07,   1.21381251884757E-07,
            5E-07,                  78.3155592640345,
            -1.16709023313781,
        ],
        'C1.55' => [
            1, 1.55, 1.138419954658165E-06, 4000, 8.44E-07, 50,
            5.46919954658165E-07, 9.15E-08, 5E-07, 81.2418837140749,
            -1.20023846000998,
        ],
        'C1.6' => [
            1, 1.6, 7.6341154151784E-07, 8000, 2.16E-07, 50,
            1.3341154151784E-07, 1.3E-07, 5E-07, 84.6819995308298,
            -1.23901139081572,
        ],
        'C1.65' => [
            1, 1.65, 1.089114785882936E-06, 4000, 8.82E-07, 50,
            5.21614785882936E-07, 6.75E-08, 5E-07, 87.8372295538845,
            -1.27498816352225,
        ],
        'C1.667' => [
            1, 1.667, 6.839999999999999E-07, 8000, 2.2E-07, 50, 1.3E-07,
            5.4E-08, 5E-07, 88.8581044049173, -1.28669743384461,
        ],
        'C1.7' => [
            1,                      1.7,
            1.0849647766184396E-06, 4000,
            9E-07,                  70,
            5.10090997292897E-07,   7.48737793255427E-08,
            5E-07,                  90.7333120627518,
            -1.30853956435519,
        ],
        'C1.8' => [
            1, 1.8, 6.80147636745969E-07, 8000, 2.34E-07, 100,
            1.22147636745969E-07, 5.8E-08, 5E-07, 95.8077523163297,
            -1.3693463552985,
        ],
        'C1.9' => [
            1,                    1.9,
            1.09351540408377E-06, 4000,
            9.6E-07,              100,
            4.69863059127107E-07, 1.23652344956663E-07,
            5E-07,                100.023132017295,
            -1.42305839115265,
        ],
        'C2' => [
            1, 2, 8E-07, 8000, 2.49E-07, 42, 1.1E-07, 1.9E-07, 5E-07,
            103.506634318208, -1.47101822583144,
        ],
        'C2.2' => [
            1,                     2.2,
            1.035102755996909E-06, 4000,
            1.05E-06,              70,
            4.23455640774506E-07,  1.11647115222403E-07,
            5E-07,                 108.654722068711,
            -1.55319919357162,
        ],
        'C2.5' => [
            1, 2.5, 9.28E-07, 4000, 1.12E-06, 61, 3.9E-07, 3.8E-08, 5E-07,
            113.003459332856, -1.65149812975875,
        ],
        'C3' => [
            1, 3, 1.5879999999999999E-06, 8000, 3E-07, 16, 8.8E-08, 1E-06,
            5E-07, 114.538489804021, -1.77183933353766,
        ],
        'C4' => [
            1, 4, 8.731488580208149E-07, 4000, 1.31E-06, 100,
            3.12148858020815E-07, 6.1E-08, 5E-07, 108.42791357809,
            -1.9262887758072,
        ],
        'C5' => [
            1, 5, 8.949999999999999E-07, 4000, 1.38E-06, 100, 2.95E-07, 1E-07,
            5E-07, 99.1523647101572, -2.02460296710399,
        ],
        'C6' => [
            1, 6, 9E-07, 4000, 1.43E-06, 20, 2.9E-07, 1.1E-07, 5E-07,
            90.2350008763815, -2.09612315913494,
        ],
        'C7' => [
            1, 7, 1.604577258877604E-06, 4000, 1.47E-06, 70,
            2.85577258877604E-07, 3.19E-07, 1E-06, 82.4658386579933,
            -2.15278923435504,
        ],
        'P1.1' => [
            0, 1.1, 9.886974120660998E-07, 8000, 2.5E-07, 24,
            1.886974120661E-07, 3E-07, 5E-07, 118.190649621423,
            -1.18989590916949,
        ],
        'P1.15' => [
            0,                     1.15,
            1.318377740495394E-06, 8000,
            6.83E-07,              100,
            4.99786750207676E-07,  3.18590990287718E-07,
            5E-07,                 308.535316078086,
            -2.00236300685825,
        ],
        'P1.17' => [
            0,                     1.17,
            1.216073200965371E-06, 8000,
            6.7E-07,               100,
            4.88491076626E-07,     2.27582124339371E-07,
            5E-07,                 441.457483711759,
            -2.39249086681047,
        ],
        'P1.2' => [
            0, 1.2, 1.0299999999999999E-06, 8000, 2.4E-07, 30, 1.8E-07,
            3.5E-07, 5E-07, 750.657922315278, -3.07511470169219,
        ],
        'P1.23' => [
            0, 1.23, 1.251706363807934E-06, 8000, 6.4E-07, 100,
            4.61706363807934E-07, 2.9E-07, 5E-07, 1280.09755807535,
            -3.91009051398649,
        ],
        'P1.25' => [
            0, 1.25, 1.274498152742862E-06, 8000, 6.31E-07, 100,
            4.54498152742862E-07, 3.2E-07, 5E-07, 1835.42724920336,
            -4.57305612734795,
        ],
        'P1.3' => [
            0, 1.3, 1.11E-06, 8000, 2.32E-07, 34, 1.7E-07, 4.4E-07, 5E-07,
            4675.19951552032, -6.75797526138565,
        ],
        'P1.35' => [
            0, 1.35, 1.409060877689594E-06, 8000, 6E-07, 100,
            4.26060877689594E-07, 4.83E-07, 5E-07, 12870.837447165,
            -10.1225394232029,
        ],
        'P1.4' => [
            0, 1.4, 1.161E-06, 8000, 2.24E-07, 20, 2.75E-07, 3.86E-07, 5E-07,
            39967.7084262245, -15.6489715936575,
        ],
        'P1.45' => [
            0, 1.45, 1.475105999856659E-06, 8000, 5.69E-07, 100,
            4.05105999856659E-07, 5.7E-07, 5E-07, 149364.356044368,
            -25.5765593976355,
        ],
        'P1.5' => [
            0, 1.5, 1.736279161337295E-06, 8000, 5.58E-07, 100,
            3.96279161337295E-07, 8.4E-07, 5E-07, 752849.002018866,
            -45.9835612809612,
        ],
        'P1.55' => [
            0, 1.55, 1.47824097699693E-06, 8000, 5.48E-07, 100,
            3.8824097699693E-07, 5.9E-07, 5E-07, 6491381.24182815,
            -98.6557492130365,
        ],
        'P1.6' => [
            0,                     1.6,
            1.673191767506391E-06, 8000,
            5.4E-07,               100,
            3.80844752179144E-07,  7.92347015327247E-07,
            5E-07,                 570498850.647285,
            -768.593235782215,
        ],
        'P1.65' => [
            0, 1.65, 1.454023583337079E-06, 8000, 5.3E-07, 100,
            3.73986487852325E-07, 5.80037095484754E-07, 5E-07,,,
        ],
        'P1.667' => [
            0,        1.667,   1.127E-06, 8000, 2.1E-07, 25,
            1.57E-07, 4.7E-07, 5E-07,,,
        ],
        'P1.7' => [
            0, 1.7, 4.43702740234804E-06, 4000, 2.1E-06, 100,
            1.47020281149324E-06, 2.4668245908548E-06, 5E-07,,,
        ],
        'P1.8' => [
            0, 1.8, 1.42592736535861E-06, 8000, 5.05E-07, 100,
            3.5592736535861E-07, 5.7E-07, 5E-07,,,
        ],
        'P1.9' => [
            0, 1.9, 1.625623111464022E-06, 8000, 4.91E-07, 100,
            3.45623111464022E-07, 7.8E-07, 5E-07,,,
        ],
        'P2' => [ 0, 2, 1.33E-06, 8000, 4.8E-07, 58, 3.3E-07, 5E-07, 5E-07,,, ],
        'P2.2' => [
            0, 2.2, 3.67252548579766E-06, 4000, 1.8E-06, 100,
            1.27993624943201E-06, 1.89258923636565E-06, 5E-07,,,
        ],
        'P2.5' =>
          [ 0, 2.5, 2.07E-06, 4000, 1.72E-06, 100, 1.2E-06, 3.7E-07, 5E-07,,, ],
        'P3' =>
          [ 0, 3, 2.67E-06, 8000, 4E-07, 100, 2.7E-07, 1.9E-06, 5E-07,,, ],
        'P4' =>
          [ 0, 4, 3.35E-06, 4000, 1.42E-06, 100, 9.7E-07, 1.88E-06, 5E-07,,, ],
        'P5' => [
            0, 5, 1.665226161410439E-06, 4000, 1.31E-06, 100,
            8.85226161410439E-07, 2.8E-07, 5E-07,,,
        ],
        'P6' =>
          [ 0, 6, 1.52E-06, 4000, 1.23E-06, 100, 8.3E-07, 1.9E-07, 5E-07,,, ],
        'P7' => [
            0, 7, 2.530622424541362E-06, 4000, 5.1E-07, 500,
            3.80622424541362E-07, 1.15E-06, 1E-06,,,
        ],
    };
}
1;
