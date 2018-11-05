package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='C', gamma=5.5
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'C5x5'} = {

    table_name  => 'C5x5',
    symmetry    => 1,
    gamma       => 5.5,
    data_source => 'C4000_G5x5/moc2_from_moc_r5000/',

    shock_table_info => [ 1, 5.5, 1.02e-06, 4000, 1.4e-06, 50, 2.9e-07, 2.3e-07, 5e-07, 94.584, -2.0629 ],

    shock_table => [
        [ -5.71714124, 12.00033763,  -1.99997921, -5.71875456, 0.9983854 ],
        [ -5.24979804, 11.06566488,  -1.99994707, -5.2523737,  0.99742108 ],
        [ -4.85075367, 10.26760345,  -1.99990159, -4.85459481, 0.99615168 ],
        [ -4.29339035, 9.15297826,   -1.99969917, -4.3001065,  0.99326236 ],
        [ -3.82739073, 8.22121111,   -1.99923589, -3.83811369, 0.98922385 ],
        [ -3.44133434, 7.44954351,   -1.9983486,  -3.45714586, 0.98407724 ],
        [ -3.14283482, 6.85321931,   -1.99699965, -3.1641973,  0.97844318 ],
        [ -2.86471299, 6.29809008,   -1.99478456, -2.89300553, 0.97138544 ],
        [ -2.61682179, 5.80397736,   -1.9914774,  -2.65318611, 0.96313954 ],
        [ -2.39266724, 5.35807052,   -1.98674685, -2.4383187,  0.95363055 ],
        [ -2.18838053, 4.95282569,   -1.98024437, -2.24456924, 0.94283129 ],
        [ -1.9987113,  4.57801511,   -1.9714955,  -2.06686811, 0.93057671 ],
        [ -1.82023859, 4.22713159,   -1.95995318, -1.90198664, 0.91670773 ],
        [ -1.64937953, 3.89347896,   -1.94489098, -1.7466675,  0.90095728 ],
        [ -1.48224538, 3.56998515,   -1.92527802, -1.5975577,  0.88289296 ],
        [ -1.3137409,  3.24764887,   -1.89946264, -1.45052503, 0.86174516 ],
        [ -1.13321078, 2.90782315,   -1.86381785, -1.29725847, 0.83560536 ],
        [ -0.95384236, 2.57736973,   -1.81924145, -1.14997712, 0.8060267 ],
        [ -0.79044178, 2.28398817,   -1.77042086, -1.02067777, 0.77613973 ],
        [ -0.63761244, 2.01733671,   -1.71811069, -0.90433931, 0.74601107 ],
        [ -0.49246337, 1.77187209,   -1.66341779, -0.79822647, 0.71591401 ],
        [ -0.35331604, 1.54426373,   -1.60761808, -0.70066895, 0.68621002 ],
        [ -0.21552997, 1.32667978,   -1.55047019, -0.60816554, 0.65649325 ],
        [ -0.07736477, 1.11645472,   -1.49265314, -0.51951061, 0.626899 ],
        [ 0.06158086,  0.91305853,   -1.43528287, -0.43443682, 0.59781131 ],
        [ 0.20288109,  0.71426792,   -1.37885149, -0.35199333, 0.56932888 ],
        [ 0.34790571,  0.51833029,   -1.32382465, -0.27145884, 0.54157467 ],
        [ 0.49793922,  0.32375409,   -1.27061145, -0.1922464,  0.51467718 ],
        [ 0.6542428,   0.12920257,   -1.21955464, -0.11385436, 0.48876146 ],
        [ 0.81832434,  -0.06684096,  -1.17086437, -0.03572924, 0.46391056 ],
        [ 0.99164664,  -0.26570382,  -1.12473768, 0.04258666,  0.44022073 ],
        [ 1.17564437,  -0.46857701,  -1.0813502,  0.12148154,  0.41779285 ],
        [ 1.37242964,  -0.67728022,  -1.04071412, 0.20157122,  0.39665402 ],
        [ 1.58380111,  -0.89316431,  -1.00292851, 0.28327277,  0.37688314 ],
        [ 1.80497132,  -1.11113341,  -0.96900903, 0.36460567,  0.35904679 ],
        [ 2.03436269,  -1.32988041,  -0.93898955, 0.44510275,  0.34320097 ],
        [ 2.27381326,  -1.55144682,  -0.91236692, 0.52554978,  0.32911287 ],
        [ 2.52517207,  -1.77772425,  -0.88873748, 0.6066574,   0.31659507 ],
        [ 2.79000875,  -2.01024052,  -0.86779821, 0.68899187,  0.30550786 ],
        [ 3.07007718,  -2.250612,    -0.84928329, 0.77314298,  0.2957268 ],
        [ 3.36700919,  -2.50029431,  -0.83297778, 0.85963917,  0.28715066 ],
        [ 3.68245419,  -2.76072628,  -0.81869437, 0.9490024,   0.27968952 ],
        [ 4.01835668,  -3.03356722,  -0.80625669, 1.04183121,  0.27325629 ],
        [ 4.37668204,  -3.32047471,  -0.79550978, 1.13872542,  0.26777227 ],
        [ 4.75934949,  -3.62306421,  -0.78631377, 1.24027493,  0.26316384 ],
        [ 5.16908717,  -3.94358629,  -0.77852253, 1.34728715,  0.25935175 ],
        [ 5.60800795,  -4.28380645,  -0.77201334, 1.46041179,  0.25626599 ],
        [ 6.07962117,  -4.6465754,   -0.76664939, 1.58066429,  0.25382763 ],
        [ 6.58762273,  -5.0348765,   -0.76230468, 1.70910651,  0.25196143 ],
        [ 7.13669776,  -5.45244103,  -0.75885428, 1.84705018,  0.25059156 ],
        [ 7.73271972,  -5.90389025,  -0.75617667, 1.99610259,  0.24964338 ],
        [ 8.38435124,  -6.39593644,  -0.7541516,  2.15856264,  0.24904363 ],
        [ 9.1018428,   -6.93646401,  -0.75266956, 2.33711809,  0.24872433 ],
        [ 9.90063554,  -7.53724018,  -0.75162529, 2.535743,    0.24862139 ],
        [ 10.80416167, -8.21600651,  -0.75092272, 2.76039475,  0.24867719 ],
        [ 11.84704569, -8.99887483,  -0.75047788, 3.01981654,  0.24884111 ],
        [ 13.08570621, -9.92828046,  -0.75021806, 3.32818747,  0.24906942 ],
        [ 14.61815771, -11.07783294, -0.75008268, 3.71007643,  0.24932481 ],
        [ 16.64141307, -12.59536912, -0.75002327, 4.21479479,  0.24957579 ],
        [ 19.63828919, -14.84305807, -0.75000374, 4.96310748,  0.24979448 ],
        [ 25.37763281, -19.14757196, -0.75000014, 6.39731285,  0.24995053 ]
    ],

    energy_table => [
        [ -5.71714124, 4.909e-05,  7.758e-05,  4.909e-05,  7.758e-05,     6.0965651,   -1.1425547 ],
        [ -5.24979804, 0.00010246, 0.00016066, 0.00010246, 0.00016066,    5.5603786,   -1.1526299 ],
        [ -4.85075367, 0.00019108, 0.00029714, 0.00019108, 0.00029714,    5.0986147,   -1.1625164 ],
        [ -4.29339035, 0.0004518,  0.00069201, 0.0004518,  0.00069201,    4.4465118,   -1.1786191 ],
        [ -3.82739073, 0.00091708, 0.00138076, 0.00091708, 0.00138076,    3.8939091,   -1.1939963 ],
        [ -3.44133434, 0.00163153, 0.00241093, 0.00163153, 0.00241093,    3.4303518,   -1.2082432 ],
        [ -3.14283482, 0.00252627, 0.00366532, 0.00252627, 0.00366532,    3.0679624,   -1.2202239 ],
        [ -2.86471299, 0.00376663, 0.00535066, 0.00376663, 0.00535066,    2.7269882,   -1.2318667 ],
        [ -2.61682179, 0.00533565, 0.00740607, 0.00533565, 0.00740607,    2.4203206,   -1.2423339 ],
        [ -2.39266724, 0.00725539, 0.00981786, 0.00725539, 0.00981786,    2.1407861,   -1.251595 ],
        [ -2.18838053, 0.00953032, 0.01254173, 0.00953032, 0.01254173,    1.8842578,   -1.2596248 ],
        [ -1.9987113,  0.01218709, 0.01555074, 0.01218709, 0.01555074,    1.6446593,   -1.266598 ],
        [ -1.82023859, 0.01524644, 0.01879734, 0.01524644, 0.01879734,    1.4180428,   -1.2726607 ],
        [ -1.64937953, 0.01874787, 0.02223668, 0.01874787, 0.02223668,    1.2001218,   -1.2780155 ],
        [ -1.48224538, 0.02276209, 0.02582653, 0.02276209, 0.02582653,    0.98610124,  -1.2828889 ],
        [ -1.3137409,  0.02742627, 0.02953372, 0.02742627, 0.02953372,    0.76952826,  -1.2875248 ],
        [ -1.13321078, 0.03310976, 0.03338951, 0.03310976, 0.03338951,    0.53665424,  -1.2923389 ],
        [ -0.95384236, 0.03941558, 0.0368338,  0.03941558, 0.0368338,     0.30442292,  -1.2972439 ],
        [ -0.79044178, 0.0456527,  0.03939912, 0.0456527,  0.03939912,    0.092076347, -1.3022117 ],
        [ -0.63761244, 0.05181589, 0.04113869, 0.05181589, 0.04113869,    -0.10732086, -1.3077643 ],
        [ -0.49246337, 0.05786592, 0.04210953, 0.05786592, 0.04210953,    -0.29756337, -1.3142665 ],
        [ -0.35331604, 0.06375222, 0.04239029, 0.06375222, 0.04239029,    -0.48091981, -1.32189 ],
        [ -0.21552997, 0.06957721, 0.04206516, 0.06957721, 0.04206516,    -0.66362755, -1.3308527 ],
        [ -0.07736477, 0.07533512, 0.04120004, 0.07533512, 0.04120004,    -0.84817397, -1.3411123 ],
        [ 0.06158086,  0.08097247, 0.03987775, 0.08097247, 0.03987775,    -1.0352728,  -1.352409 ],
        [ 0.20288109,  0.08649065, 0.0381777,  0.08649065, 0.0381777,     -1.2272088,  -1.3644918 ],
        [ 0.34790571,  0.09188477, 0.03617773, 0.09188477, 0.03617773,    -1.426007,   -1.3770468 ],
        [ 0.49793922,  0.09714708, 0.03395337, 0.09714708, 0.03395337,    -1.6335814,  -1.3897565 ],
        [ 0.6542428,   0.10226852, 0.03157649, 0.10226852, 0.03157649,    -1.8518208,  -1.4023253 ],
        [ 0.81832434,  0.10724638, 0.02911041, 0.10724638, 0.02911041,    -2.0829624,  -1.4144788 ],
        [ 0.99164664,  0.11207348, 0.02661398, 0.11207348, 0.02661398,    -2.3291812,  -1.4260243 ],
        [ 1.17564437,  0.11673977, 0.02414064, 0.11673977, 0.02414064,    -2.5926292,  -1.4368022 ],
        [ 1.37242964,  0.12124892, 0.02172913, 0.12124892, 0.02172913,    -2.8764234,  -1.4466943 ],
        [ 1.58380111,  0.12559246, 0.01941705, 0.12559246, 0.01941705,    -3.1832446,  -1.4556408 ],
        [ 1.80497132,  0.12964678, 0.01729403, 0.12964678, 0.01729403,    -3.5061303,  -1.4633921 ],
        [ 2.03436269,  0.13338946, 0.01538409, 0.13338946, 0.01538409,    -3.8426514,  -1.4699653 ],
        [ 2.27381326,  0.13686262, 0.01366945, 0.13686262, 0.01366945,    -4.1953723,  -1.4755144 ],
        [ 2.52517207,  0.14010018, 0.01213226, 0.14010018, 0.01213226,    -4.5669081,  -1.480199 ],
        [ 2.79000875,  0.14312609, 0.01075691, 0.14312609, 0.01075691,    -4.9595022,  -1.4841292 ],
        [ 3.07007718,  0.14596172, 0.00952748, 0.14596172, 0.00952748,    -5.3756718,  -1.4873934 ],
        [ 3.36700919,  0.14862304, 0.00842963, 0.14862304, 0.00842963,    -5.8177787,  -1.4900864 ],
        [ 3.68245419,  0.15112308, 0.00744988, 0.15112308, 0.00744988,    -6.2882111,  -1.4923026 ],
        [ 4.01835668,  0.15347426, 0.00657503, 0.15347426, 0.00657503,    -6.7898274,  -1.494105 ],
        [ 4.37668204,  0.15568606, 0.00579328, 0.15568606, 0.00579328,    -7.325496,   -1.4955499 ],
        [ 4.75934949,  0.15776525, 0.00509425, 0.15776525, 0.00509425,    -7.8980506,  -1.4967004 ],
        [ 5.16908717,  0.15972036, 0.0044675,  0.15972036, 0.0044675,     -8.511516,   -1.4975881 ],
        [ 5.60800795,  0.16155418, 0.00390515, 0.16160918, 0.0041363082,  -9.169014,   -1.4982805 ],
        [ 6.07962117,  0.16327301, 0.00339892, 0.16347896, 0.0038124059,  -9.8757712,  -1.4988043 ],
        [ 6.58762273,  0.16488023, 0.00294221, 0.16535976, 0.0035838227,  -10.637281,  -1.4991846 ],
        [ 7.13669776,  0.16637897, 0.00252929, 0.16720946, 0.0031678586,  -11.460537,  -1.4994615 ],
        [ 7.73271972,  0.16777164, 0.00215537, 0.16898453, 0.0027967826,  -12.35432,   -1.4996564 ],
        [ 8.38435124,  0.16906211, 0.00181607, 0.17069187, 0.0024466566,  -13.331598,  -1.4997912 ],
        [ 9.1018428,   0.17025112, 0.0015085,  0.17232345, 0.0021085216,  -14.407726,  -1.499881 ],
        [ 9.90063554,  0.1713409,  0.00122999, 0.17387762, 0.0017910553,  -15.605849,  -1.4999414 ],
        [ 10.80416167, 0.17233414, 0.00097831, 0.17534179, 0.0014504327,  -16.96111,   -1.4999828 ],
        [ 11.84704569, 0.17323146, 0.00075298, 0.17668867, 0.0011595041,  -18.525436,  -1.5000139 ],
        [ 13.08570621, 0.1740326,  0.00055148, 0.17797804, 0.00095005972, -20.383465,  -1.5000504 ],
        [ 14.61815771, 0.17473442, 0.00037572, 0.17920057, 0.00065461883, -22.682256,  -1.5001331 ],
        [ 16.64141307, 0.17533095, 0.0002265,  0.18023987, 0.00039461359, -25.717565,  -1.500249 ],
        [ 19.63828919, 0.17580863, 0.00010707, 0.18107208, 0.00018652738, -30.213797,  -1.5 ],
        [ 25.37763281, 0.1761349,  2.55e-05,   0.1816405,  4.4421409e-05, -38.822812,  -1.5 ]
    ],

    impulse_table => [
        [
            -6.6245081, 13.815063,  -6.6251589,    0.060125017, 0,          -0.16943448,
            0,          -2.2628374, -0.0097080894, -0.51850636, 0.43014202, 0.99993527,
            -5.78212
        ],
        [
            -6.5022893, 13.570625,  -6.5030247,   0.062928744, 0,          -0.169262,
            0,          -2.2628405, -0.010319846, -0.51833379, 0.43014184, 0.99992154,
            -5.6918431
        ],
        [
            -6.2146072, 12.995263,  -6.2155879,   0.069986523, 0,          -0.16876234,
            0,          -2.2628502, -0.011916329, -0.51783379, 0.43014139, 0.99987647,
            -5.4806255
        ],
        [
            -5.8091421, 12.184338,  -5.8106135,   0.081094286, 0,          -0.16776348,
            0,          -2.2628766, -0.014594453, -0.51683379, 0.43014027, 0.99976573,
            -5.1862866
        ],
        [
            -5.2983164, 11.162699, -5.3007699,   0.097176857, 0,         -0.16576798,
            0,          -2.262959, -0.018841305, -0.51483379, 0.4301368, 0.99947696,
            -4.8219492
        ],
        [
            -4.6051693, 9.7764657,  -4.6100809,   0.12301449,  0,          -0.16079604,
            0,          -2.2633298, -0.026645168, -0.50983378, 0.43012068, 0.99846124,
            -4.3419066
        ],
        [
            -3.9120221, 8.3904141, -3.9218705,   0.15362433,  0,          -0.15094861,
            0,          -2.264725, -0.037678408, -0.49983378, 0.43005634, 0.99555949,
            -3.8836588
        ],
        [
            -3.506557, 7.5798876,  -3.5213633,   0.17363901,  0,          -0.14126471,
            0,         -2.2669167, -0.046137424, -0.49020204, 0.42994922, 0.99185065,
            -3.6290552
        ],
        [
            -2.9957314, 6.5595247,  -3.020514,    0.20077407, 0,          -0.12250708,
            0,          -2.2734307, -0.059517295, -0.4713082, 0.42960722, 0.98281162,
            -3.3271629
        ],
        [
            -2.6592591, 5.8885049,  -2.6940925,   0.21961478, 0,          -0.10470217,
            0,          -2.2824617, -0.070325382, -0.4533415, 0.42909644, 0.97231424,
            -3.1430465
        ],
        [
            -2.3025842, 5.1792097, -2.3526114,   0.24026738,  0,          -0.080003656,
            0,          -2.300078, -0.083773247, -0.42836166, 0.42801991, 0.95488144,
            -2.9643827
        ],
        [
            -2.0402199, 4.6598949,  -2.1055553,   0.25586066,  0,          -0.057829368,
            0,          -2.3218816, -0.095033304, -0.40587602, 0.42658222, 0.93627249,
            -2.8464101
        ],
        [
            -1.8971191, 4.3780268,  -1.9727077,  0.26450783,  0,          -0.044429611,
            0,          -2.3384236, -0.10163696, -0.39233621, 0.42543092, 0.92352202,
            -2.7877845
        ],
        [
            -1.609437, 3.8158779, -1.7107616,  0.28220898,  0,          -0.015390277,
            0,         -2.385644, -0.11568621, -0.36350595, 0.42191896, 0.89130128,
            -2.6842666
        ],
        [
            -1.3862934, 3.3858999,   -1.5133942, 0.29624507, 0, 0.0081663491, 0, -2.4395808,
            -0.1269177, -0.34163332, 0.41758345, 0.85952322, -2.6192317
        ],
        [
            -1.2039719,  3.0402486,   -1.3567681, 0.30791729, 0, 0.027347709, 0, -2.4984837,
            -0.13590818, -0.32561803, 0.41253634, 0.82889699, -2.5774919
        ],
        [
            -1.0498212,  2.7532011,   -1.2281274, 0.31792595, 0, 0.043081617, 0, -2.560905,
            -0.14308104, -0.31423828, 0.40689449, 0.79979664, -2.551027
        ],
        [
            -0.91628981, 2.5092491,   -1.1198329, 0.3266899,  0, 0.05608876, 0, -2.6256488,
            -0.14878018, -0.30682761, 0.40077333, 0.77239357, -2.5350937
        ],
        [
            -0.79850678, 2.2982771,   -1.0269435, 0.33447807, 0, 0.066920019, 0, -2.6917382,
            -0.1532919,  -0.30202398, 0.39428215, 0.74673447, -2.5266665
        ],
        [
            -0.69314626, 2.1133007,   -0.94607929, 0.34147278, 0, 0.075998976, 0, -2.7583895,
            -0.15685291, -0.29968694, 0.38752105,  0.72279029, -2.5237106
        ],
        [
            -0.59783608, 1.9492851,   -0.87482711, 0.34780368, 0, 0.083655173, 0, -2.8249866,
            -0.15965646, -0.29909275, 0.38057929,  0.70048791, -2.5248045
        ],
        [
            -0.51082471, 1.8024805,   -0.81140719, 0.35356718, 0, 0.0901482, 0, -2.8910563,
            -0.1618581,  -0.29990905, 0.37353472,  0.67973055, -2.5289265
        ],
        [
            -0.35667403, 1.5496644,   -0.70297445, 0.36367556, 0, 0.10042872, 0, -3.0202915,
            -0.16492668, -0.30363606, 0.35939296,  0.64241889, -2.5434427
        ],
        [
            -0.22314263, 1.3384951,   -0.61316946, 0.372241,   0, 0.10804952, 0, -3.1443044,
            -0.16677241, -0.30996359, 0.34550705,  0.60999481, -2.5632309
        ],
        [
            -0.1053596,  1.1584051,  -0.53714393, 0.37957458, 0, 0.11380282, 0, -3.2622972,
            -0.16784045, -0.3174928, 0.33214991,  0.58167251, -2.5859484
        ],
        [
            9.1715208e-07, 1.0022195,  -0.47164267, 0.38590861,  0,         0.11821826,
            0,             -3.3740616, -0.16840896, -0.32561005, 0.3194883, 0.55678082,
            -2.6101867
        ],
        [
            0.095311097, 0.86487736,  -0.41438934, 0.39142247, 0, 0.12165692, 0, -3.4797185,
            -0.16865375, -0.33394336, 0.30761123,  0.5347631,  -2.635081
        ],
        [
            0.18232247, 0.74269797,  -0.36373968, 0.39625721, 0, 0.12436957, 0, -3.5795599,
            -0.1686873, -0.34229721, 0.29655276,  0.51516344, -2.660094
        ],
        [
            0.266172,    0.62777517,  -0.31635126, 0.40072359, 0, 0.12662906, 0, -3.6785009,
            -0.16857515, -0.35071445, 0.28582213,  0.49678701, -2.6861113
        ],
        [
            0.40546603,  0.4427345,   -0.24059055, 0.40770316, 0, 0.12969243, 0, -3.84798,
            -0.16813806, -0.36678069, 0.26812317,  0.46747932, -2.7330971
        ],
        [
            0.53062917,  0.28239911,  -0.17551344, 0.41349244, 0, 0.13180285, 0, -4.0048127,
            -0.16754435, -0.38210035, 0.25263761,  0.44253329, -2.7788735
        ],
        [
            0.6931481, 0.08198971, -0.094958234, 0.42032797,  0,          0.13378631,
            0,         -4.2134762, -0.16658633,  -0.40447218, 0.23347665, 0.41218079,
            -2.8425943
        ],
        [
            0.99325269,  -0.26750988, 0.04329351, 0.43102341, 0, 0.13573077, 0, -4.6090058,
            -0.16455507, -0.44996124, 0.20167488, 0.36215967, -2.9702821
        ],
        [
            1.0986132,   -0.38461173, 0.08895252, 0.43423235, 0, 0.13601657, 0, -4.7498938,
            -0.16382403, -0.4674424,  0.19169236, 0.3463766,  -3.0175049
        ],
        [
            1.2089613,   -0.50448259, 0.13533798, 0.43731749, 0, 0.13614921, 0, -4.898139,
            -0.16307044, -0.48662399, 0.18188691, 0.33077694, -3.0680039
        ],
        [
            1.3862953,   -0.69169197, 0.20706147, 0.44173368, 0, 0.13607611, 0, -5.1373053,
            -0.16190743, -0.52019887, 0.16745763, 0.30758414, -3.1509974
        ],
        [
            1.6094388,   -0.91882307, 0.29290683, 0.44645116, 0, 0.13561393, 0, -5.4389282,
            -0.16056216, -0.56663493, 0.15143682, 0.28140564, -3.2579159
        ],
        [
            1.7917604,   -1.0983195,  0.3598558,  0.44970808, 0, 0.13502767, 0, -5.6852722,
            -0.15957718, -0.60645836, 0.13991467, 0.26223811, -3.3467787
        ],
        [
            1.9459111,   -1.2463439, 0.41449166, 0.45210148, 0, 0.13443557, 0, -5.8931764,
            -0.15882828, -0.6438664, 0.13113956, 0.24742143, -3.4226787
        ],
        [
            2.0794425,   -1.3720887,  0.46051012, 0.4539406,  0, 0.13387523, 0, -6.0728623,
            -0.15824076, -0.67904305, 0.12417653, 0.23551758, -3.4888557
        ],
        [
            2.302586,    -1.5776562,  0.53499703, 0.45659317, 0, 0.13288106, 0, -6.3720961,
            -0.15737998, -0.74175021, 0.11370485, 0.21735305, -3.6000936
        ],
        [
            2.7080511,   -1.9388708,  0.66382246,  0.46032369, 0, 0.13102796, 0, -6.912032,
            -0.15616276, -0.87602944, 0.097761614, 0.18904769, -3.8034858
        ],
        [
            2.9957332,   -2.1873039,  0.75106832,  0.46231202, 0, 0.129758, 0, -7.2920377,
            -0.15552455, -0.99173061, 0.088369653, 0.17198623, -3.9482811
        ],
        [
            3.4011983,   -2.5287442, 0.86944146,  0.46443409, 0, 0.12810639, 0, -7.8234671,
            -0.15487379, -1.1797422, 0.077205205, 0.15132255, -4.1525566
        ],
        [
            3.9122261,   -2.9478071, 1.0127317,   0.46629272, 0, 0.12630613, 0, -8.4870692,
            -0.15436749, -1.4712181, 0.065765504, 0.12883086, -4.4159788
        ],
        [
            4.6052753,   -3.5016521, 1.1995982, 0.46786989, 0, 0.12436074, 0, -9.3780694,
            -0.15398038, -2.0111532, 0.0536072, 0.10590599, -4.7619993
        ],
        [
            5.2983873,   -4.0441126, 1.3807557,   0.46879055,  0, 0.1229163, 0, -10.2614,
            -0.15355919, -2.7442401, 0.044158072, 0.087675632, -5.1084784
        ],
        [
            6.214787,    -4.7501126, 1.6149346,   0.46945042,  0, 0.12159211, 0, -11.421197,
            -0.15364255, -4.2707104, 0.034528057, 0.068804388, -5.5667403
        ],
        [
            6.9077759,   -5.2785771, 1.7896287,   0.46971791, 0, 0.12091589, 0, -12.294186,
            -0.14919057, -5.7082859, 0.028810108, 0.05749535, -5.9132958
        ],
        [
            7.6010048,   -5.8042574, 1.9632099,   0.46987559,  0, 0.12043919, 0, -13.165209,
            -0.15114382, -8.1485299, 0.024102558, 0.048143019, -6.2599595
        ],
        [
            8.5174171,   -6.4962666, 2.1916964,   0.46998886,  0, 0.12000331, 0, -14.314403,
            -0.14445068, -12.173616, 0.019088606, 0.038152484, -6.7182108
        ],
        [
            9.2104756,   -7.0182192, 2.3641363,   0.4700338,  0, 0.11978373, 0, -15.182423,
            -0.27602959, -17.482818, 0.016021785, 0.03203123, -7.0647608
        ],
        [
            10.819963,   -8.2278723, 2.7643243,   0.47007301,  0, 0.1194896, 0, -17.196286,
            -0.13003945, -34.417402, 0.010690848, 0.021379237, -7.8695389
        ],
        [
            11.513023, -8.7481808, 2.9367081,   0.47007451, 0,            0.11941975,
            0,         -18.063019, -0.11123476, -41.56601,  0.0089861421, 0.017971059,
            -8.2160842
        ],
        [
            13.122508, -9.9558895, 3.3373537,   0.47005419, 0,            0.11932638,
            0,         -20.075351, -0.15372593, -128.77258, 0.0060062576, 0.012012283,
            -9.02088
        ],
        [
            13.815533, -10.475778, 3.5100123,   0.47003402, 0,            0.11930422,
            0,         -20.94173,  -0.15372699, -181.98664, 0.0050502702, 0.010100436,
            -9.3674308
        ],
        [
            16.118199, -12.202945, 4.0842275,   0.46987591, 0,            0.11926748,
            0,         -23.820173, -0.15372878, -574.88581, 0.0028395609, 0.0056791689,
            -10.519044
        ],
        [
            18.420859, -13.929979, 4.6590407,   0.46937055, 0,            0.11925588,
            0,         -26.698545, -0.15372935, -1817.3259, 0.0015966174, 0.0031934211,
            -11.671255
        ]
    ],

    tail_shock_table => [
        [ 4.5710641, -2.0629,    -0.021551008,  -0.021551008 ],
        [ 5.2782988, -2.6465974, -0.027794387,  -0.011421234 ],
        [ 6.2046828, -3.579949,  -0.024690769,  -0.0064415454 ],
        [ 6.9017708, -4.4410278, -0.021825744,  -0.0044458216 ],
        [ 7.5974377, -5.4682357, -0.019044334,  -0.0031615895 ],
        [ 8.5156258, -7.1394391, -0.015721405,  -0.0020737595 ],
        [ 9.2094117, -8.6928876, -0.013521695,  -0.0015291613 ],
        [ 10.819646, -13.573985, -0.0094139099, -0.00077562337 ],
        [ 11.512835, -16.436718, -0.0080511599, -0.00057315803 ],
        [ 13.122451, -25.236874, -0.0055159619, -0.00030411282 ],
        [ 13.815499, -30.299574, -0.0046796335, -0.00023077104 ]
    ],

    blast_info => [
        0.17076167, 0,         0.51979849, -0.26645721, 0,         0,          1.1589055,   -0.55231774,
        0.576023,   0.1608173, 0.92888643, -0.21764904, 0.1205426, 0.43014281, 0.085431983, 0,
        94.584,     -2.0629,   146.10506,  -2.3926448
    ],

};
1;
