package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='C', gamma=6.5
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'C6x5'} = {

    table_name => 'C6x5',
    symmetry   => 1,
    gamma      => 6.5,

    shock_table_info => [ 1, 6.5, 1.02e-06, 4000, 1.4e-06, 50, 2.9e-07, 2.35e-07, 5e-07, 86.198, -2.1259 ],

    shock_table => [
        [ -5.62811172, 12.00033779,  -1.99997871, -5.62974449, 0.99836591 ],
        [ -5.16950181, 11.08313152,  -1.99994672, -5.17208585, 0.99741268 ],
        [ -4.72320808, 10.19057641,  -1.99989028, -4.72724854, 0.9959516 ],
        [ -4.1602465,  9.06476705,   -1.99966226, -4.16735143, 0.99287108 ],
        [ -3.70387589, 8.15227669,   -1.99916148, -3.71511147, 0.98870625 ],
        [ -3.3200588,  7.3851258,    -1.99819891, -3.3365905,  0.98334737 ],
        [ -2.99491154, 6.7356559,    -1.99655792, -3.01785818, 0.97683192 ],
        [ -2.73492809, 6.21686213,   -1.99422333, -2.76477091, 0.96980345 ],
        [ -2.49196337, 5.73273533,   -1.99065624, -2.53013373, 0.96129152 ],
        [ -2.27302631, 5.29741879,   -1.98562619, -2.3206965,  0.95156148 ],
        [ -2.07171268, 4.89833673,   -1.97871316, -2.13021312, 0.94046176 ],
        [ -1.88473061, 4.52916694,   -1.96947175, -1.95550077, 0.92790368 ],
        [ -1.70827329, 4.18265476,   -1.95731948, -1.79298265, 0.91369397 ],
        [ -1.53845197, 3.85154143,   -1.94145725, -1.63915476, 0.89751285 ],
        [ -1.371515,   3.52908811,   -1.92078221, -1.49084086, 0.87890167 ],
        [ -1.20218142, 3.20604792,   -1.89349005, -1.34382178, 0.85701957 ],
        [ -1.01676687, 2.85838648,   -1.85501554, -1.18741891, 0.82939653 ],
        [ -0.84127403, 2.53670845,   -1.80946182, -1.04441982, 0.79973345 ],
        [ -0.6803832,  2.24948084,   -1.75977039, -0.91812782, 0.7697644 ],
        [ -0.52945396, 1.9878049,    -1.70681564, -0.80420041, 0.73962952 ],
        [ -0.38609633, 1.74701348,   -1.6518394,  -0.70030388, 0.70967083 ],
        [ -0.24784472, 1.52249417,   -1.59575302, -0.60423281, 0.68005457 ],
        [ -0.11005865, 1.30657393,   -1.53823012, -0.51257808, 0.65034834 ],
        [ 0.02799127,  1.09822159,   -1.48034494, -0.42483679, 0.62089611 ],
        [ 0.16727611,  0.89604717,   -1.42295522, -0.34038249, 0.59195194 ],
        [ 0.3090242,   0.69836546,   -1.36667554, -0.25849421, 0.56368056 ],
        [ 0.45480303,  0.50316936,   -1.31187544, -0.1783481,  0.53615989 ],
        [ 0.60586888,  0.30903838,   -1.25896556, -0.09938975, 0.50952161 ],
        [ 0.76346133,  0.11468957,   -1.20828812, -0.0211418,  0.48389215 ],
        [ 0.92916999,  -0.08146504,  -1.1600246,  0.05697549,  0.45934095 ],
        [ 1.10443442,  -0.28069739,  -1.11437881, 0.13539544,  0.43596795 ],
        [ 1.29091889,  -0.48442462,  -1.07147457, 0.2145918,   0.41384709 ],
        [ 1.49042204,  -0.69409552,  -1.03140135, 0.29503392,  0.39304503 ],
        [ 1.70528016,  -0.9115973,   -0.99415828, 0.37734135,  0.37358827 ],
        [ 1.92746136,  -1.12872243,  -0.96119545, 0.4583729,   0.3562711 ],
        [ 2.15804113,  -1.34690293,  -0.93204283, 0.53870141,  0.34088661 ],
        [ 2.39932306,  -1.56857899,  -0.90616261, 0.61925284,  0.32718383 ],
        [ 2.65258859,  -1.79509151,  -0.88322975, 0.70053275,  0.31501732 ],
        [ 2.91988851,  -2.02838229,  -0.86290797, 0.7832524,   0.3042303 ],
        [ 3.20264411,  -2.26976028,  -0.84496836, 0.86788876,  0.29471862 ],
        [ 3.50263789,  -2.5208036,   -0.8291913,  0.95501091,  0.28637932 ],
        [ 3.8215151,   -2.78294112,  -0.81539513, 1.04513464,  0.27912632 ],
        [ 4.1612525,   -3.05785407,  -0.80340633, 1.13886375,  0.27287485 ],
        [ 4.52402392,  -3.34736314,  -0.79306743, 1.23685085,  0.26754577 ],
        [ 4.91176085,  -3.65308801,  -0.78424289, 1.33968419,  0.26306865 ],
        [ 5.32718809,  -3.97727258,  -0.77679003, 1.44816662,  0.25936681 ],
        [ 5.77261711,  -4.32183522,  -0.77058452, 1.56299614,  0.25637042 ],
        [ 6.25135401,  -4.68946717,  -0.76549405, 1.68513303,  0.25400351 ],
        [ 6.76769436,  -5.08360927,  -0.761388,   1.81578876,  0.25218952 ],
        [ 7.32652099,  -5.50813705,  -0.75814337, 1.95632034,  0.25085454 ],
        [ 7.93450582,  -5.96826968,  -0.75563841, 2.1085307,   0.24992499 ],
        [ 8.60030904,  -6.4707086,   -0.75375736, 2.27471381,  0.24933037 ],
        [ 9.33578003,  -7.02453564,  -0.75239066, 2.45795247,  0.24900395 ],
        [ 10.15815836, -7.64286019,  -0.75143605, 2.66266548,  0.24888402 ],
        [ 11.09227481, -8.34446435,  -0.75080166, 2.89515727,  0.24891516 ],
        [ 12.17735301, -9.15890413,  -0.75040595, 3.16531683,  0.24904865 ],
        [ 13.47661079, -10.13370802, -0.75017954, 3.48902088,  0.24924254 ],
        [ 15.10466231, -11.35492997, -0.75006491, 3.89498503,  0.24946135 ],
        [ 17.29692078, -12.99920174, -0.75001684, 4.44211987,  0.24967497 ],
        [ 20.66480578, -15.52513983, -0.75000228, 5.28334131,  0.24985639 ]
    ],

    energy_table => [
        [ -5.62811172, 3.316e-05,  5.384e-05,  3.316e-05,  5.384e-05,     6.0343301,   -1.1337072 ],
        [ -5.16950181, 6.963e-05,  0.00011214, 6.963e-05,  0.00011214,    5.5124491,   -1.1426104 ],
        [ -4.72320808, 0.00014239, 0.00022709, 0.00014239, 0.00022709,    5.0004907,   -1.1522756 ],
        [ -4.1602465,  0.00034707, 0.00054461, 0.00034707, 0.00054461,    4.3481522,   -1.1662668 ],
        [ -3.70387589, 0.00070603, 0.0010885,  0.00070603, 0.0010885,     3.8131258,   -1.1791387 ],
        [ -3.3200588,  0.00126896, 0.00191894, 0.00126896, 0.00191894,    3.3583611,   -1.1910921 ],
        [ -2.99491154, 0.00206472, 0.00305786, 0.00206472, 0.00305786,    2.9693607,   -1.2019288 ],
        [ -2.73492809, 0.00302288, 0.00438537, 0.00302288, 0.00438537,    2.6557261,   -1.2109289 ],
        [ -2.49196337, 0.00428324, 0.0060695,  0.00428324, 0.0060695,     2.3604766,   -1.2195168 ],
        [ -2.27302631, 0.00581903, 0.00803737, 0.00581903, 0.00803737,    2.0926265,   -1.2273295 ],
        [ -2.07171268, 0.00765549, 0.0102799,  0.00765549, 0.0102799,     1.8448231,   -1.2345996 ],
        [ -1.88473061, 0.00980355, 0.01276079, 0.00980355, 0.01276079,    1.6133382,   -1.2415141 ],
        [ -1.70827329, 0.01228769, 0.01544882, 0.01228769, 0.01544882,    1.3936797,   -1.2483135 ],
        [ -1.53845197, 0.01515124, 0.01831632, 0.01515124, 0.01831632,    1.1811193,   -1.2552014 ],
        [ -1.371515,   0.01845845, 0.02132974, 0.01845845, 0.02132974,    0.97100041,  -1.262331 ],
        [ -1.20218142, 0.02233571, 0.02446578, 0.02233571, 0.02446578,    0.7566169,   -1.2699392 ],
        [ -1.01676687, 0.02718446, 0.02779923, 0.02718446, 0.02779923,    0.52036059,  -1.278672 ],
        [ -0.84127403, 0.03231804, 0.03063407, 0.03231804, 0.03063407,    0.29522096,  -1.2873676 ],
        [ -0.6803832,  0.0374255,  0.03276751, 0.0374255,  0.03276751,    0.087436207, -1.2958832 ],
        [ -0.52945396, 0.04248875, 0.034232,   0.04248875, 0.034232,      -0.10877611, -1.3045899 ],
        [ -0.38609633, 0.0474629,  0.03506993, 0.0474629,  0.03506993,    -0.29642014, -1.3137294 ],
        [ -0.24784472, 0.05233648, 0.03534675, 0.05233648, 0.03534675,    -0.47868623, -1.323472 ],
        [ -0.11005865, 0.05719688, 0.03512367, 0.05719688, 0.03512367,    -0.66174353, -1.3340779 ],
        [ 0.02799127,  0.06200417, 0.03445315, 0.06200417, 0.03445315,    -0.84667586, -1.3454509 ],
        [ 0.16727611,  0.06673349, 0.0333994,  0.06673349, 0.0333994,     -1.0348983,  -1.3574213 ],
        [ 0.3090242,   0.07137364, 0.03202875, 0.07137364, 0.03202875,    -1.2281861,  -1.3697813 ],
        [ 0.45480303,  0.07592633, 0.03040273, 0.07592633, 0.03040273,    -1.4287978,  -1.3823227 ],
        [ 0.60586888,  0.08038287, 0.02858322, 0.08038287, 0.02858322,    -1.6385878,  -1.3947993 ],
        [ 0.76346133,  0.08473365, 0.02662966, 0.08473365, 0.02662966,    -1.8593951,  -1.4069718 ],
        [ 0.92916999,  0.08897701, 0.02459368, 0.08897701, 0.02459368,    -2.0935607,  -1.4186624 ],
        [ 1.10443442,  0.09310448, 0.02252488, 0.09310448, 0.02252488,    -2.3432299,  -1.4297141 ],
        [ 1.29091889,  0.09711051, 0.02046596, 0.09711051, 0.02046596,    -2.6108783,  -1.4399943 ],
        [ 1.49042204,  0.10098942, 0.01845375, 0.10098942, 0.01845375,    -2.8991782,  -1.4494081 ],
        [ 1.70528016,  0.10474193, 0.01651564, 0.10474193, 0.01651564,    -3.2115946,  -1.4579319 ],
        [ 1.92746136,  0.10821113, 0.0147523,  0.10821113, 0.0147523,     -3.5364095,  -1.4652263 ],
        [ 2.15804113,  0.11142468, 0.0131594,  0.11142468, 0.0131594,     -3.8750482,  -1.4714201 ],
        [ 2.39932306,  0.1144219,  0.01172094, 0.1144219,  0.01172094,    -4.2307772,  -1.4766853 ],
        [ 2.65258859,  0.11722227, 0.0104271,  0.11722227, 0.0104271,     -4.605397,   -1.4811312 ],
        [ 2.91988851,  0.11984973, 0.00926369, 0.11984973, 0.00926369,    -5.0018598,  -1.4848605 ],
        [ 3.20264411,  0.12231745, 0.00822013, 0.12231745, 0.00822013,    -5.4222051,  -1.4879749 ],
        [ 3.50263789,  0.12463915, 0.0072847,  0.12463915, 0.0072847,     -5.869025,   -1.4905578 ],
        [ 3.8215151,   0.12682466, 0.00644693, 0.12682466, 0.00644693,    -6.344712,   -1.4926662 ],
        [ 4.1612525,   0.12888373, 0.00569633, 0.12888373, 0.00569633,    -6.852156,   -1.4943949 ],
        [ 4.52402392,  0.13082453, 0.00502315, 0.13082453, 0.00502315,    -7.3945737,  -1.495781 ],
        [ 4.91176085,  0.13265169, 0.00441922, 0.13265169, 0.00441922,    -7.9747819,  -1.4968759 ],
        [ 5.32718809,  0.13437147, 0.00387622, 0.13437513, 0.0040083978,  -8.5968383,  -1.4977295 ],
        [ 5.77261711,  0.13598605, 0.00338763, 0.13608305, 0.0036742827,  -9.2641328,  -1.4983823 ],
        [ 6.25135401,  0.13749928, 0.00294699, 0.13777316, 0.0034088098,  -9.9816121,  -1.4988838 ],
        [ 6.76769436,  0.13891502, 0.00254851, 0.13947692, 0.0031425568,  -10.755653,  -1.4992377 ],
        [ 7.32652099,  0.14023532, 0.00218755, 0.141121,   0.0027707308,  -11.593556,  -1.4995006 ],
        [ 7.93450582,  0.1414627,  0.00186003, 0.1427014,  0.0024404827,  -12.505298,  -1.4996852 ],
        [ 8.60030904,  0.14259897, 0.00156267, 0.14421855, 0.0021244692,  -13.503846,  -1.4998109 ],
        [ 9.33578003,  0.14364569, 0.00129284, 0.14566563, 0.0018050795,  -14.606951,  -1.4998932 ],
        [ 10.15815836, 0.14460471, 0.0010484,  0.14701787, 0.00150281,    -15.840459,  -1.4999489 ],
        [ 11.09227481, 0.14547683, 0.00082791, 0.14831177, 0.0012810479,  -17.241609,  -1.4999868 ],
        [ 12.17735301, 0.14626289, 0.00062989, 0.14954953, 0.0010083648,  -18.869229,  -1.5000161 ],
        [ 13.47661079, 0.14696156, 0.0004549,  0.15068559, 0.00079070963, -20.818158,  -1.5000527 ],
        [ 15.10466231, 0.14756978, 0.00030263, 0.15174275, 0.00052598691, -23.260362,  -1.5001457 ],
        [ 17.29692078, 0.14808047, 0.0001749,  0.15263034, 0.00030397381, -26.549269,  -1.5002067 ],
        [ 20.66480578, 0.14847862, 7.535e-05,  0.15332232, 0.00013095945, -31.601637,  -1.5 ]
    ],

    impulse_table => [
        [
            -6.5354785, 13.815062,  -6.5361372,   0.06412088,  0,          -0.15844981,
            0,          -2.2986297, -0.011053226, -0.51018154, 0.44009743, 0.9999501,
            -5.9666283
        ],
        [
            -6.5022892, 13.748684,  -6.5029701,   0.064935728, 0,          -0.15840089,
            0,          -2.2986306, -0.011238181, -0.51013258, 0.44009738, 0.99994737,
            -5.9414458
        ],
        [
            -6.2146071,  13.173321,   -6.2155151, 0.072398115, 0, -0.15790142, 0, -2.2986406,
            -0.01297673, -0.50963258, 0.44009692, 0.99991625,  -5.7241759
        ],
        [
            -5.809142, 12.362395,  -5.8105042,   0.084210335, 0,          -0.15690319,
            0,         -2.2986671, -0.015893166, -0.50863257, 0.44009595, 0.99983852,
            -5.4212912
        ],
        [
            -5.2983164, 11.340755,  -5.3005876,   0.10146382,  0,          -0.1549099,
            0,          -2.2987492, -0.020517908, -0.50663257, 0.44009296, 0.99963144,
            -5.0461169
        ],
        [
            -4.6051692, 9.9545131,  -4.6097164,   0.12958332,  0,         -0.14994963,
            0,          -2.2991147, -0.029016026, -0.50163257, 0.4400791, 0.9988822,
            -4.5510998
        ],
        [
            -3.912022, 8.5684258,  -3.9211363,   0.16362468,  0,          -0.14015385,
            0,         -2.3004727, -0.041029915, -0.49181591, 0.44002376, 0.99667427,
            -4.0772627
        ],
        [
            -3.5065569, 7.7578409,  -3.5202579,   0.18638605,  0,          -0.13056037,
            0,          -2.3025845, -0.050239113, -0.48218278, 0.43993162, 0.99378598,
            -3.8130427
        ],
        [
            -2.9957313, 6.7372926,  -3.018659,    0.21800201,  0,          -0.11209522,
            0,          -2.3087975, -0.064799405, -0.46365259, 0.43963738, 0.98658893,
            -3.4981875
        ],
        [
            -2.6592591, 6.0659953, -2.6914772,   0.24053427,  0,          -0.094714375,
            0,          -2.317335, -0.076551648, -0.44604889, 0.43919764, 0.97805924,
            -3.3047834
        ],
        [
            -2.3025841, 5.3561215,  -2.348844,    0.26580577,  0,          -0.07083484,
            0,          -2.3338601, -0.091156526, -0.42180261, 0.43826973, 0.9636242,
            -3.1153575
        ],
        [
            -2.0402199, 4.8360422,  -2.1006258,  0.28524535,  0,          -0.049600427,
            0,          -2.3541814, -0.10336663, -0.40044361, 0.43702825, 0.94793502,
            -2.9886913
        ],
        [
            -1.897119, 4.55357,    -1.9670016,  0.29612793,  0,          -0.036843848,
            0,         -2.3695387, -0.11051943, -0.38748435, 0.43603218, 0.93704643,
            -2.9249925
        ],
        [
            -1.6094369, 3.9896159, -1.7031188,  0.31852765,  0,          -0.0093134921,
            0,          -2.413209, -0.12572427, -0.36088314, 0.43298363, 0.90909707,
            -2.8104732
        ],
        [
            -1.3862934,  3.5574894,   -1.5038427, 0.33628701, 0, 0.013009788, 0, -2.4629306,
            -0.13788636, -0.34075219, 0.42919952, 0.88098882, -2.7361198
        ],
        [
            -1.2039718,  3.2094383,  -1.3453564, 0.35096551, 0, 0.031270465, 0, -2.5171653,
            -0.14765338, -0.3262212, 0.42476635, 0.853435,   -2.6862688
        ],
        [
            -1.0498212,  2.9198285, -1.2149204, 0.36344245, 0, 0.046360175, 0, -2.5746578,
            -0.15549258, -0.316431, 0.41977623, 0.82685658, -2.6526434
        ],
        [
            -0.91628977, 2.6732292,   -1.1049071, 0.37426546, 0, 0.058943954, 0, -2.634376,
            -0.1617745,  -0.31010405, 0.41432216, 0.80148918, -2.6303384
        ],
        [
            -0.79850673, 2.4595878,   -1.010382,  0.38379777, 0, 0.06951915, 0, -2.6954745,
            -0.16680114, -0.30672485, 0.40849428, 0.7774477,  -2.6162135
        ],
        [
            -0.69314621, 2.2719677,   -0.9279681, 0.39229158, 0, 0.07846448, 0, -2.7572684,
            -0.17081893, -0.30488407, 0.40237715, 0.75476834, -2.6081453
        ],
        [
            -0.59783603, 2.1053715,   -0.8552539, 0.39992864, 0, 0.086073953, 0, -2.8192108,
            -0.17402766, -0.30508027, 0.39604808, 0.7334369,  -2.6046397
        ],
        [
            -0.51082466, 1.9560725,  -0.79045732, 0.40684436, 0, 0.092579893, 0, -2.8808722,
            -0.17658822, -0.3064712, 0.38957619,  0.71340795, -2.6046136
        ],
        [
            -0.35667398, 1.6985849,   -0.67951581, 0.4189059,  0, 0.10298912, 0, -3.0021114,
            -0.18025319, -0.31166242, 0.37643859,  0.67699365, -2.6119832
        ],
        [
            -0.22314259, 1.4832021,   -0.58749985, 0.42908012, 0, 0.1107965, 0, -3.1192319,
            -0.1825603,  -0.31880752, 0.36335402,  0.64493544, -2.6258793
        ],
        [
            -0.10535955, 1.2993503,   -0.5095244, 0.43777418, 0, 0.11674185, 0, -3.2313506,
            -0.18398315, -0.32707639, 0.35059791, 0.61663202, -2.6436944
        ],
        [
            9.6629291e-07, 1.1398206,  -0.44229876, 0.44528081,  0,          0.12133109,
            0,             -3.3381285, -0.18482376, -0.33583691, 0.33835333, 0.59153579,
            -2.6638162
        ],
        [
            0.095311146, 0.99950736,  -0.38351443, 0.45181953, 0, 0.12491642, 0, -3.4395504,
            -0.18527683, -0.34515482, 0.32673229,  0.56917236, -2.6852231
        ],
        [
            0.18232252,  0.87468271,  -0.33149882, 0.45755954, 0, 0.12774699, 0, -3.535784,
            -0.18546986, -0.35443434, 0.31579413,  0.54914044, -2.7072563
        ],
        [
            0.21994287,  0.82166495,  -0.30948792, 0.45996966, 0, 0.12884698, 0, -3.5783836,
            -0.18549578, -0.35853743, 0.31099844,  0.54061325, -2.7174786
        ],
        [
            0.3364732,   0.66099709,  -0.24309478, 0.4671481,  0, 0.13181193, 0, -3.7137798,
            -0.18538636, -0.37257061, 0.29602383,  0.51478398, -2.7516352
        ],
        [
            0.40546607,  0.56833966,  -0.20502487, 0.47118778, 0, 0.13327396, 0, -3.7961642,
            -0.18520561, -0.38157084, 0.28715666,  0.49994634, -2.7735186
        ],
        [
            0.53062922,  0.40472796, -0.13821266, 0.47810805, 0, 0.13542749, 0, -3.949252,
            -0.18470418, -0.3984221, 0.27126413,  0.47397086, -2.8160674
        ],
        [
            0.69314815, 0.20041728, -0.055555522, 0.486307,    0,          0.13738888,
            0,          -4.1538184, -0.18380767,  -0.42277145, 0.25134534, 0.44216489,
            -2.8761495
        ],
        [
            0.99325274,  -0.15524399, 0.086125985, 0.49920336, 0, 0.13908777, 0, -4.5437396,
            -0.18176646, -0.47282616, 0.21776422,  0.38932132, -2.9985588
        ],
        [
            1.0986133,   -0.27420625, 0.13285546, 0.5030896,  0, 0.13923534, 0, -4.6831576,
            -0.18100657, -0.49134139, 0.20711341, 0.37255056, -3.0443071
        ],
        [
            1.2089613,   -0.39587196, 0.18029283, 0.50683355, 0, 0.13920052, 0, -4.8300909,
            -0.18021362, -0.51285603, 0.19661013, 0.35593436, -3.0934408
        ],
        [
            1.3862953,   -0.5856604,  0.25356702, 0.51220595, 0, 0.13881958, 0, -5.0675637,
            -0.17897519, -0.54937743, 0.18109066, 0.33116306, -3.174575
        ],
        [
            1.6094389,   -0.81556265, 0.34114177, 0.51796231, 0, 0.13792246, 0, -5.3676372,
            -0.17752462, -0.59986603, 0.16378931, 0.30311841, -3.2796424
        ],
        [
            1.7917604,   -0.99698519, 0.40934143, 0.52194675, 0, 0.13695705, 0, -5.6130835,
            -0.17645262, -0.6445146,  0.15131258, 0.28253653, -3.3673162
        ],
        [
            1.9459111,   -1.1464331,  0.46493381, 0.52488026, 0, 0.13603692, 0, -5.8204261,
            -0.17563259, -0.68473755, 0.14179726, 0.2666033,  -3.4423931
        ],
        [
            2.0794425,   -1.2732796,  0.51171469, 0.52713755, 0, 0.1351913, 0, -5.9997415,
            -0.17498715, -0.72286439, 0.13424125, 0.25378996, -3.5079693
        ],
        [
            2.3025861,   -1.4804483,  0.58735231, 0.53039807, 0, 0.13372553, 0, -6.2985425,
            -0.17403885, -0.79186411, 0.12287236, 0.23421999, -3.6183911
        ],
        [
            2.7080512,   -1.8439518,  0.71793759, 0.5349931,  0, 0.13105769, 0, -6.8381004,
            -0.17269247, -0.93935115, 0.10556148, 0.20369401, -3.8207493
        ],
        [
            2.9957332,   -2.093633,  0.80622252,  0.53744661, 0, 0.12925662, 0, -7.2180341,
            -0.17198678, -1.0574347, 0.095369085, 0.18528355, -3.9650541
        ],
        [
            3.4011983,   -2.4364403, 0.92582852,  0.54006818, 0, 0.12693401, 0, -7.7495096,
            -0.17127078, -1.2651517, 0.083263507, 0.16298302, -4.1688587
        ],
        [
            3.9122263,   -2.8567489, 1.070372,    0.54236823, 0, 0.12442111, 0, -8.4132849,
            -0.17071549, -1.585491,  0.070875386, 0.13878789, -4.431459
        ],
        [
            4.6052076,   -3.4116648, 1.2585291,   0.54432157, 0, 0.1217198, 0, -9.3044838,
            -0.17029919, -2.1694407, 0.057731301, 0.11402645, -4.777374
        ],
        [
            5.298432,    -3.9549286, 1.440705,    0.54546206,  0, 0.11971921, 0, -10.188219,
            -0.17010848, -2.9935166, 0.047529347, 0.094355663, -5.123886
        ],
        [
            6.2146549,  -4.661368,  1.6758085,   0.5462806,   0, 0.11788977, 0, -11.348062,
            -0.1694085, -4.5641344, 0.037149345, 0.074022323, -5.5820572
        ],
        [
            6.9079352,   -5.1903221, 1.8511285,   0.54661285,  0, 0.11695547, 0, -12.221566,
            -0.17002346, -6.4444163, 0.030988165, 0.061839283, -5.9287579
        ],
        [
            7.6010176,   -5.7160701, 2.0251117,  0.54680876,  0, 0.11629629, 0, -13.092515,
            -0.17002894, -9.0243119, 0.02592178, 0.051775394, -6.2753481
        ],
        [
            8.5174688,   -6.4082592, 2.254057,    0.5469499,   0, 0.11569762, 0, -14.241861,
            -0.24248195, -13.562515, 0.020526471, 0.041025798, -6.7336172
        ],
        [
            9.2103431,   -6.9301463, 2.426716,    0.54700629,  0, 0.11539487, 0, -15.109703,
            -0.66784588, -16.843772, 0.017228311, 0.034443075, -7.0800753
        ],
        [
            10.819847,   -8.1399059, 2.8273491,   0.54705754,  0, 0.11498931, 0, -17.123661,
            -0.17006869, -44.175276, 0.011494888, 0.022987079, -7.8848565
        ],
        [
            11.51311, -8.6603879, 2.9999187,   0.54706132, 0,            0.11489298,
            0,        -17.990665, -0.16422074, -60.001216, 0.0096612839, 0.019321227,
            -8.231498
        ],
        [
            13.122397,   -9.867977,  3.4007452,    0.54704414, 0, 0.11476432, 0, -20.002777,
            -0.17007949, -139.09489, 0.0064576532, 0.01291506, -9.0361716
        ],
        [
            13.815623, -10.388022, 3.5735255,   0.54702512, 0,            0.11473379,
            0,         -20.86941,  -0.17008104, -196.5937,  0.0054295227, 0.010858943,
            -9.3828046
        ],
        [
            16.11829, -12.1152,   4.1479036,   0.54687331, 0,            0.11468331,
            0,        -23.747888, -0.17008365, -621.02493, 0.0030526774, 0.0061054456,
            -10.534284
        ],
        [
            18.420749, -13.842087, 4.722756,    0.54638663, 0,            0.11466746,
            0,         -26.626015, -0.17008449, -1962.9782, 0.0017164719, 0.0034332722,
            -11.68597
        ]
    ],

    tail_shock_table => [
        [ 4.4810107, -2.1259,    -0.022117878,  -0.022117878 ],
        [ 4.569372,  -2.1960035, -0.025004447,  -0.018478719 ],
        [ 5.2770899, -2.8167182, -0.025839381,  -0.010283635 ],
        [ 6.2039117, -3.8096746, -0.022799338,  -0.0059058326 ],
        [ 6.901549,  -4.7265574, -0.02012594,   -0.0041027289 ],
        [ 7.5972223, -5.8198135, -0.017553624,  -0.0029293923 ],
        [ 8.5155623, -7.5989357, -0.014488942,  -0.001927196 ],
        [ 9.2092103, -9.2523331, -0.01246305,   -0.0014230113 ],
        [ 10.819509, -14.524064, -0.0087240822, -0.00070204715 ],
        [ 11.512909, -17.496067, -0.0074227134, -0.00053399853 ],
        [ 13.122337, -26.862368, -0.0050866081, -0.00028278604 ],
        [ 13.815587, -32.241966, -0.004313986,  -0.0002143571 ]
    ],

    blast_info => [
        0.15990032, 0,        0.51159738, -0.29016855, 0,          0,          1.1021165,   -0.62804833,
        0.51201194, 0.168254, 0.92976027, -0.23788105, 0.10513098, 0.44009828, 0.084134355, 0,
        86.198,     -2.1259,  128.64729,  -2.4344598
    ],

};
1;
