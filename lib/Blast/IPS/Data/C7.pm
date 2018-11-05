package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='C', gamma=7
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'C7'} = {

    table_name  => 'C7',
    symmetry    => 1,
    gamma       => 7,
    data_source => 'C4000_G7/moc2_from_moc_r5000/',

    shock_table_info => [ 1, 7, 1.1e-06, 4000, 1.47e-06, 70, 2.86e-07, 3.19e-07, 5e-07, 82.466, -2.1528 ],

    shock_table => [
        [ -5.58930901, 12.00033786,  -1.9999785,  -5.59094962, 0.99835806 ],
        [ -5.12913159, 11.07999675,  -1.99994604, -5.13173212, 0.99739615 ],
        [ -4.64342607, 10.10862299,  -1.99987784, -4.64765613, 0.99576126 ],
        [ -4.06031073, 8.9425245,    -1.99961503, -4.06790154, 0.99238192 ],
        [ -3.65969885, 8.14153629,   -1.99914244, -3.67104975, 0.98858978 ],
        [ -3.25281691, 7.32830863,   -1.99807731, -3.2699113,  0.98277686 ],
        [ -2.94692768, 6.71734109,   -1.99646118, -2.9702002,  0.97650029 ],
        [ -2.68336134, 6.19143488,   -1.99402288, -2.71373887, 0.96925758 ],
        [ -2.44554702, 5.71762175,   -1.99043182, -2.4842004,  0.96079712 ],
        [ -2.22734447, 5.28382465,   -1.98530449, -2.2755831,  0.95097879 ],
        [ -2.02715049, 4.88703656,   -1.97828822, -2.08628224, 0.93981471 ],
        [ -1.84084021, 4.51928458,   -1.96890977, -1.91232594, 0.92717203 ],
        [ -1.66460682, 4.17332743,   -1.95656487, -1.75015335, 0.91284266 ],
        [ -1.49523336, 3.84323236,   -1.94049302, -1.59688404, 0.89655802 ],
        [ -1.32829639, 3.52096452,   -1.9195136,  -1.44874242, 0.87779008 ],
        [ -1.15873888, 3.19774627,   -1.8918061,  -1.30173433, 0.85570528 ],
        [ -0.97220467, 2.84836691,   -1.85257948, -1.14466667, 0.82770436 ],
        [ -0.7978315,  2.52918782,   -1.8067913,  -1.00287586, 0.79804094 ],
        [ -0.63783638, 2.24400034,   -1.75694477, -0.8775561,  0.76810133 ],
        [ -0.4875369,  1.98385113,   -1.70387007, -0.76435207, 0.73799972 ],
        [ -0.34480668, 1.7445359,    -1.64888203, -0.66113775, 0.70811932 ],
        [ -0.20657275, 1.5204666,    -1.59262525, -0.56529487, 0.67848713 ],
        [ -0.06885542, 1.30509133,   -1.53501901, -0.47390037, 0.64880555 ],
        [ 0.06897057,  1.09751709,   -1.47717916, -0.38650847, 0.61943822 ],
        [ 0.20825542,  0.89578366,   -1.41979764, -0.30225322, 0.59055523 ],
        [ 0.35022743,  0.69824062,   -1.36348914, -0.22043138, 0.56232257 ],
        [ 0.49623019,  0.50320875,   -1.30871251, -0.14035617, 0.53486321 ],
        [ 0.64751997,  0.3092637,    -1.25587259, -0.061471,   0.50830585 ],
        [ 0.80533636,  0.11511856,   -1.20530473, 0.01670381,  0.48277338 ],
        [ 0.97126895,  -0.08081763,  -1.15718343, 0.09475007,  0.45833168 ],
        [ 1.14698124,  -0.28006948,  -1.11165087, 0.17320092,  0.43504879 ],
        [ 1.33391357,  -0.48378945,  -1.06889221, 0.25242516,  0.41302949 ],
        [ 1.53408851,  -0.69366213,  -1.02894822, 0.33298323,  0.3923149 ],
        [ 1.74961842,  -0.91133239,  -0.99186518, 0.41540122,  0.37295556 ],
        [ 1.97188604,  -1.12806063,  -0.95915276, 0.49633943,  0.35577694 ],
        [ 2.20282725,  -1.3461348,   -0.9301999,  0.5766923,   0.3404988 ],
        [ 2.44433311,  -1.56759681,  -0.90452509, 0.65723777,  0.32690021 ],
        [ 2.69804651,  -1.7941168,   -0.8817657,  0.7386002,   0.31481625 ],
        [ 2.96579428,  -2.02742916,  -0.8616103,  0.82141488,  0.3041037 ],
        [ 3.24899774,  -2.2688448,   -0.84382868, 0.90615911,  0.29465813 ],
        [ 3.54966332,  -2.52012881,  -0.82818984, 0.99346627,  0.28637101 ],
        [ 3.86943625,  -2.78270252,  -0.81451711, 1.08384734,  0.27915935 ],
        [ 4.2100746,   -3.05806615,  -0.80264749, 1.17784271,  0.27294437 ],
        [ 4.57353901,  -3.3478745,   -0.79242695, 1.27604872,  0.26764914 ],
        [ 4.96245166,  -3.65429703,  -0.78370147, 1.37923858,  0.26319404 ],
        [ 5.3788564,   -3.97903951,  -0.77634459, 1.48803269,  0.25951115 ],
        [ 5.82586889,  -4.32464617,  -0.77021857, 1.60333729,  0.2565244 ],
        [ 6.3065924,   -4.69364492,  -0.76519755, 1.72605617,  0.25416157 ],
        [ 6.82532171,  -5.08947306,  -0.76115239, 1.85739845,  0.25234728 ],
        [ 7.38673899,  -5.51585232,  -0.75796175, 1.99866922,  0.25100871 ],
        [ 7.99791648,  -5.97830459,  -0.75550208, 2.15177097,  0.25007186 ],
        [ 8.66811425,  -6.48398106,  -0.75365726, 2.31914588,  0.24946663 ],
        [ 9.40878097,  -7.04165806,  -0.75232039, 2.50377529,  0.24912758 ],
        [ 10.23755617, -7.6647443,   -0.75138926, 2.71017716,  0.24899331 ],
        [ 11.18067057, -8.37307123,  -0.75077208, 2.94500407,  0.24900889 ],
        [ 12.27794751, -9.19664201,  -0.75038871, 3.21829464,  0.24912623 ],
        [ 13.59480476, -10.18463363, -0.75017055, 3.54647468,  0.24930389 ],
        [ 15.25085643, -11.42684857, -0.75006091, 3.95950883,  0.24950692 ],
        [ 17.49271588, -13.10831729, -0.75001547, 4.51910739,  0.2497056 ],
        [ 20.97220359, -15.71795579, -0.750002,   5.3882857,   0.24987339 ]
    ],

    energy_table => [
        [ -5.58930901, 2.823e-05,  4.63e-05,   2.823e-05,  4.63e-05,      6.0084765,   -1.1296778 ],
        [ -5.12913159, 5.986e-05,  9.738e-05,  5.986e-05,  9.738e-05,     5.4867615,   -1.1381293 ],
        [ -4.64342607, 0.00013136, 0.00021135, 0.00013136, 0.00021135,    4.9317084,   -1.1480319 ],
        [ -4.06031073, 0.00033307, 0.00052656, 0.00033307, 0.00052656,    4.2585947,   -1.1615926 ],
        [ -3.65969885, 0.00062444, 0.00097163, 0.00062444, 0.00097163,    3.7912509,   -1.1721147 ],
        [ -3.25281691, 0.00116874, 0.00178107, 0.00116874, 0.00178107,    3.3120494,   -1.1838433 ],
        [ -2.94692768, 0.0018545,  0.00277031, 0.0018545,  0.00277031,    2.9485212,   -1.1932464 ],
        [ -2.68336134, 0.00273827, 0.00400487, 0.00273827, 0.00400487,    2.632928,    -1.2016874 ],
        [ -2.44554702, 0.00386227, 0.00551919, 0.00386227, 0.00551919,    2.3462286,   -1.2095824 ],
        [ -2.22734447, 0.00525516, 0.00731923, 0.00525516, 0.00731923,    2.0814893,   -1.2170989 ],
        [ -2.02715049, 0.00691884, 0.00936843, 0.00691884, 0.00936843,    1.8371307,   -1.224312 ],
        [ -1.84084021, 0.00887053, 0.01164255, 0.00887053, 0.01164255,    1.6083868,   -1.2314867 ],
        [ -1.66460682, 0.01113583, 0.01411595, 0.01113583, 0.01411595,    1.3907361,   -1.2388461 ],
        [ -1.49523336, 0.01374692, 0.01675488, 0.01374692, 0.01675488,    1.1802846,   -1.2465284 ],
        [ -1.32829639, 0.01677446, 0.01953983, 0.01677446, 0.01953983,    0.97153499,  -1.2547022 ],
        [ -1.15873888, 0.02033388, 0.02244687, 0.02033388, 0.02244687,    0.75806193,  -1.2635676 ],
        [ -0.97220467, 0.02481416, 0.02555656, 0.02481416, 0.02555656,    0.52142623,  -1.2738659 ],
        [ -0.7978315,  0.02950441, 0.02817405, 0.02950441, 0.02817405,    0.29843843,  -1.2839694 ],
        [ -0.63783638, 0.03417698, 0.03015424, 0.03417698, 0.03015424,    0.092250182, -1.2937497 ],
        [ -0.4875369,  0.0388188,  0.03152627, 0.0388188,  0.03152627,    -0.10291047, -1.3035569 ],
        [ -0.34480668, 0.04338169, 0.03232527, 0.04338169, 0.03232527,    -0.28965577, -1.3135869 ],
        [ -0.20657275, 0.04787546, 0.0326112,  0.04787546, 0.0326112,     -0.47193521, -1.3240542 ],
        [ -0.06885542, 0.05235964, 0.03243618, 0.05235964, 0.03243618,    -0.65502459, -1.3351914 ],
        [ 0.06897057,  0.05679414, 0.03184887, 0.05679414, 0.03184887,    -0.83983954, -1.3468995 ],
        [ 0.20825542,  0.06116819, 0.03090581, 0.06116819, 0.03090581,    -1.0282826,  -1.359062 ],
        [ 0.35022743,  0.0654708,  0.02966613, 0.0654708,  0.02966613,    -1.2221183,  -1.3714989 ],
        [ 0.49623019,  0.06969621, 0.02818773, 0.06969621, 0.02818773,    -1.4232902,  -1.3840129 ],
        [ 0.64751997,  0.07383627, 0.02652756, 0.07383627, 0.02652756,    -1.6336412,  -1.3963904 ],
        [ 0.80533636,  0.07788199, 0.02474031, 0.07788199, 0.02474031,    -1.8550032,  -1.4084057 ],
        [ 0.97126895,  0.08183173, 0.02287355, 0.08183173, 0.02287355,    -2.0897069,  -1.4199164 ],
        [ 1.14698124,  0.08568224, 0.02097061, 0.08568224, 0.02097061,    -2.3402204,  -1.430799 ],
        [ 1.33391357,  0.08942272, 0.01907364, 0.08942272, 0.01907364,    -2.6086964,  -1.4409008 ],
        [ 1.53408851,  0.09305163, 0.01721483, 0.09305163, 0.01721483,    -2.8981319,  -1.4501596 ],
        [ 1.74961842,  0.09656486, 0.01542212, 0.09656486, 0.01542212,    -3.2116702,  -1.45854 ],
        [ 1.97188604,  0.09980762, 0.01379277, 0.09980762, 0.01379277,    -3.5367293,  -1.4656974 ],
        [ 2.20282725,  0.10281844, 0.0123166,  0.10281844, 0.0123166,     -3.8759944,  -1.4717875 ],
        [ 2.44433311,  0.10562781, 0.01098219, 0.10562781, 0.01098219,    -4.2321301,  -1.4769719 ],
        [ 2.69804651,  0.10825753, 0.009779,   0.10825753, 0.009779,      -4.6074775,  -1.4813537 ],
        [ 2.96579428,  0.11072688, 0.00869556, 0.11072688, 0.00869556,    -5.0046544,  -1.4850259 ],
        [ 3.24899774,  0.11304788, 0.00772237, 0.11304788, 0.00772237,    -5.4257061,  -1.4881026 ],
        [ 3.54966332,  0.11523461, 0.00684823, 0.11523461, 0.00684823,    -5.8735603,  -1.4906614 ],
        [ 3.86943625,  0.1172955,  0.00606393, 0.1172955,  0.00606393,    -6.3506139,  -1.492759 ],
        [ 4.2100746,   0.11923785, 0.00536052, 0.11923785, 0.00536052,    -6.8594343,  -1.4944577 ],
        [ 4.57353901,  0.12106818, 0.00472939, 0.12106818, 0.00472939,    -7.4029002,  -1.495816 ],
        [ 4.96245166,  0.12279395, 0.00416194, 0.12279395, 0.00416194,    -7.9848823,  -1.4969048 ],
        [ 5.3788564,   0.12441765, 0.00365158, 0.12442945, 0.0038174407,  -8.6084105,  -1.4977526 ],
        [ 5.82586889,  0.12594411, 0.00319144, 0.12606949, 0.0035183825,  -9.2780895,  -1.498404 ],
        [ 6.3065924,   0.12737554, 0.00277598, 0.127681,   0.0032071273,  -9.9985537,  -1.4988997 ],
        [ 6.82532171,  0.12871511, 0.00239993, 0.12928607, 0.0029713387,  -10.776184,  -1.4992507 ],
        [ 7.38673899,  0.12996398, 0.00205923, 0.13086548, 0.0026410054,  -11.617977,  -1.499509 ],
        [ 7.99791648,  0.13112511, 0.0017499,  0.13236969, 0.0023118954,  -12.534511,  -1.4996917 ],
        [ 8.66811425,  0.13220067, 0.00146878, 0.1338017,  0.0019546134,  -13.539654,  -1.4998161 ],
        [ 9.40878097,  0.1331909,  0.00121375, 0.13514226, 0.0016829957,  -14.650555,  -1.4998977 ],
        [ 10.23755617, 0.1340976,  0.00098276, 0.13643657, 0.001462943,   -15.893662,  -1.4999532 ],
        [ 11.18067057, 0.13492217, 0.00077429, 0.13774728, 0.001281411,   -17.308312,  -1.4999917 ],
        [ 12.27794751, 0.13566425, 0.00058346, 0.13897557, 0.00097190186, -18.954236,  -1.5000236 ],
        [ 13.59480476, 0.13632322, 0.00042242, 0.14006537, 0.00069841967, -20.929578,  -1.5000673 ],
        [ 15.25085643, 0.13689581, 0.00027907, 0.14101213, 0.00046136649, -23.413819,  -1.5001845 ],
        [ 17.49271588, 0.13737465, 0.0001593,  0.14180382, 0.00026334685, -26.777282,  -1.5002437 ],
        [ 20.97220359, 0.1377448,  6.674e-05,  0.14241581, 0.00011033525, -31.997054,  -1.5 ]
    ],

    impulse_table => [
        [
            -6.4966758,  13.815062,   -6.4973376, 0.066015386, 0, -0.15357466, 0, -2.3159792,
            -0.01169993, -0.50624112, 0.44402881, 0.99995485,  -6.042821
        ],
        [
            -6.2146071, 13.250926,  -6.2154847,   0.073525115, 0,          -0.15308374,
            0,          -2.3159891, -0.013472046, -0.50574957, 0.44402835, 0.99992855,
            -5.8274978
        ],
        [
            -5.809142, 12.44,      -5.8104586,   0.085667302, 0,          -0.15208583,
            0,         -2.3160156, -0.016499798, -0.50474957, 0.44402744, 0.99986141,
            -5.5212781
        ],
        [
            -5.2983164, 11.418359,  -5.3005115,   0.10346942,  0,          -0.15009367,
            0,          -2.3160973, -0.021301051, -0.50274956, 0.44402464, 0.99968103,
            -5.1418594
        ],
        [
            -4.6051692, 10.032114,  -4.6095645,   0.13265852,  0,          -0.14513916,
            0,          -2.3164595, -0.030123445, -0.49774956, 0.44401167, 0.99902126,
            -4.6409451
        ],
        [
            -3.912022, 8.646013,   -3.9208307,   0.16830889,  0,          -0.13536815,
            0,         -2.3177977, -0.042595345, -0.48793246, 0.44395991, 0.99705342,
            -4.1609231
        ],
        [
            -3.5065569, 7.8354051,  -3.5197967,   0.19235834,  0,          -0.12581714,
            0,          -2.3198703, -0.052154876, -0.47848154, 0.44387373, 0.99445578,
            -3.8928626
        ],
        [
            -2.9957313, 6.8147836,  -3.0178848,   0.22607597,  0,          -0.10748668,
            0,          -2.3259431, -0.067266313, -0.45994829, 0.44359846, 0.98792626,
            -3.5727989
        ],
        [
            -2.6592591, 6.1433778,  -2.6903867,   0.2503404,   0,          -0.090297797,
            0,          -2.3342587, -0.079459301, -0.44252403, 0.44318701, 0.98012526,
            -3.3756412
        ],
        [
            -2.3025841, 5.4332755, -2.3472728,   0.27778169,  0,          -0.066782729,
            0,          -2.350306, -0.094604264, -0.41882683, 0.44231835, 0.96682351,
            -3.1818443
        ],
        [
            -2.0402199, 4.9128951,  -2.0985702,  0.2990326,   0,          -0.045959797,
            0,          -2.3699908, -0.10725789, -0.39784022, 0.44115529, 0.95226107,
            -3.0516311
        ],
        [
            -1.897119, 4.6301831,  -1.9646216,  0.31097104,  0,          -0.033481708,
            0,         -2.3848449, -0.11466713, -0.38526406, 0.44022146, 0.94210202,
            -2.9858557
        ],
        [
            -1.6094369, 4.0655119,  -1.699929,   0.33560194,  0,          -0.0065981137,
            0,          -2.4270241, -0.13041159, -0.35927823, 0.43735958, 0.91585898,
            -2.8668155
        ],
        [
            -1.3862934,  3.6325239,   -1.4998517, 0.35514585, 0, 0.015199286, 0, -2.4749919,
            -0.14300774, -0.34024818, 0.43379941, 0.88925582, -2.7886115
        ],
        [
            -1.2039718,  3.2835028,   -1.3405815, 0.37127969, 0, 0.033065048, 0, -2.5272936,
            -0.15313588, -0.32651002, 0.42961795, 0.86299461, -2.7353901
        ],
        [
            -1.0498211, 2.9928478,   -1.2093852, 0.38496296, 0, 0.047873758, 0, -2.5827478,
            -0.1612836, -0.31738389, 0.42489788, 0.83750457, -2.6987721
        ],
        [
            -0.91628976, 2.7451588,   -1.09864,   0.39680071, 0, 0.060267501, 0, -2.6403853,
            -0.16783408, -0.31155811, 0.41972344, 0.81303903, -2.6737889
        ],
        [
            -0.79850672, 2.5304079,   -1.0034148, 0.40719839, 0, 0.070722107, 0, -2.6994114,
            -0.17309717, -0.30847858, 0.41417703, 0.78973439, -2.6572553
        ],
        [
            -0.6931462,  2.34168,     -0.92033476, 0.41643973, 0, 0.079598352, 0, -2.7591799,
            -0.17732432, -0.30750708, 0.40833674,  0.76764898, -2.6470139
        ],
        [
            -0.59783602, 2.1739922,   -0.84698834, 0.42473033, 0, 0.087176105, 0, -2.8191715,
            -0.18071881, -0.30779224, 0.40227469,  0.74678963, -2.6415424
        ],
        [
            -0.51082465, 2.0236307,   -0.78159408, 0.43222372, 0, 0.093676614, 0, -2.8789754,
            -0.18344411, -0.30926298, 0.39605604,  0.72712999, -2.6397343
        ],
        [
            -0.35667397, 1.7641317,   -0.66955618, 0.4452654,  0, 0.10412299, 0, -2.9968167,
            -0.18738431, -0.31524176, 0.38337307,  0.69121076, -2.6440126
        ],
        [
            -0.22314257, 1.546913,    -0.5765669, 0.45624631, 0, 0.11199764, 0, -3.1109759,
            -0.18990657, -0.32286626, 0.3706649,  0.65940878, -2.6553389
        ],
        [
            -0.10535954, 1.3614063,   -0.49772785, 0.46562132, 0, 0.11801721, 0, -3.2205451,
            -0.19149685, -0.33199341, 0.35820475,  0.63119943, -2.6710007
        ],
        [
            9.7658868e-07, 1.200393,   -0.42973436, 0.47371376,  0,          0.12267623,
            0,             -3.3251394, -0.19246756, -0.34080003, 0.34617978, 0.60608841,
            -2.6893038
        ],
        [
            0.095311156, 1.058752,    -0.37026513, 0.48076367, 0, 0.12632195, 0, -3.4246926,
            -0.19302183, -0.35098603, 0.3347096,   0.58363798, -2.7091621
        ],
        [
            0.18232253,  0.9327399,  -0.31763669, 0.48695475, 0, 0.12920199, 0, -3.5193237,
            -0.19329346, -0.3601921, 0.32386272,  0.56347195, -2.7298667
        ],
        [
            0.19508658,  0.91451599,  -0.31004784, 0.4878423,  0, 0.12959001, 0, -3.5334865,
            -0.19331642, -0.36178064, 0.3222485,   0.56054549, -2.7331012
        ],
        [
            0.33647321,  0.71703103,  -0.22818407, 0.49730482, 0, 0.13333458, 0, -3.6947488,
            -0.19331864, -0.37897915, 0.30413445,  0.52876931, -2.7720975
        ],
        [
            0.40546608,  0.62350936,  -0.18966307, 0.50166928, 0, 0.13481684, 0, -3.7760981,
            -0.19317473, -0.38881816, 0.29523672,  0.51373944, -2.793105
        ],
        [
            0.53062923,  0.45840387,  -0.12206457, 0.50915253, 0, 0.13698942, 0, -3.92749,
            -0.19272341, -0.40613454, 0.27921821,  0.48736966, -2.8341971
        ],
        [
            0.69314816, 0.25230659, -0.038452304, 0.51803077,  0,          0.13894138,
            0,          -4.1301897, -0.19186456,  -0.43157361, 0.25902394, 0.45498796,
            -2.8926309
        ],
        [
            0.99325275,  -0.10619089, 0.10479228, 0.5320259,  0, 0.14052888, 0, -4.5175368,
            -0.18983266, -0.48382749, 0.22473124, 0.40098889, -3.0126384
        ],
        [
            1.0986133,   -0.22601276, 0.15201057, 0.53625109, 0, 0.14061302, 0, -4.6562767,
            -0.18906334, -0.50363056, 0.21380206, 0.38380609, -3.0577139
        ],
        [
            1.2089613,   -0.34850798, 0.19992814, 0.5403251,  0, 0.14050086, 0, -4.8026039,
            -0.18825576, -0.52562715, 0.20300422, 0.36676278, -3.1062241
        ],
        [
            1.3862953,   -0.53949006, 0.27391047, 0.54617731, 0, 0.13997674, 0, -5.0392926,
            -0.18698679, -0.56311446, 0.18701872, 0.3413233,  -3.1865087
        ],
        [
            1.6094389,   -0.77067648, 0.36227488, 0.55245584, 0, 0.13887557, 0, -5.3386478,
            -0.18549154, -0.61585933, 0.16916354, 0.31248266, -3.2907297
        ],
        [
            1.7917604,   -0.95299188, 0.43104511, 0.55680682, 0, 0.13773158, 0, -5.5836775,
            -0.18438142, -0.66152202, 0.15627092, 0.29129414, -3.3778599
        ],
        [
            1.9459111,   -1.1031007,  0.48707392, 0.56001288, 0, 0.13665643, 0, -5.7907595,
            -0.18353089, -0.70545901, 0.1464319,  0.27488044, -3.4525609
        ],
        [
            2.0794425,   -1.2304579,  0.53420202, 0.56248121, 0, 0.13567512, 0, -5.9699029,
            -0.18285962, -0.74473416, 0.13861608, 0.26167489, -3.5178624
        ],
        [
            2.3025861,   -1.4383694,  0.61036322, 0.56604894, 0, 0.13398439, 0, -6.2685033,
            -0.18187204, -0.81619739, 0.12685357, 0.24149775, -3.6279124
        ],
        [
            2.7080512,   -1.8029345,  0.74174766, 0.57108126, 0, 0.13092577, 0, -6.8078865,
            -0.18046803, -0.96844739, 0.10894241, 0.21001059, -3.8298007
        ],
        [
            2.9957333,   -2.0531943, 0.83050319,  0.57377093, 0, 0.12887021, 0, -7.1877875,
            -0.17973109, -1.0959388, 0.098399045, 0.19101572, -3.9738828
        ],
        [
            3.4011984,   -2.3966353, 0.95066635,  0.57664578, 0, 0.12622403, 0, -7.7192851,
            -0.17898732, -1.3029921, 0.085881601, 0.16800609, -4.1774755
        ],
        [
            3.912024,    -2.8173563, 1.0957179,   0.57916908, 0, 0.12336471, 0, -8.3828804,
            -0.17840734, -1.6307818, 0.073084038, 0.14397133, -4.4342189
        ],
        [
            4.2486978,   -3.0890439, 1.1883727,   0.5803553, 0, 0.12176624, 0, -8.8170554,
            -0.17828202, -1.9122188, 0.066032748, 0.1299001, -4.607551
        ],
        [
            4.6052479,   -3.372989,  1.2845291,   0.58131461, 0, 0.12029914, 0, -9.274524,
            -0.17810374, -2.2493665, 0.059505919, 0.11751772, -4.7856136
        ],
        [
            5.2984158,   -3.9165388, 1.4671322,  0.58256668, 0, 0.11802917, 0, -10.158307,
            -0.17788629, -3.0938451, 0.04897955, 0.09722763, -5.1320876
        ],
        [
            6.2147425,   -4.6233229, 1.7026936,   0.58346582,  0, 0.11595392, 0, -11.318404,
            -0.17489316, -4.6402724, 0.038273958, 0.076260394, -5.5903087
        ],
        [
            6.9078644,   -5.1522784, 1.8782182,   0.58383065, 0, 0.11489383, 0, -12.191777,
            -0.17783501, -6.6836728, 0.031924402, 0.0637062,  -5.936931
        ],
        [
            7.6009716,   -5.6781273, 2.0524018,   0.58404578,  0, 0.11414425, 0, -13.062807,
            -0.17527475, -9.1401613, 0.026702933, 0.053334939, -6.2835346
        ],
        [
            8.5174408,   -6.370399,  2.2815502,   0.58420059,  0, 0.1134665, 0, -14.212221,
            -0.20546841, -14.295249, 0.021143678, 0.042259114, -6.7418139
        ],
        [
            9.2105225,   -6.8924752, 2.4543776,   0.58426226,  0, 0.11312317, 0, -15.080347,
            -0.17787353, -20.616452, 0.017744896, 0.035475701, -7.0883764
        ],
        [
            10.819834,   -8.1021319, 2.8551568,   0.58431727,  0, 0.11266348, 0, -17.094096,
            -0.17789324, -45.771509, 0.011839675, 0.023676547, -7.8930626
        ],
        [
            11.513098, -8.6226243, 3.0277863,  0.58432046, 0,            0.11255428,
            0,         -17.96111,  -0.1778985, -64.618109, 0.0099509795, 0.019900568,
            -8.2397054
        ],
        [
            13.122387, -9.8302269, 3.4287143,   0.58429751, 0,            0.11240842,
            0,         -19.973234, -0.17790588, -144.11383, 0.0066512024, 0.013302158,
            -9.0443821
        ],
        [
            13.815613, -10.350275, 3.6015262,   0.58427356, 0,            0.11237381,
            0,         -20.839873, -0.20171226, -230.78049, 0.0055922249, 0.011184357,
            -9.3910168
        ],
        [
            16.11828, -12.077458, 4.1759759,   0.58408421, 0,            0.11231658,
            0,        -23.718351, -0.17791072, -643.41944, 0.0031441375, 0.0062884115,
            -10.542508
        ],
        [
            18.42074, -13.804346, 4.7508677,  0.58347816, 0,            0.11229818,
            0,        -26.596504, -0.1779117, -2033.7542, 0.0017677961, 0.0035360667,
            -11.694234
        ]
    ],

    tail_shock_table => [
        [ 4.4381565, -2.1528,    -0.021808543,  -0.021808543 ],
        [ 4.5684523, -2.2576922, -0.025128168,  -0.017366144 ],
        [ 5.2764955, -2.8962679, -0.024987456,  -0.0097999397 ],
        [ 6.2037058, -3.9168572, -0.021983667,  -0.0056703392 ],
        [ 7.5970709, -5.9829668, -0.016912351,  -0.0028256144 ],
        [ 8.515481,  -7.8118434, -0.013958445,  -0.0018613212 ],
        [ 9.2093582, -9.548971,  -0.012052563,  -0.001342088 ],
        [ 10.819487, -14.876178, -0.0083748297, -0.00069650525 ],
        [ 11.512892, -17.958725, -0.0071410411, -0.00052381111 ]
    ],

    blast_info => [
        0.15508241, 0,          0.50771445, -0.30124418, 0,           0,          1.0761297,   -0.66414267,
        0.48193114, 0.17172963, 0.92993572, -0.24733344, 0.098915974, 0.44402967, 0.083440649, 0,
        82.466,     -2.1528,    121.82247,  -2.4615092
    ],

};
1;