package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='P', gamma=1.667
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'P1.667'} = {

    shock_table_info => [ 0, 1.667, 1.13e-06, 8000, 2.1e-07, 25, 1.57e-07, 4.7e-07, 5e-07 ],

    shock_table => [
        [ -11.89911045, 12.0000663,  -0.99999232, -11.90095972, 0.99907452 ],
        [ -11.08146799, 11.18243389, -0.9999826,  -11.0842525,  0.99860583 ],
        [ -9.72600894,  9.82702826,  -0.99992814, -9.73150001,  0.99724714 ],
        [ -8.63033221,  8.73149503,  -0.9997845,  -8.63984739,  0.9952209 ],
        [ -7.76073329,  7.86219463,  -0.99948631, -7.7754667,   0.99258319 ],
        [ -7.04420438,  7.14620326,  -0.99895005, -7.06534645,  0.98932947 ],
        [ -6.4342758,   6.53715579,  -0.99807288, -6.46305111,  0.98543607 ],
        [ -5.90216994,  6.00640104,  -0.99673226, -5.93985311,  0.98087198 ],
        [ -5.42795523,  5.5341632,   -0.99478041, -5.47591048,  0.97558767 ],
        [ -4.99819465,  5.10719333,  -0.99204444, -5.05789509,  0.96952674 ],
        [ -4.60247451,  4.71531326,  -0.98831637, -4.67555552,  0.96261029 ],
        [ -4.23239668,  4.35042778,  -0.98334189, -4.32072675,  0.9547307 ],
        [ -3.88068277,  4.00566357,  -0.97680056, -3.98646959,  0.9457382 ],
        [ -3.53971847,  3.67399729,  -0.96825163, -3.66572136,  0.93539489 ],
        [ -3.20020985,  3.34708905,  -0.95702186, -3.3501473,   0.9232938 ],
        [ -2.84569657,  3.01039616,  -0.94180091, -3.02537553,  0.90854528 ],
        [ -2.4640853,   2.65483553,  -0.92084162, -2.68209551,  0.89011696 ],
        [ -2.11984166,  2.34172143,  -0.89761973, -2.37886945,  0.87122886 ],
        [ -1.80060054,  2.05907705,  -0.87257657, -2.10377386,  0.85194082 ],
        [ -1.50474115,  1.80467776,  -0.84679044, -1.85452671,  0.83277983 ],
        [ -1.22460225,  1.57108557,  -0.82067841, -1.62387841,  0.81377732 ],
        [ -0.94911738,  1.34865442,  -0.79405751, -1.40233137,  0.79458 ],
        [ -0.66996832,  1.13079568,  -0.76683715, -1.18326622,  0.77493611 ],
        [ -0.38127038,  0.913428,    -0.73914567, -0.96246368,  0.75475902 ],
        [ -0.08664585,  0.6996882,   -0.71202355, -0.74307203,  0.73464275 ],
        [ 0.21767968,   0.48704903,  -0.68575283, -0.52256532,  0.71466436 ],
        [ 0.53544319,   0.27320339,  -0.66059871, -0.29864114,  0.69491806 ],
        [ 0.87107568,   0.05556408,  -0.63676795, -0.06870848,  0.6754811 ],
        [ 1.22945181,   -0.1685446,  -0.61445302, 0.16990306,   0.656445 ],
        [ 1.61671631,   -0.40239017, -0.59380221, 0.42046077,   0.63788749 ],
        [ 2.02653226,   -0.64187388, -0.57551,    0.67822853,   0.6204462 ],
        [ 2.43538762,   -0.87399655, -0.5604588,  0.92870442,   0.60514144 ],
        [ 2.85182629,   -1.10469471, -0.54793085, 1.17780451,   0.59150415 ],
        [ 3.27929866,   -1.33661712, -0.53753432, 1.42799222,   0.57932934 ],
        [ 3.72801836,   -1.57579032, -0.52883033, 1.68541023,   0.56829246 ],
        [ 4.19544411,   -1.82124596, -0.52170966, 1.94868138,   0.5584391 ],
        [ 4.68884723,   -2.07716326, -0.51590547, 2.22197446,   0.54959465 ],
        [ 5.2147604,    -2.34719519, -0.51122877, 2.50886438,   0.54165625 ],
        [ 5.7745888,    -2.63231346, -0.50755731, 2.81006686,   0.53461758 ],
        [ 6.37712996,   -2.93723533, -0.50472796, 3.13025511,   0.52838554 ],
        [ 7.02816668,   -3.26510028, -0.50261877, 3.47241241,   0.5229282 ],
        [ 7.74115401,   -3.62287939, -0.50110075, 3.8434904,    0.51817155 ],
        [ 8.52783109,   -4.01664352, -0.50007217, 4.24944716,   0.51408751 ],
        [ 9.40875574,   -4.4568575,  -0.49943457, 4.70072019,   0.51062645 ],
        [ 10.41164153,  -4.95754115, -0.49910362, 5.21129662,   0.50774855 ],
        [ 11.55874472,  -5.52998864, -0.49900614, 5.79233147,   0.50544447 ],
        [ 12.95640742,  -6.22746671, -0.49907908, 6.49738,      0.50359206 ],
        [ 14.69741181,  -7.09652532, -0.49926981, 7.37280379,   0.50218856 ],
        [ 17.18460718,  -8.33865386, -0.49953858, 8.62035517,   0.50111934 ],
        [ 19.23703232,  -9.36409667, -0.49970334, 9.64834868,   0.50065893 ]
    ],

    energy_table => [
        [ -11.89911045, 0.00850866, 0.00339846, 0.00850866, 0.00339846, 11.909067,    -0.8843809 ],
        [ -11.08146799, 0.01179285, 0.00470497, 0.01179285, 0.00470497, 11.162813,    -0.93689714 ],
        [ -9.72600894,  0.02023556, 0.00804458, 0.02023556, 0.00804458, 9.8384932,    -1.0069286 ],
        [ -8.63033221,  0.03124573, 0.01234891, 0.03124573, 0.01234891, 8.7087417,    -1.0507385 ],
        [ -7.76073329,  0.04400375, 0.01724881, 0.04400375, 0.01724881, 7.7814663,    -1.0800997 ],
        [ -7.04420438,  0.0581922,  0.02256812, 0.0581922,  0.02256812, 6.9994112,    -1.1017494 ],
        [ -6.4342758,   0.07361288, 0.02817169, 0.07361288, 0.02817169, 6.3220752,    -1.118438 ],
        [ -5.90216994,  0.09009943, 0.03393229, 0.09009943, 0.03393229, 5.7232707,    -1.1314311 ],
        [ -5.42795523,  0.1075436,  0.03973958, 0.1075436,  0.03973958, 5.1841588,    -1.1413901 ],
        [ -4.99819465,  0.12584145, 0.04548075, 0.12584145, 0.04548075, 4.6918669,    -1.1486435 ],
        [ -4.60247451,  0.14493514, 0.05105334, 0.14493514, 0.05105334, 4.2361811,    -1.1533453 ],
        [ -4.23239668,  0.16481085, 0.05636054, 0.16481085, 0.05636054, 3.8087272,    -1.1555183 ],
        [ -3.88068277,  0.18550947, 0.061308,   0.18550947, 0.061308,   3.4021547,    -1.1550647 ],
        [ -3.53971847,  0.20719185, 0.06580714, 0.20719185, 0.06580714, 3.0086198,    -1.1517337 ],
        [ -3.20020985,  0.23022388, 0.06976405, 0.23022388, 0.06976405, 2.618424,     -1.1450006 ],
        [ -2.84569657,  0.25557138, 0.07307188, 0.25557138, 0.07307188, 2.214096,     -1.1337005 ],
        [ -2.4640853,   0.28394426, 0.07539018, 0.28394426, 0.07539018, 1.7842627,    -1.1162255 ],
        [ -2.11984166,  0.31006805, 0.07616411, 0.31006805, 0.07616411, 1.4031577,    -1.0953481 ],
        [ -1.80060054,  0.33433558, 0.07567236, 0.33433558, 0.07567236, 1.0569507,    -1.0715335 ],
        [ -1.50474115,  0.35652757, 0.07418361, 0.35652757, 0.07418361, 0.74347412,   -1.0444624 ],
        [ -1.22460225,  0.37701118, 0.07192617, 0.37701118, 0.07192617, 0.45488046,   -1.013792 ],
        [ -0.94911738,  0.39643843, 0.06900966, 0.39643843, 0.06900966, 0.1800359,    -0.98298082 ],
        [ -0.66996832,  0.41522229, 0.06548704, 0.41522229, 0.06548704, -0.090205647, -0.95633612 ],
        [ -0.38127038,  0.43354778, 0.06140526, 0.43354778, 0.06140526, -0.36278837,  -0.93594561 ],
        [ -0.08664585,  0.45098836, 0.05695231, 0.45098836, 0.05695231, -0.63606515,  -0.9228475 ],
        [ 0.21767968,   0.46760162, 0.05221913, 0.46760162, 0.05221913, -0.91543535,  -0.9162183 ],
        [ 0.53544319,   0.48340992, 0.04729275, 0.48340992, 0.04729275, -1.2059856,   -0.91483781 ],
        [ 0.87107568,   0.49843041, 0.04225049, 0.49843041, 0.04225049, -1.5132043,   -0.91746125 ],
        [ 1.22945181,   0.51265077, 0.03716952, 0.51265077, 0.03716952, -1.8428113,   -0.92304635 ],
        [ 1.61671631,   0.52605188, 0.03212129, 0.52605188, 0.03212129, -2.2016648,   -0.93089706 ],
        [ 2.02653226,   0.53821429, 0.02733058, 0.53821429, 0.02733058, -2.5850089,   -0.93767343 ],
        [ 2.43538762,   0.54850945, 0.02312468, 0.54850945, 0.02312468, -2.9693068,   -0.97079421 ],
        [ 2.85182629,   0.55734606, 0.01940623, 0.55734606, 0.01940623, -3.3866719,   -1.0278278 ],
        [ 3.27929866,   0.56492424, 0.01613763, 0.56492424, 0.01613763, -3.8372745,   -1 ],
        [ 3.72801836,   0.57149667, 0.01324238, 0.57149667, 0.01324238, -4.2859942,   -1 ],
        [ 4.19544411,   0.57708204, 0.01073653, 0.57708204, 0.01073653, -4.75342,     -1 ],
        [ 4.68884723,   0.58182719, 0.00857347, 0.58182719, 0.00857347, -5.2468231,   -1 ],
        [ 5.2147604,    0.58583076, 0.00672265, 0.58583076, 0.00672265, -5.7727363,   -1 ],
        [ 5.7745888,    0.58914234, 0.00517267, 0.58914234, 0.00517267, -6.3325647,   -1 ],
        [ 6.37712996,   0.59185473, 0.00388924, 0.59185473, 0.00388924, -6.9351058,   -1 ],
        [ 7.02816668,   0.59403128, 0.00284946, 0.59403128, 0.00284946, -7.5861426,   -1 ],
        [ 7.74115401,   0.59575108, 0.002021,   0.59575108, 0.002021,   -8.2991299,   -1 ],
        [ 8.52783109,   0.59707293, 0.0013796,  0.59707293, 0.0013796,  -9.085807,    -1 ],
        [ 9.40875574,   0.59806081, 0.00089726, 0.59806081, 0.00089726, -9.9667316,   -1 ],
        [ 10.41164153,  0.59877157, 0.00054839, 0.59877157, 0.00054839, -10.969617,   -1 ],
        [ 11.55874472,  0.59925209, 0.00031148, 0.59925209, 0.00031148, -12.116721,   -1 ],
        [ 12.95640742,  0.59956637, 0.00015595, 0.59956637, 0.00015595, -13.514383,   -1 ],
        [ 14.69741181,  0.59974817, 6.57e-05,   0.59974817, 6.57e-05,   -15.255388,   -1 ],
        [ 17.18460718,  0.5998419,  1.905e-05,  0.5998419,  1.905e-05,  -17.742583,   -1 ],
        [ 19.23703232,  0.59986637, 6.84e-06,   0.59986637, 6.84e-06,   -19.795008,   -1 ]
    ],

    impulse_table => [
        [
            -13.714451, 13.815401,  -13.715197,    0.66294313, 0,          0,
            0,          -7.9917105, 0.00015514751, -20.704038, 0.22720456, 0.99313346,
            -5.8983186
        ],
        [
            -2.3025851, 2.5069459,  -2.5390297,    0.49757243, 0,          0,
            0,          -2.4886302, 0.00015514879, -20.604039, 0.20643174, 0.50629093,
            -1.6219396
        ],
        [
            -1.6094379, 1.8938271,  -1.9420787,    0.46960004, 0,          0,
            0,          -2.3329436, 0.00015516027, -20.504039, 0.19033855, 0.41867316,
            -1.4586058
        ],
        [
            -1.3862943, 1.7050201,  -1.7563563,    0.46488491, 0,          0,
            0,          -2.3089323, 0.00015517296, -20.454039, 0.18351893, 0.3911035,
            -1.4122862
        ],
        [
            -1.2039728, 1.5541757,  -1.6071053,    0.46302096, 0,          0,
            0,          -2.2995976, 0.00015519115, -20.404039, 0.17734157, 0.36905231,
            -1.3767266
        ],
        [
            -0.91629069, 1.3226408,  -1.3762858,    0.46373404, 0,          0,
            0,           -2.3031606, 0.00015524475, -20.304039, 0.16655315, 0.33537143,
            -1.3247507
        ],
        [
            -0.69314714, 1.1485962,  -1.2012472,    0.46710899, 0,          0,
            0,           -2.3201919, 0.00015532149, -20.204039, 0.15741728, 0.31036089,
            -1.2878394
        ],
        [
            -0.51082558,  1.0099849,  -1.0608301, 0.47141525, 0, 0, 0, -2.342353,
            0.0001554211, -20.104039, 0.14956011, 0.29074128, -1.2598157
        ],
        [
            -0.3566749, 0.8952769,  -0.94392098,   0.47595419, 0,          0,
            0,          -2.3662562, 0.00015554327, -20.004039, 0.14271753, 0.27476571,
            -1.2375693
        ],
        [
            -0.22314351, 0.79771965, -0.84397821,   0.48043225, 0,          0,
            0,           -2.3904134, 0.00015568775, -19.904039, 0.13669594, 0.26139975,
            -1.2193349
        ],
        [
            4.1604174e-08, 0.63832825, -0.67966938,   0.48879423, 0,          0,
            0,             -2.4371503, 0.00015604317, -19.704039, 0.12656543, 0.24007307,
            -1.1909241
        ],
        [
            0.26236431, 0.45648895, -0.49069461,   0.49959448, 0,          0,
            0,          -2.5009502, 0.00015674175, -19.404039, 0.11478961, 0.21664436,
            -1.1606408
        ],
        [
            0.53062829, 0.27638498, -0.3019878,    0.51121987, 0,          0,
            0,          -2.5745246, 0.00015798383, -19.004039, 0.10319577, 0.19453559,
            -1.132916
        ],
        [
            0.69314722, 0.16994359, -0.18978927,   0.518367,   0,          0,
            0,          -2.6226191, 0.00015915206, -18.704039, 0.09649692, 0.18203884,
            -1.1176007
        ],
        [
            0.99325181, -0.021743393, 0.013408272,   0.53147931, 0,           0,
            0,          -2.7174324,   0.00016269528, -18.004039, 0.084912361, 0.16069259,
            -1.0920275
        ],
        [
            1.0986123, -0.087645713, 0.083575731,   0.53600166, 0,           0,
            0,         -2.7523702,   0.00016457924, -17.704039, 0.081105433, 0.15371496,
            -1.0838303
        ],
        [
            1.2089604, -0.15594141, 0.15644092,    0.54067028, 0,           0,
            0,         -2.7897955,  0.00016706923, -17.354039, 0.077267345, 0.14668563,
            -1.0756553
        ],
        [
            1.3862944, -0.26422103, 0.2722482,     0.54799769, 0,           0,
            0,         -2.8516041,  0.00017257945, -16.704039, 0.071418669, 0.13596961,
            -1.063359
        ],
        [
            1.609438, -0.39806697, 0.4158168,     0.5568533,  0,           0,
            0,        -2.9320476,  0.00018356203, -15.704039, 0.064608413, 0.12346175,
            -1.0492823
        ],
        [
            1.7917595, -0.50560133, 0.53143733,    0.56374012, 0,           0,
            0,         -2.999777,   0.00019813289, -14.704039, 0.059481829, 0.11400996,
            -1.0388742
        ],
        [
            1.9459102, -0.59534188, 0.62807611,    0.56928776, 0,           0,
            0,         -3.058325,   0.00021716948, -13.704039, 0.055440155, 0.10652982,
            -1.0308131
        ],
        [
            2.0794416, -0.67226757, 0.7110004,     0.57386768, 0,          0,
            0,         -3.1099252,  0.00024194312, -12.704039, 0.05214595, 0.10041155,
            -1.0243731
        ],
        [
            2.3025851, -0.7992666, 0.84802694,    0.58097793, 0,           0,
            0,         -3.1978498, 0.00031710146, -10.704039, 0.047046024, 0.090897118,
            -1.0147879
        ],
        [
            2.7080502, -1.0256292, 1.0924397,    0.59092902, 0,           0,
            0,         -3.3624908, 0.0009960365, -5.7040383, 0.038962997, 0.075698896,
            -1.0031007
        ],
        [
            2.9957323,  -1.1832721,  1.2626143,   0.57236643, 0, 0, 0, -3.4853656,
            0.02284138, -0.70403813, 0.034051358, 0.06606542, -1.0525944
        ]
    ],

    blast_info => [ 0, 0, 0, 0, 0, 0, 0.66274558, 0, 0, 0, 0, 0, 0.35315806, 0.22720457, 0, 0, 0, 0, 0, 0 ],

};
1;
