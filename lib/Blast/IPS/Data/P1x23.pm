package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='P', gamma=1.23
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'P1x23'} = {

    table_name  => 'P1x23',
    symmetry    => 0,
    gamma       => 1.23,
    data_source => 'P8000_G1x23/moc2_from_moc_r1e5/',

    shock_table_info => [ 0, 1.23, 1.25e-06, 8000, 6.4e-07, 100, 4.62e-07, 2.9e-07, 5e-07, 1280.1, -3.9101 ],

    shock_table => [
        [ -12.8889881,  12.00017792,  -0.99999322, -12.89072508, 0.99913076 ],
        [ -11.85446223, 10.96566706,  -0.99998093, -11.85737759, 0.99854023 ],
        [ -11.14676031, 10.25798907,  -0.99995206, -11.15091588, 0.997918 ],
        [ -9.88508243,  8.99643109,   -0.99983287, -9.89290534,  0.99607392 ],
        [ -8.97210008,  8.08369813,   -0.99958363, -8.98447514,  0.99377685 ],
        [ -8.22127989,  7.33334341,   -0.99911918, -8.23933944,  0.99089695 ],
        [ -7.5882538,   6.70109441,   -0.99834519, -7.61311107,  0.98743834 ],
        [ -7.03708991,  6.15114489,   -0.99713962, -7.06994307,  0.98335282 ],
        [ -6.54841626,  5.66426854,   -0.99536383, -6.59051405,  0.97861088 ],
        [ -6.1072352,   5.22564796,   -0.99285137, -6.15992579,  0.97316048 ],
        [ -5.70198601,  4.82395067,   -0.98939953, -5.76677259,  0.9669249 ],
        [ -5.32423253,  4.45102579,   -0.98476805, -5.40281323,  0.95981309 ],
        [ -4.96616348,  4.09944834,   -0.97864969, -5.06054526,  0.95168541 ],
        [ -4.61996757,  3.76196497,   -0.97062684, -4.73265025,  0.94232453 ],
        [ -4.27627127,  3.43009683,   -0.96006629, -4.41061063,  0.9313587 ],
        [ -3.92007692,  3.09055783,   -0.94580584, -4.08117733,  0.91803258 ],
        [ -3.52638844,  2.72203667,   -0.92550459, -3.72306287,  0.90079638 ],
        [ -3.1726937,   2.398564,     -0.90291225, -3.40753794,  0.88300866 ],
        [ -2.84586976,  2.10737453,   -0.87848268, -3.12189268,  0.8647142 ],
        [ -2.54185476,  1.84410313,   -0.85310091, -2.86177826,  0.84627077 ],
        [ -2.2550926,   1.60312051,   -0.82737922, -2.62172094,  0.82784563 ],
        [ -1.97388654,  1.37413188,   -0.80112515, -2.39155208,  0.80907402 ],
        [ -1.69018087,  1.15065429,   -0.77429278, -2.16475178,  0.78972478 ],
        [ -1.39790855,  0.92835638,   -0.74698982, -1.93686736,  0.76968138 ],
        [ -1.09992665,  0.70979732,   -0.72015446, -1.7105386,   0.74945706 ],
        [ -0.79318506,  0.49293515,   -0.6941244,  -1.4837778,   0.72916967 ],
        [ -0.47364402,  0.27518919,   -0.66912644, -1.25404511,  0.70888651 ],
        [ -0.13721957,  0.0541485,    -0.6453821,  -1.01898986,  0.68870217 ],
        [ 0.2209029,    -0.17289201,  -0.62307265, -0.77597786,  0.66870911 ],
        [ 0.6065538,    -0.409078,    -0.60234832, -0.52195108,  0.64900176 ],
        [ 1.02668396,   -0.65802808,  -0.58334894, -0.25341974,  0.62969414 ],
        [ 1.44267374,   -0.8973538,   -0.56778376, 0.00492391,   0.61271781 ],
        [ 1.86424805,   -1.13389587,  -0.55483651, 0.25996262,   0.59754351 ],
        [ 2.29784483,   -1.37204193,  -0.5440153,  0.51602625,   0.58388306 ],
        [ 2.74835631,   -1.61501644,  -0.53498452, 0.77622997,   0.57156503 ],
        [ 3.22007213,   -1.86553871,  -0.527489,   1.04316276,   0.56047595 ],
        [ 3.71708544,   -2.12610975,  -0.52132067, 1.31918756,   0.55053495 ],
        [ 4.24415852,   -2.39949946,  -0.51629588, 1.60695362,   0.54167046 ],
        [ 4.80690645,   -2.68885056,  -0.51225215, 1.90949707,   0.53382013 ],
        [ 5.41293297,   -2.99826456,  -0.50904096, 2.23084166,   0.52692077 ],
        [ 6.06412067,   -3.32889078,  -0.50655589, 2.57195753,   0.52098128 ],
        [ 6.77848341,   -3.69003088,  -0.50464546, 2.94222609,   0.51587837 ],
        [ 7.56596937,   -4.08683132,  -0.50321573, 3.3467068,    0.51159931 ],
        [ 8.44496717,   -4.5286627,   -0.50216769, 3.7947766,    0.50809558 ],
        [ 9.44784887,   -5.03186721,  -0.50140921, 4.30284551,   0.50530346 ],
        [ 10.60830028,  -5.61338927,  -0.50087012, 4.8879064,    0.50319413 ],
        [ 12.02924956,  -6.3247986,   -0.50048189, 5.60173164,   0.50167408 ],
        [ 13.76927909,  -7.19540354,  -0.50022577, 6.47373668,   0.50073521 ],
        [ 15.80421478,  -8.21316986,  -0.50008769, 7.49215544,   0.5002728 ],
        [ 18.33531773,  -9.47884937,  -0.50002611, 8.75810095,   0.50007784 ],
        [ 21.65553988,  -11.13900256, -0.50000493, 10.41833808,  0.5000148 ]
    ],

    energy_table => [
        [ -12.8889881,  0.10454242, 0.0195409,  0.10454242, 0.0195409,     11.523817,   -1.0017181 ],
        [ -11.85446223, 0.1268396,  0.02369691, 0.1268396,  0.02369691,    10.48708,    -1.0028245 ],
        [ -11.14676031, 0.14476344, 0.0270279,  0.14476344, 0.0270279,     9.7770465,   -1.0040786 ],
        [ -9.88508243,  0.18316931, 0.03411239, 0.18316931, 0.03411239,    8.5084584,   -1.0077643 ],
        [ -8.97210008,  0.21705171, 0.04026004, 0.21705171, 0.04026004,    7.5868763,   -1.0121382 ],
        [ -8.22127989,  0.24939064, 0.04598016, 0.24939064, 0.04598016,    6.8252637,   -1.0176862 ],
        [ -7.5882538,   0.28013849, 0.05122204, 0.28013849, 0.05122204,    6.1792742,   -1.0243909 ],
        [ -7.03708991,  0.30968293, 0.05601007, 0.30968293, 0.05601007,    5.6127894,   -1.0323907 ],
        [ -6.54841626,  0.33810606, 0.06031523, 0.33810606, 0.06031523,    5.1062971,   -1.041784 ],
        [ -6.1072352,   0.3655609,  0.06411989, 0.3655609,  0.06411989,    4.6445625,   -1.0526594 ],
        [ -5.70198601,  0.3922207,  0.06740633, 0.3922207,  0.06740633,    4.2157129,   -1.0651679 ],
        [ -5.32423253,  0.41821346, 0.07014641, 0.41821346, 0.07014641,    3.8109014,   -1.0795812 ],
        [ -4.96616348,  0.44373315, 0.07231175, 0.44373315, 0.07231175,    3.4216377,   -1.0958364 ],
        [ -4.61996757,  0.46905329, 0.07386412, 0.46905329, 0.07386412,    3.0393458,   -1.1145813 ],
        [ -4.27627127,  0.49461172, 0.07474147, 0.49461172, 0.07474147,    2.6527481,   -1.1353523 ],
        [ -3.92007692,  0.52127499, 0.07481941, 0.52127499, 0.07481941,    2.2444554,   -1.1582565 ],
        [ -3.52638844,  0.55056876, 0.07379558, 0.55056876, 0.07379558,    1.7832448,   -1.1825718 ],
        [ -3.1726937,   0.57635137, 0.07182929, 0.57635137, 0.07182929,    1.3614608,   -1.1924924 ],
        [ -2.84586976,  0.5994113,  0.06915378, 0.5994113,  0.06915378,    0.97173066,  -1.1795964 ],
        [ -2.54185476,  0.61997014, 0.06599554, 0.61997014, 0.06599554,    0.61675834,  -1.1367403 ],
        [ -2.2550926,   0.63840709, 0.0625201,  0.63840709, 0.0625201,     0.29913537,  -1.057045 ],
        [ -1.97388654,  0.65546528, 0.0587514,  0.65546528, 0.0587514,     0.015833851, -0.97323485 ],
        [ -1.69018087,  0.67156359, 0.05470321, 0.67156359, 0.05470321,    -0.25048459, -0.90857025 ],
        [ -1.39790855,  0.6869244,  0.05039678, 0.6869244,  0.05039678,    -0.50695836, -0.86270528 ],
        [ -1.09992665,  0.70128345, 0.04598423, 0.70128345, 0.04598423,    -0.75952924, -0.84546459 ],
        [ -0.79318506,  0.71470198, 0.04152955, 0.71470198, 0.04152955,    -1.0181926,  -0.84650352 ],
        [ -0.47364402,  0.7272552,  0.03708,    0.7272552,  0.03708,       -1.2897644,  -0.85581491 ],
        [ -0.13721957,  0.73898128, 0.03268462, 0.73898128, 0.03268462,    -1.5797833,  -0.86940972 ],
        [ 0.2209029,    0.74990482, 0.02838894, 0.74990482, 0.02838894,    -1.8939394,  -0.88510543 ],
        [ 0.6065538,    0.76003641, 0.0242367,  0.76003641, 0.0242367,     -2.2385514,  -0.90126879 ],
        [ 1.02668396,   0.7693664,  0.020274,   0.7693664,  0.020274,      -2.6207182,  -0.91701392 ],
        [ 1.44267374,   0.77708004, 0.01689937, 0.77708004, 0.01689937,    -3.0052243,  -0.93058912 ],
        [ 1.86424805,   0.78357403, 0.01399013, 0.78357403, 0.01399013,    -3.4002168,  -0.94233279 ],
        [ 2.29784483,   0.78907856, 0.01147549, 0.78907856, 0.01147549,    -3.8112114,  -0.95250918 ],
        [ 2.74835631,   0.79374453, 0.00930871, 0.79374453, 0.00930871,    -4.2424979,  -0.96130675 ],
        [ 3.22007213,   0.79768295, 0.00745424, 0.79768295, 0.00745424,    -4.6979275,  -0.96887099 ],
        [ 3.71708544,   0.80098244, 0.00588217, 0.80098244, 0.00588217,    -5.1812479,  -0.97532662 ],
        [ 4.24415852,   0.80372126, 0.00456404, 0.80372126, 0.00456404,    -5.6969244,  -0.98075498 ],
        [ 4.80690645,   0.80596909, 0.00347293, 0.80596909, 0.00347293,    -6.2502707,  -0.98527452 ],
        [ 5.41293297,   0.80779087, 0.00258221, 0.80779087, 0.00258221,    -6.8486661,  -0.98900447 ],
        [ 6.06412067,   0.80922984, 0.00187441, 0.80922984, 0.00187441,    -7.4938107,  -0.99197849 ],
        [ 6.77848341,   0.81035808, 0.0013167,  0.81035808, 0.0013167,     -8.2034295,  -0.99433485 ],
        [ 7.56596937,   0.81121644, 0.00089071, 0.81122099, 0.00089565653, -8.9873012,  -0.9961418 ],
        [ 8.44496717,   0.8118506,  0.00057503, 0.81185332, 0.00057503,    -9.8636253,  -0.99749253 ],
        [ 9.44784887,   0.81230438, 0.00034864, 0.81230687, 0.00034864,    -10.864604,  -0.99847707 ],
        [ 10.60830028,  0.81261141, 0.00019524, 0.81261341, 0.00019524,    -12.023791,  -0.99934571 ],
        [ 12.02924956,  0.81281001, 9.594e-05,  0.81281421, 9.657851e-05,  -13.444567,  -1 ],
        [ 13.76927909,  0.81292149, 4.019e-05,  0.81292609, 4.019e-05,     -15.184597,  -1 ],
        [ 15.80421478,  0.81297281, 1.452e-05,  0.81297516, 1.452e-05,     -17.219533,  -1 ],
        [ 18.33531773,  0.81299366, 4.1e-06,    0.81299366, 4.1e-06,       -19.750636,  -1 ],
        [ 21.65553988,  0.8130003,  7.8e-07,    0.8130003,  7.8e-07,       -23.070858,  -1 ]
    ],

    impulse_table => [
        [
            -14.703705, 13.814889,  -14.704405,    0.25583394, 0,          -0.73433652,
            0,          -8.6906861, -0.0044044976, -1.2321453, 0.13638211, 0.90840883,
            -4.0635137
        ],
        [
            -13.815512, 12.926698,  -13.816605,    0.25564593, 0,         -0.73433593,
            0,          -8.2465927, -0.0044044976, -1.2321447, 0.1363819, 0.89186326,
            -3.8974498
        ],
        [
            -13.122365, 12.233554,  -13.123911,    0.25542876, 0,          -0.73433493,
            0,          -7.9000229, -0.0044044976, -1.2321437, 0.13638146, 0.87690246,
            -3.7678687
        ],
        [
            -12.7169, 11.828091,  -12.718793,    0.25526212, 0,          -0.73433393,
            0,        -7.6972939, -0.0044044976, -1.2321427, 0.13638105, 0.86721029,
            -3.6920789
        ],
        [
            -12.206074, 11.317273, -12.208519,    0.25499787, 0,          -0.73433193,
            0,          -7.441888, -0.0044044976, -1.2321407, 0.13638024, 0.85390832,
            -3.5966114
        ],
        [
            -11.512927, 10.624141,  -11.516386,    0.2545123,  0,          -0.73432693,
            0,          -7.0953314, -0.0044044976, -1.2321357, 0.13637827, 0.83371026,
            -3.467114
        ],
        [
            -10.81978, 9.931027,   -10.824674,    0.25382564, 0,          -0.73431693,
            0,         -6.7487917, -0.0044044976, -1.2321257, 0.13637435, 0.81073537,
            -3.3376993
        ],
        [
            -10.414315, 9.5255945, -10.420311,    0.2532988,  0,          -0.73430693,
            0,          -6.546093, -0.0044044976, -1.2321157, 0.13637043, 0.79586352,
            -3.2620563
        ],
        [
            -9.9034894, 9.014835,  -9.9112404,    0.25246348, 0,          -0.73428693,
            0,          -6.290748, -0.0044044976, -1.2320957, 0.13636259, 0.77547363,
            -3.166851
        ],
        [
            -9.2103422, 8.3218517, -9.2213197,    0.25092899, 0,          -0.73423693,
            0,          -5.944344, -0.0044044976, -1.2320457, 0.13634299, 0.74457578,
            -3.0379159
        ],
        [
            -8.517195, 7.6290324,  -8.5327525,    0.24876051, 0,         -0.73413693,
            0,         -5.5981097, -0.0044044976, -1.2319457, 0.1363038, 0.70956045,
            -2.9094412
        ],
        [
            -8.1117299, 7.2238953,  -8.1308148,    0.24709829, 0,          -0.73403693,
            0,          -5.3957164, -0.0044044976, -1.2318457, 0.13626463, 0.6869954,
            -2.834614
        ],
        [
            -7.6009043, 6.7137241,  -7.6256032,    0.24446636, 0,          -0.73383693,
            0,          -5.1409824, -0.0044044976, -1.2316457, 0.13618639, 0.65622667,
            -2.7408404
        ],
        [
            -6.9077571, 6.0222065,  -6.942836,     0.23964653, 0,          -0.73333693,
            0,          -4.7961064, -0.0044044976, -1.2311457, 0.13599132, 0.61009014,
            -2.6148867
        ],
        [
            -6.2146099, 5.3322938,  -6.2644967,    0.23288148, 0,          -0.73233693,
            0,          -4.4529314, -0.0044044976, -1.2301457, 0.13560351, 0.55876923,
            -2.4911972
        ],
        [
            -5.8091448, 4.9300302,  -5.8704831,    0.22774435, 0,          -0.73133694,
            0,          -4.2536018, -0.0044044976, -1.2291457, 0.13521878, 0.52639862,
            -2.4203651
        ],
        [
            -5.2983192, 4.425512,   -5.3779483,    0.21972198, 0,          -0.72933694,
            0,          -4.0050085, -0.0044044976, -1.2271457, 0.13445839, 0.48335956,
            -2.3333206
        ],
        [
            -4.605172, 3.7476069,  -4.7187113,    0.20550021, 0,          -0.72433695,
            0,         -3.6755674, -0.0044044976, -1.2221457, 0.13260815, 0.42172033,
            -2.2204656
        ],
        [
            -3.9120248, 3.0829436,  -4.0737866,    0.18695378, 0,          -0.71433706,
            0,          -3.3636633, -0.0044044975, -1.2121458, 0.12910758, 0.35809122,
            -2.1158143
        ],
        [
            -3.5065597, 2.7036966,  -3.7052106,    0.17433692, 0,          -0.70433824,
            0,          -3.1962177, -0.0044044953, -1.2021458, 0.12584438, 0.32100192,
            -2.0593965
        ],
        [
            -2.9957341, 2.2399058,  -3.2521314,    0.15793384, 0,          -0.68436901,
            0,          -3.0135269, -0.0044044389, -1.1821458, 0.11992038, 0.27559613,
            -1.9942005
        ],
        [
            -2.3025869, 1.64252,    -2.6611129,    0.14177465, 0,          -0.63665284,
            0,          -2.8613401, -0.0044002479, -1.1348468, 0.10776892, 0.21842878,
            -1.9171607
        ],
        [
            -1.6094398, 1.0884446, -2.1012123,    0.14477167, 0,           -0.58462162,
            0,          -2.890018, -0.0043112843, -1.0807019, 0.090623011, 0.16860112,
            -1.8533889
        ],
        [
            -1.2039746, 0.7852064,  -1.7888829,   0.15328435, 0,           -0.57467776,
            0,          -2.9738603, -0.004148616, -1.0703446, 0.078974992, 0.14340211,
            -1.8217562
        ],
        [
            -0.69314902, 0.4239034,  -1.4111583,    0.16599205, 0,           -0.58191354,
            0,           -3.1163799, -0.0038232418, -1.0828147, 0.064115085, 0.11589023,
            -1.7870757
        ],
        [
            -1.8406155e-06, -0.033793195, -0.92502895,   0.18278862, 0,           -0.60830755,
            0,              -3.3491631,   -0.0032868773, -1.1210732, 0.046465188, 0.085674945,
            -1.7478309
        ],
        [
            0.26236242, -0.1986747, -0.74829939,   0.18856363, 0,           -0.6190696,
            0,          -3.4455693, -0.0030950448, -1.1355107, 0.040937706, 0.076175361,
            -1.735031
        ],
        [
            0.40546327, -0.28692153, -0.65345972,   0.19154824, 0,           -0.62462681,
            0,          -3.4997089,  -0.0029982883, -1.1429407, 0.038190475, 0.071400989,
            -1.7284895
        ],
        [
            0.53062641, -0.36319911, -0.57136918,   0.1940601,  0,           -0.62921954,
            0,          -3.5478892,  -0.0029191043, -1.1479297, 0.035934133, 0.067448424,
            -1.7230148
        ],
        [
            0.69314534, -0.46105315, -0.46593389,   0.19718306, 0,           -0.63475206,
            0,          -3.6115291,  -0.0028244663, -1.1531874, 0.033197292, 0.06261337,
            -1.7162416
        ],
        [
            0.99324993, -0.63850132, -0.27449723,   0.20254023, 0,         -0.64353529,
            0,          -3.7319946,  -0.0026749687, -1.1633774, 0.0286701, 0.054512085,
            -1.7046965
        ],
        [
            1.0986104, -0.69988128, -0.2082393,    0.20429698, 0,           -0.64616465,
            0,         -3.7751268,  -0.0026301717, -1.167017,  0.027229413, 0.051906175,
            -1.7009293
        ],
        [
            1.2089585, -0.76369376, -0.13934968,   0.20606933, 0,           -0.64867549,
            0,         -3.8207369,  -0.0025873699, -1.1686066, 0.025797015, 0.049301677,
            -1.6971377
        ],
        [
            1.3839038, -0.86392623, -0.031152042,  0.2087407,  0,           -0.6521531,
            0,         -3.8939165,  -0.0025276601, -1.1707747, 0.023677018, 0.045421825,
            -1.6914405
        ],
        [
            1.6094361, -0.99157995, 0.10657787,    0.21194222, 0,           -0.65581343,
            0,         -3.9897297,  -0.0024641579, -1.1761264, 0.021196129, 0.040843059,
            -1.6846418
        ],
        [
            1.7917576, -1.0936021, 0.21655727,    0.21433948, 0,           -0.65816749,
            0,         -4.0683114, -0.0024225699, -1.1796391, 0.019379858, 0.037464412,
            -1.6795736
        ],
        [
            1.9459083, -1.179113,  0.30864683,    0.21623935, 0,           -0.65979173,
            0,         -4.1354824, -0.0023932866, -1.1794799, 0.017964945, 0.034816748,
            -1.6755721
        ],
        [
            2.0794397, -1.2526731, 0.38778604,    0.21779495, 0,           -0.66097255,
            0,         -4.1941784, -0.0023717341, -1.1834412, 0.016822232, 0.032668415,
            -1.6723061
        ],
        [
            2.3025833, -1.3746195, 0.51879261,    0.22021672, 0,           -0.66253156,
            0,         -4.2932551, -0.0023424004, -1.1863363, 0.015071054, 0.029358707,
            -1.6672427
        ],
        [
            2.7080484, -1.5934378, 0.75317067,    0.22409057, 0,          -0.66434286,
            0,         -4.4761612, -0.0023057628, -1.1962127, 0.01233858, 0.024152062,
            -1.6592016
        ],
        [
            2.9957304, -1.7468337, 0.91686459,    0.22646641, 0,           -0.66508183,
            0,         -4.6079351, -0.0022890163, -1.2026208, 0.010703513, 0.021011759,
            -1.6543123
        ],
        [
            3.4011955, -1.9608561, 1.1443292,     0.22935301, 0,            -0.66562786,
            0,         -4.7960987, -0.0022738326, -1.2351378, 0.0087575409, 0.017250106,
            -1.6484222
        ],
        [
            3.9120212, -2.227534,  1.4261653,     0.23232662, 0,            -0.66585958,
            0,         -5.0366149, -0.0022630556, -1.258319,  0.0067982563, 0.013435907,
            -1.6424343
        ],
        [
            4.2484934, -2.4017375, 1.6093016,     0.23394167, 0,            -0.66605058,
            0,         -5.1968111, -0.0022585935, -1.2953216, 0.0057523392, 0.01138854,
            -1.6392474
        ],
        [
            4.6056894, -2.5856474, 1.8018218,    0.23540707, 0,            -0.66567215,
            0,         -5.3681741, -0.002257429, -1.3699944, 0.0048165648, 0.0091453114,
            -1.6749075
        ],
        [
            5.298356, -2.9399102, 2.1704009,     0.23763642, 0,            -0.66549421,
            0,        -5.7036041, -0.0022551597, -1.5744308, 0.0034119298, 0.0065761674,
            -1.6593553
        ],
        [
            6.2146244, -3.4050942, 2.6502773,     0.23965573, 0,            -0.66526762,
            0,         -6.1520842, -0.0022541955, -2.1888821, 0.0021606019, 0.0042208028,
            -1.6453449
        ],
        [
            6.9078833, -3.755314,  3.0089293,     0.24068806, 0,            -0.66513599,
            0,         -6.4939469, -0.0022539784, -3.213644,  0.0015284362, 0.0030063478,
            -1.6382109
        ],
        [
            7.6010326, -4.1044748, 3.3646423,     0.24142367, 0,            -0.66503833,
            0,         -6.8372655, -0.0022005277, -5.0748435, 0.0010809986, 0.0021366221,
            -1.6331387
        ],
        [
            8.5173159, -4.5649915, 3.831528,     0.24208002, 0,             -0.66495039,
            0,         -7.2926872, -0.002228853, -11.235694, 0.00068369273, 0.0013571984,
            -1.628619
        ],
        [
            9.2104952, -4.9128388, 4.1828443,     0.2424115, 0,             -0.66449514,
            0,         -7.638021,  -0.0022538988, -21.66085, 0.00048338056, 0.00096166014,
            -1.6263347
        ],
        [
            10.819842, -5.7193366, 4.9943226,     0.24284843, 0,             -0.66041281,
            0,         -8.4412641, -0.0022539169, -103.63864, 0.00021609494, 0.00043117313,
            -1.6232825
        ]
    ],

    tail_shock_table => [
        [ 7.1577434, -3.9101,    -0.0014541012,  -0.0014541012 ],
        [ 7.5864673, -5.0111695, -0.0018000105,  -0.00082400821 ],
        [ 8.5080476, -8.4637866, -0.0013786924,  -0.00034302528 ],
        [ 9.2039194, -12.47236,  -0.0010551518,  -0.00018887683 ],
        [ 10.816886, -35.671245, -0.00062876335, -3.0202965e-05 ],
        [ 11.511054, -44.091904, -0.00038994753, -2.7412098e-05 ],
        [ 13.12145,  -104.74797, -0.00018684781, -1.9538883e-05 ]
    ],

    blast_info => [
        0.73433423, 0,            1.2326143, -0.0044045104, 0,         0,          0.25580292, -0.021831962,
        3.0210714,  0.0016653099, 1.7063155, -0.0038186008, 0.4247196, 0.13638186, 0.19756671, 0,
        1280.1,     -3.9101,      1557.4461, -4.3551098
    ],

};
1;
