package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='P', gamma=1.3
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'P1x3'} = {

    table_name  => 'P1x3',
    symmetry    => 0,
    gamma       => 1.3,
    data_source => 'P8000_G1x3_inspirion/moc2_from_moc_r1e5/',

    shock_table_info => [ 0, 1.3, 1.11e-06, 8000, 2.32e-07, 34, 1.7e-07, 4.4e-07, 5e-07, 4675.2, -6.758 ],

    shock_table => [
        [ -12.64178695, 12.00006627,  -0.99999305, -12.64354541, 0.99912001 ],
        [ -11.55327663, 10.91157197,  -0.99997938, -11.55630893, 0.99848159 ],
        [ -10.65825935, 10.01659,     -0.99994115, -10.66300708, 0.99762065 ],
        [ -9.59080036,  8.94924522,   -0.99982596, -9.59890966,  0.99592965 ],
        [ -8.68066644,  8.03936958,   -0.99956787, -8.69347716,  0.99355652 ],
        [ -7.93379227,  7.29297481,   -0.99908955, -7.95245205,  0.990592 ],
        [ -7.30382014,  6.66380066,   -0.99829487, -7.32946634,  0.98703558 ],
        [ -6.75693356,  6.11815278,   -0.99706518, -6.7907593,   0.98285381 ],
        [ -6.27107344,  5.63412334,   -0.99525652, -6.31435974,  0.9779981 ],
        [ -5.83201892,  5.19767299,   -0.99270136, -5.88614365,  0.97241778 ],
        [ -5.42879645,  4.79805498,   -0.98919855, -5.49528385,  0.96604081 ],
        [ -5.05267196,  4.42682535,   -0.98450352, -5.1332558,   0.9587692 ],
        [ -4.69591917,  4.07664834,   -0.97830548, -4.79264998,  0.95046025 ],
        [ -4.35092377,  3.74047025,   -0.9701851,  -4.46634895,  0.9408963 ],
        [ -4.00839182,  3.40989748,   -0.95950522, -4.14592793,  0.92970202 ],
        [ -3.65298273,  3.07133529,   -0.94507651, -3.81785851,  0.9160964 ],
        [ -3.2617416,   2.70541667,   -0.92464052, -3.46277884,  0.89859882 ],
        [ -2.91002769,  2.38408319,   -0.9019132,  -3.1498366,   0.88056989 ],
        [ -2.58473991,  2.09461111,   -0.87733764, -2.86636424,  0.86205145 ],
        [ -2.28254156,  1.8332741,    -0.85185747, -2.60863366,  0.84344921 ],
        [ -1.99745286,  1.59406447,   -0.82605265, -2.37079919,  0.82490742 ],
        [ -1.71773896,  1.36667354,   -0.79971884, -2.14268615,  0.80605299 ],
        [ -1.43520462,  1.14452966,   -0.7727928,  -1.91768532,  0.78664439 ],
        [ -1.1437302,   0.9232927,    -0.74537858, -1.6913275,   0.76656467 ],
        [ -0.84664221,  0.70588011,   -0.71846658, -1.46660131,  0.74636561 ],
        [ -0.54041311,  0.48991048,   -0.69235534, -1.24115884,  0.72613491 ],
        [ -0.22123999,  0.2729905,    -0.66729712, -1.01264381,  0.70595787 ],
        [ 0.11499845,   0.05269475,   -0.6435139,  -0.77867875,  0.68592773 ],
        [ 0.47321343,   -0.17373041,  -0.62118424, -0.53656383,  0.66613196 ],
        [ 0.85920629,   -0.40939656,  -0.60046444, -0.28326103,  0.64666655 ],
        [ 1.27681726,   -0.65610105,  -0.58162474, -0.0172277,   0.62777585 ],
        [ 1.69095579,   -0.89367272,  -0.56617919, 0.23924942,   0.61117135 ],
        [ 2.11129741,   -1.12886858,  -0.55332639, 0.492963,     0.59632952 ],
        [ 2.5438915,    -1.36582983,  -0.54259365, 0.74797681,   0.58297627 ],
        [ 2.99347481,   -1.60768449,  -0.53365117, 1.00730178,   0.57094171 ],
        [ 3.46423937,   -1.8570952,   -0.52624595, 1.27346514,   0.56011071 ],
        [ 3.96072428,   -2.11679377,  -0.52016501, 1.54907159,   0.55039031 ],
        [ 4.48662276,   -2.38899207,  -0.51523694, 1.83617492,   0.54172591 ],
        [ 5.04836821,   -2.67726026,  -0.51129007, 2.13825821,   0.5340375 ],
        [ 5.65317286,   -2.98549919,  -0.50817929, 2.4591274,    0.52726605 ],
        [ 6.30543724,   -3.31614102,  -0.50578921, 2.80105756,   0.52139761 ],
        [ 7.01807119,   -3.67589994,  -0.50398736, 3.17075062,   0.51635438 ],
        [ 7.80467699,   -4.0717828,   -0.50266718, 3.57516273,   0.51209219 ],
        [ 8.68181627,   -4.51224958,  -0.50173312, 4.02271354,   0.50857387 ],
        [ 9.67724551,   -5.01134309,  -0.50109395, 4.52746779,   0.50574488 ],
        [ 10.83007401,  -5.5887504,   -0.50066965, 5.10915199,   0.50355856 ],
        [ 12.24501854,  -6.29694707,  -0.50038347, 5.82039808,   0.50192868 ],
        [ 14.09849459,  -7.2241919,   -0.50018548, 6.74955952,   0.50082944 ],
        [ 16.50506588,  -8.4277581,   -0.50006569, 7.95403569,   0.50026241 ],
        [ 19.33335373,  -9.84200507,  -0.50001718, 9.36858084,   0.50006512 ],
        [ 23.28159166,  -11.81615394, -0.50000244, 11.3428123,   0.5000091 ]
    ],

    energy_table => [
        [ -12.64178695, 0.0621806,  0.01434181, 0.0621806,  0.01434181,    11.50786,    -1.0016779 ],
        [ -11.55327663, 0.07992053, 0.01842119, 0.07992053, 0.01842119,    10.417073,   -1.0028252 ],
        [ -10.65825935, 0.09822151, 0.02261397, 0.09822151, 0.02261397,    9.5189874,   -1.0044378 ],
        [ -9.59080036,  0.12554896, 0.02882621, 0.12554896, 0.02882621,    8.4455065,   -1.0076091 ],
        [ -8.68066644,  0.15465917, 0.03534276, 0.15465917, 0.03534276,    7.5269203,   -1.0119463 ],
        [ -7.93379227,  0.18334407, 0.04161079, 0.18334407, 0.04161079,    6.7694938,   -1.0173547 ],
        [ -7.30382014,  0.21139452, 0.04753445, 0.21139452, 0.04753445,    6.1268745,   -1.0238979 ],
        [ -6.75693356,  0.23889252, 0.05308113, 0.23889252, 0.05308113,    5.5651039,   -1.0316522 ],
        [ -6.27107344,  0.26592497, 0.05821636, 0.26592497, 0.05821636,    5.0619503,   -1.040698 ],
        [ -5.83201892,  0.29251307, 0.06289092, 0.29251307, 0.06289092,    4.6030019,   -1.0511352 ],
        [ -5.42879645,  0.31871977, 0.06706157, 0.31871977, 0.06706157,    4.1770037,   -1.0631222 ],
        [ -5.05267196,  0.34463556, 0.07068526, 0.34463556, 0.07068526,    3.7748086,   -1.0766688 ],
        [ -4.69591917,  0.37040802, 0.07371811, 0.37040802, 0.07371811,    3.3882152,   -1.0920653 ],
        [ -4.35092377,  0.39627057, 0.07610847, 0.39627057, 0.07610847,    3.0086484,   -1.1090743 ],
        [ -4.00839182,  0.42264953, 0.07778548, 0.42264953, 0.07778548,    2.6257399,   -1.1277479 ],
        [ -3.65298273,  0.45047344, 0.07862197, 0.45047344, 0.07862197,    2.2212879,   -1.1478991 ],
        [ -3.2617416,   0.48121673, 0.07830874, 0.48121673, 0.07830874,    1.7679149,   -1.1652912 ],
        [ -2.91002769,  0.50853353, 0.07683612, 0.50853353, 0.07683612,    1.3560156,   -1.1693512 ],
        [ -2.58473991,  0.53316841, 0.07447301, 0.53316841, 0.07447301,    0.97617213,  -1.1513312 ],
        [ -2.28254156,  0.55524003, 0.07148083, 0.55524003, 0.07148083,    0.63284152,  -1.102667 ],
        [ -1.99745286,  0.57514234, 0.06805267, 0.57514234, 0.06805267,    0.32747808,  -1.0411724 ],
        [ -1.71773896,  0.59365264, 0.06423426, 0.59365264, 0.06423426,    0.044465772, -0.97746525 ],
        [ -1.43520462,  0.61121601, 0.06004918, 0.61121601, 0.06004918,    -0.22190542, -0.91820226 ],
        [ -1.1437302,   0.62806296, 0.05552558, 0.62806296, 0.05552558,    -0.48214374, -0.88071458 ],
        [ -0.84664221,  0.643863,   0.05083859, 0.643863,   0.05083859,    -0.74012354, -0.86426563 ],
        [ -0.54041311,  0.65869637, 0.04605718, 0.65869637, 0.04605718,    -1.0034939,  -0.86162416 ],
        [ -0.22123999,  0.67262241, 0.04124307, 0.67262241, 0.04124307,    -1.279027,   -0.86746598 ],
        [ 0.11499845,   0.68567571, 0.0364549,  0.68567571, 0.0364549,     -1.5721881,  -0.87824829 ],
        [ 0.47321343,   0.69787822, 0.03174642, 0.69787822, 0.03174642,    -1.8892184,  -0.89185076 ],
        [ 0.85920629,   0.70923208, 0.0271709,  0.70923208, 0.0271709,     -2.2363042,  -0.90619194 ],
        [ 1.27681726,   0.71964789, 0.02281371, 0.71964789, 0.02281371,    -2.6178977,  -0.92057012 ],
        [ 1.69095579,   0.72830317, 0.01907947, 0.72830317, 0.01907947,    -3.0019408,  -0.93322066 ],
        [ 2.11129741,   0.73562402, 0.01584158, 0.73562402, 0.01584158,    -3.3967255,  -0.94431172 ],
        [ 2.5438915,    0.74185103, 0.01303029, 0.74185103, 0.01303029,    -3.8075051,  -0.95400345 ],
        [ 2.99347481,   0.7471452,  0.01059838, 0.7471452,  0.01059838,    -4.2384793,  -0.96242993 ],
        [ 3.46423937,   0.75162599, 0.00850947, 0.75162599, 0.00850947,    -4.6934412,  -0.96970261 ],
        [ 3.96072428,   0.75539293, 0.00673094, 0.75539293, 0.00673094,    -5.1765964,  -0.97592545 ],
        [ 4.48662276,   0.75852397, 0.00523641, 0.75852397, 0.00523641,    -5.6913791,  -0.98118223 ],
        [ 5.04836821,   0.76110147, 0.0039945,  0.76110147, 0.0039945,     -6.2439453,  -0.9855764 ],
        [ 5.65317286,   0.76319516, 0.00297752, 0.76319516, 0.00297752,    -6.8412748,  -0.98920513 ],
        [ 6.30543724,   0.76485823, 0.00216416, 0.76485823, 0.00216416,    -7.4875988,  -0.99211911 ],
        [ 7.01807119,   0.76615933, 0.00152415, 0.76615933, 0.00152415,    -8.1955714,  -0.99441356 ],
        [ 7.80467699,   0.76715277, 0.00103314, 0.76715277, 0.00103314,    -8.9786116,  -0.996188 ],
        [ 8.68181627,   0.76788761, 0.00066854, 0.76788986, 0.00066854,    -9.8531124,  -0.99751076 ],
        [ 9.67724551,   0.76841246, 0.00040735, 0.76841369, 0.00040735,    -10.84665,   -0.99846732 ],
        [ 10.83007401,  0.76876963, 0.00022921, 0.76877075, 0.00022921,    -11.998202,  -0.9991255 ],
        [ 12.24501854,  0.76900218, 0.00011305, 0.76900324, 0.00011305,    -13.412313,  -0.99955718 ],
        [ 14.09849459,  0.76913879, 4.476e-05,  0.76913962, 4.476e-05,     -15.265324,  -0.99993854 ],
        [ 16.50506588,  0.76920143, 1.345e-05,  0.76920187, 1.345e-05,     -17.672339,  -1.0006453 ],
        [ 19.33335373,  0.76922177, 3.27e-06,   0.76922183, 3.27e-06,      -20.503984,  -1 ],
        [ 23.28159166,  0.7692274,  4.5e-07,    0.7692277,  1.5343151e-06, -24.452222,  -1 ]
    ],

    impulse_table => [
        [
            -14.456615, 13.814889,  -14.457324,    0.32237051, 0,          -1.2005793,
            0,          -8.4968229, -0.0021902096, -2.0252742, 0.15875046, 0.94681824,
            -4.3986941
        ],
        [
            -13.815511, 13.173785,  -13.816488,    0.32221625, 0,          -1.2005788,
            0,          -8.1762726, -0.0021902096, -2.0252737, 0.15875032, 0.93833998,
            -4.2507681
        ],
        [
            -13.122363, 12.48064,   -13.123746,    0.32198326, 0,          -1.2005778,
            0,          -7.8297023, -0.0021902096, -2.0252727, 0.15874994, 0.92764775,
            -4.0908535
        ],
        [
            -12.716898, 12.075177,  -12.718592,    0.32180448, 0,          -1.2005768,
            0,          -7.6269728, -0.0021902096, -2.0252717, 0.15874955, 0.92055419,
            -3.9973234
        ],
        [
            -12.206073, 11.564357, -12.208259,    0.32152098, 0,          -1.2005748,
            0,          -7.371566, -0.0021902096, -2.0252697, 0.15874882, 0.91062093,
            -3.8795109
        ],
        [
            -11.512925, 10.871222,  -11.51602,     0.32100002, 0,          -1.2005698,
            0,          -7.0250071, -0.0021902096, -2.0252647, 0.15874702, 0.89513325,
            -3.719706
        ],
        [
            -10.819778, 10.1781,    -10.824157,    0.32026332, 0,          -1.2005598,
            0,          -6.6784626, -0.0021902096, -2.0252547, 0.15874348, 0.87697461,
            -3.5600051
        ],
        [
            -10.414313, 9.77266,    -10.419678,    0.31969809, 0,          -1.2005498,
            0,          -6.4757593, -0.0021902096, -2.0252447, 0.15873993, 0.86493777,
            -3.4666599
        ],
        [
            -9.9034875, 9.2618853,  -9.9104183,    0.31880187, 0,          -1.2005298,
            0,          -6.2204048, -0.0021902096, -2.0252247, 0.15873284, 0.8481,
            -3.3491726
        ],
        [
            -9.2103404, 8.5688653,  -9.2201556,    0.31715541, 0,         -1.2004798,
            0,          -5.8739771, -0.0021902096, -2.0251747, 0.1587151, 0.82190059,
            -3.1900519
        ],
        [
            -8.5171932, 7.8759728,  -8.5311016,    0.31482847, 0,          -1.2003798,
            0,          -5.5276955, -0.0021902096, -2.0250747, 0.15867962, 0.79129573,
            -3.0314732
        ],
        [
            -8.1117281, 7.4707619, -8.1287868,    0.31304453, 0,          -1.2002798,
            0,          -5.325255, -0.0021902096, -2.0249747, 0.15864416, 0.77109647,
            -2.9390886
        ],
        [
            -7.6009024, 6.9604439,  -7.622974,     0.31021927, 0,         -1.2000798,
            0,          -5.0704263, -0.0021902096, -2.0247747, 0.1585733, 0.74298969,
            -2.8232689
        ],
        [
            -6.9077553, 6.268563,   -6.939092,     0.30504289, 0,          -1.1995798,
            0,          -4.7253136, -0.0021902096, -2.0242747, 0.15839652, 0.69969485,
            -2.6675765
        ],
        [
            -6.2146081, 5.5779333,  -6.2591549,    0.2977697,  0,          -1.1985798,
            0,          -4.3816641, -0.0021902096, -2.0232747, 0.15804463, 0.65000299,
            -2.5144342
        ],
        [
            -5.809143, 5.1749659,  -5.8639024,    0.29223857, 0,          -1.1975798,
            0,         -4.1818589, -0.0021902096, -2.0222747, 0.15769495, 0.61786004,
            -2.4265566
        ],
        [
            -5.2983174, 4.6690774,  -5.3693881,    0.28358209, 0,          -1.1955798,
            0,          -3.9323103, -0.0021902096, -2.0202747, 0.15700212, 0.57417969,
            -2.3182981
        ],
        [
            -4.6051702, 3.9879534,  -4.7065029,    0.26815624, 0,          -1.1905798,
            0,          -3.6004558, -0.0021902096, -2.0152747, 0.15530692, 0.50973266,
            -2.1773022
        ],
        [
            -3.912023, 3.3175998,  -4.0565001,    0.24779225, 0,          -1.1805801,
            0,         -3.2835938, -0.0021902094, -2.0052747, 0.15206441, 0.44079037,
            -2.0456172
        ],
        [
            -3.5065579, 2.933461,   -3.6841709,    0.23366361, 0,          -1.1705822,
            0,          -3.1109727, -0.0021902084, -1.9952747, 0.14900154, 0.3994124,
            -1.9741304
        ],
        [
            -2.9957323, 2.4616367,  -3.2255031,    0.21459677, 0,          -1.1506164,
            0,          -2.9170879, -0.0021901917, -1.9752747, 0.14334543, 0.34749423,
            -1.891026
        ],
        [
            -2.3025851, 1.850366,   -2.6255521,    0.19249551, 0,          -1.1019597,
            0,          -2.7312166, -0.0021895362, -1.9267345, 0.13136873, 0.28005376,
            -1.7922226
        ],
        [
            -1.6094379, 1.2806222,  -2.055791,     0.18827712, 0,          -1.0288535,
            0,          -2.6998317, -0.0021762734, -1.8502291, 0.11359966, 0.21931819,
            -1.7104064
        ],
        [
            -1.2039728, 0.96836501, -1.7376321,    0.19533738, 0,          -0.99397617,
            0,          -2.7550047, -0.0021435879, -1.8116271, 0.10088983, 0.1879035,
            -1.6701843
        ],
        [
            -0.69314717, 0.59662748, -1.3528256,    0.20906995, 0,           -0.97584779,
            0,           -2.8729451, -0.0020529193, -1.7894643, 0.083789678, 0.15307373,
            -1.6267732
        ],
        [
            1.4784578e-08, 0.12714475, -0.85794307,   0.22960097, 0,          -0.99522523,
            0,             -3.0844943, -0.0018348939, -1.8209203, 0.06200427, 0.11418913,
            -1.5789423
        ],
        [
            0.26236428, -0.041426537, -0.67821459,   0.23701023, 0,           -1.010578,
            0,          -3.1749576,   -0.0017354392, -1.8435433, 0.054877274, 0.10182757,
            -1.5636523
        ],
        [
            0.40546512, -0.13151178, -0.58181536,   0.24088942, 0,           -1.0197413,
            0,          -3.2261842,  -0.0016806087, -1.8568756, 0.051295278, 0.095589987,
            -1.5558871
        ],
        [
            0.53062827, -0.20930069, -0.49840484,   0.24417736, 0,           -1.0279454,
            0,          -3.2719783,  -0.0016332493, -1.8689639, 0.048337296, 0.090413413,
            -1.5494096
        ],
        [
            0.6931472, -0.3089865, -0.3913164,    0.24829219, 0,           -1.038601,
            0,         -3.3327188, -0.0015735212, -1.8858076, 0.044732889, 0.08406534,
            -1.5414184
        ],
        [
            0.99325179, -0.48945112, -0.19700345,  0.25541303, 0,           -1.0574284,
            0,          -3.4483377,  -0.001471518, -1.9079735, 0.038737797, 0.073389605,
            -1.5278409
        ],
        [
            1.0986123, -0.55178408, -0.12978855,   0.25776364, 0,           -1.0635738,
            0,         -3.4899041,  -0.0014389264, -1.9151878, 0.036822608, 0.069944932,
            -1.5234186
        ],
        [
            1.2089604, -0.61653857, -0.059925312,  0.26014228, 0,           -1.0696818,
            0,         -3.5339416,  -0.0014068075, -1.9219142, 0.034915229, 0.066496851,
            -1.5189703
        ],
        [
            1.3862944, -0.71953432, 0.0512469,    0.26378804, 0,           -1.0787158,
            0,         -3.6057329,  -0.001359691, -1.9323005, 0.032049565, 0.06128277,
            -1.5122008
        ],
        [
            1.6094379, -0.84740551, 0.18930197,    0.26807186, 0,          -1.088629,
            0,         -3.6977374,  -0.0013083661, -1.9444655, 0.02876815, 0.055261287,
            -1.5043168
        ],
        [
            1.7917595, -0.95057705, 0.30066852,    0.27132831, 0,           -1.0955161,
            0,         -3.7741861,  -0.0012728174, -1.9469199, 0.026333135, 0.050757019,
            -1.498372
        ],
        [
            1.9459102, -1.0369686, 0.39387805,    0.27391679, 0,          -1.1005076,
            0,         -3.8396538, -0.0012469621, -1.9527924, 0.02443296, 0.047220466,
            -1.4936758
        ],
        [
            2.0252604, -1.0811588, 0.44153325,    0.27519212, 0,           -1.1028178,
            0,         -3.873635,  -0.0012350468, -1.9523316, 0.023507982, 0.045491999,
            -1.4913715
        ],
        [
            2.3025851, -1.2342249, 0.60644123,    0.27935714, 0,           -1.1095106,
            0,         -3.9938026, -0.0012001436, -1.962889,  0.020537606, 0.039910436,
            -1.4838901
        ],
        [
            2.7080502, -1.4546102, 0.84329675,    0.28468355, 0,           -1.1161485,
            0,         -4.1730633, -0.0011646772, -1.9677065, 0.016848083, 0.032911179,
            -1.4744245
        ],
        [
            2.9957323,    -1.6088892, 1.0085906,   0.28796387, 0, -1.1192217, 0, -4.3025189,
            -0.001147913, -1.9728841, 0.014634688, 0.0286766,  -1.4686555
        ],
        [
            3.4011974, -1.8238918, 1.2381122,     0.29196421, 0,           -1.1218364,
            0,         -4.4877485, -0.0011323681, -1.9736176, 0.011994464, 0.023590202,
            -1.4616897
        ],
        [
            3.528874, -1.8910803, 1.3096237,     0.29308624, 0,          -1.1223018,
            0,        -4.5466804, -0.0011290107, -1.9882464, 0.01126494, 0.020497889,
            -1.5291672
        ],
        [
            3.9120415, -2.0914579, 1.5222556,     0.29610302, 0,           -1.1233831,
            0,         -4.7250674, -0.0011215791, -2.0033272, 0.009328887, 0.017247452,
            -1.5124418
        ],
        [
            4.6052157, -2.4500398, 1.9003156,     0.30040846, 0,            -1.1240919,
            0,         -5.0527739, -0.0011149241, -2.0517146, 0.0066261972, 0.012528748,
            -1.4887355
        ],
        [
            5.2984017, -2.8049202, 2.2714089,     0.303544,   0,            -1.1241289,
            0,         -5.3855489, -0.0011123936, -2.1495488, 0.0047016179, 0.0090344356,
            -1.4715993
        ],
        [
            6.2147824, -3.2702761, 2.7537572,     0.30639453, 0,            -1.1239125,
            0,         -5.8311151, -0.0011113037, -2.4443348, 0.0029832663, 0.0058163982,
            -1.4561195
        ],
        [
            6.9078191, -3.6203216, 3.1137833,     0.30785664, 0,            -1.1237386,
            0,         -6.1711525, -0.0011110558, -2.9361292, 0.0021132595, 0.0041507535,
            -1.4482172
        ],
        [
            7.601071, -3.9694082, 3.4707985,     0.30890177, 0,            -1.1235971,
            0,        -6.5131731, -0.0011109662, -3.9203867, 0.0014960502, 0.002953973,
            -1.4425845
        ],
        [
            8.5173829, -4.4297365, 3.9390404,     0.30983675, 0,             -1.1234623,
            0,         -6.9672824, -0.0011109291, -6.8733675, 0.00094714967, 0.0018789603,
            -1.4375569
        ],
        [
            9.210481, -4.7773931, 4.2911339,     0.3103107,  0,             -1.1232956,
            0,        -7.3118491, -0.0010857681, -11.373657, 0.00067005767, 0.0013324237,
            -1.4350123
        ],
        [
            10.819857, -5.5836349, 5.1040069,     0.31094571, 0,             -1.1217125,
            0,         -8.1140825, -0.0011109704, -51.162665, 0.00029979129, 0.00059803879,
            -1.431606
        ],
        [
            11.513263, -5.9307454, 5.4528575,     0.31109618, 0,             -1.1216711,
            0,         -8.4602798, -0.0011109165, -100.39882, 0.00021196129, 0.00042315034,
            -1.4307978
        ],
        [
            13.12244, -6.735944, 6.260512,      0.3112929,  0,             -1.1216813,
            0,        -9.2643,   -0.0011109163, -494.07733, 9.4796764e-05, 0.00018943894,
            -1.4297277
        ],
        [
            13.815583, -7.0826804, 6.6078534,     0.31133437, 0,             -1.1216974,
            0,         -9.610766,  -0.0011109163, -986.20254, 6.7028333e-05, 0.00013397939,
            -1.4294845
        ],
        [
            16.118205, -8.2343002, 7.7604948,     0.31129821, 0,             -1.1220018,
            0,         -10.76193,  -0.0011109163, -9844.8716, 2.1193399e-05, 4.2379406e-05,
            -1.429297
        ],
        [
            18.420681, -9.3856493, 8.9121711,     0.31017142, 0,             -1.125058,
            0,         -11.913044, -0.0011109163, -98420.42,  6.7009547e-06, 1.3404538e-05,
            -1.4315204
        ]
    ],

    tail_shock_table => [
        [ 8.4514717, -6.758,     -0.00073650519, -0.00073650519 ],
        [ 8.5070579, -6.9838948, -0.00082482236, -0.00062513214 ],
        [ 9.2031503, -10.483669, -0.00077888887, -0.00027967097 ],
        [ 10.816558, -26.040725, -0.00043198465, -6.8002902e-05 ]
    ],

    blast_info => [
        1.2005742, 0,            2.0254571, -0.0021902187, 0,          0,          0.3223203,  -0.014630836,
        6.0195959, 0.0004737795, 2.8316302, -0.0019055604, 0.40896362, 0.15875046, 0.23946018, 0,
        4675.2,    -6.758,       5541.1377, -7.4049291
    ],

};
1;
