package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='P', gamma=7
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'P7'} = {

    table_name  => 'P7',
    symmetry    => 0,
    gamma       => 7,
    data_source => 'P4000_G7/AF_restart_moc2/',

    shock_table_info => [ 0, 7, 2.03e-06, 4000, 5.1e-07, 500, 3.81e-07, 1.15e-06, 5e-07 ],

    shock_table => [
        [ -9.7873635,  12.00011469, -0.99998925, -9.78955181, 0.99890466 ],
        [ -8.84170448, 11.05447057, -0.99997232, -8.84521799, 0.99824021 ],
        [ -8.20260649, 10.41539396, -0.9999541,  -8.207446,   0.99757452 ],
        [ -6.98621535, 9.19911176,  -0.99984526, -6.99512354, 0.9955269 ],
        [ -6.05108108, 8.26421686,  -0.99960615, -6.06533528, 0.9928255 ],
        [ -5.29352378, 7.50710553,  -0.99916104, -5.31440495, 0.98946104 ],
        [ -4.65173862, 6.86607451,  -0.99840949, -4.68062121, 0.98537814 ],
        [ -4.09639068, 6.3119056,   -0.99723782, -4.13466671, 0.98056043 ],
        [ -3.6050869,  5.8223499,   -0.99550783, -3.65423033, 0.97496041 ],
        [ -3.16215213, 5.38190959,  -0.99305406, -3.22375912, 0.96851288 ],
        [ -2.75706164, 4.98027117,  -0.98968672, -2.83286212, 0.96115182 ],
        [ -2.38041033, 4.60830582,  -0.98516795, -2.47237104, 0.95276598 ],
        [ -2.02491117, 4.25908231,  -0.97920851, -2.13530959, 0.94321857 ],
        [ -1.68386484, 3.926387,    -0.97143791, -1.81543521, 0.93231793 ],
        [ -1.3482262,  3.60196181,  -0.96128124, -1.50457617, 0.91969445 ],
        [ -1.00713784, 3.27627686,  -0.94782467, -1.19336662, 0.9047369 ],
        [ -0.63179024, 2.9239492,   -0.92871735, -0.85727226, 0.88563905 ],
        [ -0.27762748, 2.59889724,  -0.90613718, -0.54718177, 0.86507807 ],
        [ 0.04702538,  2.30861378,  -0.88152848, -0.26966796, 0.84422927 ],
        [ 0.35193646,  2.04374986,  -0.85533768, -0.01542833, 0.82319853 ],
        [ 0.63386022,  1.80628158,  -0.82900506, 0.21379841,  0.80284957 ],
        [ 0.90655301,  1.58384567,  -0.80224395, 0.42998915,  0.78270502 ],
        [ 1.17900333,  1.36898893,  -0.7749381,  0.64048211,  0.76249606 ],
        [ 1.45907176,  1.15587608,  -0.7470067,  0.85115376,  0.74201291 ],
        [ 1.7464834,   0.94520315,  -0.71919068, 1.0614718,   0.72165848 ],
        [ 2.04175241,  0.73688931,  -0.69211538, 1.27158934,  0.7017628 ],
        [ 2.34849082,  0.52864477,  -0.66606476, 1.48384467,  0.68242959 ],
        [ 2.67072699,  0.3180821,   -0.64127324, 1.70069225,  0.66374925 ],
        [ 3.01319881,  0.10254972,  -0.61792812, 1.92487765,  0.6457973 ],
        [ 3.38134075,  -0.12083178, -0.59619807, 2.15939982,  0.62865182 ],
        [ 3.78187939,  -0.35551068, -0.57622294, 2.40786125,  0.61238138 ],
        [ 4.18432259,  -0.58395042, -0.55957389, 2.65140513,  0.59830229 ],
        [ 4.5915569,   -0.80891586, -0.54573598, 2.89250396,  0.58609556 ],
        [ 5.0100367,   -1.03479133, -0.53417958, 3.13547465,  0.5753991 ],
        [ 5.44454401,  -1.26472557, -0.52455287, 3.383386,    0.56597813 ],
        [ 5.89889924,  -1.50117849, -0.51659894, 3.63859862,  0.55766661 ],
        [ 6.37687061,  -1.74647818, -0.51010657, 3.90334182,  0.55032951 ],
        [ 6.89063129,  -2.00712708, -0.50482262, 4.18433862,  0.5437567 ],
        [ 7.42112192,  -2.27380652, -0.50079645, 4.47125609,  0.53812424 ],
        [ 7.99372891,  -2.5596243,  -0.49769324, 4.77790395,  0.53309674 ],
        [ 8.61366123,  -2.86740514, -0.49541672, 5.10695634,  0.5286279 ],
        [ 9.28144812,  -3.19767932, -0.49387597, 5.45861043,  0.52469923 ],
        [ 10.02035348, -3.56222057, -0.49294111, 5.84496537,  0.52117625 ],
        [ 10.82328667, -3.95782352, -0.49253867, 6.2621525,   0.51809362 ],
        [ 11.71884227, -4.3989026,  -0.49256584, 6.72485475,  0.51533987 ],
        [ 12.73808682, -4.90112048, -0.49295002, 7.24879364,  0.51284999 ],
        [ 13.91663575, -5.48246848, -0.49362409, 7.85182302,  0.51058445 ],
        [ 15.31029784, -6.17103668, -0.49452291, 8.56189101,  0.50850123 ],
        [ 17.09589147, -7.05505981, -0.49563014, 9.46797232,  0.5064785 ],
        [ 19.88600032, -8.43999446, -0.49704542, 10.87779978, 0.50426695 ]
    ],

    energy_table => [
        [ -9.7873635,  3.612e-05,  2.962e-05,  3.612e-05,  2.962e-05,  11.503514,   -0.99875141 ],
        [ -8.84170448, 7.819e-05,  6.357e-05,  7.819e-05,  6.357e-05,  10.558716,   -0.999497 ],
        [ -8.20260649, 0.00013123, 0.00010597, 0.00013123, 0.00010597, 9.9197637,   -1.0001316 ],
        [ -6.98621535, 0.0003477,  0.00027604, 0.0003477,  0.00027604, 8.7023814,   -1.0016587 ],
        [ -6.05108108, 0.00072601, 0.00056634, 0.00072601, 0.00056634, 7.7650894,   -1.0030879 ],
        [ -5.29352378, 0.00130416, 0.00099874, 0.00130416, 0.00099874, 7.0047139,   -1.0044409 ],
        [ -4.65173862, 0.00212209, 0.00159336, 0.00212209, 0.00159336, 6.3596867,   -1.0057275 ],
        [ -4.09639068, 0.0032064,  0.00235704, 0.0032064,  0.00235704, 5.8008332,   -1.0069477 ],
        [ -3.6050869,  0.00458291, 0.00329264, 0.00458291, 0.00329264, 5.3058398,   -1.0081131 ],
        [ -3.16215213, 0.0062761,  0.00439802, 0.0062761,  0.00439802, 4.8590706,   -1.0092354 ],
        [ -2.75706164, 0.00830526, 0.00566315, 0.00830526, 0.00566315, 4.4500246,   -1.0103234 ],
        [ -2.38041033, 0.01069729, 0.00707749, 0.01069729, 0.00707749, 4.0692892,   -1.0113901 ],
        [ -2.02491117, 0.01348249, 0.00862547, 0.01348249, 0.00862547, 3.7095574,   -1.0124472 ],
        [ -1.68386484, 0.01670331, 0.01028912, 0.01670331, 0.01028912, 3.364089,    -1.0135019 ],
        [ -1.3482262,  0.02045101, 0.01206046, 0.02045101, 0.01206046, 3.0237416,   -1.0143386 ],
        [ -1.00713784, 0.02488324, 0.01393398, 0.02488324, 0.01393398, 2.6776552,   -1.0146297 ],
        [ -0.63179024, 0.03049883, 0.01597299, 0.03049883, 0.01597299, 2.2968261,   -1.0140925 ],
        [ -0.27762748, 0.03647597, 0.01774427, 0.03647597, 0.01774427, 1.9378432,   -1.0126079 ],
        [ 0.04702538,  0.04246795, 0.01912063, 0.04246795, 0.01912063, 1.6093954,   -1.0103386 ],
        [ 0.35193646,  0.04845749, 0.02011218, 0.04845749, 0.02011218, 1.3017187,   -1.0075016 ],
        [ 0.63386022,  0.05422058, 0.02071965, 0.05422058, 0.02071965, 1.018089,    -1.0044406 ],
        [ 0.90655301,  0.05991618, 0.0210029,  0.05991618, 0.0210029,  0.74461011,  -1.001294 ],
        [ 1.17900333,  0.06564299, 0.02098789, 0.06564299, 0.02098789, 0.47223968,  -0.99819255 ],
        [ 1.45907176,  0.07148476, 0.02068249, 0.07148476, 0.02068249, 0.1931134,   -0.99525429 ],
        [ 1.7464834,   0.07735161, 0.02010158, 0.07735161, 0.02010158, -0.09252674, -0.99266658 ],
        [ 2.04175241,  0.08316993, 0.01927456, 0.08316993, 0.01927456, -0.38527557, -0.99055918 ],
        [ 2.34849082,  0.08892598, 0.01822988, 0.08892598, 0.01822988, -0.68882954, -0.98899373 ],
        [ 2.67072699,  0.09460441, 0.0169962,  0.09460441, 0.0169962,  -1.0073076,  -0.98799911 ],
        [ 3.01319881,  0.10018791, 0.015602,   0.10018791, 0.015602,   -1.3455463,  -0.98758468 ],
        [ 3.38134075,  0.10565086, 0.01407714, 0.10565086, 0.01407714, -1.709096,   -0.98775701 ],
        [ 3.78187939,  0.11096185, 0.01245276, 0.11096185, 0.01245276, -2.1048316,  -0.98856251 ],
        [ 4.18432259,  0.11565871, 0.01090636, 0.11565871, 0.01090636, -2.5028958,  -0.99005814 ],
        [ 4.5915569,   0.11980134, 0.00946056, 0.11980134, 0.00946056, -2.9064684,  -0.99290193 ],
        [ 5.0100367,   0.12347399, 0.00811639, 0.12347399, 0.00811639, -3.3227921,  -0.99731268 ],
        [ 5.44454401,  0.12672577, 0.00687789, 0.12672577, 0.00687789, -3.7572438,  -1 ],
        [ 5.89889924,  0.12958824, 0.00574999, 0.12958824, 0.00574999, -4.2115991,  -1 ],
        [ 6.37687061,  0.13208748, 0.00473589, 0.13208748, 0.00473589, -4.6895704,  -1 ],
        [ 6.89063129,  0.13427882, 0.00382353, 0.13427882, 0.00382353, -5.2033311,  -1 ],
        [ 7.42112192,  0.13609503, 0.00305017, 0.13609503, 0.00305017, -5.7338217,  -1 ],
        [ 7.99372891,  0.13764189, 0.00237846, 0.13764189, 0.00237846, -6.3064287,  -1 ],
        [ 8.61366123,  0.13893208, 0.00180831, 0.13893208, 0.00180831, -6.9263611,  -1 ],
        [ 9.28144812,  0.1399758,  0.00133984, 0.1399758,  0.00133984, -7.5941479,  -1 ],
        [ 10.02035348, 0.14081679, 0.00095706, 0.14081679, 0.00095706, -8.3330533,  -1 ],
        [ 10.82328667, 0.14145928, 0.00066094, 0.14145928, 0.00066094, -9.1359865,  -1 ],
        [ 11.71884227, 0.14194333, 0.00043532, 0.14194333, 0.00043532, -10.031542,  -1 ],
        [ 12.73808682, 0.14229583, 0.00026932, 0.14229583, 0.00026932, -11.050787,  -1 ],
        [ 13.91663575, 0.14253891, 0.00015375, 0.14253891, 0.00015375, -12.229336,  -1 ],
        [ 15.31029784, 0.14269522, 7.876e-05,  0.14269522, 7.876e-05,  -13.622998,  -1 ],
        [ 17.09589147, 0.14278942, 3.317e-05,  0.14278942, 3.317e-05,  -15.408591,  -1 ],
        [ 19.88600032, 0.14283999, 8.48e-06,   0.14283999, 8.48e-06,   -18.1987,    -1 ]
    ],

    impulse_table => [
        [
            -11.602545, 13.815287,  -11.603427,    5.5462754,  0,          0,
            0,          -7.2524588, 0.00042238529, -454.78152, 0.34056091, 1.0006244,
            -7.5186465
        ],
        [
            -2.302585, 4.5316787,  -2.3982968,    5.4085675,  0,          0,
            0,         -2.6655293, 0.00042238535, -454.68153, 0.33621142, 0.92133475,
            -2.6936996
        ],
        [
            -1.6094379, 3.8541599,  -1.7461428,    5.361321,   0,          0,
            0,          -2.3731597, 0.00042238541, -454.58153, 0.33206612, 0.87776482,
            -2.2542417
        ],
        [
            -1.3862943, 3.6385807,  -1.5396164,    5.3444767,  0,          0,
            0,          -2.2866335, 0.00042238544, -454.53153, 0.33006454, 0.86017434,
            -2.1200728
        ],
        [
            -1.2039728, 3.4636664,  -1.3723389,    5.330349,   0,          0,
            0,          -2.2194228, 0.00042238552, -454.48153, 0.32810768, 0.84441037,
            -2.0134426
        ],
        [
            -0.91629069, 3.1903577,  -1.1113709,    5.3077986,  0,          0,
            0,           -2.1207054, 0.00042238573, -454.38153, 0.32432088, 0.81693623,
            -1.8510882
        ],
        [
            -0.69314714, 2.9810389,  -0.91171436,   5.2905195,  0,          0,
            0,           -2.0511531, 0.00042238594, -454.28153, 0.32069228, 0.79343558,
            -1.7304768
        ],
        [
            -0.51082558, 2.8120376,  -0.75054676,   5.2768651,  0,          0,
            0,           -1.9994206, 0.00042238626, -454.18153, 0.31721001, 0.77285047,
            -1.6355991
        ],
        [
            -0.3566749, 2.6707412,  -0.61575438,   5.2658498,  0,         0,
            0,          -1.9595537, 0.00042238664, -454.08153, 0.3138635, 0.75451711,
            -1.5580705
        ],
        [
            -0.22314351, 2.5496324,  -0.50014059,   5.2568323,  0,          0,
            0,           -1.9280601, 0.00042238702, -453.98153, 0.31064329, 0.73798532,
            -1.4929691
        ],
        [
            4.3821648e-08, 2.3501573,  -0.30944186,   5.2431459,  0,          0,
            0,             -1.8820794, 0.00042238811, -453.78153, 0.30454864, 0.70911981,
            -1.3885677
        ],
        [
            0.26236431, 2.1207233,  -0.089447145,  5.2298064,  0,          0,
            0,          -1.8392115, 0.00042239015, -453.48153, 0.29616657, 0.67345847,
            -1.2730384
        ],
        [
            0.53062829, 1.8923713,  0.13052926,    5.219863,   0,          0,
            0,          -1.8084132, 0.00042239384, -453.08153, 0.28617602, 0.63549229,
            -1.1631582
        ],
        [
            0.69314722,   1.7573023,  0.26126786, 5.2158604,  0, 0, 0, -1.7962794,
            0.0004223973, -452.78153, 0.27942206, 0.61196656, -1.1006794
        ],
        [
            0.99325182, 1.5146676, 0.49756942,    5.2126849,  0,          0,
            0,          -1.786761, 0.00042240763, -452.08153, 0.26561881, 0.5679585,
            -0.99337034
        ],
        [
            1.0986123, 1.4316119,  0.57894509,    5.21287,    0,          0,
            0,         -1.7873169, 0.00042241311, -451.78153, 0.26039161, 0.55244863,
            -0.95814471
        ],
        [
            1.2089604, 1.3458192,  0.66329108,    5.2137711,  0,          0,
            0,         -1.7900158, 0.00042242026, -451.43153, 0.25472477, 0.53622809,
            -0.92259075
        ],
        [
            1.3862944, 1.2105034,  0.79696019,    5.2166788,  0,          0,
            0,         -1.7987709, 0.00042243564, -450.78153, 0.24524771, 0.51029964,
            -0.86827065
        ],
        [
            1.609438, 1.0446607,  0.96191608,    5.2227203,  0,          0,
            0,        -1.8172085, 0.00042246495, -449.78153, 0.23278352, 0.47810428,
            -0.80470836
        ],
        [
            1.7917595, 0.91273774, 1.0940748,     5.2294239,  0,          0,
            0,         -1.8380752, 0.00042250096, -448.78153, 0.22226849, 0.45231662,
            -0.75658153
        ],
        [
            1.9459102, 0.80363353, 1.2040284,     5.2361753,  0,        0,
            0,         -1.8595441, 0.00042254368, -447.78153, 0.213228, 0.4309724,
            -0.71844696
        ],
        [
            2.0794416, 0.7108668,  1.2979919,     5.2427224,  0,          0,
            0,         -1.8808179, 0.00042259313, -446.78153, 0.20533727, 0.41287441,
            -0.68723032
        ],
        [
            2.3025851, 0.55930688, 1.4524531,     5.2548727,  0,          0,
            0,         -1.9215575, 0.00042271234, -444.78153, 0.19214091, 0.38353843,
            -0.6386518
        ],
        [
            2.7080502, 0.29419809, 1.7254272,     5.2797188,  0,          0,
            0,         -2.0105434, 0.00042312833, -439.78153, 0.16861822, 0.33348882,
            -0.56101354
        ],
        [
            2.9957323, 0.11335249, 1.9135903,     5.2985549,  0,          0,
            0,         -2.0839178, 0.00042371353, -434.78153, 0.15266628, 0.30074179,
            -0.51345612
        ],
        [
            3.4011974, -0.1326596, 2.1718742,     5.32541,    0,          0,
            0,         -2.1995485, 0.00042539573, -424.78153, 0.13167656, 0.25864541,
            -0.45578692
        ],
        [
            3.912023, -0.43012573, 2.4872449,     5.3576446,  0,          0,
            0,        -2.3621619,  0.00043085477, -404.78153, 0.10815353, 0.21231941,
            -0.39672442
        ],
        [
            4.2484953, -0.61978351, 2.6897338,     5.3766397,  0,          0,
            0,         -2.4777904,  0.00043924674, -384.78153, 0.09452726, 0.18571446,
            -0.36505935
        ],
        [
            4.6051702, -0.81634232, 2.9004801,     5.3934715,  0,           0,
            0,         -2.6065252,  0.00045791241, -354.78153, 0.081665473, 0.1606722,
            -0.33715389
        ],
        [
            5.2983174, -1.187806, 3.3004082,     5.4048688,  0,           0,
            0,         -2.871593, 0.00059630931, -254.78152, 0.060958509, 0.12036349,
            -0.29940933
        ]
    ],

    blast_info => [ 0, 0, 0, 0, 0, 0, 5.5421057, 0, 0, 0, 0, 0, 0.23096679, 0.34056093, 0, 0, 0, 0, 0, 0 ],

};
1;