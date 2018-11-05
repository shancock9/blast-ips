package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='C', gamma=2.2
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'C2x2'} = {

    table_name  => 'C2x2',
    symmetry    => 1,
    gamma       => 2.2,
    data_source => 'C4000_G2x2/moc2_from_moc_r2000/',

    shock_table_info => [ 1, 2.2, 1.04e-06, 4000, 1.05e-06, 70, 4.23e-07, 1.12e-07, 5e-07, 108.65, -1.5532 ],

    shock_table => [
        [ -6.29450485, 12.00033668,  -1.99998311, -6.29595896, 0.99854485 ],
        [ -5.71452454, 10.84039278,  -1.99994613, -5.71712304, 0.99739818 ],
        [ -5.44535993, 10.30208066,  -1.99992257, -5.44876237, 0.99659191 ],
        [ -5.09580423, 9.60301052,   -1.99984514, -5.10063372, 0.99515925 ],
        [ -4.86918644, 9.14982256,   -1.99973953, -4.87524789, 0.99392098 ],
        [ -4.38348092, 8.17862522,   -1.99931248, -4.39335036, 0.99008532 ],
        [ -3.99115447, 7.39438165,   -1.99849533, -4.00579705, 0.98526148 ],
        [ -3.66287218, 6.73851463,   -1.99710478, -3.68325531, 0.97943938 ],
        [ -3.37959995, 6.17307167,   -1.99491472, -3.40673309, 0.97256974 ],
        [ -3.12902158, 5.67356521,   -1.99164702, -3.16398675, 0.96457452 ],
        [ -2.90419524, 5.22627557,   -1.98699523, -2.94811517, 0.95541256 ],
        [ -2.69885086, 4.81887316,   -1.9805791,  -2.75296047, 0.94497669 ],
        [ -2.5082859,  4.44221881,   -1.97193926, -2.57397367, 0.9331281 ],
        [ -2.32869354, 4.08904444,   -1.96051087, -2.40756246, 0.9196823 ],
        [ -2.15658239, 3.75284094,   -1.94556619, -2.25055582, 0.90437358 ],
        [ -1.98796263, 3.42634911,   -1.92605903, -2.09950737, 0.88675595 ],
        [ -1.8173275,  3.09980177,   -1.90026412, -1.94992266, 0.86600291 ],
        [ -1.63280867, 2.75235075,   -1.86425792, -1.7924687,  0.840014 ],
        [ -1.45246261, 2.41999161,   -1.81999585, -1.64353995, 0.81098356 ],
        [ -1.28776495, 2.12412905,   -1.7715479,  -1.51236437, 0.78149341 ],
        [ -1.13365865, 1.85503746,   -1.71977327, -1.39420538, 0.75165468 ],
        [ -0.98743196, 1.60745349,   -1.66585343, -1.28646279, 0.72176909 ],
        [ -0.84664837, 1.37677861,   -1.61074985, -1.18693464, 0.6920331 ],
        [ -0.7064678,  1.15493953,   -1.5541187,  -1.09202829, 0.66199746 ],
        [ -0.56640251, 0.94125569,   -1.497124,   -1.00140516, 0.63207005 ],
        [ -0.42528419, 0.73399373,   -1.4405287,  -0.91430553, 0.60248074 ],
        [ -0.28188218, 0.53143406,   -1.38492727, -0.8300069,  0.57341535 ],
        [ -0.13469102, 0.3316158,    -1.33069066, -0.74771676, 0.54498559 ],
        [ 0.01746653,  0.13318166,   -1.27822541, -0.66691886, 0.51735916 ],
        [ 0.17600941,  -0.06541528,  -1.22778411, -0.58704098, 0.4906537 ],
        [ 0.3420835,   -0.26525577,  -1.17965638, -0.50771713, 0.46503042 ],
        [ 0.51719725,  -0.46775677,  -1.13398345, -0.42846187, 0.4405889 ],
        [ 0.70277613,  -0.67412219,  -1.09091176, -0.34888809, 0.41744209 ],
        [ 0.90050687,  -0.88574669,  -1.05051479, -0.26854643, 0.39566919 ],
        [ 1.11249592,  -1.10435071,  -1.0127964,  -0.18687801, 0.37531304 ],
        [ 1.33906397,  -1.32978056,  -0.97804358, -0.1040224,  0.35656751 ],
        [ 1.5735003,   -1.55536579,  -0.94725277, -0.02242381, 0.3400019 ],
        [ 1.81795929,  -1.78349367,  -0.91987368, 0.05885074,  0.32534084 ],
        [ 2.0739984,   -2.01581471,  -0.89553137, 0.14044541,  0.31239619 ],
        [ 2.34356433,  -2.25421856,  -0.87388578, 0.22307319,  0.30099364 ],
        [ 2.62818015,  -2.50012811,  -0.85469258, 0.30727412,  0.2910062 ],
        [ 2.92929359,  -2.75485778,  -0.83774583, 0.39354892,  0.28232401 ],
        [ 3.24899774,  -3.020229,    -0.82283584, 0.48256903,  0.27483327 ],
        [ 3.58885118,  -3.29758497,  -0.80980876, 0.57484561,  0.26844676 ],
        [ 3.95053901,  -3.58836956,  -0.79851832, 0.6709292,   0.2630788 ],
        [ 4.33595918,  -3.89419725,  -0.788823,   0.77143271,  0.25864409 ],
        [ 4.74754581,  -4.21710378,  -0.78058135, 0.87711291,  0.25505572 ],
        [ 5.18814238,  -4.55943654,  -0.77365875, 0.98883322,  0.2522289 ],
        [ 5.65925117,  -4.92250817,  -0.76794576, 1.10712515,  0.2500874 ],
        [ 6.1649713,   -5.30964309,  -0.76330215, 1.23318047,  0.24854154 ],
        [ 6.70938542,  -5.72413576,  -0.75960204, 1.36818366,  0.24750747 ],
        [ 7.29776409,  -6.17017328,  -0.7567205,  1.51361236,  0.2469024 ],
        [ 7.93736628,  -6.65343033,  -0.75453467, 1.67143205,  0.24664626 ],
        [ 8.63883985,  -7.18211268,  -0.75292612, 1.84444026,  0.24666393 ],
        [ 9.41642243,  -7.76709676,  -0.75178655, 2.03631889,  0.24688687 ],
        [ 10.28954316, -8.42313091,  -0.7510173,  2.25203568,  0.24725275 ],
        [ 11.28942455, -9.17378869,  -0.75052813, 2.49948654,  0.24770713 ],
        [ 12.46828365, -10.05836473, -0.75024077, 2.79179807,  0.24820422 ],
        [ 13.91153353, -11.14102419, -0.75009064, 3.15039456,  0.24870345 ],
        [ 15.79338592, -12.55250934, -0.75002487, 3.61888787,  0.24917143 ],
        [ 18.5406561,  -14.61299247, -0.75000366, 4.30404741,  0.24957771 ],
        [ 23.77878727, -18.54159591, -0.75000009, 5.61234344,  0.24988537 ]
    ],

    energy_table => [
        [ -6.29450485, 0.00142607, 0.00154762, 0.00142607, 0.00154762,    6.5282002,   -1.1762513 ],
        [ -5.71452454, 0.00267323, 0.0028905,  0.00267323, 0.0028905,     5.8404104,   -1.1966312 ],
        [ -5.44535993, 0.00357506, 0.00385602, 0.00357506, 0.00385602,    5.5169772,   -1.2073352 ],
        [ -5.09580423, 0.00520823, 0.00559334, 0.00520823, 0.00559334,    5.0923512,   -1.2234316 ],
        [ -4.86918644, 0.00664049, 0.0071055,  0.00664049, 0.0071055,     4.8138257,   -1.2356005 ],
        [ -4.38348092, 0.01113693, 0.01178614, 0.01113693, 0.01178614,    4.2068727,   -1.2664156 ],
        [ -3.99115447, 0.01682473, 0.01756698, 0.01682473, 0.01756698,    3.7047056,   -1.2962106 ],
        [ -3.66287218, 0.02364235, 0.02429465, 0.02364235, 0.02429465,    3.2747221,   -1.3261016 ],
        [ -3.37959995, 0.03154982, 0.03182436, 0.03154982, 0.03182436,    2.8950892,   -1.3562059 ],
        [ -3.12902158, 0.04051841, 0.04000681, 0.04051841, 0.04000681,    2.5516988,   -1.386377 ],
        [ -2.90419524, 0.05045929, 0.04862505, 0.05045929, 0.04862505,    2.2367801,   -1.4159174 ],
        [ -2.69885086, 0.06134015, 0.05750137, 0.06134015, 0.05750137,    1.9431791,   -1.4445608 ],
        [ -2.5082859,  0.07314146, 0.06645113, 0.07314146, 0.06645113,    1.6652856,   -1.4705011 ],
        [ -2.32869354, 0.0858651,  0.0752835,  0.0858651,  0.0752835,     1.3991233,   -1.4928372 ],
        [ -2.15658239, 0.0995574,  0.08380592, 0.0995574,  0.08380592,    1.1404078,   -1.5098461 ],
        [ -1.98796263, 0.11437309, 0.09183513, 0.11437309, 0.09183513,    0.88471805,  -1.5184785 ],
        [ -1.8173275,  0.13068498, 0.09918415, 0.13068498, 0.09918415,    0.62524798,  -1.5158217 ],
        [ -1.63280867, 0.14961312, 0.1056756,  0.14961312, 0.1056756,     0.34650646,  -1.4968626 ],
        [ -1.45246261, 0.16909457, 0.10999313, 0.16909457, 0.10999313,    0.078981686, -1.4631834 ],
        [ -1.28776495, 0.18738947, 0.11181169, 0.18738947, 0.11181169,    -0.15896096, -1.4233606 ],
        [ -1.13365865, 0.20462317, 0.11152543, 0.20462317, 0.11152543,    -0.3752293,  -1.3844434 ],
        [ -0.98743196, 0.2208034,  0.10950363, 0.2208034,  0.10950363,    -0.57504565, -1.3523388 ],
        [ -0.84664837, 0.23599542, 0.10609603, 0.23599542, 0.10609603,    -0.76351555, -1.3301115 ],
        [ -0.7064678,  0.25055906, 0.10151265, 0.25055906, 0.10151265,    -0.94876967, -1.3180796 ],
        [ -0.56640251, 0.26440236, 0.09603164, 0.26440236, 0.09603164,    -1.1329027,  -1.3155772 ],
        [ -0.42528419, 0.2775272,  0.0899042,  0.2775272,  0.0899042,     -1.3186904,  -1.3207581 ],
        [ -0.28188218, 0.28995232, 0.08335593, 0.28995232, 0.08335593,    -1.5087046,  -1.3314784 ],
        [ -0.13469102, 0.30172167, 0.07657334, 0.30172167, 0.07657334,    -1.7056584,  -1.3458006 ],
        [ 0.01746653,  0.31284901, 0.0697328,  0.31284901, 0.0697328,     -1.9116463,  -1.362037 ],
        [ 0.17600941,  0.32336309, 0.06297619, 0.32336309, 0.06297619,    -2.1289515,  -1.378911 ],
        [ 0.3420835,   0.33327039, 0.05643419, 0.33327039, 0.05643419,    -2.3593917,  -1.3954812 ],
        [ 0.51719725,  0.3425968,  0.05020022, 0.3425968,  0.05020022,    -2.6052196,  -1.4111454 ],
        [ 0.70277613,  0.35135822, 0.04435046, 0.35135822, 0.04435046,    -2.8685394,  -1.4255202 ],
        [ 0.90050687,  0.35957881, 0.03893348, 0.35957881, 0.03893348,    -3.1518011,  -1.4384085 ],
        [ 1.11249592,  0.36729188, 0.03397208, 0.36729188, 0.03397208,    -3.4580565,  -1.449752 ],
        [ 1.33906397,  0.37446833, 0.02951055, 0.37446833, 0.02951055,    -3.7877532,  -1.4595062 ],
        [ 1.5735003,   0.38092189, 0.02566402, 0.38092189, 0.02566402,    -4.1309646,  -1.4675162 ],
        [ 1.81795929,  0.38677638, 0.02233857, 0.38677638, 0.02233857,    -4.4906122,  -1.4740742 ],
        [ 2.0739984,   0.39211585, 0.0194621,  0.39211585, 0.0194621,     -4.8688041,  -1.4794247 ],
        [ 2.34356433,  0.39701512, 0.01696839, 0.39701512, 0.01696839,    -5.2682702,  -1.4837849 ],
        [ 2.62818015,  0.4015266,  0.01480494, 0.4015266,  0.01480494,    -5.6911494,  -1.4873259 ],
        [ 2.92929359,  0.40569242, 0.01292637, 0.40569242, 0.01292637,    -6.1394927,  -1.4901808 ],
        [ 3.24899774,  0.40955483, 0.01129001, 0.40955483, 0.01129001,    -6.616327,   -1.4924741 ],
        [ 3.58885118,  0.41314123, 0.00986249, 0.41314123, 0.00986249,    -7.1239064,  -1.4943038 ],
        [ 3.95053901,  0.4164753,  0.00861433, 0.4164753,  0.00861433,    -7.6646779,  -1.4957286 ],
        [ 4.33595918,  0.41957767, 0.00751965, 0.41957767, 0.00751965,    -8.2414064,  -1.496843 ],
        [ 4.74754581,  0.4224679,  0.00655542, 0.4224679,  0.00655542,    -8.8576979,  -1.4977198 ],
        [ 5.18814238,  0.42516229, 0.00570212, 0.42516229, 0.00570212,    -9.5177599,  -1.4983781 ],
        [ 5.65925117,  0.42766497, 0.00494606, 0.42779963, 0.0054611164,  -10.223794,  -1.4988737 ],
        [ 6.1649713,   0.42999063, 0.00427214, 0.4304545,  0.0050226055,  -10.981918,  -1.4992386 ],
        [ 6.70938542,  0.43214726, 0.0036692,  0.4330962,  0.0046482609,  -11.798207,  -1.4994947 ],
        [ 7.29776409,  0.43414197, 0.00312795, 0.43566412, 0.0041174521,  -12.680544,  -1.4996791 ],
        [ 7.93736628,  0.43598181, 0.00264061, 0.4381314,  0.0036176106,  -13.639792,  -1.4998064 ],
        [ 8.63883985,  0.43767468, 0.00220054, 0.44049448, 0.0031285074,  -14.691905,  -1.499893 ],
        [ 9.41642243,  0.43922577, 0.00180306, 0.44274088, 0.0026584783,  -15.858223,  -1.4999518 ],
        [ 10.28954316, 0.44063756, 0.00144436, 0.4448616,  0.002210576,   -17.167885,  -1.4999962 ],
        [ 11.28942455, 0.44191426, 0.00111903, 0.44684792, 0.0017728724,  -18.667725,  -1.5000365 ],
        [ 12.46828365, 0.44305999, 0.00083525, 0.44868048, 0.0013539424,  -20.436082,  -1.500091 ],
        [ 13.91153353, 0.44407158, 0.00058187, 0.45037283, 0.0010435354,  -22.601146,  -1.500208 ],
        [ 15.79338592, 0.44494491, 0.0003634,  0.45193901, 0.00065167611, -25.424506,  -1.5005765 ],
        [ 18.5406561,  0.44566706, 0.00018284, 0.453234,   0.00032786921, -29.548068,  -1.5 ],
        [ 23.77878727, 0.44620098, 4.936e-05,  0.45419143, 8.850532e-05,  -37.405265,  -1.5 ]
    ],

    impulse_table => [
        [
            -7.2018717, 13.815063,  -7.2024583,    0.040000242, 0,          -0.21135224,
            0,          -2.1951278, -0.0046553048, -0.51463747, 0.32930183, 0.99882958,
            -4.258279
        ],
        [
            -6.9077548, 13.226831,  -6.9085421,    0.044265988, 0,          -0.21109743,
            0,          -2.1951305, -0.0053928057, -0.51438266, 0.32930156, 0.99838834,
            -4.1054124
        ],
        [
            -6.5022897, 12.415904, -6.5034709,    0.050717606, 0,          -0.21059743,
            0,          -2.195137, -0.0066048111, -0.51388266, 0.32930087, 0.99749596,
            -3.8966001
        ],
        [
            -6.2146077, 11.840544,  -6.2161827,   0.055698764, 0,          -0.21009743,
            0,          -2.1951457, -0.007626579, -0.51338266, 0.32929994, 0.99657823,
            -3.7500504
        ],
        [
            -5.8091426, 11.029624,  -5.811506,     0.063271695, 0,          -0.20909743,
            0,          -2.1951702, -0.0093406135, -0.51238266, 0.32929735, 0.99469074,
            -3.5462425
        ],
        [
            -5.2983169,  10.008009,   -5.3022592, 0.073641546, 0, -0.20709746, 0, -2.1952484,
            -0.01205868, -0.51038266, 0.32928916, 0.99078402,  -3.2951623
        ],
        [
            -4.6051698, 8.6218793,  -4.6130684,   0.088721488, 0,        -0.20209777,
            0,          -2.1956149, -0.017053544, -0.50538266, 0.329251, 0.98063551,
            -2.968246
        ],
        [
            -3.9120226, 7.2362466, -3.9278791,   0.10371087,  0,          -0.19210183,
            0,          -2.197079, -0.024117269, -0.49538265, 0.32909835, 0.95980154,
            -2.6645655
        ],
        [
            -3.5065575, 6.4264156,  -3.5304235,   0.11152024,  0,          -0.18211692,
            0,          -2.1995125, -0.029537119, -0.48538265, 0.32884451, 0.9389646,
            -2.5026019
        ],
        [
            -2.9957318,  5.4082588,   -3.0357555, 0.11910617, 0, -0.16222431, 0, -2.207255,
            -0.03812874, -0.46538264, 0.32803727, 0.89838836, -2.3222189
        ],
        [
            -2.6592596, 4.7404897,  -2.7155923,   0.12209307,  0,          -0.14253241,
            0,          -2.2187608, -0.045102585, -0.44572593, 0.32684053, 0.85990438,
            -2.2227975
        ],
        [
            -2.3025847,  4.0378833,   -2.3835786, 0.12322296, 0, -0.11369536, 0, -2.2428411,
            -0.05385383, -0.41675698, 0.32435183, 0.80656703, -2.1398215
        ],
        [
            -2.0402204, 3.5271759,  -2.1459987,   0.12292512,  0,          -0.086219369,
            0,          -2.2747004, -0.061268026, -0.38882446, 0.32109624, 0.75845368,
            -2.0971621
        ],
        [
            -1.8971196, 3.2519518, -2.0194297,   0.12256287,  0,          -0.068956189,
            0,          -2.299986, -0.065653854, -0.37124556, 0.31854241, 0.72911769,
            -2.0816127
        ],
        [
            -1.6094375, 2.7088415,  -1.7728781,   0.12213097,  0,         -0.030561986,
            0,          -2.3758587, -0.075015881, -0.33173274, 0.3110224, 0.66439842,
            -2.0690528
        ],
        [
            -1.3862939, 2.3001723,  -1.5902571,   0.12283397, 0,          0.00015472281,
            0,          -2.4666959, -0.082340856, -0.3001615, 0.30224208, 0.61028232,
            -2.0781065
        ],
        [
            -1.2039724, 1.9768239, -1.4475466,   0.12451702,  0,          0.023547431,
            0,          -2.568098, -0.087831455, -0.27696832, 0.29262755, 0.56480847,
            -2.0983683
        ],
        [
            -1.0498217, 1.7121226,  -1.3318977,   0.12680299,  0,          0.040968727,
            0,          -2.6756486, -0.091745688, -0.26150549, 0.28254641, 0.52632617,
            -2.1245576
        ],
        [
            -0.9162903, 1.4899177,  -1.2356453,   0.12936478,  0,          0.053925783,
            0,          -2.7855497, -0.094411009, -0.25191719, 0.27229581, 0.49349826,
            -2.1537629
        ],
        [
            -0.79850727, 1.2997,     -1.1538672,   0.13198139, 0,          0.063670879,
            0,           -2.8949597, -0.096154206, -0.2464925, 0.26210419, 0.46526043,
            -2.1843052
        ],
        [
            -0.69314675, 1.1342732,  -1.0832288,   0.13452274, 0,          0.071127901,
            0,           -3.0019919, -0.097247983, -0.2440189, 0.25213969, 0.44077121,
            -2.2151938
        ],
        [
            -0.59783657, 0.98851706, -1.0213788,   0.1369213,   0,          0.076946017,
            0,           -3.1055269, -0.097896428, -0.24346586, 0.24252057, 0.41936558,
            -2.2458406
        ],
        [
            -0.51082519, 0.85867455, -0.9666034,   0.13914814,  0,          0.081574503,
            0,           -3.2049876, -0.098243496, -0.24415272, 0.23332545, 0.40051629,
            -2.2758996
        ],
        [
            -0.35667451, 0.63608537, -0.87345326,  0.14307182,  0,          0.08841504,
            0,           -3.3910082, -0.098396094, -0.24788803, 0.21637555, 0.36888864,
            -2.3335579
        ],
        [
            -0.22314312, 0.45073416, -0.79666497,  0.14634979,  0,          0.093192685,
            0,           -3.5603569, -0.098171119, -0.25288235, 0.20142824, 0.34341646,
            -2.3875068
        ],
        [
            -0.10536009, 0.29273848, -0.73181225,  0.14908925,  0,          0.09670672,
            0,           -3.7146124, -0.097783366, -0.25848847, 0.18841176, 0.32246443,
            -2.4377296
        ],
        [
            4.2974742e-07, 0.15555808, -0.67598193,  0.15139364,  0,          0.099399638,
            0,             -3.8556054, -0.097334018, -0.26421509, 0.17715153, 0.30491825,
            -2.4844531
        ],
        [
            0.09531061, 0.034671581, -0.62716986,  0.15334909,  0,          0.10153181,
            0,          -3.9850573,  -0.096871605, -0.27013865, 0.16743239, 0.28999724,
            -2.5279817
        ],
        [
            0.18232199, -0.073159707, -0.58394691,  0.15502386,  0,          0.10326476,
            0,          -4.1044761,   -0.096419364, -0.27613794, 0.15902083, 0.27714032,
            -2.5686287
        ],
        [
            0.26997834, -0.17946973, -0.54163591,  0.15660435,  0,          0.10483374,
            0,          -4.2257081,  -0.095946498, -0.28206779, 0.15101636, 0.26489928,
            -2.6103388
        ],
        [
            0.40546554,  -0.33947882, -0.47853422, 0.1588413,  0, 0.1069638, 0, -4.4144155,
            -0.09520204, -0.29244293, 0.13957291,  0.24731334, -2.6760786
        ],
        [
            0.53062868, -0.48296564, -0.42255602,  0.16069572,  0,          0.10866633,
            0,          -4.5896637,  -0.094518432, -0.30259185, 0.12996858, 0.23241646,
            -2.7379284
        ],
        [
            0.69314761, -0.6636082, -0.35291289,  0.16282295,  0,          0.11056952,
            0,          -4.8178369, -0.093659383, -0.31663552, 0.11878313, 0.21484285,
            -2.8194688
        ],
        [
            0.9932522,   -0.98237289, -0.23228388, 0.16602553, 0, 0.11338625, 0, -5.2389928,
            -0.09221138, -0.34516038, 0.10145064,  0.18698315, -2.972492
        ],
        [
            1.0986127, -1.0902739, -0.19209718,  0.16695856,  0,           0.11420779,
            0,         -5.3863764, -0.091753619, -0.35652664, 0.096241267, 0.17842997,
            -3.0266929
        ],
        [
            1.2089608, -1.2012971, -0.15107977,  0.16784428,  0,           0.11499316,
            0,         -5.5403293, -0.091303907, -0.36889252, 0.091204568, 0.17007056,
            -3.0836192
        ],
        [
            1.3862948, -1.3758179, -0.087265719, 0.16909351,  0,           0.1161156,
            0,         -5.7867177, -0.090644375, -0.39148766, 0.083909845, 0.15779529,
            -3.1753155
        ],
        [
            1.6094383, -1.5893306, -0.010246458, 0.17040507,  0,          0.11732488,
            0,         -6.0947558, -0.089920916, -0.42149155, 0.07593293, 0.14412602,
            -3.2908514
        ],
        [
            1.7917599, -1.759358,  0.050308023,  0.17129734,  0,           0.11817437,
            0,         -6.3446876, -0.089411535, -0.44741547, 0.070250181, 0.13421965,
            -3.3852274
        ],
        [
            1.9459106, -1.9003713, 0.10004062,   0.17194654,  0,           0.11881177,
            0,         -6.5547546, -0.089033281, -0.47295927, 0.065941252, 0.12660955,
            -3.4649369
        ],
        [
            2.079442,   -2.0206883,  0.14214528,  0.17244172, 0, 0.11931182, 0, -6.7358096,
            -0.0887408, -0.49674098, 0.062528768, 0.12052018, -3.5338943
        ],
        [
            2.3025855, -2.2183455, 0.21070622,   0.1731505,   0,           0.12005422,
            0,         -7.036532,  -0.088317688, -0.53822998, 0.057400108, 0.11126149,
            -3.6489039
        ],
        [
            2.7080506, -2.5681997, 0.33041733,   0.1741368,   0,           0.12116007,
            0,         -7.5774585, -0.087725231, -0.62858117, 0.049577847, 0.096886494,
            -3.8570811
        ],
        [
            2.9957327, -2.8104054, 0.41224981,   0.17465778,  0,           0.1217942,
            0,         -7.9573683, -0.087413274, -0.70350284, 0.044949794, 0.088235609,
            -4.0041532
        ],
        [
            3.4011978, -3.1449895, 0.52416362,   0.17520979, 0,           0.12252392,
            0,         -8.488098,  -0.087092697, -0.8353435, 0.039418172, 0.077755295,
            -4.2106215
        ],
        [
            3.9120234, -3.5575933, 0.66078681,   0.17568988, 0,           0.12323445,
            0,         -9.150171,  -0.086825728, -1.0312968, 0.033709456, 0.066785993,
            -4.4695857
        ],
        [
            4.2486968, -3.8252759, 0.74882394,   0.17591608, 0,           0.12360887,
            0,         -9.5832619, -0.086728206, -1.1984985, 0.030537526, 0.060187165,
            -4.6461181
        ],
        [
            4.6053434, -4.1059191, 0.84076509,   0.17609744, 0,           0.12393544,
            0,         -10.039743, -0.086633917, -1.405627,  0.027582974, 0.054553619,
            -4.8243827
        ],
        [
            5.2983453,   -4.644614, 1.0165977,   0.17633295,  0, 0.12441881, 0, -10.92141,
            -0.08633586, -1.906563, 0.022783881, 0.045268049, -5.1708835
        ],
        [
            6.2146602, -5.3475611, 1.2455273,    0.17650154, 0,           0.12485191,
            0,         -12.079581, -0.086442855, -2.9770813, 0.017857784, 0.035597778,
            -5.6291147
        ],
        [
            6.9078093, -5.8747497, 1.4172691,    0.17656955, 0,           0.12508119,
            0,         -12.951948, -0.086415055, -4.1369431, 0.014917188, 0.029775957,
            -5.9757426
        ],
        [
            7.6011782,   -6.3995953, 1.5884994,   0.17660927,  0, 0.12522473, 0, -13.822538,
            -0.21243703, -5.1872633, 0.012489346, 0.024949643, -6.3224674
        ],
        [
            8.5172016, -7.0905142, 1.8144378,    0.17663704, 0,            0.12535745,
            0,         -14.970688, -0.077684598, -8.0193676, 0.0098991492, 0.019786727,
            -6.780513
        ],
        [
            9.2103829, -7.612174,  1.9854582,    0.17664721, 0,            0.12542419,
            0,         -15.838573, -0.086382097, -12.744318, 0.0083113487, 0.016616925,
            -7.1271236
        ],
        [
            10.819792, -8.8212737, 2.3832051,    0.17665168, 0,            0.12551302,
            0,         -17.851943, -0.060226669, -19.545469, 0.0055485344, 0.011095924,
            -7.9318745
        ],
        [
            11.513038, -9.3416088, 2.5548884,    0.1766479,  0,            0.12553402,
            0,         -18.718813, -0.073717428, -27.453705, 0.0046640681, 0.009327571,
            -8.2785254
        ],
        [
            13.122508, -10.549161, 2.9542601,    0.17662515, 0,            0.12556197,
            0,         -20.730999, -0.080158287, -82.529757, 0.0031178614, 0.0062356315,
            -9.0833705
        ],
        [
            13.815531, -11.069013, 3.1265197,    0.17660684, 0,           0.12556856,
            0,         -21.59734,  -0.014401085, -20.826014, 0.002621699, 0.0052433742,
            -9.429963
        ],
        [
            16.118194, -12.796123, 3.699831,      0.1764693, 0,            0.12557933,
            0,         -24.47575,  -0.0082690397, -37.91275, 0.0014740914, 0.0029482838,
            -10.581892
        ],
        [
            18.420853, -14.52314, 4.274148,     0.17603298, 0,           0.12558202,
            0,         -27.35358, -0.086373426, -1256.6987, 0.000829197, 0.001658731,
            -11.735111
        ]
    ],

    tail_shock_table => [
        [ 4.7023259, -1.5532,    -0.029685238,  -0.029685238 ],
        [ 5.2844304, -1.8878628, -0.039353158,  -0.015307474 ],
        [ 6.2076867, -2.5261406, -0.035349878,  -0.0080172029 ],
        [ 6.9036751, -3.1160288, -0.031210663,  -0.0054049757 ],
        [ 7.5987276, -3.8204946, -0.027159441,  -0.0038041933 ],
        [ 8.5159733, -4.9668682, -0.02233609,   -0.0024831908 ],
        [ 9.2096545, -6.0335894, -0.019159164,  -0.001828982 ],
        [ 10.819575, -9.3816195, -0.013267178,  -0.00092946012 ],
        [ 11.512909, -11.309564, -0.011289522,  -0.00070065258 ],
        [ 13.12247,  -17.386867, -0.0077366559, -0.00036431208 ],
        [ 13.815508, -20.847101, -0.0065527139, -0.00027857317 ],
        [ 16.11819,  -38.09856,  -0.0037770338, -0.000114222 ]
    ],

    blast_info => [
        0.21209743, 0,          0.51538234, -0.17053545, 0,          0,          1.3357134,   -0.24717811,
        0.72803675, 0.14440277, 0.84902842, -0.13523853, 0.25490451, 0.32930212, 0.080213317, 0,
        108.65,     -1.5532,    178.08495,  -1.8283343
    ],

};
1;
