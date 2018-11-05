package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='P', gamma=1.25
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'P1x25'} = {

    table_name  => 'P1x25',
    symmetry    => 0,
    gamma       => 1.25,
    data_source => 'P8000_G1x25/moc2_from_moc_r1e5/',

    shock_table_info => [ 0, 1.25, 1.27e-06, 8000, 6.31e-07, 100, 4.54e-07, 3.2e-07, 5e-07, 1835.4, -4.5731 ],

    shock_table => [
        [ -12.81135743, 12.00017781, -0.99999317, -12.81310068, 0.99912762 ],
        [ -11.73808762, 10.92692396, -0.99998004, -11.74107082, 0.99850621 ],
        [ -10.88322119, 10.07209015, -0.99994461, -10.88779895, 0.99770601 ],
        [ -9.84932023,  9.03829229,  -0.99983992, -9.85700821,  0.99614187 ],
        [ -8.91446629,  8.10368575,  -0.99959289, -8.92676199,  0.99381696 ],
        [ -8.16083404,  7.35051107,  -0.99913624, -8.17880298,  0.9909429 ],
        [ -7.52468343,  6.71512759,  -0.99837213, -7.54945481,  0.98748201 ],
        [ -6.97101994,  6.16266576,  -0.99717901, -7.00380087,  0.98338952 ],
        [ -6.48109649,  5.67452069,  -0.99542151, -6.52312837,  0.97864418 ],
        [ -6.03866561,  5.23462746,  -0.99293066, -6.09130741,  0.97318453 ],
        [ -5.63279153,  4.83227276,  -0.9895088,  -5.69753929,  0.96694293 ],
        [ -5.2547256,   4.45899094,  -0.98491699, -5.33327278,  0.95982693 ],
        [ -4.89634409,  4.10704514,  -0.97884592, -4.99070242,  0.95169174 ],
        [ -4.55014818,  3.76948334,  -0.97088579, -4.66280563,  0.94232815 ],
        [ -4.20707679,  3.43811395,  -0.96042101, -4.3413476,   0.93137856 ],
        [ -3.85213224,  3.09961511,  -0.94630863, -4.01305557,  0.91809687 ],
        [ -3.45719395,  2.72971503,  -0.92605411, -3.65379064,  0.90079682 ],
        [ -3.10287431,  2.40546865,  -0.9034984,  -3.33771703,  0.88295596 ],
        [ -2.77573792,  2.11380496,  -0.87908825, -3.05182339,  0.86461435 ],
        [ -2.47141047,  1.85008032,  -0.8536902,  -2.79147994,  0.84611775 ],
        [ -2.18479109,  1.60904896,  -0.82796318, -2.55158973,  0.82766637 ],
        [ -1.90420993,  1.38040245,  -0.80172639, -2.32198171,  0.80890202 ],
        [ -1.62112916,  1.15724692,  -0.77489059, -2.09572832,  0.7895628 ],
        [ -1.32911589,  0.93497962,  -0.74752907, -1.86809507,  0.76950713 ],
        [ -1.0318123,   0.71676283,  -0.7206573,  -1.64233038,  0.74930344 ],
        [ -0.72569561,  0.50019681,  -0.69457156, -1.41607542,  0.72903745 ],
        [ -0.40677948,  0.28274464,  -0.66950669, -1.18683016,  0.70878031 ],
        [ -0.07097992,  0.06199941,  -0.64568758, -0.95224216,  0.68862669 ],
        [ 0.28611046,   -0.16448601, -0.62332261, -0.70994728,  0.66869082 ],
        [ 0.67060629,   -0.40004957, -0.60254056, -0.45667612,  0.6490469 ],
        [ 1.08979909,   -0.64850904, -0.5834675,  -0.18871406,  0.62979115 ],
        [ 1.50502572,   -0.88743018, -0.5678333,  0.06920484,   0.61285747 ],
        [ 1.92580093,   -1.1235323,  -0.5548289,  0.32382771,   0.59772449 ],
        [ 2.35851376,   -1.36117928, -0.54396181, 0.57945642,   0.58410363 ],
        [ 2.80803448,   -1.60358698, -0.53489516, 0.8391955,    0.57182191 ],
        [ 3.2785005,    -1.85339727, -0.5273746,  1.10555037,   0.56076715 ],
        [ 3.774264,     -2.11325122, -0.52118686, 1.38103191,   0.55085063 ],
        [ 4.29973714,   -2.38573773, -0.51615163, 1.66809608,   0.54200609 ],
        [ 4.86073072,   -2.67410429, -0.51210325, 1.96988803,   0.53416595 ],
        [ 5.46298327,   -2.9815046,  -0.5089013,  2.28944642,   0.52728556 ],
        [ 6.11591065,   -3.31291948, -0.50640811, 2.63169662,   0.52130077 ],
        [ 6.82570334,   -3.67165099, -0.50451637, 2.9998234,    0.51619563 ],
        [ 7.61044829,   -4.06697547, -0.50310368, 3.40313422,   0.51188753 ],
        [ 8.48589587,   -4.50693376, -0.50207774, 3.84963114,   0.5083476 ],
        [ 9.47178144,   -5.00153857, -0.50135204, 4.34933312,   0.50554081 ],
        [ 10.62908783,  -5.58143212, -0.50083601, 4.93304476,   0.50336705 ],
        [ 12.04507471,  -6.29032501, -0.50046912, 5.64457578,   0.50178458 ],
        [ 13.80713903,  -7.17194168, -0.50022214, 6.52775543,   0.50078225 ],
        [ 15.90569283,  -8.22152021, -0.50008523, 7.57806504,   0.50028289 ],
        [ 18.38258847,  -9.46009166, -0.50002563, 8.81691744,   0.50008315 ],
        [ 21.91502323,  -11.2263521, -0.50000448, 10.58327305,  0.50001429 ]
    ],

    energy_table => [
        [ -12.81135743, 0.08956153, 0.01790464, 0.08956153, 0.01790464,    11.518789,   -1.0017115 ],
        [ -11.73808762, 0.11098956, 0.02217602, 0.11098956, 0.02217602,    10.443225,   -1.0028734 ],
        [ -10.88322119, 0.13165319, 0.02628063, 0.13165319, 0.02628063,    9.5854012,   -1.0044222 ],
        [ -9.84932023,  0.16179907, 0.03222551, 0.16179907, 0.03222551,    8.5457238,   -1.0074731 ],
        [ -8.91446629,  0.19484496, 0.03864489, 0.19484496, 0.03864489,    7.6022894,   -1.0118958 ],
        [ -8.16083404,  0.2261634,  0.04458054, 0.2261634,  0.04458054,    6.8380411,   -1.0173398 ],
        [ -7.52468343,  0.25625234, 0.05008525, 0.25625234, 0.05008525,    6.1891128,   -1.0239521 ],
        [ -6.97101994,  0.28537864, 0.05516228, 0.28537864, 0.05516228,    5.6203249,   -1.0318347 ],
        [ -6.48109649,  0.31352991, 0.05976414, 0.31352991, 0.05976414,    5.1128464,   -1.0410666 ],
        [ -6.03866561,  0.34088539, 0.06387591, 0.34088539, 0.06387591,    4.6501551,   -1.0517661 ],
        [ -5.63279153,  0.36754839, 0.06746751, 0.36754839, 0.06746751,    4.2210466,   -1.0640535 ],
        [ -5.2547256,   0.39364216, 0.07050824, 0.39364216, 0.07050824,    3.8163666,   -1.078143 ],
        [ -4.89634409,  0.4193664,  0.07296793, 0.4193664,  0.07296793,    3.427348,    -1.0940728 ],
        [ -4.55014818,  0.44496279, 0.07480317, 0.44496279, 0.07480317,    3.0457133,   -1.1122335 ],
        [ -4.20707679,  0.47084355, 0.07595034, 0.47084355, 0.07595034,    2.6607832,   -1.1323178 ],
        [ -3.85213224,  0.49788928, 0.07628901, 0.49788928, 0.07628901,    2.2550885,   -1.1543033 ],
        [ -3.45719395,  0.52790679, 0.07550999, 0.52790679, 0.07550999,    1.7942343,   -1.176761 ],
        [ -3.10287431,  0.55437337, 0.07370891, 0.55437337, 0.07370891,    1.3741507,   -1.1856015 ],
        [ -2.77573792,  0.57808822, 0.07113457, 0.57808822, 0.07113457,    0.9862988,   -1.1724151 ],
        [ -2.47141047,  0.5992792,  0.06802251, 0.5992792,  0.06802251,    0.6332326,   -1.12693 ],
        [ -2.18479109,  0.61828998, 0.0645552,  0.61828998, 0.0645552,     0.31920071,  -1.0533228 ],
        [ -1.90420993,  0.63587927, 0.0607671,  0.63587927, 0.0607671,     0.035281673, -0.97909373 ],
        [ -1.62112916,  0.65250648, 0.05667032, 0.65250648, 0.05667032,    -0.23251243, -0.91435543 ],
        [ -1.32911589,  0.66841683, 0.05228294, 0.66841683, 0.05228294,    -0.48998446, -0.86855781 ],
        [ -1.0318123,   0.68329008, 0.04777417, 0.68329008, 0.04777417,    -0.74422774, -0.85148603 ],
        [ -0.72569561,  0.69721196, 0.04320491, 0.69721196, 0.04320491,    -1.0037226,  -0.85084173 ],
        [ -0.40677948,  0.71025456, 0.03862663, 0.71025456, 0.03862663,    -1.2761139,  -0.85910116 ],
        [ -0.07097992,  0.7224548,  0.03409165, 0.7224548,  0.03409165,    -1.5663624,  -0.871464 ],
        [ 0.28611046,   0.73382374, 0.02965338, 0.73382374, 0.02965338,    -1.880255,   -0.88655796 ],
        [ 0.67060629,   0.74438259, 0.02535392, 0.74438259, 0.02535392,    -2.2242501,  -0.9021446 ],
        [ 1.08979909,   0.75412764, 0.02123864, 0.75412764, 0.02123864,    -2.6058408,  -0.91754023 ],
        [ 1.50502572,   0.76219847, 0.01772583, 0.76219847, 0.01772583,    -2.9898063,  -0.93090521 ],
        [ 1.92580093,   0.76900122, 0.01469235, 0.76900122, 0.01469235,    -3.3841484,  -0.94251317 ],
        [ 2.35851376,   0.77477372, 0.01206634, 0.77477372, 0.01206634,    -3.7943592,  -0.95259796 ],
        [ 2.80803448,   0.77967213, 0.00980029, 0.77967213, 0.00980029,    -4.22472,    -0.96133392 ],
        [ 3.2785005,    0.78381027, 0.00785852, 0.78381027, 0.00785852,    -4.6789433,  -0.96885608 ],
        [ 3.774264,     0.78728219, 0.00620946, 0.78728219, 0.00620946,    -5.161033,   -0.97527924 ],
        [ 4.29973714,   0.79016669, 0.00482512, 0.79016669, 0.00482512,    -5.6751092,  -0.98068997 ],
        [ 4.86073072,   0.79253744, 0.00367721, 0.79253744, 0.00367721,    -6.2266947,  -0.98520755 ],
        [ 5.46298327,   0.79445664, 0.00274092, 0.79445664, 0.00274092,    -6.8213173,  -0.98892769 ],
        [ 6.11591065,   0.79598797, 0.00198912, 0.79598797, 0.00198912,    -7.4681455,  -0.99192575 ],
        [ 6.82570334,   0.7971792,  0.00140122, 0.7971792,  0.00140122,    -8.1731848,  -0.99428068 ],
        [ 7.61044829,   0.79809027, 0.00094967, 0.7980952,  0.00094967,    -8.9542897,  -0.99609845 ],
        [ 8.48589587,   0.79876438, 0.00061445, 0.79876676, 0.00061445,    -9.8270415,  -0.99746137 ],
        [ 9.47178144,   0.79924299, 0.00037586, 0.7992449,  0.00037622239, -10.811025,  -0.99844704 ],
        [ 10.62908783,  0.79957337, 0.00021087, 0.79957592, 0.00021114767, -11.967046,  -0.99933819 ],
        [ 12.04507471,  0.79978738, 0.0001039,  0.79979034, 0.00010468505, -13.382873,  -1 ],
        [ 13.80713903,  0.79990907, 4.305e-05,  0.79991322, 4.3408051e-05, -15.144937,  -1 ],
        [ 15.90569283,  0.79996501, 1.507e-05,  0.79996962, 1.519236e-05,  -17.243491,  -1 ],
        [ 18.38258847,  0.79998642, 4.37e-06,   0.7999912,  4.4021895e-06, -19.720387,  -1 ],
        [ 21.91502323,  0.79999366, 7.5e-07,    0.79999849, 7.5261558e-07, -23.252821,  -1 ]
    ],

    impulse_table => [
        [
            -14.626074, 13.814889,  -14.626777,    0.2750935,  0,          -0.84626714,
            0,          -8.6285865, -0.0036465601, -1.4199675, 0.14328994, 0.92211732,
            -4.1590617
        ],
        [
            -13.815512, 13.004328,  -13.816567,    0.27491468, 0,          -0.84626658,
            0,          -8.2233081, -0.0036465601, -1.419967,  0.14328976, 0.90841279,
            -3.9969705
        ],
        [
            -13.122365, 12.311184,  -13.123857,    0.2746924, 0,          -0.84626558,
            0,          -7.8767381, -0.0036465601, -1.419966, 0.14328933, 0.8947977,
            -3.8583758
        ],
        [
            -12.7169, 11.905721, -12.718728,    0.27452183, 0,          -0.84626458,
            0,        -7.674009, -0.0036465601, -1.419965,  0.14328893, 0.88591471,
            -3.7773143
        ],
        [
            -12.206074, 11.394902,  -12.208434,    0.27425136, 0,          -0.84626258,
            0,          -7.4186027, -0.0036465601, -1.419963,  0.14328814, 0.87365025,
            -3.6752065
        ],
        [
            -11.512927, 10.701769,  -11.516266,    0.27375435, 0,          -0.84625758,
            0,          -7.0720453, -0.0036465601, -1.419958,  0.14328622, 0.85488046,
            -3.5367024
        ],
        [
            -10.81978, 10.008653, -10.824506,    0.27305151, 0,          -0.84624758,
            0,         -6.725504, -0.0036465601, -1.419948,  0.14328242, 0.8333371,
            -3.3982867
        ],
        [
            -10.414315, 9.6032174,  -10.420104,    0.27251225, 0,          -0.84623758,
            0,          -6.5228038, -0.0036465601, -1.419938,  0.14327862, 0.81929294,
            -3.3173823
        ],
        [
            -9.9034892, 9.0924528,  -9.910971,     0.27165724, 0,          -0.84621758,
            0,          -6.2674555, -0.0036465601, -1.419918,  0.14327101, 0.7999225,
            -3.215554
        ],
        [
            -9.2103421, 8.3994568,  -9.2209375,    0.27008654, 0,          -0.84616758,
            0,          -5.9210434, -0.0036465601, -1.419868,  0.14325198, 0.7703375,
            -3.0776458
        ],
        [
            -8.5171949, 7.7066124,  -8.5322103,    0.26786682, 0,          -0.84606758,
            0,          -5.5747929, -0.0036465601, -1.419768,  0.14321393, 0.7365058,
            -2.9402213
        ],
        [
            -8.1117298, 7.3014502,  -8.130149,     0.26616524, 0,          -0.84596758,
            0,          -5.3723835, -0.0036465601, -1.419668,  0.14317591, 0.71454769,
            -2.8601738
        ],
        [
            -7.6009041, 6.7912287, -7.6247396,    0.26347078, 0,          -0.84576758,
            0,          -5.117617, -0.0036465601, -1.419468,  0.14309993, 0.6844249,
            -2.7598448
        ],
        [
            -6.907757, 6.0995869,  -6.9416058,    0.25853565, 0,          -0.84526758,
            0,         -4.7726601, -0.0036465601, -1.418968,  0.14291049, 0.63889202,
            -2.6250478
        ],
        [
            -6.2146098, 5.4094286,  -6.2627407,    0.25160623, 0,          -0.84426758,
            0,          -4.4293229, -0.0036465601, -1.417968,  0.14253369, 0.58776465,
            -2.4926006
        ],
        [
            -5.8091447, 5.0069247,  -5.8683197,    0.24634164, 0,          -0.84326759,
            0,          -4.2298305, -0.0036465601, -1.416968,  0.14215967, 0.55527045,
            -2.4167014
        ],
        [
            -5.2983191, 4.5019403,  -5.3751344,    0.23811409, 0,          -0.84126759,
            0,          -3.9809104, -0.0036465601, -1.414968,  0.14141983, 0.51178121,
            -2.3233537
        ],
        [
            -4.6051719, 3.8229444,  -4.7147001,    0.22350254, 0,         -0.8362676,
            0,          -3.6506433, -0.0036465601, -1.409968,  0.1396162, 0.44893361,
            -2.2021477
        ],
        [
            -3.9120247, 3.1563713,  -4.0681144,    0.20436739, 0,          -0.82626775,
            0,          -3.3370404, -0.0036465599, -1.399968,  0.13619119, 0.38335088,
            -2.0894963
        ],
        [
            -3.5065596, 2.7754999,  -3.6983163,    0.19126064, 0,          -0.8162692,
            0,          -3.1678176, -0.0036465581, -1.389968,  0.13298423, 0.34478174,
            -2.0286357
        ],
        [
            -2.995734, 2.3090733,  -3.2434269,    0.17399417, 0,          -0.79630119,
            0,         -2.9812666, -0.0036465174, -1.3699681, 0.12712916, 0.29721123,
            -1.9581842
        ],
        [
            -2.3025868, 1.7072148,  -2.6495395,   0.15591616, 0,          -0.74821545,
            0,          -2.8173579, -0.003644082, -1.3220354, 0.11499252, 0.23675276,
            -1.8748131
        ],
        [
            -1.6094396,   1.1481953,  -2.0865034,  0.15682367, 0, -0.6889107, 0, -2.8263769,
            -0.003591665, -1.2608022, 0.097583376, 0.18353927, -1.8058687
        ],
        [
            -1.2039745, 0.84215266, -1.7723335,    0.16508865, 0,          -0.67229498,
            0,          -2.9017105, -0.0034849142, -1.2424981, 0.08555519, 0.15644756,
            -1.771812
        ],
        [
            -0.69314886, 0.47763417, -1.392382,     0.17827007, 0,           -0.67445734,
            0,           -3.0371273, -0.0032488697, -1.2478961, 0.069943724, 0.12673454,
            -1.7346984
        ],
        [
            -1.6839784e-06, 0.016336689, -0.90351018,   0.19626638, 0,           -0.70129178,
            0,              -3.2637673,  -0.0028169284, -1.2880388, 0.050991569, 0.093943569,
            -1.6930711
        ],
        [
            0.26236258, -0.1496669, -0.72584241,   0.2025343,  0,           -0.71369106,
            0,          -3.3584386, -0.0026520201, -1.3048489, 0.044980584, 0.083600355,
            -1.6795725
        ],
        [
            0.40546342, -0.23847212, -0.63051496,   0.2057849,  0,           -0.72031426,
            0,          -3.4117244,  -0.0025668818, -1.3135453, 0.041984174, 0.078395854,
            -1.6726849
        ],
        [
            0.53062657, -0.31520851, -0.54801164,   0.20852582, 0,           -0.72589969,
            0,          -3.4592038,  -0.0024962154, -1.320788,  0.039519719, 0.074084078,
            -1.6669251
        ],
        [
            0.6931455, -0.41361777, -0.44205944,   0.21193961, 0,           -0.7327533,
            0,         -3.5219898,  -0.0024105643, -1.3283188, 0.036526928, 0.06880574,
            -1.6598035
        ],
        [
            0.99325009, -0.59198041, -0.24972339,   0.21780961, 0,           -0.74394031,
            0,          -3.6410235,  -0.0022724889, -1.3416627, 0.031569458, 0.05995205,
            -1.6476725
        ],
        [
            1.0986106, -0.65364867, -0.18316632,   0.219738,   0,          -0.74736885,
            0,         -3.6836919,  -0.0022304187, -1.3445976, 0.02999026, 0.057101504,
            -1.6437151
        ],
        [
            1.2089587, -0.71774664, -0.11397223,   0.22168515, 0,           -0.75067063,
            0,         -3.7288358,  -0.0021898925, -1.3485081, 0.028419439, 0.05425121,
            -1.6397325
        ],
        [
            1.3862927, -0.8197649, -0.0038348583, 0.22466175, 0,          -0.75537817,
            0,         -3.8023108, -0.0021320935, -1.3531954, 0.02606288, 0.04994698,
            -1.6336691
        ],
        [
            1.5668613, -0.92247799, 0.10702844,    0.22750466, 0,           -0.75943663,
            0,         -3.8782301,  -0.0020819724, -1.3574666, 0.023861142, 0.045894555,
            -1.6279062
        ],
        [
            1.7917578, -1.0489057, 0.24340018,    0.23079067, 0,           -0.76355328,
            0,         -3.9742317, -0.0020306386, -1.3584222, 0.021373814, 0.041280058,
            -1.6212803
        ],
        [
            1.9459085, -1.1346829, 0.33583969,    0.23288679, 0,           -0.76583159,
            0,         -4.0408956, -0.0020017513, -1.3627064, 0.019818525, 0.038374874,
            -1.6170745
        ],
        [
            2.0794399, -1.2084551, 0.41527049,    0.23460422, 0,           -0.76750839,
            0,         -4.0991728, -0.0019803326, -1.36439,   0.018561962, 0.036016491,
            -1.6136411
        ],
        [
            2.3025834, -1.3307198, 0.54674115,    0.23728003, 0,           -0.76973886,
            0,         -4.1975891, -0.0019510027, -1.3715763, 0.016635459, 0.032381204,
            -1.6083161
        ],
        [
            2.7080485, -1.5500138, 0.78189326,    0.24156546, 0,           -0.77244205,
            0,         -4.3794101, -0.0019139924, -1.374225,  0.013627282, 0.026657314,
            -1.5998552
        ],
        [
            2.9957306, -1.7036784, 0.9460867,     0.2441971,  0,           -0.77354792,
            0,         -4.5104947, -0.0018970109, -1.3867854, 0.011825922, 0.023201879,
            -1.5947061
        ],
        [
            3.4011957, -1.9179996, 1.1741928,     0.24739757, 0,            -0.77459706,
            0,         -4.6977876, -0.0018815291, -1.3979325, 0.0096806417, 0.019059387,
            -1.5884986
        ],
        [
            3.9120213, -2.1849474, 1.4567441,     0.2506989,  0,            -0.77518417,
            0,         -4.9373549, -0.0018703544, -1.4205509, 0.0075190011, 0.014855029,
            -1.5821776
        ],
        [
            4.2484936, -2.3592773, 1.6403015,     0.25249503, 0,            -0.77533129,
            0,         -5.0970056, -0.0018660753, -1.4522195, 0.0063642899, 0.012596477,
            -1.5787943
        ],
        [
            4.6056896, -2.5432862, 1.8332298,     0.25412463, 0,            -0.77486079,
            0,         -5.2678508, -0.0018644592, -1.524576,  0.0053306723, 0.010107953,
            -1.6152376
        ],
        [
            5.2984015, -2.8976861, 2.2025234,     0.25660718, 0,            -0.77470661,
            0,         -5.6024562, -0.0018620851, -1.6919388, 0.0037780741, 0.0072748268,
            -1.5991964
        ],
        [
            6.2146542, -3.3629088, 2.6831326,     0.25885857, 0,            -0.77446807,
            0,         -6.0500735, -0.0018610725, -2.195189,  0.0023938587, 0.0046735378,
            -1.5847317
        ],
        [
            6.9078118, -3.7130688, 3.0421865,     0.26001066, 0,            -0.77432191,
            0,         -6.3913996, -0.0018608439, -3.0344627, 0.0016941092, 0.0033307115,
            -1.5773608
        ],
        [
            7.6009756, -4.0622097, 3.3982851,     0.26083231, 0,            -0.77421163,
            0,         -6.7343494, -0.0018607617, -4.7133238, 0.0011984888, 0.0023680793,
            -1.5721162
        ],
        [
            8.5172508, -4.5226759, 3.8655687,     0.26156573, 0,             -0.77411172,
            0,         -7.189406,  -0.0018475825, -9.6252418, 0.00075820673, 0.0015048047,
            -1.5674408
        ],
        [
            9.2105278, -4.8705384, 4.2171767,    0.2619362,  0,             -0.77400971,
            0,         -7.5345966, -0.001717192, -16.612197, 0.00053611433, 0.0010664136,
            -1.565077
        ],
        [
            10.819805, -5.6769442, 5.0290198,    0.26242284, 0,             -0.77338061,
            0,         -8.3375448, -0.001860774, -85.298429, 0.00023972735, 0.00047828931,
            -1.5619219
        ]
    ],

    tail_shock_table => [
        [ 7.5175062, -4.5731,    -0.0012777516,  -0.0012777516 ],
        [ 7.5859074, -4.7595026, -0.0014463482,  -0.0010603612 ],
        [ 8.5076594, -8.0956423, -0.0012653974,  -0.00038698435 ],
        [ 9.2037218, -11.988205, -0.00098838827, -0.0002087572 ],
        [ 10.816745, -29.078485, -0.00050351655, -5.5171137e-05 ],
        [ 11.510974, -43.488121, -0.00038004115, -2.7866016e-05 ],
        [ 13.121452, -98.653266, -0.00017413812, -9.0927334e-06 ]
    ],

    blast_info => [
        0.84626432, 0,            1.4204978, -0.0036465714, 0,          0,         0.2750574,  -0.020886579,
        3.6772206,  0.0011750898, 1.9708764, -0.0031657852, 0.41997341, 0.1432897, 0.21005863, 0,
        1835.4,     -4.5731,      2152.3607, -4.9687952
    ],

};
1;