package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='P', gamma=3
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'P3'} = {

    table_name  => 'P3',
    symmetry    => 0,
    gamma       => 3,
    data_source => 'P8000_G3_inspirion_ZR1_NLD0/moc3p_from_AD_r100/',

    shock_table_info => [ 0, 3, 7.89e-07, 8000, 4e-07, 100, 6.9e-08, 2.2e-07, 5e-07 ],

    shock_table => [
        [ -10.85863198, 12.00017826, -0.99999079, -10.86065774, 0.9989861 ],
        [ -9.09234545,  10.23393348, -0.99994949, -9.09725166,  0.99754102 ],
        [ -7.89815703,  9.03986173,  -0.99983293, -7.90708808,  0.9955154 ],
        [ -6.98361242,  8.12556728,  -0.99958303, -6.99775589,  0.9928817 ],
        [ -6.23872881,  7.38114484,  -0.99912252, -6.25931501,  0.98961152 ],
        [ -5.60445291,  6.74764523,  -0.99834905, -5.63281778,  0.98564393 ],
        [ -5.05391392,  6.19831624,  -0.99714692, -5.09140906,  0.98096411 ],
        [ -4.56586518,  5.71205785,  -0.99537584, -4.61392122,  0.97552658 ],
        [ -4.12593391,  5.27467256,  -0.99287287, -4.18607789,  0.96927998 ],
        [ -3.72255944,  4.87482145,  -0.98943893, -3.79648427,  0.96214199 ],
        [ -3.34730557,  4.50434387,  -0.98483848, -3.43691218,  0.95401632 ],
        [ -2.99267348,  4.15610667,  -0.9787755,  -3.1001791,   0.94476291 ],
        [ -2.65147679,  3.82343637,  -0.97085898, -2.7795852,   0.93417388 ],
        [ -2.31559177,  3.49899674,  -0.96052255, -2.46781565,  0.92191125 ],
        [ -1.9723243,   3.17154127,  -0.94676484, -2.15379939,  0.90729577 ],
        [ -1.58912869,  2.81238453,  -0.92693591, -1.80967053,  0.88833184 ],
        [ -1.2375677,   2.49037598,  -0.90420626, -1.50079217,  0.8684658 ],
        [ -0.91441613,  2.20207798,  -0.8794965,  -1.223362,    0.84827535 ],
        [ -0.61159753,  1.93963898,  -0.8533799,  -0.96953629,  0.82794631 ],
        [ -0.33042915,  1.70334991,  -0.82711715, -0.73950812,  0.80816882 ],
        [ -0.05754453,  1.48126524,  -0.80042489, -0.52165027,  0.78848263 ],
        [ 0.21621611,   1.26586954,  -0.77315859, -0.30851751,  0.7685976 ],
        [ 0.49857725,   1.0515125,   -0.74525544, -0.09436939,  0.74830699 ],
        [ 0.78750733,   0.84020696,  -0.71762832, 0.11890506,   0.72812158 ],
        [ 1.08479987,   0.63089915,  -0.69077067, 0.33239272,   0.70826765 ],
        [ 1.39410554,   0.42129131,  -0.66496362, 0.54842505,   0.68884329 ],
        [ 1.71973308,   0.20883136,  -0.64042442, 0.76960691,   0.6699235 ],
        [ 2.06620375,   -0.00897192, -0.61736509, 0.99848757,   0.65160153 ],
        [ 2.43943861,   -0.23528806, -0.59593337, 1.23832606,   0.63393898 ],
        [ 2.84615352,   -0.47354356, -0.57628572, 1.49263858,   0.61701935 ],
        [ 3.24944253,   -0.70260066, -0.56018013, 1.73846562,   0.60243633 ],
        [ 3.6583446,    -0.92883059, -0.54679867, 1.98212232,   0.58963673 ],
        [ 4.07887183,   -1.15634454, -0.53564449, 2.22763602,   0.57829695 ],
        [ 4.5156409,    -1.38819584, -0.52637664, 2.47795657,   0.56820499 ],
        [ 4.9727892,    -1.62700948, -0.51873516, 2.73559761,   0.55920584 ],
        [ 5.45410897,   -1.87512257, -0.5125114,  3.0027691,    0.55118292 ],
        [ 5.96348947,   -2.13485449, -0.50752527, 3.28165841,   0.54404066 ],
        [ 6.50597906,   -2.40906315, -0.5036104,  3.57501899,   0.53768871 ],
        [ 7.0877743,    -2.70113866, -0.50061755, 3.88614998,   0.53204794 ],
        [ 7.71561374,   -3.01470761, -0.49841575, 4.21857065,   0.52705696 ],
        [ 8.39481396,   -3.35267148, -0.49688948, 4.5750076,    0.52267864 ],
        [ 9.14481241,   -3.72493392, -0.49591334, 4.96550989,   0.51881231 ],
        [ 9.97967552,   -4.1387063,  -0.4954004,  5.39717655,   0.51543051 ],
        [ 10.92269053,  -4.60578612, -0.49526948, 5.88178557,   0.51248944 ],
        [ 11.99870649,  -5.13877832, -0.49544956, 6.43180622,   0.50996279 ],
        [ 13.25318819,  -5.76056337, -0.49587678, 7.07011413,   0.50779638 ],
        [ 14.75837313,  -6.50740916, -0.49649382, 7.83295807,   0.5059324 ],
        [ 16.71309072,  -7.47868947, -0.49727192, 8.82017078,   0.5042595 ],
        [ 17.97592944,  -8.1069484,  -0.49771609, 9.45645104,   0.50346907 ]
    ],

    energy_table => [
        [ -10.85863198, 0.00038429, 0.00025299, 0.00038429, 0.00025299, 11.566591,    -0.97595506 ],
        [ -9.09234545,  0.00122252, 0.0007963,  0.00122252, 0.0007963,  9.8199026,    -1.0001003 ],
        [ -7.89815703,  0.00264948, 0.00170472, 0.00264948, 0.00170472, 8.6165548,    -1.0136882 ],
        [ -6.98361242,  0.00475377, 0.00301611, 0.00475377, 0.00301611, 7.6852766,    -1.0221594 ],
        [ -6.23872881,  0.00759891, 0.00474547, 0.00759891, 0.00474547, 6.9215435,    -1.0279789 ],
        [ -5.60445291,  0.01125411, 0.0069033,  0.01125411, 0.0069033,  6.2680772,    -1.0321867 ],
        [ -5.05391392,  0.01572436, 0.00945415, 0.01572436, 0.00945415, 5.6988954,    -1.035266 ],
        [ -4.56586518,  0.02102086, 0.0123598,  0.02102086, 0.0123598,  5.1930281,    -1.0375264 ],
        [ -4.12593391,  0.02714274, 0.01556828, 0.02714274, 0.01556828, 4.7361847,    -1.0391623 ],
        [ -3.72255944,  0.03410362, 0.0190274,  0.03410362, 0.0190274,  4.316747,     -1.0403053 ],
        [ -3.34730557,  0.04191575, 0.02267453, 0.04191575, 0.02267453, 3.9261997,    -1.041042 ],
        [ -2.99267348,  0.05061766, 0.02644712, 0.05061766, 0.02644712, 3.5569166,    -1.04143 ],
        [ -2.65147679,  0.0602917,  0.03028313, 0.0602917,  0.03028313, 3.2015458,    -1.0415031 ],
        [ -2.31559177,  0.07110892, 0.03412457, 0.07110892, 0.03412457, 2.8517335,    -1.0412693 ],
        [ -1.9723243,   0.0834832,  0.03793566, 0.0834832,  0.03793566, 2.4943681,    -1.0407085 ],
        [ -1.58912869,  0.09877531, 0.04178013, 0.09877531, 0.04178013, 2.0957278,    -1.0393167 ],
        [ -1.2375677,   0.11399088, 0.04465381, 0.11399088, 0.04465381, 1.7306633,    -1.0354675 ],
        [ -0.91441613,  0.12874711, 0.04653912, 0.12874711, 0.04653912, 1.3969244,    -1.0285526 ],
        [ -0.61159753,  0.14300998, 0.04753053, 0.14300998, 0.04753053, 1.0866543,    -1.0196407 ],
        [ -0.33042915,  0.15641921, 0.04773643, 0.15641921, 0.04773643, 0.80126017,   -1.0099447 ],
        [ -0.05754453,  0.1693987,  0.04728833, 0.1693987,  0.04728833, 0.5270084,    -1.0001074 ],
        [ 0.21621611,   0.18221463, 0.04624761, 0.18221463, 0.04624761, 0.25456479,   -0.99071431 ],
        [ 0.49857725,   0.19505762, 0.04463871, 0.19505762, 0.04463871, -0.023871232, -0.98225242 ],
        [ 0.78750733,   0.2076616,  0.0425412,  0.2076616,  0.0425412,  -0.30653627,  -0.97534839 ],
        [ 1.08479987,   0.21994341, 0.04003529, 0.21994341, 0.04003529, -0.59559216,  -0.97027071 ],
        [ 1.39410554,   0.23189112, 0.03719024, 0.23189112, 0.03719024, -0.89505114,  -0.9670553 ],
        [ 1.71973308,   0.24349503, 0.03406988, 0.24349503, 0.03406988, -1.2095705,   -0.96561993 ],
        [ 2.06620375,   0.254721,   0.03073988, 0.254721,   0.03073988, -1.5440311,   -0.96584086 ],
        [ 2.43943861,   0.26553984, 0.02726046, 0.26553984, 0.02726046, -1.9047196,   -0.96763929 ],
        [ 2.84615352,   0.27589313, 0.0236975,  0.27589313, 0.0236975,  -2.2988289,   -0.97113321 ],
        [ 3.24944253,   0.28478725, 0.0204634,  0.28478725, 0.0204634,  -2.6913268,   -0.97600032 ],
        [ 3.6583446,    0.2925409,  0.01751772, 0.2925409,  0.01751772, -3.0915588,   -0.99140047 ],
        [ 4.07887183,   0.29933263, 0.01484243, 0.29933263, 0.01484243, -3.5139206,   -1.0100246 ],
        [ 4.5156409,    0.30527516, 0.01242887, 0.30527516, 0.01242887, -3.9576367,   -1 ],
        [ 4.9727892,    0.31045009, 0.0102709,  0.31045009, 0.0102709,  -4.414785,    -1 ],
        [ 5.45410897,   0.31492054, 0.00836333, 0.31492054, 0.00836333, -4.8961048,   -1 ],
        [ 5.96348947,   0.31874249, 0.00669917, 0.31874249, 0.00669917, -5.4054853,   -1 ],
        [ 6.50597906,   0.32197373, 0.00526685, 0.32197373, 0.00526685, -5.9479749,   -1 ],
        [ 7.0877743,    0.32467017, 0.00405249, 0.32467017, 0.00405249, -6.5297701,   -1 ],
        [ 7.71561374,   0.32688279, 0.00304189, 0.32688279, 0.00304189, -7.1576095,   -1 ],
        [ 8.39481396,   0.32865628, 0.0022217,  0.32865628, 0.0022217,  -7.8368098,   -1 ],
        [ 9.14481241,   0.33006208, 0.00156433, 0.33006208, 0.00156433, -8.5868082,   -1 ],
        [ 9.97967552,   0.33114165, 0.00105452, 0.33114165, 0.00105452, -9.4216713,   -1 ],
        [ 10.92269053,  0.33194291, 0.00067281, 0.33194291, 0.00067281, -10.364686,   -1 ],
        [ 11.99870649,  0.33250843, 0.00040128, 0.33250843, 0.00040128, -11.440702,   -1 ],
        [ 13.25318819,  0.33288594, 0.00021873, 0.33288594, 0.00021873, -12.695184,   -1 ],
        [ 14.75837313,  0.33311939, 0.00010511, 0.33311939, 0.00010511, -14.200369,   -1 ],
        [ 16.71309072,  0.33325165, 4.034e-05,  0.33325165, 4.034e-05,  -16.155087,   -1 ],
        [ 17.97592944,  0.33328959, 2.166e-05,  0.33328959, 2.166e-05,  -17.417925,   -1 ]
    ],

    impulse_table => [
        [
            -12.673661, 13.815199,  -12.674478,    1.8809556,  0,          0,
            0,          -7.5088769, 0.00059693965, -87.491292, 0.30121118, 0.99965847,
            -8.3713421
        ],
        [
            -2.3025851, 3.4865065,  -2.455828,     1.7084805,  0,          0,
            0,          -2.4443124, 0.00059694008, -87.391295, 0.29058061, 0.78533953,
            -1.9444942
        ],
        [
            -1.6094379, 2.8312218,  -1.8277228,    1.6609959,  0,          0,
            0,          -2.2024863, 0.00059694205, -87.291295, 0.28112065, 0.70621459,
            -1.6307684
        ],
        [
            -1.3862944, 2.6256205,  -1.6306072,    1.646458,   0,          0,
            0,          -2.1387929, 0.00059694379, -87.241295, 0.27675935, 0.6777067,
            -1.5381848
        ],
        [
            -1.2039728, 2.4600396,  -1.47165,      1.6354122,  0,          0,
            0,          -2.0929742, 0.00059694608, -87.191295, 0.27261419, 0.65351581,
            -1.4658368
        ],
        [
            -0.91629073, 2.2037268, -1.2249523,    1.6202211,  0,         0,
            0,           -2.033205, 0.00059695245, -87.091295, 0.2649004, 0.61402104,
            -1.3579424
        ],
        [
            -0.69314718, 2.0095303,  -1.0372836,    1.6108894,  0,          0,
            0,           -1.9981865, 0.00059696131, -86.991295, 0.25785642, 0.58257705,
            -1.2796738
        ],
        [
            -0.51082562, 1.8541055,  -0.88645495,   1.6051594,  0,          0,
            0,           -1.9772768, 0.00059697276, -86.891295, 0.25138597, 0.55657607,
            -1.219294
        ],
        [
            -0.35667494, 1.7250914,  -0.76074378,   1.6017649,  0,          0,
            0,           -1.9650938, 0.00059698686, -86.791295, 0.24541208, 0.5345021,
            -1.1707578
        ],
        [
            -0.22314355, 1.6151686,  -0.65321602,   1.5999414,  0,          0,
            0,           -1.9586112, 0.00059700363, -86.691295, 0.23987214, 0.51539253,
            -1.1305701
        ],
        [
            -9.3259e-10, 1.4353697,  -0.47639777,   1.5992274,  0,          0,
            0,           -1.9560886, 0.00059704529, -86.491295, 0.22989608, 0.48365548,
            -1.0672206
        ],
        [
            0.26236426, 1.2302959,  -0.27312534,   1.6019992,  0,          0,
            0,          -1.9659482, 0.00059712819, -86.191295, 0.21713791, 0.4469697,
            -0.99870581
        ],
        [
            0.53062825, 1.0276765, -0.070421897,  1.6084934,  0,          0,
            0,          -1.989431, 0.00059727692, -85.791295, 0.20317594, 0.41054584,
            -0.93508108
        ],
        [
            0.69314718, 0.90834065, 0.049892847,   1.6139874,  0,          0,
            0,          -2.0097385, 0.00059741712, -85.491295, 0.19439764, 0.38916235,
            -0.89956352
        ],
        [
            0.99325177, 0.69450619, 0.26727824,    1.6266371,  0,          0,
            0,          -2.0581401, 0.00059783968, -84.791295, 0.17786145, 0.35127648,
            -0.83965379
        ],
        [
            1.0986123, 0.62136625, 0.34216942,    1.6317032,  0,          0,
            0,         -2.0782092, 0.00059806168, -84.491295, 0.17203061, 0.33851402,
            -0.8202834
        ],
        [
            1.2089603,    0.5457972,  0.41983518, 1.6372853,  0, 0, 0, -2.1008051,
            0.0005983517, -84.141295, 0.16594624, 0.32546558, -0.80087871
        ],
        [
            1.3862944, 0.42648789, 0.54304252,    1.6467355,  0,          0,
            0,         -2.140284,  0.00059897909, -83.491295, 0.15627147, 0.30519964,
            -0.77151409
        ],
        [
            1.6094379, 0.27990542, 0.69537592,    1.6592291,  0,          0,
            0,         -2.1950489, 0.00060017043, -82.491295, 0.14438873, 0.28096482,
            -0.73758819
        ],
        [
            1.7917595, 0.16288662, 0.81771586,    1.6697313,  0,          0,
            0,         -2.2436034, 0.00060163757, -81.491295, 0.13500794, 0.26222887,
            -0.71222034
        ],
        [
            1.9459101, 0.065749594, 0.91973652,    1.6786947,  0,          0,
            0,         -2.2870812,  0.00060338313, -80.491295, 0.12735252, 0.24713828,
            -0.69232509
        ],
        [
            2.0794415, -0.017139046, 1.0071089,     1.6864536,  0,          0,
            0,         -2.3264004,   0.00060541038, -79.491295, 0.12094645, 0.2346193,
            -0.6761841
        ],
        [
            2.3025851, -0.15322722, 1.1511462,     1.699272,   0,          0,
            0,         -2.3952547,  0.00061032626, -77.491295, 0.11074038, 0.21483012,
            -0.65135614
        ],
        [
            2.7080502, -0.39352909, 1.4070506,     1.7215129,  0,           0,
            0,         -2.5292563,  0.00062790051, -72.491295, 0.093868233, 0.18238195,
            -0.61261328
        ],
        [
            2.9957323, -0.55926429, 1.5845055,     1.7359435,  0,           0,
            0,         -2.630325,   0.00065381927, -67.491295, 0.083211838, 0.16196575,
            -0.58974294
        ],
        [
            3.4011974, -0.78720496, 1.829508,     1.7532274,  0,           0,
            0,         -2.7798827,  0.0007372961, -57.491295, 0.069964225, 0.13658094,
            -0.56375271
        ],
        [
            3.912023,     -1.066633,  2.1307934,   1.7645103,  0, 0, 0, -2.9782192,
            0.0011467852, -37.491295, 0.055970275, 0.10967975, -0.5431713
        ],
        [
            4.2484952, -1.2468711, 2.325377,     1.747182,   0,           0,
            0,         -3.1140779, 0.0029996227, -17.491295, 0.048208221, 0.094675192,
            -0.54722857
        ]
    ],

    blast_info => [ 0, 0, 0, 0, 0, 0, 1.8799196, 0, 0, 0, 0, 0, 0.27795619, 0.30121117, 0, 0, 0, 0, 0, 0 ],

};
1;
