package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='C', gamma=2
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'C2'} = {

    table_name  => 'C2',
    symmetry    => 1,
    gamma       => 2,
    data_source => 'C8000_G2_inspirion/moc2_from_moc_r2000/',

    shock_table_info => [ 1, 2, 8e-07, 8000, 2.49e-07, 42, 1.1e-07, 1.9e-07, 5e-07, 103.51, -1.471 ],

    shock_table => [
        [ -6.37381146, 12.00017718,  -1.99998362, -6.37524347, 0.99856698 ],
        [ -5.77060139, 10.79377453,  -1.99994526, -5.7732206,  0.99737743 ],
        [ -4.96034929, 9.17337044,   -1.99974978, -4.96624785, 0.99408479 ],
        [ -4.46967424, 8.19222894,   -1.9993329,  -4.47932607, 0.99030484 ],
        [ -4.07507385, 7.40342872,   -1.99853308, -4.08942575, 0.98555578 ],
        [ -3.74463522, 6.74323763,   -1.99716516, -3.7646563,  0.97980736 ],
        [ -3.45965849, 6.17437125,   -1.99500372, -3.4863545,  0.97301576 ],
        [ -3.20838628, 5.67345397,   -1.99178145, -3.24281091, 0.96512801 ],
        [ -2.9826444,  5.224307,     -1.98718094, -3.02592384, 0.95607033 ],
        [ -2.77649817, 4.81526941,   -1.98082561, -2.82986025, 0.94574593 ],
        [ -2.58524463, 4.43720029,   -1.97225748, -2.65006828, 0.93401789 ],
        [ -2.40518861, 4.08304867,   -1.9609229,  -2.48305475, 0.92071351 ],
        [ -2.23260064, 3.74583369,   -1.94608479, -2.32542256, 0.90555318 ],
        [ -2.06351818, 3.41834811,   -1.92670091, -2.17374693, 0.88809418 ],
        [ -1.89253214, 3.09100591,   -1.90106934, -2.02361012, 0.86752631 ],
        [ -1.70777316, 2.74293259,   -1.86530196, -1.86564871, 0.84176699 ],
        [ -1.52668829, 2.40900978,   -1.82118257, -1.7157792,  0.81287845 ],
        [ -1.36132391, 2.11174293,   -1.77285823, -1.58375022, 0.78349238 ],
        [ -1.20668624, 1.84151009,   -1.72122265, -1.46486787, 0.75373621 ],
        [ -1.05997278, 1.59287936,   -1.66743768, -1.35645735, 0.7238974 ],
        [ -0.91866016, 1.36110658,   -1.61243437, -1.25625433, 0.69415635 ],
        [ -0.77788946, 1.13809186,   -1.5558668,  -1.16065315, 0.66406334 ],
        [ -0.63738357, 0.92348229,   -1.49898216, -1.06945867, 0.63407316 ],
        [ -0.49593293, 0.71545951,   -1.44252317, -0.98187547, 0.60440931 ],
        [ -0.35208397, 0.51197625,   -1.38699158, -0.89704563, 0.57521492 ],
        [ -0.20455568, 0.31138899,   -1.33284141, -0.81431108, 0.54665171 ],
        [ -0.05204554, 0.11216263,   -1.28042723, -0.73308373, 0.51886692 ],
        [ 0.1067123,   -0.08705862,  -1.2300475,  -0.65287011, 0.49200908 ],
        [ 0.27299538,  -0.28752952,  -1.1819434,  -0.57323479, 0.46621988 ],
        [ 0.4481245,   -0.49045241,  -1.136303,   -0.49377728, 0.44162911 ],
        [ 0.63375735,  -0.69730696,  -1.09320848, -0.41400246, 0.41831921 ],
        [ 0.83143647,  -0.90932665,  -1.05276873, -0.33352382, 0.39639132 ],
        [ 1.04319268,  -1.12816295,  -1.01499761, -0.25180713, 0.37589537 ],
        [ 1.27054667,  -1.35484777,  -0.98000603, -0.16855728, 0.35693136 ],
        [ 1.50565985,  -1.58152069,  -0.94900515, -0.08665922, 0.34018706 ],
        [ 1.75063856,  -1.81054033,  -0.92144436, -0.00518453, 0.32538373 ],
        [ 2.00727252,  -2.04378071,  -0.89692172, 0.07659433,  0.31231601 ],
        [ 2.27730161,  -2.28294776,  -0.87511666, 0.15932895,  0.30081748 ],
        [ 2.56225543,  -2.52947936,  -0.85578009, 0.24356871,  0.29075626 ],
        [ 2.86380761,  -2.78488633,  -0.83869135, 0.32988415,  0.28201283 ],
        [ 3.18362655,  -3.05063593,  -0.82366252, 0.41883002,  0.27448257 ],
        [ 3.52347423,  -3.32824919,  -0.81052572, 0.510981,    0.26806931 ],
        [ 3.88518448,  -3.61929164,  -0.79912959, 0.60693046,  0.26268343 ],
        [ 4.27068258,  -3.92539705,  -0.78933423, 0.70729972,  0.2582391 ],
        [ 4.68206754,  -4.24833743,  -0.78100591, 0.81276151,  0.25465173 ],
        [ 5.12165964,  -4.59006107,  -0.77401482, 0.92405206,  0.25183692 ],
        [ 5.5921672,   -4.95282,     -0.76823211, 1.04201172,  0.24970942 ],
        [ 6.09708645,  -5.33947057,  -0.7635276,  1.16768132,  0.24818289 ],
        [ 6.64029991,  -5.75315717,  -0.75977622, 1.30219839,  0.24717279 ],
        [ 7.22727849,  -6.19822218,  -0.75685089, 1.44709254,  0.24659492 ],
        [ 7.86388025,  -6.67928454,  -0.75463261, 1.60398576,  0.24636875 ],
        [ 8.56215412,  -7.20561347,  -0.7529959,  1.77602199,  0.24641748 ],
        [ 9.33513707,  -7.78718232,  -0.75183458, 1.96658752,  0.24667208 ],
        [ 10.20265827, -8.43904298,  -0.75104825, 2.18074863,  0.24706984 ],
        [ 11.19554006, -9.18446982,  -0.75054642, 2.42630178,  0.24755599 ],
        [ 12.36259929, -10.06020774, -0.75025079, 2.71552915,  0.24808325 ],
        [ 13.78884905, -11.13012455, -0.75009529, 3.06975071,  0.24861153 ],
        [ 15.64150089, -12.51971326, -0.75002655, 3.53083053,  0.24910664 ],
        [ 18.3259696,  -14.53309695, -0.75000403, 4.20018889,  0.24953758 ],
        [ 23.34909541, -18.3004468,  -0.75000011, 5.45464354,  0.24986758 ]
    ],

    energy_table => [
        [ -6.37381146, 0.0024492,  0.00244011, 0.0024492,  0.00244011,    6.5809005,   -1.1792938 ],
        [ -5.77060139, 0.00446334, 0.00443304, 0.00446334, 0.00443304,    5.8639674,   -1.2008854 ],
        [ -4.96034929, 0.00995104, 0.00979928, 0.00995104, 0.00979928,    4.8775001,   -1.2397237 ],
        [ -4.46967424, 0.01609806, 0.0156976,  0.01609806, 0.0156976,     4.2625896,   -1.2704483 ],
        [ -4.07507385, 0.02360068, 0.02273121, 0.02360068, 0.02273121,    3.7557945,   -1.3016433 ],
        [ -3.74463522, 0.03237239, 0.030717,   0.03237239, 0.030717,      3.3208891,   -1.3333037 ],
        [ -3.45965849, 0.04232802, 0.03945944, 0.04232802, 0.03945944,    2.9367117,   -1.3653682 ],
        [ -3.20838628, 0.05337823, 0.04874743, 0.05337823, 0.04874743,    2.589807,    -1.397788 ],
        [ -2.9826444,  0.06544652, 0.05836986, 0.06544652, 0.05836986,    2.2707805,   -1.4305766 ],
        [ -2.77649817, 0.07846947, 0.06811444, 0.07846947, 0.06811444,    1.9726081,   -1.4627477 ],
        [ -2.58524463, 0.09241314, 0.07777641, 0.09241314, 0.07777641,    1.6899539,   -1.4931947 ],
        [ -2.40518861, 0.10725902, 0.08714169, 0.10725902, 0.08714169,    1.418503,    -1.5204156 ],
        [ -2.23260064, 0.12306869, 0.09601501, 0.12306869, 0.09601501,    1.153976,    -1.5424952 ],
        [ -2.06351818, 0.14000514, 0.10419598, 0.14000514, 0.10419598,    0.89154596,  -1.5567165 ],
        [ -1.89253214, 0.15846101, 0.11146956, 0.15846101, 0.11146956,    0.62456804,  -1.5581919 ],
        [ -1.70777316, 0.17965497, 0.11760371, 0.17965497, 0.11760371,    0.33731893,  -1.540167 ],
        [ -1.52668829, 0.20132517, 0.12131289, 0.20132517, 0.12131289,    0.061001525, -1.5022222 ],
        [ -1.36132391, 0.22150473, 0.1223589,  0.22150473, 0.1223589,     -0.18383696, -1.4535704 ],
        [ -1.20668624, 0.24036304, 0.12120142, 0.24036304, 0.12120142,    -0.40470504, -1.4035783 ],
        [ -1.05997278, 0.25795035, 0.118264,   0.25795035, 0.118264,      -0.60718853, -1.3612909 ],
        [ -0.91866016, 0.27437232, 0.11393146, 0.27437232, 0.11393146,    -0.79699206, -1.3312411 ],
        [ -0.77788946, 0.29003483, 0.10842049, 0.29003483, 0.10842049,    -0.98272233, -1.3139806 ],
        [ -0.63738357, 0.3048291,  0.10204748, 0.3048291,  0.10204748,    -1.1665869,  -1.3086795 ],
        [ -0.49593293, 0.31877574, 0.09508183, 0.31877574, 0.09508183,    -1.3517132,  -1.3128723 ],
        [ -0.35208397, 0.33192708, 0.08775037, 0.33192708, 0.08775037,    -1.5411687,  -1.3237808 ],
        [ -0.20455568, 0.34431819, 0.08025758, 0.34431819, 0.08025758,    -1.737483,   -1.3388829 ],
        [ -0.05204554, 0.35598327, 0.07277847, 0.35598327, 0.07277847,    -1.942968,   -1.3561469 ],
        [ 0.1067123,   0.3669495,  0.06546319, 0.3669495,  0.06546319,    -2.1597203,  -1.3740903 ],
        [ 0.27299538,  0.37724112, 0.05843556, 0.37724112, 0.05843556,    -2.389739,   -1.3916676 ],
        [ 0.4481245,   0.38688172, 0.05179253, 0.38688172, 0.05179253,    -2.6350046,  -1.4081936 ],
        [ 0.63375735,  0.39590785, 0.04559634, 0.39590785, 0.04559634,    -2.8979285,  -1.4232878 ],
        [ 0.83143647,  0.4043432,  0.03989442, 0.4043432,  0.03989442,    -3.18074,    -1.436756 ],
        [ 1.04319268,  0.41222587, 0.03470364, 0.41222587, 0.03470364,    -3.486365,   -1.4485481 ],
        [ 1.27054667,  0.41956926, 0.03003895, 0.41956926, 0.03003895,    -3.8169825,  -1.4586826 ],
        [ 1.50565985,  0.4261469,  0.02604059, 0.4261469,  0.02604059,    -4.1610286,  -1.4669631 ],
        [ 1.75063856,  0.43209157, 0.02260274, 0.43209157, 0.02260274,    -4.5213329,  -1.4737087 ],
        [ 2.00727252,  0.43749972, 0.01964127, 0.43749972, 0.01964127,    -4.90033,    -1.4791912 ],
        [ 2.27730161,  0.44244693, 0.01708559, 0.44244693, 0.01708559,    -5.300434,   -1.4836356 ],
        [ 2.56225543,  0.44699042, 0.01487747, 0.44699042, 0.01487747,    -5.7237815,  -1.4872338 ],
        [ 2.86380761,  0.45117888, 0.01296593, 0.45117888, 0.01296593,    -6.1727589,  -1.4901301 ],
        [ 3.18362655,  0.45505159, 0.01130763, 0.45505159, 0.01130763,    -6.6497521,  -1.4924404 ],
        [ 3.52347423,  0.45864122, 0.0098653,  0.45864122, 0.0098653,     -7.1573131,  -1.4942777 ],
        [ 3.88518448,  0.46197454, 0.00860709, 0.46197454, 0.00860709,    -7.6981097,  -1.4957173 ],
        [ 4.27068258,  0.46507338, 0.00750587, 0.46507338, 0.00750587,    -8.2749554,  -1.4968461 ],
        [ 4.68206754,  0.46795578, 0.00653846, 0.46795578, 0.00653846,    -8.8909468,  -1.4977175 ],
        [ 5.12165964,  0.47063654, 0.00568517, 0.47063654, 0.00568517,    -9.5494995,  -1.4983735 ],
        [ 5.5921672,   0.47312807, 0.00492933, 0.47324976, 0.0054838119,  -10.254634,  -1.4988706 ],
        [ 6.09708645,  0.47544188, 0.00425661, 0.47592032, 0.0050594927,  -11.011553,  -1.4992327 ],
        [ 6.64029991,  0.47758585, 0.00365564, 0.47855416, 0.0046528874,  -11.82604,   -1.49949 ],
        [ 7.22727849,  0.4795685,  0.00311659, 0.48112124, 0.0041240502,  -12.706274,  -1.4996721 ],
        [ 7.86388025,  0.48139354, 0.00263246, 0.48358222, 0.0036239734,  -13.661016,  -1.4997969 ],
        [ 8.56215412,  0.48307399, 0.00219517, 0.48594082, 0.0031356662,  -14.708322,  -1.4998802 ],
        [ 9.33513707,  0.48461289, 0.00180039, 0.48818133, 0.0026703649,  -15.86773,   -1.4999333 ],
        [ 10.20265827, 0.48601455, 0.00144458, 0.49030108, 0.0022205868,  -17.168973,  -1.4999677 ],
        [ 11.19554006, 0.48728329, 0.00112405, 0.49225888, 0.0017366266,  -18.658278,  -1.4999898 ],
        [ 12.36259929, 0.48842097, 0.00083901, 0.49407303, 0.0013871205,  -20.408865,  -1.5000056 ],
        [ 13.78884905, 0.48942714, 0.00058697, 0.49582762, 0.0010434915,  -22.54826,   -1.5000252 ],
        [ 15.64150089, 0.49029738, 0.00036926, 0.49737463, 0.00065641279, -25.327311,  -1.500079 ],
        [ 18.3259696,  0.49101941, 0.00018873, 0.49865812, 0.00033547601, -29.35438,   -1.5000903 ],
        [ 23.34909541, 0.49155928, 5.376e-05,  0.49961776, 9.5559378e-05, -36.88909,   -1.5 ]
    ],

    impulse_table => [
        [
            -7.281258, 13.815063,  -7.2818356,    0.03749185,  0,          -0.21070872,
            0,         -2.2107159, -0.0043007923, -0.50309314, 0.31068878, 0.99801856,
            -4.0232116
        ],
        [
            -6.9077549, 13.068059,  -6.9085943,    0.042600456, 0,          -0.21039704,
            0,          -2.2107192, -0.0051838642, -0.50278146, 0.31068844, 0.99712347,
            -3.8443865
        ],
        [
            -6.5022898, 12.257132,  -6.5035491,    0.048731845, 0,          -0.20989704,
            0,          -2.2107257, -0.0063489111, -0.50228146, 0.31068768, 0.99569005,
            -3.6522087
        ],
        [
            -6.2146078, 11.681773,  -6.2162867,    0.053450119, 0,          -0.20939704,
            0,          -2.2107345, -0.0073310911, -0.50178146, 0.31068666, 0.99425974,
            -3.5173527
        ],
        [
            -5.8091427, 10.870855,  -5.8116627,    0.060594174, 0,          -0.20839704,
            0,          -2.2107594, -0.0089787162, -0.50078146, 0.31068379, 0.99140846,
            -3.3298739
        ],
        [
            -5.298317, 9.849244,   -5.3025162,   0.070309298, 0,          -0.20639705,
            0,         -2.2108388, -0.011591473, -0.49878146, 0.31067478, 0.9857431,
            -3.0991266
        ],
        [
            -4.6051699, 8.4631407,  -4.613593,    0.084247739, 0,          -0.20139715,
            0,          -2.2112113, -0.016392816, -0.49378145, 0.31063259, 0.97179461,
            -2.7994622
        ],
        [
            -3.9120227, 7.0776085,  -3.9289347,   0.097721853, 0,          -0.19139881,
            0,          -2.2127011, -0.023182908, -0.48378145, 0.31046391, 0.94479872,
            -2.5229149
        ],
        [
            -3.5065576, 6.2679455, -3.5320177,   0.10442788,  0,          -0.18140603,
            0,          -2.215183, -0.028392964, -0.47378145, 0.31018354, 0.9189652,
            -2.3769944
        ],
        [
            -2.9957319, 5.2503164,  -3.0384403,   0.11034963,  0,          -0.16146645,
            0,          -2.2231135, -0.036653163, -0.45378144, 0.30929271, 0.87059942,
            -2.2173734
        ],
        [
            -2.6592597, 4.583316,   -2.7193797,   0.11205728,  0,          -0.14166337,
            0,          -2.2349716, -0.043360909, -0.43394826, 0.30797432, 0.8263325,
            -2.1322145
        ],
        [
            -2.3025848, 3.8822622, -2.3890251,   0.11168637,  0,          -0.11250022,
            0,          -2.260001, -0.051786999, -0.40461605, 0.30524107, 0.76685609,
            -2.0651028
        ],
        [
            -2.0402205, 3.3734963, -2.1530868,   0.11023957,  0,         -0.084506352,
            0,          -2.293453, -0.058938397, -0.37645641, 0.3016822, 0.71471612,
            -2.0346887
        ],
        [
            -1.8971197,  3.0997289,   -2.0275913, 0.10926848, 0, -0.06682498, 0, -2.3202186,
            -0.06317476, -0.35838751, 0.29890316, 0.68353687, -2.0260199
        ],
        [
            -1.6094376, 2.5606063,  -1.783614,    0.10776205,  0,          -0.027416286,
            0,          -2.4013721, -0.072216373, -0.31770731, 0.29078101, 0.61622403,
            -2.0276066
        ],
        [
            -1.386294, 2.1561089,  -1.6033718,   0.10786245,  0,          0.0038168626,
            0,         -2.4995666, -0.079230171, -0.28538289, 0.28140482, 0.56136559,
            -2.047702
        ],
        [
            -1.2039725, 1.8368404,  -1.4628231,   0.1092112,   0,         0.027061262,
            0,          -2.6096798, -0.084362342, -0.26232356, 0.2712591, 0.51619144,
            -2.0767987
        ],
        [
            -1.0498218, 1.5759728, -1.3491198,   0.1112793,   0,          0.043884663,
            0,          -2.726314, -0.087875466, -0.24727858, 0.26074687, 0.47857878,
            -2.1101632
        ],
        [
            -0.9162904, 1.3572866,  -1.2546099,   0.11365622,  0,          0.05606552,
            0,          -2.8448598, -0.090139768, -0.23866441, 0.25018171, 0.44691355,
            -2.1452622
        ],
        [
            -0.79850737, 1.1702567,  -1.1743903,   0.11608688,  0,          0.065030582,
            0,           -2.9620238, -0.091520894, -0.23414894, 0.23979508, 0.4199695,
            -2.1806939
        ],
        [
            -0.69314685, 1.007699,   -1.1051473,   0.11843286,  0,          0.071783515,
            0,           -3.0757765, -0.092309633, -0.23250923, 0.22974952, 0.39681138,
            -2.2156767
        ],
        [
            -0.59783667, 0.86451739, -1.0445487,   0.12062863, 0,          0.076997573,
            0,           -3.1850351, -0.092711683, -0.2325786, 0.22015281, 0.3767209,
            -2.2497821
        ],
        [
            -0.51082529, 0.73698585, -0.99089953, 0.12265007, 0, 0.081120127, 0, -3.2893367,
            -0.0928644,  -0.2336196, 0.21107061,  0.35914164, -2.2827873
        ],
        [
            -0.35667461, 0.51834742, -0.89968833,  0.12617231,  0,          0.087194545,
            0,           -3.4829147, -0.092748794, -0.23792559, 0.19456326, 0.32986783,
            -2.3451632
        ],
        [
            -0.22314322, 0.33622494, -0.82450461,  0.12907857,  0,          0.091447165,
            0,           -3.6577403, -0.092365147, -0.24328476, 0.18026141, 0.30648161,
            -2.4026883
        ],
        [
            -0.10536018, 0.18090292, -0.76099915,  0.13148425,  0,          0.094596708,
            0,           -3.8160538, -0.091885439, -0.24885758, 0.16800567, 0.28736334,
            -2.4556904
        ],
        [
            3.3099852e-07, 0.04596597, -0.70631519, 0.13349294,  0,          0.097031978,
            0,             -3.9601235, -0.09138493, -0.25479553, 0.15755464, 0.27142922,
            -2.5046208
        ],
        [
            0.095310511, -0.073014098, -0.65849047,  0.13518757,  0,          0.098978993,
            0,           -4.0919582,   -0.090896631, -0.26076234, 0.14864164, 0.25793033,
            -2.5499347
        ],
        [
            0.18232189, -0.17920668, -0.61612722,  0.13663225,  0,          0.10057719,
            0,          -4.2132577,  -0.090434363, -0.26613835, 0.14099258, 0.24633432,
            -2.5920493
        ],
        [
            0.23696096, -0.24475898, -0.59013143,  0.13749006,  0,          0.10150432,
            0,          -4.2897827,  -0.090140209, -0.26983732, 0.13642902, 0.2393934,
            -2.6188133
        ],
        [
            0.33647257, -0.362009,  -0.54393468,  0.13895697,  0,          0.10306007,
            0,          -4.4296748, -0.089603355, -0.27666525, 0.12858334, 0.22739764,
            -2.6681033
        ],
        [
            0.40546544, -0.4417526, -0.51273888,  0.13990407,  0,          0.1040491,
            0,          -4.526948,  -0.089234049, -0.28155098, 0.12348544, 0.21954817,
            -2.7026366
        ],
        [
            0.53062858, -0.58337827, -0.45778617,  0.14148381,  0,         0.10568161,
            0,          -4.7037519,  -0.088578007, -0.29123068, 0.1149116, 0.20622446,
            -2.7659049
        ],
        [
            0.69314751, -0.761853,  -0.38936422,  0.143289,    0,          0.10753815,
            0,          -4.9334915, -0.087765741, -0.30520922, 0.10497197, 0.19054738,
            -2.8489837
        ],
        [
            0.9932521,    -1.0772643, -0.27069301, 0.14599417, 0, 0.11035852, 0, -5.3565426,
            -0.086418367, -0.3326221, 0.089642321, 0.16577057, -3.0041286
        ],
        [
            1.0986126, -1.1841626, -0.23111144,  0.14677973,  0,           0.11119973,
            0,         -5.5043654, -0.085996742, -0.34356075, 0.085047279, 0.15817941,
            -3.0588993
        ],
        [
            1.2089607, -1.2942195, -0.19068739,  0.14752445,  0,           0.11201217,
            0,         -5.6586839, -0.085584062, -0.35638267, 0.080608326, 0.15076597,
            -3.1163443
        ],
        [
            1.3862947, -1.4673518, -0.12774605,  0.14857335,  0,           0.11318934,
            0,         -5.9054967, -0.084981434, -0.37643159, 0.074184114, 0.13988851,
            -3.2087328
        ],
        [
            1.6094382, -1.6793661, -0.051700118, 0.1496728,   0,           0.11447882,
            0,         -6.2138534, -0.084323291, -0.40585572, 0.067162804, 0.12778529,
            -3.3249407
        ],
        [
            1.7917598, -1.8483445, 0.0081492407, 0.15041986,  0,           0.11539916,
            0,         -6.4639181, -0.083861261, -0.43211222, 0.062161353, 0.11901819,
            -3.4197379
        ],
        [
            1.9459105, -1.9885764, 0.057341161,  0.15096294,  0,           0.1160975,
            0,         -6.6740338, -0.083518761, -0.45639282, 0.058368448, 0.11228478,
            -3.4997343
        ],
        [
            2.0794419,   -2.1082859,  0.099014747, 0.15137692, 0, 0.1166502, 0, -6.8550957,
            -0.08325416, -0.47803703, 0.055363898, 0.10689736, -3.5688978
        ],
        [
            2.3025854, -2.3050506, 0.16692246,   0.15196921,  0,           0.11747912,
            0,         -7.1557769, -0.082871613, -0.52015624, 0.050846501, 0.098705986,
            -3.6841813
        ],
        [
            2.7080505, -2.6536082, 0.28563027,   0.15279272,  0,           0.11873157,
            0,         -7.6965283, -0.082336025, -0.60652095, 0.043950621, 0.085985754,
            -3.8926977
        ],
        [
            2.9957326, -2.8950934, 0.36686788,   0.15322743,  0,          0.11945968,
            0,         -8.0762775, -0.082055541, -0.68007938, 0.03986634, 0.078327882,
            -4.0399236
        ],
        [
            3.4011977, -3.2288739, 0.4780736,    0.15368833,  0,           0.12030751,
            0,         -8.6067701, -0.081764039, -0.80563926, 0.034979632, 0.069046373,
            -4.2465359
        ],
        [
            3.7519938, -3.5125976, 0.57182364,   0.15398113,  0,           0.12090524,
            0,         -9.0619311, -0.081581373, -0.93347045, 0.031399619, 0.061492025,
            -4.4341503
        ],
        [
            3.9120478, -3.6407488, 0.61398233,   0.15408957, 0,           0.12114234,
            0,         -9.2685941, -0.081515432, -1.0002912, 0.029929972, 0.058756455,
            -4.5140733
        ],
        [
            4.605245,    -4.188285,  0.79317609,  0.15442834,  0, 0.12196779, 0, -10.157693,
            -0.08131696, -1.3598829, 0.024504386, 0.048474421, -4.8604894
        ],
        [
            5.2984903, -4.7267192, 0.968504,     0.15462508, 0,           0.12254566,
            0,         -11.039431, -0.081208707, -1.867635,  0.020246934, 0.040232351,
            -5.2071204
        ],
        [
            6.2146916, -5.4292112, 1.1968531,    0.15476604, 0,           0.12307176,
            0,         -12.197242, -0.080828872, -2.8349033, 0.015874671, 0.031646537,
            -5.6652906
        ],
        [
            6.9077947, -5.9561939, 1.3682713,    0.15482319, 0,           0.12336122,
            0,         -13.069442, -0.080160879, -3.9106981, 0.013262692, 0.026474392,
            -6.0118878
        ],
        [
            7.6011365, -6.4809056, 1.5392481,    0.15485693, 0,           0.12352474,
            0,         -13.939923, -0.078967607, -5.3820961, 0.011105267, 0.022185178,
            -6.3585917
        ],
        [
            8.5173403, -7.171867,  1.7649793,   0.15488128, 0,            0.12368533,
            0,         -15.088234, -0.07680636, -8.2084489, 0.0088024148, 0.017594729,
            -6.816721
        ],
        [
            9.2105136, -7.6934769, 1.9358496,    0.15489109, 0,            0.123766,
            0,         -15.956078, -0.073962439, -11.133949, 0.0073908191, 0.014776598,
            -7.16332
        ],
        [
            10.819914, -8.9025177, 2.3333476,    0.15490042, 0,            0.12387339,
            0,         -17.969399, -0.059336706, -19.870022, 0.0049342358, 0.0098674672,
            -7.9680352
        ],
        [
            11.512959, -9.4226907, 2.5049048,   0.15490147, 0,            0.12389878,
            0,         -18.83601,  -0.05898569, -27.947706, 0.0041479382, 0.0082953758,
            -8.3145615
        ],
        [
            13.122428, -10.630229, 2.9041456,    0.15490019, 0,            0.12393265,
            0,         -20.848185, -0.081067632, -86.189431, 0.0027728719, 0.0055456482,
            -9.119305
        ],
        [
            13.81565, -11.150228, 3.0764138,    0.15489817, 0,            0.12394067,
            0,        -21.714777, -0.013302414, -19.835258, 0.0023314939, 0.0046629446,
            -9.4659215
        ],
        [
            16.118113, -12.877184, 3.6495812,    0.15488134, 0,           0.1239539,
            0,         -24.592941, -0.081066699, -384.76012, 0.001311009, 0.0026220349,
            -10.617191
        ],
        [
            18.420772, -14.604199, 4.2238461,    0.15482697, 0,             0.12395796,
            0,         -27.471289, -0.081066521, -1216.2939, 0.00073717393, 0.0014744167,
            -11.768638
        ]
    ],

    tail_shock_table => [
        [ 4.6537794, -1.471,     -0.030241375,  -0.030241375 ],
        [ 5.2852351, -1.8143216, -0.040558728,  -0.014856056 ],
        [ 6.2080507, -2.4211798, -0.036213647,  -0.0078029429 ],
        [ 6.9038586, -2.9824634, -0.031909802,  -0.0052705925 ],
        [ 7.5988038, -3.6531598, -0.027732103,  -0.0037171816 ],
        [ 8.5161715, -4.7449438, -0.022777889,  -0.0024318603 ],
        [ 9.2098204, -5.7603994, -0.019523201,  -0.001793828 ],
        [ 10.819708, -8.9491052, -0.013504476,  -0.00091359382 ],
        [ 11.512836, -10.81781,  -0.01152298,   -0.00067744725 ],
        [ 13.122391, -16.535918, -0.0078509431, -0.00036302874 ],
        [ 13.815628, -19.852529, -0.006657135,  -0.00027681128 ]
    ],

    blast_info => [
        0.211397,   0,          0.50379894, -0.16392816, 0,          0,          1.3247718,   -0.22118894,
        0.71299196, 0.14739202, 0.82394895, -0.12941628, 0.27590941, 0.31068898, 0.077440671, 0,
        103.51,     -1.471,     167.95147,  -1.7219333
    ],

};
1;
