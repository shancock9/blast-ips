package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='P', gamma=5.5
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'P5x5'} = {

    table_name  => 'P5x5',
    symmetry    => 0,
    gamma       => 5.5,
    data_source => 'P4000_G5x5/moc3p_restart_r100/',

    shock_table_info => [ 0, 5.5, 1.63e-06, 4000, 1.27e-06, 100, 8.55e-07, 2.8e-07, 5e-07 ],

    shock_table => [
        [ -10.07047387, 12.00046259, -0.99998961, -10.0726254, 0.99892309 ],
        [ -9.10889933,  11.03890298, -0.99997282, -9.11238138, 0.99825599 ],
        [ -8.56969281,  10.49971342, -0.99995908, -8.57425478, 0.99771392 ],
        [ -7.3750661,   9.30518189,  -0.99986327, -7.38337125, 0.99583086 ],
        [ -6.40536909,  8.33570829,  -0.99964028, -6.4188893,  0.99319712 ],
        [ -5.62748716,  7.55824935,  -0.99921821, -5.64749432, 0.98990576 ],
        [ -4.97644059,  6.90791944,  -0.99850386, -5.00424061, 0.98593192 ],
        [ -4.4147411,   6.34734572,  -0.9973847,  -4.45169682, 0.98123941 ],
        [ -3.91802131,  5.85230437,  -0.99572279, -3.96559492, 0.97577167 ],
        [ -3.47066109,  5.40734898,  -0.99335546, -3.53042851, 0.96946746 ],
        [ -3.06203878,  5.00206453,  -0.99009627, -3.1357018,  0.96226413 ],
        [ -2.68278228,  4.62734621,  -0.98571408, -2.77226202, 0.95405834 ],
        [ -2.32539395,  4.27604224,  -0.97992583, -2.4329114,  0.944717 ],
        [ -1.98268539,  3.9414456,   -0.97235988, -2.11092659, 0.934041 ],
        [ -1.64716523,  3.61677767,  -0.96249946, -1.7995482,  0.92172742 ],
        [ -1.3075693,   3.29203441,  -0.94946832, -1.48894113, 0.90718358 ],
        [ -0.93905592,  2.94536482,  -0.93122053, -1.15792404, 0.88886886 ],
        [ -0.58068266,  2.61549825,  -0.90892076, -0.84295292, 0.86850075 ],
        [ -0.25295769,  2.32151789,  -0.88452497, -0.56166358, 0.84781007 ],
        [ 0.05373287,   2.05415776,  -0.85852275, -0.30481632, 0.82693608 ],
        [ 0.3373405,    1.8143549,   -0.83226488, -0.07314758, 0.80666726 ],
        [ 0.61060292,   1.59055026,  -0.80558682, 0.14455218,  0.78661656 ],
        [ 0.88246269,   1.37523299,  -0.77839587, 0.35566982,  0.76652636 ],
        [ 1.16109243,   1.16223318,  -0.75058119, 0.56640013,  0.74616022 ],
        [ 1.44752489,   0.95125426,  -0.72275109, 0.77719323,  0.72581994 ],
        [ 1.74127138,   0.74297321,  -0.69563293, 0.98744869,  0.70590577 ],
        [ 2.04650338,   0.5346943,   -0.66946397, 1.199914,    0.68648259 ],
        [ 2.36693965,   0.32424284,  -0.64451327, 1.41682802,  0.66766174 ],
        [ 2.70683344,   0.10925832,  -0.62100479, 1.64062876,  0.64954295 ],
        [ 3.07133331,   -0.11300682, -0.59911557, 1.87416327,  0.63220916 ],
        [ 3.46742635,   -0.34620163, -0.57895892, 2.12123149,  0.61570661 ],
        [ 3.86980062,   -0.57562662, -0.56193924, 2.36598778,  0.60121432 ],
        [ 4.27598631,   -0.80091391, -0.54781493, 2.60757396,  0.58864324 ],
        [ 4.6939307,    -1.02731243, -0.53599508, 2.85122137,  0.57758233 ],
        [ 5.12223785,   -1.25471964, -0.5262569,  3.09648149,  0.56793259 ],
        [ 5.57477132,   -1.49094996, -0.51810731, 3.35148117,  0.55929887 ],
        [ 6.04808347,   -1.7345396,  -0.51147877, 3.61435294,  0.55169771 ],
        [ 6.55135651,   -1.99053891, -0.50611456, 3.89025377,  0.54493264 ],
        [ 7.08278325,   -2.25832773, -0.50191586, 4.17821687,  0.53898733 ],
        [ 7.65548092,   -2.54479089, -0.49867558, 4.48532718,  0.53368783 ],
        [ 8.26437678,   -2.84766621, -0.49631858, 4.80883077,  0.52905598 ],
        [ 8.92981497,   -3.17734244, -0.4946718,  5.15946221,  0.52492308 ],
        [ 9.65599747,   -3.53615315, -0.49365005, 5.53927912,  0.52127262 ],
        [ 10.45663841,  -3.93115522, -0.49315287, 5.95529008,  0.51804304 ],
        [ 11.34781806,  -4.37058614, -0.49309217, 6.41563636,  0.51518504 ],
        [ 12.35501198,  -4.86735139, -0.49338907, 6.93319492,  0.51264215 ],
        [ 13.51466689,  -5.4398357,  -0.49397438, 7.52630456,  0.51036206 ],
        [ 14.88580733,  -6.11769533, -0.49478682, 8.2245982,   0.50829035 ],
        [ 15.83854747,  -6.5893702,  -0.495358,   8.70830172,  0.50713713 ]
    ],

    energy_table => [
        [ -10.07047387, 6.111e-05,  4.829e-05,  6.111e-05,  4.829e-05,  11.45123,    -0.9963711 ],
        [ -9.10889933,  0.00013027, 0.00010211, 0.00013027, 0.00010211, 10.49206,    -0.99864547 ],
        [ -8.56969281,  0.00019857, 0.0001548,  0.00019857, 0.0001548,  9.9532373,   -0.99994825 ],
        [ -7.3750661,   0.00050056, 0.00038426, 0.00050056, 0.00038426, 8.7569251,   -1.0029217 ],
        [ -6.40536909,  0.00104742, 0.00079039, 0.00104742, 0.00079039, 7.7832057,   -1.0054253 ],
        [ -5.62748716,  0.00187434, 0.00138895, 0.00187434, 0.00138895, 7.0003064,   -1.0075142 ],
        [ -4.97644059,  0.00302295, 0.00219706, 0.00302295, 0.00219706, 6.3437879,   -1.0093271 ],
        [ -4.4147411,   0.00452848, 0.00322316, 0.00452848, 0.00322316, 5.7764024,   -1.0109447 ],
        [ -3.91802131,  0.00642421, 0.00446974, 0.00642421, 0.00446974, 5.2738852,   -1.0124208 ],
        [ -3.47066109,  0.00873745, 0.00593013, 0.00873745, 0.00593013, 4.8206665,   -1.0137902 ],
        [ -3.06203878,  0.01148819, 0.00758762, 0.01148819, 0.00758762, 4.4061499,   -1.0150771 ],
        [ -2.68278228,  0.01470492, 0.00942429, 0.01470492, 0.00942429, 4.0209457,   -1.0163045 ],
        [ -2.32539395,  0.01842185, 0.01141727, 0.01842185, 0.01141727, 3.6575209,   -1.0174924 ],
        [ -1.98268539,  0.02269351, 0.01354321, 0.02269351, 0.01354321, 3.3086198,   -1.0186619 ],
        [ -1.64716523,  0.02760923, 0.01577858, 0.02760923, 0.01577858, 2.9666435,   -1.0198026 ],
        [ -1.3075693,   0.03336365, 0.01811523, 0.03336365, 0.01811523, 2.62013,     -1.0204414 ],
        [ -0.93905592,  0.04050174, 0.02060308, 0.04050174, 0.02060308, 2.2440552,   -1.0199231 ],
        [ -0.58068266,  0.04828977, 0.02281052, 0.04828977, 0.02281052, 1.8787496,   -1.0180088 ],
        [ -0.25295769,  0.0560524,  0.02449852, 0.0560524,  0.02449852, 1.5455228,   -1.0149134 ],
        [ 0.05373287,   0.06375859, 0.02568503, 0.06375859, 0.02568503, 1.2347964,   -1.0109436 ],
        [ 0.3373405,    0.07115159, 0.02638329, 0.07115159, 0.02638329, 0.94866622,  -1.0065985 ],
        [ 0.61060292,   0.0784089,  0.02666902, 0.0784089,  0.02666902, 0.67420524,  -1.0021136 ],
        [ 0.88246269,   0.08565526, 0.02658009, 0.08565526, 0.02658009, 0.40238548,  -0.99770746 ],
        [ 1.16109243,   0.0930066,  0.02613066, 0.0930066,  0.02613066, 0.12500716,  -0.99356874 ],
        [ 1.44752489,   0.10038528, 0.02534006, 0.10038528, 0.02534006, -0.15901196, -0.98997028 ],
        [ 1.74127138,   0.1076747,  0.02424931, 0.1076747,  0.02424931, -0.449329,   -0.98715831 ],
        [ 2.04650338,   0.11487415, 0.02289299, 0.11487415, 0.02289299, -0.75027123, -0.985288 ],
        [ 2.36693965,   0.12195935, 0.02130864, 0.12195935, 0.02130864, -1.0657719,  -0.98458691 ],
        [ 2.70683344,   0.12890214, 0.01953468, 0.12890214, 0.01953468, -1.4004231,  -0.98542466 ],
        [ 3.07133331,   0.13567129, 0.01760977, 0.13567129, 0.01760977, -1.7599421,  -0.98852261 ],
        [ 3.46742635,   0.14223959, 0.01557024, 0.14223959, 0.01557024, -2.1524307,  -1.0004165 ],
        [ 3.86980062,   0.14810639, 0.01361381, 0.14810639, 0.01361381, -2.5588624,  -1.048802 ],
        [ 4.27598631,   0.15326094, 0.01179458, 0.15326094, 0.01179458, -3.0007465,  -1 ],
        [ 4.6939307,    0.15783081, 0.0101054,  0.15783081, 0.0101054,  -3.4186909,  -1 ],
        [ 5.12223785,   0.16182347, 0.00857159, 0.16182347, 0.00857159, -3.846998,   -1 ],
        [ 5.57477132,   0.16537542, 0.0071615,  0.16537542, 0.0071615,  -4.2995315,  -1 ],
        [ 6.04808347,   0.16845871, 0.00590206, 0.16845871, 0.00590206, -4.7728436,  -1 ],
        [ 6.55135651,   0.17113791, 0.00478011, 0.17113791, 0.00478011, -5.2761167,  -1 ],
        [ 7.08278325,   0.17341078, 0.00380728, 0.17341078, 0.00380728, -5.8075434,  -1 ],
        [ 7.65548092,   0.17534075, 0.00296513, 0.17534075, 0.00296513, -6.3802411,  -1 ],
        [ 8.26437678,   0.17692329, 0.00226267, 0.17692329, 0.00226267, -6.9891369,  -1 ],
        [ 8.92981497,   0.17822456, 0.00167624, 0.17822456, 0.00167624, -7.6545751,  -1 ],
        [ 9.65599747,   0.17926081, 0.00120286, 0.17926081, 0.00120286, -8.3807576,  -1 ],
        [ 10.45663841,  0.18006594, 0.00083056, 0.18006594, 0.00083056, -9.1813986,  -1 ],
        [ 11.34781806,  0.18067147, 0.00054749, 0.18067147, 0.00054749, -10.072578,  -1 ],
        [ 12.35501198,  0.18111044, 0.00034023, 0.18111044, 0.00034023, -11.079772,  -1 ],
        [ 13.51466689,  0.18141365, 0.00019575, 0.18141365, 0.00019575, -12.239427,  -1 ],
        [ 14.88580733,  0.18161029, 0.00010125, 0.18161029, 0.00010125, -13.610567,  -1 ],
        [ 15.83854747,  0.18168759, 6.385e-05,  0.18168759, 6.385e-05,  -14.563308,  -1 ]
    ],

    impulse_table => [
        [
            -11.884907,   13.814887,  -11.885774, 3.9367976, 0, 0, 0, -7.3014346,
            0.0027163799, -83.614233, 0.33313834, 1.0000291, -9.6007894
        ],
        [
            -2.302585,    4.2536961, -2.4113709, 3.7875796, 0, 0, 0, -2.5864922,
            0.0027163831, -83.51424, 0.32755761, 0.8972727, -2.4769059
        ],
        [
            -1.6094379,   3.5804891, -1.7648022, 3.7386611,  0, 0, 0, -2.3051612,
            0.0027163947, -83.41424, 0.32230877, 0.84487975, -2.0651753
        ],
        [
            -1.3862943,   3.3669132, -1.5605003, 3.721709,  0, 0, 0, -2.2236016,
            0.0027164041, -83.36424, 0.31979765, 0.8242567, -1.9404753
        ],
        [
            -1.2039728,   3.1939115, -1.3952109, 3.7077387,  0, 0, 0, -2.1610545,
            0.0027164161, -83.31424, 0.31735678, 0.80600546, -1.8417743
        ],
        [
            -0.91629068,  2.9241799, -1.1377026, 3.685997,   0, 0, 0, -2.0708872,
            0.0027164479, -83.21424, 0.31267176, 0.77468082, -1.692263
        ],
        [
            -0.69314713,  2.7181453, -0.94100511, 3.6699011,  0, 0, 0, -2.0090124,
            0.0027164906, -83.11424, 0.30822819,  0.74834257, -1.581864
        ],
        [
            -0.51082557,  2.5521734, -0.78243014, 3.6576141, 0, 0, 0, -1.9642315,
            0.0027165445, -83.01424, 0.3040041,   0.7255951, -1.4954651
        ],
        [
            -0.35667489,  2.4136829, -0.64994759, 3.6480462,  0, 0, 0, -1.9307016,
            0.0027166099, -82.91424, 0.29998043,  0.70557678, -1.4251789
        ],
        [
            -0.2231435,   2.2951823, -0.53641619, 3.640496,   0, 0, 0, -1.905017,
            0.0027166869, -82.81424, 0.2961405,   0.68771162, -1.3663905
        ],
        [
            5.3377405e-08, 2.1004169,  -0.34935085,  3.6296775, 0,          0,
            0,             -1.8693363, 0.0027168763, -82.61424, 0.28895487, 0.65692036,
            -1.2725813
        ],
        [
            0.26236432,   1.8770226, -0.13383232, 3.6202733,  0, 0, 0, -1.8393455,
            0.0027172501, -82.31424, 0.27923761,  0.61954435, -1.1694899
        ],
        [
            0.5306283,    1.6552937, 0.081406627, 3.6147738,  0, 0, 0, -1.822265,
            0.0027179178, -81.91424, 0.2678877,   0.58050035, -1.0721907
        ],
        [
            0.69314723,   1.5243929, 0.20923095, 3.6135507,  0, 0, 0, -1.8185586,
            0.0027185467, -81.61424, 0.26034789, 0.55666661, -1.0172136
        ],
        [
            0.99325183,   1.2896111,  0.44014074, 3.6154081,  0, 0, 0, -1.8245095,
            0.0027204438, -80.914239, 0.245247,   0.51276567, -0.923421
        ],
        [
            1.0986123,    1.2093232,  0.51963851, 3.6172656,  0, 0, 0, -1.8303959,
            0.0027214418, -80.614239, 0.23962983, 0.49749294, -0.89281812
        ],
        [
            1.2089604,    1.1264176,  0.60203458, 3.6198357,  0, 0, 0, -1.8385896,
            0.0027227467, -80.264239, 0.23359959, 0.48162534, -0.86202998
        ],
        [
            1.3862944,    0.99568771, 0.73261997, 3.6252063,  0, 0, 0, -1.855942,
            0.0027255732, -79.614239, 0.22364655, 0.45647442, -0.81519851
        ],
        [
            1.609438,     0.8354655,  0.89380892, 3.6338784,  0, 0, 0, -1.8847094,
            0.0027309499, -78.614239, 0.21079238, 0.42558986, -0.76075295
        ],
        [
            1.7917595,    0.70796541, 1.0230049,  3.6422782,  0, 0, 0, -1.9135568,
            0.0027375826, -77.614239, 0.20014237, 0.40111112, -0.7198248
        ],
        [
            1.9459102,   0.60245716, 1.1305456, 3.6501073,  0, 0, 0, -1.941428,
            0.002745484, -76.614239, 0.191118,  0.38101443, -0.68761142
        ],
        [
            2.0794416,    0.5126877,  1.2224923,  3.6573012, 0, 0, 0, -1.9679846,
            0.0027546699, -75.614239, 0.18333524, 0.3640852, -0.66141307
        ],
        [
            2.3025851,    0.36587391, 1.373744,   3.6698787,  0, 0, 0, -2.0170053,
            0.0027769719, -73.614239, 0.17050192, 0.33684946, -0.62102893
        ],
        [
            2.7080503,    0.10850272, 1.6414191, 3.6929602,  0, 0, 0, -2.1193387,
            0.0028568669, -68.614239, 0.1481363, 0.29093138, -0.55794122
        ],
        [
            2.9957323, -0.06755279, 1.8262389,    3.7077391,  0,          0,
            0,         -2.2009499,  0.0029749671, -63.614239, 0.13329553, 0.26123406,
            -0.52084055
        ],
        [
            3.4011974, -0.30775505, 2.0803683,    3.7215162,  0,          0,
            0,         -2.3266706,  0.0033568628, -53.614238, 0.11411211, 0.22342807,
            -0.47934543
        ],
        [
            3.9120231, -0.59931918, 2.391343,     3.7022658,  0,           0,
            0,         -2.4999963,  0.0052597973, -33.614237, 0.093010515, 0.18224841,
            -0.44992843
        ],
        [
            4.2484953,   -0.78584204, 2.5913807,   3.5795598,  0, 0, 0, -2.6229865,
            0.014342098, -13.614236,  0.080952086, 0.15857034, -0.46941211
        ]
    ],

    blast_info => [ 0, 0, 0, 0, 0, 0, 3.9338399, 0, 0, 0, 0, 0, 0.2402124, 0.33313836, 0, 0, 0, 0, 0, 0 ],

};
1;