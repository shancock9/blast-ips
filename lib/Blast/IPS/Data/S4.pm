package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=4
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S4'} = {

    table_name => 'S4',
    symmetry   => 2,
    gamma      => 4,

    shock_table_info => [ 2, 4, 9.1e-07, 4000, 9.36e-07, 10, 1.8e-07, 2.3e-07, 5e-07, 83.486, -0.80308 ],

    shock_table => [
        [ -4.06012486,  12.0003368,    -2.99997052, -4.0613796,  0.99811673 ],
        [ -3.77304581,  11.13911052,   -2.99993025, -3.7749765,  0.99710123 ],
        [ -3.50813588,  10.34440347,   -2.99987489, -3.51100987, 0.995683 ],
        [ -3.08513094,  9.07549441,    -2.99955969, -3.09055818, 0.99183824 ],
        [ -2.79827581,  8.21512963,    -2.99895525, -2.80663259, 0.98741677 ],
        [ -2.52888727,  7.40739589,    -2.99766587, -2.54142854, 0.98108469 ],
        [ -2.32152752,  6.78598772,    -2.99565644, -2.3386779,  0.97409093 ],
        [ -2.13822007,  6.2371249,     -2.99249414, -2.16084765, 0.96575947 ],
        [ -1.96556957,  5.72086622,    -2.98747014, -1.99495931, 0.95545244 ],
        [ -1.81889501,  5.28314372,    -2.98067178, -1.85560762, 0.94427819 ],
        [ -1.68050588,  4.87127605,    -2.97100652, -1.72580441, 0.93118 ],
        [ -1.55771,     4.50716656,    -2.95860831, -1.61230133, 0.91703323 ],
        [ -1.43768318,  4.15301525,    -2.94167706, -1.50319698, 0.90048164 ],
        [ -1.32292205,  3.81663105,    -2.91955875, -1.40090298, 0.88173974 ],
        [ -1.21073281,  3.49063081,    -2.89072625, -1.30315349, 0.86030134 ],
        [ -1.09738189,  3.16502998,    -2.85263427, -1.20702615, 0.83520318 ],
        [ -0.97556363,  2.82061231,    -2.79977954, -1.10713114, 0.8041482 ],
        [ -0.85470987,  2.48610007,    -2.73369977, -1.01202092, 0.76914089 ],
        [ -0.744536,    2.18880319,    -2.66126054, -0.92919903, 0.73383032 ],
        [ -0.64164319,  1.918888,      -2.58380067, -0.85549919, 0.69837278 ],
        [ -0.54423337,  1.67107746,    -2.50316548, -0.78917581, 0.66314924 ],
        [ -0.45058834,  1.44049688,    -2.42077141, -0.72869935, 0.62836142 ],
        [ -0.35765714,  1.21944389,    -2.33629758, -0.67192274, 0.59354592 ],
        [ -0.26441559,  1.00558891,    -2.25087209, -0.61819887, 0.55890895 ],
        [ -0.17051991,  0.79823954,    -2.16606501, -0.56732649, 0.52486509 ],
        [ -0.07490154,  0.5951437,     -2.08258549, -0.51874612, 0.4915208 ],
        [ 0.02318007,   0.39490837,    -2.00125956, -0.47214397, 0.4590725 ],
        [ 0.12462065,   0.19593731,    -1.92262967, -0.42718861, 0.42764194 ],
        [ 0.2303159,    -0.00322678,   -1.84713152, -0.38361253, 0.3973408 ],
        [ 0.34118779,   -0.20396148,   -1.77510253, -0.34119576, 0.36826724 ],
        [ 0.45807957,   -0.40739317,   -1.70686255, -0.29979777, 0.34053215 ],
        [ 0.58213224,   -0.61506035,   -1.64251067, -0.25922035, 0.31417067 ],
        [ 0.71442053,   -0.8282653,    -1.58217242, -0.21934269, 0.2892375 ],
        [ 0.85639255,   -1.04879255,   -1.52580717, -0.17998542, 0.26572373 ],
        [ 1.00623969,   -1.27348596,   -1.47444548, -0.14182708, 0.24407622 ],
        [ 1.16176067,   -1.49913996,   -1.42863115, -0.10542187, 0.22455513 ],
        [ 1.32406063,   -1.72759452,   -1.38764221, -0.07044451, 0.20688769 ],
        [ 1.49467672,   -1.9611226,    -1.35079491, -0.03655092, 0.19080932 ],
        [ 1.67460167,   -2.20110296,   -1.31765083, -0.00357006, 0.17615556 ],
        [ 1.86520866,   -2.44933201,   -1.28777574, 0.02869822,  0.16275954 ],
        [ 2.06772835,   -2.70732876,   -1.26083395, 0.06038734,  0.15049448 ],
        [ 2.28344228,   -2.97661414,   -1.23653369, 0.09160795,  0.1392508 ],
        [ 2.51360846,   -3.25862987,   -1.21462501, 0.12244177,  0.12893617 ],
        [ 2.76012872,   -3.55555099,   -1.19483919, 0.15302778,  0.11944698 ],
        [ 3.02453568,   -3.86904572,   -1.17698179, 0.1834262,   0.11071281 ],
        [ 3.30870812,   -4.20115414,   -1.16086142, 0.21371463,  0.1026631 ],
        [ 3.61471813,   -4.55409652,   -1.14630509, 0.24396468,  0.09523458 ],
        [ 3.94483085,   -4.93027181,   -1.13315716, 0.27424065,  0.08837095 ],
        [ 4.30170472,   -5.33248173,   -1.12127153, 0.30461511,  0.08201896 ],
        [ 4.68819061,   -5.76369724,   -1.11052009, 0.33514725,  0.07613288 ],
        [ 5.10953361,   -6.22948481,   -1.10074304, 0.36603923,  0.07064667 ],
        [ 5.56096555,   -6.72436338,   -1.09200319, 0.39676733,  0.06561973 ],
        [ 6.05632219,   -7.2632526,    -1.08400211, 0.42807349,  0.06090237 ],
        [ 6.59622469,   -7.8464864,    -1.0767339,  0.45973799,  0.05651022 ],
        [ 7.18529096,   -8.47874838,   -1.07012494, 0.49178934,  0.05241854 ],
        [ 7.82893411,   -9.16553124,   -1.06410616, 0.52426783,  0.04860325 ],
        [ 8.53296325,   -9.91270256,   -1.05861812, 0.55719977,  0.04504411 ],
        [ 9.30358436,   -10.72650771,  -1.05360937, 0.59059823,  0.04172384 ],
        [ 10.14780117,  -11.61399239,  -1.04903312, 0.62447955,  0.03862615 ],
        [ 11.07281584,  -12.58237039,  -1.04485019, 0.65883744,  0.03573787 ],
        [ 12.08642955,  -13.63944808,  -1.04102553, 0.69366103,  0.03304679 ],
        [ 13.19744286,  -14.79403916,  -1.03752668, 0.72894642,  0.03054074 ],
        [ 14.41525604,  -16.05554479,  -1.03432521, 0.76468057,  0.02820883 ],
        [ 15.74946919,  -17.4335428,   -1.03139677, 0.80083172,  0.0260417 ],
        [ 17.21088235,  -18.93882359,  -1.02871817, 0.83737854,  0.02402956 ],
        [ 18.72257513,  -20.49212194,  -1.02639287, 0.87233021,  0.02225858 ],
        [ 20.35280199,  -22.16360133,  -1.02427615, 0.90725135,  0.02062603 ],
        [ 22.19989372,  -24.05361364,  -1.02225644, 0.94385379,  0.01904931 ],
        [ 24.08012356,  -25.97401434,  -1.02052059, 0.97834959,  0.01767876 ],
        [ 26.20299469,  -28.13863721,  -1.01886146, 1.01443682,  0.01635486 ],
        [ 28.61136895,  -30.59047432,  -1.01727817, 1.05224566,  0.01507825 ],
        [ 31.03416766,  -33.05346181,  -1.0159338,  1.08741734,  0.01398371 ],
        [ 33.76931364,  -35.8303813,   -1.01464814, 1.12418183,  0.01292748 ],
        [ 36.87207648,  -38.97663369,  -1.0134206,  1.1626703,   0.01190998 ],
        [ 39.9416289,   -42.08576341,  -1.01239376, 1.19787809,  0.01105178 ],
        [ 43.40012307,  -45.58537428,  -1.01141063, 1.23463114,  0.01022384 ],
        [ 47.31530294,  -49.54333624,  -1.01047085, 1.27305539,  0.00942641 ],
        [ 51.77026365,  -54.04288681,  -1.00957405, 1.31329304,  0.00865977 ],
        [ 56.09570516,  -58.40810889,  -1.00883931, 1.34934755,  0.00802734 ],
        [ 60.97223854,  -63.32598428,  -1.00813567, 1.38696746,  0.00741786 ],
        [ 66.49663213,  -68.89341063,  -1.0074629,  1.426282,    0.0068315 ],
        [ 71.68016321,  -74.11419882,  -1.00692568, 1.46044347,  0.00636062 ],
        [ 77.47923077,  -79.95189387,  -1.00640963, 1.49597789,  0.00590597 ],
        [ 83.9944203,   -86.50718747,  -1.0059146,  1.53299138,  0.00546763 ],
        [ 91.34766729,  -93.90213419,  -1.00544048, 1.57160297,  0.00504569 ],
        [ 99.6879264,   -102.28582292, -1.00498713, 1.61194678,  0.00464022 ],
        [ 109.19867586, -111.84188386, -1.00455444, 1.65417475,  0.00425129 ],
        [ 117.80145766, -120.48238529, -1.00422308, 1.68942915,  0.0039521 ],
        [ 127.4409947,  -130.16105662, -1.00390478, 1.72609958,  0.00366359 ],
        [ 138.2896973,  -141.05042006, -1.00359949, 1.76429789,  0.00338578 ],
        [ 150.55727236, -153.3603083,  -1.00330714, 1.80414981,  0.00311871 ],
        [ 164.50082685, -167.34797036, -1.00302767, 1.84579738,  0.00286241 ],
        [ 176.24686013, -179.12835812, -1.00282648, 1.87830762,  0.00267727 ],
        [ 189.27511817, -192.19214653, -1.00263247, 1.91199441,  0.00249823 ],
        [ 203.77825776, -206.73207693, -1.00244561, 1.94694228,  0.00232529 ]
    ],

    energy_table => [
        [ -4.06012486,  0.00010795, 0.00023751, 0.00010795, 0.00023751,    5.1472686,    -1.5613105 ],
        [ -3.77304581,  0.00020265, 0.00044326, 0.00020265, 0.00044326,    4.6974973,    -1.5727127 ],
        [ -3.50813588,  0.00036104, 0.00078424, 0.00036104, 0.00078424,    4.2794042,    -1.5847448 ],
        [ -3.08513094,  0.00089898, 0.00192284, 0.00089898, 0.00192284,    3.60466,      -1.6078351 ],
        [ -2.79827581,  0.00165329, 0.00348431, 0.00165329, 0.00348431,    3.1409711,    -1.627457 ],
        [ -2.52888727,  0.00290229, 0.00600366, 0.00290229, 0.00600366,    2.69977,      -1.6499527 ],
        [ -2.32152752,  0.00443982, 0.00901392, 0.00443982, 0.00901392,    2.3556944,    -1.6685357 ],
        [ -2.13822007,  0.00641774, 0.01276096, 0.00641774, 0.01276096,    2.0483455,    -1.6865449 ],
        [ -1.96556957,  0.00901012, 0.01748265, 0.00901012, 0.01748265,    1.7555594,    -1.705641 ],
        [ -1.81889501,  0.0119348,  0.02257353, 0.0119348,  0.02257353,    1.5041629,    -1.7227693 ],
        [ -1.68050588,  0.01544901, 0.02838021, 0.01544901, 0.02838021,    1.264602,     -1.7399087 ],
        [ -1.55771,     0.0192934,  0.0343609,  0.0192934,  0.0343609,     1.049985,     -1.7557282 ],
        [ -1.43768318,  0.02380362, 0.04089663, 0.02380362, 0.04089663,    0.83831502,   -1.7713127 ],
        [ -1.32292205,  0.02888037, 0.04764305, 0.02888037, 0.04764305,    0.6341824,    -1.7857459 ],
        [ -1.21073281,  0.03460827, 0.05448581, 0.03460827, 0.05448581,    0.43307484,   -1.7987768 ],
        [ -1.09738189,  0.04117464, 0.06133055, 0.04117464, 0.06133055,    0.22847151,   -1.8099816 ],
        [ -0.97556363,  0.04906936, 0.068147,   0.04906936, 0.068147,      0.0073362096, -1.8190678 ],
        [ -0.85470987,  0.05766145, 0.07381489, 0.05766145, 0.07381489,    -0.21295808,  -1.8251797 ],
        [ -0.744536,    0.06601567, 0.07758797, 0.06601567, 0.07758797,    -0.41428249,  -1.8288677 ],
        [ -0.64164319,  0.07411726, 0.07963597, 0.07411726, 0.07963597,    -0.60260725,  -1.8319488 ],
        [ -0.54423337,  0.08191093, 0.0801486,  0.08191093, 0.0801486,     -0.78120906,  -1.8359243 ],
        [ -0.45058834,  0.08938854, 0.07934821, 0.08938854, 0.07934821,    -0.95335166,  -1.8417758 ],
        [ -0.35765714,  0.09668066, 0.07741122, 0.09668066, 0.07741122,    -1.1248356,   -1.8500818 ],
        [ -0.26441559,  0.10376979, 0.07450523, 0.10376979, 0.07450523,    -1.2977895,   -1.860861 ],
        [ -0.17051991,  0.11059794, 0.07083095, 0.11059794, 0.07083095,    -1.4730798,   -1.8736712 ],
        [ -0.07490154,  0.11716998, 0.06656519, 0.11716998, 0.06656519,    -1.6529002,   -1.8879243 ],
        [ 0.02318007,   0.12347092, 0.06188657, 0.12347092, 0.06188657,    -1.8388072,   -1.9028657 ],
        [ 0.12462065,   0.12949848, 0.05695338, 0.12949848, 0.05695338,    -2.0326149,   -1.9178174 ],
        [ 0.2303159,    0.13524998, 0.0519082,  0.13524998, 0.0519082,     -2.236119,    -1.9322085 ],
        [ 0.34118779,   0.14072316, 0.046876,   0.14072316, 0.046876,      -2.4511401,   -1.9455504 ],
        [ 0.45807957,   0.1459114,  0.04196788, 0.1459114,  0.04196788,    -2.6793212,   -1.9575034 ],
        [ 0.58213224,   0.15082039, 0.0372648,  0.15082039, 0.0372648,     -2.9228699,   -1.9679055 ],
        [ 0.71442053,   0.15545049, 0.03283453, 0.15545049, 0.03283453,    -3.1838549,   -1.976642 ],
        [ 0.85639255,   0.15981237, 0.02871784, 0.15981237, 0.02871784,    -3.4650604,   -1.9838277 ],
        [ 1.00623969,   0.16383071, 0.02501825, 0.16383071, 0.02501825,    -3.7628243,   -1.9894021 ],
        [ 1.16176067,   0.16746284, 0.02178586, 0.16746284, 0.02178586,    -4.0725868,   -1.9934468 ],
        [ 1.32406063,   0.17076302, 0.01896781, 0.17076302, 0.01896781,    -4.3964066,   -1.9963647 ],
        [ 1.49467672,   0.17378287, 0.01650913, 0.17378287, 0.01650913,    -4.7372287,   -1.9983895 ],
        [ 1.67460167,   0.17655454, 0.01436938, 0.17655454, 0.01436938,    -5.0969393,   -1.9997225 ],
        [ 1.86520866,   0.17911014, 0.01250774, 0.17911014, 0.01250774,    -5.4782006,   -2.0005343 ],
        [ 2.06772835,   0.18147386, 0.0108901,  0.18147386, 0.0108901,     -5.8834098,   -2.0009737 ],
        [ 2.28344228,   0.18366636, 0.009486,   0.18366636, 0.009486,      -6.3150788,   -2.0011331 ],
        [ 2.51360846,   0.18570475, 0.00826869, 0.18570475, 0.00826869,    -6.7756757,   -2.0011345 ],
        [ 2.76012872,   0.18760831, 0.00721197, 0.18760831, 0.00721197,    -7.268992,    -2.0010767 ],
        [ 3.02453568,   0.18938972, 0.00629535, 0.18938972, 0.00629535,    -7.7980787,   -2.0009261 ],
        [ 3.30870812,   0.19106163, 0.00549996, 0.19106163, 0.00549996,    -8.3666544,   -2.0007413 ],
        [ 3.61471813,   0.19263521, 0.00480937, 0.19263521, 0.00480937,    -8.9788779,   -2.0005943 ],
        [ 3.94483085,   0.19412022, 0.0042093,  0.19412022, 0.0042093,     -9.6392744,   -2.0004472 ],
        [ 4.30170472,   0.19552587, 0.00368718, 0.19552587, 0.00368718,    -10.353154,   -2.0003167 ],
        [ 4.68819061,   0.19685983, 0.00323233, 0.19685983, 0.00323233,    -11.126226,   -2.0002189 ],
        [ 5.10953361,   0.19813472, 0.00283381, 0.19813472, 0.00283381,    -11.968986,   -2.0001458 ],
        [ 5.56096555,   0.19933368, 0.00249038, 0.19934935, 0.0025547706,  -12.871901,   -2.0000928 ],
        [ 6.05632219,   0.20048947, 0.00218715, 0.20058027, 0.0024771825,  -13.862649,   -2.0000563 ],
        [ 6.59622469,   0.20159596, 0.00192144, 0.20196685, 0.0026071894,  -14.942477,   -2.0000327 ],
        [ 7.18529096,   0.20265672, 0.00168787, 0.20341327, 0.0023177068,  -16.120623,   -2.0000181 ],
        [ 7.82893411,   0.20367517, 0.0014858,  0.20483853, 0.0021207868,  -17.407917,   -2.0000095 ],
        [ 8.53296325,   0.20465412, 0.0013037,  0.20626519, 0.0019389029,  -18.81598,    -2.0000047 ],
        [ 9.30358436,   0.20559552, 0.00114529, 0.20769522, 0.0017743868,  -20.357224,   -2.0000023 ],
        [ 10.14780117,  0.20650139, 0.00100587, 0.209128,   0.0016229758,  -22.045659,   -2.0000011 ],
        [ 11.07281584,  0.20737302, 0.00088317, 0.2105594,  0.0014764682,  -23.895689,   -2.0000005 ],
        [ 12.08642955,  0.20821149, 0.0007752,  0.21198412, 0.0013384668,  -25.922917,   -2.0000003 ],
        [ 13.19744286,  0.20901805, 0.0006802,  0.21339719, 0.0012081941,  -28.144944,   -2.0000002 ],
        [ 14.41525604,  0.20979366, 0.00059664, 0.21479238, 0.0010862105,  -30.58057,    -2.0000002 ],
        [ 15.74946919,  0.21053891, 0.0005232,  0.21616285, 0.00097088682, -33.248997,   -2.0000002 ],
        [ 17.21088235,  0.21125465, 0.00045868, 0.21758375, 0.0009309355,  -36.171823,   -2.0000002 ],
        [ 18.72257513,  0.21190594, 0.00040488, 0.21890477, 0.00082069234, -39.195209,   -2 ],
        [ 20.35280199,  0.21252622, 0.00035775, 0.22016143, 0.00072436562, -42.455663,   -2 ],
        [ 22.19989372,  0.21314564, 0.00031453, 0.22141502, 0.00063619982, -46.149846,   -2 ],
        [ 24.08012356,  0.2137023,  0.00027882, 0.22254055, 0.00056348415, -49.910306,   -2 ],
        [ 26.20299469,  0.21425812, 0.000246,   0.22366342, 0.00049674627, -54.156048,   -2 ],
        [ 28.61136895,  0.21481301, 0.00021593, 0.22478357, 0.0004356938,  -58.972797,   -2 ],
        [ 31.03416766,  0.21530542, 0.0001914,  0.22577689, 0.00038596616, -63.818394,   -2 ],
        [ 33.76931364,  0.215797,   0.00016886, 0.22676793, 0.00034030658, -69.288686,   -2 ],
        [ 36.87207648,  0.21628769, 0.00014821, 0.22775663, 0.00029851938, -75.494212,   -2 ],
        [ 39.9416289,   0.21671627, 0.00013161, 0.22861977, 0.00026498145, -81.633317,   -2 ],
        [ 43.40012307,  0.21714411, 0.00011634, 0.22948103, 0.00023413442, -88.550305,   -2 ],
        [ 47.31530294,  0.21757116, 0.00010233, 0.23034036, 0.00020585316, -96.380665,   -2 ],
        [ 51.77026365,  0.21799739, 8.952e-05,  0.23119772, 0.00018001395, -105.29059,   -2 ],
        [ 56.09570516,  0.21836206, 7.945e-05,  0.23193101, 0.00015971745, -113.94147,   -2 ],
        [ 60.97223854,  0.21872607, 7.018e-05,  0.23266278, 0.00014104988, -123.69454,   -2 ],
        [ 66.49663213,  0.21908942, 6.168e-05,  0.23339303, 0.00012393618, -134.74332,   -2 ],
        [ 71.68016321,  0.21939169, 5.516e-05,  0.23400037, 0.00011080832, -145.11039,   -2 ],
        [ 77.47923077,  0.21969346, 4.913e-05,  0.23460661, 9.8665926e-05, -156.70852,   -2 ],
        [ 83.9944203,   0.21999473, 4.356e-05,  0.23521173, 8.746668e-05,  -169.7389,    -2 ],
        [ 91.34766729,  0.22029548, 3.844e-05,  0.23581572, 7.7169292e-05, -184.44539,   -2 ],
        [ 99.6879264,   0.22059571, 3.374e-05,  0.23641855, 6.7732444e-05, -201.12591,   -2 ],
        [ 109.19867586, 0.22089539, 2.945e-05,  0.23702021, 5.9115656e-05, -220.14741,   -2 ],
        [ 117.80145766, 0.22113474, 2.63e-05,   0.23750069, 5.2785594e-05, -237.35297,   -2 ],
        [ 127.4409947,  0.22137373, 2.339e-05,  0.23798039, 4.6933924e-05, -256.63205,   -2 ],
        [ 138.2896973,  0.22161236, 2.07e-05,   0.23845932, 4.1540574e-05, -278.32945,   -2 ],
        [ 150.55727236, 0.22185061, 1.824e-05,  0.23893745, 3.6585392e-05, -302.8646,    -2 ],
        [ 164.50082685, 0.22208849, 1.598e-05,  0.23941479, 3.2048522e-05, -330.75171,   -2 ],
        [ 176.24686013, 0.22226664, 1.441e-05,  0.23977226, 2.8908491e-05, -354.24378,   -2 ],
        [ 189.27511817, 0.22244458, 1.296e-05,  0.24012928, 2.5984376e-05, -380.3003,    -2 ],
        [ 203.77825776, 0.22262229, 1.16e-05,   0.24048583, 2.3267944e-05, -409.30657,   -2 ]
    ],

    impulse_table => [
        [
            -4.6649614, 13.814838,   -4.6654677,    0.039016636, -0.0038296439, -0.067115865,
            -2.6460003, -0.28598157, -0.0049773817, -0.28711853, 0.44710028,    0.99988522,
            -3.9111115
        ],
        [
            -4.6051695, 13.635463,   -4.6057233,   0.040093032, -0.0040656064, -0.066536047,
            -2.64542,   -0.31587925, -0.005284063, -0.28653814, 0.44710019,    0.99986914,
            -3.8733803
        ],
        [
            -4.1997044, 12.419072,   -4.200722,     0.048097909, -0.0060983333, -0.061544387,
            -2.640422,  -0.51863325, -0.0079260226, -0.28153813, 0.44709917,    0.99968119,
            -3.621975
        ],
        [
            -3.9120223, 11.556033,   -3.9135893,   0.054560233, -0.0081309129, -0.056560634,
            -2.6354258, -0.66251398, -0.010567844, -0.27653813, 0.44709728,    0.99940077,
            -3.4489777
        ],
        [
            -3.5065572, 10.339668,   -3.509438,    0.064815411, -0.012195149, -0.04662733,
            -2.6254415, -0.86540354, -0.015850616, -0.26666768, 0.4470896,    0.99854892,
            -3.2144358
        ],
        [
            -2.9957316, 8.8073412,  -3.0019395,   0.079589914, -0.020314782, -0.026970272,
            -2.6055224, -1.1215787, -0.026407824, -0.24692692, 0.44705025,   0.9956404,
            -2.9387841
        ],
        [
            -2.3025844, 6.7292429,  -2.3202327,   0.10217077,  -0.040446184, 0.020045688,
            -2.5562455, -1.4742092, -0.052640831, -0.19979101, 0.44669985,   0.98141878,
            -2.6169309
        ],
        [
            -1.8971193, 5.5164657,  -1.929723,  0.11582212,  -0.059965219, 0.062542295,
            -2.5082027, -1.6914862, -0.0782682, -0.15718372, 0.44575516,   0.95824898,
            -2.4724464
        ],
        [
            -1.6094372, 4.6603588,  -1.6599012,  0.12515101,  -0.078313602, 0.099718304,
            -2.4619888, -1.8606417, -0.10265411, -0.12092833, 0.44394148,   0.92793424,
            -2.4005353
        ],
        [
            -1.3862937, 4.0020717,  -1.457125,   0.13193324,   -0.094991609, 0.1317752,
            -2.4181533, -2.0097021, -0.12510458, -0.091771712, 0.44102337,   0.89255255,
            -2.3690293
        ],
        [
            -1.2039721, 3.4710942,  -1.297342,   0.13710156,   -0.10967401, 0.1592749,
            -2.3771605, -2.1511366, -0.14503867, -0.069619611, 0.4368266,   0.85415746,
            -2.363456
        ],
        [
            -1.0498214, 3.0298027, -1.1675771,  0.14122958,   -0.12223232, 0.18276143,
            -2.3393462, -2.291109, -0.16211164, -0.053835804, 0.43125621,  0.81456152,
            -2.3757521
        ],
        [
            -0.91629002, 2.6555511,  -1.0599552, 0.14468427,   -0.13270756, 0.20269671,
            -2.3048833,  -2.4324948, -0.1762533, -0.043771708, 0.42430517,  0.77521536,
            -2.4007521
        ],
        [
            -0.79850699, 2.3334426, -0.96928484, 0.14769127,   -0.14126534, 0.2194985,
            -2.2737706,  -2.576196, -0.1876323,  -0.037901833, 0.41605082,  0.73717429,
            -2.4348385
        ],
        [
            -0.69314647, 2.0529986,  -0.89193381, 0.15037828,   -0.14814682, 0.23357217,
            -2.2458545,  -2.7219001, -0.19657842, -0.035847967, 0.40664056,  0.70112774,
            -2.4753446
        ],
        [
            -0.59783629, 1.8064726,  -0.82524812, 0.15281042,   -0.15362294, 0.24531541,
            -2.2208725,  -2.8686243, -0.20349584, -0.036435196, 0.39627109,  0.66746314,
            -2.5202565
        ],
        [
            -0.51082491, 1.587933,   -0.76722697, 0.15501829,   -0.15795866, 0.25510629,
            -2.1985104,  -3.0151248, -0.20879105, -0.038916971, 0.38516613,  0.63634058,
            -2.5680426
        ],
        [
            -0.43078221, 1.3927275,  -0.71632743, 0.15701698,   -0.16139036, 0.26328773,
            -2.1784435,  -3.1601794, -0.21282752, -0.042627097, 0.37355669,  0.60776249,
            -2.6175425
        ],
        [
            -0.35667423, 1.217148,   -0.67133952, 0.15881708,   -0.16411597, 0.27015704,
            -2.1603685,  -3.3027469, -0.21590575, -0.047397311, 0.36166581,  0.581631,
            -2.6678843
        ],
        [
            -0.28768136, 1.0582048,  -0.63130216, 0.16042968,   -0.16629422, 0.27596292,
            -2.1440155,  -3.4420301, -0.21826196, -0.052670439, 0.34969808,  0.55779077,
            -2.7184212
        ],
        [
            -0.22314284, 0.91346518, -0.59544344, 0.1618677,    -0.16804888, 0.28090806,
            -2.1291526,  -3.5774746, -0.220076,   -0.058687884, 0.33783342,  0.53605891,
            -2.7686807
        ],
        [
            -0.16251822, 0.78093588, -0.56313809, 0.16314568,   -0.16947493, 0.28515484,
            -2.1155842,  -3.7087368, -0.2214821,  -0.064566938, 0.32622413,  0.5162446,
            -2.8183239
        ],
        [
            -0.10535981, 0.65897373,  -0.53387582, 0.1642789,  -0.17064469, 0.28883229, -2.1031471, -3.8356404,
            -0.22257992, -0.07059534, 0.31499414,  0.49816107, -2.8671136
        ],
        [
            -0.051292584, 0.54621291, -0.50723632, 0.16528259,   -0.17161313, 0.29204257,
            -2.0917039,   -3.9581347, -0.22344315, -0.076541649, 0.30423776,  0.48163236,
            -2.9148889
        ],
        [
            -0.040991776, 0.52501139, -0.50227331, 0.16546607,   -0.17178376, 0.29262645,
            -2.0895588,   -3.9817157, -0.22359034, -0.077562495, 0.30218327,  0.47854954,
            -2.9241585
        ],
        [
            0.09531089, 0.25261201, -0.43985193, 0.16765659,   -0.17368302, 0.29960086,
            -2.0622748, -4.2998323, -0.22510582, -0.093725848, 0.27528348,  0.43983407,
            -3.0512951
        ],
        [
            0.18232227,  0.086221052, -0.40300291, 0.16882744, -0.17460629, 0.3034117, -2.045934, -4.5075469,
            -0.22573833, -0.10494026, 0.25881298,  0.41717992, -3.1361962
        ],
        [
            0.26236497, -0.062078377, -0.37101794, 0.16975667,  -0.17530157, 0.3065465,
            -2.0316574, -4.7008179,   -0.22615182, -0.11494379, 0.24442798,  0.39775246,
            -3.2163599
        ],
        [
            0.40546582,  -0.31681255, -0.31803076, 0.17109936, -0.1762612, 0.31140641, -2.0079745, -5.0496353,
            -0.22660002, -0.13400098, 0.22090348,  0.3662697,  -3.3635088
        ],
        [
            0.53062896,  -0.52980838, -0.27567127, 0.17198584, -0.1768788, 0.31500871, -1.9892217, -5.3563911,
            -0.22678193, -0.15119686, 0.20276686,  0.34192726, -3.495147
        ],
        [
            0.69314789,  -0.79451128, -0.22553581, 0.17281765, -0.17746688, 0.31896239, -1.9676094, -5.7547796,
            -0.22683884, -0.17371543, 0.18247222,  0.31430269, -3.6686554
        ],
        [
            0.99325248,  -1.2543102,  -0.14500833, 0.17369011, -0.17815963, 0.32462516, -1.9355466, -6.485224,
            -0.22667902, -0.21567552, 0.15321386,  0.27316565, -3.9925041
        ],
        [
            1.098613,    -1.4083706,  -0.11983936, 0.17385896, -0.17832748, 0.32621986, -1.9265349, -6.7391438,
            -0.22658784, -0.23090861, 0.14499651,  0.26122316, -4.1064233
        ],
        [
            1.2089611,  -1.566273,  -0.094951039, 0.17398432,  -0.17847574, 0.32771265,
            -1.9182274, -7.0034448, -0.22648608,  -0.24651568, 0.13731615,  0.2498698,
            -4.2255883
        ],
        [
            1.3862951,  -1.8135105, -0.057761485, 0.17410479,  -0.17866983, 0.32978377,
            -1.9070728, -7.4245432, -0.22632131,  -0.27137298, 0.12665931,  0.23377133,
            -4.416513
        ],
        [
            1.6094386,  -2.1148735, -0.015212193, 0.17416331,  -0.17885941, 0.33191283,
            -1.8963375, -7.9480343, -0.22612678,  -0.30409781, 0.1156401,   0.21664518,
            -4.6553869
        ],
        [
            1.7917602,   -2.3543496,  0.01656464, 0.17416569, -0.17898271, 0.33333123, -1.8897931, -8.3706697,
            -0.22598539, -0.33036386, 0.10819193, 0.20475518, -4.8492672
        ],
        [
            1.9459109,   -2.5527989,  0.041625235, 0.17414952, -0.17907039, 0.33434446, -1.8855166, -8.7246452,
            -0.22587991, -0.35054478, 0.1027647,   0.19591429, -5.0122578
        ],
        [
            2.0794423,   -2.7220897,  0.062146382, 0.17412791, -0.17913654, 0.33510505, -1.8825649, -9.0289412,
            -0.22579929, -0.37050464, 0.098599659, 0.18902037, -5.1527704
        ],
        [
            2.3027882,   -3.0005167,  0.094292894, 0.17408396, -0.17923652, 0.33617473, -1.8797283, -9.5334369,
            -0.22568477, -0.40283896, 0.092552729, 0.17722861, -5.3940615
        ],
        [
            2.9957458,   -3.8351349,  0.18022608,  0.17395101, -0.1794296, 0.33831307, -1.8738564, -11.068858,
            -0.22544349, -0.50119067, 0.078933077, 0.15439589, -6.1059558
        ],
        [
            3.9121548,   -4.8932253,  0.27134274,  0.17384378, -0.17955246, 0.33960586, -1.871702, -13.049337,
            -0.22529632, -0.63076272, 0.067774263, 0.13429602, -7.0364037
        ],
        [
            4.6053029,   -5.6715606,  0.32878795,  0.17380418, -0.17960077, 0.34004212, -1.8719528, -14.522549,
            -0.22525017, -0.72855218, 0.062049745, 0.12350286, -7.7350496
        ],
        [
            5.2986443,   -6.4372784,  0.37918899,  0.17378393, -0.17962209, 0.34025808, -1.8717965, -15.981943,
            -0.22522573, -0.82626833, 0.057663868, 0.11503936, -8.4314387
        ],
        [
            6.214758,    -7.4348161,  0.43761403,  0.1737718,  -0.17959899, 0.34038817, -1.8717918, -17.895003,
            -0.22519734, -0.95845339, 0.053169141, 0.10622642, -9.3495624
        ],
        [
            6.9078641,   -8.1814622, 0.47699423,  0.17376781, -0.1795514, 0.34043157, -1.8717919, -19.334001,
            -0.22495502, -1.0349814, 0.050430397, 0.1008058,  -10.043398
        ],
        [
            7.6009213,  -8.9226752, 0.51304092,  0.17376585,  -0.17979893, 0.34045328, -1.8717923, -20.767465,
            -0.2251956, -1.1465508, 0.048095994, 0.096164854, -10.736843
        ],
        [
            8.5173628,   -9.8961868, 0.55649649,  0.17376468,  -0.17967332, 0.34046631, -1.8717927, -22.656408,
            -0.22509722, -1.2676684, 0.045472727, 0.090934753, -11.653532
        ],
        [
            9.2103822,   -10.628283, 0.58669214,  0.1737643,   -0.1796442, 0.34047065, -1.8717928, -24.080842,
            -0.21435015, -1.2078324, 0.043759552, 0.087513801, -12.346639
        ],
        [
            10.819812,   -12.317885, 0.6497023,   0.17376397,  -0.17964423, 0.34047413, -1.8717929, -27.378578,
            -0.20484253, -1.3203289, 0.040440298, 0.080879552, -13.956143
        ],
        [
            11.513022,   -13.041931, 0.67429738,  0.17376392,  -0.17964456, 0.34047456, -1.8717929, -28.795385,
            -0.20079021, -1.3663789, 0.039229909, 0.078459299, -14.649363
        ],
        [
            13.122442,   -14.716216, 0.72664999,  0.17376385,  -0.17964427, 0.34047491, -1.8717929, -32.078242,
            -0.19181467, -1.4685638, 0.036797196, 0.073594291, -16.258791
        ],
        [
            13.81565,    -15.434913, 0.74744114, 0.17376383, -0.17964451, 0.34047496, -1.871793, -33.489848,
            -0.18817619, -1.5107677, 0.03588184, 0.07176363, -16.952
        ],
        [
            16.118273,   -17.813792, 0.8103359,  0.17376377,  -0.1796452, 0.34047499, -1.8717932, -38.170561,
            -0.17711979, -1.6442957, 0.03327357, 0.066547136, -19.254624
        ],
        [
            18.420692,   -20.182207, 0.86556098,  0.17376369,  -0.1796446, 0.34047554, -1.8717946, -42.840828,
            -0.16753497, -1.7691179, 0.031166265, 0.062332524, -21.557044
        ]
    ],

    tail_shock_table => [
        [ 4.4342523, -0.80308,    -0.0405023,   -0.0405023 ],
        [ 4.5913146, -0.81906896, -0.048142977, -0.031931802 ],
        [ 5.2913144, -0.887203,   -0.054450553, -0.020470343 ],
        [ 6.2116557, -0.97044571, -0.056235229, -0.013343132 ],
        [ 6.9062517, -1.0294719,  -0.056203193, -0.010075415 ],
        [ 7.6000858, -1.0857374,  -0.055660291, -0.0077755477 ],
        [ 8.5170139, -1.1565708,  -0.05455871,  -0.0056486156 ],
        [ 9.2102024, -1.2078102,  -0.053588323, -0.0044916954 ],
        [ 10.819774, -1.3203244,  -0.051210791, -0.0027076226 ],
        [ 11.513002, -1.3663766,  -0.050197629, -0.0021910572 ],
        [ 13.122438, -1.4685633,  -0.047953682, -0.0013423861 ],
        [ 13.815648, -1.5107675,  -0.047044055, -0.0010825111 ],
        [ 16.118272, -1.6443883,  -0.044282551, -0.00049718133 ]
    ],

    blast_info => [
        0.07653254,  2.6554192,   0.29657278,  -0.52840794,  3.3050029, 0.00062328477,
        1.7425461,   -0.40656279, 0.35462092,  0.63044949,   0.6352101, -0.37319608,
        0.088497324, 0.44710124,  0.043440943, -0.044911029, 83.486,    -0.80308,
        660.73062,   -0.99451734
    ],

};
1;
