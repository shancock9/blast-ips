package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=1.12
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S1x12'} = {

    table_name  => 'S1x12',
    symmetry    => 2,
    gamma       => 1.12,
    data_source => 'S8000_G1x12/moc2_from_moc_r50/',

    shock_table_info => [ 2, 1.12, 1.86e-06, 8000, 2e-07, 10, 9.21e-07, 4.4e-07, 5e-07, 1692.7, -0.45529 ],

    shock_table => [
        [ -4.54674842,  10.70016224,   -2.99992855, -4.54870245, 0.99706616 ],
        [ -4.21810139,  9.71426079,    -2.99981332, -4.2213024,  0.99519107 ],
        [ -3.8586612,   8.63606313,    -2.99944791, -3.86415553, 0.99173722 ],
        [ -3.61242787,  7.89756674,    -2.99883092, -3.62038608, 0.98801931 ],
        [ -3.36272331,  7.14888729,    -2.99753039, -3.37431603, 0.98252287 ],
        [ -3.14862684,  6.50733831,    -2.99531964, -3.16463966, 0.97582197 ],
        [ -2.95468591,  5.92674435,    -2.99166647, -2.97615161, 0.96753618 ],
        [ -2.79971231,  5.46346411,    -2.98680279, -2.82685183, 0.95889978 ],
        [ -2.65190513,  5.0224871,     -2.97961251, -2.68585783, 0.94852031 ],
        [ -2.51652719,  4.619738,      -2.96975684, -2.55821968, 0.93673378 ],
        [ -2.38809171,  4.2391317,     -2.95626093, -2.43875745, 0.9230986 ],
        [ -2.2677608,   3.88440947,    -2.93858283, -2.32857678, 0.90774584 ],
        [ -2.15085875,  3.54217958,    -2.91526146, -2.22346513, 0.89005157 ],
        [ -2.03440461,  3.20440025,    -2.88440326, -2.12098978, 0.86933056 ],
        [ -1.91447923,  2.86086574,    -2.84296455, -2.0181911,  0.84441752 ],
        [ -1.78055697,  2.48395911,    -2.78329245, -1.90721992, 0.8120072 ],
        [ -1.66040765,  2.15341784,    -2.71682746, -1.81161704, 0.7787431 ],
        [ -1.54953118,  1.8560854,     -2.64485257, -1.72713113, 0.74471684 ],
        [ -1.44483443,  1.5831073,     -2.56856803, -1.65095823, 0.71002265 ],
        [ -1.34430366,  1.32883087,    -2.48927487, -1.5813339,  0.67485671 ],
        [ -1.24610146,  1.08834953,    -2.40792026, -1.51679922, 0.63932132 ],
        [ -1.14910496,  0.8587715,     -2.325656,   -1.45651408, 0.60367684 ],
        [ -1.05213449,  0.63725026,    -2.24330394, -1.39970287, 0.56810409 ],
        [ -0.95426822,  0.42171871,    -2.16167677, -1.34583766, 0.53283884 ],
        [ -0.85470438,  0.21051363,    -2.08150639, -1.29452565, 0.49813149 ],
        [ -0.75258305,  0.00197741,    -2.00333438, -1.24540472, 0.46418894 ],
        [ -0.64701339,  -0.20546771,   -1.9275795,  -1.19816196, 0.43119254 ],
        [ -0.53727776,  -0.41293784,   -1.85470742, -1.15261487, 0.39936118 ],
        [ -0.42258919,  -0.62159434,   -1.78505816, -1.10858844, 0.36886809 ],
        [ -0.30199196,  -0.83280296,   -1.71882118, -1.0658859,  0.33982247 ],
        [ -0.17443434,  -1.04797528,   -1.65611858, -1.02432797, 0.3123014 ],
        [ -0.03894442,  -1.26828072,   -1.59709992, -0.98380651, 0.28638809 ],
        [ 0.10570456,   -1.49520926,   -1.5417875,  -0.94417663, 0.26210404 ],
        [ 0.26079013,   -1.7302217,    -1.49019266, -0.90532615, 0.23945907 ],
        [ 0.42293664,   -1.9679777,    -1.44356185, -0.86819834, 0.21899883 ],
        [ 0.59232658,   -2.20885907,   -1.40159573, -0.83270064, 0.20058641 ],
        [ 0.76940228,   -2.45362473,   -1.36389403, -0.79868387, 0.18403819 ],
        [ 0.95617681,   -2.70509897,   -1.32979853, -0.76574514, 0.1690566 ],
        [ 1.15325619,   -2.96406297,   -1.29902541, -0.73379702, 0.15550859 ],
        [ 1.36225434,   -3.23256889,   -1.27117295, -0.70261465, 0.14320969 ],
        [ 1.58460236,   -3.51233104,   -1.2459431,  -0.67204864, 0.13202157 ],
        [ 1.82109698,   -3.80422027,   -1.22314381, -0.64206009, 0.12185399 ],
        [ 2.07456961,   -4.1115559,    -1.20243243, -0.61238374, 0.11255017 ],
        [ 2.34482851,   -4.43392949,   -1.18375644, -0.58314017, 0.10408487 ],
        [ 2.63501342,   -4.7749122,    -1.16683162, -0.55409138, 0.09632931 ],
        [ 2.94646609,   -5.13587415,   -1.1515337,  -0.52522442, 0.08922847 ],
        [ 3.28125136,   -5.51900969,   -1.13771193, -0.49647092, 0.08271654 ],
        [ 3.64162538,   -5.92669605,   -1.12522828, -0.46776831, 0.07673476 ],
        [ 4.03043632,   -6.36193873,   -1.1139459,  -0.4390322,  0.07122571 ],
        [ 4.45032376,   -6.8274671,    -1.10375443, -0.41022002, 0.06614546 ],
        [ 4.90452037,   -7.32663411,   -1.09454392, -0.38127104, 0.06145092 ],
        [ 5.39665212,   -7.86318231,   -1.08621366, -0.35212672, 0.05710387 ],
        [ 5.93093911,   -8.44145587,   -1.07866983, -0.32272278, 0.05306944 ],
        [ 6.5121963,    -9.06639154,   -1.07182735, -0.29299486, 0.04931686 ],
        [ 7.14503444,   -9.74266028,   -1.06561759, -0.26291944, 0.04582365 ],
        [ 7.83566102,   -10.47659435,  -1.05996883, -0.23242723, 0.04256457 ],
        [ 8.59008092,   -11.27425682,  -1.05482496, -0.20149345, 0.03952198 ],
        [ 9.4152972,    -12.14272127,  -1.05013362, -0.17008317, 0.03667897 ],
        [ 10.31871166,  -13.08943392,  -1.0458504,  -0.13817878, 0.03402186 ],
        [ 11.30792526,  -14.12200984,  -1.04193794, -0.10578415, 0.03153971 ],
        [ 12.39213852,  -15.24969456,  -1.03835979, -0.07287915, 0.02922052 ],
        [ 13.58015166,  -16.48127741,  -1.03508759, -0.03948501, 0.02705597 ],
        [ 14.88216478,  -17.82696755,  -1.03209398, -0.00560757, 0.02503696 ],
        [ 16.30857791,  -19.29714836,  -1.02935576, 0.02872761,  0.02315592 ],
        [ 17.87159106,  -20.90402808,  -1.02685026, 0.06351458,  0.0214044 ],
        [ 19.58340423,  -22.65978354,  -1.02455854, 0.09872257,  0.01977555 ],
        [ 21.41554062,  -24.53497887,  -1.02250602, 0.13355922,  0.01829359 ],
        [ 23.3853882,   -26.54728106,  -1.02065219, 0.16822209,  0.01693586 ],
        [ 25.40526142,  -28.60720257,  -1.01904581, 0.20119569,  0.01574325 ],
        [ 27.68057066,  -30.92405362,  -1.01751271, 0.23567129,  0.01459084 ],
        [ 30.25524401,  -33.54188147,  -1.01605205, 0.27176792,  0.01347914 ],
        [ 33.18299445,  -36.51453332,  -1.01466301, 0.30962006,  0.0124086 ],
        [ 36.17427061,  -39.54783951,  -1.01347267, 0.34531297,  0.01148069 ],
        [ 39.56428539,  -42.98154726,  -1.01233829, 0.3826757,   0.01058683 ],
        [ 42.97082858,  -46.42843849,  -1.01137637, 0.41740261,  0.00982113 ],
        [ 46.81204981,  -50.31154473,  -1.01045764, 0.45367197,  0.00908279 ],
        [ 51.16386599,  -54.70690841,  -1.00958146, 0.49160929,  0.00837203 ],
        [ 55.463258,    -59.04587917,  -1.00884933, 0.52628348,  0.0077729 ],
        [ 60.30059046,  -63.92427821,  -1.00814904, 0.56244967,  0.00719518 ],
        [ 65.76812312,  -69.43448508,  -1.00748015, 0.60022756,  0.00663898 ],
        [ 71.97866215,  -75.68944094,  -1.00684251, 0.63975217,  0.00610445 ],
        [ 77.99791841,  -81.74827112,  -1.00632044, 0.67513591,  0.0056636 ],
        [ 84.77173066,  -88.56315812,  -1.0058208,  0.7120235,   0.00523883 ],
        [ 92.43009144,  -96.26421631,  -1.00534336, 0.75053665,  0.00483019 ],
        [ 101.13210402, -105.01068639, -1.00488792, 0.79081279,  0.00443778 ],
        [ 109.32133769, -113.23842425, -1.00452503, 0.82583469,  0.00412319 ],
        [ 118.50949281, -122.46651516, -1.00417715, 0.86229014,  0.00381994 ],
        [ 128.86417527, -132.86267867, -1.00384413, 0.90029284,  0.00352807 ],
        [ 140.58953091, -144.63118852, -1.00352588, 0.93997059,  0.00324761 ],
        [ 153.93621915, -158.0228486,  -1.00322227, 0.98146791,  0.00297859 ],
        [ 165.98712608, -170.11115129, -1.00298984, 1.01608543,  0.00277164 ],
        [ 179.46665593, -183.62943947, -1.00276663, 1.05206667,  0.00257204 ],
        [ 194.60797075, -198.81098087, -1.0025526,  1.08951738,  0.00237981 ],
        [ 211.69365635, -215.9384799,  -1.00234768, 1.12855607,  0.00219497 ]
    ],

    energy_table => [
        [ -4.54674842,  0.30972354, 0.09948411, 0.30972354, 0.09948411,    5.0589581,     -1.6316005 ],
        [ -4.21810139,  0.34419229, 0.11045412, 0.34419229, 0.11045412,    4.518254,      -1.6628617 ],
        [ -3.8586612,   0.38622998, 0.12365696, 0.38622998, 0.12365696,    3.9136283,     -1.707307 ],
        [ -3.61242787,  0.41786298, 0.13335521, 0.41786298, 0.13335521,    3.4889858,     -1.746523 ],
        [ -3.36272331,  0.45243058, 0.14355055, 0.45243058, 0.14355055,    3.0473078,     -1.7981055 ],
        [ -3.14862684,  0.48410517, 0.15231334, 0.48410517, 0.15231334,    2.6569607,     -1.8560957 ],
        [ -2.95468591,  0.51439087, 0.15991614, 0.51439087, 0.15991614,    2.2912138,     -1.9204724 ],
        [ -2.79971231,  0.53961012, 0.16543747, 0.53961012, 0.16543747,    1.989306,      -1.9814959 ],
        [ -2.65190513,  0.56440518, 0.16990815, 0.56440518, 0.16990815,    1.6917222,     -2.0502384 ],
        [ -2.51652719,  0.58762853, 0.17299,    0.58762853, 0.17299,       1.409588,      -2.1233423 ],
        [ -2.38809171,  0.60997019, 0.17468735, 0.60997019, 0.17468735,    1.1320876,     -2.210349 ],
        [ -2.2677608,   0.63101757, 0.17488871, 0.63101757, 0.17488871,    0.86050811,    -2.3153359 ],
        [ -2.15085875,  0.65139877, 0.17352176, 0.65139877, 0.17352176,    0.58320838,    -2.4153678 ],
        [ -2.03440461,  0.67144194, 0.17039118, 0.67144194, 0.17039118,    0.29690579,    -2.5049933 ],
        [ -1.91447923,  0.69158245, 0.16513991, 0.69158245, 0.16513991,    -0.0092483214, -2.6369342 ],
        [ -1.78055697,  0.71316676, 0.15676245, 0.71316676, 0.15676245,    -0.37496428,   -2.7580146 ],
        [ -1.66040765,  0.73144035, 0.14710077, 0.73144035, 0.14710077,    -0.70927054,   -2.7648375 ],
        [ -1.54953118,  0.74718286, 0.1366468,  0.74718286, 0.1366468,     -1.0140266,    -2.6622048 ],
        [ -1.44483443,  0.76092636, 0.1257595,  0.76092636, 0.1257595,     -1.2842087,    -2.3721311 ],
        [ -1.34430366,  0.77301827, 0.11474098, 0.77301827, 0.11474098,    -1.5025557,    -1.9563836 ],
        [ -1.24610146,  0.7837489,  0.10380013, 0.7837489,  0.10380013,    -1.6739989,    -1.6644584 ],
        [ -1.14910496,  0.793298,   0.09314547, 0.793298,   0.09314547,    -1.8276508,    -1.5638817 ],
        [ -1.05213449,  0.80183018, 0.0829202,  0.80183018, 0.0829202,     -1.9773414,    -1.573251 ],
        [ -0.95426822,  0.80946601, 0.07324912, 0.80946601, 0.07324912,    -2.1342312,    -1.632393 ],
        [ -0.85470438,  0.81630267, 0.06422955, 0.81630267, 0.06422955,    -2.2997252,    -1.6920475 ],
        [ -0.75258305,  0.82242931, 0.05591999, 0.82242931, 0.05591999,    -2.4756467,    -1.7493419 ],
        [ -0.64701339,  0.82792417, 0.04835046, 0.82792417, 0.04835046,    -2.6632348,    -1.7998958 ],
        [ -0.53727776,  0.83284675, 0.04154029, 0.83284675, 0.04154029,    -2.8633691,    -1.8430992 ],
        [ -0.42258919,  0.83725374, 0.03548185, 0.83725374, 0.03548185,    -3.0770665,    -1.8791858 ],
        [ -0.30199196,  0.84120095, 0.03014309, 0.84120095, 0.03014309,    -3.3057076,    -1.9087816 ],
        [ -0.17443434,  0.84473866, 0.02547952, 0.84473866, 0.02547952,    -3.5509243,    -1.9326383 ],
        [ -0.03894442,  0.84790785, 0.02144374, 0.84790785, 0.02144374,    -3.8142508,    -1.9514944 ],
        [ 0.10570456,   0.85074961, 0.01797718, 0.85074961, 0.01797718,    -4.0977649,    -1.9661054 ],
        [ 0.26079013,   0.85329945, 0.01502128, 0.85329945, 0.01502128,    -4.4036922,    -1.97718 ],
        [ 0.42293664,   0.85552943, 0.01258108, 0.85552943, 0.01258108,    -4.7250555,    -1.9851534 ],
        [ 0.59232658,   0.85748366, 0.01057259, 0.85748366, 0.01057259,    -5.0618907,    -1.9907499 ],
        [ 0.76940228,   0.85920359, 0.00891911, 0.85920359, 0.00891911,    -5.414817,     -1.9945979 ],
        [ 0.95617681,   0.86073611, 0.00754596, 0.86073611, 0.00754596,    -5.7876558,    -1.997175 ],
        [ 1.15325619,   0.86210639, 0.00640536, 0.86210639, 0.00640536,    -6.1814611,    -1.9988079 ],
        [ 1.36225434,   0.86334161, 0.00545294, 0.86334161, 0.00545294,    -6.5993415,    -1.9997966 ],
        [ 1.58460236,   0.86446183, 0.00465499, 0.86446183, 0.00465499,    -7.0440754,    -2.0003552 ],
        [ 1.82109698,   0.86548049, 0.00398618, 0.86548049, 0.00398618,    -7.5171952,    -2.0006165 ],
        [ 2.07456961,   0.86641642, 0.00342103, 0.86641642, 0.00342103,    -8.0243143,    -2.0006731 ],
        [ 2.34482851,   0.8672742,  0.00294553, 0.8672742,  0.00294553,    -8.5650103,    -2.0006284 ],
        [ 2.63501342,   0.86806813, 0.00254211, 0.86806813, 0.00254211,    -9.1455528,    -2.0005567 ],
        [ 2.94646609,   0.86880446, 0.00219961, 0.86880446, 0.00219961,    -9.7686187,    -2.0004577 ],
        [ 3.28125136,   0.86949014, 0.00190794, 0.86949014, 0.00190794,    -10.438322,    -2.0003466 ],
        [ 3.64162538,   0.8701311,  0.0016588,  0.8701311,  0.0016588,     -11.159176,    -2.0002535 ],
        [ 4.03043632,   0.87073294, 0.00144514, 0.87073294, 0.00144514,    -11.936879,    -2.0001771 ],
        [ 4.45032376,   0.87129971, 0.00126145, 0.87129971, 0.00126145,    -12.776714,    -2.000118 ],
        [ 4.90452037,   0.87183531, 0.00110295, 0.87183531, 0.00110295,    -13.685149,    -2.0000754 ],
        [ 5.39665212,   0.87234308, 0.00096567, 0.87234308, 0.00096567,    -14.66944,     -2.0000463 ],
        [ 5.93093911,   0.87282604, 0.0008464,  0.87282604, 0.0008464,     -15.738033,    -2.0000272 ],
        [ 6.5121963,    0.87328681, 0.00073813, 0.87328681, 0.00073813,    -16.900558,    -2.000015 ],
        [ 7.14503444,   0.8737271,  0.00065209, 0.8737271,  0.00065209,    -18.166241,    -2.000008 ],
        [ 7.83566102,   0.87414901, 0.00057264, 0.87416755, 0.00060938592, -19.547498,    -2.000004 ],
        [ 8.59008092,   0.87455378, 0.00050297, 0.87460989, 0.00058535514, -21.05634,     -2.000002 ],
        [ 9.4152972,    0.87494266, 0.00044177, 0.87508725, 0.00057124613, -22.706773,    -2.0000009 ],
        [ 10.31871166,  0.87531657, 0.00038796, 0.8755983,  0.00055717154, -24.513603,    -2.0000004 ],
        [ 11.30792526,  0.87567608, 0.00034063, 0.87613652, 0.0005327895,  -26.49203,     -2.0000002 ],
        [ 12.39213852,  0.87602199, 0.00029898, 0.87669727, 0.00050294816, -28.660457,    -2.0000001 ],
        [ 13.58015166,  0.87635463, 0.00026235, 0.8773117,  0.00054545211, -31.036483,    -2.0000001 ],
        [ 14.88216478,  0.87667447, 0.00023013, 0.87797606, 0.00047759302, -33.64051,     -2.0000001 ],
        [ 16.30857791,  0.8769818,  0.00020181, 0.87861333, 0.00041813803, -36.493336,    -2.0000001 ],
        [ 17.87159106,  0.87727707, 0.00017693, 0.87922465, 0.00036603602, -39.619363,    -2.0000001 ],
        [ 19.58340423,  0.87756055, 0.00015508, 0.87981072, 0.00032040108, -43.042989,    -2 ],
        [ 21.41554062,  0.87782683, 0.00013629, 0.88036056, 0.00028124828, -46.707262,    -2 ],
        [ 23.3853882,   0.87807867, 0.00011999, 0.88088,    0.00024735431, -50.646957,    -2 ],
        [ 25.40526142,  0.87830685, 0.00010641, 0.88135018, 0.00021915938, -54.686703,    -2 ],
        [ 27.68057066,  0.87853428, 9.394e-05,  0.88181841, 0.00019331813, -59.237322,    -2 ],
        [ 30.25524401,  0.87876092, 8.253e-05,  0.88228465, 0.00016971021, -64.386669,    -2 ],
        [ 33.18299445,  0.87898675, 7.213e-05,  0.88274887, 0.00014821774, -70.242169,    -2 ],
        [ 36.17427061,  0.87918928, 6.359e-05,  0.88316491, 0.00013058019, -76.224722,    -2 ],
        [ 39.56428539,  0.87939111, 5.578e-05,  0.88357928, 0.00011447298, -83.004751,    -2 ],
        [ 42.97082858,  0.8795699,  4.942e-05,  0.88394617, 0.00010137452, -89.817838,    -2 ],
        [ 46.81204981,  0.87974812, 4.358e-05,  0.8843117,  8.9365182e-05, -97.50028,     -2 ],
        [ 51.16386599,  0.87992573, 3.825e-05,  0.88467584, 7.8388781e-05, -106.20391,    -2 ],
        [ 55.463258,    0.88008064, 3.396e-05,  0.88499333, 6.958951e-05,  -114.8027,     -2 ],
        [ 60.30059046,  0.88023508, 3.003e-05,  0.88530972, 6.150421e-05,  -124.47736,    -2 ],
        [ 65.76812312,  0.88038902, 2.642e-05,  0.88562503, 5.4097293e-05, -135.41243,    -2 ],
        [ 71.97866215,  0.88054247, 2.312e-05,  0.88593923, 4.7335904e-05, -147.8335,     -2 ],
        [ 77.99791841,  0.88067359, 2.054e-05,  0.88620765, 4.2027719e-05, -159.87202,    -2 ],
        [ 84.77173066,  0.88080434, 1.815e-05,  0.88647524, 3.7147761e-05, -173.41964,    -2 ],
        [ 92.43009144,  0.8809347,  1.597e-05,  0.886742,   3.2675418e-05, -188.73636,    -2 ],
        [ 101.13210402, 0.88106468, 1.398e-05,  0.88700791, 2.8590593e-05, -206.14039,    -2 ],
        [ 109.32133769, 0.88117269, 1.245e-05,  0.88722885, 2.5468207e-05, -222.51886,    -2 ],
        [ 118.50949281, 0.88128042, 1.105e-05,  0.88744919, 2.2589837e-05, -240.89517,    -2 ],
        [ 128.86417527, 0.88138787, 9.75e-06,   0.88766893, 1.9944176e-05, -261.60453,    -2 ],
        [ 140.58953091, 0.88149503, 8.57e-06,   0.88788806, 1.7520396e-05, -285.05524,    -2 ],
        [ 153.93621915, 0.88160191, 7.49e-06,   0.88810658, 1.5307607e-05, -311.74862,    -2 ],
        [ 165.98712608, 0.8816872,  6.69e-06,   0.88828095, 1.3682093e-05, -335.85043,    -2 ],
        [ 179.46665593, 0.88177231, 5.96e-06,   0.88845492, 1.2179323e-05, -362.80949,    -2 ],
        [ 194.60797075, 0.88185722, 5.28e-06,   0.88862849, 1.0794056e-05, -393.09212,    -2 ],
        [ 211.69365635, 0.88194195, 4.66e-06,   0.88880166, 9.5209755e-06, -427.26349,    -2 ]
    ],

    impulse_table => [
        [
            -5.5849819,  13.81484,    -5.5853933,    0.011640591, -0.00028396791, -0.12994161,
            -0.98413055, -0.52861571, -0.0007278716, -0.28225447, 0.095622386,    0.75151449,
            -1.6411506
        ],
        [
            -4.9618452,  11.945436,   -4.9628931,    0.015079925, -0.00052953432, -0.12669543,
            -0.98088436, -0.84019107, -0.0013573118, -0.27900829, 0.095621829,    0.69641645,
            -1.5889091
        ],
        [
            -4.6051702,  10.875424,  -4.6069602,    0.017312731, -0.0007564776, -0.12369543,
            -0.97788436, -1.0185397, -0.0019390168, -0.27600829, 0.095620813,   0.65955363,
            -1.5636285
        ],
        [
            -4.1997051,  9.6590755, -4.2029958,    0.020005521, -0.0011347164, -0.11869543,
            -0.97288436, -1.221311, -0.0029085253, -0.27100829, 0.095617374,   0.61221628,
            -1.5412118
        ],
        [
            -3.912023,   8.7961213,  -3.9170935,    0.021934812, -0.0015129552, -0.11369543,
            -0.96788436, -1.3652271, -0.0038780337, -0.26600829, 0.095610733,   0.57474484,
            -1.5310596
        ],
        [
            -3.5065579,  7.5801022, -3.5158911,    0.024473944, -0.0022694328, -0.10369543,
            -0.95788437, -1.568268, -0.0058170505, -0.25600829, 0.095583492,   0.51591622,
            -1.5285349
        ],
        [
            -2.9957323,  6.0495612,  -3.0159061,    0.026733609, -0.003782388, -0.083695433,
            -0.93788437, -1.8252713, -0.0096950842, -0.23600829, 0.09544344,   0.43116778,
            -1.5557135
        ],
        [
            -2.3025851,  3.9868442,  -2.3602713,   0.025631459, -0.007564776, -0.033695435,
            -0.88788437, -2.1860639, -0.019390168, -0.1860083,  0.094224895,  0.29986803,
            -1.6928729
        ],
        [
            -1.89712,    2.8115736,  -2.0035666,   0.021304804, -0.011347164, 0.016304561,
            -0.83788437, -2.4276632, -0.029085253, -0.1360083,  0.09117464,   0.2215516,
            -1.8693898
        ],
        [
            -1.609438,   2.0157455,  -1.7723111,   0.017100986, -0.01512948, 0.066300766,
            -0.78788661, -2.6492737, -0.038780259, -0.0860083,  0.086079652, 0.17115449,
            -2.0587625
        ],
        [
            -1.3862944,  1.4340686,  -1.6099849,   0.015100848,  -0.018876921, 0.11490696,
            -0.73874471, -2.9002679, -0.048435684, -0.036332096, 0.079363243,  0.13756682,
            -2.2459659
        ],
        [
            -1.2039728,  0.98765666, -1.4901908,   0.015772019,    -0.021364378, 0.14109733,
            -0.70838187, -3.2369493, -0.055514609, -0.00086984822, 0.071789086,  0.11447607,
            -2.4223404
        ],
        [
            -1.0498222,  0.63206528, -1.3983902,   0.016897699,  -0.021881922, 0.15049157,
            -0.69626697, -3.6190127, -0.056704331, 0.0044297844, 0.06411674,   0.098090312,
            -2.5843332
        ],
        [
            -0.91629078, 0.34021378, -1.325857,    0.01778982,   -0.022013383, 0.15632313,
            -0.68882407, -3.9719045, -0.056653437, 0.0034228534, 0.056914737,  0.086094087,
            -2.731543
        ],
        [
            -0.79850774, 0.094770811, -1.2670662,  0.018447826,  -0.02207693, 0.1607172,
            -0.68343354, -4.2880584,  -0.05645574, 0.0013759678, 0.050521645, 0.077051276,
            -2.8650258
        ],
        [
            -0.69314723, -0.11579673, -1.2183791,   0.018932106,   -0.022120154, 0.16424695,
            -0.67928013, -4.5716767,  -0.056253066, -0.0011434513, 0.045080989,  0.070052415,
            -2.9863134
        ],
        [
            -0.59783705, -0.29943253, -1.1773182,   0.019293085,   -0.022153918, 0.16717014,
            -0.67597246, -4.8277211,  -0.056068051, -0.0037945565, 0.040597409,  0.064507296,
            -3.0969624
        ],
        [
            -0.51082567, -0.46177772, -1.1421476,   0.019566529,   -0.022181797, 0.16963795,
            -0.67328059, -5.0604842,  -0.055902828, -0.0063006671, 0.036989861,  0.060022883,
            -3.1983844
        ],
        [
            -0.43078296, -0.60694852, -1.1116194,   0.019777051,   -0.022205424, 0.17175112,
            -0.67105426, -5.2734815,  -0.055755546, -0.0088728992, 0.034119607,  0.056330675,
            -3.2918011
        ],
        [
            -0.40906169, -0.645689,  -1.1036217,   0.019828121,   -0.022211533, 0.17229986,
            -0.67048884, -5.3310711, -0.055716462, -0.0096321618, 0.03340734,   0.055394249,
            -3.3171488
        ],
        [
            -0.28768212, -0.85734599, -1.0610464,   0.02007208,   -0.022243456, 0.17518228,
            -0.66760956, -5.6510481,  -0.055505438, -0.013792927, 0.029866528,  0.050624523,
            -3.4586023
        ],
        [
            -0.2231436,  -0.96674905, -1.0397849,   0.02017682,   -0.022258996, 0.17659405,
            -0.66625928, -5.8198342,  -0.055398578, -0.016187803, 0.028246022,  0.048377678,
            -3.5335958
        ],
        [
            -0.16251897, -1.0676758, -1.0206212,   0.020261876,  -0.022272744, 0.17784823,
            -0.6650964,  -5.9774897, -0.055301631, -0.018445997, 0.02686387,   0.046428932,
            -3.6038508
        ],
        [
            -0.10536056, -1.1612874, -1.0032313,   0.020331668, -0.022284996, 0.17896975,
            -0.66408709, -6.125321,  -0.055213358, -0.02078581, 0.025671133,  0.044722748,
            -3.6698935
        ],
        [
            -0.051293339, -1.248527,  -0.98735682,  0.020389467,  -0.022295975, 0.1799784,
            -0.66320557,  -6.2644255, -0.055132635, -0.022973527, 0.024631328,  0.043216363,
            -3.732172
        ],
        [
            -4.4479394e-08, -1.3301728, -0.97278767,  0.020437744,  -0.022305874, 0.18089037,
            -0.66243115,    -6.3957325, -0.055058574, -0.024873974, 0.023716713,  0.041876382,
            -3.7910709
        ],
        [
            0.095310135, -1.4791639, -0.94690954,  0.020512819,  -0.022323014, 0.18247483,
            -0.66114066, -6.6380141, -0.054927423, -0.028748369, 0.022181834,  0.03959519,
            -3.9000085
        ],
        [
            0.18232151,  -1.6123128, -0.92454427,  0.020567399,  -0.022337337, 0.18380389,
            -0.66011587, -6.8572816, -0.054814966, -0.032470125, 0.020943484,  0.037723881,
            -3.9988631
        ],
        [
            0.26236422,  -1.732567,  -0.90494939,  0.020607976,  -0.022349491, 0.18493446,
            -0.65928918, -7.0574004, -0.054717548, -0.036167963, 0.019922295,  0.036159054,
            -4.0892845
        ],
        [
            0.40546506,  -1.9427155, -0.87204253, 0.020662381,  -0.022369017, 0.18675447,
            -0.65805347, -7.4114954, -0.05455731, -0.042468574, 0.018333901,  0.033683952,
            -4.2497024
        ],
        [
            0.53062821,  -2.1219382, -0.84527161,  0.020695427, -0.022384042, 0.18815515,
            -0.65719016, -7.7174973, -0.054431075, -0.04820994, 0.017151654,  0.031807148,
            -4.3887304
        ],
        [
            0.69314714,  -2.3490354, -0.81297473, 0.020723982,  -0.02240108, 0.18973982,
            -0.65631668, -8.1099748, -0.05428534, -0.055856095, 0.015847314, 0.029699518,
            -4.5675421
        ],
        [
            0.99325173,  -2.7542864, -0.75952782,  0.020749935, -0.022426433, 0.19208379,
            -0.65524981, -8.8216295, -0.054064177, -0.06984527, 0.01395559,   0.026566453,
            -4.8930709
        ],
        [
            1.0986122,   -2.8928609, -0.74239078,  0.020753978,  -0.022433765, 0.19275678,
            -0.65499932, -9.0678403, -0.053999483, -0.074850525, 0.013413394,  0.025650089,
            -5.0060629
        ],
        [
            1.2089603,   -3.0362037, -0.72523171,  0.020756457, -0.022440702, 0.19339087,
            -0.65478631, -9.3238697, -0.053938144, -0.0807398,  0.012899385,  0.024773116,
            -5.1237554
        ],
        [
            1.3862943,   -3.2630926, -0.69918749, 0.020757845, -0.022450456, 0.19427629,
            -0.65453449, -9.7316739, -0.05385169, -0.08912172, 0.012171369,  0.023516438,
            -5.3116185
        ],
        [
            1.6094379,   -3.5432427, -0.66878405,  0.020756797,  -0.022460603, 0.19519123,
            -0.65432753, -10.239051, -0.053761527, -0.099917886, 0.01139548,   0.022156942,
            -5.5460287
        ],
        [
            1.7917594,   -3.7682979, -0.64565215, 0.020754696, -0.022467374, 0.19580212,
            -0.65422202, -10.649363, -0.05370083, -0.10978784, 0.01085433,   0.021195466,
            -5.7361302
        ],
        [
            1.9459101,   -3.9562126, -0.62715225,  0.020752584, -0.022472273, 0.19623928,
            -0.65416288, -10.993618, -0.053657041, -0.11628067, 0.010449831,  0.020469182,
            -5.8959944
        ],
        [
            2.0794415,   -4.1174131, -0.61183581,  0.020750634, -0.022475933, 0.19656738,
            -0.65413028, -11.290036, -0.053624311, -0.1243969,  0.010132774,  0.019895129,
            -6.0339098
        ],
        [
            2.3027863,    -4.3841053, -0.58754184,  0.02074756, -0.022481325, 0.19702858, -0.65409177, -11.782476,
            -0.053578342, -0.1345868, 0.0096608706, 0.01887754, -6.2706049
        ],
        [
            2.9958012,   -5.192631,  -0.52084769,  0.02073962,  -0.022492856, 0.19795516,
            -0.65412333, -13.288679, -0.053486717, -0.16947652, 0.0085371024, 0.01687531,
            -6.9730713
        ],
        [
            3.9123293,   -6.2301859, -0.44753705,  0.020733913, -0.022500041, 0.19851389,
            -0.65419931, -15.24464,  -0.053431672, -0.21574116, 0.007536404,  0.015000645,
            -7.8962953
        ],
        [
            4.6052014,   -6.9981539, -0.40010706,  0.020731892, -0.022502168, 0.1986987,
            -0.65421323, -16.705253, -0.053412992, -0.25074906, 0.0069891927, 0.013944087,
            -8.5917165
        ],
        [
            5.2984307,   -7.7564176, -0.35777525,  0.020730872, -0.022503225, 0.19879121,
            -0.65422069, -18.155921, -0.053403638, -0.28578202, 0.0065529977, 0.0130895,
            -9.2863402
        ],
        [
            6.2147706,   -8.7471131, -0.30793378,  0.020730264, -0.022503857, 0.19884671,
            -0.65422532, -20.061577, -0.053398026, -0.33209187, 0.0060910409, 0.012175736,
            -10.203589
        ],
        [
            6.9078222,   -9.4896255, -0.27393555,  0.020730065, -0.022504066, 0.19886521,
            -0.65422689, -21.496103, -0.053396156, -0.36711668, 0.0058028924, 0.011602685,
            -10.896967
        ],
        [
            7.6010532,   -10.227708, -0.2425352,   0.020729967, -0.022504169, 0.19887447,
            -0.65422768, -22.926462, -0.053395221, -0.40214986, 0.0055535686, 0.011105617,
            -11.590371
        ],
        [
            8.5172793,   -11.197447, -0.20438062, 0.020729908, -0.022504232, 0.19888002,
            -0.65422815, -24.811326, -0.05339466, -0.44845151, 0.0052699377, 0.01053928,
            -12.506707
        ],
        [
            9.2104935,   -11.92754,  -0.17766265,  0.020729888, -0.022504249, 0.19888187,
            -0.65422831, -26.233917, -0.053394473, -0.48348282, 0.0050827509, 0.010165209,
            -13.19996
        ],
        [
            10.819919,   -13.613095, -0.12145856,  0.020729871, -0.022504265, 0.19888335,
            -0.65422844, -29.527565, -0.053349437, -0.5570087,  0.0047164335, 0.0094328101,
            -14.809418
        ],
        [
            11.512928,   -14.335535, -0.099366517, 0.020729867, -0.022504266, 0.19888354,
            -0.65422845, -30.942554, -0.053056682, -0.57578218, 0.0045817519, 0.0091634758,
            -15.502431
        ],
        [
            13.122547,   -16.007347, -0.052045318, 0.020729863, -0.02250427,  0.19888368,
            -0.65422847, -34.22311,  -0.051793139, -0.61738878, 0.0043092444, 0.0086184836,
            -17.112054
        ],
        [
            13.815554,   -16.724871, -0.033161986, 0.020729861, -0.022504276, 0.1988837,
            -0.65422847, -35.633331, -0.051125953, -0.63454596, 0.0042061864, 0.0084123704,
            -17.805061
        ],
        [
            16.118176,   -19.101125, 0.024296548,  0.020729857, -0.022491092, 0.19888371,
            -0.68881342, -40.311381, -0.048769355, -0.68881342, 0.0039108351, 0.0078216703,
            -20.107684
        ],
        [
            18.420795,   -21.467763, 0.075117193,  0.020729852, -0.022504253, 0.19888452,
            -0.65422767, -44.980031, -0.053394291, -0.94891498, 0.0036705402, 0.0073410807,
            -22.410303
        ]
    ],

    tail_shock_table => [
        [ 7.4343491, -0.45529,    -0.03116194,  -0.03116194 ],
        [ 7.6006608, -0.46076971, -0.038041529, -0.023824503 ],
        [ 8.5171162, -0.48995542, -0.045210147, -0.013189651 ],
        [ 9.2104098, -0.51099629, -0.046887477, -0.0092773221 ],
        [ 10.819901, -0.55700788, -0.047633455, -0.0042613213 ],
        [ 11.512919, -0.57578175, -0.047372055, -0.002984262 ],
        [ 13.122545, -0.61738867, -0.046243877, -0.0010475442 ],
        [ 13.815553, -0.63455682, -0.045650654, -0.00049959304 ]
    ],

    blast_info => [
        0.13369536, 0.98788457, 0.28606275, -0.1939018,  1.2080796,  0.0037986787, 3.0477731, -0.075647774,
        0.33252495, 0.47938916, 0.46827428, -0.13131982, 0.45149296, 0.095622371,  0.0185088, -0.020093086,
        1692.7,     -0.45529,   34734.53,   -0.54639768
    ],

};
1;