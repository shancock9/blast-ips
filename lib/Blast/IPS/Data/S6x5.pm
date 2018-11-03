package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=6.5
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S6x5'} = {

    table_name => 'S6x5',
    symmetry   => 2,
    gamma      => 6.5,

    shock_table_info => [ 2, 6.5, 8.27e-07, 8000, 3.78e-07, 4, 3.64e-08, 2.9e-07, 5e-07, 0.59028, -0.83858 ],

    shock_table => [
        [ -3.90224117,  12.00017786,   -2.99996806, -3.90354728, 0.99803958 ],
        [ -3.62566258,  11.17045316,   -2.99992677, -3.6276409,  0.99702964 ],
        [ -3.44997359,  10.64339974,   -2.99987595, -3.45254917, 0.99613178 ],
        [ -3.31302799,  10.23257975,   -2.99984886, -3.31619181, 0.995247 ],
        [ -2.88942605,  8.96190149,    -2.99947046, -2.89540666, 0.99100383 ],
        [ -2.62494077,  8.16866028,    -2.99882094, -2.63384567, 0.98658838 ],
        [ -2.35396093,  7.35621181,    -2.99735612, -2.36735836, 0.97978693 ],
        [ -2.14727559,  6.73691245,    -2.99509463, -2.16558014, 0.97233682 ],
        [ -1.96296451,  6.18518483,    -2.99150227, -1.98715429, 0.96337971 ],
        [ -1.80104846,  5.70120424,    -2.98625597, -1.83196354, 0.95312485 ],
        [ -1.65906399,  5.27767572,    -2.97909976, -1.69741068, 0.94178324 ],
        [ -1.523553,    4.87461991,    -2.96892868, -1.57066343, 0.9284163 ],
        [ -1.39892467,  4.50539943,    -2.95543943, -1.45585892, 0.91347106 ],
        [ -1.27969404,  4.15403753,    -2.93741766, -1.34793582, 0.89635935 ],
        [ -1.16825275,  3.82789965,    -2.91457611, -1.24906794, 0.87749788 ],
        [ -1.05706144,  3.50542455,    -2.88442315, -1.15268638, 0.8555769 ],
        [ -0.941951,    3.1756403,     -2.84368954, -1.05567568, 0.82932379 ],
        [ -0.81898435,  2.82926981,    -2.78756742, -0.95563417, 0.79707985 ],
        [ -0.70073861,  2.50351093,    -2.7200169,  -0.86341765, 0.76202053 ],
        [ -0.59245865,  2.2128788,     -2.64629853, -0.78279148, 0.72672643 ],
        [ -0.49089406,  1.94802937,    -2.5676801,  -0.71076317, 0.69132618 ],
        [ -0.39493137,  1.70547792,    -2.48647387, -0.6460884,  0.65640075 ],
        [ -0.30135924,  1.4767054,     -2.40270078, -0.58629382, 0.6215647 ],
        [ -0.20881075,  1.25827375,    -2.31744577, -0.53037159, 0.58695474 ],
        [ -0.11587138,  1.04689132,    -2.23148332, -0.47742125, 0.55262077 ],
        [ -0.02214819,  0.84175877,    -2.14632583, -0.42721524, 0.51894399 ],
        [ 0.0733437,    0.64082382,    -2.06273906, -0.37924387, 0.48604587 ],
        [ 0.1715459,    0.44229595,    -1.98136244, -0.33310079, 0.45404125 ],
        [ 0.27322584,   0.24487477,    -1.90284586, -0.28852715, 0.42308533 ],
        [ 0.37937827,   0.04694057,    -1.82755008, -0.24522167, 0.39325198 ],
        [ 0.49090551,   -0.1528124,    -1.75583561, -0.20298452, 0.36464088 ],
        [ 0.60867491,   -0.35552523,   -1.6880144,  -0.16167645, 0.33735284 ],
        [ 0.73390612,   -0.56283402,   -1.62415527, -0.12108561, 0.31140619 ],
        [ 0.86771641,   -0.77607192,   -1.56438231, -0.0810936,  0.28685136 ],
        [ 1.01139173,   -0.99673893,   -1.50873149, -0.04157989, 0.26371008 ],
        [ 1.16098243,   -1.21860556,   -1.45884573, -0.00373962, 0.24269122 ],
        [ 1.31668041,   -1.44218774,   -1.41430031, 0.03253061,  0.22365851 ],
        [ 1.47956964,   -1.66922959,   -1.3744291,  0.0675207,   0.20636896 ],
        [ 1.6504208,    -1.9009223,    -1.33873663, 0.1014039,   0.19064665 ],
        [ 1.83151063,   -2.1403592,    -1.30653574, 0.13459044,  0.17622299 ],
        [ 2.02205323,   -2.38649356,   -1.27776328, 0.16688812,  0.1631031 ],
        [ 2.13184589,   -2.52596929,   -1.26314705, 0.18442015,  0.15635394 ],
        [ 2.19354419,   -2.6036673,    -1.25554863, 0.19395621,  0.15279124 ],
        [ 2.2425893,    -2.66510189,   -1.24976155, 0.20138273,  0.15006894 ],
        [ 2.28032492,   -2.71217943,   -1.24539622, 0.20700717,  0.14803743 ],
        [ 2.56158811,   -3.05832244,   -1.21695731, 0.24667161,  0.13443515 ],
        [ 2.76861594,   -3.3084136,    -1.19949781, 0.27360076,  0.12590013 ],
        [ 3.06017847,   -3.65503089,   -1.17884845, 0.30875873,  0.11556541 ],
        [ 3.25768371,   -3.88665845,   -1.16694616, 0.33097207,  0.10948489 ],
        [ 3.59051445,   -4.27212851,   -1.14994554, 0.36588887,  0.10058911 ],
        [ 3.78541994,   -4.49541194,   -1.14141184, 0.38504337,  0.0960349 ],
        [ 4.17604473,   -4.938337,     -1.12686884, 0.42095645,  0.08807758 ],
        [ 4.54227611,   -5.34891296,   -1.11564325, 0.45202737,  0.08176452 ],
        [ 4.81480666,   -5.65196184,   -1.10846828, 0.47373928,  0.0776452 ],
        [ 5.28198991,   -6.16730396,   -1.09804211, 0.50853891,  0.07150748 ],
        [ 5.72441416,   -6.65124103,   -1.0898559,  0.53905397,  0.0665612 ],
        [ 6.22837164,   -7.19845798,   -1.08206075, 0.5713491,   0.06173203 ],
        [ 6.77187304,   -7.78459715,   -1.07505196, 0.60365981,  0.05728265 ],
        [ 7.34893708,   -8.40312182,   -1.06882598, 0.63551856,  0.05323742 ],
        [ 8.01958126,   -9.11783742,   -1.06279232, 0.66984002,  0.04922581 ],
        [ 8.71440933,   -9.85443049,   -1.05757853, 0.70278088,  0.04568349 ],
        [ 9.48162997,   -10.66392208,  -1.05275703, 0.73651401,  0.04234056 ],
        [ 10.3302466,   -11.55536053,  -1.0482977,  0.77107145,  0.03918806 ],
        [ 11.26686122,  -12.53522893,  -1.04418791, 0.8063519,   0.0362283 ],
        [ 12.28267474,  -13.59398113,  -1.04046372, 0.84172962,  0.03349932 ],
        [ 13.3872878,   -14.74136344,  -1.03707416, 0.87730172,  0.03097512 ],
        [ 14.63370112,  -16.03193492,  -1.03387953, 0.91436171,  0.02855883 ],
        [ 15.94091386,  -17.38154989,  -1.03107644, 0.95025082,  0.02640821 ],
        [ 17.36412656,  -18.84711402,  -1.02851297, 0.98637882,  0.02441495 ],
        [ 18.91261531,  -20.43788702,  -1.02616823, 1.02271771,  0.02256889 ],
        [ 20.67381188,  -22.24314377,  -1.02393326, 1.06085376,  0.020788 ],
        [ 22.50932702,  -24.12074415,  -1.02197929, 1.09752679,  0.01921359 ],
        [ 24.595349,    -26.25061338,  -1.02011531, 1.1359792,   0.01769564 ],
        [ 26.86326667,  -28.56217035,  -1.01841889, 1.17448486,  0.01630013 ],
        [ 29.32120783,  -31.06343465,  -1.01687778, 1.21293309,  0.01502029 ],
        [ 31.97521731,  -33.76033003,  -1.01548043, 1.25120068,  0.01384949 ],
        [ 34.82855812,  -36.65598667,  -1.01421615, 1.28915257,  0.01278127 ],
        [ 37.88073024,  -39.74975807,  -1.01307464, 1.32664031,  0.01180937 ],
        [ 41.34273492,  -43.25508413,  -1.01198387, 1.36585968,  0.01087376 ],
        [ 45.04325759,  -46.99808729,  -1.01100327, 1.40448791,  0.01002663 ],
        [ 49.25186016,  -51.25096924,  -1.0100669,  1.4449232,   0.00921221 ],
        [ 53.74437843,  -55.78677627,  -1.00922892, 1.48461527,  0.00847863 ],
        [ 58.50307905,  -60.58756581,  -1.00848141, 1.52335342,  0.00782029 ],
        [ 63.9084393,   -66.03678379,  -1.00776708, 1.56386621,  0.00718756 ],
        [ 69.61264136,  -71.7834318,   -1.00713328, 1.60321118,  0.00662308 ],
        [ 76.09928302,  -78.31432236,  -1.00652774, 1.64436383,  0.00608094 ],
        [ 82.90941256,  -85.16703636,  -1.0059937,  1.68409433,  0.00560043 ],
        [ 90.6571672,   -92.95919333,  -1.00548345, 1.72564653,  0.00513916 ],
        [ 98.73516534,  -101.07963031, -1.0050365,  1.76547601,  0.00473329 ],
        [ 107.92162301, -110.31033549, -1.00460933, 1.80711718,  0.00434372 ],
        [ 117.41143478, -119.84207804, -1.00423807, 1.84668112,  0.00400377 ],
        [ 128.18691368, -130.66125534, -1.003883,   1.88801507,  0.0036774 ],
        [ 139.18283376, -141.69814533, -1.0035772,  1.92685891,  0.0033953 ],
        [ 151.63203763, -154.19000803, -1.0032844,  1.96739347,  0.00312426 ],
        [ 165.80103976, -168.40350514, -1.00300452, 2.00976592,  0.00286432 ],
        [ 180.10495002, -182.74864305, -1.00276655, 2.04910826,  0.00264259 ],
        [ 196.31380846, -199.00044477, -1.00253872, 2.09016751,  0.00242968 ],
        [ 212.33408199, -215.05981878, -1.00234764, 2.1276205,   0.00225062 ]
    ],

    energy_table => [
        [ -3.90224117,  2.818e-05,  6.863e-05,  2.818e-05,  6.863e-05,     5.0221225,    -1.511239 ],
        [ -3.62566258,  5.514e-05,  0.00013332, 5.514e-05,  0.00013332,    4.6038762,    -1.5133418 ],
        [ -3.44997359,  8.423e-05,  0.00020256, 8.423e-05,  0.00020256,    4.337873,     -1.5148979 ],
        [ -3.31302799,  0.00011699, 0.00028002, 0.00011699, 0.00028002,    4.1303247,    -1.5167932 ],
        [ -2.88942605,  0.00031964, 0.00075103, 0.00031964, 0.00075103,    3.4861831,    -1.5260122 ],
        [ -2.62494077,  0.00059231, 0.00137025, 0.00059231, 0.00137025,    3.0816864,    -1.5369661 ],
        [ -2.35396093,  0.00110192, 0.00249741, 0.00110192, 0.00249741,    2.6630915,    -1.5522132 ],
        [ -2.14727559,  0.00175244, 0.0038932,  0.00175244, 0.0038932,     2.3410946,    -1.5676224 ],
        [ -1.96296451,  0.00262777, 0.00571012, 0.00262777, 0.00571012,    2.0505679,    -1.5863365 ],
        [ -1.80104846,  0.00372087, 0.00789518, 0.00372087, 0.00789518,    1.7922847,    -1.6058707 ],
        [ -1.65906399,  0.00501071, 0.01036742, 0.00501071, 0.01036742,    1.5629434,    -1.6262179 ],
        [ -1.523553,    0.00660665, 0.01328146, 0.00660665, 0.01328146,    1.3411553,    -1.6482942 ],
        [ -1.39892467,  0.00845552, 0.01647044, 0.00845552, 0.01647044,    1.1343999,    -1.6706443 ],
        [ -1.27969404,  0.01062364, 0.01996766, 0.01062364, 0.01996766,    0.93387697,   -1.6934411 ],
        [ -1.16825275,  0.01304741, 0.02357921, 0.01304741, 0.02357921,    0.74394591,   -1.7154181 ],
        [ -1.05706144,  0.01588086, 0.0274137,  0.01588086, 0.0274137,     0.55197451,   -1.7374436 ],
        [ -0.941951,    0.01926975, 0.03146342, 0.01926975, 0.03146342,    0.35067207,   -1.7595352 ],
        [ -0.81898435,  0.02339775, 0.03562253, 0.02339775, 0.03562253,    0.13289469,   -1.7816605 ],
        [ -0.70073861,  0.02782501, 0.0391579,  0.02782501, 0.0391579,     -0.078988049, -1.801228 ],
        [ -0.59245865,  0.03221147, 0.04173983, 0.03221147, 0.04173983,    -0.27495137,  -1.8177154 ],
        [ -0.49089406,  0.03654263, 0.04341875, 0.03654263, 0.04341875,    -0.46032249,  -1.8323067 ],
        [ -0.39493137,  0.04075538, 0.04425595, 0.04075538, 0.04425595,    -0.63680401,  -1.8457771 ],
        [ -0.30135924,  0.04490649, 0.04435255, 0.04490649, 0.04435255,    -0.81012983,  -1.8589629 ],
        [ -0.20881075,  0.04898987, 0.04378641, 0.04898987, 0.04378641,    -0.98278181,  -1.8722173 ],
        [ -0.11587138,  0.05301004, 0.04263695, 0.05301004, 0.04263695,    -1.1574086,   -1.8856985 ],
        [ -0.02214819,  0.05693268, 0.04100026, 0.05693268, 0.04100026,    -1.3347818,   -1.8992506 ],
        [ 0.0733437,    0.06075329, 0.03896928, 0.06075329, 0.03896928,    -1.5167995,   -1.9126859 ],
        [ 0.1715459,    0.06446702, 0.03663415, 0.06446702, 0.03663415,    -1.7052939,   -1.9257605 ],
        [ 0.27322584,   0.068063,   0.03408495, 0.068063,   0.03408495,    -1.9017693,   -1.9381946 ],
        [ 0.37937827,   0.07153848, 0.03140081, 0.07153848, 0.03140081,    -2.1081674,   -1.9497454 ],
        [ 0.49090551,   0.07488641, 0.02865674, 0.07488641, 0.02865674,    -2.3262487,   -1.9602058 ],
        [ 0.60867491,   0.07809829, 0.02592116, 0.07809829, 0.02592116,    -2.5576976,   -1.9694361 ],
        [ 0.73390612,   0.08117427, 0.02324639, 0.08117427, 0.02324639,    -2.804887,    -1.9773654 ],
        [ 0.86771641,   0.08410974, 0.02067953, 0.08410974, 0.02067953,    -3.0699794,   -1.984023 ],
        [ 1.01139173,   0.08690278, 0.01825633, 0.08690278, 0.01825633,    -3.3554844,   -1.9893979 ],
        [ 1.16098243,   0.08946648, 0.01607504, 0.08946648, 0.01607504,    -3.6534293,   -1.9934133 ],
        [ 1.31668041,   0.09181419, 0.01413463, 0.09181419, 0.01413463,    -3.9640719,   -1.9963621 ],
        [ 1.47956964,   0.09397257, 0.01241527, 0.09397257, 0.01241527,    -4.2894625,   -1.9984487 ],
        [ 1.6504208,    0.09596041, 0.01089919, 0.09596041, 0.01089919,    -4.6310485,   -1.9998455 ],
        [ 1.83151063,   0.09780892, 0.00955728, 0.09780892, 0.00955728,    -4.9933011,   -2.0007126 ],
        [ 2.02205323,   0.09951477, 0.00838454, 0.09951477, 0.00838454,    -5.3745842,   -2.0012097 ],
        [ 2.13184589,   0.10040274, 0.00780119, 0.10040274, 0.00780119,    -5.5943132,   -2.0013813 ],
        [ 2.19354419,   0.10087466, 0.00749913, 0.10087466, 0.00749913,    -5.7177975,   -2.0014528 ],
        [ 2.2425893,    0.10123682, 0.00727115, 0.10123682, 0.00727115,    -5.8159602,   -2.0014959 ],
        [ 2.28032492,   0.101508,   0.00710267, 0.101508,   0.00710267,    -5.8914884,   -2.0014971 ],
        [ 2.56158811,   0.10334602, 0.00601259, 0.10334602, 0.00601259,    -6.454408,    -2.0013708 ],
        [ 2.76861594,   0.10452175, 0.00536453, 0.10452175, 0.00536453,    -6.8687433,   -2.0012978 ],
        [ 3.06017847,   0.10597308, 0.00461989, 0.10597308, 0.00461989,    -7.4522251,   -2.0010771 ],
        [ 3.25768371,   0.10684339, 0.0042034,  0.10684339, 0.0042034,     -7.8474287,   -2.0009169 ],
        [ 3.59051445,   0.10814226, 0.00362445, 0.10814226, 0.00362445,    -8.5133611,   -2.0007189 ],
        [ 3.78541994,   0.10882061, 0.0033425,  0.10882061, 0.0033425,     -8.9033013,   -2.0006115 ],
        [ 4.17604473,   0.11003106, 0.00287419, 0.11006391, 0.0030087716,  -9.6847493,   -2.0004257 ],
        [ 4.54227611,   0.11101751, 0.00252525, 0.11118261, 0.0032345525,  -10.41734,    -2.0002969 ],
        [ 4.81480666,   0.11167547, 0.00230861, 0.11214286, 0.0037689468,  -10.962471,   -2.0002286 ],
        [ 5.28198991,   0.11267962, 0.00200227, 0.11387402, 0.0033292823,  -11.896921,   -2.0001456 ],
        [ 5.72441416,   0.11351232, 0.00176994, 0.115226,   0.0028409246,  -12.78182,    -2.0000935 ],
        [ 6.22837164,   0.11434839, 0.00155582, 0.11657044, 0.0025269806,  -13.78977,    -2.0000555 ],
        [ 6.77187304,   0.11514158, 0.00136958, 0.11785487, 0.0022264072,  -14.876795,   -2.0000326 ],
        [ 7.34893708,   0.11588423, 0.00121071, 0.11907168, 0.0019965099,  -16.030937,   -2.0000181 ],
        [ 8.01958126,   0.11664352, 0.00105463, 0.12047956, 0.0020951658,  -17.372233,   -2.0000092 ],
        [ 8.71440933,   0.11733535, 0.00093548, 0.12184442, 0.001842613,   -18.761893,   -2.0000045 ],
        [ 9.48162997,   0.1180088,  0.00082406, 0.12316897, 0.0016184372,  -20.296337,   -2.0000022 ],
        [ 10.3302466,   0.11866444, 0.00072472, 0.12445489, 0.0014195303,  -21.993571,   -2.0000009 ],
        [ 11.26686122,  0.11930041, 0.00063652, 0.12569907, 0.0012437607,  -23.866801,   -2.0000003 ],
        [ 12.28267474,  0.11990652, 0.00055959, 0.12688209, 0.0010910663,  -25.898428,   -2.0000001 ],
        [ 13.3872878,   0.12048609, 0.00049221, 0.12801101, 0.00095781623, -28.107654,   -2.0000001 ],
        [ 14.63370112,  0.12106012, 0.00043114, 0.12912702, 0.0008374726,  -30.600481,   -2 ],
        [ 15.94091386,  0.12158884, 0.00037964, 0.13015326, 0.0007362909,  -33.214907,   -2.0000001 ],
        [ 17.36412656,  0.12209577, 0.00033434, 0.13113573, 0.00064754155, -36.061332,   -2.0000001 ],
        [ 18.91261531,  0.12258156, 0.00029451, 0.13207601, 0.00056968622, -39.15831,    -2 ],
        [ 20.67381188,  0.12306696, 0.00025805, 0.13301438, 0.00049857228, -42.680703,   -2 ],
        [ 22.50932702,  0.12351152, 0.00022744, 0.1338729,  0.000438998,   -46.351733,   -2 ],
        [ 24.595349,    0.12395565, 0.00019941, 0.13472975, 0.00038453245, -50.523777,   -2 ],
        [ 26.86326667,  0.12437912, 0.00017494, 0.13554603, 0.00033706791, -55.059612,   -2 ],
        [ 29.32120783,  0.12478195, 0.00015362, 0.13632192, 0.00029577011, -59.975495,   -2 ],
        [ 31.97521731,  0.12516415, 0.00013507, 0.13705758, 0.00025989437, -65.283514,   -2 ],
        [ 34.82855812,  0.12552577, 0.00011897, 0.13775318, 0.00022878014, -70.990195,   -2 ],
        [ 37.88073024,  0.12586683, 0.00010501, 0.1384089,  0.00020183651, -77.09454,    -2 ],
        [ 41.34273492,  0.12620742, 9.221e-05,  0.13906339, 0.00017715476, -84.018549,   -2 ],
        [ 45.04325759,  0.12652751, 8.118e-05,  0.13967821, 0.00015589367, -91.419594,   -2 ],
        [ 49.25186016,  0.12684712, 7.108e-05,  0.14029188, 0.00013644957, -99.836799,   -2 ],
        [ 53.74437843,  0.1271463,  6.242e-05,  0.1408661,  0.00011979184, -108.82184,   -2 ],
        [ 58.50307905,  0.12742512, 5.502e-05,  0.14140108, 0.00010555174, -118.33924,   -2 ],
        [ 63.9084393,   0.12770351, 4.824e-05,  0.1419351,  9.2514107e-05, -129.14996,   -2 ],
        [ 69.61264136,  0.12796162, 4.247e-05,  0.1424301,  8.1433934e-05, -140.55836,   -2 ],
        [ 76.09928302,  0.12821935, 3.719e-05,  0.14292423, 7.129547e-05,  -153.53165,   -2 ],
        [ 82.90941256,  0.12845688, 3.273e-05,  0.14337957, 6.2733554e-05, -167.1519,    -2 ],
        [ 90.6571672,   0.12869406, 2.865e-05,  0.14383415, 5.4900348e-05, -182.64741,   -2 ],
        [ 98.73516534,  0.12891115, 2.522e-05,  0.14425015, 4.8330255e-05, -198.80341,   -2 ],
        [ 107.92162301, 0.12912792, 2.209e-05,  0.1446655,  4.2316682e-05, -217.17633,   -2 ],
        [ 117.41143478, 0.12932471, 1.948e-05,  0.14504249, 3.7310547e-05, -236.15595,   -2 ],
        [ 128.18691368, 0.12952121, 1.708e-05,  0.14541891, 3.2722985e-05, -257.70691,   -2 ],
        [ 139.18283376, 0.12969783, 1.511e-05,  0.1457572,  2.8935816e-05, -279.69875,   -2 ],
        [ 151.63203763, 0.12987421, 1.329e-05,  0.14609501, 2.545787e-05,  -304.59715,   -2 ],
        [ 165.80103976, 0.13005035, 1.163e-05,  0.14643233, 2.2274974e-05, -332.93516,   -2 ],
        [ 180.10495002, 0.13020672, 1.028e-05,  0.14673177, 1.9682292e-05, -361.54298,   -2 ],
        [ 196.31380846, 0.13036289, 9.04e-06,   0.14703082, 1.7302224e-05, -393.9607,    -2 ],
        [ 212.33408199, 0.13049938, 8.04e-06,   0.14729216, 1.5386417e-05, -426.00124,   -2 ]
    ],

    impulse_table => [
        [
            -4.5071309, 13.814838, -4.5076578,    0.046000237, -0.0077082523, -0.00053243264,
            -2.685246,  -0.348137, -0.0097996095, -0.14925345, 0.49383194,    0.9999567,
            -4.3530681
        ],
        [
            -4.1997043, 12.892561,   -4.2005401,   0.053501392, -0.010480015, 0.003254705,
            -2.6812837, -0.50187879, -0.013325557, -0.14544044, 0.49383126,   0.99990942,
            -4.1355469
        ],
        [
            -3.9120222, 12.029521,  -3.9133093,   0.061597036, -0.013966561, 0.0079003479,
            -2.6762988, -0.6457792, -0.017764398, -0.14075473, 0.49382995,   0.99981775,
            -3.9358683
        ],
        [
            -3.5065571, 10.813145,   -3.5089228,   0.075052392, -0.020915813, 0.01674991,
            -2.6663518, -0.84873385, -0.026630932, -0.13185899, 0.49382463,   0.99950998,
            -3.662173
        ],
        [
            -2.9957315, 9.2807691, -3.0008273,   0.095986011, -0.034652912, 0.032855676,
            -2.6465722, -1.105146, -0.044278223, -0.11641913, 0.4937974,    0.99830372,
            -3.334657
        ],
        [
            -2.3025843, 7.2022285,  -2.3170615,   0.13239181,   -0.067446031, 0.06751514,
            -2.5980343, -1.4588599, -0.087246241, -0.092074423, 0.4935548,    0.99119188,
            -2.9415035
        ],
        [
            -1.8971192, 5.9882666,  -1.9238451,  0.1574894,   -0.097367252, 0.097811098,
            -2.5512108, -1.6774013, -0.12722372, -0.08337717, 0.49289919,   0.97790928,
            -2.7587489
        ],
        [
            -1.6094371,  5.1299112,  -1.650785,  0.17607951, -0.12406985, 0.12501134, -2.5065049, -1.8472037,
            -0.16297348, -0.0829564, 0.49163435, 0.95895347, -2.6617162
        ],
        [
            -1.3862935, 4.4680792,  -1.4443313,  0.19027202,   -0.1474683, 0.14946354,
            -2.4642194, -1.9954726, -0.19389115, -0.086638109, 0.48958236, 0.93532312,
            -2.6108423
        ],
        [
            -1.203972, 3.9321516,  -1.2805261,  0.20132274,   -0.16763809, 0.17132128,
            -2.424569, -2.1339913, -0.21991579, -0.092466952, 0.48659503,  0.90820151,
            -2.5885646
        ],
        [
            -1.0498213, 3.4845491,  -1.1464975,  0.21008304,   -0.18477133, 0.1906996,
            -2.3876812, -2.2685164, -0.24136042, -0.099885132, 0.48256497,  0.87877513,
            -2.5855538
        ],
        [
            -0.91628992, 3.1028024,   -1.0344756, 0.21715912, -0.19914106, 0.20772655, -2.3535968, -2.4018785,
            -0.25874346, -0.10849691, 0.47743396, 0.848123,   -2.5961958
        ],
        [
            -0.79850688, 2.7722965,   -0.93937113, 0.22298621, -0.21106759, 0.22256056, -2.3222763, -2.5352929,
            -0.27266159, -0.11807516, 0.47119665,  0.81715529, -2.6167808
        ],
        [
            -0.75080966, 2.6404732,   -0.90195917, 0.22523626, -0.21562505, 0.22843648, -2.3093548, -2.5940612,
            -0.27783855, -0.12235665, 0.46811414,  0.80363123, -2.6282752
        ],
        [
            -0.59783619, 2.2271199,   -0.78670437, 0.23203886, -0.22893334, 0.24642833, -2.2674471, -2.8028137,
            -0.29242764, -0.13909829, 0.45562945,  0.75695216, -2.677952
        ],
        [
            -0.51082481, 1.9993661,   -0.72461269, 0.23563292, -0.23550719, 0.25589227, -2.2435885, -2.9360929,
            -0.29929313, -0.14994616, 0.44651293,  0.72860182, -2.7150994
        ],
        [
            -0.4307821,  1.7951782,   -0.66985757, 0.23876242, -0.24087906, 0.26399782, -2.221831, -3.0682148,
            -0.30469962, -0.16107727, 0.43669546,  0.70175552, -2.7549701
        ],
        [
            -0.35667413, 1.6109975,  -0.62124745, 0.24150418, -0.24527831, 0.27094557, -2.2019691, -3.1985439,
            -0.30896758, -0.1722726, 0.4263354,   0.67652135, -2.7966638
        ],
        [
            -0.28768126, 1.4439269,   -0.57782712, 0.24391587, -0.24889519, 0.27691556, -2.1838062, -3.3265277,
            -0.31235058, -0.18308721, 0.41559376,  0.65292713, -2.8394799
        ],
        [
            -0.22314274, 1.2915826,   -0.53882215, 0.24604275, -0.25188437, 0.28206455, -2.1671608, -3.4517253,
            -0.31504672, -0.19405888, 0.40462667,  0.63094505, -2.882876
        ],
        [
            -0.16251812, 1.1519866,  -0.50359842, 0.24792187, -0.25436987, 0.28652609, -2.1518707, -3.573813,
            -0.31720867, -0.2048431, 0.39357934,  0.61051104, -2.9264353
        ],
        [
            -0.1053597,  1.0234855,   -0.47163245, 0.24958445, -0.25645027, 0.29041213, -2.1377887, -3.6925743,
            -0.31895381, -0.21508923, 0.3825818,   0.59153938, -2.9698413
        ],
        [
            -0.05129248, 0.90469352,  -0.44248998, 0.25105735, -0.25820355, 0.29381556, -2.1247882, -3.8078848,
            -0.32037208, -0.22474426, 0.37174602,  0.57393324, -3.0128565
        ],
        [
            8.1448994e-07, 0.7944387,  -0.4158074,  0.25236394,  -0.25969142, 0.29681313,
            -2.1127576,    -3.9196935, -0.32153246, -0.23410557, 0.36116438,  0.55759187,
            -3.0553056
        ],
        [
            0.095310994, 0.59571632,  -0.36864747, 0.25455724, -0.26205608, 0.30183214, -2.0912222, -4.132874,
            -0.32327999, -0.25239762, 0.34103401,  0.52830765, -3.1380364
        ],
        [
            0.18232237,  0.42099032,  -0.32822611, 0.25629881, -0.26382673, 0.3058521, -2.0725362, -4.3326074,
            -0.32449353, -0.26954781, 0.32255146,  0.50294297, -3.2174252
        ],
        [
            0.26236508,  0.26558502,  -0.29313952, 0.25769127, -0.26518456, 0.30913477, -2.0561987, -4.5197372,
            -0.32535457, -0.28505972, 0.30583873,  0.48085238, -3.2932351
        ],
        [
            0.40546592, -0.00050747096, -0.23505351, 0.25972267,  -0.26709793, 0.3141615,
            -2.0290866, -4.8601262,     -0.32643461, -0.31375371, 0.27748708,  0.44442753,
            -3.434157
        ],
        [
            0.53062907,  -0.22208564, -0.18869007, 0.26107804, -0.2683574, 0.31782536, -2.0076179, -5.161764,
            -0.32703151, -0.3393879,  0.25491168,  0.41577132, -3.5617689
        ],
        [
            0.693148,    -0.49623402, -0.13394217, 0.26236046, -0.26958129, 0.32177539, -1.982871, -5.5559257,
            -0.32749048, -0.37163197, 0.22909986,  0.38278415, -3.7316433
        ],
        [
            0.99325259, -0.96931242, -0.046388227, 0.26371269,  -0.27106008, 0.32729498,
            -1.9461113, -6.2833592,  -0.32781496,  -0.43213073, 0.19125825,  0.33288175,
            -4.052139
        ],
        [
            1.0986131,  -1.1270052, -0.019135765, 0.26397299,  -0.27142492, 0.3288209,
            -1.9357592, -6.5371609, -0.32784329,  -0.45266997, 0.18056631,  0.31825186,
            -4.1655836
        ],
        [
            1.2089612,   -1.2882496,  0.007755862, 0.26416436, -0.27174931, 0.33023892, -1.926203, -6.8016665,
            -0.32784846, -0.47493646, 0.17056889,  0.30429636, -4.2845145
        ],
        [
            1.3862952,   -1.5400174,  0.047830289, 0.26434341, -0.27217635, 0.3321915, -1.9133441, -7.2235785,
            -0.32782491, -0.51075596, 0.1567101,   0.28444326, -4.4754752
        ],
        [
            1.6094387,   -1.8458944,  0.09351818, 0.26442022, -0.27259468, 0.33418291, -1.9009324, -7.7485877,
            -0.32777065, -0.55380331, 0.14242301, 0.26326189, -4.7148591
        ],
        [
            1.7917603,   -2.0882932,  0.12752649, 0.26441063, 0, 0.33550254, 0, -8.1726482,
            -0.35551907, -0.72337845, 0.13280644, 0.24853495, -4.9093684
        ],
        [
            1.945911,    -2.2887905,  0.15428015, 0.26437411, 0, 0.33644293, 0, -8.5278628,
            -0.14055643, 0.074527136, 0.12582649, 0.23758073, -5.0729629
        ],
        [
            2.3027887,   -2.7401277,  0.21031937, 0.2642465, -0.27041852, 0.33813912, -1.6651364, -9.3394155,
            -0.41386264, -0.82437623, 0.11277095, 0.2148842, -5.4546075
        ],
        [
            2.9958854,   -3.5791053,  0.30126079,  0.26400171, -0.26390222, 0.34011845, -1.4456715, -10.879478,
            -0.40912255, -0.99389241, 0.095557429, 0.18641992, -6.169971
        ],
        [
            3.9122026,   -4.6397996, 0.39704373,  0.26380772, -0.26700416, 0.34131708, -1.4765902, -12.863361,
            -0.40979965, -1.2131985, 0.081673832, 0.16164878, -7.1028145
        ],
        [
            4.6055618,   -5.4194617, 0.45717025,  0.26373641, -0.27417569, 0.34171992, -1.4765684, -14.338622,
            -0.40978037, -1.3780982, 0.074633787, 0.14845624, -7.8026531
        ],
        [
            5.2984002,   -6.1853205, 0.50971076,  0.26370004, -0.27906222, 0.34192166, -1.5210118, -15.797975,
            -0.40977131, -1.5424048, 0.069283596, 0.13817367, -8.4990834
        ],
        [
            6.2147684,   -7.1837371, 0.57050852,  0.26367825, -0.29795467, 0.34204305, -1.7445166, -17.712387,
            -0.40976614, -1.7593639, 0.063825282, 0.12749737, -9.4178234
        ],
        [
            6.9078913,   -7.9307166, 0.61138226,  0.26367111, -0.30006277, 0.34208356, -1.765886, -19.151812,
            -0.40976448, -1.9233348, 0.060512736, 0.12094997, -10.111808
        ],
        [
            7.6009565,   -8.672181,  0.6487338,   0.26366761, -0.27834107, 0.34210383, -1.6920319, -20.585571,
            -0.40976366, -2.0872368, 0.057695041, 0.11535274, -10.805331
        ],
        [
            8.5172026,   -9.6457334, 0.69367892,  0.26366555, -0.29482897, 0.34211599, -2.0045178, -22.474375,
            -0.40976317, -2.3038768, 0.054534222, 0.10905368, -11.721868
        ],
        [
            9.2104235,   -10.378192, 0.72488042, 0.26366488, -0.28417735, 0.34212005, -1.9396913, -23.899373,
            -0.40976301, -2.4677687, 0.05247112, 0.10493491, -12.415193
        ],
        [
            10.819855,   -12.068059, 0.78985941, 0.26366435,  -0.26895948, 0.34212329, -1.4765648, -27.197365,
            -0.40976289, -2.848252,  0.04847874, 0.096956031, -14.024711
        ],
        [
            11.513065,   -12.792192, 0.81518448,  0.26366429,  -0.26895958, 0.3421237, -1.4765649, -28.614255,
            -0.40976288, -3.0121295, 0.047023841, 0.094046961, -14.717933
        ],
        [
            13.122485,   -14.466644, 0.86902474,  0.26366423,  -0.26895965, 0.34212402, -1.4765649, -31.897268,
            -0.40976287, -3.3926004, 0.044100945, 0.088201748, -16.327363
        ],
        [
            13.815693,   -15.185399, 0.89038225,  0.26366422,  -0.26895964, 0.34212407, -1.4765649, -33.308929,
            -0.40976287, -3.556476,  0.043001535, 0.086002999, -17.020572
        ],
        [
            16.118116,   -17.564228, 0.95490679,  0.26366419,  -0.26895964, 0.34212417, -1.4765648, -37.989384,
            -0.40976288, -4.1007727, 0.039870094, 0.079740181, -19.322996
        ],
        [
            18.420735,  -19.932965, 1.0114814,   0.26366415,  -0.26895964, 0.34212503, -1.4765643, -42.660166,
            -0.4097629, -4.6451151, 0.037340697, 0.074681387, -21.625616
        ]
    ],

    tail_shock_table => [
        [ 0.35687692, -0.83858,    -0.050591015, -0.050591015 ],
        [ 2.9259266,  -0.82146116, -0.051058427, -0.050123602 ],
        [ 3.8820082,  -0.98587469, -0.050239748, -0.0082887387 ],
        [ 4.589646,   -1.0630202,  -0.050101507, -0.0057116485 ],
        [ 5.290042,   -1.1352031,  -0.049333336, -0.0039357775 ],
        [ 6.2112243,  -1.2249909,  -0.048010115, -0.0023423081 ],
        [ 6.9060468,  -1.2894248,  -0.046902873, -0.0015058032 ],
        [ 7.5999995,  -1.3513121,  -0.045768528, -0.00088113699 ],
        [ 8.5168023,  -1.4297544,  -0.044285068, -0.00027996641 ]
    ],

    blast_info => [
        0.011405989,  2.6962759,   0.1601353,   -0.88850875,  3.4104636,  0.0010126545,
        0.18331545,   -0.69898776, 0.066205506, 1.0619088,    0.58206332, -0.59934466,
        0.0085707929, 0.49383314,  0.040563723, -0.041378405, 0.59028,    -0.83858,
        12.051547,    -0.86942239
    ],

};
1;
