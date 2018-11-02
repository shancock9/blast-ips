package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='C', gamma=1.7
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'C1.7'} = {

    shock_table_info => [ 1, 1.7, 6.79e-07, 8000, 2.2e-07, 70, 1.27e-07, 5.2e-08, 5e-07, 90.735, -1.3086 ],

    shock_table => [
        [ -6.52959346, 12.00017706,  -1.99998453, -6.5309851,  0.99860741 ],
        [ -5.91977686, 10.78056107,  -1.99994762, -5.92233907, 0.99743457 ],
        [ -5.62763373, 10.19629405,  -1.99991227, -5.63106674, 0.99656124 ],
        [ -5.00941899, 8.95997168,   -1.99969842, -5.01579836, 0.99360124 ],
        [ -4.54662574, 8.03461498,   -1.99923861, -4.5567773,  0.98980076 ],
        [ -4.16008756, 7.26198217,   -1.99835373, -4.17506196, 0.9849258 ],
        [ -3.84095842, 6.62445836,   -1.99688991, -3.86161268, 0.97916462 ],
        [ -3.5618044,  6.06731099,   -1.99458364, -3.58918468, 0.97231948 ],
        [ -3.31478724, 5.57500244,   -1.99116949, -3.34994382, 0.96438237 ],
        [ -3.09229264, 5.13247863,   -1.98632289, -3.13634784, 0.95528082 ],
        [ -2.88849789, 4.72831242,   -1.97965251, -2.9426869,  0.9449077 ],
        [ -2.69903596, 4.35404289,   -1.97069223, -2.76474298, 0.93312976 ],
        [ -2.51998771, 4.00219738,   -1.95885114, -2.59883062, 0.91974431 ],
        [ -2.34776995, 3.66611687,   -1.94335954, -2.44171408, 0.90446046 ],
        [ -2.17806875, 3.33796925,   -1.92307379, -2.28968792, 0.88677217 ],
        [ -2.00484321, 3.00708296,   -1.89606431, -2.13785286, 0.86573925 ],
        [ -1.81246992, 2.64590116,   -1.85731559, -1.97385228, 0.83861143 ],
        [ -1.63442945, 2.31908655,   -1.81245256, -1.82704945, 0.80991719 ],
        [ -1.47092468, 2.02663884,   -1.76360189, -1.69697778, 0.78069398 ],
        [ -1.31720571, 1.75946638,   -1.71161677, -1.57922505, 0.75104167 ],
        [ -1.17127608, 1.51356824,   -1.65786393, -1.47177606, 0.7213642 ],
        [ -1.02940358, 1.28224765,   -1.60274125, -1.37154021, 0.69156318 ],
        [ -0.88784906, 1.05935991,   -1.54626909, -1.27577888, 0.66140488 ],
        [ -0.74653643, 0.8448525,    -1.48973433, -1.18443922, 0.63138154 ],
        [ -0.60390301, 0.6363818,    -1.43369153, -1.09651356, 0.60163909 ],
        [ -0.45864522, 0.43215371,   -1.37866123, -1.01126226, 0.57235423 ],
        [ -0.30953342, 0.23061417,   -1.32507949, -0.92807322, 0.54370266 ],
        [ -0.15528663, 0.0302698,    -1.27327586, -0.84638221, 0.51584004 ],
        [ 0.00537479,  -0.17024117,  -1.2235175,  -0.76569875, 0.48891724 ],
        [ 0.17367343,  -0.37209627,  -1.17604071, -0.68562222, 0.46308913 ],
        [ 0.35104205,  -0.57662,     -1.13098852, -0.60570666, 0.43847293 ],
        [ 0.53902637,  -0.78515045,  -1.0884647,  -0.52551435, 0.4151707 ],
        [ 0.73928093,  -0.99903598,  -1.04854811, -0.44461524, 0.39327207 ],
        [ 0.95382529,  -1.21990107,  -1.01125822, -0.3624863,  0.37283178 ],
        [ 1.18416128,  -1.4487502,   -0.97670937, -0.27884088, 0.35395188 ],
        [ 1.42257429,  -1.67786271,  -0.94606964, -0.19649332, 0.33729808 ],
        [ 1.67091225,  -1.90933567,  -0.9188371,  -0.11460481, 0.32260981 ],
        [ 1.9311685,   -2.14522728,  -0.89459653, -0.03237746, 0.30966765 ],
        [ 2.20485368,  -2.38703329,  -0.8730551,  0.05077233,  0.29831409 ],
        [ 2.49372664,  -2.63639437,  -0.85394973, 0.13546946,  0.2884048 ],
        [ 2.79930993,  -2.89468979,  -0.83707456, 0.22224508,  0.27982354 ],
        [ 3.12349589,  -3.16357526,  -0.82223372, 0.31172256,  0.2724578 ],
        [ 3.46782258,  -3.4443865,   -0.80927312, 0.40442045,  0.26621486 ],
        [ 3.83420523,  -3.73876012,  -0.79803937, 0.50096175,  0.26100087 ],
        [ 4.2245456,   -4.04831516,  -0.78839447, 0.60196892,  0.25672805 ],
        [ 4.64057971,  -4.37454535,  -0.7802116,  0.70803017,  0.25331203 ],
        [ 5.08525978,  -4.71990028,  -0.77334856, 0.82005014,  0.2506608 ],
        [ 5.56147782,  -5.08677213,  -0.76767791, 0.93891898,  0.2486882 ],
        [ 6.07192529,  -5.47740103,  -0.76307926, 1.06548186,  0.24730959 ],
        [ 6.62128396,  -5.89554735,  -0.75941984, 1.20107952,  0.24643649 ],
        [ 7.21562207,  -6.34600458,  -0.756572,   1.34739106,  0.24598413 ],
        [ 7.86079443,  -6.83338569,  -0.75441989, 1.50604042,  0.24587209 ],
        [ 8.56804717,  -7.36635349,  -0.75284119, 1.67997511,  0.24602338 ],
        [ 9.35141577,  -7.95563453,  -0.75172728, 1.87282916,  0.24636833 ],
        [ 10.2321278,  -8.61732949,  -0.75097768, 2.09001537,  0.24684476 ],
        [ 11.24260414, -9.37590925,  -0.75050316, 2.33972847,  0.2473985 ],
        [ 12.43426056, -10.27006531, -0.7502268,  2.63490018,  0.24798267 ],
        [ 13.89850945, -11.36846256, -0.75008378, 2.99844912,  0.24855778 ],
        [ 15.81856214, -12.80859127, -0.7500222,  3.47624006,  0.24908984 ],
        [ 18.65223427, -14.93387268, -0.75000301, 4.18279918,  0.24954685 ],
        [ 24.25257445, -19.13413193, -0.75000005, 5.5815154,   0.24988769 ]
    ],

    energy_table => [
        [ -6.52959346, 0.00699301, 0.00574785, 0.00699301, 0.00574785,    6.6740158,    -1.1812889 ],
        [ -5.91977686, 0.01153862, 0.00946502, 0.01153862, 0.00946502,    5.9473741,    -1.2035715 ],
        [ -5.62763373, 0.01466057, 0.01200668, 0.01466057, 0.01200668,    5.5940799,    -1.2164491 ],
        [ -5.00941899, 0.02429245, 0.01977928, 0.02429245, 0.01977928,    4.8327255,    -1.250701 ],
        [ -4.54662574, 0.03536058, 0.02855996, 0.03536058, 0.02855996,    4.2472702,    -1.2825821 ],
        [ -4.16008756, 0.04824167, 0.03854088, 0.04824167, 0.03854088,    3.7458421,    -1.3163236 ],
        [ -3.84095842, 0.06215165, 0.04899123, 0.06215165, 0.04899123,    3.3207336,    -1.3503009 ],
        [ -3.5618044,  0.07731696, 0.05994959, 0.07731696, 0.05994959,    2.9393447,    -1.3855104 ],
        [ -3.31478724, 0.09347235, 0.07106969, 0.09347235, 0.07106969,    2.5928858,    -1.4220134 ],
        [ -3.09229264, 0.11049692, 0.08210777, 0.11049692, 0.08210777,    2.2725989,    -1.4596859 ],
        [ -2.88849789, 0.12831636, 0.09284288, 0.12831636, 0.09284288,    1.9713595,    -1.4977808 ],
        [ -2.69903596, 0.14687319, 0.10305283, 0.14687319, 0.10305283,    1.6841295,    -1.536248 ],
        [ -2.51998771, 0.16617909, 0.11253523, 0.16617909, 0.11253523,    1.405646,     -1.572583 ],
        [ -2.34776995, 0.18630701, 0.12108216, 0.18630701, 0.12108216,    1.1319654,    -1.6049565 ],
        [ -2.17806875, 0.20750089, 0.12848745, 0.20750089, 0.12848745,    0.85695989,   -1.630134 ],
        [ -2.00484321, 0.2303046,  0.13449056, 0.2303046,  0.13449056,    0.5728805,    -1.6402979 ],
        [ -1.81246992, 0.25662827, 0.13870278, 0.25662827, 0.13870278,    0.2572544,    -1.6243025 ],
        [ -1.63442945, 0.28146577, 0.1398357,  0.28146577, 0.1398357,     -0.029235915, -1.5762131 ],
        [ -1.47092468, 0.30424014, 0.1383318,  0.30424014, 0.1383318,     -0.28201156,  -1.5073223 ],
        [ -1.31720571, 0.32525346, 0.13472578, 0.32525346, 0.13472578,    -0.50812836,  -1.4335563 ],
        [ -1.17127608, 0.34455377, 0.12952272, 0.34455377, 0.12952272,    -0.71214286,  -1.3694102 ],
        [ -1.02940358, 0.36248692, 0.12308434, 0.36248692, 0.12308434,    -0.90247657,  -1.3240169 ],
        [ -0.88784906, 0.37939335, 0.11564488, 0.37939335, 0.11564488,    -1.0874162,   -1.2988096 ],
        [ -0.74653643, 0.3951703,  0.10756816, 0.3951703,  0.10756816,    -1.2698708,   -1.2911552 ],
        [ -0.60390301, 0.40991164, 0.0991107,  0.40991164, 0.0991107,     -1.4540346,   -1.2963845 ],
        [ -0.45864522, 0.42368131, 0.09050349, 0.42368131, 0.09050349,    -1.6431162,   -1.3098892 ],
        [ -0.30953342, 0.43653352, 0.08194592, 0.43653352, 0.08194592,    -1.8396902,   -1.3279958 ],
        [ -0.15528663, 0.44852218, 0.07360136, 0.44852218, 0.07360136,    -2.0460764,   -1.3481458 ],
        [ 0.00537479,  0.45969446, 0.0656031,  0.45969446, 0.0656031,     -2.2643649,   -1.3685518 ],
        [ 0.17367343,  0.47008833, 0.05805863, 0.47008833, 0.05805863,    -2.4964296,   -1.3880807 ],
        [ 0.35104205,  0.47974972, 0.05103948, 0.47974972, 0.05103948,    -2.7443502,   -1.4060859 ],
        [ 0.53902637,  0.4887229,  0.0445908,  0.4887229,  0.0445908,     -3.010329,    -1.4222445 ],
        [ 0.73928093,  0.49704972, 0.03873527, 0.49704972, 0.03873527,    -3.296707,    -1.4364428 ],
        [ 0.95382529,  0.50477821, 0.03347064, 0.50477821, 0.03347064,    -3.6063526,   -1.4487098 ],
        [ 1.18416128,  0.51193152, 0.02879415, 0.51193152, 0.02879415,    -3.9413876,   -1.4591252 ],
        [ 1.42257429,  0.51830737, 0.02482378, 0.51830737, 0.02482378,    -4.2903913,   -1.4675424 ],
        [ 1.67091225,  0.52403805, 0.02144266, 0.52403805, 0.02144266,    -4.6557898,   -1.4743335 ],
        [ 1.9311685,   0.52922975, 0.01855286, 0.52922975, 0.01855286,    -5.040301,    -1.4798067 ],
        [ 2.20485368,  0.53395719, 0.01607864, 0.53395719, 0.01607864,    -5.4459853,   -1.4842041 ],
        [ 2.49372664,  0.53828447, 0.01395425, 0.53828447, 0.01395425,    -5.8753114,   -1.4877414 ],
        [ 2.79930993,  0.54225978, 0.01212638, 0.54225978, 0.01212638,    -6.3304351,   -1.490569 ],
        [ 3.12349589,  0.54592641, 0.0105482,  0.54592641, 0.0105482,     -6.8140723,   -1.4928109 ],
        [ 3.46782258,  0.54931528, 0.00918223, 0.54931528, 0.00918223,    -7.3284379,   -1.4945843 ],
        [ 3.83420523,  0.55245478, 0.00799543, 0.55245478, 0.00799543,    -7.8763216,   -1.4959678 ],
        [ 4.2245456,   0.55536702, 0.00696039, 0.55536702, 0.00696039,    -8.4605,      -1.4970391 ],
        [ 4.64057971,  0.55806823, 0.00605459, 0.55806823, 0.00605459,    -9.0835197,   -1.4978611 ],
        [ 5.08525978,  0.56057759, 0.00525712, 0.56057759, 0.00525712,    -9.7497505,   -1.4984834 ],
        [ 5.56147782,  0.56290783, 0.00455164, 0.56305115, 0.0052679812,  -10.463487,   -1.4989517 ],
        [ 6.07192529,  0.56506639, 0.00392545, 0.56566216, 0.0048750225,  -11.228728,   -1.4992921 ],
        [ 6.62128396,  0.56706457, 0.00336657, 0.56818487, 0.0043702312,  -12.052457,   -1.4995328 ],
        [ 7.21562207,  0.56891183, 0.00286536, 0.5706288,  0.0038762902,  -12.943748,   -1.4997017 ],
        [ 7.86079443,  0.57061072, 0.00241551, 0.5729697,  0.0033947942,  -13.91136,    -1.4998174 ],
        [ 8.56804717,  0.57217095, 0.00201009, 0.57520351, 0.0029270638,  -14.972144,   -1.4998945 ],
        [ 9.35141577,  0.57359737, 0.00164439, 0.57731815, 0.002478877,   -16.147139,   -1.4999441 ],
        [ 10.2321278,  0.5748952,  0.00131523, 0.57930771, 0.0020513422,  -17.468176,   -1.4999775 ],
        [ 11.24260414, 0.57606846, 0.00101977, 0.58116694, 0.0016416069,  -18.983883,   -1.5000015 ],
        [ 12.43426056, 0.57711867, 0.0007561,  0.58287691, 0.001246419,   -20.771382,   -1.5000244 ],
        [ 13.89850945, 0.57804545, 0.00052399, 0.58444525, 0.00090870258, -22.967813,   -1.5000649 ],
        [ 15.81856214, 0.57884436, 0.00032415, 0.58584283, 0.00056393111, -25.848081,   -1.500179 ],
        [ 18.65223427, 0.57950245, 0.0001596,  0.58705166, 0.00029591245, -30.099433,   -1.5 ],
        [ 24.25257445, 0.57998344, 3.935e-05,  0.58794344, 7.2963742e-05, -38.499943,   -1.5 ]
    ],

    impulse_table => [
        [
            -7.43704, 13.815064,  -7.4376014,   0.032642155, 0,          -0.20466814,
            0,        -2.2588586, -0.003741227, -0.47347796, 0.27200436, 0.9943642,
            -3.5767029
        ],
        [
            -6.9077551, 12.756496, -6.9087083,    0.039034511, 0,          -0.20425717,
            0,          -2.258863, -0.0048746878, -0.47306699, 0.27200389, 0.99128886,
            -3.3664819
        ],
        [
            -6.50229, 11.945571,  -6.5037202,    0.0445263,   0,          -0.20375717,
            0,        -2.2588698, -0.0059702489, -0.47256699, 0.27200299, 0.98784208,
            -3.207358
        ],
        [
            -6.2146079, 11.370212,  -6.2165146,    0.048726557, 0,          -0.20325717,
            0,          -2.2588792, -0.0068938496, -0.47206699, 0.27200176, 0.98460053,
            -3.0957554
        ],
        [
            -5.8091428, 10.559299,  -5.8120052,    0.055037879, 0,          -0.20225717,
            0,          -2.2589056, -0.0084432069, -0.47106699, 0.27199836, 0.97852146,
            -2.9407673
        ],
        [
            -5.2983172, 9.5377014,  -5.3030904,   0.063508408, 0,          -0.20025717,
            0,          -2.2589902, -0.010900133, -0.46906699, 0.27198759, 0.96737221,
            -2.7504994
        ],
        [
            -4.60517, 8.1516614,  -4.6147416,   0.075344305, 0,          -0.19525717,
            0,        -2.2593868, -0.015415116, -0.46406699, 0.27193705, 0.94266473,
            -2.5051067
        ],
        [
            -3.9120228, 6.7663809,  -3.9312484,   0.086141754, 0,          -0.18525736,
            0,          -2.2609742, -0.021800262, -0.45406698, 0.27173517, 0.90010466,
            -2.2825988
        ],
        [
            -3.5065577, 5.9571337,  -3.535511,    0.090971436, 0,          -0.17525853,
            0,          -2.2636228, -0.026699729, -0.44406698, 0.27139989, 0.86280956,
            -2.1686349
        ],
        [
            -2.9957321, 4.9408107,  -3.0443258,   0.094167709, 0,          -0.15527353,
            0,          -2.2721191, -0.034468693, -0.42406698, 0.27033702, 0.79822299,
            -2.0503826
        ],
        [
            -2.6592599, 4.2757016,  -2.7276814,   0.093859773, 0,          -0.13534103,
            0,          -2.2849124, -0.040781203, -0.40422257, 0.26877061, 0.74316851,
            -1.9936386
        ],
        [
            -2.3025849, 3.5784144,  -2.4009456,   0.091212003, 0,          -0.10572966,
            0,          -2.3122406, -0.048724216, -0.37453389, 0.26554743, 0.67362791,
            -1.9580733
        ],
        [
            -2.0402207, 3.0742693,  -2.1685614,   0.088104303, 0,          -0.076928277,
            0,          -2.3493909, -0.055488214, -0.3454689,  0.26139735, 0.61603016,
            -1.9518874
        ],
        [
            -1.8971198, 2.8039088,  -2.0453725,   0.086301513, 0,          -0.05855566,
            0,          -2.3795778, -0.059506784, -0.32671808, 0.25819139, 0.58288489,
            -1.9567465
        ],
        [
            -1.6094378, 2.2738775,  -1.8068618,   0.083506798, 0,          -0.017543745,
            0,          -2.4731449, -0.068072545, -0.28427291, 0.24898281, 0.5142731,
            -1.9858418
        ],
        [
            -1.3862942, 1.8785558,  -1.6315841,   0.083111595, 0,         0.013924808,
            0,          -2.5890047, -0.074539974, -0.25132433, 0.2386235, 0.46104777,
            -2.0270056
        ],
        [
            -1.2039726,  1.5679773,   -1.4954728, 0.08431104, 0, 0.035792791, 0, -2.7198762,
            -0.07892445, -0.22941324, 0.22770634, 0.41885564, -2.0726252
        ],
        [
            -1.049822, 1.3150555,  -1.3857051,   0.086276772, 0,          0.050466602,
            0,         -2.8573424, -0.081565426, -0.21688551, 0.21668561, 0.38476236,
            -2.1191499
        ],
        [
            -0.91629057, 1.1035005,  -1.2946766,   0.08848277,  0,          0.060457794,
            0,           -2.9946769, -0.082992758, -0.21097576, 0.20588528, 0.3567379,
            -2.1648656
        ],
        [
            -0.79850753, 0.92281428, -1.2175386,   0.090661291, 0,         0.06751129,
            0,           -3.1278015, -0.083668165, -0.20895912, 0.1955219, 0.33334801,
            -2.2089637
        ],
        [
            -0.69314702, 0.76588226, -1.1510303,   0.092696164, 0,         0.072697279,
            0,           -3.254739,  -0.083902801, -0.20915227, 0.1857295, 0.31355967,
            -2.2511
        ],
        [
            -0.59783684, 0.62769191, -1.0928677,   0.094548081, 0,          0.076657098,
            0,           -3.3747877, -0.083885112, -0.21051552, 0.17658131, 0.29661526,
            -2.2911727
        ],
        [
            -0.51082546, 0.50459952, -1.041398,    0.096213603, 0,          0.07978155,
            0,           -3.4879231, -0.083725327, -0.21259387, 0.16810715, 0.28194956,
            -2.3292051
        ],
        [
            -0.35667478, 0.29346997, -0.95391298,  0.099038423, 0,          0.084419545,
            0,           -3.694799,  -0.083209614, -0.21806612, 0.15315972, 0.25783942,
            -2.3995206
        ],
        [
            -0.22314339, 0.11742303, -0.88179126,  0.10130651,  0,          0.087730318,
            0,           -3.8789092, -0.082612714, -0.22372889, 0.14068453, 0.23883776,
            -2.462971
        ],
        [
            -0.10536035, -0.032900755, -0.82084384,  0.10314768,  0,          0.090240857,
            0,           -4.0439164,   -0.082024572, -0.22951477, 0.13034486, 0.22346056,
            -2.5205202
        ],
        [
            1.6138916e-07, -0.16366094, -0.76832882,  0.10466305,  0,          0.092228904,
            0,             -4.1929533,  -0.081475568, -0.23508605, 0.12177462, 0.21074359,
            -2.5730238
        ],
        [
            0.095310341, -0.27910158, -0.72236719,  0.10592759,  0,         0.093854499,
            0,           -4.3285665,  -0.080973733, -0.24074422, 0.1146067, 0.20003546,
            -2.6212026
        ],
        [
            0.15196479, -0.34650257, -0.69570986,  0.10663264,  0,         0.094753572,
            0,          -4.4094226,  -0.080676751, -0.24406998, 0.1106022, 0.19400984,
            -2.6501002
        ],
        [
            0.26236443, -0.4753647, -0.64511566,  0.1079106,   0,          0.096380168,
            0,          -4.5673033, -0.080106747, -0.25101834, 0.10332069, 0.18295386,
            -2.7068726
        ],
        [
            0.3364724, -0.56011634, -0.61210904,  0.10870055,  0,           0.097389331,
            0,         -4.6734048,  -0.079733634, -0.25583541, 0.098798206, 0.1760126,
            -2.7452667
        ],
        [
            0.40546527, -0.63782065, -0.58203595,  0.10938964,  0,           0.098275934,
            0,          -4.772195,   -0.079394822, -0.26051952, 0.094834479, 0.16987619,
            -2.7811774
        ],
        [
            0.53062841,  -0.77600206, -0.52900505, 0.110533,   0, 0.099767956, 0, -4.9512851,
            -0.07880459, -0.26980437, 0.088205639, 0.15949249, -2.8466472
        ],
        [
            0.69314734,  -0.9504624,  -0.46286818, 0.11183108, 0, 0.10151079, 0, -5.1832539,
            -0.07808994, -0.28209155, 0.080568711, 0.14732092, -2.9320881
        ],
        [
            0.99325193, -1.2596471, -0.34785492,  0.11376165,  0,           0.10425901,
            0,          -5.6088273, -0.076933587, -0.30894499, 0.068860977, 0.12816654,
            -3.0904368
        ],
        [
            1.0986125, -1.3646765, -0.30940331,  0.11431938,  0,           0.10510313,
            0,         -5.7571812, -0.076577384, -0.31923779, 0.065362313, 0.12231352,
            -3.1460541
        ],
        [
            1.2089605,    -1.4729294, -0.27008626, 0.11484708, 0, 0.10592954, 0, -5.9119107,
            -0.076230939, -0.3303318, 0.06198512,  0.11660274, -3.204263
        ],
        [
            1.3862945, -1.6434609, -0.20877302,  0.11558872,  0,           0.10714608,
            0,         -6.1591369, -0.075728083, -0.34985018, 0.057099825, 0.10823062,
            -3.297659
        ],
        [
            1.6094381, -1.8526576, -0.13454082,  0.11636428,  0,          0.10850488,
            0,         -6.4676953, -0.075182378, -0.37784447, 0.05176003, 0.098921198,
            -3.4148273
        ],
        [
            1.7917596, -2.0196557, -0.076003589, 0.11689041,  0,           0.10949162,
            0,         -6.7177428, -0.07480135,  -0.40246793, 0.047953746, 0.09217918,
            -3.5102151
        ],
        [
            1.9459103, -2.158406,  -0.027817309, 0.1172724,   0,           0.11024923,
            0,         -6.9277565, -0.074519435, -0.42473958, 0.045064644, 0.087000532,
            -3.5906056
        ],
        [
            2.0794417, -2.2769595, 0.013054939,  0.11756343,  0,           0.11085462,
            0,         -7.1086841, -0.074301882, -0.44681254, 0.042773968, 0.082856085,
            -3.6600478
        ],
        [
            2.3025853, -2.4720208, 0.079750852,  0.11797951, 0,           0.11177145,
            0,         -7.4090775, -0.073988003, -0.4846551, 0.039325657, 0.076551839,
            -3.7756962
        ],
        [
            2.7080504, -2.8180857, 0.1966012,    0.11855772,  0,          0.11317762,
            0,         -7.9492057, -0.073547627, -0.56887281, 0.03404972, 0.066752305,
            -3.9846402
        ],
        [
            2.9957324, -3.0581766, 0.2767419,    0.11886263,  0,         0.11400474,
            0,         -8.3284959, -0.073316523, -0.63839893, 0.0309168, 0.060845117,
            -4.1320525
        ],
        [
            3.4011975, -3.3903924, 0.38664772,   0.11918591,  0,           0.11497741,
            0,         -8.8583742, -0.073075206, -0.75671059, 0.027159421, 0.053676016,
            -4.3388287
        ],
        [
            3.9120232,   -3.8007799,  0.52123522,  0.11946692,  0, 0.11594322, 0, -9.5194895,
            -0.07287616, -0.93588415, 0.023265978, 0.046155159, -4.5979974
        ],
        [
            4.2486966, -4.0673493, 0.60816642,   0.11959885, 0,           0.11645502,
            0,         -9.9520373, -0.072780326, -1.0924714, 0.021095522, 0.041614271,
            -4.7747324
        ],
        [
            4.6052008, -4.3469313, 0.69906385,   0.1197049,  0,           0.11690845,
            0,         -10.407838, -0.072706009, -1.2821278, 0.019070022, 0.037742065,
            -4.9529458
        ],
        [
            5.298358, -4.8844036, 0.87335853,   0.11984295, 0,          0.11758732,
            0,        -11.288941, -0.072478855, -1.7414575, 0.01576899, 0.031343101,
            -5.299542
        ],
        [
            6.2147534, -5.5863134, 1.1007843,    0.11994177, 0,           0.11819935,
            0,         -12.446536, -0.071688514, -2.6377954, 0.012370765, 0.024664909,
            -5.757807
        ],
        [
            6.9079728, -6.1130446, 1.2716891,    0.11998177, 0,           0.11850353,
            0,         -13.318659, -0.075283906, -3.6694474, 0.010337982, 0.020637957,
            -6.1044615
        ],
        [
            7.601065, -6.6373438, 1.4421803,    0.12000531, 0,            0.11871757,
            0,        -14.188675, -0.069702159, -5.0116875, 0.0086585064, 0.017298106,
            -6.4510402
        ],
        [
            8.5172329, -7.3280961, 1.6674741,    0.12002217, 0,            0.11890702,
            0,         -15.336815, -0.080281718, -7.742472,  0.0068641387, 0.013720745,
            -6.90915
        ],
        [
            9.2103915, -7.8496116, 1.8380903,    0.12002881, 0,            0.11900235,
            0,         -16.204582, -0.061669149, -9.7954917, 0.0057637951, 0.011523831,
            -7.2557421
        ],
        [
            10.819977, -9.0586939, 2.2352191,    0.12003433, 0,            0.1191294,
            0,         -18.218066, -0.098998638, -20.76652,  0.0038481374, 0.0076955256,
            -8.0605557
        ],
        [
            11.513019, -9.5788453, 2.4066479,   0.12003429, 0,            0.11915946,
            0,         -19.084659, -0.16954613, -29.113549, 0.0032349688, 0.006469564,
            -8.4070863
        ],
        [
            13.122485, -10.786359, 2.8056699,    0.12003016, 0,            0.11919956,
            0,         -21.096814, -0.013494784, -15.041991, 0.0021626007, 0.0043251295,
            -9.2118525
        ],
        [
            13.815707, -11.306354, 2.9778691,    0.12002632, 0,            0.11920904,
            0,         -21.963403, -0.011429835, -18.039476, 0.0018183696, 0.003636706,
            -9.5584871
        ],
        [
            16.11817, -13.033303, 3.5508788,     0.11999662, 0,            0.11922467,
            0,        -24.841557, -0.0065571547, -32.816234, 0.0010224882, 0.0020449862,
            -10.709891
        ],
        [
            18.420828, -14.760318, 4.1250555,     0.11990187, 0,             0.11922923,
            0,         -27.719885, -0.0037381596, -59.256387, 0.00057495639, 0.0011499551,
            -11.861766
        ]
    ],

    tail_shock_table => [
        [ 4.5222624, -1.3086,    -0.034304741,  -0.034304741 ],
        [ 4.5848776, -1.3360858, -0.039387989,  -0.028328211 ],
        [ 5.2863115, -1.6787191, -0.042318235,  -0.013417608 ],
        [ 6.208723,  -2.2285802, -0.037229913,  -0.0071326795 ],
        [ 6.9044004, -2.7379864, -0.032654666,  -0.0048536488 ],
        [ 7.5989481, -3.3468308, -0.028297088,  -0.0034451 ],
        [ 8.5161726, -4.3384448, -0.023179917,  -0.0022681532 ],
        [ 9.2097629, -5.2608678, -0.019838403,  -0.0016788173 ],
        [ 10.81979,  -8.157282,  -0.013687978,  -0.0008597608 ],
        [ 11.512908, -9.8224804, -0.011633121,  -0.0006500074 ],
        [ 13.122452, -15.065949, -0.0079509675, -0.00034076642 ],
        [ 13.815687, -18.052903, -0.006728562,  -0.00026117688 ],
        [ 16.118166, -32.81611,  -0.0038571516, -0.00010815058 ]
    ],

    blast_info => [
        0.20525711, 0,          0.47407847, -0.15415117, 0,          0,          1.2782794,   -0.18370548,
        0.66539115, 0.15731156, 0.76644388, -0.12061669, 0.31617329, 0.27200445, 0.070586245, 0,
        90.735,     -1.3086,    144.50898,  -1.5161767
    ],

};
1;
