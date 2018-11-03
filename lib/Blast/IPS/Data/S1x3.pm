package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=1.3
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S1x3'} = {

    table_name => 'S1x3',
    symmetry   => 2,
    gamma      => 1.3,

    shock_table_info => [ 2, 1.3, 8.9e-07, 8000, 1.99e-07, 10, 2.2e-07, 1.7e-07, 5e-07, 2103.2, -0.58791 ],

    shock_table => [
        [ -4.70223143,  12.00017664,   -2.99997917, -4.70328608, 0.99841721 ],
        [ -4.45655798,  11.263163,     -2.99995647, -4.45808291, 0.99771089 ],
        [ -3.96633083,  9.79252728,    -2.99982473, -3.96951469, 0.99521686 ],
        [ -3.59782065,  8.6871166,     -2.9994631,  -3.60336057, 0.99166845 ],
        [ -3.29907098,  7.79112733,    -2.99868401, -3.30775561, 0.98692171 ],
        [ -3.0503742,   7.04552293,    -2.99723068, -3.06300758, 0.98094626 ],
        [ -2.84019686,  6.41580123,    -2.99481177, -2.85754665, 0.97379094 ],
        [ -2.66316417,  5.88591324,    -2.99120708, -2.68583804, 0.96569558 ],
        [ -2.50091219,  5.40098954,    -2.98578028, -2.52989989, 0.95607931 ],
        [ -2.35433674,  4.96386969,    -2.97811598, -2.39053789, 0.94508483 ],
        [ -2.21961152,  4.56330574,    -2.96760331, -2.26402517, 0.93257617 ],
        [ -2.09375147,  4.19063983,    -2.9534948,  -2.14751722, 0.91836908 ],
        [ -1.97460992,  3.83980427,    -2.93492173, -2.03903302, 0.90226365 ],
        [ -1.85837973,  3.50003452,    -2.91038224, -1.93521383, 0.88367154 ],
        [ -1.74293336,  3.16582162,    -2.87808721, -1.83441451, 0.86202626 ],
        [ -1.62267206,  2.82222056,    -2.83428678, -1.73228513, 0.83578019 ],
        [ -1.491976,    2.45562314,    -2.77317689, -1.62515954, 0.80273422 ],
        [ -1.37431303,  2.13319451,    -2.70533846, -1.53265795, 0.76896018 ],
        [ -1.26536122,  1.8423499,     -2.63201884, -1.45072793, 0.7345355 ],
        [ -1.16231806,  1.57506888,    -2.5545468,  -1.37682103, 0.69960832 ],
        [ -1.06377994,  1.32724346,    -2.47471494, -1.30959833, 0.6645733 ],
        [ -0.96680947,  1.09123175,    -2.39255768, -1.24686805, 0.62911758 ],
        [ -0.87060488,  0.8650488,     -2.30946051, -1.18805265, 0.59358972 ],
        [ -0.77430626,  0.64665301,    -2.22652368, -1.13259559, 0.55827113 ],
        [ -0.67691079,  0.43381506,    -2.14449697, -1.0799316,  0.52335149 ],
        [ -0.57758884,  0.22484403,    -2.0640931,  -1.0296672,  0.48905549 ],
        [ -0.47544149,  0.01804122,    -1.98582576, -0.98143875, 0.45556359 ],
        [ -0.36973696,  -0.18782376,   -1.91022771, -0.9350202,  0.42309443 ],
        [ -0.25959251,  -0.39416833,   -1.83763286, -0.89016627, 0.39180135 ],
        [ -0.14420109,  -0.60215298,   -1.76836555, -0.8467118,  0.36184462 ],
        [ -0.02259609,  -0.81312368,   -1.7026069,  -0.80447462, 0.33332678 ],
        [ 0.10617532,   -1.02829361,   -1.64052132, -0.76332318, 0.30634216 ],
        [ 0.24323289,   -1.24905499,   -1.58217498, -0.72311458, 0.28093846 ],
        [ 0.38979197,   -1.47684534,   -1.52759502, -0.68372418, 0.25714082 ],
        [ 0.54720056,   -1.71320357,   -1.47676662, -0.64503688, 0.23495021 ],
        [ 0.71050332,   -1.95055887,   -1.43130134, -0.60833182, 0.21507382 ],
        [ 0.88081751,   -2.19076871,   -1.39051473, -0.57326062, 0.19721338 ],
        [ 1.05954105,   -2.4359198,    -1.35378355, -0.53949077, 0.18109427 ],
        [ 1.24792769,   -2.68775069,   -1.32063391, -0.50678425, 0.16650555 ],
        [ 1.44720883,   -2.94786448,   -1.29067892, -0.4749552,  0.15327361 ],
        [ 1.65881894,   -3.21804003,   -1.26357071, -0.44382702, 0.14124171 ],
        [ 1.883915,     -3.49963018,   -1.23905171, -0.41329843, 0.13029332 ],
        [ 2.12400764,   -3.79438033,   -1.21686365, -0.38324561, 0.12031188 ],
        [ 2.38037768,   -4.10370479,   -1.19680555, -0.35359899, 0.11120726 ],
        [ 2.65508639,   -4.4299132,    -1.17864849, -0.32422374, 0.10287711 ],
        [ 2.94960391,   -4.77456009,   -1.16223795, -0.29507684, 0.09525378 ],
        [ 3.26600246,   -5.13987596,   -1.1474046,  -0.2660737,  0.08826365 ],
        [ 3.60634462,   -5.52803912,   -1.13400423, -0.23715483, 0.08184555 ],
        [ 3.97328522,   -5.94186341,   -1.12189148, -0.20823474, 0.0759383 ],
        [ 4.36927013,   -6.38388462,   -1.11094941, -0.17927038, 0.07049503 ],
        [ 4.79733801,   -6.85726647,   -1.10106061, -0.15019768, 0.06546905 ],
        [ 5.26072023,   -7.36534479,   -1.09212084, -0.12096572, 0.06082052 ],
        [ 5.76344177,   -7.91228173,   -1.0840278,  -0.09150167, 0.05651039 ],
        [ 6.30952168,   -8.50218516,   -1.07669644, -0.06176375, 0.05250831 ],
        [ 6.90397417,   -9.14019253,   -1.07004314, -0.03168595, 0.04878399 ],
        [ 7.55200936,   -9.83159984,   -1.06399759, -0.00122611, 0.04531353 ],
        [ 8.25923415,   -10.58208161,  -1.05849797, 0.02964619,  0.04207675 ],
        [ 9.03245304,   -11.39853555,  -1.05348521, 0.06098122,  0.03905358 ],
        [ 9.87846879,   -12.28780751,  -1.0489118,  0.09279497,  0.03622933 ],
        [ 10.80468299,  -13.25733209,  -1.04473565, 0.12509689,  0.03359102 ],
        [ 11.81949648,  -14.31554853,  -1.04091827, 0.15790157,  0.03112626 ],
        [ 12.9315097,   -15.47106412,  -1.03742765, 0.19120075,  0.02882521 ],
        [ 14.15032284,  -16.73349001,  -1.03423441, 0.22498974,  0.02667826 ],
        [ 15.48573596,  -18.11261135,  -1.03131364, 0.2592441,   0.02467736 ],
        [ 16.9489491,   -19.61962937,  -1.02864154, 0.29395207,  0.02281398 ],
        [ 18.46398782,  -21.1762526,   -1.02632006, 0.32723756,  0.02116827 ],
        [ 20.16438713,  -22.91948902,  -1.02412869, 0.36185492,  0.01959058 ],
        [ 21.94999295,  -24.74640072,  -1.02219122, 0.39553974,  0.0181752 ],
        [ 23.96945697,  -26.8087615,   -1.02034592, 0.43082668,  0.01680833 ],
        [ 26.07699099,  -28.95741732,  -1.01872299, 0.46493412,  0.01559024 ],
        [ 28.459083,    -31.38221348,  -1.01717606, 0.50063369,  0.0144146 ],
        [ 31.16445423,  -34.13200129,  -1.01570398, 0.53805551,  0.01328191 ],
        [ 33.97931731,  -36.98920907,  -1.01441927, 0.57399703,  0.01228178 ],
        [ 37.17273303,  -40.22665963,  -1.01319558, 0.61163703,  0.01131854 ],
        [ 40.46265958,  -43.55820572,  -1.01213533, 0.64745032,  0.01047519 ],
        [ 44.18697233,  -47.32577652,  -1.01112416, 0.68490902,  0.00966289 ],
        [ 47.97483625,  -51.15408741,  -1.01025553, 0.72014251,  0.00895858 ],
        [ 52.24786766,  -55.46911805,  -1.00942568, 0.75693384,  0.00827982 ],
        [ 57.09122872,  -60.35615806,  -1.00863415, 0.79541048,  0.00762678 ],
        [ 61.95925407,  -65.2645353,   -1.00796237, 0.83114276,  0.00706803 ],
        [ 67.44996457,  -70.79715241,  -1.00732021, 0.86843414,  0.00652985 ],
        [ 73.67273999,  -77.06351736,  -1.00670735, 0.90741278,  0.00601235 ],
        [ 79.82307523,  -83.25348737,  -1.0061949,  0.94301661,  0.00557659 ],
        [ 86.74613536,  -90.21769071,  -1.00570444, 0.98013215,  0.00515681 ],
        [ 94.57538072,  -98.08971042,  -1.00523578, 1.01888214,  0.00475309 ],
        [ 103.47423954, -107.03311322, -1.0047887,  1.05940501,  0.00436548 ],
        [ 112.10488614, -115.70346731, -1.00442252, 1.0956684,   0.00404611 ],
        [ 121.82443781, -125.46425318, -1.00407193, 1.13346087,  0.00373867 ],
        [ 132.8221657,  -136.50486741, -1.0037368,  1.17290781,  0.00344318 ],
        [ 143.12921188, -146.84901586, -1.00346923, 1.2071457,   0.0032061 ],
        [ 154.646881,   -158.40512454, -1.00321225, 1.24272245,  0.00297737 ],
        [ 167.57073757, -171.36886052, -1.00296576, 1.27974066,  0.002757 ],
        [ 182.13734906, -185.97690621, -1.0027297,  1.31831522,  0.00254501 ],
        [ 198.63503301, -202.51770847, -1.00250399, 1.35857544,  0.00234141 ],
        [ 217.41803654, -221.34566066, -1.00228857, 1.40066755,  0.00214623 ]
    ],

    energy_table => [
        [ -4.70223143,  0.06019761, 0.04165331, 0.06019761, 0.04165331,    5.6793741,    -1.6004213 ],
        [ -4.45655798,  0.07134978, 0.0493506,  0.07134978, 0.0493506,     5.2843477,    -1.6174258 ],
        [ -3.96633083,  0.10012459, 0.06912695, 0.10012459, 0.06912695,    4.4821554,    -1.6612462 ],
        [ -3.59782065,  0.12908362, 0.08881618, 0.12908362, 0.08881618,    3.8630778,    -1.7051851 ],
        [ -3.29907098,  0.15845979, 0.10841937, 0.15845979, 0.10841937,    3.3475418,    -1.7549217 ],
        [ -3.0503742,   0.18773331, 0.12739569, 0.18773331, 0.12739569,    2.9050373,    -1.8078087 ],
        [ -2.84019686,  0.2163539,  0.14519624, 0.2163539,  0.14519624,    2.5200118,    -1.86313 ],
        [ -2.66316417,  0.24346079, 0.16115172, 0.24346079, 0.16115172,    2.1855212,    -1.9219545 ],
        [ -2.50091219,  0.27082203, 0.17612903, 0.27082203, 0.17612903,    1.8688435,    -1.9865302 ],
        [ -2.35433674,  0.29762048, 0.18944914, 0.29762048, 0.18944914,    1.5730636,    -2.0580563 ],
        [ -2.21961152,  0.32393101, 0.20096383, 0.32393101, 0.20096383,    1.2908233,    -2.1331698 ],
        [ -2.09375147,  0.34984079, 0.2105138,  0.34984079, 0.2105138,     1.0178466,    -2.2055626 ],
        [ -1.97460992,  0.37538074, 0.21790144, 0.37538074, 0.21790144,    0.75093658,   -2.2847567 ],
        [ -1.85837973,  0.40102601, 0.22298556, 0.40102601, 0.22298556,    0.48033528,   -2.376205 ],
        [ -1.74293336,  0.42693638, 0.22540733, 0.42693638, 0.22540733,    0.20050066,   -2.4536998 ],
        [ -1.62267206,  0.45403313, 0.22462454, 0.45403313, 0.22462454,    -0.098313168, -2.5181925 ],
        [ -1.491976,    0.483107,   0.21952512, 0.483107,   0.21952512,    -0.43218746,  -2.5353 ],
        [ -1.37431303,  0.50847432, 0.21106091, 0.50847432, 0.21106091,    -0.72845545,  -2.4491669 ],
        [ -1.26536122,  0.53090166, 0.20017667, 0.53090166, 0.20017667,    -0.98835868,  -2.2746858 ],
        [ -1.16231806,  0.55089783, 0.18761204, 0.55089783, 0.18761204,    -1.2119513,   -2.0510322 ],
        [ -1.06377994,  0.56872678, 0.17405244, 0.56872678, 0.17405244,    -1.4028559,   -1.8561756 ],
        [ -0.96680947,  0.58491837, 0.15979189, 0.58491837, 0.15979189,    -1.5751021,   -1.7374471 ],
        [ -0.87060488,  0.59959428, 0.14529074, 0.59959428, 0.14529074,    -1.7385466,   -1.6915088 ],
        [ -0.77430626,  0.61289061, 0.13091814, 0.61289061, 0.13091814,    -1.9007216,   -1.6938237 ],
        [ -0.67691079,  0.62495494, 0.11694442, 0.62495494, 0.11694442,    -2.0666517,   -1.7207991 ],
        [ -0.57758884,  0.63589874, 0.10359898, 0.63589874, 0.10359898,    -2.2392993,   -1.7577492 ],
        [ -0.47544149,  0.64582919, 0.09104476, 0.64582919, 0.09104476,    -2.4208959,   -1.7970119 ],
        [ -0.36973696,  0.65482599, 0.07941511, 0.65482599, 0.07941511,    -2.6129508,   -1.834528 ],
        [ -0.25959251,  0.66297379, 0.06878068, 0.66297379, 0.06878068,    -2.8170367,   -1.8683395 ],
        [ -0.14420109,  0.67034212, 0.05918201, 0.67034212, 0.05918201,    -3.0344943,   -1.8976219 ],
        [ -0.02259609,  0.67700288, 0.05061426, 0.67700288, 0.05061426,    -3.2669335,   -1.9222531 ],
        [ 0.10617532,   0.68301818, 0.0430504,  0.68301818, 0.0430504,     -3.5159409,   -1.9424407 ],
        [ 0.24323289,   0.68844992, 0.03643606, 0.68844992, 0.03643606,    -3.7834397,   -1.9585896 ],
        [ 0.38979197,   0.69335469, 0.03070305, 0.69335469, 0.03070305,    -4.0715644,   -1.971192 ],
        [ 0.54720056,   0.6977849,  0.02577382, 0.6977849,  0.02577382,    -4.3827368,   -1.9807823 ],
        [ 0.71050332,   0.70164919, 0.0217088,  0.70164919, 0.0217088,     -4.7068711,   -1.9876192 ],
        [ 0.88081751,   0.70504892, 0.0183441,  0.70504892, 0.0183441,     -5.0458797,   -1.9924411 ],
        [ 1.05954105,   0.70806763, 0.01554526, 0.70806763, 0.01554526,    -5.4023419,   -1.9957442 ],
        [ 1.24792769,   0.71076753, 0.01320908, 0.71076753, 0.01320908,    -5.7785626,   -1.9978919 ],
        [ 1.44720883,   0.7131974,  0.01125356, 0.7131974,  0.01125356,    -6.1768787,   -1.9992698 ],
        [ 1.65881894,   0.71539821, 0.00961128, 0.71539821, 0.00961128,    -6.6000579,   -2.0000843 ],
        [ 1.883915,     0.71740013, 0.00822996, 0.71740013, 0.00822996,    -7.0503356,   -2.0005149 ],
        [ 2.12400764,   0.71923075, 0.00706486, 0.71923075, 0.00706486,    -7.530679,    -2.0006893 ],
        [ 2.38037768,   0.72091089, 0.00608062, 0.72091089, 0.00608062,    -8.0436042,   -2.0007031 ],
        [ 2.65508639,   0.72246216, 0.00524577, 0.72246216, 0.00524577,    -8.5932092,   -2.0006532 ],
        [ 2.94960391,   0.72389867, 0.00453671, 0.72389867, 0.00453671,    -9.1824273,   -2.0005581 ],
        [ 3.26600246,   0.72523483, 0.00393261, 0.72523483, 0.00393261,    -9.8153794,   -2.0004371 ],
        [ 3.60634462,   0.72648209, 0.00341663, 0.72648209, 0.00341663,    -10.496193,   -2.0003335 ],
        [ 3.97328522,   0.72765155, 0.00297433, 0.72765155, 0.00297433,    -11.230178,   -2.0002412 ],
        [ 4.36927013,   0.72875125, 0.00259429, 0.72875125, 0.00259429,    -12.022226,   -2.0001647 ],
        [ 4.79733801,   0.72978903, 0.00226671, 0.72978903, 0.00226671,    -12.878417,   -2.0001084 ],
        [ 5.26072023,   0.73077132, 0.00198345, 0.73077132, 0.00198345,    -13.805221,   -2.000069 ],
        [ 5.76344177,   0.73170442, 0.0017378,  0.73170442, 0.0017378,     -14.81069,    -2.0000419 ],
        [ 6.30952168,   0.73259289, 0.00152436, 0.73259289, 0.00152436,    -15.902867,   -2.0000237 ],
        [ 6.90397417,   0.73344146, 0.00133802, 0.73344146, 0.00133802,    -17.091782,   -2.0000132 ],
        [ 7.55200936,   0.73425359, 0.00117459, 0.73425359, 0.00117459,    -18.387858,   -2.000007 ],
        [ 8.25923415,   0.73503189, 0.00103164, 0.73504965, 0.001095559,   -19.802311,   -2.0000034 ],
        [ 9.03245304,   0.73577928, 0.00090614, 0.73588869, 0.0010714523,  -21.34875,    -2.0000016 ],
        [ 9.87846879,   0.73649754, 0.00079586, 0.73678042, 0.0010419344,  -23.040783,   -2.0000007 ],
        [ 10.80468299,  0.73718813, 0.00069889, 0.73772601, 0.00099454736, -24.893211,   -2.0000003 ],
        [ 11.81949648,  0.7378525,  0.00061357, 0.73871499, 0.00094904301, -26.922839,   -2.0000002 ],
        [ 12.9315097,   0.73849153, 0.0005385,  0.73973985, 0.00089090566, -29.146865,   -2.0000001 ],
        [ 14.15032284,  0.73910613, 0.00047246, 0.74083785, 0.00097201489, -31.584491,   -2 ],
        [ 15.48573596,  0.73969687, 0.00041439, 0.74205212, 0.00085106236, -34.255318,   -2 ],
        [ 16.9489491,   0.74026451, 0.00036336, 0.74321698, 0.00074507498, -37.181744,   -2 ],
        [ 18.46398782,  0.74078161, 0.00032075, 0.74427657, 0.0006567961,  -40.211822,   -2 ],
        [ 20.16438713,  0.74129294, 0.00028207, 0.74532296, 0.0005768641,  -43.61262,    -2 ],
        [ 21.94999295,  0.74176623, 0.0002492,  0.74629038, 0.00050907531, -47.183832,   -2 ],
        [ 23.96945697,  0.742238,   0.00021911, 0.74725364, 0.00044714498, -51.22276,    -2 ],
        [ 26.07699099,  0.74267204, 0.00019368, 0.74813905, 0.00039490012, -55.437828,   -2 ],
        [ 28.459083,    0.74310466, 0.00017039, 0.74902079, 0.00034712783, -60.202012,   -2 ],
        [ 31.16445423,  0.7435358,  0.00014913, 0.7498988,  0.0003035888,  -65.612754,   -2 ],
        [ 33.97931731,  0.74392965, 0.00013134, 0.75070032, 0.00026719554, -71.242481,   -2 ],
        [ 37.17273303,  0.74432216, 0.00011509, 0.7514986,  0.00023398625, -77.629312,   -2 ],
        [ 40.46265958,  0.74467779, 0.00010158, 0.75222146, 0.0002064136,  -84.209165,   -2 ],
        [ 44.18697233,  0.74503224, 8.922e-05,  0.75294156, 0.00018120376, -91.657791,   -2 ],
        [ 47.97483625,  0.74535021, 7.902e-05,  0.75358726, 0.00016043244, -99.233518,   -2 ],
        [ 52.24786766,  0.74566718, 6.967e-05,  0.75423067, 0.00014138611, -107.77958,   -2 ],
        [ 57.09122872,  0.74598312, 6.111e-05,  0.75487176, 0.00012397785, -117.4663,    -2 ],
        [ 61.95925407,  0.74626308, 5.415e-05,  0.75543964, 0.0001098101,  -127.20235,   -2 ],
        [ 67.44996457,  0.74654218, 4.775e-05,  0.75600563, 9.6810764e-05, -138.18378,   -2 ],
        [ 73.67273999,  0.74682042, 4.19e-05,   0.75656973, 8.4922034e-05, -150.62933,   -2 ],
        [ 79.82307523,  0.74706316, 3.72e-05,   0.75706173, 7.5385975e-05, -162.93,      -2 ],
        [ 86.74613536,  0.74730521, 3.288e-05,  0.75755224, 6.6619444e-05, -176.77612,   -2 ],
        [ 94.57538072,  0.74754657, 2.892e-05,  0.75804125, 5.8585847e-05, -192.43461,   -2 ],
        [ 103.47423954, 0.74778722, 2.53e-05,   0.75852874, 5.1248805e-05, -210.23232,   -2 ],
        [ 112.10488614, 0.74799291, 2.246e-05,  0.75894536, 4.5487407e-05, -227.49362,   -2 ],
        [ 121.82443781, 0.74819807, 1.985e-05,  0.75936085, 4.0189619e-05, -246.93272,   -2 ],
        [ 132.8221657,  0.74840269, 1.745e-05,  0.75977519, 3.5333621e-05, -268.92818,   -2 ],
        [ 143.12921188, 0.74857279, 1.561e-05,  0.76011958, 3.1608804e-05, -289.54227,   -2 ],
        [ 154.646881,   0.7487425,  1.391e-05,  0.76046317, 2.8163733e-05, -312.57761,   -2 ],
        [ 167.57073757, 0.74891182, 1.235e-05,  0.76080593, 2.4986043e-05, -338.42532,   -2 ],
        [ 182.13734906, 0.74908075, 1.09e-05,   0.76114788, 2.2063873e-05, -367.55854,   -2 ],
        [ 198.63503301, 0.74924928, 9.58e-06,   0.76148899, 1.9385338e-05, -400.55391,   -2 ],
        [ 217.41803654, 0.74941741, 8.37e-06,   0.76182927, 1.6938854e-05, -438.11992,   -2 ]
    ],

    impulse_table => [
        [
            -5.3071211, 13.81484,    -5.3075467,    0.01874984,  -0.00053490862, -0.15243899,
            -1.3237827, -0.32947529, -0.0010537376, -0.34429109, 0.18444439,     0.9485123,
            -1.9272291
        ],
        [
            -4.9618451, 12.779014,   -4.9625595,    0.021691793, -0.00075549406, -0.15039516,
            -1.3217389, -0.50211597, -0.0014882776, -0.34224726, 0.18444414,     0.93461248,
            -1.8579618
        ],
        [
            -4.6051702, 11.708995,   -4.6063902,    0.025056685, -0.0010792772, -0.14739516,
            -1.3187389, -0.68045966, -0.0021261108, -0.33924726, 0.18444329,    0.91630744,
            -1.7899171
        ],
        [
            -4.1997051, 10.49262,    -4.2019458,    0.029222297, -0.0016189158, -0.14239516,
            -1.3137389, -0.88321355, -0.0031891662, -0.33424726, 0.18444047,    0.88921872,
            -1.718501
        ],
        [
            -3.912023,  9.6296142,  -3.9154773,    0.032316759, -0.0021585544, -0.13739516,
            -1.3087389, -1.0270959, -0.0042522216, -0.32924726, 0.18443505,    0.86486557,
            -1.673025
        ],
        [
            -3.5065579, 8.4133845,  -3.5129125,    0.036652687, -0.0032378317, -0.12739516,
            -1.2987389, -1.2299979, -0.0063783325, -0.31924726, 0.18441281,    0.82132121,
            -1.6191434
        ],
        [
            -2.9957323, 6.8817614,  -3.0094509,   0.041350661, -0.0053963861, -0.10739516,
            -1.2787389, -1.4862849, -0.010630554, -0.29924726, 0.18429827,    0.74672387,
            -1.5764117
        ],
        [
            -2.3025851, 4.8098371,  -2.3417431,   0.043325659, -0.010792772, -0.057395183,
            -1.2287389, -1.8406763, -0.021261108, -0.24924726, 0.18328714,   0.60003647,
            -1.5973744
        ],
        [
            -1.89712,   3.6129599, -1.9695752,   0.039882532, -0.016189092, -0.0073976519,
            -1.1787403, -2.064753, -0.031891576, -0.19924726, 0.18064278,   0.48927731,
            -1.6845287
        ],
        [
            -1.6094379, 2.7847476,  -1.721245,    0.035068768, -0.021582738, 0.042526435,
            -1.1287834, -2.2508392, -0.042518588, -0.14938051, 0.17586626,   0.40453469,
            -1.7996374
        ],
        [
            -1.3862944, 2.165653,   -1.5418928,   0.031071141, -0.026931667, 0.091491407,
            -1.0793792, -2.4347021, -0.053087193, -0.0997805,  0.16887458,   0.33985294,
            -1.9275066
        ],
        [
            -1.2039728, 1.682152,   -1.4062634,   0.02920476,   -0.031902385, 0.13431655,
            -1.0338942, -2.6408175, -0.063124052, -0.053800401, 0.15999686,   0.29048445,
            -2.0590784
        ],
        [
            -1.0498221, 1.2927832,  -1.3003577,  0.02935832,   -0.035511136, 0.16244662,
            -1.0004802, -2.8857408, -0.07080002, -0.020667241, 0.14982432,   0.25259622,
            -2.1886277
        ],
        [
            -0.91629075, 0.97146227, -1.2155569,   0.03038344,    -0.037367408, 0.17815999,
            -0.98026365, -3.1590453, -0.074678226, -0.0055938422, 0.13902177,   0.22322036,
            -2.3128572
        ],
        [
            -0.79850771, 0.70078896, -1.1462133,   0.031529084,   -0.038226475, 0.18781961,
            -0.96736825, -3.4338201, -0.076182457, -0.0011041754, 0.12818899,   0.20014281,
            -2.4301848
        ],
        [
            -0.6931472,  0.46874345, -1.0884756,   0.032547225,    -0.038668483, 0.19458209,
            -0.95819097, -3.695119,  -0.076732056, -0.00061764547, 0.11779306,   0.18174826,
            -2.5401096
        ],
        [
            -0.59783702, 0.26680069, -1.0396391,   0.033391314,   -0.038929129, 0.19974508,
            -0.95116149, -3.938688,  -0.076909871, -0.0018053015, 0.10815377,   0.16686895,
            -2.6427381
        ],
        [
            -0.51082564, 0.088776176, -0.99775911,  0.034075499,   -0.039100498, 0.20389964,
            -0.94553051, -4.1644584,  -0.076931733, -0.0038174981, 0.09945763,   0.15466039,
            -2.738479
        ],
        [
            -0.43078293, -0.06991038, -0.96140814,  0.034627384,   -0.039222396, 0.20735507,
            -0.9408879,  -4.3737075,  -0.076883978, -0.0062727092, 0.09178353,   0.1445081,
            -2.8278661
        ],
        [
            -0.35667496, -0.21271666, -0.92951892,  0.035073722,   -0.039314124, 0.2102938,
            -0.93698395, -4.5680455,  -0.076804447, -0.0088778981, 0.085130359,  0.13596055,
            -2.9114628
        ],
        [
            -0.28768209, -0.34229944, -0.90127997,  0.035436822,  -0.039386024, 0.21283353,
            -0.93365238, -4.7490604,  -0.076711139, -0.011226458, 0.079442194,  0.12868204,
            -2.989814
        ],
        [
            -0.22314357, -0.46073377, -0.87606451,  0.035734351,  -0.039444131, 0.21505564,
            -0.93077743, -4.9182006,  -0.076613053, -0.013988841, 0.074620707,  0.12241994,
            -3.0634245
        ],
        [
            -0.16251895, -0.56966406, -0.8533817,  0.035980015,  -0.039492207, 0.21701904,
            -0.92827354, -5.0767472,  -0.07651487, -0.016510592, 0.070509502,  0.11698164,
            -3.1327512
        ],
        [
            -0.15177028, -0.58875151, -0.84945776,  0.036020387,  -0.039500158, 0.21735528,
            -0.92784827, -5.1048134,  -0.076497039, -0.016985455, 0.069819612,  0.11605978,
            -3.1450547
        ],
        [
            -0.051293311, -0.76405122, -0.81413237,  0.036355703,  -0.039567416, 0.22033699,
            -0.92413484,  -5.3663854,  -0.076326808, -0.021623116, 0.063885077,  0.10801451,
            -3.260141
        ],
        [
            -1.69756e-08, -0.85146632, -0.79699902,  0.036500274,  -0.03959748, 0.22175286,
            -0.92241079,  -5.4992924,  -0.076238759, -0.024010427, 0.061181576, 0.10427794,
            -3.3188877
        ],
        [
            0.095310163, -1.0104422, -0.76666332,  0.036728048, -0.039647099, 0.22420815,
            -0.91949041, -5.7449685, -0.076075853, -0.02886267, 0.05666469,   0.09792861,
            -3.4279097
        ],
        [
            0.18232154,  -1.1519287, -0.74055532,  0.036896257,  -0.039686486, 0.22626478,
            -0.91712249, -5.9676776, -0.075930026, -0.033169311, 0.053043742,  0.09273687,
            -3.5271741
        ],
        [
            0.26236425,  -1.2792517, -0.71777142,  0.037023046,  -0.039718581, 0.22801295,
            -0.91517487, -6.1711481, -0.075799721, -0.037488328, 0.050077079,  0.088411343,
            -3.6181855
        ],
        [
            0.40546509,  -1.5007451, -0.67971247,  0.037196207, -0.039767977, 0.2308264,
            -0.91218774, -6.5314641, -0.075578397, -0.04548062, 0.045504401,  0.081608117,
            -3.7800169
        ],
        [
            0.53062823,  -1.6886888, -0.6489486,   0.037304027,  -0.039804393, 0.23299176,
            -0.91003367, -6.8429638, -0.075398837, -0.052418799, 0.042139879,  0.076488597,
            -3.9205221
        ],
        [
            0.69314716,  -1.925678,  -0.61208175,  0.037400184,  -0.039844315, 0.23544247,
            -0.90777857, -7.2424447, -0.075186531, -0.062021162, 0.038473244,  0.070788907,
            -4.1014145
        ],
        [
            0.99325176,  -2.3457531, -0.55168224,  0.037494408,  -0.039901777, 0.23907133,
            -0.90487621, -7.9661539, -0.074855595, -0.079999756, 0.033254939,  0.062433935,
            -4.4308423
        ],
        [
            1.0986123,   -2.4886704, -0.5324782,   0.037511357,  -0.039918103, 0.24011428,
            -0.90415678, -8.2162549, -0.074757005, -0.086107964, 0.031784073,  0.060021215,
            -4.5451429
        ],
        [
            1.2089603,   -2.6361645, -0.51332745,  0.037523383,  -0.039933486, 0.24109747,
            -0.90353262, -8.476154,  -0.074662866, -0.092954353, 0.030400948,  0.057726846,
            -4.6641499
        ],
        [
            1.3862943,   -2.8689822, -0.4844074,   0.037533915, -0.039954977, 0.24247063,
            -0.90275714, -8.8897375, -0.074529435, -0.10471362, 0.028461909,  0.054465735,
            -4.853989
        ],
        [
            1.6094379,   -3.155498, -0.45086645,  0.037537455, -0.039977352, 0.2438909,
            -0.90208333, -9.403641, -0.074388875, -0.11919991, 0.026422944,  0.050975766,
            -5.0906255
        ],
        [
            1.7917595,   -3.3850096, -0.42550047,  0.037535768, -0.039992387, 0.24484007,
            -0.90171316, -9.8187063, -0.074293945, -0.12978325, 0.025018902,  0.048533076,
            -5.2823359
        ],
        [
            1.9459101,   -3.5762556, -0.40530597, 0.037532734, -0.040003298, 0.24551979,
            -0.90149019, -10.166606, -0.07422511, -0.14099983, 0.023979675,  0.046702735,
            -5.4434127
        ],
        [
            2.0794415,   -3.7400643, -0.38864588,  0.03752934,  -0.040011438, 0.24602975,
            -0.90134953, -10.465926, -0.074173636, -0.14820622, 0.023171502,  0.045265451,
            -5.5822736
        ],
        [
            2.3027867,   -4.0106237, -0.36232785,  0.037523278, -0.040023404, 0.24674619,
            -0.90120581, -10.962735, -0.074101001, -0.16326519, 0.021979463,  0.042684111,
            -5.8220714
        ],
        [
            2.995773,    -4.8281656, -0.29070429, 0.037506315, -0.040046523, 0.24817994,
            -0.90095422, -12.479296, -0.07395341, -0.20828244, 0.019196215,  0.037828312,
            -6.5266364
        ],
        [
            3.9122097,   -5.8732866, -0.21290063,  0.037493304, -0.040066424, 0.2490466,
            -0.90130718, -14.443834, -0.073868828, -0.26802474, 0.016787511,  0.033372045,
            -7.4514111
        ],
        [
            4.6054328,   -6.6455686, -0.1629658,   0.037488586, -0.040068616, 0.24933299,
            -0.90114001, -15.909525, -0.073836793, -0.3132709,  0.015497845,  0.030899568,
            -8.1478441
        ],
        [
            5.2983316,   -7.4064086, -0.11868469,  0.03748619,  -0.04007103, 0.24947699,
            -0.90114857, -17.362639, -0.073821881, -0.35851041, 0.014484414, 0.028922665,
            -8.8425055
        ],
        [
            6.2147105,   -8.4000468, -0.06677278,  0.037484758, -0.040072473, 0.24956346,
            -0.90115417, -19.271368, -0.073812925, -0.4183456,  0.013422493,  0.026827249,
            -9.7600396
        ],
        [
            6.9077744,   -9.1442589, -0.031500599, 0.037484289, -0.040072953, 0.24959228,
            -0.90115611, -20.707614, -0.07380994,  -0.46359919, 0.0127657,    0.025522737,
            -10.45352
        ],
        [
            7.6010114,   -9.8837278, 0.00098844088, 0.037484058, -0.040073192, 0.2496067,
            -0.90115709, -22.13935,  -0.073808448,  -0.50886326, 0.012200472,  0.024396696,
            -11.146978
        ],
        [
            8.5172412,   -10.85495, 0.040364172,  0.037483922, -0.040073335, 0.24961535,
            -0.90115769, -24.02567, -0.073807553, -0.5686863,  0.011560574,  0.023119484,
            -12.063348
        ],
        [
            9.2104567,   -11.585967, 0.067876414,  0.037483877, -0.040073382, 0.24961823,
            -0.90115789, -25.449163, -0.073807255, -0.61394759, 0.011139907,  0.022278994,
            -12.756613
        ],
        [
            10.819883,   -13.273212, 0.12560718,   0.037483841, -0.040073417, 0.24962054,
            -0.90115804, -28.744455, -0.073791568, -0.71392618, 0.010320078,  0.020639996,
            -14.366081
        ],
        [
            11.513093,   -13.996443, 0.14825763,   0.037483835, -0.040073424, 0.24962083,
            -0.90115806, -30.16042,  -0.073478126, -0.73886693, 0.010019616,  0.020039153,
            -15.059296
        ],
        [
            13.122512,   -15.669163, 0.19667197,   0.037483829, -0.040073423, 0.24962106,
            -0.90115808, -33.441658, -0.071843509, -0.79404399, 0.0094135,    0.018826984,
            -16.66872
        ],
        [
            13.815519,   -16.387089, 0.21596555,   0.037483828, -0.040073419, 0.24962108,
            -0.90115809, -34.852272, -0.070947134, -0.81677281, 0.0091847657, 0.018369524,
            -17.361728
        ],
        [
            16.118142,  -18.764432, 0.27458079,   0.037483824, -0.040073422, 0.24962109,
            -0.9011581, -39.531389, -0.067733552, -0.88857173, 0.0085307193, 0.017061438,
            -19.664351
        ],
        [
            18.420761,   -21.131887, 0.32632159,   0.03748382,  -0.04005021,  0.24962156,
            -0.95559409, -44.200847, -0.064565911, -0.95559409, 0.0080001098, 0.01600022,
            -21.966971
        ]
    ],

    tail_shock_table => [
        [ 7.6514948, -0.58791,    -0.03688871,  -0.03688871 ],
        [ 8.5170329, -0.62466015, -0.052972475, -0.018838099 ],
        [ 9.2103497, -0.65271086, -0.055478299, -0.013486481 ],
        [ 10.819861, -0.71392481, -0.056762786, -0.0067950438 ],
        [ 11.513081, -0.73886625, -0.056521664, -0.0050962012 ],
        [ 13.122509, -0.79404385, -0.055264245, -0.0025013359 ],
        [ 13.815518, -0.81677274, -0.054574722, -0.0017573121 ],
        [ 16.118142, -0.88858654, -0.052104336, -0.00018622054 ]
    ],

    blast_info => [
        0.15739511, 1.3287381,  0.34930477, -0.21261118, 1.6151814,  0.0027799119, 3.6048765,  -0.10792771,
        0.43148504, 0.40179621, 0.58110032, -0.14573456, 0.39270413, 0.18444438,   0.02883371, -0.030825705,
        2103.2,     -0.58791,   40181.464,  -0.70474841
    ],

};
1;
