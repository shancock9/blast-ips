package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=1.667
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S1.667'} = {

    shock_table_info => [ 2, 1.667, 8.2e-07, 8000, 2.14e-07, 10, 1.2e-07, 2e-07, 5e-07, 1719, -0.70345 ],

    shock_table => [
        [ -4.47131967,  12.00017656,   -2.99997696, -4.47242876, 0.99833545 ],
        [ -4.34086863,  11.60882679,   -2.99995356, -4.3422176,  0.99797521 ],
        [ -4.18432739,  11.13920907,   -2.99992945, -4.18603368, 0.99743842 ],
        [ -4.00830247,  10.61114536,   -2.99990805, -4.01052493, 0.9966627 ],
        [ -3.90394164,  10.29807267,   -2.99990634, -3.90654119, 0.99609575 ],
        [ -3.80898673,  10.01321978,   -2.99989199, -3.81198478, 0.99549639 ],
        [ -3.75423088,  9.84896131,    -2.99982618, -3.75748598, 0.99510968 ],
        [ -3.36690888,  8.68712263,    -2.99944452, -3.37273544, 0.99123623 ],
        [ -3.07655735,  7.8163253,     -2.99867172, -3.08557716, 0.98641492 ],
        [ -2.82203786,  7.05327322,    -2.99715821, -2.83527547, 0.98003001 ],
        [ -2.61566768,  6.43497396,    -2.99473421, -2.63374459, 0.97268482 ],
        [ -2.43594758,  5.89706555,    -2.99100619, -2.45966996, 0.96409622 ],
        [ -2.27481527,  5.41552665,    -2.9855003,  -2.3050949,  0.9541035 ],
        [ -2.12902364,  4.98079254,    -2.97772943, -2.16679676, 0.94267606 ],
        [ -1.99510137,  4.58267559,    -2.96709304, -2.04139114, 0.92969842 ],
        [ -1.87013712,  4.21273373,    -2.95285423, -1.92610224, 0.9149965 ],
        [ -1.75144347,  3.86330586,    -2.93405429, -1.81846069, 0.89829516 ],
        [ -1.63666896,  3.5278949,     -2.90945043, -1.71642463, 0.87922099 ],
        [ -1.52265719,  3.19794423,    -2.87710678, -1.61741379, 0.85707375 ],
        [ -1.40441145,  2.86020889,    -2.83346983, -1.51760582, 0.83042664 ],
        [ -1.27454318,  2.49606791,    -2.77187453, -1.41190242, 0.79663364 ],
        [ -1.15783906,  2.1764507,     -2.70344383, -1.3209015,  0.76227052 ],
        [ -1.04989503,  1.88853752,    -2.62940478, -1.24047575, 0.72741128 ],
        [ -0.94788567,  1.6242454,     -2.55110531, -1.16805218, 0.69220861 ],
        [ -0.85100924,  1.38094967,    -2.4709114,  -1.10267575, 0.65728555 ],
        [ -0.75583037,  1.14967237,    -2.38850323, -1.0417842,  0.62214644 ],
        [ -0.6606744,   0.92638273,    -2.30452529, -0.98426559, 0.58680099 ],
        [ -0.56538356,  0.71078418,    -2.22073602, -0.93002273, 0.55177668 ],
        [ -0.46894693,  0.50063867,    -2.13792894, -0.87848562, 0.51724522 ],
        [ -0.37052079,  0.29423388,    -2.05684762, -0.82925362, 0.48340904 ],
        [ -0.26922237,  0.08991267,    -1.97804738, -0.78197219, 0.45043537 ],
        [ -0.16413266,  -0.1139088,    -1.90194577, -0.73633655, 0.41846772 ],
        [ -0.054509,    -0.31834833,   -1.82900029, -0.69217375, 0.38769107 ],
        [ 0.06046656,   -0.52457702,   -1.75952995, -0.64931888, 0.35824852 ],
        [ 0.18184761,   -0.73407809,   -1.69366519, -0.60756663, 0.33021201 ],
        [ 0.31061902,   -0.94809318,   -1.63155988, -0.56678778, 0.30366628 ],
        [ 0.44790054,   -1.1679879,    -1.57328453, -0.52685325, 0.27865976 ],
        [ 0.59490197,   -1.39516909,   -1.51886836, -0.4876522,  0.25522014 ],
        [ 0.75092634,   -1.62815515,   -1.46889949, -0.44955775, 0.23361331 ],
        [ 0.91295437,   -1.86244672,   -1.42420978, -0.413316,   0.21421169 ],
        [ 1.08238397,   -2.10026209,   -1.38406074, -0.37854189, 0.19670579 ],
        [ 1.26035692,   -2.34328759,   -1.34791333, -0.34497817, 0.18086815 ],
        [ 1.4482166,    -2.59336101,   -1.31528679, -0.31238468, 0.16649363 ],
        [ 1.64717352,   -2.85203589,   -1.28580978, -0.28059304, 0.15342303 ],
        [ 1.85838466,   -3.12072444,   -1.25917527, -0.24947724, 0.14152449 ],
        [ 2.08336066,   -3.40122415,   -1.23508204, -0.21889045, 0.13066789 ],
        [ 2.32339214,   -3.69499592,   -1.21329977, -0.18874674, 0.12075446 ],
        [ 2.58016892,   -4.0039423,    -1.19359658, -0.15893468, 0.11168451 ],
        [ 2.85527744,   -4.32979096,   -1.17578408, -0.12938179, 0.1033782 ],
        [ 3.15059686,   -4.6745767,    -1.15968168, -0.10000671, 0.09575934 ],
        [ 3.46799715,   -5.04028422,   -1.14513448, -0.07075197, 0.08876397 ],
        [ 3.80994126,   -5.42954181,   -1.13198398, -0.04152921, 0.08232638 ],
        [ 4.17868276,   -5.84469602,   -1.12010355, -0.01229407, 0.07639637 ],
        [ 4.57706811,   -6.28872605,   -1.10936515, 0.01702229,  0.07092308 ],
        [ 5.00813597,   -6.76478221,   -1.09965573, 0.04647574,  0.06586293 ],
        [ 5.4755179,    -7.27662574,   -1.09086735, 0.07613458,  0.06117462 ],
        [ 5.98283886,   -7.82796598,   -1.08290955, 0.10603769,  0.05682553 ],
        [ 6.53451816,   -8.42333411,   -1.07569399, 0.13624303,  0.05278351 ],
        [ 7.13537008,   -9.0676391,    -1.06914293, 0.16679822,  0.04902121 ],
        [ 7.79080481,   -9.7663812,    -1.06318628, 0.19774905,  0.04551444 ],
        [ 8.50682928,   -10.52564772,  -1.05776188, 0.22913618,  0.04224216 ],
        [ 9.28964796,   -11.3516896,   -1.05281748, 0.26097706,  0.03918773 ],
        [ 10.14666359,  -12.25197783,  -1.04830341, 0.2933065,   0.03633444 ],
        [ 11.08527773,  -13.23393648,  -1.04417965, 0.32612673,  0.03366986 ],
        [ 12.1136912,   -14.30578802,  -1.04040991, 0.35943996,  0.03118236 ],
        [ 13.24050442,  -15.47613513,  -1.03696295, 0.39323418,  0.02886196 ],
        [ 14.47551756,  -16.75479211,  -1.03380948, 0.42750705,  0.02669851 ],
        [ 15.82873069,  -18.15174622,  -1.03092482, 0.46223551,  0.02468359 ],
        [ 17.31134383,  -19.67819384,  -1.02828581, 0.49740384,  0.02280856 ],
        [ 18.88832591,  -21.29787507,  -1.02593576, 0.53199868,  0.021112 ],
        [ 20.62829918,  -23.08102816,  -1.02376103, 0.56730907,  0.01951828 ],
        [ 22.47666448,  -24.97146735,  -1.02182006, 0.60201764,  0.01807571 ],
        [ 24.57274049,  -27.11128586,  -1.01997232, 0.63840684,  0.01668379 ],
        [ 26.79231697,  -29.37331798,  -1.01833054, 0.67401031,  0.015431 ],
        [ 29.12104014,  -31.7429879,   -1.01687645, 0.70860383,  0.01430827 ],
        [ 31.75254986,  -34.41703052,  -1.01548928, 0.7447927,   0.0132251 ],
        [ 34.74059872,  -37.44932923,  -1.01416806, 0.78270772,  0.0121819 ],
        [ 37.87220238,  -40.62343037,  -1.01300624, 0.81936585,  0.01125481 ],
        [ 41.42775859,  -44.22320424,  -1.01189914, 0.85775203,  0.01036248 ],
        [ 45.12637165,  -47.96398797,  -1.01093179, 0.89458407,  0.00957532 ],
        [ 49.32169678,  -52.20318341,  -1.0100093,  0.93312234,  0.00881787 ],
        [ 53.64336503,  -56.56633145,  -1.00920917, 0.96976052,  0.00815526 ],
        [ 58.53542848,  -61.50152425,  -1.0084453,  1.00805366,  0.00751755 ],
        [ 63.5110249,   -66.5174629,   -1.00778849, 1.04404789,  0.00696504 ],
        [ 69.1242819,   -72.1726266,   -1.00716047, 1.08161094,  0.00643295 ],
        [ 75.4875309,   -78.57947626,  -1.00656093, 1.12087224,  0.00592139 ],
        [ 81.88429598,  -85.01653747,  -1.00605171, 1.15731476,  0.00548399 ],
        [ 89.10090155,  -92.27500893,  -1.00556459, 1.19533038,  0.00506296 ],
        [ 97.28202476,  -100.4996983,  -1.00509937, 1.23504959,  0.00465837 ],
        [ 105.37058788, -108.62789211, -1.00471013, 1.27131691,  0.0043179 ],
        [ 114.47879249, -117.77725333, -1.00433738, 1.30911264,  0.0039901 ],
        [ 124.78374849, -128.12501709, -1.003981,   1.34856201,  0.00367502 ],
        [ 136.5026042,  -139.88847317, -1.00364086, 1.38980643,  0.0033727 ],
        [ 147.87455236, -151.30019806, -1.00336213, 1.42670857,  0.00312375 ],
        [ 160.68879541, -164.1557666,  -1.00309516, 1.46516124,  0.00288423 ],
        [ 175.19786073, -178.70783464, -1.00283987, 1.50529288,  0.00265415 ],
        [ 188.80505695, -192.35225188, -1.00263598, 1.54012265,  0.00246964 ],
        [ 204.02114087, -207.60691643, -1.0024401,  1.57631322,  0.00229172 ]
    ],

    energy_table => [
        [ -4.47131967,  0.00783589, 0.00938925, 0.00783589, 0.00938925,    5.5495044,   -1.5982706 ],
        [ -4.34086863,  0.00916132, 0.01097231, 0.00916132, 0.01097231,    5.3404549,   -1.6068789 ],
        [ -4.18432739,  0.01104979, 0.01322452, 0.01104979, 0.01322452,    5.088092,    -1.6175956 ],
        [ -4.00830247,  0.01363952, 0.01630609, 0.01363952, 0.01630609,    4.8022706,   -1.631178 ],
        [ -3.90394164,  0.01545117, 0.01845662, 0.01545117, 0.01845662,    4.6315802,   -1.6395697 ],
        [ -3.80898673,  0.01730612, 0.02065379, 0.01730612, 0.02065379,    4.4755503,   -1.6479429 ],
        [ -3.75423088,  0.01847445, 0.02203509, 0.01847445, 0.02203509,    4.385166,    -1.6536803 ],
        [ -3.36690888,  0.02928493, 0.03470996, 0.02928493, 0.03470996,    3.736432,    -1.7017113 ],
        [ -3.07655735,  0.04125946, 0.04849012, 0.04125946, 0.04849012,    3.2365065,   -1.7470365 ],
        [ -2.82203786,  0.05555425, 0.06451712, 0.05555425, 0.06451712,    2.7862198,   -1.7978776 ],
        [ -2.61566768,  0.07048605, 0.08069164, 0.07048605, 0.08069164,    2.4103868,   -1.8502625 ],
        [ -2.43594758,  0.08644149, 0.09725027, 0.08644149, 0.09725027,    2.0733025,   -1.9052315 ],
        [ -2.27481527,  0.10343628, 0.11396718, 0.10343628, 0.11396718,    1.7620284,   -1.9628432 ],
        [ -2.12902364,  0.12123363, 0.13034523, 0.12123363, 0.13034523,    1.4717663,   -2.0187607 ],
        [ -1.99510137,  0.13973575, 0.14602297, 0.13973575, 0.14602297,    1.1979866,   -2.0740378 ],
        [ -1.87013712,  0.15890176, 0.16067254, 0.15890176, 0.16067254,    0.9353405,   -2.1340003 ],
        [ -1.75144347,  0.1787722,  0.17399263, 0.1787722,  0.17399263,    0.67841463,  -2.1857218 ],
        [ -1.63666896,  0.19942571, 0.18564192, 0.19942571, 0.18564192,    0.42520678,  -2.2281595 ],
        [ -1.52265719,  0.22116225, 0.19527618, 0.22116225, 0.19527618,    0.16867579,  -2.2570623 ],
        [ -1.40441145,  0.24470767, 0.20241964, 0.24470767, 0.20241964,    -0.09907346, -2.2551166 ],
        [ -1.27454318,  0.27128253, 0.20603498, 0.27128253, 0.20603498,    -0.39062516, -2.2088258 ],
        [ -1.15783906,  0.29530606, 0.20495654, 0.29530606, 0.20495654,    -0.64461257, -2.1252044 ],
        [ -1.04989503,  0.31720441, 0.20018285, 0.31720441, 0.20018285,    -0.86891062, -2.0215428 ],
        [ -0.94788567,  0.33725925, 0.19253617, 0.33725925, 0.19253617,    -1.069693,   -1.9225398 ],
        [ -0.85100924,  0.3554595,  0.18285397, 0.3554595,  0.18285397,    -1.2517342,  -1.8474715 ],
        [ -0.75583037,  0.37233771, 0.1715623,  0.37233771, 0.1715623,     -1.4246162,  -1.7994575 ],
        [ -0.6606744,   0.38807642, 0.15908726, 0.38807642, 0.15908726,    -1.5942344,  -1.778075 ],
        [ -0.56538356,  0.40261437, 0.14598514, 0.40261437, 0.14598514,    -1.7632437,  -1.7780326 ],
        [ -0.46894693,  0.41604799, 0.13264072, 0.41604799, 0.13264072,    -1.9351422,  -1.7922892 ],
        [ -0.37052079,  0.42844665, 0.11939271, 0.42844665, 0.11939271,    -2.1125334,  -1.8147549 ],
        [ -0.26922237,  0.43988091, 0.10651225, 0.43988091, 0.10651225,    -2.2976661,  -1.8408689 ],
        [ -0.16413266,  0.45041766, 0.09421084, 0.45041766, 0.09421084,    -2.4925694,  -1.8675344 ],
        [ -0.054509,    0.46010045, 0.08266919, 0.46010045, 0.08266919,    -2.6987702,  -1.8927883 ],
        [ 0.06046656,   0.46897852, 0.07200846, 0.46897852, 0.07200846,    -2.9178176,  -1.9155147 ],
        [ 0.18184761,   0.47711365, 0.06228636, 0.47711365, 0.06228636,    -3.151651,   -1.9352096 ],
        [ 0.31061902,   0.48455462, 0.05353421, 0.48455462, 0.05353421,    -3.4020486,  -1.9517323 ],
        [ 0.44790054,   0.49135247, 0.04574546, 0.49135247, 0.04574546,    -3.671042,   -1.9651827 ],
        [ 0.59490197,   0.49755603, 0.03888779, 0.49755603, 0.03888779,    -3.9608323,  -1.9758112 ],
        [ 0.75092634,   0.503146,   0.03297744, 0.503146,   0.03297744,    -4.2698469,  -1.9838926 ],
        [ 0.91295437,   0.50807386, 0.02802867, 0.50807386, 0.02802867,    -4.591855,   -1.9897486 ],
        [ 1.08238397,   0.51245753, 0.02387002, 0.51245753, 0.02387002,    -4.9294004,  -1.9938533 ],
        [ 1.26035692,   0.51638267, 0.02036905, 0.51638267, 0.02036905,    -5.2845529,  -1.9966357 ],
        [ 1.4482166,    0.51992124, 0.01741403, 0.51992124, 0.01741403,    -5.6598571,  -1.9984974 ],
        [ 1.64717352,   0.52312803, 0.01491615, 0.52312803, 0.01491615,    -6.057621,   -1.9996684 ],
        [ 1.85838466,   0.52604683, 0.01280255, 0.52604683, 0.01280255,    -6.4800679,  -2.000339 ],
        [ 2.08336066,   0.52871778, 0.01100988, 0.52871778, 0.01100988,    -6.9301495,  -2.0006711 ],
        [ 2.32339214,   0.53117088, 0.00948793, 0.53117088, 0.00948793,    -7.4103979,  -2.0007759 ],
        [ 2.58016892,   0.53343454, 0.0081928,  0.53343454, 0.0081928,     -7.9241516,  -2.0007662 ],
        [ 2.85527744,   0.53553084, 0.00708907, 0.53553084, 0.00708907,    -8.4745755,  -2.0006992 ],
        [ 3.15059686,   0.53747989, 0.0061464,  0.53747989, 0.0061464,     -9.0654041,  -2.0005755 ],
        [ 3.46799715,   0.53929791, 0.00533989, 0.53929791, 0.00533989,    -9.7003644,  -2.0004511 ],
        [ 3.80994126,   0.54100104, 0.00464774, 0.54100104, 0.00464774,    -10.384388,  -2.0003427 ],
        [ 4.17868276,   0.54260101, 0.00405266, 0.54260101, 0.00405266,    -11.121976,  -2.0002427 ],
        [ 4.57706811,   0.54410949, 0.00353952, 0.54410949, 0.00353952,    -11.918825,  -2.0001654 ],
        [ 5.00813597,   0.54553608, 0.00309591, 0.54553608, 0.00309591,    -12.781018,  -2.0001093 ],
        [ 5.4755179,    0.54688982, 0.00271118, 0.54688982, 0.00271118,    -13.715821,  -2.0000685 ],
        [ 5.98283886,   0.54817732, 0.00237688, 0.54817732, 0.00237688,    -14.730489,  -2.0000412 ],
        [ 6.53451816,   0.54940523, 0.00208555, 0.54940523, 0.00208555,    -15.833864,  -2.0000234 ],
        [ 7.13537008,   0.550579,   0.00183423, 0.550579,   0.00183423,    -17.035578,  -2.0000128 ],
        [ 7.79080481,   0.55170346, 0.00160805, 0.55170723, 0.0016382504,  -18.346453,  -2.0000066 ],
        [ 8.50682928,   0.55278226, 0.00141241, 0.55284899, 0.0015728222,  -19.778505,  -2.0000032 ],
        [ 9.28964796,   0.55381823, 0.00124064, 0.55406757, 0.0015201651,  -21.344144,  -2.0000015 ],
        [ 10.14666359,  0.55481439, 0.0010896,  0.55534851, 0.0014650338,  -23.058176,  -2.0000006 ],
        [ 11.08527773,  0.55577247, 0.00095674, 0.55668518, 0.0013893827,  -24.935405,  -2.0000002 ],
        [ 12.1136912,   0.55669409, 0.00083984, 0.55806868, 0.0013048369,  -26.992232,  -2.0000001 ],
        [ 13.24050442,  0.55758038, 0.00073702, 0.55948688, 0.0012121471,  -29.245858,  -2 ],
        [ 14.47551756,  0.55843271, 0.00064657, 0.5609249,  0.0011181483,  -31.715885,  -2 ],
        [ 15.82873069,  0.55925188, 0.00056706, 0.56237029, 0.0010166663,  -34.422311,  -2 ],
        [ 17.31134383,  0.5600389,  0.00049718, 0.56394009, 0.0010165902,  -37.387537,  -2 ],
        [ 18.88832591,  0.56077413, 0.00043741, 0.56544239, 0.00089316621, -40.541501,  -2 ],
        [ 20.62829918,  0.5614873,  0.00038427, 0.56689775, 0.00078369272, -44.021448,  -2 ],
        [ 22.47666448,  0.56215396, 0.0003387,  0.56825666, 0.00069002374, -47.718178,  -2 ],
        [ 24.57274049,  0.56281864, 0.00029705, 0.56961013, 0.0006045601,  -51.91033,   -2 ],
        [ 26.79231697,  0.56343712, 0.00026153, 0.57086836, 0.00053179956, -56.349483,  -2 ],
        [ 29.12104014,  0.56400972, 0.00023129, 0.5720323,  0.00046997079, -61.00693,   -2 ],
        [ 31.75254986,  0.56458061, 0.00020359, 0.57319191, 0.00041337974, -66.269949,  -2 ],
        [ 34.74059872,  0.56514972, 0.00017828, 0.5743471,  0.00036175195, -72.246047,  -2 ],
        [ 37.87220238,  0.5656734,  0.00015694, 0.57540945, 0.0003182701,  -78.509254,  -2 ],
        [ 41.42775859,  0.56619545, 0.00013744, 0.57646791, 0.000278584,   -85.620367,  -2 ],
        [ 45.12637165,  0.56667251, 0.00012111, 0.57743469, 0.00024536769, -93.017593,  -2 ],
        [ 49.32169678,  0.56714812, 0.00010618, 0.57839809, 0.00021501917, -101.40824,  -2 ],
        [ 53.64336503,  0.56757918, 9.375e-05,  0.57927092, 0.0001897927,  -110.05158,  -2 ],
        [ 58.53542848,  0.56800897, 8.237e-05,  0.58014089, 0.00016670294, -119.83571,  -2 ],
        [ 63.5110249,   0.56839467, 7.299e-05,  0.58092136, 0.00014765938, -129.7869,   -2 ],
        [ 69.1242819,   0.56877929, 6.436e-05,  0.58169944, 0.00013018156, -141.01341,  -2 ],
        [ 75.4875309,   0.5691628,  5.647e-05,  0.58247509, 0.00011419256, -153.73991,  -2 ],
        [ 81.88429598,  0.56950276, 5.004e-05,  0.5831625,  0.00010116962, -166.53344,  -2 ],
        [ 89.10090155,  0.56984182, 4.414e-05,  0.58384794, 8.9211956e-05, -180.96665,  -2 ],
        [ 97.28202476,  0.57017995, 3.873e-05,  0.5845314,  7.8268295e-05, -197.3289,   -2 ],
        [ 105.37058788, 0.57047505, 3.439e-05,  0.58512778, 6.9484966e-05, -213.50603,  -2 ],
        [ 114.47879249, 0.57076942, 3.039e-05,  0.58572261, 6.1406075e-05, -231.72243,  -2 ],
        [ 124.78374849, 0.57106305, 2.673e-05,  0.58631587, 5.3998895e-05, -252.33235,  -2 ],
        [ 136.5026042,  0.57135593, 2.338e-05,  0.58690754, 4.7230982e-05, -275.77006,  -2 ],
        [ 147.87455236, 0.57160637, 2.075e-05,  0.58741342, 4.1914345e-05, -298.51395,  -2 ],
        [ 160.68879541, 0.57185624, 1.833e-05,  0.58791811, 3.7024044e-05, -324.14244,  -2 ],
        [ 175.19786073, 0.57210554, 1.611e-05,  0.5884216,  3.2540436e-05, -353.16057,  -2 ],
        [ 188.80505695, 0.57231285, 1.441e-05,  0.58884025, 2.9100613e-05, -380.37496,  -2 ],
        [ 204.02114087, 0.57251975, 1.284e-05,  0.58925806, 2.5918645e-05, -410.80713,  -2 ]
    ],

    impulse_table => [
        [
            -5.0762094, 13.814839,   -5.0766569,    0.026109094, -0.00097428934, -0.15334894,
            -1.7072497, -0.23126674, -0.0015786156, -0.37471572, 0.28281827,     0.99367503,
            -2.4843639
        ],
        [
            -4.9618449, 13.471746,   -4.9623761,    0.027420191, -0.0010923349, -0.15259247,
            -1.7064932, -0.28845017, -0.0017698818, -0.37395925, 0.28281819,    0.99274508,
            -2.443936
        ],
        [
            -4.6051699, 12.401725,   -4.6060772,    0.031827131, -0.0015604785, -0.14959247,
            -1.7034932, -0.46679241, -0.0025284025, -0.37095925, 0.2828176,     0.98887324,
            -2.3206454
        ],
        [
            -4.1997048, 11.18534,    -4.2013722,    0.037385765, -0.0023407177, -0.14459246,
            -1.6984932, -0.66954089, -0.0037926038, -0.36595925, 0.28281543,    0.98191556,
            -2.1868451
        ],
        [
            -3.9120227, 10.322315,   -3.9145909,    0.041617881, -0.003120957, -0.13959246,
            -1.6934932, -0.81341253, -0.0050568051, -0.36095925, 0.28281132,   0.97448985,
            -2.0972174
        ],
        [
            -3.5065576, 9.1060051,  -3.5112794,    0.047786615, -0.0046814354, -0.12959247,
            -1.6834932, -1.0162705, -0.0075852075, -0.35095925, 0.28279466,    0.95863751,
            -1.9808298
        ],
        [
            -2.995732,  7.5739703,  -3.0059193,   0.055195151, -0.0078023906, -0.10959259,
            -1.6634933, -1.2723309, -0.012642011, -0.33095924, 0.28270872,    0.92436975,
            -1.8574034
        ],
        [
            -2.3025848, 5.4984486,  -2.3316165,  0.061499837, -0.015604465, -0.059604632,
            -1.6134995, -1.6246998, -0.02528364, -0.28095923, 0.28194585,   0.8326692,
            -1.7581799
        ],
        [
            -1.8971197, 4.292457,   -1.9508375,   0.060733006, -0.023400494, -0.0097621409,
            -1.563581,  -1.8432601, -0.037917989, -0.23110789, 0.27991529,   0.74180385,
            -1.7622381
        ],
        [
            -1.6094377, 3.44876,    -1.6925497,   0.057592895, -0.03115216, 0.039316549,
            -1.5140613, -2.0182247, -0.050498696, -0.18170302, 0.27611895,  0.65786681,
            -1.8094624
        ],
        [
            -1.3862941, 2.8089431,  -1.5026007,   0.054335357, -0.0387137, 0.085907996,
            -1.4658892, -2.1816355, -0.062845347, -0.13409099, 0.27026779, 0.58358505,
            -1.8799587
        ],
        [
            -1.2039725, 2.3018358,  -1.3563933,   0.052153586,  -0.045724783, 0.12694902,
            -1.4210538, -2.3503659, -0.074478409, -0.091022673, 0.2623394,    0.51967491,
            -1.9635698
        ],
        [
            -1.0498219, 1.8883451,  -1.2404225,   0.051369916,  -0.051637227, 0.15956142,
            -1.382356,  -2.5336374, -0.084538136, -0.055965432, 0.25256577,   0.46564478,
            -2.0539407
        ],
        [
            -0.91629047, 1.5440464,  -1.1463593,   0.051687187,  -0.05605233, 0.18334191,
            -1.3517295,  -2.7332664, -0.092168531, -0.031642045, 0.24135921,  0.42038867,
            -2.1468271
        ],
        [
            -0.79850744, 1.2524038,  -1.0686735,   0.052619184, -0.059030169, 0.20010617,
            -1.3287492,  -2.9440554, -0.097221066, -0.01729933, 0.22921505,   0.38259567,
            -2.2394252
        ],
        [
            -0.69314692, 1.0016831,  -1.0035162,  0.053781016,  -0.060946187, 0.21207323,
            -1.3115389,  -3.1578368, -0.10028024, -0.010218137, 0.21662726,   0.35098983,
            -2.3299704
        ],
        [
            -0.59783674, 0.78331463, -0.94812181, 0.054950483,   -0.062183932, 0.22092279,
            -1.2982839,  -3.3678442, -0.10207189, -0.0075924871, 0.20403327,   0.32444275,
            -2.417432
        ],
        [
            -0.51082536, 0.59091592, -0.90045729, 0.056025343,   -0.063009095, 0.2277324,
            -1.2877218,  -3.5700644, -0.10312051, -0.0073065108, 0.1917871,    0.30200752,
            -2.5012759
        ],
        [
            -0.43078266, 0.41965714, -0.85900012, 0.056970878,   -0.063581557, 0.23316207,
            -1.2790511,  -3.7626546, -0.10373877, -0.0084446297, 0.18015281,   0.28291244,
            -2.5812896
        ],
        [
            -0.35667468, 0.26583131, -0.82259236, 0.057784928, -0.063994414, 0.2376195,
            -1.2717637,  -3.9450763, -0.10410312, -0.01045509, 0.16930994,   0.26653818,
            -2.6574605
        ],
        [
            -0.28768181, 0.12655509, -0.79034094, 0.058478838,  -0.06430258, 0.24136435,
            -1.2655272,  -4.1174825, -0.10431349, -0.012953209, 0.15936472,  0.25239111,
            -2.7298937
        ],
        [
            -0.22314329, -0.00044382611, -0.76154796, 0.059068238,  -0.064539488, 0.24456859,
            -1.2601139,  -4.2803604,     -0.10442794, -0.015379542, 0.15036346,   0.24007834,
            -2.7987602
        ],
        [
            -0.16251867, -0.11697762, -0.73566153, 0.059568881,  -0.064726221, 0.24735057,
            -1.2553627,  -4.434337,   -0.10448122, -0.018353456, 0.14230585,   0.22928654,
            -2.8642626
        ],
        [
            -0.10536025, -0.22450799, -0.71224051, 0.059995024,  -0.064876568, 0.24979475,
            -1.251156,   -4.5800779,  -0.10449483, -0.021058611, 0.135155,     0.21976481,
            -2.9266152
        ],
        [
            -0.051293034, -0.32422705, -0.69092834, 0.060358921,  -0.064999821, 0.2519631,
            -1.2474034,   -4.7182355,  -0.10448244, -0.023778435, 0.12882507,   0.21131114,
            -2.9860315
        ],
        [
            0.0064942487, -0.42876194, -0.66901485, 0.060708162,  -0.065114574, 0.25414052,
            -1.2436223,   -4.866036,   -0.10444821, -0.026980441, 0.12252026,   0.20283753,
            -3.0499129
        ],
        [
            0.18232182, -0.73488118, -0.60741007, 0.061548769, -0.065382618, 0.25999483,
            -1.2334504, -5.3151042,  -0.10425768, -0.03785593, 0.1060238,    0.18027102,
            -3.2458433
        ],
        [
            0.26236453, -0.86882745, -0.58167035, 0.061837627,  -0.065474774, 0.26232683,
            -1.229429,  -5.5186118,  -0.10414563, -0.042972955, 0.099695888,  0.17141622,
            -3.3354436
        ],
        [
            0.40546537, -1.1008632, -0.53883386, 0.062239091,  -0.065607338, 0.26605529,
            -1.2230951, -5.8802056, -0.10392859, -0.052570089, 0.089914963,  0.1574587,
            -3.4957035
        ],
        [
            0.53062851, -1.2968199, -0.50436938, 0.062494013,  -0.065698125, 0.26890921,
            -1.2183758, -6.1936888, -0.10373375, -0.061045698, 0.082718615,  0.14694761,
            -3.6355893
        ],
        [
            0.69314744, -1.542779,  -0.46327427, 0.062725813,  -0.065791255, 0.27212738,
            -1.2132628, -6.5964538, -0.10348716, -0.073055903, 0.074899284,  0.1352579,
            -3.8164022
        ],
        [
            0.99325203, -1.9760035, -0.39646531, 0.062960774,  -0.065915116, 0.2768778,
            -1.206331,  -7.3270986, -0.10307732, -0.095300109, 0.063858905,  0.1181941,
            -4.1469368
        ],
        [
            1.0986125,  -2.1226947, -0.37536218, 0.063005132, -0.065948516, 0.27824105,
            -1.2045307, -7.579678,  -0.10295065, -0.10304942, 0.06077455,   0.11329233,
            -4.2618321
        ],
        [
            1.2089606,  -2.2737584, -0.35438454, 0.063037807, -0.065979424, 0.27952575,
            -1.2029307, -7.8421293, -0.10282815, -0.11180172, 0.057888225,  0.10864533,
            -4.381517
        ],
        [
            1.3862946,  -2.5116034, -0.32283231, 0.063069009, -0.066022003, 0.28132003,
            -1.2008811, -8.259659,  -0.10265182, -0.12551113, 0.053867847,  0.10206824,
            -4.5724975
        ],
        [
            1.6094382, -2.8034173, -0.28642616, 0.063084077, -0.066065692, 0.28317552,
            -1.199025, -8.7781749, -0.10246389, -0.14288395, 0.049677712,  0.095072011,
            -4.8105709
        ],
        [
            1.7917597,  -3.03657,   -0.25902372, 0.063084697, -0.066094906, 0.28441575,
            -1.1979597, -9.1966822, -0.1023354,  -0.1583744,  0.046817635,  0.090205335,
            -5.0033965
        ],
        [
            1.9459104,  -3.2304995, -0.23728567, 0.063080658, -0.066115926, 0.28530339,
            -1.1972931, -9.5472592, -0.10224237, -0.16998764, 0.044715359,  0.086576832,
            -5.1653558
        ],
        [
            2.0794418,   -3.3963833, -0.21940286, 0.063075335, -0.06613192, 0.28597058, -1.1968525, -9.8487287,
            -0.10217178, -0.1808336, 0.043089731, 0.083739291, -5.3049318
        ],
        [
            2.3027871,  -3.6699781, -0.19124301, 0.063064485, -0.066155653, 0.28690745,
            -1.1964264, -10.348801, -0.1020722,  -0.19855661, 0.04070772,   0.078646232,
            -5.5465775
        ],
        [
            2.995761,   -4.4943995, -0.11512711, 0.06303157,  -0.066201054, 0.28878216,
            -1.1953855, -11.87319,  -0.10186894, -0.25465258, 0.035226308,  0.069234534,
            -6.2532092
        ],
        [
            3.912114,   -5.5450193, -0.03320703, 0.063005154, -0.066233587, 0.28992375,
            -1.1957032, -13.844148, -0.10174797, -0.32908136, 0.030583536,  0.060730156,
            -7.179522
        ],
        [
            4.6054898,   -6.3202462,  0.019032947, 0.062995424, -0.06624275, 0.29030138, -1.1956463, -15.31338,
            -0.10170632, -0.38545969, 0.028137648, 0.056068518, -7.8767626
        ],
        [
            5.298459,   -7.0832032, 0.065155067, 0.062990459, -0.066247657, 0.29049002,
            -1.1956519, -16.76893,  -0.10168568, -0.44182019, 0.026235274,  0.052371213,
            -8.5718627
        ],
        [
            6.2146772,  -8.0786514, 0.11900399,  0.062987488, -0.066250594, 0.29060329,
            -1.1956567, -18.679466, -0.10167329, -0.51634076, 0.024259102,  0.048479986,
            -9.4894809
        ],
        [
            6.9079534,  -8.8242354, 0.15549757,  0.062986514, -0.066251573, 0.29064107,
            -1.1956585, -20.117346, -0.10166915, -0.57272692, 0.023043858,  0.046069036,
            -10.183264
        ],
        [
            7.6009964,  -9.5644278, 0.18901977,  0.062986037, -0.066252062, 0.29065995,
            -1.1956594, -21.549634, -0.10166709, -0.62909255, 0.022002727,  0.043996283,
            -10.876575
        ],
        [
            8.5172296, -10.536648, 0.22957528,  0.062985757, -0.066252364, 0.29067129,
            -1.19566,  -23.436964, -0.10166585, -0.7036086,  0.020827736,  0.041651875,
            -11.792979
        ],
        [
            9.2104463,  -11.268287, 0.25786199,  0.062985666, -0.066252455, 0.29067507,
            -1.1956602, -24.861079, -0.10166544, -0.75998615, 0.020057401,  0.040113026,
            -12.486256
        ],
        [
            10.819874,  -12.956664, 0.31709712,  0.062985596, -0.066252535, 0.29067809,
            -1.1956603, -28.1575,   -0.10148003, -0.87288188, 0.018560361,  0.037120375,
            -14.095735
        ],
        [
            11.513084,  -13.680285, 0.34029672,  0.062985587, -0.066252546, 0.29067847,
            -1.1956603, -29.573854, -0.10080138, -0.90417025, 0.018012993,  0.036025815,
            -14.788951
        ],
        [
            13.122503,  -15.353753, 0.38981517,   0.062985581, -0.066252572, 0.29067877,
            -1.1956604, -32.855838, -0.098226639, -0.97335261, 0.016910724,  0.033821414,
            -16.398377
        ],
        [
            13.815511,  -16.071947, 0.40952421,   0.06298558, -0.066252555, 0.29067881,
            -1.1956604, -34.26672,  -0.096914179, -1.0018354, 0.016495386,  0.032990756,
            -17.091385
        ],
        [
            16.118133,  -18.450018, 0.46932228,   0.06298558, -0.06625255, 0.29067881,
            -1.1956604, -38.946569, -0.092359069, -1.0917541, 0.01530956,  0.030619118,
            -19.394008
        ],
        [
            18.420753, -20.818024, 0.52201732,  0.062985579, -0.066252577, 0.29067821,
            -1.195661, -43.616579, -0.08796574, -1.1756137,  0.014349371,  0.028698742,
            -21.696628
        ]
    ],

    tail_shock_table => [
        [ 7.4499071, -0.70345,    -0.041632777, -0.041632777 ],
        [ 7.6003922, -0.7117947,  -0.049585103, -0.033094862 ],
        [ 8.5169779, -0.76079838, -0.058307882, -0.019504097 ],
        [ 9.2103169, -0.79603669, -0.060221734, -0.014471729 ],
        [ 10.819846, -0.8728798,  -0.060875909, -0.0078980991 ],
        [ 11.513069, -0.90416921, -0.060468774, -0.0061826255 ],
        [ 13.1225,   -0.9733524,  -0.058924208, -0.0035136695 ],
        [ 13.815509, -1.0018353,  -0.058136885, -0.0027328713 ],
        [ 16.118133, -1.0917984,  -0.055408089, -0.0010418925 ]
    ],

    blast_info => [
        0.15959245, 1.7134927,  0.38092943, -0.25284024, 2.0707224,  0.0017293655, 3.6934624,   -0.15604781,
        0.49420116, 0.38980434, 0.65561346, -0.1756547,  0.30613069, 0.28281849,   0.037783791, -0.039743581,
        1719,       -0.70345,   29264.445,  -0.84730506
    ],

};
1;
