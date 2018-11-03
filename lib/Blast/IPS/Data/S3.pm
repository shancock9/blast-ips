package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=3
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S3'} = {

    table_name => 'S3',
    symmetry   => 2,
    gamma      => 3,

    shock_table_info => [ 2, 3, 1.44e-06, 8000, 2.3e-07, 10, 6.36e-08, 8.8e-07, 5e-07, 289.43, -0.80301 ],

    shock_table => [
        [ -4.16830247,  12.00017702,   -2.99997236, -4.16951744, 0.99817646 ],
        [ -3.53922615,  10.11299026,   -2.99985977, -3.54235061, 0.99530621 ],
        [ -3.20990726,  9.12511712,    -2.99960569, -3.21503262, 0.99229327 ],
        [ -2.91037377,  8.22671013,    -2.99902372, -2.91841723, 0.98789017 ],
        [ -2.66570809,  7.49306634,    -2.99796156, -2.6773367,  0.98246735 ],
        [ -2.41298021,  6.73565222,    -2.99567051, -2.43000828, 0.97427697 ],
        [ -2.25152444,  6.25218504,    -2.99297595, -2.27325971, 0.9671185 ],
        [ -2.08132223,  5.74313484,    -2.98834991, -2.10944644, 0.95738364 ],
        [ -1.92903604,  5.28851771,    -2.98172485, -1.96446587, 0.9462379 ],
        [ -1.7905228,   4.8761019,     -2.97256371, -1.83424481, 0.93358738 ],
        [ -1.66175139,  4.49407299,    -2.96013271, -1.71492203, 0.91919919 ],
        [ -1.54113309,  4.13796331,    -2.94370731, -1.60499893, 0.90298259 ],
        [ -1.42557475,  3.79897641,    -2.92214824, -1.50168744, 0.88455289 ],
        [ -1.31283868,  3.47105646,    -2.89402558, -1.40312409, 0.86347707 ],
        [ -1.19938183,  3.14472174,    -2.85694198, -1.3065196,  0.8388605 ],
        [ -1.07833671,  2.8018683,     -2.80581161, -1.20676951, 0.80858385 ],
        [ -0.95594788,  2.46232193,    -2.74046994, -1.10989917, 0.77370527 ],
        [ -0.84462665,  2.16113718,    -2.66868789, -1.02570133, 0.73847054 ],
        [ -0.7408257,   1.88803782,    -2.59179654, -0.95086769, 0.70302645 ],
        [ -0.6426464,   1.6374611,     -2.51161681, -0.88356564, 0.66774203 ],
        [ -0.54847531,  1.40476576,    -2.42970142, -0.82231982, 0.63287738 ],
        [ -0.45506202,  1.1817132,     -2.34560409, -0.76483321, 0.59791109 ],
        [ -0.36133883,  0.96586884,    -2.26042244, -0.71043299, 0.56304191 ],
        [ -0.26698287,  0.75659444,    -2.17575236, -0.65893464, 0.52870409 ],
        [ -0.17107512,  0.55194569,    -2.09244293, -0.60985265, 0.49506811 ],
        [ -0.07276095,  0.35026054,    -2.01121665, -0.56280675, 0.46229837 ],
        [ 0.02888957,   0.14986436,    -1.93258981, -0.51744816, 0.43051573 ],
        [ 0.13470608,   -0.05058085,   -1.85704776, -0.47353658, 0.39986161 ],
        [ 0.24556147,   -0.2523823,    -1.78495153, -0.43086527, 0.37045262 ],
        [ 0.36246351,   -0.45697263,   -1.71651467, -0.38922897, 0.34236405 ],
        [ 0.48630801,   -0.6654761,    -1.65197486, -0.34851208, 0.31569356 ],
        [ 0.61830477,   -0.87944408,   -1.59137459, -0.30854098, 0.29046604 ],
        [ 0.7596604,    -1.10029854,   -1.53476566, -0.26919802, 0.26671365 ],
        [ 0.91015499,   -1.3272533,    -1.48264801, -0.2307573,  0.2446608 ],
        [ 1.06614758,   -1.55481885,   -1.43617185, -0.19417636, 0.22481915 ],
        [ 1.22900388,   -1.78523075,   -1.39453213, -0.15905951, 0.20687463 ],
        [ 1.39984334,   -2.02019433,   -1.35713782, -0.12514141, 0.19059751 ],
        [ 1.58001134,   -2.26159283,   -1.32345919, -0.09216964, 0.17577899 ],
        [ 1.7706821,    -2.51096515,   -1.29309492, -0.05997442, 0.16226213 ],
        [ 1.9731201,    -2.76988765,   -1.2656952,  -0.02840806, 0.14990976 ],
        [ 2.18872927,   -3.0400417,    -1.24095032, 0.00266378,  0.13860012 ],
        [ 2.41885565,   -3.32297399,   -1.21860581, 0.03333627,  0.12823457 ],
        [ 2.66497462,   -3.62034278,   -1.19843059, 0.06369687,  0.11872377 ],
        [ 2.92877813,   -3.9340207,    -1.18021025, 0.09383428,  0.10998475 ],
        [ 3.21214664,   -4.26605549,   -1.16375123, 0.12383168,  0.101943 ],
        [ 3.51695251,   -4.61844225,   -1.14888857, 0.15374652,  0.09453664 ],
        [ 3.84566235,   -4.99382003,   -1.13545448, 0.18366898,  0.08770101 ],
        [ 4.2005339,    -5.39454262,   -1.12331599, 0.21364222,  0.08138773 ],
        [ 4.58441845,   -5.82359631,   -1.11233944, 0.24373536,  0.07554682 ],
        [ 5.0003596,    -6.28413638,   -1.10240767, 0.27400363,  0.07013558 ],
        [ 5.45179388,   -6.77970946,   -1.0934133,  0.3045024,   0.06511513 ],
        [ 5.94235081,   -7.31403113,   -1.08526232, 0.33527149,  0.06045252 ],
        [ 6.47645386,   -7.89163735,   -1.07786385, 0.36637062,  0.05611493 ],
        [ 7.05892076,   -8.51743936,   -1.07113879, 0.39784806,  0.05207455 ],
        [ 7.69476446,   -9.19651187,   -1.06501991, 0.42973055,  0.04830883 ],
        [ 8.38979402,   -9.93473595,   -1.0594446,  0.46205362,  0.0447962 ],
        [ 9.15001541,   -10.73815725,  -1.05436018, 0.49483,     0.04151938 ],
        [ 9.98243239,   -11.61383374,  -1.0497173,  0.52808466,  0.0384613 ],
        [ 10.89424716,  -12.56898896,  -1.04547494, 0.56181869,  0.03560837 ],
        [ 11.8934609,   -13.61164548,  -1.04159594, 0.59603399,  0.03294778 ],
        [ 12.98847423,  -14.75020606,  -1.03804829, 0.63071806,  0.03046849 ],
        [ 14.1884874,   -15.99387191,  -1.03480298, 0.6658578,   0.02815995 ],
        [ 15.50350055,  -17.35264054,  -1.03183385, 0.70143773,  0.02601219 ],
        [ 16.9443137,   -18.8373046,   -1.02911739, 0.73743875,  0.02401583 ],
        [ 18.45007586,  -20.38506711,  -1.02673698, 0.77222955,  0.02224066 ],
        [ 20.18212331,  -22.16137713,  -1.02444223, 0.8092067,   0.02050528 ],
        [ 22.05632948,  -24.07939121,  -1.02236778, 0.84610676,  0.01891536 ],
        [ 24.07542792,  -26.14170748,  -1.02049612, 0.88279224,  0.01746277 ],
        [ 26.23948437,  -28.3482429,   -1.01881098, 0.91911239,  0.01613954 ],
        [ 28.6989312,   -30.85192037,  -1.01720499, 0.95719449,  0.01486425 ],
        [ 31.33289533,  -33.52924973,  -1.01576488, 0.99478217,  0.01370837 ],
        [ 34.13469912,  -36.37336791,  -1.01447698, 1.03168828,  0.01266429 ],
        [ 37.31698478,  -39.59971084,  -1.01324855, 1.07034504,  0.01165892 ],
        [ 40.69282789,  -43.01838019,  -1.01215513, 1.10813778,  0.01075588 ],
        [ 44.24158952,  -46.60850878,  -1.01118529, 1.14483669,  0.00994815 ],
        [ 48.26046879,  -50.67042484,  -1.01025879, 1.18321108,  0.00917035 ],
        [ 52.46128663,  -54.9125742,   -1.0094417,  1.22024378,  0.00847916 ],
        [ 57.21818094,  -59.71247107,  -1.00866091, 1.25895217,  0.00781387 ],
        [ 62.15340664,  -64.68870575,  -1.00797678, 1.29602884,  0.00722694 ],
        [ 67.73661891,  -70.31457579,  -1.00732264, 1.33475809,  0.00666206 ],
        [ 74.08516567,  -76.70756971,  -1.00669825, 1.37528041,  0.00611937 ],
        [ 80.64212037,  -83.30661749,  -1.00615623, 1.4138086,   0.00564538 ],
        [ 88.08544877,  -90.79378519,  -1.00563843, 1.45408593,  0.00518993 ],
        [ 95.67848133,  -98.42789288,  -1.00519296, 1.49195774,  0.00479596 ],
        [ 104.27165039, -107.06380276, -1.00476684, 1.53149751,  0.00441714 ],
        [ 114.04710498, -116.88380722, -1.00435992, 1.57284855,  0.00405355 ],
        [ 123.91000396, -126.78795594, -1.00401423, 1.61125667,  0.00374318 ],
        [ 135.07945489, -138.00034383, -1.00368354, 1.65135369,  0.00344492 ],
        [ 146.11157231, -149.07152985, -1.0034064,  1.68793866,  0.00319391 ],
        [ 158.52261288, -161.52315325, -1.0031406,  1.72603907,  0.00295224 ],
        [ 172.55075943, -175.5935211,  -1.00288609, 1.76577913,  0.00271992 ],
        [ 188.48779585, -191.57455755, -1.00264279, 1.80729907,  0.00249698 ],
        [ 203.93845363, -207.06446545, -1.00244314, 1.84442428,  0.00231337 ]
    ],

    energy_table => [
        [ -4.16830247,  0.00031029, 0.00061283, 0.00031029, 0.00061283,    5.2652537,    -1.5830404 ],
        [ -3.53922615,  0.001068,   0.00208485, 0.001068,   0.00208485,    4.2605899,    -1.6176089 ],
        [ -3.20990726,  0.00202504, 0.00391302, 0.00202504, 0.00391302,    3.7243358,    -1.6436441 ],
        [ -2.91037377,  0.00359944, 0.00686351, 0.00359944, 0.00686351,    3.227849,     -1.6733097 ],
        [ -2.66570809,  0.00572024, 0.01074549, 0.00572024, 0.01074549,    2.8152947,    -1.7027366 ],
        [ -2.41298021,  0.00915285, 0.01683572, 0.00915285, 0.01683572,    2.3806477,    -1.7386678 ],
        [ -2.25152444,  0.01228717, 0.02220923, 0.01228717, 0.02220923,    2.0979859,    -1.7634741 ],
        [ -2.08132223,  0.01665689, 0.02942365, 0.01665689, 0.02942365,    1.7955481,    -1.7905043 ],
        [ -1.92903604,  0.02172653, 0.03740836, 0.02172653, 0.03740836,    1.5210297,    -1.816107 ],
        [ -1.7905228,   0.02748831, 0.04599648, 0.02748831, 0.04599648,    1.2677792,    -1.838922 ],
        [ -1.66175139,  0.03398488, 0.05506904, 0.03398488, 0.05506904,    1.0297129,    -1.8579881 ],
        [ -1.54113309,  0.04118133, 0.06436693, 0.04118133, 0.06436693,    0.80456192,   -1.8735356 ],
        [ -1.42557475,  0.04915853, 0.07374605, 0.04915853, 0.07374605,    0.58729586,   -1.8847471 ],
        [ -1.31283868,  0.0579943,  0.08298433, 0.0579943,  0.08298433,    0.37431037,   -1.8911661 ],
        [ -1.19938183,  0.06792135, 0.09189986, 0.06792135, 0.09189986,    0.15952469,   -1.8922138 ],
        [ -1.07833671,  0.07957158, 0.10035367, 0.07957158, 0.10035367,    -0.069403538, -1.8869575 ],
        [ -0.95594788,  0.09228694, 0.10706497, 0.09228694, 0.10706497,    -0.29981319,  -1.8756996 ],
        [ -0.84462665,  0.10444598, 0.1110041,  0.10444598, 0.1110041,     -0.50791931,  -1.8622184 ],
        [ -0.7408257,   0.11606518, 0.11250743, 0.11606518, 0.11250743,    -0.70052233,  -1.8495331 ],
        [ -0.6426464,   0.12709803, 0.11191877, 0.12709803, 0.11191877,    -0.88155378,  -1.8402347 ],
        [ -0.54847531,  0.13754127, 0.10960253, 0.13754127, 0.10960253,    -1.0545209,   -1.8358963 ],
        [ -0.45506202,  0.14761328, 0.10581655, 0.14761328, 0.10581655,    -1.2259405,   -1.8370657 ],
        [ -0.36133883,  0.15730465, 0.10081879, 0.15730465, 0.10081879,    -1.3983041,   -1.8436171 ],
        [ -0.26698287,  0.16654435, 0.09491134, 0.16654435, 0.09491134,    -1.5726922,   -1.854663 ],
        [ -0.17107512,  0.17533583, 0.08835614, 0.17533583, 0.08835614,    -1.7511994,   -1.868962 ],
        [ -0.07276095,  0.18368105, 0.08139464, 0.18368105, 0.08139464,    -1.9357232,   -1.8851937 ],
        [ 0.02888957,   0.19158963, 0.07423675, 0.19158963, 0.07423675,    -2.128228,    -1.9021598 ],
        [ 0.13470608,   0.19906263, 0.06707276, 0.19906263, 0.06707276,    -2.3304298,   -1.9188409 ],
        [ 0.24556147,   0.20610417, 0.06006267, 0.20610417, 0.06006267,    -2.5440691,   -1.9344964 ],
        [ 0.36246351,   0.21272524, 0.05333133, 0.21272524, 0.05333133,    -2.7711165,   -1.9486456 ],
        [ 0.48630801,   0.21892868, 0.04698421, 0.21892868, 0.04698421,    -3.0132909,   -1.96096 ],
        [ 0.61830477,   0.22473164, 0.04108604, 0.22473164, 0.04108604,    -3.2729034,   -1.9713451 ],
        [ 0.7596604,    0.23014695, 0.03568219, 0.23014695, 0.03568219,    -3.5522507,   -1.9798472 ],
        [ 0.91015499,   0.23514169, 0.0308405,  0.23514169, 0.0308405,     -3.8507911,   -1.9865415 ],
        [ 1.06614758,   0.23961665, 0.02666378, 0.23961665, 0.02666378,    -4.1611317,   -1.991447 ],
        [ 1.22900388,   0.24365603, 0.02305905, 0.24365603, 0.02305905,    -4.4857896,   -1.9949453 ],
        [ 1.39984334,   0.24732113, 0.01995057, 0.24732113, 0.01995057,    -4.8268598,   -1.9973991 ],
        [ 1.58001134,   0.25066581, 0.01726863, 0.25066581, 0.01726863,    -5.1869101,   -1.9990409 ],
        [ 1.7706821,    0.25373038, 0.01495618, 0.25373038, 0.01495618,    -5.5681951,   -2.0000762 ],
        [ 1.9731201,    0.25654926, 0.01296282, 0.25654926, 0.01296282,    -5.9731666,   -2.0006656 ],
        [ 2.18872927,   0.25915233, 0.01124421, 0.25915233, 0.01124421,    -6.404573,    -2.000941 ],
        [ 2.41885565,   0.26156338, 0.00976298, 0.26156338, 0.00976298,    -6.8650592,   -2.0010085 ],
        [ 2.66497462,   0.26380345, 0.00848626, 0.26380345, 0.00848626,    -7.3575438,   -2.0009821 ],
        [ 2.92877813,   0.26589164, 0.00738507, 0.26589164, 0.00738507,    -7.8854042,   -2.0008783 ],
        [ 3.21214664,   0.26784472, 0.00643437, 0.26784472, 0.00643437,    -8.4523651,   -2.0007162 ],
        [ 3.51695251,   0.26967622, 0.00561315, 0.26967622, 0.00561315,    -9.0621709,   -2.0005731 ],
        [ 3.84566235,   0.27140021, 0.00490234, 0.27140021, 0.00490234,    -9.7197564,   -2.0004385 ],
        [ 4.2005339,    0.27302665, 0.00428664, 0.27302665, 0.00428664,    -10.42963,    -2.0003133 ],
        [ 4.58441845,   0.2745659,  0.00375229, 0.2745659,  0.00375229,    -11.197497,   -2.000217 ],
        [ 5.0003596,    0.27602648, 0.00328775, 0.27602648, 0.00328775,    -12.029452,   -2.0001456 ],
        [ 5.45179388,   0.27741601, 0.00288314, 0.27741601, 0.00288314,    -12.932372,   -2.0000928 ],
        [ 5.94235081,   0.27874061, 0.00253018, 0.27874061, 0.00253018,    -13.91352,    -2.0000568 ],
        [ 6.47645386,   0.28000657, 0.00222145, 0.28000787, 0.0022257993,  -14.981749,   -2.0000328 ],
        [ 7.05892076,   0.28121896, 0.00195087, 0.28130469, 0.0023031204,  -16.146696,   -2.0000184 ],
        [ 7.69476446,   0.28238168, 0.00171207, 0.28274831, 0.0021656396,  -17.418392,   -2.0000097 ],
        [ 8.38979402,   0.28349821, 0.00150607, 0.28420268, 0.0020172183,  -18.808455,   -2.0000049 ],
        [ 9.15001541,   0.28457103, 0.001323,   0.28568216, 0.0018738716,  -20.3289,     -2.0000023 ],
        [ 9.98243239,   0.28560286, 0.00116198, 0.28718553, 0.0017419864,  -21.993736,   -2.000001 ],
        [ 10.89424716,  0.28659542, 0.00102031, 0.28870966, 0.0016028057,  -23.817366,   -2.0000004 ],
        [ 11.8934609,   0.28755036, 0.00089563, 0.29024511, 0.0014713244,  -25.815793,   -2.0000002 ],
        [ 12.98847423,  0.28846884, 0.00078595, 0.29178457, 0.0013415891,  -28.00582,    -2.0000001 ],
        [ 14.1884874,   0.28935198, 0.00068948, 0.29331761, 0.0012154296,  -30.405847,   -2.0000001 ],
        [ 15.50350055,  0.29020084, 0.00060466, 0.29483558, 0.0010966484,  -33.035873,   -2 ],
        [ 16.9443137,   0.29101637, 0.00053013, 0.29633385, 0.00098566517, -35.917499,   -2 ],
        [ 18.45007586,  0.2917657,  0.00046737, 0.29774289, 0.00095158101, -38.929024,   -2 ],
        [ 20.18212331,  0.292523,   0.00040926, 0.29928386, 0.00083224511, -42.393118,   -2 ],
        [ 22.05632948,  0.29324106, 0.00035885, 0.30074329, 0.00072896684, -46.141531,   -2 ],
        [ 24.07542792,  0.29391996, 0.00031522, 0.30212176, 0.0006397161,  -50.179728,   -2 ],
        [ 26.23948437,  0.29455985, 0.0002775,  0.30341983, 0.00056270227, -54.507841,   -2 ],
        [ 28.6989312,   0.29519836, 0.00024302, 0.30471409, 0.00049239377, -59.426734,   -2 ],
        [ 31.33289533,  0.295798,   0.00021337, 0.30592865, 0.0004320316,  -64.694663,   -2 ],
        [ 34.13469912,  0.29635892, 0.00018794, 0.30706407, 0.00038029964, -70.29827,    -2 ],
        [ 37.31698478,  0.2969186,  0.00016467, 0.3081963,  0.00033302931, -76.662841,   -2 ],
        [ 40.69282789,  0.29743977, 0.00014482, 0.30925009, 0.00029273491, -83.414528,   -2 ],
        [ 44.24158952,  0.29792265, 0.00012791, 0.31022602, 0.00025846201, -90.512051,   -2 ],
        [ 48.26046879,  0.29840445, 0.00011242, 0.31119939, 0.00022706749, -98.549809,   -2 ],
        [ 52.46128663,  0.29884822, 9.931e-05,  0.31209557, 0.00020052059, -106.95145,   -2 ],
        [ 57.21818094,  0.299291,   8.729e-05,  0.31298949, 0.00017619619, -116.46523,   -2 ],
        [ 62.15340664,  0.29969601, 7.719e-05,  0.3138069,  0.0001557581,  -126.33569,   -2 ],
        [ 67.73661891,  0.30010014, 6.792e-05,  0.31462233, 0.00013701391, -137.50211,   -2 ],
        [ 74.08516567,  0.30050338, 5.944e-05,  0.31543577, 0.00011988233, -150.1992,    -2 ],
        [ 80.64212037,  0.30086916, 5.239e-05,  0.3161735,  0.00010563988, -163.31311,   -2 ],
        [ 88.08544877,  0.30123417, 4.593e-05,  0.31690952, 9.2604495e-05, -178.19977,   -2 ],
        [ 95.67848133,  0.301562,   4.061e-05,  0.31757047, 8.1855931e-05, -193.38583,   -2 ],
        [ 104.27165039, 0.30188917, 3.572e-05,  0.31822999, 7.1995217e-05, -210.57217,   -2 ],
        [ 114.04710498, 0.30221568, 3.125e-05,  0.31888808, 6.298017e-05,  -230.12308,   -2 ],
        [ 123.91000396, 0.30250533, 2.761e-05,  0.31947182, 5.5643319e-05, -249.84888,   -2 ],
        [ 135.07945489, 0.30279443, 2.428e-05,  0.32005439, 4.8913516e-05, -272.18778,   -2 ],
        [ 146.11157231, 0.30304694, 2.159e-05,  0.32056318, 4.350026e-05,  -294.25202,   -2 ],
        [ 158.52261288, 0.30329902, 1.912e-05,  0.32107105, 3.8510924e-05, -319.0741,    -2 ],
        [ 172.55075943, 0.30355066, 1.684e-05,  0.321578,   3.3926941e-05, -347.13039,   -2 ],
        [ 188.48779585, 0.30380185, 1.476e-05,  0.32208402, 2.9729563e-05, -379.00446,   -2 ],
        [ 203.93845363, 0.3040168,  1.312e-05,  0.32251701, 2.6425798e-05, -409.90578,   -2 ]
    ],

    impulse_table => [
        [
            -4.7731922, 13.814838,   -4.7736824,    0.035693817, -0.0025496796, -0.10209685,
            -2.3970216, -0.24351436, -0.0034374282, -0.33340822, 0.40958802,    0.99971801,
            -3.5294346
        ],
        [
            -4.6051696, 13.310772,   -4.6058003,    0.038434733, -0.0030161767, -0.10055036,
            -2.395475,  -0.32752902, -0.0040663509, -0.33186157, 0.40958778,    0.99960686,
            -3.4341673
        ],
        [
            -4.1997045, 12.094382,   -4.2008635,   0.045749898, -0.0045242569, -0.09555151,
            -2.3904754, -0.53027908, -0.006099518, -0.32686156, 0.40958653,    0.99912399,
            -3.2097678
        ],
        [
            -3.9120224, 11.231346,  -3.9138036,    0.051540363, -0.0060323178, -0.090554153,
            -2.3854763, -0.6741533, -0.0081326652, -0.32186156, 0.40958414,    0.99845601,
            -3.0561007
        ],
        [
            -3.5065573, 10.014989,   -3.5098388,  0.060482892, -0.0090482929, -0.08056721,
            -2.3754808, -0.87702038, -0.01219881, -0.31186156, 0.40957445,    0.99658348,
            -2.8492642
        ],
        [
            -2.9957317, 8.4827099, -3.0028048,   0.072702044, -0.015078478, -0.060652934,
            -2.3555099, -1.133112, -0.020329292, -0.29186154, 0.40952469,   0.99082093,
            -2.6099017
        ],
        [
            -2.3025845, 6.4050321,  -2.3227043,   0.089159144, -0.030107441, -0.011707353,
            -2.3058707, -1.4854312, -0.040607762, -0.24286432, 0.40908193,   0.96630585,
            -2.3421997
        ],
        [
            -1.8971194, 5.193379,   -1.9343075,   0.097176807, -0.044911537, 0.034767592,
            -2.2571085, -1.7027383, -0.060652528, -0.19603142, 0.40789066,   0.93069453,
            -2.2341259
        ],
        [
            -1.6094373, 4.339381,   -1.6670074,   0.1016189,   -0.05915966, 0.077194593,
            -2.2099274, -1.8729994, -0.080100664, -0.15317032, 0.40561391,  0.88779427,
            -2.1912896
        ],
        [
            -1.3862937, 3.6843633,  -1.4670771,   0.10441401,  -0.072421179, 0.11451167,
            -2.1651122, -2.0249388, -0.098428991, -0.11583524, 0.40197827,   0.84092476,
            -2.1845904
        ],
        [
            -1.2039722, 3.15784,    -1.3103727,  0.10648743,  -0.084285306, 0.14641185,
            -2.1234322, -2.1715457, -0.11505077, -0.08503368, 0.39680618,   0.7928252,
            -2.2007965
        ],
        [
            -1.0498215, 2.7220571,  -1.1838225, 0.10829942,   -0.094474895, 0.17314844,
            -2.0855152, -2.3192086, -0.129467,  -0.061402496, 0.39003767,   0.74557753,
            -2.2324176
        ],
        [
            -0.91629012, 2.3541168,  -1.0794562, 0.11006104,   -0.10289996, 0.19524111,
            -2.0517222,  -2.4706478, -0.141405,  -0.044280433, 0.38173376,  0.70060554,
            -2.2745823
        ],
        [
            -0.79850708, 2.0388139,  -0.99199851, 0.11183967,   -0.10964602, 0.2133055,
            -2.0220816,  -2.6262322, -0.15087559, -0.032978978, 0.3720601,   0.65875338,
            -2.3238556
        ],
        [
            -0.69314656, 1.7653666,  -0.91775099, 0.11362457, -0.11492024, 0.22797959, -1.9963247, -2.784863,
            -0.1581291,  -0.0261218, 0.36125728,  0.62041082, -2.3777213
        ],
        [
            -0.59783638, 1.5257746,   -0.85401328, 0.11537281, -0.1189843, 0.2398782, -1.9740004, -2.944723,
            -0.1635485,  -0.02316644, 0.34960674,  0.5856476,  -2.4343172
        ],
        [
            -0.51082501, 1.3139198,  -0.79875688, 0.11703867,   -0.12209854, 0.24955561,
            -1.9545946,  -3.1038717, -0.16754079, -0.022741677, 0.3374001,   0.55433085,
            -2.4922723
        ],
        [
            -0.4307823, 1.1250307,  -0.75042641, 0.11858798,   -0.12448927, 0.25748218,
            -1.9376192, -3.2606204, -0.17046698, -0.024033767, 0.32491588,  0.526215,
            -2.5505916
        ],
        [
            -0.35667433, 0.95533497, -0.7078107,  0.12000166,   -0.12633754, 0.26403834,
            -1.9226537,  -3.4136892, -0.17261476, -0.026404539, 0.31240435,  0.50100379,
            -2.6085665
        ],
        [
            -0.28768146, 0.80181977,  -0.66995504, 0.1212737,  -0.12778113, 0.26952204, -1.9093532, -3.5622167,
            -0.17419963, -0.03000892, 0.30007939,  0.47838954, -2.6657046
        ],
        [
            -0.22314294, 0.66205643, -0.63609848, 0.12240734,   -0.12892207, 0.27416209,
            -1.8974419,  -3.7056945, -0.17537763, -0.033797667, 0.28811552,  0.4580759,
            -2.7216744
        ],
        [
            -0.16251831, 0.53407214, -0.60562898, 0.12341138,   -0.12983503, 0.27813284,
            -1.8867018,  -3.8438832, -0.17625994, -0.038024447, 0.27664836,  0.43978948,
            -2.7762626
        ],
        [
            -0.1053599, 0.41625312, -0.57805033, 0.12429738,   -0.13057453, 0.28156701,
            -1.8769586, -3.9767338, -0.17692551, -0.042447784, 0.26577714,  0.42328469,
            -2.8293421
        ],
        [
            -0.051292678, 0.30726705, -0.55295623, 0.12507784,   -0.1311806, 0.28456628,
            -1.8680739,   -4.1043251, -0.17743076, -0.046786891, 0.25556771, 0.40834449,
            -2.8808471
        ],
        [
            6.1591356e-07, 0.20600804, -0.53001206, 0.12576505,   -0.13168279, 0.28720895,
            -1.8599351,    -4.2268165, -0.17781597, -0.051167006, 0.246049,    0.39477908,
            -2.9307553
        ],
        [
            0.04478674,  0.11923539,  -0.51064214, 0.12632255, -0.13207064, 0.28936922, -1.8530536, -4.3347297,
            -0.17808877, -0.05518108, 0.23793117,  0.38341743, -2.9750804
        ],
        [
            0.18232217, -0.13824279, -0.45480752, 0.1277941,    -0.13302182, 0.29525746,
            -1.833209,  -4.6702331,  -0.17863905, -0.068275388, 0.21450336,  0.35127885,
            -3.1148518
        ],
        [
            0.26236488, -0.28228878, -0.42467592, 0.12849434,   -0.13344387, 0.29823998,
            -1.8225333, -4.8673646,  -0.17881067, -0.076147495, 0.2020775,   0.33440814,
            -3.1982303
        ],
        [
            0.40546572, -0.53028221, -0.37471443, 0.12949668,   -0.13402441, 0.30290577,
            -1.805019,  -5.2213099,  -0.17893873, -0.091158009, 0.18219455,  0.30734618,
            -3.3499395
        ],
        [
            0.53062887,  -0.73821907, -0.33471703, 0.13015218, -0.1343979, 0.3064011, -1.7913238, -5.5310651,
            -0.17892726, -0.10454193, 0.16715133,  0.28663079, -3.4845245
        ],
        [
            0.6931478,   -0.99737444, -0.28729253, 0.13076273, -0.13475482, 0.31027537, -1.7757323, -5.9318309,
            -0.17881437, -0.12269777, 0.15052794,  0.26331136, -3.6607296
        ],
        [
            0.99325239,  -1.4493786,  -0.21088596, 0.13139908, -0.13518005, 0.31589118, -1.7530066, -6.663787,
            -0.17848293, -0.15718425, 0.12680217,  0.22888506, -3.9872391
        ],
        [
            1.0986129,   -1.6013009,  -0.18693924, 0.13152211, -0.13528454, 0.31748545, -1.7467143, -6.9176926,
            -0.17835708, -0.16914383, 0.12016464,  0.21894293, -4.1016135
        ],
        [
            1.208961,   -1.7572325,  -0.16322655, 0.13161382, -0.13537762, 0.31898249, -1.7409565, -7.1817932,
            -0.1782276, -0.18198804, 0.11396373,  0.20950767, -4.2210754
        ],
        [
            1.386295,    -2.0017887,  -0.12773183, 0.13170275, -0.13550052, 0.32106564, -1.7332925, -7.6023003,
            -0.17803076, -0.20302354, 0.10535719,  0.19615046, -4.4121924
        ],
        [
            1.6094385, -2.3004646, -0.087029708, 0.13174773,  -0.13562177, 0.32321315,
            -1.726005, -8.1247805, -0.17781035,  -0.22917139, 0.096445173, 0.18195902,
            -4.6509802
        ],
        [
            1.7917601,  -2.5381886, -0.056568765, 0.13175181,  -0.13570117, 0.32464633,
            -1.7216141, -8.546496,  -0.177655,    -0.25033414, 0.090408042, 0.17211152,
            -4.84463
        ],
        [
            1.9459108,  -2.7354026, -0.032508045, 0.13174229,  -0.1357579,  0.32567141,
            -1.7187664, -8.8996813, -0.17754049,  -0.26848254, 0.085999647, 0.16478903,
            -5.007364
        ],
        [
            2.0794422,  -2.9037741, -0.012780934, 0.13172846,  -0.13580075, 0.32644116,
            -1.7168176, -9.2033043, -0.17745253,  -0.28463032, 0.08261002,  0.15907757,
            -5.1476309
        ],
        [
            2.3027878,   -3.1809156,  0.018164341, 0.13169964, -0.13586597, 0.32752346, -1.7150159, -9.7067258,
            -0.17732799, -0.31012085, 0.077677136, 0.1491248,  -5.3894078
        ],
        [
            2.9957966,   -4.0129766,  0.1011374,   0.13161045, -0.13599041, 0.32968801, -1.7111815, -11.239511,
            -0.17707007, -0.39205114, 0.066502878, 0.13026484, -6.0996288
        ],
        [
            3.9122092,   -5.0692993,  0.18946312,  0.13153791, -0.13607019, 0.33099543, -1.709907, -13.217658,
            -0.17691288, -0.50030445, 0.057265942, 0.11354513, -7.0287596
        ],
        [
            4.6052394,   -5.8467506,  0.24530528,  0.13151107, -0.1361034, 0.33143567, -1.7103152, -14.689535,
            -0.17686309, -0.58210057, 0.052494309, 0.10451956, -7.72676
        ],
        [
            5.2985264,   -6.6119071,  0.29439967,  0.13149733,  -0.13611488, 0.33165283, -1.7100144, -16.148096,
            -0.17683569, -0.66388588, 0.048822393, 0.097418559, -8.4228011
        ],
        [
            6.2146099,   -7.6089595,  0.35141361,  0.13148911,  -0.13612345, 0.33178414, -1.7100133, -18.060472,
            -0.17682005, -0.77191519, 0.045046384, 0.090005244, -9.3406989
        ],
        [
            6.9079066,   -8.3555605,  0.38991013,  0.13148641,  -0.13612633, 0.33182796, -1.710014, -19.499546,
            -0.17681484, -0.85365502, 0.042739555, 0.085436167, -10.034654
        ],
        [
            7.6009593,   -9.0965686,  0.42517474,  0.13148509,  -0.13612777, 0.33184986, -1.7100146, -20.932757,
            -0.17681223, -0.93535788, 0.040771287, 0.081521421, -10.728056
        ],
        [
            8.5171982,   -10.069655, 0.46772338, 0.13148431,  -0.13612861, 0.331863, -1.710015, -22.821038,
            -0.17634309, -1.0177229, 0.03855756, 0.077106791, -11.644518
        ],
        [
            9.2104168,   -10.801831, 0.49733062,  0.13148405,  -0.13612889, 0.33186739, -1.7100152, -24.245733,
            -0.17474513, -1.0642402, 0.037110008, 0.074215892, -12.337815
        ],
        [
            10.819846,   -12.491193, 0.55916139,  0.13148385,  -0.13612909, 0.33187089, -1.7100153, -27.543204,
            -0.16901945, -1.1660602, 0.034304248, 0.068607686, -13.947311
        ],
        [
            11.513056,   -13.215155, 0.58331946,  0.13148382,  -0.13612911, 0.33187133, -1.7100153, -28.959919,
            -0.16620877, -1.2076364, 0.033280543, 0.066560683, -14.64053
        ],
        [
            13.122476,   -14.88928,  0.63478231,  0.1314838,   -0.13612923, 0.33187168, -1.7100153, -32.242601,
            -0.15959158, -1.2997413, 0.031222238, 0.062444397, -16.249958
        ],
        [
            13.815683,   -15.607919, 0.65523463, 0.13148379, -0.13612914, 0.33187172, -1.7100153, -33.654143,
            -0.15680484, -1.3377274, 0.0304475,  0.06089496, -16.943166
        ],
        [
            16.118106,  -17.986433, 0.71714785,  0.13148378,  -0.1361292, 0.33187175, -1.7100153, -38.334276,
            -0.1481074, -1.4577545, 0.028239355, 0.056478707, -19.24559
        ],
        [
            18.420726,   -20.354932, 0.77157632,  0.13148377,  -0.13612933, 0.33187216, -1.7100149, -43.004814,
            -0.14037876, -1.5698341, 0.026454263, 0.052908526, -21.54821
        ]
    ],

    tail_shock_table => [
        [ 5.6706841, -0.80301,    -0.04127305,  -0.04127305 ],
        [ 6.2117637, -0.84750697, -0.054606993, -0.02554784 ],
        [ 6.9064289, -0.9017431,  -0.057842512, -0.018089506 ],
        [ 7.6001941, -0.95319011, -0.058832775, -0.013543481 ],
        [ 8.5168789, -1.0176869,  -0.058782007, -0.0096344016 ],
        [ 9.2102524, -1.0642222,  -0.058248972, -0.0075955259 ],
        [ 10.819811, -1.1660566,  -0.056339955, -0.0045503158 ],
        [ 11.513038, -1.2076346,  -0.055402994, -0.0036880099 ],
        [ 13.122472, -1.2997409,  -0.053197206, -0.0022872159 ],
        [ 13.815682, -1.3377272,  -0.052268287, -0.0018619694 ],
        [ 16.118106, -1.457832,   -0.049372163, -0.00090960286 ]
    ],

    blast_info => [
        0.11055004, 2.4054787,  0.34179936, -0.40663495, 2.9410003,  0.0007202227, 2.5808454,   -0.30161774,
        0.44558731, 0.50893787, 0.66456498, -0.28708537, 0.14790029, 0.40958878,   0.043827926, -0.045376369,
        289.43,     -0.80301,   2736.0762,  -0.9722916
    ],

};
1;
