package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=5
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S5'} = {

    shock_table_info => [ 2, 5, 1.13e-06, 4000, 9.2e-07, 10, 1.55e-07, 4.7e-07, 5e-07, 27.612, -0.78157 ],

    shock_table => [
        [ -3.98450494,  12.00033703,   -2.99996929, -3.98578557, 0.99807785 ],
        [ -3.70145664,  11.15120316,   -2.99992821, -3.70341528, 0.99705921 ],
        [ -3.52253607,  10.61445509,   -2.99987722, -3.52509843, 0.99615165 ],
        [ -3.38929735,  10.21475503,   -2.99985422, -3.39242745, 0.99529773 ],
        [ -2.99764269,  9.03990028,    -2.99952665, -3.0032819,  0.99151866 ],
        [ -2.69466456,  8.13119986,    -2.99882574, -2.7035616,  0.98660028 ],
        [ -2.43468111,  7.35171193,    -2.99744377, -2.44784644, 0.9801388 ],
        [ -2.23182306,  6.74385193,    -2.99530554, -2.24970625, 0.97297731 ],
        [ -2.04663244,  6.18944071,    -2.99184598, -2.07029591, 0.96418158 ],
        [ -1.88159557,  5.69606729,    -2.9866893,  -1.91198037, 0.95393406 ],
        [ -1.73783211,  5.26716151,    -2.97964881, -1.77562211, 0.94263314 ],
        [ -1.60048073,  4.85854908,    -2.96957856, -1.64703658, 0.92926203 ],
        [ -1.47395993,  4.48363676,    -2.9561273,  -1.53038577, 0.91424392 ],
        [ -1.35561591,  4.13478297,    -2.93852713, -1.42315755, 0.89741666 ],
        [ -1.24320274,  3.80566694,    -2.91584556, -1.32330806, 0.8785593 ],
        [ -1.13029967,  3.47808834,    -2.88559276, -1.22533277, 0.85644888 ],
        [ -1.01631899,  3.15136268,    -2.84569759, -1.12915223, 0.83060526 ],
        [ -0.89501737,  2.80935319,    -2.79106699, -1.03027158, 0.79901212 ],
        [ -0.77611021,  2.48132506,    -2.72402925, -0.93730766, 0.76397262 ],
        [ -0.66717127,  2.18846688,    -2.6506683,  -0.85598094, 0.72861583 ],
        [ -0.56528284,  1.92230298,    -2.57252583, -0.78353091, 0.69319794 ],
        [ -0.46872671,  1.6777761,     -2.49145916, -0.71828334, 0.65809999 ],
        [ -0.37557157,  1.44951439,    -2.40863275, -0.65858944, 0.62341852 ],
        [ -0.28277788,  1.22993995,    -2.32365845, -0.60235294, 0.58867424 ],
        [ -0.18977793,  1.01782425,    -2.23806225, -0.54921277, 0.55423508 ],
        [ -0.095951,    0.81183544,    -2.15314439, -0.49880679, 0.52040265 ],
        [ -0.00033263,  0.60997526,    -2.06969186, -0.45064139, 0.48731311 ],
        [ 0.09797292,   0.41055061,    -1.98838948, -0.40433512, 0.45510406 ],
        [ 0.19963742,   0.21243991,    -1.9099584,  -0.35966944, 0.42396496 ],
        [ 0.30580676,   0.01371952,    -1.83464044, -0.31627454, 0.39392719 ],
        [ 0.41710028,   -0.186408,     -1.76298195, -0.27405891, 0.36516911 ],
        [ 0.53465896,   -0.3895957,    -1.69510355, -0.23277224, 0.33771888 ],
        [ 0.65961225,   -0.59732396,   -1.63113163, -0.19223536, 0.31162011 ],
        [ 0.79302019,   -0.81084246,   -1.571213,   -0.15234391, 0.28693428 ],
        [ 0.93611186,   -1.03157941,   -1.51539594, -0.11298641, 0.26368894 ],
        [ 1.08595901,   -1.25478695,   -1.46500781, -0.07510085, 0.24245839 ],
        [ 1.24165417,   -1.4792922,    -1.42005056, -0.03887928, 0.22328022 ],
        [ 1.4044327,    -1.70708425,   -1.37979438, -0.00398384, 0.20588091 ],
        [ 1.57551569,   -1.93996999,   -1.34365025, 0.02985119,  0.1900395 ],
        [ 1.75633636,   -2.17990783,   -1.3111019,  0.0628733,   0.17556 ],
        [ 1.94814941,   -2.42850356,   -1.28176327, 0.09524486,  0.16229869 ],
        [ 2.15087537,   -2.68561143,   -1.25546911, 0.12689144,  0.15021009 ],
        [ 2.36781407,   -2.95531638,   -1.23166186, 0.15823874,  0.13906512 ],
        [ 2.59948803,   -3.23810414,   -1.21020461, 0.18924052,  0.12882519 ],
        [ 2.84777348,   -3.53610689,   -1.19083794, 0.22002536,  0.119393 ],
        [ 3.11394983,   -3.85068957,   -1.17338575, 0.25061989,  0.11070943 ],
        [ 3.40029786,   -4.1843644,    -1.15763124, 0.28114427,  0.1026932 ],
        [ 3.7090886,    -4.5395677,    -1.14339909, 0.31168118,  0.09528107 ],
        [ 4.0421847,    -4.91822668,   -1.13055591, 0.34224845,  0.08842973 ],
        [ 4.40284549,   -5.32381565,   -1.1189367,  0.3729663,   0.08207682 ],
        [ 4.79352046,   -5.75884212,   -1.10843083, 0.40385128,  0.07618702 ],
        [ 5.21725317,   -6.22644573,   -1.09892534, 0.43494619,  0.07072187 ],
        [ 5.6778804,    -6.7305961,    -1.09030911, 0.46632197,  0.06564148 ],
        [ 6.18003184,   -7.27607095,   -1.08247985, 0.49806434,  0.06090819 ],
        [ 6.7281306,    -7.86736599,   -1.07535887, 0.530207,    0.05649588 ],
        [ 7.32639422,   -8.50871777,   -1.06888185, 0.56274428,  0.05238519 ],
        [ 7.98043561,   -9.20582223,   -1.06298018, 0.59571949,  0.04855179 ],
        [ 8.69666366,   -9.96517066,   -1.05759277, 0.62917806,  0.04497363 ],
        [ 9.48088412,   -10.79256965,  -1.05267441, 0.66310381,  0.04163678 ],
        [ 10.33950057,  -11.69442676,  -1.04818294, 0.69748301,  0.03852739 ],
        [ 11.28131508,  -12.6796267,   -1.0440727,  0.73236586,  0.03562699 ],
        [ 12.31292871,  -13.75471011,  -1.0403158,  0.76768891,  0.03292762 ],
        [ 13.44394199,  -14.92931841,  -1.03687786, 0.80347072,  0.03041498 ],
        [ 14.68355516,  -16.21263756,  -1.03373224, 0.83968575,  0.02807878 ],
        [ 16.04176832,  -17.61465328,  -1.0308545,  0.87630852,  0.02590896 ],
        [ 17.52998148,  -19.14676634,  -1.02822133, 0.9133265,   0.02389505 ],
        [ 19.12804778,  -20.78798642,  -1.02585581, 0.95000742,  0.02206155 ],
        [ 20.88048492,  -22.58377347,  -1.02368218, 0.98713313,  0.0203557 ],
        [ 22.68648315,  -24.43079128,  -1.02179628, 1.02250834,  0.01885854 ],
        [ 24.73164266,  -26.51863266,  -1.01999521, 1.05955905,  0.01741329 ],
        [ 26.8124472,   -28.63938799,  -1.01844603, 1.09445298,  0.01615757 ],
        [ 29.16065497,  -31.02912218,  -1.01696415, 1.13093299,  0.01494508 ],
        [ 31.82339595,  -33.73509511,  -1.01554891, 1.16912888,  0.01377636 ],
        [ 34.50087179,  -36.45254951,  -1.01434635, 1.20463945,  0.01277468 ],
        [ 37.52222925,  -39.51546525,  -1.01319551, 1.24173764,  0.01180837 ],
        [ 40.94818852,  -42.9846935,   -1.01209596, 1.28055429,  0.01087779 ],
        [ 44.33609434,  -46.41197895,  -1.01117557, 1.31604435,  0.01009312 ],
        [ 48.15180648,  -50.26860541,  -1.01029383, 1.35307481,  0.00933629 ],
        [ 52.46968569,  -54.62905935,  -1.00945046, 1.39177125,  0.00860753 ],
        [ 57.38091733,  -59.58466806,  -1.00864517, 1.43227546,  0.00790706 ],
        [ 62.1475857,   -64.39092998,  -1.00798502, 1.46855364,  0.00732933 ],
        [ 67.51967218,  -69.80416682,  -1.00735247, 1.50639242,  0.00677267 ],
        [ 73.60327287,  -75.93060361,  -1.00674736, 1.54592063,  0.0062372 ],
        [ 79.30959975,  -81.67401958,  -1.00626392, 1.5802559,   0.00580725 ],
        [ 85.69155461,  -88.09442985,  -1.00579932, 1.61595991,  0.00539217 ],
        [ 92.85932525,  -95.30212807,  -1.00535345, 1.65313845,  0.00499202 ],
        [ 100.94644968, -103.43077096, -1.0049262,  1.69191016,  0.00460688 ],
        [ 110.11600528, -112.64356966, -1.00451749, 1.73240872,  0.0042368 ],
        [ 120.56880314, -123.14148675, -1.00412721, 1.77478554,  0.00388186 ],
        [ 130.02065903, -132.63090389, -1.0038282,  1.81015436,  0.00360883 ],
        [ 140.60844988, -143.2576662,  -1.00354086, 1.84693466,  0.00334555 ],
        [ 152.52069403, -155.21040349, -1.00326516, 1.88523788,  0.00309205 ],
        [ 165.98661661, -168.71846618, -1.00300104, 1.92518925,  0.00284836 ],
        [ 181.28715455, -184.06293266, -1.00274845, 1.96693021,  0.0026145 ],
        [ 194.17243501, -196.98242901, -1.00256655, 1.99950612,  0.00244558 ],
        [ 208.46045186, -211.30583359, -1.00239108, 2.03325447,  0.00228221 ]
    ],

    energy_table => [
        [ -3.98450494,  5.518e-05,  0.00012839, 5.518e-05,  0.00012839,    5.0728292,   -1.5409864 ],
        [ -3.70145664,  0.00010638, 0.00024587, 0.00010638, 0.00024587,    4.63555,     -1.5486423 ],
        [ -3.52253607,  0.0001607,  0.00036952, 0.0001607,  0.00036952,    4.358042,    -1.5532553 ],
        [ -3.38929735,  0.00021817, 0.00049951, 0.00021817, 0.00049951,    4.1508657,   -1.5573375 ],
        [ -2.99764269,  0.00053127, 0.0011971,  0.00053127, 0.0011971,     3.5381502,   -1.5712218 ],
        [ -2.69466456,  0.00104574, 0.00231669, 0.00104574, 0.00231669,    3.060512,    -1.5853824 ],
        [ -2.43468111,  0.00185043, 0.00402112, 0.00185043, 0.00402112,    2.6463523,   -1.6010205 ],
        [ -2.23182306,  0.00286391, 0.00610493, 0.00286391, 0.00610493,    2.3203066,   -1.6145296 ],
        [ -2.04663244,  0.0042326,  0.0088271,  0.0042326,  0.0088271,     2.0200819,   -1.6305186 ],
        [ -1.88159557,  0.00594784, 0.01210893, 0.00594784, 0.01210893,    1.7496106,   -1.6475004 ],
        [ -1.73783211,  0.00794166, 0.01576093, 0.00794166, 0.01576093,    1.5116774,   -1.6646252 ],
        [ -1.60048073,  0.0103903,  0.02002468, 0.0103903,  0.02002468,    1.2817801,   -1.6838058 ],
        [ -1.47395993,  0.01320935, 0.02464826, 0.01320935, 0.02464826,    1.0675763,   -1.7030842 ],
        [ -1.35561591,  0.01641094, 0.02954412, 0.01641094, 0.02954412,    0.86491373,  -1.7223589 ],
        [ -1.24320274,  0.02001387, 0.0346147,  0.02001387, 0.0346147,     0.67024379,  -1.741121 ],
        [ -1.13029967,  0.02422232, 0.03996005, 0.02422232, 0.03996005,    0.47260168,  -1.7595682 ],
        [ -1.01631899,  0.02908673, 0.0453736,  0.02908673, 0.0453736,     0.27100682,  -1.7770245 ],
        [ -0.89501737,  0.03492527, 0.05079954, 0.03492527, 0.05079954,    0.054373111, -1.7936406 ],
        [ -0.77611021,  0.04124693, 0.05537124, 0.04124693, 0.05537124,    -0.15980473, -1.8077307 ],
        [ -0.66717127,  0.04746303, 0.05856719, 0.04746303, 0.05856719,    -0.35738691, -1.8190806 ],
        [ -0.56528284,  0.05353711, 0.06047609, 0.05353711, 0.06047609,    -0.54324324, -1.8290971 ],
        [ -0.46872671,  0.05942086, 0.06122029, 0.05942086, 0.06122029,    -0.72030956, -1.8388835 ],
        [ -0.37557157,  0.06511867, 0.06095163, 0.06511867, 0.06095163,    -0.89206631, -1.8492206 ],
        [ -0.28277788,  0.07072722, 0.05979116, 0.07072722, 0.05979116,    -1.064166,   -1.8606897 ],
        [ -0.18977793,  0.07620342, 0.05786141, 0.07620342, 0.05786141,    -1.2377731,  -1.8732965 ],
        [ -0.095951,    0.08151642, 0.05530217, 0.08151642, 0.05530217,    -1.4141591,  -1.8867902 ],
        [ -0.00033263,  0.0866611,  0.05224663, 0.0866611,  0.05224663,    -1.5952416,  -1.9008083 ],
        [ 0.09797292,   0.09163059, 0.04882353, 0.09163059, 0.04882353,    -1.7828094,  -1.914907 ],
        [ 0.19963742,   0.09640845, 0.04516188, 0.09640845, 0.04516188,    -1.9782127,  -1.9286099 ],
        [ 0.30580676,   0.10100083, 0.0413641,  0.10100083, 0.0413641,     -2.1837003,  -1.941551 ],
        [ 0.41710028,   0.10538964, 0.03754029, 0.10538964, 0.03754029,    -2.4004918,  -1.9533861 ],
        [ 0.53465896,   0.10957819, 0.0337702,  0.10957819, 0.0337702,     -2.6308074,  -1.9639116 ],
        [ 0.65961225,   0.113566,   0.03012299, 0.113566,   0.03012299,    -2.8768364,  -1.9729851 ],
        [ 0.79302019,   0.11734862, 0.02665819, 0.11734862, 0.02665819,    -3.1406205,  -1.9806271 ],
        [ 0.93611186,   0.12092573, 0.02341865, 0.12092573, 0.02341865,    -3.4245464,  -1.9868593 ],
        [ 1.08595901,   0.12421119, 0.02050977, 0.12421119, 0.02050977,    -3.7226851,  -1.9915743 ],
        [ 1.24165417,   0.12719933, 0.01794703, 0.12719933, 0.01794703,    -4.0330799,  -1.9950222 ],
        [ 1.4044327,    0.12993194, 0.01569398, 0.12993194, 0.01569398,    -4.3580677,  -1.9974828 ],
        [ 1.57551569,   0.1324426,  0.01371671, 0.1324426,  0.01371671,    -4.6999793,  -1.9991644 ],
        [ 1.75633636,   0.13476112, 0.01198246, 0.13476112, 0.01198246,    -5.061594,   -2.0002497 ],
        [ 1.94814941,   0.13690918, 0.01046413, 0.13690918, 0.01046413,    -5.4453487,  -2.0008835 ],
        [ 2.15087537,   0.13889243, 0.00914518, 0.13889243, 0.00914518,    -5.8510254,  -2.0011937 ],
        [ 2.36781407,   0.14074679, 0.00798944, 0.14074679, 0.00798944,    -6.2851815,  -2.0012697 ],
        [ 2.59948803,   0.14247712, 0.00698255, 0.14247712, 0.00698255,    -6.7488199,  -2.0012272 ],
        [ 2.84777348,   0.14409804, 0.00610478, 0.14409804, 0.00610478,    -7.2456884,  -2.0011243 ],
        [ 3.11394983,   0.14561781, 0.00534114, 0.14561781, 0.00534114,    -7.7783192,  -2.0009426 ],
        [ 3.40029786,   0.14704861, 0.00467577, 0.14704861, 0.00467577,    -8.3512537,  -2.000756 ],
        [ 3.7090886,    0.14839969, 0.00409565, 0.14839969, 0.00409565,    -8.9690431,  -2.0005962 ],
        [ 4.0421847,    0.14967678, 0.00359036, 0.14967678, 0.00359036,    -9.6354062,  -2.0004419 ],
        [ 4.40284549,   0.15088924, 0.00314908, 0.15088924, 0.00314908,    -10.356859,  -2.0003118 ],
        [ 4.79352046,   0.15204155, 0.00276386, 0.15204155, 0.00276386,    -11.138309,  -2.0002138 ],
        [ 5.21725317,   0.15313881, 0.0024273,  0.15313881, 0.0024273,     -11.985846,  -2.0001416 ],
        [ 5.6778804,    0.15418659, 0.00213269, 0.15444864, 0.0030585485,  -12.907152,  -2.0000894 ],
        [ 6.18003184,   0.15519026, 0.00187414, 0.1559365,  0.0027222284,  -13.911488,  -2.0000534 ],
        [ 6.7281306,    0.15615302, 0.001647,   0.1573323,  0.0024146385,  -15.007707,  -2.0000313 ],
        [ 7.32639422,   0.15707666, 0.0014481,  0.15870142, 0.0021712875,  -16.204247,  -2.0000171 ],
        [ 7.98043561,   0.15796423, 0.00127201, 0.16005248, 0.0019650245,  -17.512338,  -2.0000088 ],
        [ 8.69666366,   0.15881849, 0.00111829, 0.16139241, 0.0017820989,  -18.944798,  -2.0000043 ],
        [ 9.48088412,   0.15964025, 0.00098242, 0.16272257, 0.0016154911,  -20.513241,  -2.000002 ],
        [ 10.33950057,  0.16043059, 0.0008629,  0.1640421,  0.0014624643,  -22.230475,  -2.0000009 ],
        [ 11.28131508,  0.16119188, 0.00075759, 0.16535074, 0.0013203544,  -24.114104,  -2.0000004 ],
        [ 12.31292871,  0.16192389, 0.00066496, 0.16664145, 0.0011859841,  -26.177332,  -2.0000002 ],
        [ 13.44394199,  0.16262818, 0.00058344, 0.16791373, 0.0010670995,  -28.439359,  -2.0000001 ],
        [ 14.68355516,  0.16330536, 0.00051175, 0.16924648, 0.0010373851,  -30.918585,  -2.0000001 ],
        [ 16.04176832,  0.16395605, 0.00044873, 0.17056449, 0.00090821802, -33.635012,  -2.0000001 ],
        [ 17.52998148,  0.16458113, 0.00039336, 0.17182875, 0.00079503014, -36.611438,  -2.0000001 ],
        [ 19.12804778,  0.1651702,  0.00034561, 0.17301859, 0.00069762625, -39.807571,  -2 ],
        [ 20.88048492,  0.16573757, 0.00030346, 0.17416324, 0.00061187079, -43.312445,  -2 ],
        [ 22.68648315,  0.16625278, 0.00026833, 0.17520161, 0.00054051263, -46.924442,  -2 ],
        [ 24.73164266,  0.16676737, 0.00023608, 0.17623778, 0.00047513053, -51.014761,  -2 ],
        [ 26.8124472,   0.16722991, 0.00020941, 0.17716839, 0.00042115258, -55.17637,   -2 ],
        [ 29.16065497,  0.16769184, 0.00018489, 0.17809709, 0.00037156429, -59.872785,  -2 ],
        [ 31.82339595,  0.16815309, 0.0001624,  0.17902382, 0.00032615678, -65.198267,  -2 ],
        [ 34.50087179,  0.16856249, 0.00014404, 0.17984587, 0.00028913816, -70.553219,  -2 ],
        [ 37.52222925,  0.16897128, 0.00012716, 0.18066628, 0.00025511832, -76.595934,  -2 ],
        [ 40.94818852,  0.16937943, 0.00011168, 0.18148499, 0.00022395684, -83.447852,  -2 ],
        [ 44.33609434,  0.16973599, 9.924e-05,  0.18219993, 0.00019892589, -90.223664,  -2 ],
        [ 48.15180648,  0.17009201, 8.778e-05,  0.18291351, 0.00017588539, -97.855088,  -2 ],
        [ 52.46968569,  0.17044745, 7.725e-05,  0.1836257,  0.00015474469, -106.49085,  -2 ],
        [ 57.38091733,  0.1708023,  6.762e-05,  0.18433645, 0.00013541402, -116.31331,  -2 ],
        [ 62.1475857,   0.17110596, 6.005e-05,  0.18494451, 0.00012021837, -125.84665,  -2 ],
        [ 67.51967218,  0.17140914, 5.308e-05,  0.18555148, 0.00010623233, -136.59082,  -2 ],
        [ 73.60327287,  0.17171184, 4.668e-05,  0.18615733, 9.3401499e-05, -148.75802,  -2 ],
        [ 79.30959975,  0.1719637,  4.176e-05,  0.18666133, 8.3552332e-05, -160.17067,  -2 ],
        [ 85.69155461,  0.1722152,  3.721e-05,  0.18716454, 7.4436934e-05, -172.93458,  -2 ],
        [ 92.85932525,  0.17246633, 3.301e-05,  0.18766693, 6.6024455e-05, -187.27013,  -2 ],
        [ 100.94644968, 0.17271709, 2.915e-05,  0.1881685,  5.8284339e-05, -203.44437,  -2 ],
        [ 110.11600528, 0.17296746, 2.56e-05,   0.18866923, 5.1186661e-05, -221.78349,  -2 ],
        [ 120.56880314, 0.17321743, 2.236e-05,  0.1891691,  4.4701401e-05, -242.68908,  -2 ],
        [ 130.02065903, 0.17341712, 1.998e-05,  0.18956838, 3.9934343e-05, -261.59279,  -2 ],
        [ 140.60844988, 0.17361655, 1.777e-05,  0.18996709, 3.5525123e-05, -282.76838,  -2 ],
        [ 152.52069403, 0.1738157,  1.574e-05,  0.19036524, 3.1459032e-05, -306.59286,  -2 ],
        [ 165.98661661, 0.17401459, 1.387e-05,  0.19076281, 2.7721163e-05, -333.52471,  -2 ],
        [ 181.28715455, 0.1742132,  1.216e-05,  0.19115981, 2.4296828e-05, -364.12578,  -2 ],
        [ 194.17243501, 0.17436197, 1.097e-05,  0.19145717, 2.1925611e-05, -389.89635,  -2 ],
        [ 208.46045186, 0.17451058, 9.87e-06,   0.19175419, 1.9716389e-05, -418.47238,  -2 ]
    ],

    impulse_table => [
        [
            -4.5893415, 13.814838,   -4.5898582,    0.041706048, -0.005229232, -0.037765624,
            -2.7466396, -0.32101906, -0.0066887939, -0.24071724, 0.47045576,   0.99993147,
            -4.1581466
        ],
        [
            -4.1997043, 12.64593,   -4.2006315,    0.050020165, -0.0077203562, -0.032953604,
            -2.7418036, -0.5158627, -0.0098753766, -0.23587677, 0.47045487,    0.99983149,
            -3.9005582
        ],
        [
            -3.9120222, 11.782891,   -3.91345,     0.057096276, -0.010293033, -0.02800571,
            -2.7368114, -0.65974999, -0.013166533, -0.23087677, 0.47045328,   0.99967188,
            -3.7153608
        ],
        [
            -3.5065571, 10.56652,    -3.5091817,  0.068569917, -0.015435196, -0.018200518,
            -2.7268412, -0.86266097, -0.01974623, -0.22110483, 0.47044683,   0.99916385,
            -3.463046
        ],
        [
            -2.9957315, 9.0341676,  -3.0013869,   0.085750097, -0.025692413, 0.00091968545,
            -2.7069776, -1.1189107, -0.032883231, -0.20190444, 0.47041375,   0.99731922,
            -3.1636071
        ],
        [
            -2.3025843, 6.9558359,  -2.3186542,   0.11406529,  -0.05091169, 0.044754693,
            -2.6580121, -1.4718096, -0.065361477, -0.15802799, 0.47011907,  0.98747652,
            -2.8056785
        ],
        [
            -1.8971192, 5.7424362,  -1.9267973,   0.13277661,  -0.074872735, 0.082361095,
            -2.610495,  -1.6891611, -0.096651331, -0.12174249, 0.46932366,   0.97021474,
            -2.637428
        ],
        [
            -1.6094371, 4.8851493,  -1.6553637,  0.14631079,   -0.096880051, 0.11464523,
            -2.5649212, -1.8577952, -0.12582652, -0.094159271, 0.4677927,    0.94643407,
            -2.547371
        ],
        [
            -1.3862936, 4.2250094,  -1.4507599,  0.15646726,   -0.11651442, 0.14267813,
            -2.5217029, -2.0053933, -0.15209502, -0.074652116, 0.46531887,  0.91755933,
            -2.5010392
        ],
        [
            -1.203972,  3.6914606,  -1.2889835,  0.16428909,   -0.13360717, 0.16713222,
            -2.4811602, -2.1441756, -0.17496907, -0.061897276, 0.46173858,  0.8851804,
            -2.4828206
        ],
        [
            -1.0498214, 3.246917,   -1.157114,   0.1704745,    -0.14817546, 0.18841197,
            -2.4435079, -2.2801894, -0.19429904, -0.054661556, 0.45694668,  0.85081497,
            -2.4840277
        ],
        [
            -0.91628996, 2.8688393,  -1.0473311,  0.17550783,   -0.16036796, 0.20681267,
            -2.4088444,  -2.4163567, -0.21021869, -0.051695798, 0.45090567,  0.81577027,
            -2.4991556
        ],
        [
            -0.79850692, 2.5424884,  -0.95449566, 0.17972035,   -0.17041943, 0.22260523,
            -2.3771531,  -2.5537902, -0.22305476, -0.052082104, 0.44364688,  0.78108053,
            -2.5244039
        ],
        [
            -0.69314641, 2.2575601,  -0.87502009, 0.18333124,   -0.17861059, 0.23606797,
            -2.3483164,  -2.6925114, -0.23323645, -0.055014488, 0.43526349,  0.74749909,
            -2.5570013
        ],
        [
            -0.59783623, 2.0064712,  -0.8062851,  0.18648083,   -0.18523412, 0.24748884,
            -2.3221464,  -2.8319214, -0.24122063, -0.059911821, 0.42589738,  0.71552564,
            -2.5948615
        ],
        [
            -0.51082485, 1.7834251,  -0.74631408, 0.18925786,   -0.19056884, 0.25715442,
            -2.2984084,  -2.9711354, -0.2474395,  -0.065895208, 0.4157231,   0.6854511,
            -2.6363889
        ],
        [
            -0.43078214, 1.5838689,  -0.69357863, 0.19171986,   -0.19486344, 0.26533586,
            -2.2768555,  -3.1092126, -0.25227124, -0.072648047, 0.40493192,  0.6574073,
            -2.680358
        ],
        [
            -0.35667417, 1.4041598,  -0.64687539, 0.19390686,   -0.19832913, 0.27227791,
            -2.2572442,  -3.2452979, -0.25602881, -0.080034929, 0.3937183,   0.63141278,
            -2.7258297
        ],
        [
            -0.2876813, 1.2413449,  -0.60524395, 0.19584977,   -0.20113902, 0.27819351,
            -2.2393525, -3.3786917, -0.25896128, -0.087937982, 0.38226934,  0.60741089,
            -2.7720891
        ],
        [
            -0.22314278, 1.0930079,  -0.56790936, 0.19757464,   -0.20343143, 0.28326234,
            -2.2229768,  -3.5088705, -0.26126185, -0.096009636, 0.37075763,  0.58529902,
            -2.8185971
        ],
        [
            -0.20634798, 1.0550351,  -0.55844688, 0.1980076,    -0.20397834, 0.28451019,
            -2.2187774,  -3.5435432, -0.26179668, -0.098086182, 0.36764783,  0.5796162,
            -2.8311879
        ],
        [
            -0.10535974, 0.83213338,  -0.50371883, 0.20046151, -0.20687525, 0.2914274, -2.1940959, -3.7582918,
            -0.26452228, -0.11150669, 0.34813783,  0.5462249,  -2.9108589
        ],
        [
            -0.051292521, 0.7165632,  -0.47591719, 0.20166441,  -0.20817741, 0.29474367,
            -2.181306,    -3.8772132, -0.26567819, -0.11932465, 0.33726595,  0.52898592,
            -2.956105
        ],
        [
            7.7316702e-07, 0.60928527, -0.45047894, 0.20273131,  -0.20927288, 0.29766162,
            -2.1694611,    -3.992222,  -0.26660967, -0.12669139, 0.32679895,  0.51309869,
            -3.0005405
        ],
        [
            0.095310953, 0.41584648,  -0.40554772, 0.20451987, -0.21099461, 0.302548, -2.148223, -4.2107213,
            -0.26798196, -0.14125887, 0.3072704,   0.4848855,  -3.0866107
        ],
        [
            0.18232233,  0.24562272,  -0.36705485, 0.20593609, -0.212267, 0.30647003, -2.129749, -4.4145824,
            -0.26890527, -0.15489038, 0.28974157,  0.46069706, -3.1686426
        ],
        [
            0.26236504,  0.094066794, -0.33364611, 0.20706485, -0.21323193, 0.30968406, -2.1135555, -4.6049137,
            -0.26953845, -0.16729869, 0.27417806,  0.43980403, -3.2465496
        ],
        [
            0.40546588,  -0.16585528, -0.27832414, 0.20870438, -0.21457355, 0.31463693, -2.0865858, -4.9497332,
            -0.27028875, -0.19085343, 0.24828319,  0.40567319, -3.3904778
        ],
        [
            0.53062902,  -0.38276012, -0.23413502, 0.20979292, -0.21544316, 0.31827954, -2.0651399, -5.2540825,
            -0.27066246, -0.21153324, 0.22801322,  0.37907149, -3.5200286
        ],
        [
            0.69314795,  -0.65175825, -0.18189446, 0.21081886, -0.2162757, 0.3222459, -2.0403231, -5.6504962,
            -0.27089835, -0.23841815, 0.20509821,  0.34868504, -3.6916335
        ],
        [
            0.99325255, -1.1175903, -0.098162317, 0.21189822,  -0.21726112, 0.32786709,
            -2.0032797, -6.3795404, -0.27093996,  -0.28829101, 0.17179101,  0.30310955,
            -4.0136487
        ],
        [
            1.0986131,  -1.2733004, -0.072043288, 0.21210674,  -0.21750028, 0.32943804,
            -1.9928112, -6.6334042, -0.27089518,  -0.30585398, 0.16240716,  0.28981966,
            -4.127275
        ],
        [
            1.2089611, -1.4327223, -0.046240832, 0.21226094,  -0.21771154, 0.33090425,
            -1.983136, -6.897798,  -0.27083341,  -0.32401335, 0.15363344,  0.2771661,
            -4.2462648
        ],
        [
            1.3862951,  -1.6820204, -0.0077344394, 0.21240764,  -0.21798805, 0.33293263,
            -1.9701001, -7.3192687, -0.27071911,   -0.35450634, 0.14146254,  0.25919775,
            -4.4371154
        ],
        [
            1.6094387,   -1.9854397,  0.036248881, 0.21247549, -0.21825754, 0.33501121, -1.9574977, -7.8434557,
            -0.27057268, -0.39060344, 0.12889246,  0.24005803, -4.6761341
        ],
        [
            1.7917602,   -2.2262493,  0.069046131, 0.21247389, -0.21843231, 0.3363929, -1.9497821, -8.2667449,
            -0.27046172, -0.42096222, 0.12041093,  0.22676175, -4.8702454
        ],
        [
            1.9459109, -2.425634,   0.094881396, 0.21245017, -0.21855667, 0.33737935, -1.9447182, -8.6212894,
            -0.270377, -0.44766009, 0.11424109,  0.2168739,  -5.0334728
        ],
        [
            2.0794423,   -2.5956195, 0.11601796, 0.21242009, -0.21865035, 0.33811929, -1.9412167, -8.9260731,
            -0.27031208, -0.4688904, 0.10951328, 0.20916425, -5.1742052
        ],
        [
            2.3027884,   -2.8750104,  0.14909384, 0.21236029, -0.21879051, 0.33915991, -1.9377334, -9.4313509,
            -0.27022027, -0.50537372, 0.10266209, 0.19616892, -5.4151256
        ],
        [
            2.9957468,   -3.7115619,  0.23731803,  0.21218276, -0.21906389, 0.34123927, -1.9307561, -10.968777,
            -0.27002429, -0.61901724, 0.087298823, 0.1705623,  -6.1284667
        ],
        [
            3.9120738,   -4.7708242,  0.33057863,  0.21204072, -0.21920595, 0.34249727, -1.9273197, -12.950766,
            -0.26990573, -0.76834396, 0.074799132, 0.1481395,  -7.0599244
        ],
        [
            4.6053167,   -5.5497843,  0.38925972,  0.21198837, -0.21928328, 0.34292212, -1.9274388, -14.424941,
            -0.26986995, -0.88093963, 0.068420107, 0.13614438, -7.7590981
        ],
        [
            5.2985016,   -6.315665,   0.44065331,  0.21196163, -0.21924315, 0.34313344, -1.9273634, -15.884499,
            -0.26985131, -0.99336061, 0.063550359, 0.12676384, -8.4555731
        ],
        [
            6.2148398,   -7.3137414, 0.50017919,  0.21194562, -0.21927778, 0.34326036, -1.9273568, -17.798435,
            -0.26941638, -1.1609031, 0.058570133, 0.11700937, -9.3740818
        ],
        [
            6.9079534,   -8.0605523, 0.54024811,  0.21194036, -0.21954945, 0.3433027, -1.9273567, -19.237649,
            -0.26797913, -1.2890392, 0.055541578, 0.11101872, -10.067983
        ],
        [
            7.6010142,   -8.8018916, 0.57689662,  0.21193777, -0.2193586, 0.34332388, -1.927357, -20.671267,
            -0.25831242, -1.2028803, 0.052962621, 0.10589342, -10.761462
        ],
        [
            8.5172578,   -9.7753204, 0.62103449, 0.21193625, -0.2190933, 0.34333659, -1.9273573, -22.559945,
            -0.25206987, -1.2856356, 0.05006725, 0.10012193, -11.677973
        ],
        [
            9.2104779,   -10.507707, 0.65169904,  0.21193575,  -0.21936074, 0.34334083, -1.9273575, -23.984875,
            -0.24583029, -1.3331918, 0.048176378, 0.096346524, -12.371287
        ],
        [
            10.819908,   -12.197448, 0.71561744,  0.21193534,  -0.21936091, 0.34334422, -1.9273576, -27.282758,
            -0.23357556, -1.453882,  0.044515621, 0.089030014, -13.980798
        ],
        [
            11.513118,   -12.921541, 0.74054913,  0.21193528,  -0.21936124, 0.34334464, -1.9273576, -28.699613,
            -0.22858326, -1.5033899, 0.043181147, 0.086361684, -14.674018
        ],
        [
            13.122538,   -14.595917, 0.79358812,  0.21193521,  -0.21936092, 0.34334498, -1.9273576, -31.982563,
            -0.21778548, -1.6134144, 0.040499652, 0.080999184, -16.283447
        ],
        [
            13.815546,   -15.314439, 0.81463435,  0.2119352,   -0.21936101, 0.34334502, -1.9273576, -33.393795,
            -0.21347983, -1.6588972, 0.039491158, 0.078982257, -16.976456
        ],
        [
            16.118169,   -17.693406, 0.87828372,  0.21193515, -0.21936103, 0.34334503, -1.9273578, -38.0746,
            -0.20055351, -1.8029999, 0.036617172, 0.07323434, -19.279081
        ],
        [
            18.420789,   -20.062094, 0.93413378,  0.21193508,  -0.21936112, 0.34334606, -1.9378916, -42.745342,
            -0.18948235, -1.9378916, 0.034295546, 0.068591088, -21.5817
        ]
    ],

    tail_shock_table => [
        [ 3.3461627, -0.78157,    -0.038392847, -0.038392847 ],
        [ 3.883845,  -0.84450574, -0.050816521, -0.023310252 ],
        [ 4.59045,   -0.92140593, -0.05345766,  -0.015863492 ],
        [ 5.2907041, -0.99261887, -0.053992524, -0.011463332 ],
        [ 6.2115371, -1.0803508,  -0.053470701, -0.0078157697 ],
        [ 6.9062359, -1.1429045,  -0.052652551, -0.0059835124 ],
        [ 7.6001237, -1.2027519,  -0.051666674, -0.004643598 ],
        [ 8.5168856, -1.2783263,  -0.050254755, -0.0033709947 ],
        [ 9.210286,  -1.3331658,  -0.049166913, -0.0026661248 ],
        [ 10.819867, -1.4538767,  -0.046715275, -0.0015648908 ],
        [ 11.513097, -1.5033872,  -0.04571673,  -0.0012434561 ],
        [ 13.122534, -1.613532,   -0.043560642, -0.0007135144 ],
        [ 13.815544, -1.6589912,  -0.042698578, -0.000551626 ],
        [ 16.118169, -1.8030744,  -0.040112363, -0.00018836016 ]
    ],

    blast_info => [
        0.047910341, 2.7567969,   0.25076242,  -0.65838426,  3.4694193,  0.00068138954,
        1.0285974,   -0.51472259, 0.24143548,  0.77110034,   0.61073091, -0.45960471,
        0.048659559, 0.47045686,  0.042387029, -0.043872134, 27.612,     -0.78157,
        201.33416,   -0.99440769
    ],

};
1;
