package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=5.5
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S5x5'} = {

    table_name => 'S5x5',
    symmetry   => 2,
    gamma      => 5.5,

    shock_table_info => [ 2, 5.5, 2.75e-06, 4000, 9.3e-07, 10, 1.47e-07, 2.1e-06, 5e-07, 16.14, -0.76349 ],

    shock_table => [
        [ -3.95386944,  12.00033713,   -2.99996882, -3.95515989, 0.9980631 ],
        [ -3.67731512,  11.17068499,   -2.99992852, -3.67926965, 0.99706541 ],
        [ -3.500186,    10.639311,     -2.99987839, -3.5027361,  0.99617011 ],
        [ -3.35485503,  10.20333579,   -2.99984692, -3.35802724, 0.99523438 ],
        [ -2.94662951,  8.97877955,    -2.99949035, -2.95248891, 0.99118663 ],
        [ -2.67925635,  8.17686954,    -2.99885679, -2.68801857, 0.98680404 ],
        [ -2.41076353,  7.37186061,    -2.997456,   -2.42389655, 0.98018777 ],
        [ -2.19599415,  6.72831902,    -2.9951656,  -2.21415739, 0.97255163 ],
        [ -2.0105796,   6.17327281,    -2.99159691, -2.03462244, 0.9636036 ],
        [ -1.85024527,  5.6939962,     -2.98647139, -1.88089865, 0.95352419 ],
        [ -1.70401858,  5.25779068,    -2.9791704,  -1.74228591, 0.94190449 ],
        [ -1.57047402,  4.8605618,     -2.96920926, -1.61734657, 0.92877903 ],
        [ -1.4452968,   4.48967686,    -2.95577016, -1.50199075, 0.91383627 ],
        [ -1.32650492,  4.13955713,    -2.93795497, -1.39441354, 0.89686228 ],
        [ -1.21364389,  3.80921051,    -2.91497048, -1.29423847, 0.8778278 ],
        [ -1.10186048,  3.48497611,    -2.88477387, -1.19731054, 0.85583462 ],
        [ -0.98765586,  3.15772082,    -2.84453272, -1.1010195,  0.82984212 ],
        [ -0.86344314,  2.80775496,    -2.788063,   -0.99991588, 0.79732486 ],
        [ -0.74498384,  2.48134568,    -2.72057231, -0.90750599, 0.76222549 ],
        [ -0.63649277,  2.1900842,     -2.64690198, -0.8267034,  0.72687275 ],
        [ -0.53482827,  1.92490624,    -2.56841828, -0.75459106, 0.69143724 ],
        [ -0.43849607,  1.68135379,    -2.48712625, -0.68966342, 0.65636665 ],
        [ -0.34570238,  1.45438023,    -2.40430881, -0.63035679, 0.62180339 ],
        [ -0.25277117,  1.23489935,    -2.31897939, -0.5741889,  0.58702437 ],
        [ -0.15983996,  1.02337973,    -2.23330094, -0.52123816, 0.55265816 ],
        [ -0.06601303,  0.81784131,    -2.14831815, -0.47097683, 0.51890045 ],
        [ 0.02955881,   0.61653887,    -2.06491318, -0.42297314, 0.48592468 ],
        [ 0.12786435,   0.41758102,    -1.98368097, -0.37679776, 0.45383164 ],
        [ 0.22975278,   0.21951641,    -1.9052084,  -0.33216022, 0.42275642 ],
        [ 0.3360374,    0.02107992,    -1.82999585, -0.2888407,  0.39283166 ],
        [ 0.44755485,   -0.17893688,   -1.75842699, -0.24665698, 0.36416987 ],
        [ 0.56533746,   -0.38198333,   -1.69069395, -0.20540304, 0.33682695 ],
        [ 0.69051468,   -0.58954365,   -1.6269135,  -0.16489796, 0.31084224 ],
        [ 0.82414656,   -0.80287197,   -1.56722245, -0.1250356,  0.28627295 ],
        [ 0.96768609,   -1.02373919,   -1.51157889, -0.08564394, 0.26310911 ],
        [ 1.11764636,   -1.24656626,   -1.4614991,  -0.04780703, 0.2420046 ],
        [ 1.27345233,   -1.47070803,   -1.41683795, -0.01162146, 0.22293761 ],
        [ 1.43645479,   -1.69831093,   -1.3768378,  0.02327334,  0.20562403 ],
        [ 1.60776171,   -1.93101665,   -1.34094224, 0.05711533,  0.18985813 ],
        [ 1.78835844,   -2.17019734,   -1.30870721, 0.0900731,   0.17547765 ],
        [ 1.97972364,   -2.41778511,   -1.27968026, 0.12236282,  0.16231151 ],
        [ 2.18356925,   -2.67589701,   -1.25347284, 0.15418556,  0.15021166 ],
        [ 2.40116197,   -2.9459981,    -1.22982376, 0.18562978,  0.13908498 ],
        [ 2.63342448,   -3.22909681,   -1.20853136, 0.21671733,  0.12886479 ],
        [ 2.88229961,   -3.52741133,   -1.18932582, 0.24758741,  0.11944993 ],
        [ 3.14926816,   -3.84254604,   -1.17201652, 0.27828941,  0.110775 ],
        [ 3.43660988,   -4.17700424,   -1.15639102, 0.30893896,  0.10276136 ],
        [ 3.74619393,   -4.53275751,   -1.14229347, 0.33957666,  0.09535664 ],
        [ 4.0806857,    -4.91265124,   -1.12955582, 0.37029636,  0.08849985 ],
        [ 4.442542,     -5.31924305,   -1.11804584, 0.40114138,  0.0821466 ],
        [ 4.83501377,   -5.75593838,   -1.10762814, 0.43219409,  0.07624831 ],
        [ 5.26634716,   -6.23157542,   -1.09809123, 0.46385619,  0.07070806 ],
        [ 5.72316855,   -6.7312256,    -1.08966897, 0.49498145,  0.06569008 ],
        [ 6.22791833,   -7.27921654,   -1.081904,   0.52691,     0.06094679 ],
        [ 6.77901585,   -7.87344546,   -1.0748399,  0.55924687,  0.05652425 ],
        [ 7.38067859,   -8.5181446,    -1.06841344, 0.59198302,  0.05240371 ],
        [ 8.03891941,   -9.21943086,   -1.06255372, 0.6251783,   0.04855909 ],
        [ 8.75954708,   -9.98315127,   -1.05720615, 0.6588444,   0.04497226 ],
        [ 9.54856732,   -10.81532398,  -1.05232396, 0.69297372,  0.04162811 ],
        [ 10.41318366,  -11.72319309,  -1.04786162, 0.72758175,  0.03851009 ],
        [ 11.3609981,   -12.71437906,  -1.0437805,  0.76266771,  0.03560425 ],
        [ 12.39981171,  -13.79667528,  -1.04004788, 0.79821044,  0.03289894 ],
        [ 13.53882498,  -14.97929968,  -1.03663196, 0.83420972,  0.0303814 ],
        [ 14.78663815,  -16.27081497,  -1.0335079,  0.87062055,  0.0280424 ],
        [ 16.1542513,   -17.68224029,  -1.03064897, 0.90744509,  0.02586992 ],
        [ 17.65246447,  -19.2243541,   -1.02803353, 0.94465204,  0.02385451 ],
        [ 19.29347766,  -20.90934862,  -1.02564088, 0.98222224,  0.02198635 ],
        [ 21.09129087,  -22.75123212,  -1.02345133, 1.02014973,  0.02025545 ],
        [ 23.02570387,  -24.72904715,  -1.02148032, 1.05776613,  0.01867911 ],
        [ 25.30610418,  -27.05615599,  -1.01954697, 1.09851511,  0.01711488 ],
        [ 27.16021725,  -28.94524551,  -1.01821415, 1.12921724,  0.01602691 ],
        [ 29.4945692,   -31.32039017,  -1.01677557, 1.16521101,  0.01484233 ],
        [ 31.82558489,  -33.68904586,  -1.01554998, 1.19859427,  0.01382474 ],
        [ 34.43742356,  -36.33992413,  -1.01437466, 1.23338512,  0.01284082 ],
        [ 37.37645116,  -39.31949964,  -1.01324897, 1.26969221,  0.01189091 ],
        [ 40.6989698,   -42.6841965,   -1.01217169, 1.30763796,  0.01097533 ],
        [ 43.97454916,  -45.99812691,  -1.0112684,  1.34229135,  0.01020263 ],
        [ 47.09991115,  -49.15751653,  -1.01052421, 1.37315502,  0.0095616 ],
        [ 50.56272958,  -52.65550521,  -1.00980632, 1.40516492,  0.00894035 ],
        [ 54.41288516,  -56.54205377,  -1.00911552, 1.43840169,  0.00833901 ],
        [ 58.70979411,  -60.87666739,  -1.00845122, 1.47295423,  0.00775771 ],
        [ 63.52467859,  -65.73066732,  -1.00781331, 1.50892161,  0.00719657 ],
        [ 68.94349267,  -71.1901185,   -1.00720189, 1.54641462,  0.00665571 ],
        [ 75.0707301,   -77.35963925,  -1.00661652, 1.58555735,  0.00613524 ],
        [ 80.80944034,  -83.13494402,  -1.00614854, 1.61953771,  0.00571717 ],
        [ 87.21824922,  -89.58167923,  -1.00569863, 1.65485238,  0.0053134 ],
        [ 94.40490419,  -96.80769415,  -1.00526641, 1.6916029,   0.004924 ],
        [ 102.49965296, -104.94335049, -1.00485207, 1.7299034,   0.00454903 ],
        [ 111.66111854, -114.14739953, -1.0044554,  1.76988187,  0.00418854 ],
        [ 119.88773473, -122.40938593, -1.00415073, 1.80316999,  0.00391061 ],
        [ 131.4908258,  -134.05845558, -1.0037857,  1.84654749,  0.00357633 ],
        [ 142.00866217, -144.61460174, -1.00350624, 1.88277803,  0.00331944 ],
        [ 153.81867206, -156.46439284, -1.00323777, 1.9204813,   0.00307195 ],
        [ 167.14036363, -169.82745609, -1.00298046, 1.95977619,  0.00283388 ],
        [ 182.24181376, -184.9720025,  -1.0027342,  2.00079672,  0.00260525 ],
        [ 194.93165485, -197.6953882,  -1.00255665, 2.03278552,  0.00244 ],
        [ 208.97448816, -211.77289292, -1.00238526, 2.06590143,  0.00228008 ]
    ],

    energy_table => [
        [ -3.95386944,  4.269e-05,  0.0001012,  4.269e-05,  0.0001012,     5.0480621,   -1.5309912 ],
        [ -3.67731512,  8.206e-05,  0.0001932,  8.206e-05,  0.0001932,     4.6238746,   -1.5365279 ],
        [ -3.500186,    0.00012439, 0.00029134, 0.00012439, 0.00029134,    4.3514048,   -1.5398627 ],
        [ -3.35485503,  0.00017468, 0.00040717, 0.00017468, 0.00040717,    4.1274234,   -1.5431597 ],
        [ -2.94662951,  0.00044891, 0.00102845, 0.00044891, 0.00102845,    3.4951978,   -1.5542948 ],
        [ -2.67925635,  0.00082461, 0.00186054, 0.00082461, 0.00186054,    3.0786443,   -1.5652223 ],
        [ -2.41076353,  0.00150264, 0.00332345, 0.00150264, 0.00332345,    2.6564322,   -1.5798927 ],
        [ -2.19599415,  0.00240508, 0.00521133, 0.00240508, 0.00521133,    2.3158541,   -1.5940096 ],
        [ -2.0105796,   0.00357896, 0.00758471, 0.00357896, 0.00758471,    2.018985,    -1.6106572 ],
        [ -1.85024527,  0.0050079,  0.01036568, 0.0050079,  0.01036568,    1.7594177,   -1.6281036 ],
        [ -1.70401858,  0.0067524,  0.01361701, 0.0067524,  0.01361701,    1.5201198,   -1.6469407 ],
        [ -1.57047402,  0.00880665, 0.01725925, 0.00880665, 0.01725925,    1.2989047,   -1.6669678 ],
        [ -1.4452968,   0.01121291, 0.02128473, 0.01121291, 0.02128473,    1.0890088,   -1.6875805 ],
        [ -1.32650492,  0.01399463, 0.02562883, 0.01399463, 0.02562883,    0.88732152,  -1.7084656 ],
        [ -1.21364389,  0.01713931, 0.0301529,  0.01713931, 0.0301529,     0.6933609,   -1.728909 ],
        [ -1.10186048,  0.02077267, 0.03488042, 0.02077267, 0.03488042,    0.49895389,  -1.7489813 ],
        [ -0.98765586,  0.02503533, 0.03975544, 0.02503533, 0.03975544,    0.29806389,  -1.7684276 ],
        [ -0.86344314,  0.03029062, 0.04478219, 0.03029062, 0.04478219,    0.077133868, -1.7877976 ],
        [ -0.74498384,  0.03584954, 0.0489353,  0.03584954, 0.0489353,     -0.13568167, -1.804265 ],
        [ -0.63649277,  0.04132695, 0.05188052, 0.04132695, 0.05188052,    -0.332197,   -1.8178471 ],
        [ -0.53482827,  0.04670192, 0.05369464, 0.04670192, 0.05369464,    -0.51762645, -1.8298743 ],
        [ -0.43849607,  0.05191948, 0.05447412, 0.05191948, 0.05447412,    -0.6944444,  -1.8413123 ],
        [ -0.34570238,  0.05697512, 0.05435118, 0.05697512, 0.05435118,    -0.86582591, -1.8528729 ],
        [ -0.25277117,  0.06198894, 0.05342666, 0.06198894, 0.05342666,    -1.0385707,  -1.8652177 ],
        [ -0.15983996,  0.06688365, 0.05180957, 0.06688365, 0.05180957,    -1.2124998,  -1.8782679 ],
        [ -0.06601303,  0.07164577, 0.0496188,  0.07164577, 0.0496188,     -1.3893645,  -1.8918637 ],
        [ 0.02955881,   0.07626422, 0.04697401, 0.07626422, 0.04697401,    -1.5708406,  -1.9056913 ],
        [ 0.12786435,   0.08073674, 0.04398637, 0.08073674, 0.04398637,    -1.7588729,  -1.9194013 ],
        [ 0.22975278,   0.08505477, 0.04076439, 0.08505477, 0.04076439,    -1.9551418,  -1.9326248 ],
        [ 0.3360374,    0.08920855, 0.0374103,  0.08920855, 0.0374103,     -2.16125,    -1.9450097 ],
        [ 0.44755485,   0.09318959, 0.03401631, 0.09318959, 0.03401631,    -2.3788311,  -1.9562825 ],
        [ 0.56533746,   0.09699577, 0.03065835, 0.09699577, 0.03065835,    -2.6098926,  -1.966285 ],
        [ 0.69051468,   0.10062604, 0.02739952, 0.10062604, 0.02739952,    -2.8566284,  -1.9748671 ],
        [ 0.82414656,   0.10407572, 0.02429448, 0.10407572, 0.02429448,    -3.1210721,  -1.982081 ],
        [ 0.96768609,   0.10734865, 0.02137884, 0.10734865, 0.02137884,    -3.40607,    -1.9879606 ],
        [ 1.11764636,   0.1103531,  0.01875983, 0.1103531,  0.01875983,    -3.7045705,  -1.9923672 ],
        [ 1.27345233,   0.11309075, 0.01644619, 0.11309075, 0.01644619,    -4.0152906,  -1.9956 ],
        [ 1.43645479,   0.11560037, 0.0144056,  0.11560037, 0.0144056,     -4.3408041,  -1.9979017 ],
        [ 1.60776171,   0.11790976, 0.01261068, 0.11790976, 0.01261068,    -4.6832235,  -1.9994573 ],
        [ 1.78835844,   0.12004064, 0.01103651, 0.12004064, 0.01103651,    -5.0444317,  -2.0004416 ],
        [ 1.97972364,   0.12201635, 0.00965614, 0.12201635, 0.00965614,    -5.4273195,  -2.0010121 ],
        [ 2.18356925,   0.12385725, 0.00844527, 0.12385725, 0.00844527,    -5.8352582,  -2.0012762 ],
        [ 2.40116197,   0.12557574, 0.00738553, 0.12557574, 0.00738553,    -6.2707358,  -2.0013216 ],
        [ 2.63342448,   0.12718015, 0.00646135, 0.12718015, 0.00646135,    -6.7355625,  -2.0012661 ],
        [ 2.88229961,   0.12868438, 0.00565463, 0.12868438, 0.00565463,    -7.2336191,  -2.0011432 ],
        [ 3.14926816,   0.13009687, 0.00495145, 0.13009687, 0.00495145,    -7.7678363,  -2.0009493 ],
        [ 3.43660988,   0.13142836, 0.00433775, 0.13142836, 0.00433775,    -8.3427616,  -2.0007625 ],
        [ 3.74619393,   0.13268548, 0.00380256, 0.13268548, 0.00380256,    -8.9621394,  -2.0005975 ],
        [ 4.0806857,    0.13387646, 0.0033352,  0.13387646, 0.0033352,     -9.631294,   -2.0004402 ],
        [ 4.442542,     0.13500684, 0.00292705, 0.13500684, 0.00292705,    -10.355138,  -2.0003104 ],
        [ 4.83501377,   0.13608305, 0.00257003, 0.13613775, 0.0028629377,  -11.140181,  -2.0002123 ],
        [ 5.26634716,   0.137121,   0.00225421, 0.13736353, 0.0028580304,  -12.00292,   -2.0001394 ],
        [ 5.72316855,   0.13808702, 0.0019847,  0.13872833, 0.0030454585,  -12.916613,  -2.0000885 ],
        [ 6.22791833,   0.13902601, 0.00174454, 0.14018611, 0.0026241186,  -13.926145,  -2.0000525 ],
        [ 6.77901585,   0.13992717, 0.00153341, 0.14154268, 0.0023243804,  -15.028362,  -2.0000308 ],
        [ 7.38067859,   0.14079205, 0.00134704, 0.1428683,  0.002086088,   -16.2317,    -2.0000168 ],
        [ 8.03891941,   0.14162388, 0.00118518, 0.14417317, 0.0018818134,  -17.548189,  -2.0000087 ],
        [ 8.75954708,   0.14242428, 0.00104145, 0.14546194, 0.0017005008,  -18.989449,  -2.0000042 ],
        [ 9.54856732,   0.1431943,  0.00091499, 0.14673723, 0.0015368554,  -20.567491,  -2.000002 ],
        [ 10.41318366,  0.14393553, 0.00080365, 0.14799904, 0.0013867191,  -22.296725,  -2.0000009 ],
        [ 11.3609981,   0.14464907, 0.0007056,  0.14924636, 0.0012487897,  -24.192355,  -2.0000004 ],
        [ 12.39981171,  0.14533559, 0.0006193,  0.15047388, 0.0011187902,  -26.269982,  -2.0000002 ],
        [ 13.53882498,  0.14599615, 0.00054335, 0.15168175, 0.0010051084,  -28.548009,  -2.0000002 ],
        [ 14.78663815,  0.14663097, 0.00047659, 0.15294522, 0.00096747224, -31.043635,  -2.0000002 ],
        [ 16.1542513,   0.14724115, 0.00041789, 0.1541829,  0.00084699297, -33.778862,  -2.0000001 ],
        [ 17.65246447,  0.14782719, 0.00036632, 0.15536986, 0.00074143981, -36.775288,  -2.0000001 ],
        [ 19.29347766,  0.14838981, 0.00032104, 0.15650788, 0.00064897832, -40.057315,  -2 ],
        [ 21.09129087,  0.14892992, 0.00028128, 0.15759909, 0.00056796982, -43.652941,  -2 ],
        [ 23.02570387,  0.14943962, 0.00024694, 0.15862777, 0.00049813753, -47.521767,  -2 ],
        [ 25.30610418,  0.1499645,  0.00021467, 0.15968605, 0.00043261613, -52.082568,  -2 ],
        [ 27.16021725,  0.15034209, 0.00019328, 0.16044677, 0.00038926207, -55.790794,  -2 ],
        [ 29.4945692,   0.1507664,  0.00017102, 0.16130103, 0.00034418697, -60.459498,  -2 ],
        [ 31.82558489,  0.15114308, 0.00015276, 0.16205894, 0.00030725582, -65.121529,  -2 ],
        [ 34.43742356,  0.15151927, 0.00013588, 0.16281547, 0.00027315754, -70.345207,  -2 ],
        [ 37.37645116,  0.15189495, 0.00012032, 0.16357058, 0.00024176197, -76.223262,  -2 ],
        [ 40.6989698,   0.15227009, 0.00010602, 0.16432422, 0.00021293105, -82.868299,  -2 ],
        [ 43.97454916,  0.15259786, 9.449e-05,  0.16498243, 0.00018971196, -89.419458,  -2 ],
        [ 47.09991115,  0.15287843, 8.533e-05,  0.16554568, 0.00017125527, -95.670182,  -2 ],
        [ 50.56272958,  0.15315864, 7.678e-05,  0.16610804, 0.00015405871, -102.59582,  -2 ],
        [ 54.41288516,  0.15343849, 6.884e-05,  0.1666695,  0.00013808814, -110.29613,  -2 ],
        [ 58.70979411,  0.15371795, 6.148e-05,  0.16723004, 0.0001232876,  -118.88995,  -2 ],
        [ 63.52467859,  0.15399701, 5.467e-05,  0.16778965, 0.00010960931, -128.51972,  -2 ],
        [ 68.94349267,  0.15427566, 4.84e-05,   0.1683483,  9.7008634e-05, -139.35734,  -2 ],
        [ 75.0707301,   0.15455388, 4.263e-05,  0.16890597, 8.5433878e-05, -151.61182,  -2 ],
        [ 80.80944034,  0.1547854,  3.82e-05,   0.16936995, 7.6538455e-05, -163.08924,  -2 ],
        [ 87.21824922,  0.15501661, 3.409e-05,  0.16983322, 6.8297838e-05, -175.90686,  -2 ],
        [ 94.40490419,  0.15524749, 3.03e-05,   0.17029579, 6.0681646e-05, -190.28017,  -2 ],
        [ 102.49965296, 0.15547805, 2.68e-05,   0.17075763, 5.3666567e-05, -206.46967,  -2 ],
        [ 111.66111854, 0.15570827, 2.358e-05,  0.17121874, 4.7224753e-05, -224.7926,   -2 ],
        [ 119.88773473, 0.15589219, 2.121e-05,  0.17158709, 4.2466492e-05, -241.24583,  -2 ],
        [ 131.4908258,  0.15612179, 1.848e-05,  0.17204685, 3.6991322e-05, -264.45201,  -2 ],
        [ 142.00866217, 0.15630521, 1.647e-05,  0.17241411, 3.2972932e-05, -285.48768,  -2 ],
        [ 153.81867206, 0.15648839, 1.462e-05,  0.17278087, 2.9260731e-05, -309.1077,   -2 ],
        [ 167.14036363, 0.15667134, 1.291e-05,  0.17314713, 2.5843717e-05, -335.75109,  -2 ],
        [ 182.24181376, 0.15685404, 1.135e-05,  0.17351288, 2.2708345e-05, -365.95399,  -2 ],
        [ 194.93165485, 0.15699091, 1.026e-05,  0.17378686, 2.0533276e-05, -391.33367,  -2 ],
        [ 208.97448816, 0.15712764, 9.25e-06,   0.17406055, 1.8503977e-05, -419.41934,  -2 ]
    ],

    impulse_table => [
        [
            -4.558706,  13.814838,  -4.5592267,    0.043031401, -0.0059895052, -0.024646076,
            -2.7524332, -0.3341855, -0.0076297267, -0.21559899, 0.47934624,    0.99994294,
            -4.2479606
        ],
        [
            -4.1997043, 12.737836,   -4.2005966,   0.05102953,  -0.008575791, -0.020169199,
            -2.7479143, -0.51371223, -0.010924579, -0.21117947, 0.47934543,   0.99986757,
            -4.0049736
        ],
        [
            -3.9120222, 11.874797,   -3.9133963,   0.058420855, -0.011432949, -0.015257989,
            -2.7429245, -0.65760296, -0.014565056, -0.20617946, 0.47934395,   0.99973901,
            -3.8148732
        ],
        [
            -3.5065571, 10.658424,  -3.5090829,   0.070516969, -0.017141607, -0.0055787995,
            -2.7329618, -0.8605251, -0.021841864, -0.19649435, 0.47933794,   0.99932271,
            -3.5553169
        ],
        [
            -2.9957315, 9.1260621, -3.0011735,   0.088912374, -0.028513046, 0.013048653,
            -2.7131256, -1.116814, -0.036361358, -0.17779,    0.47930716,   0.99777441,
            -3.2460635
        ],
        [
            -2.3025843, 7.0476475,  -2.318046,    0.12001583,  -0.05629074, 0.054501767,
            -2.6643038, -1.4698765, -0.072137562, -0.13669265, 0.47903295,  0.98920682,
            -2.8735262
        ],
        [
            -1.8971192, 5.8340255,  -1.92567,    0.14104125,  -0.082352735, 0.089411478,
            -2.6170198, -1.6873832, -0.10631169, -0.10549474, 0.47829244,   0.9737407,
            -2.6964945
        ],
        [
            -1.6094371,  4.9763167,   -1.6536154, 0.15642881, -0.10602637, 0.11957885, -2.5717303, -1.8560303,
            -0.13781411, -0.08402039, 0.47686583, 0.95204824, -2.6004609
        ],
        [
            -1.3862936, 4.31551,    -1.4483058,  0.16805995,   -0.12699168, 0.14608153,
            -2.5288027, -2.0033729, -0.16584923, -0.070501629, 0.47455702,  0.9253653,
            -2.5496663
        ],
        [
            -1.203972,  3.7810283,   -1.2857568, 0.177052,   -0.14516919, 0.16943186, -2.4885122, -2.1415235,
            -0.1900197, -0.06326062, 0.47120788, 0.89512042, -2.5278056
        ],
        [
            -1.0498213, 3.3352893,   -1.1530666, 0.18416129, -0.16064192, 0.18991469, -2.4510375, -2.2764599,
            -0.2102975, -0.06063053, 0.46671156, 0.86270834, -2.5258315
        ],
        [
            -0.91628994, 2.9557867,  -1.0424355,  0.18991964,   -0.17360433, 0.20774914,
            -2.4164543,  -2.4110919, -0.22692751, -0.061845483, 0.46102152,  0.82935793,
            -2.5380694
        ],
        [
            -0.79850691, 2.6278328,  -0.94874342, 0.19469898,   -0.18432145, 0.22315376,
            -2.3847386,  -2.5465809, -0.24032003, -0.065566664, 0.45415395,  0.79606734,
            -2.5606506
        ],
        [
            -0.69314639, 2.3411838, -0.86842,    0.19875354,   -0.19309311, 0.23636594,
            -2.3557801,  -2.683047, -0.25095798, -0.071117759, 0.44618293,  0.76358776,
            -2.590783
        ],
        [
            -0.59783621, 2.0883177,  -0.7988601,  0.20225316,   -0.20022396, 0.24763915,
            -2.3294073,  -2.8200157, -0.25932904, -0.078047168, 0.43722992,  0.73243983,
            -2.6263775
        ],
        [
            -0.51082483, 1.8634912,  -0.73809764, 0.20531023,   -0.20600096, 0.25723068,
            -2.3054085,  -2.9567267, -0.26588149, -0.086159473, 0.42745005,  0.70294839,
            -2.6658389
        ],
        [
            -0.43078213, 1.6621943,  -0.68461126, 0.20800038,   -0.21067938, 0.26538777,
            -2.2835564,  -3.0923448, -0.27100212, -0.094589226, 0.41701794,  0.67528411,
            -2.7079391
        ],
        [
            -0.35667415, 1.4808145,   -0.63720163, 0.21037699, -0.21447668, 0.272337, -2.2636264, -3.2260924,
            -0.27500929, -0.10335276, 0.40611504,  0.64950391, -2.7517319
        ],
        [
            -0.33293751, 1.4237638,   -0.62245009, 0.21111309, -0.21558809, 0.27443831, -2.257319, -3.2705316,
            -0.27615281, -0.10662084, 0.40237415,  0.64125554, -2.7667068
        ],
        [
            -0.22314276, 1.1665975,   -0.55695989, 0.21434329, -0.22010976, 0.2833811, -2.2287089, -3.4855239,
            -0.28064154, -0.12168103, 0.38359867,  0.60345756, -2.8416652
        ],
        [
            -0.16251814, 1.0293642,   -0.52271959, 0.21599362, -0.22220435, 0.28778953, -2.2133605, -3.6103564,
            -0.2826158,  -0.13048672, 0.37230446,  0.58301597, -2.8868364
        ],
        [
            -0.10535972, 0.90306461, -0.49166886, 0.21745602, -0.22394575, 0.2916211, -2.1992131, -3.7315937,
            -0.28419497, -0.139353,  0.36116828,  0.56414324, -2.931695
        ],
        [
            -0.051292504, 0.78631345, -0.46337654, 0.21875246,  -0.22540435, 0.29497238,
            -2.1861378,   -3.849119,  -0.28546673, -0.14761382, 0.35029791,  0.54671627,
            -2.9760135
        ],
        [
            7.9085384e-07, 0.67794785, -0.43748392, 0.2199026,   -0.22663523, 0.29792205,
            -2.174021,     -3.9628979, -0.28649793, -0.15600088, 0.33977688,  0.5306133,
            -3.0196285
        ],
        [
            0.095310971, 0.48258367,  -0.39174033, 0.22183206, -0.22857725, 0.3028608, -2.1522906, -4.1793673,
            -0.28803063, -0.17195384, 0.32000441,  0.50192072, -3.1043294
        ],
        [
            0.18232235,  0.31072436,  -0.35254604, 0.22336179, -0.23001883, 0.30682139, -2.1333824, -4.3816652,
            -0.28907531, -0.18660793, 0.3021056,   0.47722828, -3.1852847
        ],
        [
            0.26236506,  0.15777323,  -0.318528,  0.22458266, -0.23111604, 0.31006258, -2.1168062, -4.5707909,
            -0.28980189, -0.20073389, 0.28610445, 0.45583506, -3.262341
        ],
        [
            0.4054659,   -0.10438018, -0.26220393, 0.22635916, -0.23264802, 0.31504554, -2.0891986, -4.9139513,
            -0.29068418, -0.22573611, 0.25928718,  0.42076765, -3.4050518
        ],
        [
            0.53062904, -0.32297171, -0.2172276, 0.22754097, -0.2336456, 0.31869835, -2.0672473, -5.2172902,
            -0.2911447, -0.24820683, 0.23815986, 0.39334249, -3.5338154
        ],
        [
            0.69314797,  -0.59382613, -0.1640801, 0.22865655, -0.23460471, 0.3226618, -2.0418452, -5.6128664,
            -0.29146419, -0.27702422, 0.21417269, 0.36192715, -3.7047098
        ],
        [
            0.99325256,  -1.0622685,  -0.07896599, 0.22983137, -0.23574616, 0.32825149, -2.0039233, -6.3413081,
            -0.29160656, -0.33145397, 0.17919008,  0.31466158, -4.0260717
        ],
        [
            1.0986131,  -1.2186928, -0.052437045, 0.23005805,  -0.23602428, 0.32980783,
            -1.9932034, -6.5951457, -0.29158663,  -0.35050119, 0.16932234,  0.30085208,
            -4.1396063
        ],
        [
            1.2089611,  -1.3787716, -0.026240504, 0.23022542,  -0.23627042, 0.33125849,
            -1.9832923, -6.859577,  -0.29154701,  -0.37021292, 0.16009538,  0.2876949,
            -4.2585515
        ],
        [
            1.3862952,   -1.6289598,  0.012833258, 0.23038363, -0.23659285, 0.33326194, -1.9699341, -7.2812038,
            -0.29146262, -0.40108059, 0.14729806,  0.26899923, -4.4494119
        ],
        [
            1.6094387,   -1.9332651,  0.057433601, 0.23045507, -0.2369075, 0.33531198, -1.957015, -7.8056855,
            -0.29134426, -0.44034926, 0.1340888,   0.24907328, -4.6885335
        ],
        [
            1.7917603,   -2.1746484,  0.090669616, 0.23045108, -0.23711168, 0.33667336, -1.9491011, -8.2292524,
            -0.29125274, -0.47389273, 0.12518333,  0.2352267,  -4.8827703
        ],
        [
            1.9459109,   -2.3744352,  0.11683785, 0.23042324, -0.23725687, 0.3376447, -1.9439031, -8.5840387,
            -0.29118127, -0.49954659, 0.11871008, 0.22492882, -5.0461173
        ],
        [
            2.0794423,   -2.5447181,  0.13823843, 0.23038869, -0.23736592, 0.33837297, -1.9403084, -8.8890298,
            -0.29112687, -0.52333189, 0.11375314, 0.21689967, -5.1869577
        ],
        [
            2.3027886,   -2.8245213,  0.17171303, 0.23032084, -0.23752882, 0.3393974, -1.936665, -9.3946397,
            -0.29105008, -0.56312626, 0.1065759,  0.2034492,  -5.4277408
        ],
        [
            2.9958113,   -3.6619699, 0.26092561,  0.23012085, -0.23779357, 0.34144397, -1.9238166, -10.933054,
            -0.29088492, -0.6839731, 0.090511611, 0.176747,   -6.1418059
        ],
        [
            3.9120778,   -4.7216901,  0.35509815,  0.2299614,  -0.23794143, 0.34268224, -1.9172562, -12.915615,
            -0.29078616, -0.84255834, 0.077483122, 0.15341959, -7.0736938
        ],
        [
            4.6053624,   -5.5009088,  0.41430609,  0.22990269, -0.23812963, 0.34309949, -1.9173244, -14.39019,
            -0.29075757, -0.96208873, 0.070849002, 0.14095989, -7.7731034
        ],
        [
            5.2983664,   -6.2667252, 0.46611416,  0.2298727,  -0.23832117, 0.34330785, -1.9172742, -15.849566,
            -0.29071897, -1.0774637, 0.065793108, 0.13122859, -8.469504
        ],
        [
            6.2149154,   -7.2651474, 0.52611678,  0.22985473, -0.23860558, 0.34343298, -1.9172677, -17.7641,
            -0.29069952, -1.2326955, 0.060625621, 0.12111221, -9.3882953
        ],
        [
            6.9078323,   -8.011807,  0.56646743,  0.22984884, -0.23795713, 0.34347471, -1.9172676, -19.202981,
            -0.28707644, -1.2795086, 0.057487111, 0.11490574, -10.082026
        ],
        [
            7.6010947,   -8.7534086, 0.60338297,  0.22984593, -0.23885993, 0.3434956, -1.9172679, -20.637069,
            -0.28758883, -1.4041905, 0.054814049, 0.10959427, -10.775721
        ],
        [
            8.5173392,   -9.7268853, 0.64781485, 0.22984421, -0.23812944, 0.34350813, -1.9172683, -22.525798,
            -0.26533578, -1.3351002, 0.05181486, 0.10361635, -11.692242
        ],
        [
            9.2103596,   -10.459089, 0.67866631,  0.22984364,  -0.23813012, 0.34351231, -1.9172684, -23.950346,
            -0.25931109, -1.391379,  0.049857118, 0.099707603, -12.38536
        ],
        [
            10.81979,    -12.148882, 0.74297204,  0.22984316,  -0.23813069, 0.34351565, -1.9172685, -27.248277,
            -0.24597369, -1.5154931, 0.046066408, 0.092131509, -13.994874
        ],
        [
            11.513,      -12.872991, 0.76804733,  0.22984309,  -0.23813038, 0.34351607, -1.9172685, -28.665149,
            -0.24060011, -1.5664528, 0.044684738, 0.089368827, -14.688095
        ],
        [
            13.12242,    -14.547399, 0.82137906,  0.229843,    -0.23813035, 0.3435164, -1.9172686, -31.948128,
            -0.22905172, -1.6797919, 0.041908634, 0.083817141, -16.297524
        ],
        [
            13.815628,   -15.266139, 0.84254252,  0.22984297,  -0.23813042, 0.34351645, -1.9172686, -33.359778,
            -0.22446571, -1.7266841, 0.040864339, 0.081728617, -16.990733
        ],
        [
            16.118251,   -17.645135, 0.90651282,  0.22984289,  -0.23813023, 0.34351648, -1.9172686, -38.04061,
            -0.21074966, -1.8752864, 0.037889405, 0.075778805, -19.293358
        ],
        [
            18.420871,  -20.013845, 0.96262667,  0.22984281,  -0.2381035, 0.34351686, -2.0144764, -42.711371,
            -0.1990435, -2.0144764, 0.035486412, 0.070972825, -21.595978
        ]
    ],

    tail_shock_table => [
        [ 2.8275201, -0.76349,    -0.038793378, -0.038793378 ],
        [ 2.928708,  -0.77680937, -0.045025039, -0.031904097 ],
        [ 3.8831381, -0.89253364, -0.052028648, -0.016962296 ],
        [ 4.5901164, -0.9703998,  -0.052880668, -0.01198502 ],
        [ 5.2903659, -1.0427656,  -0.052701055, -0.0087847783 ],
        [ 6.211526,  -1.1322624,  -0.051748877, -0.006027805 ],
        [ 6.9060688, -1.1962058,  -0.050774443, -0.0046143374 ],
        [ 7.6001803, -1.2575216,  -0.049704576, -0.0035691158 ],
        [ 8.5169568, -1.335045,   -0.048244625, -0.0025693995 ],
        [ 9.2101624, -1.3913512,  -0.047148337, -0.0020131117 ],
        [ 10.819748, -1.5154874,  -0.044722651, -0.0011413463 ],
        [ 11.512979, -1.5664499,  -0.043745554, -0.00088674199 ],
        [ 13.122416, -1.6799172,  -0.041649163, -0.00046769802 ],
        [ 13.815626, -1.7267797,  -0.040814346, -0.00034008362 ],
        [ 16.118251, -1.8753957,  -0.038320341, -5.5481352e-05 ]
    ],

    blast_info => [
        0.035092025, 2.7629053,   0.22607424,  -0.72835023,  3.4874457,  0.00075333274,
        0.71583369,  -0.57178185, 0.18014581,  0.85189596,   0.60211877, -0.50255111,
        0.033166143, 0.47934738,  0.041789615, -0.043296352, 16.14,      -0.76349,
        114.21911,   -0.98491737
    ],

};
1;
