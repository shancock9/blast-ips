package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=1.1
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S1x1'} = {

    table_name => 'S1x1',
    symmetry   => 2,
    gamma      => 1.1,

    shock_table_info => [ 2, 1.1, 1.15e-06, 8000, 2.03e-07, 10, 3e-07, 3.55e-07, 5e-07, 1590.9, -0.431 ],

    shock_table => [
        [ -5.03701586,  12.00017661,   -2.99998069, -5.03803112, 0.99847635 ],
        [ -4.72942015,  11.0773991,    -2.99995142, -4.73103111, 0.99758165 ],
        [ -4.72572497,  11.06631376,   -2.99995088, -4.7273449,  0.99756818 ],
        [ -4.71833461,  11.04414281,   -2.99994978, -4.71997261, 0.99754103 ],
        [ -4.71654301,  11.03876834,   -2.99994951, -4.71818542, 0.9975344 ],
        [ -4.71430351,  11.03204968,   -2.99994917, -4.71595145, 0.99752609 ],
        [ -4.71262388,  11.02701128,   -2.99994891, -4.71427598, 0.99751984 ],
        [ -4.7102724,   11.0199565,    -2.99994855, -4.71193034, 0.99751107 ],
        [ -4.70870475,  11.01525435,   -2.99994831, -4.7103666,  0.9975052 ],
        [ -4.70646525,  11.00853496,   -2.99994796, -4.70813269, 0.99749679 ],
        [ -4.70433772,  11.00215353,   -2.99994762, -4.7060105,  0.99748878 ],
        [ -4.70232217,  10.99610577,   -2.99994731, -4.70400002, 0.99748116 ],
        [ -4.69985872,  10.98871706,   -2.99994692, -4.70154278, 0.99747182 ],
        [ -4.69750724,  10.98166182,   -2.99994654, -4.69919726, 0.99746287 ],
        [ -4.68407022,  10.94135182,   -2.99994434, -4.68579468, 0.99741113 ],
        [ -4.33470779,  9.89329863,    -2.9998476,  -4.33762181, 0.9956228 ],
        [ -4.06999856,  9.09923512,    -2.99966182, -4.07433599, 0.99348043 ],
        [ -3.85052729,  8.44093043,    -2.9993235,  -3.85656054, 0.99092465 ],
        [ -3.60373409,  7.70079879,    -2.99857754, -3.61248114, 0.98682754 ],
        [ -3.3672426,   6.99181339,    -2.99711491, -3.37973475, 0.98116073 ],
        [ -3.16411969,  6.38325067,    -2.99470809, -3.1810932,  0.97436368 ],
        [ -2.97992059,  5.83195159,    -2.99084864, -3.00234335, 0.96608059 ],
        [ -2.80635913,  5.3133438,     -2.98471954, -2.83551992, 0.9558217 ],
        [ -2.6803871,   4.937759,      -2.97786812, -2.71568299, 0.94647531 ],
        [ -2.54512113,  4.53562958,    -2.96721809, -2.58845668, 0.93423635 ],
        [ -2.41836527,  4.16037368,    -2.9528851,  -2.47089354, 0.92027888 ],
        [ -2.30135356,  3.81586894,    -2.93456462, -2.3640864,  0.90486432 ],
        [ -2.18534732,  3.47679124,    -2.91011656, -2.26013532, 0.88680808 ],
        [ -2.06844528,  3.13841085,    -2.87754703, -2.15768124, 0.86545424 ],
        [ -1.9458325,   2.7881964,     -2.83308769, -2.05312789, 0.83930428 ],
        [ -1.81459764,  2.42023651,    -2.77216607, -1.94506333, 0.80679285 ],
        [ -1.69623992,  2.09600373,    -2.7047208,  -1.85151002, 0.77343938 ],
        [ -1.58648321,  1.80304975,    -2.63196824, -1.76846631, 0.73931003 ],
        [ -1.4824583,   1.53319737,    -2.55510693, -1.69335102, 0.70451277 ],
        [ -1.3824874,   1.28169655,    -2.47560509, -1.62466705, 0.66932939 ],
        [ -1.28450916,  1.0431105,     -2.39415433, -1.56082334, 0.63376493 ],
        [ -1.18740068,  0.81460935,    -2.31183456, -1.50101162, 0.59806501 ],
        [ -1.09020626,  0.5939159,     -2.22960677, -1.44461566, 0.56248613 ],
        [ -0.99200406,  0.37897813,    -2.1482485,  -1.39111598, 0.52726039 ],
        [ -0.89188035,  0.16791668,    -2.06837169, -1.34007241, 0.492597 ],
        [ -0.78908717,  -0.04066009,   -1.99058935, -1.29119333, 0.45874209 ],
        [ -0.68284566,  -0.2480995,    -1.91538141, -1.24422044, 0.42590871 ],
        [ -0.5723262,   -0.45573591,   -1.8431094,  -1.19892153, 0.39427402 ],
        [ -0.45662985,  -0.66491684,   -1.77403702, -1.15508574, 0.36397894 ],
        [ -0.33491288,  -0.87678249,   -1.70842061, -1.11256818, 0.33516184 ],
        [ -0.20612353,  -1.0927351,    -1.64637301, -1.07119301, 0.30789491 ],
        [ -0.06917793,  -1.31411618,   -1.58798747, -1.03082236, 0.28223553 ],
        [ 0.0770387,    -1.54222016,   -1.53333439, -0.99135021, 0.25822405 ],
        [ 0.23391587,   -1.77867198,   -1.48237837, -0.95263833, 0.23584874 ],
        [ 0.39874979,   -2.01911448,   -1.43616341, -0.91547594, 0.21556435 ],
        [ 0.56847565,   -2.25930463,   -1.39519114, -0.88045299, 0.19758481 ],
        [ 0.7480148,    -2.50636759,   -1.35795275, -0.8464836,  0.18123825 ],
        [ 0.93691686,   -2.75963501,   -1.32438752, -0.81367623, 0.16648876 ],
        [ 1.13466809,   -3.01848454,   -1.2943284,  -0.78209497, 0.15325498 ],
        [ 1.34500994,   -3.28778829,   -1.26702576, -0.7511599,  0.14119857 ],
        [ 1.56982142,   -3.56976185,   -1.24218006, -0.72068787, 0.13017998 ],
        [ 1.80911542,   -3.86424718,   -1.21971981, -0.69076674, 0.12016194 ],
        [ 2.06561138,   -4.17441165,   -1.19932278, -0.66115224, 0.11099662 ],
        [ 2.33898244,   -4.49968905,   -1.18094365, -0.63197828, 0.10266206 ],
        [ 2.63250467,   -4.84380883,   -1.16429206, -0.60299491, 0.0950268 ],
        [ 2.94770694,   -5.20835518,   -1.14923669, -0.57417374, 0.08803253 ],
        [ 3.28665223,   -5.5955108,    -1.13563249, -0.54545171, 0.08161608 ],
        [ 3.65159506,   -6.00764378,   -1.12334559, -0.51677065, 0.07572065 ],
        [ 4.04538216,   -6.44775174,   -1.11224235, -0.48805017, 0.07029041 ],
        [ 4.47065186,   -6.91855868,   -1.10221434, -0.45925062, 0.0652825 ],
        [ 4.9308357,    -7.42363151,   -1.09314963, -0.43030182, 0.06065286 ],
        [ 5.42955854,   -7.96670395,   -1.08495052, -0.40115008, 0.05636478 ],
        [ 5.97123953,   -8.55232843,   -1.07752293, -0.37172488, 0.0523832 ],
        [ 6.56029281,   -9.18500336,   -1.07078888, -0.34198719, 0.04868112 ],
        [ 7.20212856,   -9.87025158,   -1.0646726,  -0.31187767, 0.04523204 ],
        [ 7.90255373,   -10.61396753,  -1.059109,   -0.28135238, 0.04201426 ],
        [ 8.66797285,   -11.4226313,   -1.05404059, -0.2503741,  0.03900914 ],
        [ 9.50518872,   -12.30309725,  -1.04941809, -0.21892129, 0.03620139 ],
        [ 10.42160298,  -13.26280687,  -1.04519803, -0.18697976, 0.03357772 ],
        [ 11.42581651,  -14.31041336,  -1.0413402,  -0.15452472, 0.03112513 ],
        [ 12.52582974,  -15.45390346,  -1.03781397, -0.12157978, 0.02883512 ],
        [ 13.73124287,  -16.7028957,   -1.03458879, -0.08814429, 0.02669776 ],
        [ 15.05265599,  -18.06800622,  -1.03163737, -0.05421847, 0.02470378 ],
        [ 16.50006912,  -19.55919957,  -1.02893815, -0.01984236, 0.02284655 ],
        [ 18.08628227,  -21.18929689,  -1.026468,   0.01498821,  0.02111714 ],
        [ 19.68231065,  -22.82585465,  -1.02437764, 0.04747578,  0.01963018 ],
        [ 21.49965911,  -24.68562774,  -1.02236886, 0.08179808,  0.01817999 ],
        [ 23.56170659,  -26.79177244,  -1.02045866, 0.11780347,  0.01678052 ],
        [ 25.68615592,  -28.9578802,   -1.01880611, 0.15211487,  0.01555293 ],
        [ 28.09096023,  -31.4059616,   -1.01723221, 0.18805371,  0.01436842 ],
        [ 30.53642375,  -33.89186955,  -1.01588185, 0.22190128,  0.01333965 ],
        [ 33.29444525,  -36.69186614,  -1.01459324, 0.2572867,   0.0123466 ],
        [ 36.41950658,  -39.86055737,  -1.01336546, 0.29433412,  0.01138963 ],
        [ 39.55913824,  -43.04047321,  -1.01232449, 0.32877324,  0.01056958 ],
        [ 43.09857555,  -46.62172331,  -1.01133028, 0.36474667,  0.00977857 ],
        [ 47.1074426,   -50.67405614,  -1.01038223, 0.40237875,  0.00901686 ],
        [ 51.06692048,  -54.67303242,  -1.00959009, 0.43677729,  0.00837456 ],
        [ 55.52051278,  -59.16760156,  -1.00883245, 0.47265859,  0.00775502 ],
        [ 60.55267651,  -64.24233915,  -1.0081089,  0.51014117,  0.00715838 ],
        [ 66.26660711,  -70.00057414,  -1.00741908, 0.54935856,  0.00658479 ],
        [ 71.80248083,  -75.57591447,  -1.00685435, 0.58446826,  0.00611159 ],
        [ 78.02992812,  -81.84431777,  -1.00631386, 0.62107042,  0.0056555 ],
        [ 85.06760126,  -88.92455666,  -1.00579739, 0.6592854,   0.00521661 ],
        [ 93.06060496,  -96.96187152,  -1.0053047,  0.6992488,   0.004795 ],
        [ 100.57907587, -104.51871195, -1.00491211, 0.73399749,  0.0044569 ],
        [ 109.01068314, -112.9901073,  -1.00453574, 0.77016678,  0.00413091 ],
        [ 118.50788348, -122.52862629, -1.00417544, 0.80786871,  0.00381705 ],
        [ 129.25619515, -133.31991313, -1.00383108, 0.84722917,  0.00351536 ],
        [ 141.48317094, -145.59166287, -1.00350254, 0.88839039,  0.00322588 ],
        [ 152.51638153, -156.6620955,  -1.00325099, 0.92272369,  0.00300312 ],
        [ 164.85028883, -169.03457202, -1.00300941, 0.95840539,  0.00278821 ],
        [ 178.69587545, -182.92017885, -1.00277772, 0.99553933,  0.00258118 ],
        [ 194.30855527, -198.57444627, -1.00255585, 1.03424186,  0.00238204 ],
        [ 211.99987327, -216.30905027, -1.00234375, 1.07464401,  0.0021908 ]
    ],

    energy_table => [
        [ -5.03701586,  0.32842029, 0.08954992, 0.32842029, 0.08954992,    5.7680088,   -1.5981387 ],
        [ -4.72942015,  0.35715012, 0.09735788, 0.35715012, 0.09735788,    5.2730329,   -1.6208789 ],
        [ -4.72572497,  0.35751005, 0.09745555, 0.35751005, 0.09745555,    5.2670429,   -1.6211682 ],
        [ -4.71833461,  0.35823101, 0.09765115, 0.35823101, 0.09765115,    5.2550597,   -1.6217637 ],
        [ -4.71654301,  0.358406,   0.09769865, 0.358406,   0.09769865,    5.2521541,   -1.6219084 ],
        [ -4.71430351,  0.35862486, 0.09775801, 0.35862486, 0.09775801,    5.2485216,   -1.6220894 ],
        [ -4.71262388,  0.3587891,  0.0978026,  0.3587891,  0.0978026,     5.245797,    -1.6222252 ],
        [ -4.7102724,   0.35901915, 0.09786498, 0.35901915, 0.09786498,    5.2419821,   -1.6224156 ],
        [ -4.70870475,  0.3591726,  0.09790667, 0.3591726,  0.09790667,    5.2394386,   -1.6225426 ],
        [ -4.70646525,  0.35939193, 0.09796608, 0.35939193, 0.09796608,    5.2358047,   -1.6227242 ],
        [ -4.70433772,  0.35960042, 0.09802273, 0.35960042, 0.09802273,    5.2323522,   -1.6228968 ],
        [ -4.70232217,  0.35979804, 0.09807623, 0.35979804, 0.09807623,    5.229081,    -1.6230605 ],
        [ -4.69985872,  0.36003973, 0.09814193, 0.36003973, 0.09814193,    5.2250824,   -1.6232608 ],
        [ -4.69750724,  0.36027058, 0.09820446, 0.36027058, 0.09820446,    5.2212651,   -1.6234523 ],
        [ -4.68407022,  0.36159256, 0.09856304, 0.36159256, 0.09856304,    5.1994434,   -1.6245414 ],
        [ -4.33470779,  0.39770714, 0.10832352, 0.39770714, 0.10832352,    4.626972,    -1.6560659 ],
        [ -4.06999856,  0.42742344, 0.11627409, 0.42742344, 0.11627409,    4.1850965,   -1.6857272 ],
        [ -3.85052729,  0.45369685, 0.12319584, 0.45369685, 0.12319584,    3.8121361,   -1.7162706 ],
        [ -3.60373409,  0.48508982, 0.13124666, 0.48508982, 0.13124666,    3.3838788,   -1.758589 ],
        [ -3.3672426,   0.51705191, 0.13904733, 0.51705191, 0.13904733,    2.962707,    -1.8101008 ],
        [ -3.16411969,  0.54596308, 0.14556461, 0.54596308, 0.14556461,    2.5899426,   -1.8679364 ],
        [ -2.97992059,  0.57328826, 0.15102596, 0.57328826, 0.15102596,    2.2403966,   -1.9331307 ],
        [ -2.80635913,  0.59989556, 0.15542132, 0.59989556, 0.15542132,    1.8990777,   -2.0053823 ],
        [ -2.6803871,   0.61963786, 0.15789709, 0.61963786, 0.15789709,    1.642906,    -2.0660977 ],
        [ -2.54512113,  0.64112453, 0.15961138, 0.64112453, 0.15961138,    1.3587073,   -2.141288 ],
        [ -2.41836527,  0.66139799, 0.16006007, 0.66139799, 0.16006007,    1.0825054,   -2.2232958 ],
        [ -2.30135356,  0.68009135, 0.15923481, 0.68009135, 0.15923481,    0.81756985,  -2.33551 ],
        [ -2.18534732,  0.69844972, 0.1570242,  0.69844972, 0.1570242,     0.53843242,  -2.4500006 ],
        [ -2.06844528,  0.71659949, 0.1532104,  0.71659949, 0.1532104,     0.2468657,   -2.5400469 ],
        [ -1.9458325,   0.73504654, 0.14736912, 0.73504654, 0.14736912,    -0.07048393, -2.6706123 ],
        [ -1.81459764,  0.7538614,  0.13901186, 0.7538614,  0.13901186,    -0.43253356, -2.8003503 ],
        [ -1.69623992,  0.7697815,  0.12975376, 0.7697815,  0.12975376,    -0.76841278, -2.796197 ],
        [ -1.58648321,  0.78349503, 0.11997034, 0.78349503, 0.11997034,    -1.0710766,  -2.6497135 ],
        [ -1.4824583,   0.7954584,  0.10994575, 0.7954584,  0.10994575,    -1.3360772,  -2.3760266 ],
        [ -1.3824874,   0.80595149, 0.09994429, 0.80595149, 0.09994429,    -1.5571391,  -1.9406023 ],
        [ -1.28450916,  0.81526074, 0.09010003, 0.81526074, 0.09010003,    -1.7212856,  -1.5951024 ],
        [ -1.18740068,  0.82354492, 0.08057626, 0.82354492, 0.08057626,    -1.8684613,  -1.5276568 ],
        [ -1.09020626,  0.8309311,  0.07150441, 0.8309311,  0.07150441,    -2.0181159,  -1.5658193 ],
        [ -0.99200406,  0.83752859, 0.06297945, 0.83752859, 0.06297945,    -2.17447,    -1.6271056 ],
        [ -0.89188035,  0.84343127, 0.05506546, 0.84343127, 0.05506546,    -2.3409487,  -1.6939411 ],
        [ -0.78908717,  0.84871121, 0.0478127,  0.84871121, 0.0478127,     -2.5183681,  -1.7531197 ],
        [ -0.68284566,  0.85343378, 0.04124411, 0.85343378, 0.04124411,    -2.7076018,  -1.8041924 ],
        [ -0.5723262,   0.85765833, 0.03535951, 0.85765833, 0.03535951,    -2.9096475,  -1.8473053 ],
        [ -0.45662985,  0.86143863, 0.03014013, 0.86143863, 0.03014013,    -3.1256945,  -1.8830314 ],
        [ -0.33491288,  0.86481959, 0.02555813, 0.86481959, 0.02555813,    -3.3568979,  -1.9121196 ],
        [ -0.20612353,  0.86784572, 0.02156952, 0.86784572, 0.02156952,    -3.6048759,  -1.9354182 ],
        [ -0.06917793,  0.87055536, 0.01812621, 0.87055536, 0.01812621,    -3.8713738,  -1.9537428 ],
        [ 0.0770387,    0.87298202, 0.01517762, 0.87298202, 0.01517762,    -4.1582503,  -1.9678613 ],
        [ 0.23391587,   0.87515846, 0.01266865, 0.87515846, 0.01266865,    -4.46795,    -1.9785076 ],
        [ 0.39874979,   0.87706878, 0.01059334, 0.87706878, 0.01059334,    -4.7948291,  -1.9861657 ],
        [ 0.56847565,   0.87871841, 0.00891217, 0.87871841, 0.00891217,    -5.1324718,  -1.9914319 ],
        [ 0.7480148,    0.88018782, 0.00751267, 0.88018782, 0.00751267,    -5.4904089,  -1.9950672 ],
        [ 0.93691686,   0.88149327, 0.00635511, 0.88149327, 0.00635511,    -5.8675641,  -1.9974825 ],
        [ 1.13466809,   0.88265194, 0.00540121, 0.88265194, 0.00540121,    -6.2627587,  -1.998993 ],
        [ 1.34500994,   0.88370062, 0.0046017,  0.88370062, 0.0046017,     -6.6833536,  -1.9998969 ],
        [ 1.56982142,   0.88465648, 0.00392875, 0.88465648, 0.00392875,    -7.13303,    -2.0004014 ],
        [ 1.80911542,   0.88552642, 0.00336456, 0.88552642, 0.00336456,    -7.6117558,  -2.0006315 ],
        [ 2.06561138,   0.88632587, 0.00288797, 0.88632587, 0.00288797,    -8.124925,   -2.0006727 ],
        [ 2.33898244,   0.88705843, 0.00248723, 0.88705843, 0.00248723,    -8.6718456,  -2.0006169 ],
        [ 2.63250467,   0.88773664, 0.00214726, 0.88773664, 0.00214726,    -9.2590598,  -2.0005374 ],
        [ 2.94770694,   0.88836618, 0.00185846, 0.88836618, 0.00185846,    -9.8896199,  -2.0004372 ],
        [ 3.28665223,   0.88895278, 0.00161243, 0.88895278, 0.00161243,    -10.567638,  -2.000329 ],
        [ 3.65159506,   0.8895014,  0.00140221, 0.8895014,  0.00140221,    -11.297625,  -2.0002387 ],
        [ 4.04538216,   0.89001671, 0.00122189, 0.89001671, 0.00122189,    -12.085277,  -2.0001655 ],
        [ 4.47065186,   0.89050212, 0.00106682, 0.89050212, 0.00106682,    -12.935872,  -2.0001096 ],
        [ 4.9308357,    0.8909611,  0.00093295, 0.8909611,  0.00093295,    -13.856279,  -2.0000697 ],
        [ 5.42955854,   0.8913964,  0.00081698, 0.8913964,  0.00081698,    -14.853751,  -2.0000424 ],
        [ 5.97123953,   0.89181067, 0.00071638, 0.89181067, 0.00071638,    -15.93713,   -2.0000247 ],
        [ 6.56029281,   0.89220579, 0.0006346,  0.89220579, 0.0006346,     -17.115247,  -2.0000135 ],
        [ 7.20212856,   0.89258377, 0.00055184, 0.89258377, 0.00055184,    -18.398924,  -2.0000072 ],
        [ 7.90255373,   0.89294589, 0.00048463, 0.89296022, 0.00051922608, -19.799778,  -2.0000036 ],
        [ 8.66797285,   0.89329345, 0.00042567, 0.89335014, 0.0005068091,  -21.330618,  -2.0000017 ],
        [ 9.50518872,   0.89362735, 0.00037387, 0.89376894, 0.00049473963, -23.005051,  -2.0000008 ],
        [ 10.42160298,  0.89394834, 0.00032833, 0.89421714, 0.00048069609, -24.83788,   -2.0000004 ],
        [ 11.42581651,  0.8942572,  0.00028825, 0.89468838, 0.00045793765, -26.846307,  -2.0000002 ],
        [ 12.52582974,  0.89455419, 0.00025301, 0.89517813, 0.00043218506, -29.046334,  -2.0000001 ],
        [ 13.73124287,  0.89483979, 0.000222,   0.89572156, 0.0004624296,  -31.45716,   -2.0000001 ],
        [ 15.05265599,  0.89511446, 0.00019472, 0.89629319, 0.0004048863,  -34.099987,  -2.0000001 ],
        [ 16.50006912,  0.89537833, 0.00017076, 0.89684139, 0.00035447856, -36.994813,  -2.0000001 ],
        [ 18.08628227,  0.89563187, 0.00014969, 0.89736732, 0.00031029993, -40.167239,  -2.0000001 ],
        [ 19.68231065,  0.89585664, 0.00013256, 0.89783297, 0.00027444197, -43.359296,  -2 ],
        [ 21.49965911,  0.89608264, 0.00011672, 0.8983006,  0.00024136702, -46.993993,  -2 ],
        [ 23.56170659,  0.89630785, 0.00010226, 0.89876608, 0.00021123218, -51.118088,  -2 ],
        [ 25.68615592,  0.89651187, 9.024e-05,  0.89918736, 0.00018624724, -55.366987,  -2 ],
        [ 28.09096023,  0.89671519, 7.926e-05,  0.89960681, 0.00016343191, -60.176595,  -2 ],
        [ 30.53642375,  0.89689756, 7.02e-05,   0.89998274, 0.00014465458, -65.067522,  -2 ],
        [ 33.29444525,  0.89707931, 6.19e-05,   0.90035714, 0.00012745394, -70.583565,  -2 ],
        [ 36.41950658,  0.89726042, 5.43e-05,   0.90073,    0.00011174735, -76.833688,  -2 ],
        [ 39.55913824,  0.89742087, 4.812e-05,  0.90106012, 9.8975355e-05, -83.112951,  -2 ],
        [ 43.09857555,  0.89758079, 4.245e-05,  0.90138899, 8.7264919e-05, -90.191826,  -2 ],
        [ 47.1074426,   0.89774017, 3.726e-05,  0.90171659, 7.6561908e-05, -98.20956,   -2 ],
        [ 51.06692048,  0.89787917, 3.31e-05,   0.90200218, 6.7981516e-05, -106.12852,  -2 ],
        [ 55.52051278,  0.89801773, 2.927e-05,  0.90228679, 6.0097075e-05, -115.0357,   -2 ],
        [ 60.55267651,  0.89815586, 2.576e-05,  0.90257038, 5.2874284e-05, -125.10003,  -2 ],
        [ 66.26660711,  0.89829353, 2.255e-05,  0.90285297, 4.6279798e-05, -136.52789,  -2 ],
        [ 71.80248083,  0.89841117, 2.003e-05,  0.90309436, 4.1102435e-05, -147.59964,  -2 ],
        [ 78.02992812,  0.89852846, 1.772e-05,  0.903335,   3.6342036e-05, -160.05453,  -2 ],
        [ 85.06760126,  0.8986454,  1.559e-05,  0.90357486, 3.1978696e-05, -174.12988,  -2 ],
        [ 93.06060496,  0.89876199, 1.365e-05,  0.90381394, 2.7992711e-05, -190.11588,  -2 ],
        [ 100.57907587, 0.89885887, 1.217e-05,  0.90401258, 2.494529e-05,  -205.15283,  -2 ],
        [ 109.01068314, 0.8989555,  1.08e-05,   0.90421067, 2.2135477e-05, -222.01604,  -2 ],
        [ 118.50788348, 0.89905187, 9.54e-06,   0.90440821, 1.955233e-05,  -241.01044,  -2 ],
        [ 129.25619515, 0.89914798, 8.39e-06,   0.90460519, 1.7185161e-05, -262.50707,  -2 ],
        [ 141.48317094, 0.89924383, 7.33e-06,   0.9048016,  1.5023475e-05, -286.96102,  -2 ],
        [ 152.51638153, 0.89932031, 6.56e-06,   0.90495832, 1.3434998e-05, -309.02744,  -2 ],
        [ 164.85028883, 0.89939663, 5.84e-06,   0.90511468, 1.1966109e-05, -333.69525,  -2 ],
        [ 178.69587545, 0.89947277, 5.18e-06,   0.90527067, 1.0611566e-05, -361.38643,  -2 ],
        [ 194.30855527, 0.89954873, 4.57e-06,   0.90542628, 9.3662365e-06, -392.61179,  -2 ],
        [ 211.99987327, 0.89962452, 4.02e-06,   0.90558153, 8.2251187e-06, -427.99442,  -2 ]
    ],

    impulse_table => [
        [
            -5.6419056,  13.81484,   -5.6423152,    0.010487848, -0.00025113065, -0.12441001,
            -0.92721557, -0.5760752, -0.0006805186, -0.26907232, 0.082446019,    0.69367601,
            -1.6243953
        ],
        [
            -4.9618453,  11.774667,   -4.9629817,    0.013897969, -0.00049573106, -0.12095612,
            -0.92376167, -0.91611404, -0.0013433415, -0.26561843, 0.082445388,    0.63125841,
            -1.5780887
        ],
        [
            -4.6051704,  10.704657,  -4.6071115,    0.015933504, -0.00070818723, -0.11795612,
            -0.92076167, -1.0944644, -0.0019190592, -0.26261843, 0.082444344,    0.59360266,
            -1.5587312
        ],
        [
            -4.1997053,  9.4883176,  -4.2032743,    0.018372765, -0.0010622808, -0.11295612,
            -0.91576167, -1.2972416, -0.0028785888, -0.25761843, 0.082440812,   0.54614037,
            -1.5432475
        ],
        [
            -3.9120232,  8.6253798,  -3.9175233,    0.020104469, -0.0014163745, -0.10795612,
            -0.91076168, -1.4411692, -0.0038381184, -0.25261843, 0.082433991,   0.50919907,
            -1.5382302
        ],
        [
            -3.5065581,  7.4094307,  -3.5166835,    0.022345525, -0.0021245617, -0.097956125,
            -0.90076168, -1.6442572, -0.0057571776, -0.24261843, 0.08240602,    0.45219473,
            -1.5434601
        ],
        [
            -2.9957324,  5.8792457,  -3.0176252,    0.024219932, -0.0035409362, -0.077956128,
            -0.88076168, -1.9015035, -0.0095952961, -0.22261843, 0.082262299,   0.37197953,
            -1.5819559
        ],
        [
            -2.3025853,  3.8194836,  -2.365201,    0.022715892, -0.0070818723, -0.027956137,
            -0.83076169, -2.2644695, -0.019190592, -0.17261844, 0.081017648,   0.2518985,
            -1.7402723
        ],
        [
            -1.8971202, 2.6506918,  -2.0125212,   0.018421452, -0.010622808, 0.022043854,
            -0.7807617, -2.5120588, -0.028785888, -0.12261845, 0.077944105,  0.18295753,
            -1.9343205
        ],
        [
            -1.6094381,  1.8636505,  -1.785522,    0.01457931,   -0.014163697, 0.072040861,
            -0.73076347, -2.7461928, -0.038381138, -0.072618459, 0.072924981,  0.13988701,
            -2.1379013
        ],
        [
            -1.3862945,  1.2911274,  -1.6272179,   0.013263503, -0.017640563, 0.11925851,
            -0.68257931, -3.0233562, -0.047906157, -0.02313238, 0.066491293,  0.11182703,
            -2.3354994
        ],
        [
            -1.203973,   0.85303859, -1.5109735,   0.014251196,  -0.01935328, 0.13819279,
            -0.65967555, -3.405867,  -0.053202003, 0.0035300445, 0.059449744, 0.092851523,
            -2.518856
        ],
        [
            -1.0498223,  0.50455838, -1.4221957,  0.015284827,  -0.019586818, 0.14556633,
            -0.65013514, -3.8059015, -0.05351098, 0.0048646945, 0.052528133,  0.079538009,
            -2.6853237
        ],
        [
            -0.91629091, 0.21863993, -1.3521983,   0.016045761,  -0.019652217, 0.150595,
            -0.64386869, -4.1645799, -0.053295999, 0.0030672185, 0.046220964,  0.069864083,
            -2.8353486
        ],
        [
            -0.79850787, -0.021874713, -1.2955292,   0.016591335,   -0.019690545, 0.15450386,
            -0.63923811, -4.4830027,   -0.053057666, 0.00053729514, 0.040783588,  0.062606908,
            -2.9705906
        ],
        [
            -0.69314735, -0.22833155, -1.2486239,   0.01698684,   -0.019720301, 0.15768047,
            -0.63564785, -4.7674522,  -0.052845437, -0.002045071, 0.03628689,   0.05700673,
            -3.0929687
        ],
        [
            -0.59783717, -0.40851128, -1.2090697,   0.017278865,   -0.019745282, 0.16032404,
            -0.63278565, -5.0236308,  -0.052660994, -0.0046798554, 0.03268176,   0.052577394,
            -3.2042812
        ],
        [
            -0.5108258,  -0.5679236, -1.1751839,   0.017498633,   -0.019766749, 0.16256067,
            -0.63045863, -5.2561625, -0.052500142, -0.0072649135, 0.029850928,  0.048998484,
            -3.3060901
        ],
        [
            -0.46603967, -0.64819824, -1.1585218,   0.01759652,    -0.019777308, 0.16364947,
            -0.62935893, -5.3752685,  -0.052420251, -0.0085048923, 0.028557803,  0.047311878,
            -3.3584879
        ],
        [
            -0.35667512, -0.83948181, -1.1199155,   0.017798163,  -0.01980167, 0.16613915,
            -0.62693168, -5.6642795,  -0.052233402, -0.012314355, 0.025783051, 0.043589009,
            -3.4862455
        ],
        [
            -0.28768225, -0.95691107, -1.0969849,   0.017901812, -0.019815997, 0.16759232,
            -0.62557514, -5.8451789,  -0.052121521, -0.01461658, 0.024269759,  0.041499403,
            -3.5666033
        ],
        [
            -0.22314373, -1.0646475, -1.0764625,   0.017984817, -0.019828678, 0.16887391,
            -0.62441872, -6.0133593, -0.052021018, -0.01688979, 0.022994359,  0.039705361,
            -3.6415481
        ],
        [
            -0.1625191,  -1.1640987, -1.0579545,   0.018052078,  -0.019839967, 0.17001238,
            -0.62342454, -6.1704075, -0.051930244, -0.019081392, 0.021904847,  0.03814846,
            -3.7117184
        ],
        [
            -0.10536069, -1.2563939, -1.0411506,   0.018107162,  -0.019850074, 0.17103028,
            -0.62256398, -6.3176411, -0.051847816, -0.021225158, 0.020963279,  0.036784488,
            -3.7776532
        ],
        [
            -0.051293468, -1.342453,  -1.0258026,   0.018152701,  -0.019859168, 0.17194565,
            -0.62181417,  -6.4561642, -0.051772678, -0.023230435, 0.020141322,  0.035579438,
            -3.8398089
        ],
        [
            -1.7392799e-07, -1.4230336, -1.0117092,   0.018190678,  -0.019867397, 0.17277319,
            -0.62115672,    -6.5869098, -0.051703894, -0.025263102, 0.019417403,  0.034506771,
            -3.8985759
        ],
        [
            0.095310006, -1.5701783, -0.98665752,  0.018249595,  -0.019881689, 0.17421055,
            -0.62006452, -6.8281346, -0.051582432, -0.029197112, 0.018200426,  0.0326788,
            -4.0072377
        ],
        [
            0.18232138,  -1.7017815, -0.9649859,   0.018292306,  -0.019893687, 0.17541597,
            -0.61920055, -7.046435,  -0.051478625, -0.032744789, 0.017216436,  0.031177247,
            -4.1058156
        ],
        [
            0.26236409,  -1.8207221, -0.94598195,  0.018323969,  -0.019903897, 0.17644111,
            -0.61850584, -7.2456715, -0.051388904, -0.035963985, 0.016403423,  0.029919992,
            -4.1959695
        ],
        [
            0.40546493,  -2.0287527, -0.91403096,  0.018366249, -0.019920363, 0.17809102,
            -0.61747272, -7.5982263, -0.051241743, -0.04223426, 0.015135683,  0.027927893,
            -4.3558947
        ],
        [
            0.53062808,  -2.2063382, -0.88800213,  0.01839176,   -0.019933064, 0.17936032,
            -0.61675542, -7.9029359, -0.051126177, -0.047789358, 0.014189331,  0.026414084,
            -4.4944915
        ],
        [
            0.69314701,  -2.4315666, -0.85655645,  0.018413607,  -0.019947505, 0.18079596,
            -0.61603314, -8.2938282, -0.050993149, -0.054780298, 0.013142158,  0.024710243,
            -4.6727585
        ],
        [
            0.9932516,   -2.8339875, -0.80440995,  0.018432972,  -0.019969018, 0.18291859,
            -0.61516098, -9.0028575, -0.050791768, -0.068585638, 0.011616824,  0.022168697,
            -4.9973675
        ],
        [
            1.0986121,  -2.9717243, -0.78766128,  0.01843584,   -0.01997526, 0.18352802,
            -0.6149584, -9.2482403, -0.050733055, -0.073345336, 0.011178021, 0.021423109,
            -5.1100683
        ],
        [
            1.2089602,   -3.1142629, -0.77087693,  0.018437466,  -0.019981144, 0.18410195,
            -0.61478762, -9.5034529, -0.050677332, -0.078103077, 0.010761295,  0.020708528,
            -5.2274751
        ],
        [
            1.3862942, -3.3399956, -0.74537528, 0.018438119,  -0.019989449, 0.18490352,
            -0.614586, -9.9100478, -0.05059892, -0.086875091, 0.010169761,  0.019682667,
            -5.4149174
        ],
        [
            1.6094377,   -3.6188934, -0.71556578,  0.018436804,  -0.019998068, 0.18573158,
            -0.61442286, -10.416072, -0.050517153, -0.096508367, 0.0095375196, 0.018570191,
            -5.6488635
        ],
        [
            1.7917593,   -3.8430645, -0.69285812,  0.018434827, -0.020003835, 0.18628453,
            -0.61434155, -10.825404, -0.050462087, -0.10496274, 0.009095363,  0.017781608,
            -5.8386413
        ],
        [
            1.94591,     -4.0303141, -0.67468068,  0.018432902, -0.020007972, 0.18668005,
            -0.61429897, -11.16891,  -0.050422924, -0.11254954, 0.0087641757, 0.017184871,
            -5.9982596
        ],
        [
            2.0794414,   -4.1909913, -0.6596203,  0.018431172, -0.020011071, 0.18697675,
            -0.61427438, -11.464733, -0.05039356, -0.11919473, 0.0085041511, 0.01671255,
            -6.135988
        ],
        [
            2.302786,   -4.4569025, -0.63571281,  0.018428476, -0.02001562,  0.18739386,
            -0.6142462, -11.956277, -0.050351678, -0.12986753, 0.0081164036, 0.015878217,
            -6.3720397
        ],
        [
            2.9958308,  -5.2636108, -0.56996077, 0.018421677, -0.020025672, 0.18823004,
            -0.6143033, -13.460414, -0.05026925, -0.16295358, 0.0071891956, 0.014219202,
            -7.0741424
        ],
        [
            3.9120691,   -6.2992421, -0.49753584,  0.018416852, -0.020031596, 0.18873296,
            -0.61436259, -15.413956, -0.050219383, -0.20680212, 0.0063588866, 0.012659881,
            -7.9967744
        ],
        [
            4.6055151,   -7.0670147, -0.45054369,  0.018415149, -0.020033393, 0.18890025,
            -0.61437578, -16.874876, -0.050202506, -0.24001715, 0.0059023115, 0.0117771,
            -8.6926528
        ],
        [
            5.2983318,   -7.8241985, -0.40861589, 0.018414292, -0.020034286, 0.18898385,
            -0.61438278, -18.324021, -0.05019407, -0.27320792, 0.0055378911, 0.011062531,
            -9.3867954
        ],
        [
            6.2146648,  -8.8142637, -0.35917025,  0.018413781, -0.020034819, 0.18903406,
            -0.6143871, -20.22903,  -0.050189004, -0.31710796, 0.0051507945, 0.010296488,
            -10.303994
        ],
        [
            6.9079142,   -9.5566233, -0.32540552,  0.018413614, -0.020034996, 0.1890508,
            -0.61438856, -21.663603, -0.050187316, -0.35031981, 0.0049088475, 0.0098152017,
            -10.997554
        ],
        [
            7.6009441,  -10.294192, -0.29422156,  0.018413531, -0.020035081, 0.18905917,
            -0.6143893, -23.093254, -0.050186472, -0.38352047, 0.0046993813, 0.0093975388,
            -11.690748
        ],
        [
            8.5173695,   -11.26382,  -0.25629049,  0.018413482, -0.020035136, 0.18906419,
            -0.61438974, -24.978214, -0.050185966, -0.42742246, 0.0044607246, 0.008920969,
            -12.607277
        ],
        [
            9.2103836,   -11.9935,   -0.22973058,  0.018413465, -0.020035155, 0.18906586,
            -0.61438988, -26.400199, -0.050185797, -0.46062138, 0.004303167,  0.0086060962,
            -13.300328
        ],
        [
            10.819809, -13.678686, -0.17381391, 0.01841345,  -0.020035165, 0.1890672,
            -0.61439,  -29.693489, -0.05012834, -0.52876854, 0.003994462,  0.0079888764,
            -14.909784
        ],
        [
            11.513018,   -14.401206, -0.15181909,  0.018413447, -0.020035164, 0.18906737,
            -0.61439002, -31.108763, -0.049833563, -0.5464737,  0.0038808454, 0.0077616665,
            -15.602997
        ],
        [
            13.122436,   -16.072565, -0.10470913,  0.018413443, -0.020035169, 0.1890675,
            -0.61439003, -34.388673, -0.048622221, -0.58570788, 0.0036509241, 0.0073018423,
            -17.212419
        ],
        [
            13.815644,   -16.790207, -0.085896766, 0.018413441, -0.020035167, 0.18906751,
            -0.61439004, -35.799214, -0.047989575, -0.60189737, 0.0035638921, 0.0071277806,
            -17.905627
        ],
        [
            16.118266,   -19.166222, -0.028652515, 0.018413437, -0.02001834,  0.18906754,
            -0.65310308, -40.477029, -0.045766497, -0.65310308, 0.0033144196, 0.0066288378,
            -20.20825
        ],
        [
            18.420685,   -21.532472, 0.021994418,  0.018413432, -0.020035151, 0.18906842,
            -0.61438915, -45.145097, -0.050185633, -0.90183568, 0.0031113392, 0.0062226772,
            -22.510669
        ]
    ],

    tail_shock_table => [
        [ 7.3723261, -0.431,      -0.029687701, -0.029687701 ],
        [ 7.6005717, -0.43809805, -0.037332324, -0.021464569 ],
        [ 8.517215,  -0.46558708, -0.043459453, -0.012071524 ],
        [ 9.2103043, -0.48540305, -0.044952799, -0.0084697955 ],
        [ 10.819792, -0.52876792, -0.045571255, -0.0038177456 ],
        [ 11.513009, -0.54647344, -0.045303265, -0.0026304292 ],
        [ 13.122435, -0.58570792, -0.044202032, -0.00083153935 ],
        [ 13.815643, -0.60190783, -0.043629221, -0.0003233208 ]
    ],

    blast_info => [
        0.12795606, 0.93076261, 0.2726629,  -0.19190606, 1.1392837,  0.0039531505, 2.9158309,   -0.070818764,
        0.31346722, 0.50196174, 0.44556117, -0.12974821, 0.45897934, 0.082445974,  0.016739491, -0.018213781,
        1590.9,     -0.431,     32622.619,  -0.51726287
    ],

};
1;
