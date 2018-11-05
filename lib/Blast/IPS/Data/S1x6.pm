package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=1.6
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S1x6'} = {

    table_name  => 'S1x6',
    symmetry    => 2,
    gamma       => 1.6,
    data_source => 'S8000_G1x6/moc2_from_moc_r100/',

    shock_table_info => [ 2, 1.6, 1.09e-06, 8000, 2.1e-07, 10, 1.28e-07, 4.6e-07, 5e-07, 1834.8, -0.68914 ],

    shock_table => [
        [ -4.50135843,  12.00017658,   -2.99997732, -4.50245891, 0.99834838 ],
        [ -4.13094467,  10.88894903,   -2.99993109, -4.1328636,  0.99711889 ],
        [ -3.84809547,  10.0404284,    -2.99985509, -3.85102997, 0.99559199 ],
        [ -3.43277968,  8.79459884,    -2.99950364, -3.4382576,  0.99176189 ],
        [ -3.12115287,  7.85997402,    -2.99873872, -3.12990834, 0.9868145 ],
        [ -2.8760393,   7.12508992,    -2.99737048, -2.88870732, 0.98089332 ],
        [ -2.66373443,  6.48895493,    -2.99504138, -2.68118806, 0.97363204 ],
        [ -2.47998323,  5.93891456,    -2.99142757, -2.50302645, 0.96513074 ],
        [ -2.31717138,  5.45227347,    -2.98610766, -2.34665769, 0.95531391 ],
        [ -2.17011852,  5.01367237,    -2.97857703, -2.20697105, 0.94408148 ],
        [ -2.0350765,   4.61209339,    -2.96823142, -2.08031379, 0.9313038 ],
        [ -1.90932842,  4.23966554,    -2.95436379, -1.96408551, 0.91683376 ],
        [ -1.79018687,  3.88871121,    -2.93605299, -1.85580235, 0.9004149 ],
        [ -1.67496445,  3.55172852,    -2.91203228, -1.75310696, 0.88163763 ],
        [ -1.56095269,  3.22143789,    -2.88051801, -1.6537979,  0.85989517 ],
        [ -1.44304287,  2.88419262,    -2.83806176, -1.5539097,  0.83378434 ],
        [ -1.31183091,  2.51563614,    -2.77713199, -1.44666016, 0.80015737 ],
        [ -1.194119,    2.19260637,    -2.70924471, -1.35445082, 0.76591238 ],
        [ -1.08550312,  1.9022424,     -2.63570731, -1.27312233, 0.7311615 ],
        [ -0.98293389,  1.63583322,    -2.55776968, -1.19991337, 0.69600873 ],
        [ -0.88560956,  1.39075087,    -2.47781848, -1.13386413, 0.66108767 ],
        [ -0.79020674,  1.15825668,    -2.39567317, -1.07246555, 0.6259588 ],
        [ -0.69491589,  0.93395713,    -2.31188382, -1.01450276, 0.59058924 ],
        [ -0.59953598,  0.71745014,    -2.22817673, -0.95985075, 0.55549446 ],
        [ -0.50309935,  0.5065849,     -2.14539254, -0.90795931, 0.52086686 ],
        [ -0.4046732,   0.29944959,    -2.06420886, -0.85837776, 0.48688264 ],
        [ -0.30348676,  0.09461492,    -1.98528186, -0.8108047,  0.45375287 ],
        [ -0.19862101,  -0.10952373,   -1.90902858, -0.76492667, 0.42162273 ],
        [ -0.08933327,  -0.31410187,   -1.83590478, -0.72056333, 0.39067856 ],
        [ 0.02532927,   -0.52054499,   -1.76614957, -0.67749421, 0.36103121 ],
        [ 0.14626242,   -0.73005665,   -1.69999763, -0.63557134, 0.33279809 ],
        [ 0.27445757,   -0.94390715,   -1.63759488, -0.59465565, 0.3060617 ],
        [ 0.41095526,   -1.16335152,   -1.57904084, -0.55463385, 0.28088368 ],
        [ 0.55705437,   -1.38995726,   -1.52432006, -0.51536107, 0.25727283 ],
        [ 0.71296677,   -1.62357894,   -1.47375823, -0.47699667, 0.23538139 ],
        [ 0.87465887,   -1.85813044,   -1.42858365, -0.44056301, 0.21575249 ],
        [ 1.04364057,   -2.09601936,   -1.38800652, -0.40563725, 0.19805296 ],
        [ 1.22105365,   -2.33894591,   -1.35147824, -0.37195519, 0.18204987 ],
        [ 1.40824148,   -2.58875965,   -1.31850979, -0.33927072, 0.16753326 ],
        [ 1.6063255,    -2.84690808,   -1.28873595, -0.30742453, 0.15434564 ],
        [ 1.81672991,   -3.11515449,   -1.26180806, -0.27624553, 0.14233558 ],
        [ 2.0406411,    -3.39488738,   -1.2374635,  -0.24563224, 0.13138929 ],
        [ 2.2795411,    -3.68781449,   -1.21544577, -0.21546817, 0.12139543 ],
        [ 2.53477084,   -3.99542188,   -1.19554767, -0.18568006, 0.11226493 ],
        [ 2.80831169,   -4.3199209,    -1.17754543, -0.15614491, 0.10390106 ],
        [ 3.1018576,    -4.66312722,   -1.16127075, -0.12680024, 0.09623266 ],
        [ 3.41747994,   -5.02726091,   -1.14655786, -0.09756767, 0.08919006 ],
        [ 3.75724114,   -5.41449234,   -1.1332649,  -0.06839321, 0.08271513 ],
        [ 4.12359611,   -5.82740438,   -1.12125405, -0.03921113, 0.07675193 ],
        [ 4.51919168,   -6.26875693,   -1.11040138, -0.00996494, 0.07125142 ],
        [ 4.94726728,   -6.74192892,   -1.10058651, 0.01941908,  0.06616615 ],
        [ 5.41105477,   -7.25024633,   -1.09170826, 0.04898581,  0.0614583 ],
        [ 5.91437968,   -7.79764612,   -1.08366961, 0.07879167,  0.05709179 ],
        [ 6.46166169,   -8.38866409,   -1.076381,   0.10889751,  0.05303384 ],
        [ 7.05751537,   -9.02799905,   -1.06976576, 0.1393435,   0.04925783 ],
        [ 7.7073512,    -9.72115794,   -1.06375207, 0.17017962,  0.04573875 ],
        [ 8.41717632,   -10.47423359,  -1.05827645, 0.20145005,  0.04245498 ],
        [ 9.19319536,   -11.29348117,  -1.05328572, 0.23317515,  0.03938954 ],
        [ 10.04261118,  -12.18616553,  -1.04873042, 0.26538492,  0.03652607 ],
        [ 10.9726254,   -13.15950685,  -1.04457034, 0.29807832,  0.03385233 ],
        [ 11.99163891,  -14.22194204,  -1.04076736, 0.33126828,  0.03135573 ],
        [ 13.10825214,  -15.38207659,  -1.03728977, 0.36494499,  0.02902621 ],
        [ 14.33206528,  -16.64952018,  -1.03410843, 0.39910253,  0.02685391 ],
        [ 15.67307841,  -18.03426316,  -1.0311982,  0.43372032,  0.02483027 ],
        [ 17.14209155,  -19.54709292,  -1.02853624, 0.46877549,  0.02294711 ],
        [ 18.73195696,  -21.18035577,  -1.0261278,  0.50384392,  0.02121537 ],
        [ 20.48510389,  -22.97730324,  -1.02390644, 0.53957689,  0.01959339 ],
        [ 22.35075809,  -24.88565708,  -1.02192556, 0.57472591,  0.01812594 ],
        [ 24.32024994,  -26.8965477,   -1.02016401, 0.60908862,  0.01680342 ],
        [ 26.54886044,  -29.16817087,  -1.01848536, 0.64507614,  0.01552703 ],
        [ 28.89074385,  -31.55155644,  -1.01699977, 0.68006477,  0.01438364 ],
        [ 31.54144878,  -34.24539101,  -1.01558353, 0.71669041,  0.01328102 ],
        [ 34.31025135,  -37.05556739,  -1.01433698, 0.75206899,  0.01229981 ],
        [ 37.44292978,  -40.23124254,  -1.01314792, 0.78907863,  0.01135408 ],
        [ 40.68983394,  -43.51910171,  -1.01210787, 0.82455267,  0.01051863 ],
        [ 44.35886536,  -47.23068354,  -1.0111149,  0.86162896,  0.00971355 ],
        [ 48.52521418,  -51.44131079,  -1.01016852, 0.9004401,   0.00893907 ],
        [ 52.82187555,  -55.77985071,  -1.00934808, 0.93735431,  0.00826177 ],
        [ 57.69127659,  -60.69280963,  -1.00856524, 0.9759535,   0.00761014 ],
        [ 62.64968944,  -65.69197905,  -1.00789249, 1.01225201,  0.00704577 ],
        [ 68.25033816,  -71.33497936,  -1.00724956, 1.05015007,  0.00650244 ],
        [ 74.60759083,  -77.73631198,  -1.00663615, 1.08978097,  0.00598028 ],
        [ 81.00681719,  -84.17629306,  -1.00611544, 1.12658472,  0.00553398 ],
        [ 88.23596624,  -91.44780111,  -1.00561761, 1.16499625,  0.00510458 ],
        [ 96.44328491,  -99.69921797,  -1.00514246, 1.20515044,  0.00469212 ],
        [ 104.56971958, -107.86578484, -1.00474515, 1.24183448,  0.00434519 ],
        [ 113.73421604, -117.07197726, -1.0043649,  1.28008488,  0.00401133 ],
        [ 124.11956013, -127.5007111,  -1.00400159, 1.32003164,  0.00369059 ],
        [ 134.16193224, -137.58173345, -1.00370354, 1.35573282,  0.00342614 ],
        [ 145.43789741, -148.89780522, -1.00341775, 1.3928922,   0.00317139 ],
        [ 158.15625452, -161.65784209, -1.00314413, 1.4316266,   0.00292635 ],
        [ 172.57176478, -176.11673682, -1.00288258, 1.47206755,  0.00269104 ],
        [ 186.105099,   -189.6876338,  -1.0026738,  1.50717875,  0.00250241 ],
        [ 201.25384632, -204.87532757, -1.00247329, 1.54367489,  0.00232057 ]
    ],

    energy_table => [
        [ -4.50135843,  0.01059753, 0.0119047,  0.01059753, 0.0119047,     5.5706673,    -1.5997437 ],
        [ -4.13094467,  0.01606197, 0.0180168,  0.01606197, 0.0180168,     4.9739061,    -1.6250773 ],
        [ -3.84809547,  0.02205271, 0.02468676, 0.02205271, 0.02468676,    4.5112282,    -1.6498321 ],
        [ -3.43277968,  0.03506502, 0.03903089, 0.03506502, 0.03903089,    3.8174552,    -1.6978677 ],
        [ -3.12115287,  0.0495408,  0.05469588, 0.0495408,  0.05469588,    3.281948,     -1.7461686 ],
        [ -2.8760393,   0.06485148, 0.07084997, 0.06485148, 0.07084997,    2.8485895,    -1.7947396 ],
        [ -2.66373443,  0.0816598,  0.08800131, 0.0816598,  0.08800131,    2.4626388,    -1.8480313 ],
        [ -2.47998323,  0.09938691, 0.10532434, 0.09938691, 0.10532434,    2.11827,      -1.9043617 ],
        [ -2.31717138,  0.11790949, 0.12246681, 0.11790949, 0.12246681,    1.8038524,    -1.9640939 ],
        [ -2.17011852,  0.13713023, 0.13908934, 0.13713023, 0.13908934,    1.5106541,    -2.0232366 ],
        [ -2.0350765,   0.1569757,  0.15486029, 0.1569757,  0.15486029,    1.2337852,    -2.079498 ],
        [ -1.90932842,  0.1773698,  0.16943269, 0.1773698,  0.16943269,    0.96886668,   -2.1436099 ],
        [ -1.79018687,  0.19834597, 0.18251362, 0.19834597, 0.18251362,    0.70931174,   -2.2034193 ],
        [ -1.67496445,  0.22004263, 0.19380814, 0.22004263, 0.19380814,    0.45265636,   -2.2507084 ],
        [ -1.56095269,  0.24268252, 0.20294081, 0.24268252, 0.20294081,    0.19342841,   -2.288372 ],
        [ -1.44304287,  0.26702871, 0.20946308, 0.26702871, 0.20946308,    -0.078182876, -2.2973516 ],
        [ -1.31183091,  0.29475715, 0.21235764, 0.29475715, 0.21235764,    -0.37871768,  -2.2574532 ],
        [ -1.194119,    0.31968826, 0.21051441, 0.31968826, 0.21051441,    -0.64096201,  -2.1706683 ],
        [ -1.08550312,  0.34228578, 0.20498415, 0.34228578, 0.20498415,    -0.87100028,  -2.0561066 ],
        [ -0.98293389,  0.36290597, 0.1966121,  0.36290597, 0.1966121,     -1.0759077,   -1.9436982 ],
        [ -0.88560956,  0.38155431, 0.18626105, 0.38155431, 0.18626105,    -1.2600855,   -1.8543293 ],
        [ -0.79020674,  0.39876843, 0.17437199, 0.39876843, 0.17437199,    -1.4334322,   -1.7964806 ],
        [ -0.69491589,  0.41477124, 0.16136164, 0.41477124, 0.16136164,    -1.6026676,   -1.7699109 ],
        [ -0.59953598,  0.42951672, 0.14778739, 0.42951672, 0.14778739,    -1.7709014,   -1.7678586 ],
        [ -0.50309935,  0.44310416, 0.13403878, 0.44310416, 0.13403878,    -1.9417807,   -1.7820169 ],
        [ -0.4046732,   0.45562265, 0.1204395,  0.45562265, 0.1204395,     -2.1181906,   -1.805366 ],
        [ -0.30348676,  0.46713519, 0.10727203, 0.46713519, 0.10727203,    -2.3022275,   -1.8327561 ],
        [ -0.19862101,  0.47771669, 0.09474121, 0.47771669, 0.09474121,    -2.4959387,   -1.8607874 ],
        [ -0.08933327,  0.48741758, 0.08302046, 0.48741758, 0.08302046,    -2.7008441,   -1.887308 ],
        [ 0.02532927,   0.49630287, 0.07221122, 0.49630287, 0.07221122,    -2.9187383,   -1.9111693 ],
        [ 0.14626242,   0.50442558, 0.06238033, 0.50442558, 0.06238033,    -3.1512473,   -1.9318236 ],
        [ 0.27445757,   0.51184006, 0.05355061, 0.51184006, 0.05355061,    -3.4001479,   -1.9491369 ],
        [ 0.41095526,   0.51859776, 0.04571265, 0.51859776, 0.04571265,    -3.6672995,   -1.9632153 ],
        [ 0.55705437,   0.52475588, 0.03882213, 0.52475588, 0.03882213,    -3.9550661,   -1.9743547 ],
        [ 0.71296677,   0.5303273,  0.03286051, 0.5303273,  0.03286051,    -4.2636724,   -1.982847 ],
        [ 0.87465887,   0.53522372, 0.02788558, 0.53522372, 0.02788558,    -4.5848683,   -1.9890145 ],
        [ 1.04364057,   0.53957045, 0.02371464, 0.53957045, 0.02371464,    -4.9214253,   -1.9933741 ],
        [ 1.22105365,   0.54345532, 0.02021067, 0.54345532, 0.02021067,    -5.275392,    -1.996314 ],
        [ 1.40824148,   0.54695182, 0.01725865, 0.54695182, 0.01725865,    -5.6493065,   -1.9982839 ],
        [ 1.6063255,    0.55011451, 0.01476848, 0.55011451, 0.01476848,    -6.0452913,   -1.9995282 ],
        [ 1.81672991,   0.5529919,  0.01266275, 0.5529919,  0.01266275,    -6.4661022,   -2.0002601 ],
        [ 2.0406411,    0.55562007, 0.0108804,  0.55562007, 0.0108804,     -6.9140425,   -2.0006307 ],
        [ 2.2795411,    0.55803191, 0.0093687,  0.55803191, 0.0093687,     -7.3920197,   -2.0007552 ],
        [ 2.53477084,   0.56025301, 0.00808505, 0.56025301, 0.00808505,    -7.9026755,   -2.0007571 ],
        [ 2.80831169,   0.56230932, 0.00699151, 0.56230932, 0.00699151,    -8.4499607,   -2.0007008 ],
        [ 3.1018576,    0.56421948, 0.00605849, 0.56421948, 0.00605849,    -9.0372447,   -2.0005839 ],
        [ 3.41747994,   0.56600094, 0.00526047, 0.56600094, 0.00526047,    -9.6686497,   -2.0004563 ],
        [ 3.75724114,   0.56766767, 0.00457659, 0.56766767, 0.00457659,    -10.348308,   -2.0003494 ],
        [ 4.12359611,   0.56923262, 0.00398898, 0.56923262, 0.00398898,    -11.081126,   -2.0002496 ],
        [ 4.51919168,   0.57070677, 0.0034828,  0.57070677, 0.0034828,     -11.872397,   -2.0001704 ],
        [ 4.94726728,   0.57210054, 0.00304536, 0.57210054, 0.00304536,    -12.728606,   -2.0001127 ],
        [ 5.41105477,   0.5734218,  0.00266641, 0.5734218,  0.00266641,    -13.656222,   -2.0000711 ],
        [ 5.91437968,   0.57467796, 0.00233714, 0.57467796, 0.00233714,    -14.662899,   -2.0000431 ],
        [ 6.46166169,   0.5758757,  0.00205019, 0.5758757,  0.00205019,    -15.75748,    -2.0000245 ],
        [ 7.05751537,   0.57702006, 0.00180382, 0.57702006, 0.00180382,    -16.949198,   -2.0000136 ],
        [ 7.7073512,    0.5781159,  0.00158084, 0.57811689, 0.0015943637,  -18.248875,   -2.0000071 ],
        [ 8.41717632,   0.57916724, 0.00138847, 0.57921141, 0.001513041,   -19.668529,   -2.0000035 ],
        [ 9.19319536,   0.5801768,  0.00121959, 0.58037759, 0.0014695178,  -21.220569,   -2.0000016 ],
        [ 10.04261118,  0.58114738, 0.00107113, 0.58160596, 0.0014208268,  -22.919401,   -2.0000007 ],
        [ 10.9726254,   0.58208062, 0.00094057, 0.58289253, 0.0013489305,  -24.77943,    -2.0000003 ],
        [ 11.99163891,  0.5829784,  0.00082569, 0.58422788, 0.0012720867,  -26.817457,   -2.0000001 ],
        [ 13.10825214,  0.58384189, 0.00072463, 0.58560079, 0.0011832591,  -29.050684,   -2.0000001 ],
        [ 14.33206528,  0.5846723,  0.00063572, 0.58699527, 0.0010951545,  -31.49831,    -2 ],
        [ 15.67307841,  0.58547048, 0.00055756, 0.58840063, 0.00099940267, -34.180336,   -2 ],
        [ 17.14209155,  0.58623725, 0.00048888, 0.58991101, 0.0010013077,  -37.118363,   -2 ],
        [ 18.73195696,  0.58696534, 0.00042921, 0.59140122, 0.00087786228, -40.298094,   -2 ],
        [ 20.48510389,  0.58766977, 0.00037636, 0.5928411,  0.00076880486, -43.804387,   -2 ],
        [ 22.35075809,  0.5883282,  0.00033111, 0.59418537, 0.00067562807, -47.535696,   -2 ],
        [ 24.32024994,  0.58894091, 0.00029244, 0.59543501, 0.00059614875, -51.474679,   -2 ],
        [ 26.54886044,  0.58955177, 0.00025704, 0.59667972, 0.0005235106,  -55.9319,     -2 ],
        [ 28.89074385,  0.59011725, 0.00022695, 0.59783102, 0.00046186649, -60.615667,   -2 ],
        [ 31.54144878,  0.59068099, 0.0001994,  0.59897789, 0.00040551627, -65.917077,   -2 ],
        [ 34.31025135,  0.59119975, 0.00017613, 0.60003256, 0.0003579577,  -71.454682,   -2 ],
        [ 37.44292978,  0.5917169,  0.00015481, 0.60108333, 0.00031445452, -77.720039,   -2 ],
        [ 40.68983394,  0.5921895,  0.00013692, 0.60204305, 0.00027795981, -84.213847,   -2 ],
        [ 44.35886536,  0.59266067, 0.00012051, 0.60299938, 0.00024453465, -91.55191,    -2 ],
        [ 48.52521418,  0.59313035, 0.00010552, 0.60395228, 0.00021402412, -99.884608,   -2 ],
        [ 52.82187555,  0.59355601, 9.306e-05,  0.60481553, 0.00018868604, -108.47793,   -2 ],
        [ 57.69127659,  0.59398038, 8.166e-05,  0.60567585, 0.00016551598, -118.21673,   -2 ],
        [ 62.64968944,  0.59436118, 7.226e-05,  0.6064476,  0.00014642461, -128.13356,   -2 ],
        [ 68.25033816,  0.59474089, 6.364e-05,  0.60721692, 0.00012891951, -139.33486,   -2 ],
        [ 74.60759083,  0.59511948, 5.576e-05,  0.60798377, 0.00011292224, -152.04936,   -2 ],
        [ 81.00681719,  0.59545504, 4.934e-05,  0.60866332, 9.9905976e-05, -164.84781,   -2 ],
        [ 88.23596624,  0.59578969, 4.345e-05,  0.60934088, 8.7966854e-05, -179.30611,   -2 ],
        [ 96.44328491,  0.5961234,  3.807e-05,  0.61001641, 7.7052511e-05, -195.72075,   -2 ],
        [ 104.56971958, 0.59641461, 3.375e-05,  0.61060583, 6.8302504e-05, -211.97362,   -2 ],
        [ 113.73421604, 0.59670509, 2.978e-05,  0.61119367, 6.0263265e-05, -230.30261,   -2 ],
        [ 124.11956013, 0.59699482, 2.615e-05,  0.61177991, 5.2901382e-05, -251.0733,    -2 ],
        [ 134.16193224, 0.59724255, 2.329e-05,  0.61228112, 4.7104878e-05, -271.15804,   -2 ],
        [ 145.43789741, 0.59748971, 2.065e-05,  0.61278114, 4.1761016e-05, -293.70997,   -2 ],
        [ 158.15625452, 0.59773631, 1.822e-05,  0.61327995, 3.6849362e-05, -319.14669,   -2 ],
        [ 172.57176478, 0.59798232, 1.6e-05,    0.61377755, 3.2349622e-05, -347.97771,   -2 ],
        [ 186.105099,   0.59818688, 1.429e-05,  0.61419127, 2.8900374e-05, -375.04438,   -2 ],
        [ 201.25384632, 0.59839104, 1.272e-05,  0.61460414, 2.571197e-05,  -405.34187,   -2 ]
    ],

    impulse_table => [
        [
            -5.1062481, 13.814839,   -5.1066921,    0.025118799, -0.00089631373, -0.15481658,
            -1.6509538, -0.23932539, -0.0014866937, -0.37299636, 0.26945578,     0.9914086,
            -2.3957224
        ],
        [
            -4.9618449, 13.38163,    -4.9623963,    0.026718124, -0.001035556, -0.15387535,
            -1.6500126, -0.31152838, -0.0017176513, -0.37205513, 0.26945568,   0.98989413,
            -2.3478719
        ],
        [
            -4.60517,   12.311609,  -4.6061117,    0.030991737, -0.0014793656, -0.15087535,
            -1.6470126, -0.4898707, -0.0024537876, -0.36905513, 0.26945506,    0.98491037,
            -2.2325772
        ],
        [
            -4.1997049, 11.095225,   -4.2014354,    0.036368919, -0.0022190485, -0.14587535,
            -1.6420126, -0.69261958, -0.0036806814, -0.36405513, 0.26945279,    0.97620855,
            -2.1077828
        ],
        [
            -3.9120228, 10.232202,   -3.9146885,    0.040449861, -0.0029587313, -0.14087535,
            -1.6370126, -0.83649197, -0.0049075752, -0.35905512, 0.26944852,    0.96715272,
            -2.0245096
        ],
        [
            -3.5065577, 9.0158998,  -3.5114598,    0.046368695, -0.0044380969, -0.13087535,
            -1.6270126, -1.0393531, -0.0073613628, -0.34905512, 0.2694311,     0.94832075,
            -1.9170796
        ],
        [
            -2.9957321, 7.4839032,  -3.0063079,   0.053390473, -0.0073968276, -0.11087539,
            -1.6070126, -1.2954299, -0.012268937, -0.32905512, 0.26934128,    0.90898267,
            -1.8051139
        ],
        [
            -2.3025849, 5.408721,   -2.3327303,   0.05901554,  -0.014793498, -0.060881665,
            -1.5570159, -1.6479448, -0.024537683, -0.27905511, 0.26854444,   0.80835364,
            -1.7225955
        ],
        [
            -1.8971198, 4.2036065,  -1.9529016,   0.057758518, -0.022186433, -0.010983376,
            -1.5070694, -1.8669137, -0.036801867, -0.2290551,  0.26642699,   0.71240717,
            -1.7378921
        ],
        [
            -1.6094377, 3.3614563,  -1.6957263,   0.054243772, -0.029547044, 0.038320838,
            -1.4574381, -2.0427593, -0.049026451, -0.17949721, 0.26248145,   0.62608909,
            -1.794196
        ],
        [
            -1.3862941, 2.7238176,  -1.5069849,   0.050730636, -0.036749663, 0.085409385,
            -1.40902,   -2.2079254, -0.061054337, -0.13156585, 0.25643198,   0.55126998,
            -1.8724227
        ],
        [
            -1.2039726, 2.2193326, -1.3620127,   0.048438232, -0.043443743, 0.12699527,
            -1.3639292, -2.379955, -0.072415752, -0.08782564, 0.24829091,   0.48798481,
            -1.962643
        ],
        [
            -1.0498219, 1.8086631,  -1.247247,    0.047668677,  -0.049047705, 0.15967662,
            -1.3253728, -2.5687026, -0.082189396, -0.052285925, 0.23833527,   0.43523005,
            -2.0585871
        ],
        [
            -0.91629051, 1.4671687,  -1.1543182,   0.048045923,  -0.053131418, 0.1829537,
            -1.2955278,  -2.7756795, -0.089428703, -0.028294298, 0.22701861,   0.39154845,
            -2.156092
        ],
        [
            -0.79850748, 1.1781727,  -1.0776743,   0.04901833,   -0.055790352, 0.19899776,
            -1.2736707,  -2.9941996, -0.094028331, -0.014757661, 0.21486421,   0.35541001,
            -2.2524635
        ],
        [
            -0.69314696, 0.92986896, -1.0134586,   0.050182751,   -0.057447223, 0.21030115,
            -1.2575605,  -3.2147261, -0.096690249, -0.0085559132, 0.20237816,   0.32541511,
            -2.3460635
        ],
        [
            -0.59783678, 0.71366529, -0.95890738,  0.051325195,   -0.058494621, 0.218624,
            -1.2452487,  -3.4300442, -0.098187485, -0.0063569612, 0.18999633,   0.30037358,
            -2.4359876
        ],
        [
            -0.51082541, 0.52318556, -0.91199411,  0.052356525,   -0.059184565, 0.22503421,
            -1.2354715,  -3.6362735, -0.099031799, -0.0065098571, 0.17806136,   0.27931377,
            -2.5218155
        ],
        [
            -0.4307827, 0.35362061, -0.87120547,  0.053252026,   -0.059660671, 0.23016177,
            -1.227458,  -3.831857,  -0.099509943, -0.0079693454, 0.16682002,   0.26145971,
            -2.6034314
        ],
        [
            -0.35667473, 0.20128867, -0.83539301,  0.054015561,   -0.060003599, 0.2343872,
            -1.2207307,  -4.0165188, -0.099777181, -0.0098179221, 0.15643184,   0.24619829,
            -2.6809009
        ],
        [
            -0.28768186, 0.06333169, -0.8036727,   0.054661624,  -0.060259846, 0.23795008,
            -1.2149787,  -4.1906135, -0.099918972, -0.012266883, 0.14698264,   0.23304674,
            -2.7543922
        ],
        [
            -0.22314333, -0.062498374, -0.77535513,  0.055207244,  -0.06045732, 0.24100848,
            -1.2099921,  -4.3547746,   -0.099984058, -0.014931916, 0.13849951,  0.22162454,
            -2.824125
        ],
        [
            -0.16251871, -0.17799282, -0.74989576, 0.055668596,  -0.060613466, 0.24367115,
            -1.2056208,  -4.5097339,  -0.10000116, -0.017618706, 0.1309651,    0.21163062,
            -2.8903409
        ],
        [
            -0.1053603, -0.2845953, -0.72685972,  0.056059851, -0.060739623, 0.24601581,
            -1.2017552, -4.6562323, -0.099987801, -0.02061315, 0.12432444,   0.20282544,
            -2.953284
        ],
        [
            -0.051293077, -0.38348259, -0.70589595,  0.056392963,  -0.060843441, 0.24809998,
            -1.1983121,   -4.7949763,  -0.099955078, -0.023381426, 0.11847139,   0.19501718,
            -3.0131918
        ],
        [
            -0.0094180511, -0.45882355, -0.69018915,  0.056627681,  -0.06091505, 0.24963194,
            -1.1957802,    -4.9024525,  -0.099919218, -0.025593038, 0.11421545,  0.18929911,
            -3.0597883
        ],
        [
            0.095310397, -0.64275677, -0.6528194,  0.057133628,  -0.061066752, 0.25317903,
            -1.1899317,  -5.1709415,  -0.09980043, -0.031580041, 0.10456039,   0.1761659,
            -3.1768697
        ],
        [
            0.18232177, -0.79102535, -0.62371294,  0.057475974,  -0.061169063, 0.25584805,
            -1.1855623, -5.3933018,  -0.099680679, -0.036803968, 0.097513966,  0.16641124,
            -3.2745069
        ],
        [
            0.26236448, -0.92407,   -0.59837129,  0.057737496,  -0.061248465, 0.2581046,
            -1.1819067, -5.5970212, -0.099560748, -0.041779447, 0.091721112,  0.15826613,
            -3.3644382
        ],
        [
            0.40546533, -1.1546766, -0.55617849,  0.058100034, -0.061363668, 0.26171744,
            -1.176171,  -5.9587591, -0.099336231, -0.05074782, 0.082778724,  0.14543732,
            -3.5251188
        ],
        [
            0.53062847, -1.349553,  -0.52221277,  0.058329604,  -0.061443352, 0.2644861,
            -1.1719174, -6.2721935, -0.099139457, -0.059665897, 0.076204779,  0.13578221,
            -3.6652335
        ],
        [
            0.6931474, -1.5943106, -0.48168754,  0.058537868,  -0.061525923, 0.26761101,
            -1.16733,  -6.6747396, -0.098894154, -0.070797616, 0.069063455,  0.12504784,
            -3.8462051
        ],
        [
            0.99325199, -2.0257935, -0.41574192,  0.0587482,    -0.061637122, 0.27222726,
            -1.1611583, -7.4047415, -0.098491614, -0.092617626, 0.058976627,  0.10937885,
            -4.1767761
        ],
        [
            1.0986125,  -2.1719905, -0.39489431,  0.058787734,  -0.061667381, 0.27355259,
            -1.1595646, -7.6570636, -0.098368024, -0.099950221, 0.056156596,  0.10487652,
            -4.2916372
        ],
        [
            1.2089606,  -2.3225884, -0.37416286, 0.05881677,  -0.061695477, 0.27480166,
            -1.1581537, -7.9192409, -0.09824884, -0.10821055, 0.053516374,  0.10060723,
            -4.4112713
        ],
        [
            1.3862946,  -2.5597831, -0.34296493,  0.058844266, -0.061734239, 0.27654608,
            -1.1563521, -8.3363332, -0.098077787, -0.12152683, 0.049836332,  0.094562577,
            -4.6021494
        ],
        [
            1.6094381,  -2.8509188, -0.30694441,  0.05885731, -0.061774201, 0.27835049,
            -1.1547298, -8.8543227, -0.097896006, -0.1392543, 0.04599718,   0.088129119,
            -4.8400756
        ],
        [
            1.7917597,  -3.0836099, -0.27981621,  0.05885756,  -0.061800983, 0.27955685,
            -1.1538044, -9.2724291, -0.097771967, -0.15332775, 0.043374064,  0.0836512,
            -5.0327783
        ],
        [
            1.9459104,  -3.2772025, -0.25828596,  0.058853631, -0.06182018, 0.28041982,
            -1.1532284, -9.6226909, -0.097681857, -0.16567298, 0.041444411, 0.080310828,
            -5.194638
        ],
        [
            2.0794418,  -3.4428273, -0.24056793,  0.058848616, -0.061834781, 0.28106842,
            -1.1528443, -9.9239053, -0.097613714, -0.176126,   0.039951273,  0.077697462,
            -5.3341295
        ],
        [
            2.3027871,  -3.716046,  -0.21265659,  0.058838663, -0.061856493, 0.28197962,
            -1.1524788, -10.423588, -0.097517806, -0.19271887, 0.037761684,  0.073002089,
            -5.5756081
        ],
        [
            2.9958229, -4.539701,  -0.13714231,  0.058808567, -0.061896272, 0.28380264,
            -1.151599, -11.947165, -0.097321977, -0.24714354, 0.032713901,  0.064318322,
            -6.2820113
        ],
        [
            3.9121033,  -5.589574, -0.055789995, 0.058784531, -0.061928301, 0.28491148,
            -1.1519273, -13.91717, -0.097205891, -0.31935155, 0.028427923,  0.056457801,
            -7.2080232
        ],
        [
            4.6054578, -6.3644549, -0.0038655713, 0.058775689, -0.061936281, 0.28527784,
            -1.151853, -15.385971, -0.097165558,  -0.37405383, 0.02616526,   0.052142155,
            -7.9051507
        ],
        [
            5.2984171,  -7.1271666, 0.042003136,  0.058771179, -0.06194074, 0.28546123,
            -1.1518594, -16.841229, -0.097145686, -0.42874109, 0.024403177, 0.04871588,
            -8.600189
        ],
        [
            6.2146298,  -8.1223831, 0.095582631,  0.05876848,  -0.061943409, 0.28557135,
            -1.1518644, -18.751501, -0.097133751, -0.50105089, 0.022570845,  0.045106892,
            -9.5177675
        ],
        [
            6.9079043,  -8.8678354, 0.1319077,    0.058767594, -0.061944297, 0.28560807,
            -1.1518662, -20.189237, -0.097129771, -0.55576483, 0.021443206,  0.042869407,
            -10.211536
        ],
        [
            7.6009464,  -9.6079215, 0.16528429,   0.058767159, -0.061944739, 0.28562643,
            -1.1518672, -21.621413, -0.097127783, -0.61045906, 0.020476678,  0.040945009,
            -10.90484
        ],
        [
            8.5173792,  -10.580241, 0.2056829,    0.058766903, -0.061945009, 0.28563745,
            -1.1518678, -23.509037, -0.097126589, -0.68278146, 0.019385222,  0.038767173,
            -11.821439
        ],
        [
            9.2103957,  -11.311597, 0.23385213,   0.058766819, -0.061945086, 0.28564112,
            -1.1518679, -24.932669, -0.097126192, -0.73747165, 0.01866979,   0.037337966,
            -12.514515
        ],
        [
            10.819823,  -12.999846, 0.29287443,   0.058766751, -0.061945164, 0.28564406,
            -1.1518681, -28.22896,  -0.097006648, -0.85023122, 0.017278569,  0.034556822,
            -14.123992
        ],
        [
            11.513033,  -13.723424, 0.31599608,   0.058766742, -0.061945178, 0.28564443,
            -1.1518681, -29.645269, -0.096419425, -0.88062918, 0.016769755,  0.033539353,
            -14.817208
        ],
        [
            13.122452,  -15.396806, 0.36535697,   0.058766733, -0.061945189, 0.28564472,
            -1.1518681, -32.927166, -0.094039239, -0.94784301, 0.01574492,   0.031489809,
            -16.426633
        ],
        [
            13.81566,   -16.115177, 0.38501212,   0.058766731, -0.061945167, 0.28564476,
            -1.1518681, -34.338425, -0.092804324, -0.97552363, 0.015358586,  0.030717157,
            -17.119842
        ],
        [
            16.118283,  -18.493164, 0.44463954,   0.058766725, -0.061945169, 0.28564482,
            -1.1518681, -39.018187, -0.088485153, -1.0628874,  0.014255705,  0.028511408,
            -19.422465
        ],
        [
            18.420702,  -20.860901, 0.49719155,   0.058766719, -0.061945255, 0.28564486,
            -1.1518681, -43.687728, -0.084296995, -1.1443628,  0.013362552,  0.026725103,
            -21.724885
        ]
    ],

    tail_shock_table => [
        [ 7.5150663, -0.68914,    -0.041369999, -0.041369999 ],
        [ 7.6003564, -0.69372713, -0.047446316, -0.034954904 ],
        [ 8.5171335, -0.74134728, -0.057683411, -0.01986323 ],
        [ 9.2102694, -0.77557293, -0.059790188, -0.014648565 ],
        [ 10.819796, -0.85022926, -0.060629217, -0.0079091622 ],
        [ 11.513019, -0.8806282,  -0.060262179, -0.0061611827 ],
        [ 13.122449, -0.94784281, -0.058774533, -0.003450548 ],
        [ 13.815658, -0.97552353, -0.058002707, -0.0026599697 ],
        [ 16.118282, -1.062921,   -0.055306188, -0.0009542178 ]
    ],

    blast_info => [
        0.16087533, 1.6570119,  0.37902484, -0.24537877, 2.0036232,  0.0018645358, 3.7166874,   -0.14793653,
        0.48897016, 0.38780118, 0.64817471, -0.17018313, 0.31949081, 0.26945596,   0.036729203, -0.038715727,
        1834.8,     -0.68914,   32268.334,  -0.82991671
    ],

};
1;
