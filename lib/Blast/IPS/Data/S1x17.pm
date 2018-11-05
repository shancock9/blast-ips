package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=1.17
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S1x17'} = {

    table_name  => 'S1x17',
    symmetry    => 2,
    gamma       => 1.17,
    data_source => 'S8000_G1x17/moc2_from_moc_r100/',

    shock_table_info => [ 2, 1.17, 1.28e-06, 8000, 2e-07, 10, 3.31e-07, 4.54e-07, 5e-07, 1885.2, -0.50399 ],

    shock_table => [
        [ -4.85499195,  11.94676413,   -2.99997904, -4.85604989, 0.99841226 ],
        [ -4.80516301,  11.7972789,    -2.99997566, -4.8063031,  0.9982889 ],
        [ -4.36207739,  10.46804393,   -2.99991046, -4.36429463, 0.99667054 ],
        [ -3.86367605,  8.97294338,    -2.99960788, -3.86836427, 0.99295202 ],
        [ -3.67197461,  8.39794274,    -2.99927863, -3.67822935, 0.99059053 ],
        [ -3.35474904,  7.44664266,    -2.99815987, -3.36483271, 0.98480656 ],
        [ -3.12710359,  6.76430629,    -2.99636542, -3.14131676, 0.97855197 ],
        [ -2.92834772,  6.1690245,     -2.99342721, -2.94753657, 0.97099754 ],
        [ -2.75736168,  5.65753128,    -2.98907576, -2.78221407, 0.96238136 ],
        [ -2.60160426,  5.19241989,    -2.98269918, -2.63306942, 0.95230807 ],
        [ -2.46085152,  4.77317611,    -2.97387974, -2.49980361, 0.94090128 ],
        [ -2.33062443,  4.38662736,    -2.96193402, -2.37808851, 0.92795133 ],
        [ -2.2075182,   4.02292355,    -2.94596491, -2.26473432, 0.91316826 ],
        [ -2.08916048,  3.6754311,     -2.92486578, -2.15762818, 0.89622356 ],
        [ -1.97337819,  3.33830986,    -2.8972052,  -2.05496284, 0.87667676 ],
        [ -1.85558034,  2.99910632,    -2.86025342, -1.95302716, 0.85341823 ],
        [ -1.72692091,  2.63433304,    -2.80788825, -1.84508368, 0.8237997 ],
        [ -1.60215953,  2.28787356,    -2.74382599, -1.74432421, 0.79072062 ],
        [ -1.48828082,  1.9792992,     -2.67371457, -1.65617297, 0.75688202 ],
        [ -1.38170338,  1.69825967,    -2.59879771, -1.57732284, 0.72237479 ],
        [ -1.27991798,  1.43768632,    -2.52026384, -1.50556479, 0.6873149 ],
        [ -1.1814658,   1.19350291,    -2.43956645, -1.43962697, 0.65199303 ],
        [ -1.08438336,  0.9606423,     -2.35732352, -1.37805399, 0.61639966 ],
        [ -0.98774882,  0.7368372,     -2.27468433, -1.3202086,  0.58082486 ],
        [ -0.89053645,  0.5197196,     -2.19246989, -1.26546918, 0.54547703 ],
        [ -0.79188636,  0.30745397,    -2.11144268, -1.2133886,  0.51059348 ],
        [ -0.69101966,  0.09850789,    -2.03226956, -1.16362476, 0.47641476 ],
        [ -0.58710673,  -0.10863545,   -1.95544868, -1.1158661,  0.44314214 ],
        [ -0.47924864,  -0.31549686,   -1.88134064, -1.0698284,  0.41093976 ],
        [ -0.36671362,  -0.52315853,   -1.81035385, -1.02534962, 0.38000559 ],
        [ -0.24850696,  -0.73308615,   -1.74267542, -0.98220671, 0.35044666 ],
        [ -0.12376666,  -0.94639292,   -1.67854337, -0.94027419, 0.32239175 ],
        [ 0.00851011,   -1.16434506,   -1.61808701, -0.89941568, 0.29591859 ],
        [ 0.14948681,   -1.38837065,   -1.56135287, -0.85948854, 0.2710623 ],
        [ 0.30043581,   -1.61996097,   -1.5083484,  -0.82036612, 0.24783496 ],
        [ 0.46112013,   -1.85830844,   -1.45950971, -0.78230497, 0.22642975 ],
        [ 0.62748674,   -2.09740707,   -1.41593581, -0.74626323, 0.20732629 ],
        [ 0.80153911,   -2.34035612,   -1.37673578, -0.71171194, 0.19012862 ],
        [ 0.98450649,   -2.58893663,   -1.34138594, -0.67838137, 0.17460034 ],
        [ 1.17721883,   -2.84428726,   -1.30951564, -0.64612015, 0.16057185 ],
        [ 1.38207391,   -3.10950645,   -1.28058762, -0.61456837, 0.14779969 ],
        [ 1.59893515,   -3.38431149,   -1.25449898, -0.58380321, 0.13623287 ],
        [ 1.83039089,   -3.67185857,   -1.23082996, -0.55352457, 0.12568054 ],
        [ 2.07692106,   -3.97258689,   -1.20944817, -0.52375499, 0.11608054 ],
        [ 2.34063229,   -4.28890809,   -1.1900902,  -0.49432977, 0.10731314 ],
        [ 2.62336263,   -4.6228349,    -1.17256436, -0.46515306, 0.09929152 ],
        [ 2.926734,     -4.97608475,   -1.15671826, -0.43617425, 0.0919479 ],
        [ 3.25261405,   -5.35063596,   -1.1424029,  -0.40733594, 0.08521732 ],
        [ 3.60326245,   -5.74888416,   -1.12947175, -0.37856678, 0.07903701 ],
        [ 3.98133037,   -6.17362705,   -1.11778683, -0.34978898, 0.07334899 ],
        [ 4.38926036,   -6.62738968,   -1.10723571, -0.32096473, 0.06810834 ],
        [ 4.83048804,   -7.11376525,   -1.09769689, -0.2920097,  0.06326619 ],
        [ 5.30844173,   -7.63628967,   -1.0890695,  -0.26287054, 0.0587841 ],
        [ 5.82714365,   -8.19910495,   -1.08125809, -0.23348559, 0.05462641 ],
        [ 6.39081057,   -8.80651906,   -1.07417987, -0.20381152, 0.05076379 ],
        [ 7.00465484,   -9.46386618,   -1.06775414, -0.17378311, 0.04716736 ],
        [ 7.67388505,   -10.17642632,  -1.06191483, -0.14336823, 0.04381538 ],
        [ 8.40490699,   -10.95070704,  -1.05659785, -0.11251201, 0.040686 ],
        [ 9.20412435,   -11.79316309,  -1.05175149, -0.08119374, 0.03776316 ],
        [ 10.07873933,  -12.71104811,  -1.04732875, -0.04939133, 0.0350321 ],
        [ 11.03655317,  -13.71220082,  -1.04328862, -0.01709166, 0.0324801 ],
        [ 12.08596652,  -14.80504552,  -1.03959537, 0.01570918,  0.03009611 ],
        [ 13.23617969,  -15.99880047,  -1.0362172,  0.04901161,  0.02787003 ],
        [ 14.49659282,  -17.30285636,  -1.03312731, 0.08279555,  0.02579357 ],
        [ 15.87740594,  -18.72740232,  -1.03030124, 0.11703968,  0.02385868 ],
        [ 17.39061908,  -20.2844523,   -1.02771518, 0.15174241,  0.02205653 ],
        [ 18.98971805,  -21.9259872,   -1.02542601, 0.18568214,  0.02043467 ],
        [ 20.6772994,   -23.65471108,  -1.02339026, 0.2189,      0.01897026 ],
        [ 22.58478882,  -25.60491442,  -1.02145155, 0.25370013,  0.01755535 ],
        [ 24.5595727,   -27.6203424,   -1.01975825, 0.28709959,  0.01630244 ],
        [ 26.78818752,  -29.89114027,  -1.01814374, 0.32204732,  0.01509229 ],
        [ 29.31500748,  -32.46180683,  -1.01660682, 0.358668,    0.01392543 ],
        [ 31.91645623,  -35.10467768,  -1.01527584, 0.39352926,  0.01290266 ],
        [ 34.85921278,  -38.0904662,   -1.01400702, 0.43000823,  0.01191648 ],
        [ 37.84978674,  -41.12125582,  -1.01291745, 0.46433291,  0.01106052 ],
        [ 41.22049446,  -44.53370148,  -1.01187686, 0.50018623,  0.01023477 ],
        [ 45.03737332,  -48.39396396,  -1.01088459, 0.53769161,  0.00943949 ],
        [ 48.86973501,  -52.26638255,  -1.01004262, 0.57252874,  0.00875835 ],
        [ 53.18760984,  -56.62583432,  -1.00923789, 0.60889093,  0.00810166 ],
        [ 58.07525824,  -61.55670325,  -1.00846998, 0.6469017,   0.00746957 ],
        [ 62.89994972,  -66.42067009,  -1.00782791, 0.68162262,  0.00693681 ],
        [ 68.32376635,  -71.8852319,   -1.00721341, 0.71781715,  0.0064231 ],
        [ 74.44875417,  -78.05255246,  -1.00662619, 0.75560301,  0.00592857 ],
        [ 81.39952893,  -85.04738005,  -1.00606598, 0.79511272,  0.0054533 ],
        [ 88.12995694,  -91.81704988,  -1.00560708, 0.83046341,  0.00506132 ],
        [ 95.69714134,  -99.42495618,  -1.00516765, 0.86729687,  0.00468363 ],
        [ 104.24410632, -108.0142425,  -1.00474749, 0.90573248,  0.00432026 ],
        [ 113.94565134, -117.75984283, -1.00434644, 0.94590475,  0.00397128 ],
        [ 123.06617936, -126.91851749, -1.00402671, 0.98081753,  0.00369148 ],
        [ 133.28901446, -137.18090803, -1.00372003, 1.01714094,  0.00342174 ],
        [ 144.79745071, -148.73041951, -1.00342631, 1.05498595,  0.00316208 ],
        [ 157.81437729, -161.79006586, -1.00314543, 1.09447728,  0.00291252 ],
        [ 172.61298358, -176.63317728, -1.00287731, 1.13575574,  0.00267309 ],
        [ 185.95900362, -190.01619355, -1.00267192, 1.1701719,   0.00248885 ],
        [ 200.87022737, -204.96575006, -1.00247459, 1.20592569,  0.00231113 ]
    ],

    energy_table => [
        [ -4.85499195,  0.17076127, 0.07441105, 0.17076127, 0.07441105,    5.7068013,   -1.6010375 ],
        [ -4.80516301,  0.17450956, 0.07604127, 0.17450956, 0.07604127,    5.6269484,   -1.6042902 ],
        [ -4.36207739,  0.2116629,  0.09216583, 0.2116629,  0.09216583,    4.9092214,   -1.6397556 ],
        [ -3.86367605,  0.26290766, 0.11420023, 0.26290766, 0.11420023,    4.0807988,   -1.691877 ],
        [ -3.67197461,  0.28571525, 0.12385869, 0.28571525, 0.12385869,    3.7542723,   -1.7187243 ],
        [ -3.35474904,  0.32771703, 0.141215,   0.32771703, 0.141215,      3.2009559,   -1.7773981 ],
        [ -3.12710359,  0.36136623, 0.15449633, 0.36136623, 0.15449633,    2.790922,    -1.8324344 ],
        [ -2.92834772,  0.39324721, 0.16630493, 0.39324721, 0.16630493,    2.4212939,   -1.8913487 ],
        [ -2.75736168,  0.42253801, 0.1762308,  0.42253801, 0.1762308,     2.0932456,   -1.9529287 ],
        [ -2.60160426,  0.45065413, 0.18465475, 0.45065413, 0.18465475,    1.7841872,   -2.0209011 ],
        [ -2.46085152,  0.47712734, 0.19132207, 0.47712734, 0.19132207,    1.4950754,   -2.092647 ],
        [ -2.33062443,  0.50237941, 0.1962585,  0.50237941, 0.1962585,     1.2179044,   -2.1795474 ],
        [ -2.2075182,   0.52675169, 0.19941697, 0.52675169, 0.19941697,    0.9436324,   -2.2712821 ],
        [ -2.08916048,  0.5504478,  0.20067542, 0.5504478,  0.20067542,    0.66987541,  -2.3551503 ],
        [ -1.97337819,  0.57365639, 0.19986006, 0.57365639, 0.19986006,    0.39241241,  -2.4528229 ],
        [ -1.85558034,  0.59703481, 0.19664124, 0.59703481, 0.19664124,    0.096715339, -2.5703852 ],
        [ -1.72692091,  0.62194997, 0.19013051, 0.62194997, 0.19013051,    -0.24244525, -2.6622352 ],
        [ -1.60215953,  0.64511881, 0.18079537, 0.64511881, 0.18079537,    -0.57774988, -2.702085 ],
        [ -1.48828082,  0.66510598, 0.16987705, 0.66510598, 0.16987705,    -0.88696885, -2.6291348 ],
        [ -1.38170338,  0.68258779, 0.1579409,  0.68258779, 0.1579409,     -1.1585772,  -2.3644722 ],
        [ -1.27991798,  0.69803395, 0.14541972, 0.69803395, 0.14541972,    -1.3813607,  -2.0419294 ],
        [ -1.1814658,   0.71172826, 0.1327119,  0.71172826, 0.1327119,     -1.568411,   -1.7770278 ],
        [ -1.08438336,  0.72399642, 0.12003491, 0.72399642, 0.12003491,    -1.7291655,  -1.6282642 ],
        [ -0.98774882,  0.73499452, 0.1076563,  0.73499452, 0.1076563,     -1.8838581,  -1.6096171 ],
        [ -0.89053645,  0.74487642, 0.09576557, 0.74487642, 0.09576557,    -2.041195,   -1.6421264 ],
        [ -0.79188636,  0.75376142, 0.08451952, 0.75376142, 0.08451952,    -2.2055575,  -1.6921014 ],
        [ -0.69101966,  0.7617489,  0.07403741, 0.7617489,  0.07403741,    -2.3789141,  -1.7436002 ],
        [ -0.58710673,  0.76893098, 0.0643927,  0.76893098, 0.0643927,     -2.5627655,  -1.791688 ],
        [ -0.47924864,  0.77539204, 0.0556206,  0.77539204, 0.0556206,     -2.7585218,  -1.8343242 ],
        [ -0.36671362,  0.78119631, 0.04774218, 0.78119631, 0.04774218,    -2.9672231,  -1.8708363 ],
        [ -0.24850696,  0.78641363, 0.04073611, 0.78641363, 0.04073611,    -3.1903912,  -1.9013587 ],
        [ -0.12376666,  0.79109855, 0.03457326, 0.79109855, 0.03457326,    -3.4293325,  -1.9263136 ],
        [ 0.00851011,   0.7953046,  0.02920328, 0.7953046,  0.02920328,    -3.6856544,  -1.9462961 ],
        [ 0.14948681,   0.79908269, 0.02456255, 0.79908269, 0.02456255,    -3.9613179,  -1.9619679 ],
        [ 0.30043581,   0.8024786,  0.02058252, 0.8024786,  0.02058252,    -4.2585392,  -1.9739871 ],
        [ 0.46112013,   0.80550538, 0.0172233,  0.80550538, 0.0172233,     -4.5765779,  -1.9829006 ],
        [ 0.62748674,   0.80813336, 0.01447743, 0.80813336, 0.01447743,    -4.9070906,  -1.9891599 ],
        [ 0.80153911,   0.81044809, 0.01221015, 0.81044809, 0.01221015,    -5.2537651,  -1.9935244 ],
        [ 0.98450649,   0.81250337, 0.01033025, 0.81250337, 0.01033025,    -5.6188501,  -1.9964669 ],
        [ 1.17721883,   0.81433773, 0.00876877, 0.81433773, 0.00876877,    -6.0038195,  -1.998354 ],
        [ 1.38207391,   0.81599471, 0.00746024, 0.81599471, 0.00746024,    -6.4133485,  -1.9995407 ],
        [ 1.59893515,   0.81748946, 0.00636824, 0.81748946, 0.00636824,    -6.8470689,  -2.0002214 ],
        [ 1.83039089,   0.81885285, 0.00544908, 0.81885285, 0.00544908,    -7.3100887,  -2.0005672 ],
        [ 2.07692106,   0.82009724, 0.00467662, 0.82009724, 0.00467662,    -7.803315,   -2.0006777 ],
        [ 2.34063229,   0.8212411,  0.00402409, 0.8212411,  0.00402409,    -8.3309174,  -2.0006578 ],
        [ 2.62336263,   0.82229762, 0.00347118, 0.82229762, 0.00347118,    -8.8965566,  -2.0005975 ],
        [ 2.926734,     0.82327671, 0.00300181, 0.82327671, 0.00300181,    -9.5034696,  -2.0005025 ],
        [ 3.25261405,   0.82418733, 0.00260236, 0.82418733, 0.00260236,    -10.155373,  -2.0003876 ],
        [ 3.60326245,   0.82503776, 0.00226132, 0.82503776, 0.00226132,    -10.856786,  -2.0002883 ],
        [ 3.98133037,   0.82583535, 0.0019691,  0.82583535, 0.0019691,     -11.613013,  -2.0002051 ],
        [ 4.38926036,   0.82658547, 0.0017181,  0.82658547, 0.0017181,     -12.428941,  -2.0001384 ],
        [ 4.83048804,   0.82729398, 0.00150162, 0.82729398, 0.00150162,    -13.311444,  -2.0000895 ],
        [ 5.30844173,   0.82796526, 0.00131432, 0.82796526, 0.00131432,    -14.267384,  -2.0000561 ],
        [ 5.82714365,   0.82860327, 0.00115142, 0.82860327, 0.00115142,    -15.30481,   -2.0000336 ],
        [ 6.39081057,   0.82921113, 0.00101292, 0.82921113, 0.00101292,    -16.432158,  -2.0000187 ],
        [ 7.00465484,   0.82979202, 0.00088688, 0.82979202, 0.00088688,    -17.659854,  -2.0000102 ],
        [ 7.67388505,   0.83034804, 0.00077877, 0.83034854, 0.00078948866, -18.998319,  -2.0000052 ],
        [ 8.40490699,   0.83088143, 0.00068399, 0.83091569, 0.00076218175, -20.460366,  -2.0000025 ],
        [ 9.20412435,   0.83139362, 0.00060078, 0.83151694, 0.00073898437, -22.058802,  -2.0000012 ],
        [ 10.07873933,  0.83188592, 0.00052764, 0.83215867, 0.00072557547, -23.808032,  -2.0000005 ],
        [ 11.03655317,  0.83235937, 0.00046331, 0.83283826, 0.00069570183, -25.72366,   -2.0000002 ],
        [ 12.08596652,  0.83281479, 0.00040671, 0.83354899, 0.00065989252, -27.822487,  -2.0000001 ],
        [ 13.23617969,  0.83325291, 0.00035691, 0.83429669, 0.00073898177, -30.122913,  -2.0000001 ],
        [ 14.49659282,  0.83367415, 0.00031311, 0.83516803, 0.00064706992, -32.64374,   -2.0000001 ],
        [ 15.87740594,  0.83407894, 0.00027461, 0.83600385, 0.00056654252, -35.405366,  -2.0000001 ],
        [ 17.39061908,  0.83446794, 0.00024077, 0.83680576, 0.00049596578, -38.431793,  -2.0000001 ],
        [ 18.98971805,  0.83482912, 0.000212,   0.83754927, 0.00043609643, -41.629991,  -2 ],
        [ 20.6772994,   0.83516539, 0.00018738, 0.83824061, 0.00038500946, -45.005153,  -2 ],
        [ 22.58478882,  0.83550054, 0.00016485, 0.83892889, 0.00033833319, -48.820132,  -2 ],
        [ 24.5595727,   0.83580675, 0.00014592, 0.83955707, 0.00029921468, -52.7697,    -2 ],
        [ 26.78818752,  0.83611193, 0.00012858, 0.8401826,  0.00026341553, -57.22693,   -2 ],
        [ 29.31500748,  0.83641606, 0.00011273, 0.84080542, 0.00023075903, -62.28057,   -2 ],
        [ 31.91645623,  0.83669158, 9.956e-05,  0.84136924, 0.00020365165, -67.483467,  -2 ],
        [ 34.85921278,  0.83696616, 8.75e-05,   0.84193076, 0.00017887104, -73.36898,   -2 ],
        [ 37.84978674,  0.83721246, 7.756e-05,  0.84243412, 0.00015845488, -79.350128,  -2 ],
        [ 41.22049446,  0.83745794, 6.843e-05,  0.84293557, 0.00013973459, -86.091543,  -2 ],
        [ 45.03737332,  0.8377026,  6.008e-05,  0.84343508, 0.00012262333, -93.725301,  -2 ],
        [ 48.86973501,  0.83791935, 5.328e-05,  0.84387743, 0.0001086955,  -101.39002,  -2 ],
        [ 53.18760984,  0.83813542, 4.703e-05,  0.84431821, 9.5913708e-05, -110.02577,  -2 ],
        [ 58.07525824,  0.83835079, 4.131e-05,  0.8447574,  8.4220893e-05, -119.80107,  -2 ],
        [ 62.89994972,  0.83853864, 3.672e-05,  0.84514036, 7.4838695e-05, -129.45045,  -2 ],
        [ 68.32376635,  0.83872593, 3.249e-05,  0.84552207, 6.6210306e-05, -140.29809,  -2 ],
        [ 74.44875417,  0.83891265, 2.862e-05,  0.84590253, 5.8299387e-05, -152.54806,  -2 ],
        [ 81.39952893,  0.83909879, 2.508e-05,  0.84628171, 5.107039e-05,  -166.44961,  -2 ],
        [ 88.12995694,  0.83925787, 2.229e-05,  0.84660569, 4.5390237e-05, -179.91047,  -2 ],
        [ 95.69714134,  0.83941651, 1.973e-05,  0.84692872, 4.0163679e-05, -195.04484,  -2 ],
        [ 104.24410632, 0.8395747,  1.738e-05,  0.84725078, 3.5369143e-05, -212.13877,  -2 ],
        [ 113.94565134, 0.83973243, 1.522e-05,  0.84757187, 3.0985691e-05, -231.54186,  -2 ],
        [ 123.06617936, 0.83986353, 1.358e-05,  0.8478387,  2.7631914e-05, -249.78291,  -2 ],
        [ 133.28901446, 0.8399943,  1.206e-05,  0.84810484, 2.4537389e-05, -270.22858,  -2 ],
        [ 144.79745071, 0.84012475, 1.066e-05,  0.84837028, 2.1690513e-05, -293.24546,  -2 ],
        [ 157.81437729, 0.84025486, 9.38e-06,   0.84863503, 1.9079642e-05, -319.27931,  -2 ],
        [ 172.61298358, 0.84038464, 8.21e-06,   0.84889907, 1.6693543e-05, -348.87652,  -2 ],
        [ 185.95900362, 0.84048822, 7.35e-06,   0.84910979, 1.4938881e-05, -375.56856,  -2 ],
        [ 200.87022737, 0.84059159, 6.55e-06,   0.84932005, 1.3315276e-05, -405.39101,  -2 ]
    ],

    impulse_table => [
        [
            -5.4776857, 13.814839,   -5.4781013,    0.014096226, -0.00035984991, -0.13991812,
            -1.1022062, -0.44430966, -0.0008315299, -0.30723309, 0.12506608,     0.84769566,
            -1.7073112
        ],
        [
            -4.961845,  12.267321,   -4.9627462,    0.017488799, -0.0006027652, -0.13709711,
            -1.0993852, -0.70223469, -0.0013928509, -0.30441208, 0.12506565,    0.80929936,
            -1.6452084
        ],
        [
            -4.6051701, 11.197306,   -4.6067085,   0.020128226, -0.00086109314, -0.13409711,
            -1.0963852, -0.88058088, -0.001989787, -0.30141208, 0.1250647,      0.77723453,
            -1.6063918
        ],
        [
            -4.199705,  9.980945,   -4.2025325,    0.023346463, -0.0012916397, -0.12909711,
            -1.0913852, -1.0833435, -0.0029846806, -0.29641208, 0.12506146,    0.73421533,
            -1.5683198
        ],
        [
            -3.9120229, 9.1179661,  -3.9163821,    0.025688043, -0.0017221863, -0.12409711,
            -1.0863852, -1.2272429, -0.0039795741, -0.29141207, 0.12505522,    0.69879053,
            -1.5467466
        ],
        [
            -3.5065578, 7.901844,   -3.5145794,    0.028855435, -0.0025832794, -0.11409711,
            -1.0763852, -1.4302149, -0.0059693611, -0.28141207, 0.12502963,    0.64088855,
            -1.527334
        ],
        [
            -2.9957321, 6.3707761,  -3.0130631,    0.031945944, -0.0043054657, -0.094097107,
            -1.0563852, -1.6868629, -0.0099489352, -0.26141207, 0.12489793,    0.55283349,
            -1.5307802
        ],
        [
            -2.302585,  4.3036202,  -2.3521128,  0.031799408, -0.0086109314, -0.044097101,
            -1.0063852, -2.0444796, -0.01989787, -0.21141206, 0.12374402,    0.40556955,
            -1.6261758
        ],
        [
            -1.8971199,  3.1182177,  -1.9886581,   0.027589229, -0.012916396, 0.0059028708,
            -0.95638517, -2.2773639, -0.029846805, -0.16141206, 0.1207947,    0.31041118,
            -1.769191
        ],
        [
            -1.6094378,  2.3078588,  -1.7500867,   0.022899838, -0.017221613, 0.055892995,
            -0.90639093, -2.4811164, -0.039795435, -0.1115299,  0.11568719,   0.24535772,
            -1.9308453
        ],
        [
            -1.3862942,  1.7101982,  -1.5806427,  0.019832686,  -0.021502858, 0.10513773,
            -0.85683898, -2.6972912, -0.04971418, -0.061647776, 0.10863321,   0.19990932,
            -2.0966959
        ],
        [
            -1.2039727,  1.2486216,  -1.4543934,   0.019511555,  -0.025161731, 0.14219435,
            -0.81620157, -2.9648711, -0.058640837, -0.017826694, 0.10026061,   0.16753918,
            -2.2577972
        ],
        [
            -1.049822,   0.87968099, -1.3569706,   0.020572304,   -0.026726154, 0.15813733,
            -0.79622331, -3.2946548, -0.062630153, 0.00019907578, 0.091325257,  0.14397474,
            -2.4093735
        ],
        [
            -0.9162906,  0.5764629,  -1.2796371,   0.021674901,  -0.027181031, 0.16644419,
            -0.78541256, -3.6271865, -0.063412323, 0.0030248464, 0.082492088,  0.12641218,
            -2.5496305
        ],
        [
            -0.79850757, 0.32145195, -1.216777,    0.022564572,  -0.027364723, 0.17214582,
            -0.77804843, -3.9344521, -0.063468896, 0.0019948434, 0.074240376,  0.11301095,
            -2.6785005
        ],
        [
            -0.69314705, 0.10283306, -1.164639,    0.023248897,   -0.027466729, 0.17654008,
            -0.77250798, -4.2138533, -0.063359266, 0.00012014942, 0.066856183,  0.10255376,
            -2.7967289
        ],
        [
            -0.59783687, -0.087611968, -1.120639,    0.023772515,   -0.027535322, 0.18010873,
            -0.76813305, -4.4679504,   -0.063210723, -0.0023203118, 0.060465041,  0.094224343,
            -2.9053518
        ],
        [
            -0.5108255,  -0.25575845, -1.0829489,   0.024176086,   -0.027586577, 0.18309149,
            -0.76457829, -4.6999899,  -0.063058365, -0.0049273978, 0.055075546,  0.087465577,
            -3.0054414
        ],
        [
            -0.43078279, -0.40591267, -1.0502449,   0.024490647,   -0.027627221, 0.18563188,
            -0.76163349, -4.912963,   -0.062912875, -0.0073933616, 0.050620306,  0.081889755,
            -3.0979973
        ],
        [
            -0.35667482, -0.54130207, -1.021548,    0.024738831,  -0.027660627, 0.18782557,
            -0.75915807, -5.1094326,  -0.062777366, -0.010183374, 0.046982084,  0.077221697,
            -3.183909
        ],
        [
            -0.30472291, -0.63424693, -1.002289,    0.024890792,  -0.027682046, 0.18927864,
            -0.75755652, -5.2466454,  -0.062683031, -0.012102531, 0.04469031,   0.074207578,
            -3.2441815
        ],
        [
            -0.22314342, -0.77711313, -0.97339397,  0.025097046,  -0.02771284, 0.19142796,
            -0.7552473,  -5.46111,    -0.062537247, -0.015143311, 0.04146048,  0.069864881,
            -3.3387925
        ],
        [
            -0.1625188,  -0.88097698, -0.95292897,  0.02522765,   -0.027733745, 0.19292595,
            -0.75368373, -5.6196193,  -0.062431398, -0.017645707, 0.039316042,  0.066919474,
            -3.4090066
        ],
        [
            -0.10536039, -0.97720692, -0.93437593,  0.025335268,  -0.027752086, 0.19426503,
            -0.75232114, -5.7683431,  -0.062333946, -0.019944324, 0.037469357,  0.064342369,
            -3.4750873
        ],
        [
            -0.051293167, -1.0667958, -0.91745524, 0.025424731,  -0.027768323, 0.1954692,
            -0.7511255,   -5.908353,  -0.06224406, -0.022155172, 0.035862731,  0.062068879,
            -3.5374599
        ],
        [
            1.2701063e-07, -1.1505596, -0.90194077,  0.025499706,  -0.02778281, 0.19655792,
            -0.75007101,   -6.0405615, -0.062160991, -0.024326439, 0.034452336, 0.060048359,
            -3.5964917
        ],
        [
            0.095310307, -1.3032204, -0.87441984,  0.025616812, -0.02780758, 0.19844944,
            -0.74830301, -6.2845955, -0.062012589, -0.02865246, 0.032092044, 0.056613581,
            -3.7057681
        ],
        [
            0.18232168,  -1.4394368, -0.85067618, 0.025702415,  -0.027828006, 0.20003629,
            -0.74688819, -6.5055157, -0.06188419, -0.032691478, 0.030194532,  0.053801717,
            -3.8050113
        ],
        [
            0.26236439,  -1.5622964, -0.82990654,  0.025766388,  -0.02784517, 0.2013865,
            -0.74573886, -6.7071731, -0.061772175, -0.036182428, 0.02863492,  0.05145519,
            -3.8958376
        ],
        [
            0.40546524,  -1.7766356, -0.7951015, 0.025852796,  -0.027872464, 0.20356087,
            -0.74400496, -7.0640081, -0.0615865, -0.043074754, 0.02621949,   0.047754333,
            -4.05705
        ],
        [
            0.53062838, -1.9590895, -0.76685855,  0.025905833,  -0.027893254, 0.20523495,
            -0.7427797, -7.3723499, -0.061439208, -0.049655604, 0.024430981,  0.044958288,
            -4.1968102
        ],
        [
            0.69314731,  -2.1898651, -0.73287511,  0.025952319,  -0.027916653, 0.2071298,
            -0.74152432, -7.7677358, -0.061268042, -0.057970549, 0.022468234,  0.041830672,
            -4.3765769
        ],
        [
            0.9932519,   -2.6006608, -0.67685741,  0.025996103,  -0.027951219, 0.20993425,
            -0.73996133, -8.4842742, -0.061006493, -0.073353045, 0.019644021,  0.037209672,
            -4.7037764
        ],
        [
            1.0986124,   -2.7408691, -0.65895442,  0.026003473,  -0.027961212, 0.2107401,
            -0.73958631, -8.7320382, -0.060929544, -0.078825407, 0.018840073,  0.035865404,
            -4.8173087
        ],
        [
            1.2089605,   -2.8857772, -0.64105694,  0.026008348,  -0.027970632, 0.21149929,
            -0.73926544, -8.9896045, -0.060856292, -0.085234312, 0.018080415,  0.034582329,
            -4.9355356
        ],
        [
            1.3862945,   -3.1149101, -0.61394508,  0.026011951,  -0.027983862, 0.21255952,
            -0.73887727, -9.3996914, -0.060753173, -0.094511097, 0.017008885,  0.0327499,
            -5.124191
        ],
        [
            1.609438,    -3.3974813, -0.58237508,  0.026011882, -0.027997647, 0.21365556,
            -0.73855033, -9.9096322, -0.060645203, -0.10633737, 0.015872998,  0.030776318,
            -5.3594791
        ],
        [
            1.7917596,   -3.6242393, -0.55841141,  0.026009722, -0.028006877, 0.21438776,
            -0.73838016, -10.321803, -0.060572224, -0.11729806, 0.015084767,  0.029386426,
            -5.5502079
        ],
        [
            1.9459103,   -3.8134316, -0.53927995,  0.026007251, -0.028013576, 0.214912,
            -0.73828065, -10.667481, -0.060519954, -0.12449749, 0.014497866,  0.028339943,
            -5.7105432
        ],
        [
            2.0794417,   -3.9756352, -0.52346251,  0.026004786, -0.028018543, 0.21530501,
            -0.73821788, -10.96503,  -0.060480483, -0.13166253, 0.014039265,  0.02751498,
            -5.8488194
        ],
        [
            2.3027864,   -4.2438195, -0.49841323,  0.026000807, -0.028025895, 0.21585771,
            -0.73815709, -11.459171, -0.060425269, -0.14498831, 0.013359124,  0.026043712,
            -6.0867442
        ],
        [
            2.9957484,   -5.0558025, -0.42988106, 0.025990176, -0.028040907, 0.21696658,
            -0.73811889, -12.969313, -0.06031417, -0.1835103,  0.011752027,  0.023203476,
            -6.7899373
        ],
        [
            3.9120585,   -6.0961274, -0.3549036,  0.025982345, -0.028054633, 0.21764827,
            -0.73846608, -14.928233, -0.06025176, -0.23460123, 0.01033696,   0.020565229,
            -7.7135402
        ],
        [
            4.605385,    -6.8661551, -0.30651472,  0.025979545, -0.028054016, 0.21785934,
            -0.73824399, -16.391506, -0.060225593, -0.27330028, 0.0095690216, 0.019086485,
            -8.4096582
        ],
        [
            5.2984392,   -7.6253954, -0.26345896,  0.025978129, -0.028055469, 0.21797026,
            -0.73825218, -17.843043, -0.060214313, -0.31199378, 0.0089606401, 0.017896501,
            -9.1042401
        ],
        [
            6.2147932,   -8.6172642, -0.2128464,   0.025977284, -0.028056337, 0.21803684,
            -0.73825733, -19.749912, -0.060207541, -0.36315699, 0.0083189631, 0.0166284,
            -10.021591
        ],
        [
            6.9078492,   -9.3604565, -0.1783747,   0.025977007, -0.028056625, 0.21805903,
            -0.73825908, -21.185121, -0.060205284, -0.40185229, 0.0079200437, 0.015835434,
            -10.715007
        ],
        [
            7.6010824,   -10.099095, -0.14657042,  0.025976871, -0.028056768, 0.21807013,
            -0.73825996, -22.616029, -0.060204156, -0.44055675, 0.0075756123, 0.015148943,
            -11.40843
        ],
        [
            8.5173098,   -11.06943,  -0.10796354,  0.02597679,  -0.028056859, 0.21807679,
            -0.73826049, -24.501475, -0.060203479, -0.49171046, 0.0071845399, 0.014368186,
            -12.324778
        ],
        [
            9.2105245,   -11.799894, -0.08095212,  0.025976763, -0.028056882, 0.21807901,
            -0.73826067, -25.924427, -0.060203254, -0.53041274, 0.0069268467, 0.013853253,
            -13.018036
        ],
        [
            10.81995,    -13.48613,  -0.024185228, 0.02597674,  -0.028056901, 0.21808079,
            -0.73826081, -29.218735, -0.060177013, -0.61376989, 0.0064233866, 0.012846687,
            -14.627497
        ],
        [
            11.512959,   -14.208805, -0.0018901008, 0.025976736, -0.028056904, 0.21808101,
            -0.73826083, -30.633952, -0.059888051,  -0.63473627, 0.0062385401, 0.012477038,
            -15.320511
        ],
        [
            13.122378,   -15.88086, 0.045828374,  0.025976731, -0.028056901, 0.21808119,
            -0.73826085, -33.91454, -0.058514314, -0.68117163, 0.005864964,  0.01172992,
            -16.929934
        ],
        [
            13.815586,   -16.598754, 0.064867663,  0.025976729, -0.028056919, 0.21808121,
            -0.73826085, -35.325326, -0.057773511, -0.70031902, 0.0057237548, 0.011447506,
            -17.623142
        ],
        [
            16.118208,  -18.975448, 0.12274786, 0.025976724, -0.028051537, 0.21808121,
            -0.7608321, -40.003806, -0.055135,  -0.7608321,  0.0053195561, 0.010639112,
            -19.925765
        ],
        [
            18.420827,   -21.342416, 0.1739028,    0.025976719, -0.02805689,  0.2180807,
            -0.73826139, -44.672781, -0.060203034, -1.0446189,  0.0049910829, 0.0099821659,
            -22.228385
        ]
    ],

    tail_shock_table => [
        [ 7.5420565, -0.50399,    -0.033938669, -0.033938669 ],
        [ 7.6006505, -0.50614174, -0.038423903, -0.029266759 ],
        [ 8.5171302, -0.53881232, -0.048392768, -0.015441115 ],
        [ 9.2104323, -0.56234708, -0.050435696, -0.010919956 ],
        [ 10.819931, -0.6137689,  -0.051433378, -0.0051981939 ],
        [ 11.512949, -0.63473577, -0.05118639,  -0.0037461634 ],
        [ 13.122376, -0.68117152, -0.050012239, -0.0015404454 ],
        [ 13.815585, -0.70033513, -0.049382776, -0.00091329578 ]
    ],

    blast_info => [
        0.14409705, 1.1063848,  0.31148224, -0.19897879, 1.3504024,  0.0034561824, 3.2885701,   -0.086109318,
        0.37025691, 0.44273759, 0.51214479, -0.13528123, 0.43369209, 0.12506613,   0.022202328, -0.023980252,
        1885.2,     -0.50399,   38098.076,  -0.60507904
    ],

};
1;
