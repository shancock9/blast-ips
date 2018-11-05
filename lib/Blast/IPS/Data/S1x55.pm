package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=1.55
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S1x55'} = {

    table_name  => 'S1x55',
    symmetry    => 2,
    gamma       => 1.55,
    data_source => 'S8000_G1x55/moc2_from_moc_r100/',

    shock_table_info => [ 2, 1.55, 8.93e-07, 8000, 2.1e-07, 10, 1.3e-07, 2.63e-07, 5e-07, 1915.4, -0.67705 ],

    shock_table => [
        [ -4.52616832,  12.00017656,   -2.9999776,  -4.52726203, 0.99835855 ],
        [ -4.14802827,  10.86577056,   -2.99993034, -4.14995765, 0.99710319 ],
        [ -3.982641,    10.3696215,    -2.99989719, -3.98511428, 0.99628561 ],
        [ -3.57829877,  9.15667568,    -2.99965655, -3.58283934, 0.99317441 ],
        [ -3.24259731,  8.14977028,    -2.99906139, -3.25012056, 0.98867598 ],
        [ -2.99110115,  7.3956347,     -2.99800131, -3.00208913, 0.98343777 ],
        [ -2.76222396,  6.70965962,    -2.99603977, -2.77774292, 0.97656987 ],
        [ -2.5697387,   6.13323387,    -2.99296736, -2.59049614, 0.96860965 ],
        [ -2.39964847,  5.62451996,    -2.98834309, -2.4265,     0.95933044 ],
        [ -2.24726939,  5.16962665,    -2.98171742, -2.28109695, 0.94869539 ],
        [ -2.10797231,  4.75488115,    -2.97250756, -2.14976208, 0.93655713 ],
        [ -1.97908893,  4.37252771,    -2.96007259, -2.02991285, 0.9228089 ],
        [ -1.85759591,  4.01384855,    -2.94354711, -1.91871788, 0.9072016 ],
        [ -1.74114176,  3.67226022,    -2.92185061, -1.81407677, 0.88942778 ],
        [ -1.62692712,  3.34009095,    -2.89340775, -1.71363019, 0.86894929 ],
        [ -1.51112375,  3.00711755,    -2.85561706, -1.6143683,  0.84477137 ],
        [ -1.38517462,  2.65066455,    -2.80240923, -1.50984296, 0.81428871 ],
        [ -1.26267382,  2.31122092,    -2.73722125, -1.41212699, 0.78036857 ],
        [ -1.15078463,  2.0088441,     -2.66587605, -1.32671211, 0.74587778 ],
        [ -1.0459759,   1.73335959,    -2.58961736, -1.25034987, 0.7109077 ],
        [ -0.94665395,  1.48003945,    -2.510395,   -1.18146731, 0.6758918 ],
        [ -0.85067331,  1.24295743,    -2.42921696, -1.11826806, 0.64088287 ],
        [ -0.75493457,  1.01437293,    -2.3456989,  -1.0586065,  0.60541873 ],
        [ -0.65955465,  0.79463865,    -2.2619188,  -1.00254649, 0.57014993 ],
        [ -0.56347998,  0.58133782,    -2.1787286,  -0.94945295, 0.53525499 ],
        [ -0.46581163,  0.37256863,    -2.09690207, -0.89886275, 0.50093472 ],
        [ -0.36570586,  0.16668969,    -2.01708763, -0.85041075, 0.46738464 ],
        [ -0.26227991,  -0.03788363,   -1.93977066, -0.80377663, 0.43476972 ],
        [ -0.15483873,  -0.24224736,   -1.86546893, -0.75877772, 0.40329685 ],
        [ -0.04241569,  -0.44791164,   -1.7944247,  -0.71516264, 0.37307531 ],
        [ 0.07594203,   -0.65622561,   -1.72685469, -0.67274364, 0.34421377 ],
        [ 0.20113024,   -0.86833063,   -1.66297289, -0.63139905, 0.3168242 ],
        [ 0.33414032,   -1.08544357,   -1.60290698, -0.59101233, 0.29098156 ],
        [ 0.47623677,   -1.30912348,   -1.54665312, -0.55142823, 0.26670126 ],
        [ 0.62871846,   -1.54086488,   -1.4942114,  -0.51253321, 0.24399579 ],
        [ 0.78810209,   -1.77515348,   -1.44689712, -0.47532179, 0.22344523 ],
        [ 0.95404924,   -2.01165607,   -1.40450829, -0.43981227, 0.20497156 ],
        [ 1.12806174,   -2.25265313,   -1.36635213, -0.40563324, 0.18827922 ],
        [ 1.31119364,   -2.49964499,   -1.33196304, -0.37257184, 0.17316844 ],
        [ 1.504988,     -2.75467969,   -1.30087245, -0.34037736, 0.15943579 ],
        [ 1.7103275,    -3.01884011,   -1.27279382, -0.308953,   0.14695685 ],
        [ 1.92870505,   -3.2939414,    -1.24739799, -0.27813441, 0.13558772 ],
        [ 2.16151931,   -3.58160543,   -1.22442411, -0.2478066,  0.1252144 ],
        [ 2.41035053,   -3.8836217,    -1.20363533, -0.21785973, 0.11573362 ],
        [ 2.67647009,   -4.20136261,   -1.18485142, -0.18824401, 0.10706797 ],
        [ 2.96176342,   -4.53690149,   -1.1678716,  -0.15886013, 0.09913105 ],
        [ 3.26830602,   -4.89248261,   -1.15251813, -0.12961801, 0.09184709 ],
        [ 3.5979633,    -5.27006582,   -1.13864888, -0.10047159, 0.08515742 ],
        [ 3.9531938,    -5.67225701,   -1.126116,   -0.0713433,  0.07900095 ],
        [ 4.3366475,    -6.10183406,   -1.11478852, -0.04216693, 0.07332486 ],
        [ 4.7511665,    -6.56174946,   -1.10454944, -0.01288725, 0.06808331 ],
        [ 5.1997855,    -7.05513269,   -1.09529337, 0.01654053,  0.06323636 ],
        [ 5.68653268,   -7.58616228,   -1.08691212, 0.04619826,  0.05874221 ],
        [ 6.21522977,   -8.15874171,   -1.07931957, 0.07612337,  0.05457017 ],
        [ 6.7904933,    -8.77759377,   -1.07243211, 0.10637042,  0.05069024 ],
        [ 7.41753526,   -9.44803226,   -1.0661744,  0.13699247,  0.0470759 ],
        [ 8.10176399,   -10.17553315,  -1.06048267, 0.16802004,  0.04370606 ],
        [ 8.84958505,   -10.96658662,  -1.05529684, 0.19949717,  0.04056031 ],
        [ 9.66760191,   -11.82784413,  -1.05056697, 0.23144298,  0.03762278 ],
        [ 10.56321664,  -12.76675486,  -1.04624799, 0.26387718,  0.0348789 ],
        [ 11.54423036,  -13.79114426,  -1.04230141, 0.29680355,  0.03231661 ],
        [ 12.61944367,  -14.90984172,  -1.0386921,  0.33023022,  0.0299243 ],
        [ 13.79745683,  -16.13143222,  -1.03539174, 0.36413202,  0.02769334 ],
        [ 15.08846996,  -17.46612891,  -1.03237238, 0.39850588,  0.02561388 ],
        [ 16.5032831,   -18.92472987,  -1.02960974, 0.43333749,  0.02367726 ],
        [ 18.05289625,  -20.51821016,  -1.02708299, 0.46859438,  0.02187602 ],
        [ 19.65134946,  -22.1581562,   -1.02489521, 0.50226391,  0.02029204 ],
        [ 21.4355891,   -23.98492718,  -1.02283884, 0.53708639,  0.01878155 ],
        [ 23.31815668,  -25.90872083,  -1.02101011, 0.57112901,  0.01741978 ],
        [ 25.44722287,  -28.08060999,  -1.01926747, 0.60677965,  0.01610505 ],
        [ 27.68318259,  -30.35787849,  -1.01772514, 0.64143851,  0.0149269 ],
        [ 30.21240917,  -32.93002257,  -1.01625478, 0.67771585,  0.01379036 ],
        [ 32.85257448,  -35.61134707,  -1.01496043, 0.71275435,  0.01277861 ],
        [ 35.83762492,  -38.63915921,  -1.01372573, 0.74940374,  0.01180309 ],
        [ 38.92925603,  -41.77150966,  -1.01264556, 0.78452752,  0.01094101 ],
        [ 42.42014835,  -45.30469485,  -1.01161419, 0.82123188,  0.01010996 ],
        [ 46.38088493,  -49.30942759,  -1.01063103, 0.85964658,  0.0093102 ],
        [ 50.46192921,  -53.43206952,  -1.00977858, 0.8961763,   0.0086105 ],
        [ 55.08274662,  -58.09613773,  -1.00896502, 0.93436497,  0.00793706 ],
        [ 59.78358765,  -62.83743363,  -1.00826571, 0.97026887,  0.00735356 ],
        [ 65.08811086,  -68.18397966,  -1.00759726, 1.00774539,  0.00679159 ],
        [ 71.10281451,  -74.24240388,  -1.00695932, 1.0469243,   0.00625128 ],
        [ 77.1506031,   -80.3305987,   -1.00641765, 1.08329737,  0.00578927 ],
        [ 83.97505033,  -87.19702555,  -1.00589962, 1.12124746,  0.00534455 ],
        [ 91.71341549,  -94.97907458,  -1.00540504, 1.16090542,  0.00491718 ],
        [ 99.36596863,  -102.67136504, -1.00499133, 1.19712298,  0.00455753 ],
        [ 107.98510444, -111.33176782, -1.00459526, 1.23487291,  0.00421126 ],
        [ 117.73897518, -121.12856082, -1.00421667, 1.27428083,  0.00387844 ],
        [ 128.83380023, -132.26810533, -1.00385544, 1.3154886,   0.00355909 ],
        [ 139.60269772, -143.07688569, -1.0035595,  1.35236357,  0.00329613 ],
        [ 151.74014786, -155.25577166, -1.00327611, 1.39079341,  0.00304312 ],
        [ 165.48616617, -169.04490766, -1.00300517, 1.43090704,  0.00280009 ],
        [ 178.38067785, -181.97674012, -1.00278885, 1.46572585,  0.0026052 ],
        [ 192.80295454, -196.43770119, -1.00258105, 1.50190958,  0.00241728 ],
        [ 209.00233555, -212.67723603, -1.00238174, 1.53956403,  0.00223633 ]
    ],

    energy_table => [
        [ -4.52616832,  0.01350413, 0.01435713, 0.01350413, 0.01435713,    5.5871225,   -1.5998671 ],
        [ -4.14802827,  0.02018181, 0.02142759, 0.02018181, 0.02142759,    4.9776742,   -1.6257209 ],
        [ -3.982641,    0.02405389, 0.02551384, 0.02405389, 0.02551384,    4.7077864,   -1.639654 ],
        [ -3.57829877,  0.03690707, 0.03898503, 0.03690707, 0.03898503,    4.0370936,   -1.6825975 ],
        [ -3.24259731,  0.0525583,  0.05513996, 0.0525583,  0.05513996,    3.4655897,   -1.7313824 ],
        [ -2.99110115,  0.06835484, 0.07108918, 0.06835484, 0.07108918,    3.0246962,   -1.7768847 ],
        [ -2.76222396,  0.08660761, 0.08897088, 0.08660761, 0.08897088,    2.61305,     -1.8293269 ],
        [ -2.5697387,   0.10539595, 0.10664812, 0.10539595, 0.10664812,    2.2559488,   -1.885849 ],
        [ -2.39964847,  0.12500286, 0.12417579, 0.12500286, 0.12417579,    1.9305798,   -1.9452493 ],
        [ -2.24726939,  0.14520264, 0.14110763, 0.14520264, 0.14110763,    1.6297512,   -2.0067524 ],
        [ -2.10797231,  0.16597388, 0.15716888, 0.16597388, 0.15716888,    1.3460728,   -2.0665985 ],
        [ -1.97908893,  0.18718912, 0.17198793, 0.18718912, 0.17198793,    1.0761333,   -2.1283229 ],
        [ -1.85759591,  0.20890337, 0.18530631, 0.20890337, 0.18530631,    0.81367728,  -2.195042 ],
        [ -1.74114176,  0.23116837, 0.19681009, 0.23116837, 0.19681009,    0.55417254,  -2.2516313 ],
        [ -1.62692712,  0.25420368, 0.20618241, 0.25420368, 0.20618241,    0.29439932,  -2.2998341 ],
        [ -1.51112375,  0.27850503, 0.21300714, 0.27850503, 0.21300714,    0.025087999, -2.3293382 ],
        [ -1.38517462,  0.30560279, 0.21655478, 0.30560279, 0.21655478,    -0.26880337, -2.3146805 ],
        [ -1.26267382,  0.33211457, 0.21551595, 0.33211457, 0.21551595,    -0.55012163, -2.2451522 ],
        [ -1.15078463,  0.35598331, 0.21048778, 0.35598331, 0.21048778,    -0.79608634, -2.1327312 ],
        [ -1.0459759,   0.37764727, 0.20239514, 0.37764727, 0.20239514,    -1.013179,   -2.0038185 ],
        [ -0.94665395,  0.397257,   0.19209144, 0.397257,   0.19209144,    -1.2058485,  -1.8906962 ],
        [ -0.85067331,  0.41513681, 0.18021271, 0.41513681, 0.18021271,    -1.3827598,  -1.8129217 ],
        [ -0.75493457,  0.43176882, 0.1670677,  0.43176882, 0.1670677,     -1.5534354,  -1.7700322 ],
        [ -0.65955465,  0.44704938, 0.15327979, 0.44704938, 0.15327979,    -1.7210551,  -1.7578921 ],
        [ -0.56347998,  0.46110086, 0.13925066, 0.46110086, 0.13925066,    -1.8899926,  -1.7670069 ],
        [ -0.46581163,  0.47401686, 0.125328,   0.47401686, 0.125328,      -2.063428,   -1.7885993 ],
        [ -0.36570586,  0.48587825, 0.1117998,  0.48587825, 0.1117998,     -2.2437947,  -1.8162185 ],
        [ -0.26227991,  0.49676339, 0.09888841, 0.49676339, 0.09888841,    -2.4331837,  -1.845538 ],
        [ -0.15483873,  0.50672548, 0.08678322, 0.50672548, 0.08678322,    -2.6330757,  -1.873834 ],
        [ -0.04241569,  0.51583874, 0.07559082, 0.51583874, 0.07559082,    -2.8453089,  -1.8996478 ],
        [ 0.07594203,   0.52416587, 0.06537999, 0.52416587, 0.06537999,    -3.0716221,  -1.922243 ],
        [ 0.20113024,   0.53175902, 0.0561877,  0.53175902, 0.0561877,     -3.3136057,  -1.9413565 ],
        [ 0.33414032,   0.53867194, 0.04801089, 0.53867194, 0.04801089,    -3.5730125,  -1.9570558 ],
        [ 0.47623677,   0.54496528, 0.04080779, 0.54496528, 0.04080779,    -3.8521319,  -1.9695944 ],
        [ 0.62871846,   0.55069152, 0.03452317, 0.55069152, 0.03452317,    -4.153326,   -1.9793202 ],
        [ 0.78810209,   0.555757,   0.02923346, 0.555757,   0.02923346,    -4.4694703,  -1.9865046 ],
        [ 0.95404924,   0.56022821, 0.02481706, 0.56022821, 0.02481706,    -4.7996368,  -1.9916216 ],
        [ 1.12806174,   0.56421236, 0.02111298, 0.56421236, 0.02111298,    -5.1465749,  -1.9951348 ],
        [ 1.31119364,   0.5677832,  0.01800173, 0.5677832,  0.01800173,    -5.5122139,  -1.9975024 ],
        [ 1.504988,     0.57100801, 0.01537866, 0.57100801, 0.01537866,    -5.899506,   -1.9990461 ],
        [ 1.7103275,    0.57393016, 0.01316738, 0.57393016, 0.01316738,    -6.3101146,  -1.9999817 ],
        [ 1.92870505,   0.57659369, 0.01129798, 0.57659369, 0.01129798,    -6.7469411,  -2.0004921 ],
        [ 2.16151931,   0.57903265, 0.00971474, 0.57903265, 0.00971474,    -7.2127252,  -2.0007134 ],
        [ 2.41035053,   0.58127637, 0.00837097, 0.58127637, 0.00837097,    -7.7105772,  -2.0007548 ],
        [ 2.67647009,   0.58334636, 0.00722966, 0.58334636, 0.00722966,    -8.2430152,  -2.0007259 ],
        [ 2.96176342,   0.58526496, 0.00625759, 0.58526496, 0.00625759,    -8.8138023,  -2.0006317 ],
        [ 3.26830602,   0.58705107, 0.00542735, 0.58705107, 0.00542735,    -9.4270578,  -2.0005 ],
        [ 3.5979633,    0.58871871, 0.0047171,  0.58871871, 0.0047171,     -10.086518,  -2.0003902 ],
        [ 3.9531938,    0.59028203, 0.00410768, 0.59028203, 0.00410768,    -10.797098,  -2.0002873 ],
        [ 4.3366475,    0.59175281, 0.00358331, 0.59175281, 0.00358331,    -11.564095,  -2.0001985 ],
        [ 4.7511665,    0.59314089, 0.00313094, 0.59314089, 0.00313094,    -12.393199,  -2.000133 ],
        [ 5.1997855,    0.59445447, 0.00273974, 0.59445447, 0.00273974,    -13.290484,  -2.0000861 ],
        [ 5.68653268,   0.59570233, 0.00240018, 0.59570233, 0.00240018,    -14.264011,  -2.0000528 ],
        [ 6.21522977,   0.59689035, 0.00210457, 0.59689035, 0.00210457,    -15.321425,  -2.0000308 ],
        [ 6.7904933,    0.59802435, 0.0018459,  0.59802435, 0.0018459,     -16.471965,  -2.0000176 ],
        [ 7.41753526,   0.59910945, 0.00162208, 0.59910945, 0.00162208,    -17.726057,  -2.0000093 ],
        [ 8.10176399,   0.60014929, 0.0014246,  0.60016162, 0.0014799203,  -19.094519,  -2.0000047 ],
        [ 8.84958505,   0.60114748, 0.00125132, 0.60125654, 0.0014554163,  -20.590163,  -2.0000022 ],
        [ 9.66760191,   0.60210654, 0.00109908, 0.6024198,  0.0013976597,  -22.226198,  -2.000001 ],
        [ 10.56321664,  0.60302875, 0.00096521, 0.60364651, 0.0013364684,  -24.017428,  -2.0000004 ],
        [ 11.54423036,  0.60391576, 0.00084744, 0.60492352, 0.0012635798,  -25.979455,  -2.0000002 ],
        [ 12.61944367,  0.60476917, 0.00074379, 0.60624423, 0.0011900697,  -28.129882,  -2.0000001 ],
        [ 13.79745683,  0.60558972, 0.00065264, 0.60759339, 0.0011014157,  -30.485908,  -2 ],
        [ 15.08846996,  0.60637864, 0.00057248, 0.60919519, 0.0011643081,  -33.067935,  -2 ],
        [ 16.5032831,   0.60713691, 0.00050201, 0.61073609, 0.0010193213,  -35.897561,  -2 ],
        [ 18.05289625,  0.6078651,  0.00044011, 0.61221358, 0.00089232252, -38.996787,  -2 ],
        [ 19.65134946,  0.60852599, 0.0003886,  0.6135527,  0.00078690251, -42.193694,  -2 ],
        [ 21.4355891,   0.60917633, 0.00034205, 0.61486887, 0.0006918438,  -45.762173,  -2 ],
        [ 23.31815668,  0.60978151, 0.00030226, 0.61609231, 0.00061074241, -49.527308,  -2 ],
        [ 25.44722287,  0.61038485, 0.00026582, 0.61731086, 0.00053660741, -53.785441,  -2 ],
        [ 27.68318259,  0.61094337, 0.00023484, 0.61843788, 0.0004736761,  -58.25736,   -2 ],
        [ 30.21240917,  0.61150015, 0.00020647, 0.61956052, 0.00041613482, -63.315813,  -2 ],
        [ 32.85257448,  0.6120125,  0.00018249, 0.62059284, 0.00036755458, -68.596144,  -2 ],
        [ 35.83762492,  0.61252325, 0.00016052, 0.62162128, 0.00032310378, -74.566245,  -2 ],
        [ 38.92925603,  0.61298999, 0.00014206, 0.62256056, 0.00028579884, -80.749507,  -2 ],
        [ 42.42014835,  0.61345529, 0.00012513, 0.62349648, 0.00025161919, -87.731292,  -2 ],
        [ 46.38088493,  0.61391912, 0.00010966, 0.62442899, 0.00022040574, -95.652765,  -2 ],
        [ 50.46192921,  0.61433947, 9.679e-05,  0.62527371, 0.00019447228, -103.81485,  -2 ],
        [ 55.08274662,  0.61475852, 8.502e-05,  0.62611553, 0.00017074563, -113.05649,  -2 ],
        [ 59.78358765,  0.61513454, 7.53e-05,   0.62687065, 0.00015118516, -122.45817,  -2 ],
        [ 65.08811086,  0.61550947, 6.638e-05,  0.62762334, 0.00013324034, -133.06722,  -2 ],
        [ 71.10281451,  0.61588328, 5.822e-05,  0.62837359, 0.00011683096, -145.09662,  -2 ],
        [ 77.1506031,   0.6162146,  5.158e-05,  0.6290384,  0.00010347107, -157.1922,   -2 ],
        [ 83.97505033,  0.616545,   4.547e-05,  0.62970122, 9.1208529e-05, -170.8411,   -2 ],
        [ 91.71341549,  0.61687446, 3.989e-05,  0.63036204, 7.9990546e-05, -186.31783,  -2 ],
        [ 99.36596863,  0.61716196, 3.541e-05,  0.63093859, 7.0990189e-05, -201.62293,  -2 ],
        [ 107.98510444, 0.61744873, 3.128e-05,  0.63151357, 6.2714811e-05, -218.8612,   -2 ],
        [ 117.73897518, 0.61773474, 2.75e-05,   0.63208697, 5.5130095e-05, -238.36895,  -2 ],
        [ 128.83380023, 0.61801999, 2.405e-05,  0.63265877, 4.8202733e-05, -260.5586,   -2 ],
        [ 139.60269772, 0.61826388, 2.134e-05,  0.63314759, 4.2762767e-05, -282.09639,  -2 ],
        [ 151.74014786, 0.61850719, 1.884e-05,  0.63363522, 3.7760691e-05, -306.37129,  -2 ],
        [ 165.48616617, 0.61874992, 1.656e-05,  0.63412163, 3.3176053e-05, -333.86333,  -2 ],
        [ 178.38067785, 0.61895174, 1.48e-05,   0.63452605, 2.9660101e-05, -359.65235,  -2 ],
        [ 192.80295454, 0.61915316, 1.318e-05,  0.63492961, 2.6408484e-05, -388.4969,   -2 ],
        [ 209.00233555, 0.61935415, 1.169e-05,  0.63533231, 2.3410076e-05, -420.89567,  -2 ]
    ],

    impulse_table => [
        [
            -5.131058,  13.814839,  -5.1314993,    0.024303419, -0.00083778412, -0.15561218,
            -1.6060242, -0.2471438, -0.0014175098, -0.37087645, 0.25846554,     0.98900057,
            -2.3258322
        ],
        [
            -4.9618449, 13.307201,  -4.9624138,    0.026123122, -0.00099224846, -0.15452248,
            -1.6049345, -0.3317519, -0.0016788596, -0.36978675, 0.25846543,     0.98683079,
            -2.272791
        ],
        [
            -4.60517,   12.23718,   -4.6061414,    0.030285423, -0.0014174978, -0.15152248,
            -1.6019345, -0.5100943, -0.0023983709, -0.36678675, 0.25846477,    0.98075484,
            -2.163964
        ],
        [
            -4.1997049, 11.020797,   -4.2014901,    0.035512036, -0.0021262467, -0.14652248,
            -1.5969345, -0.71284357, -0.0035975564, -0.36178675, 0.25846242,    0.97038828,
            -2.0464577
        ],
        [
            -3.9120228, 10.157775,   -3.9147724,    0.039468374, -0.0028349956, -0.14152248,
            -1.5919345, -0.85671669, -0.0047967418, -0.35678675, 0.25845802,    0.95981688,
            -1.9683296
        ],
        [
            -3.5065577, 8.9414798, -3.5116149,    0.045182984, -0.0042524934, -0.13152247,
            -1.5819345, -1.059581, -0.0071951127, -0.34678675, 0.25843997,    0.93829155,
            -1.8681575
        ],
        [
            -2.9957321, 7.4095183, -3.0066436,   0.051893308, -0.0070874887, -0.11152249,
            -1.5619345, -1.315674, -0.011991854, -0.32678674, 0.25834697,    0.89454958,
            -1.7654926
        ],
        [
            -2.3025849, 5.3346395,  -2.3336907,   0.056986805, -0.014174893, -0.061526013,
            -1.5119363, -1.6683336, -0.023983604, -0.27678673, 0.2575222,    0.7866464,
            -1.6965867
        ],
        [
            -1.8971198, 4.1303088,  -1.9546814,   0.055361939, -0.02125984, -0.011595352,
            -1.4619733, -1.8877003, -0.035972321, -0.22678672, 0.25533382,  0.68693284,
            -1.721093
        ],
        [
            -1.6094377, 3.289531,   -1.6984626,   0.051577824, -0.028320185, 0.03785885,
            -1.4122661, -2.0643764, -0.047930589, -0.17722487, 0.25126828,   0.59913458,
            -1.7848569
        ],
        [
            -1.3862942, 2.6538023,  -1.5107548,   0.047892914, -0.035245312, 0.085308303,
            -1.3636655, -2.2311465, -0.059719006, -0.12883488, 0.2450635,    0.52430651,
            -1.8694261
        ],
        [
            -1.2039726, 2.1515862,  -1.3668336,   0.04554609,   -0.041692217, 0.12728461,
            -1.3183906, -2.4061588, -0.070875007, -0.084695455, 0.2367634,    0.46188128,
            -1.9650417
        ],
        [
            -1.0498219, 1.743325,   -1.2530866,   0.044818462,  -0.047048496, 0.15991146,
            -1.28003,   -2.5998624, -0.080418219, -0.049001427, 0.22668362,   0.41042938,
            -2.0655006
        ],
        [
            -0.91629054, 1.4041984,  -1.1611115,   0.045261477,  -0.050859345, 0.18264374,
            -1.2509535,  -2.8134061, -0.087321035, -0.025372091, 0.21531096,   0.36821735,
            -2.166714
        ],
        [
            -0.7985075, 1.1174154,  -1.0853388,   0.046270787,  -0.053260527, 0.19801416,
            -1.2301011, -3.0386035, -0.091536879, -0.012681463, 0.20318932,   0.33355314,
            -2.2660874
        ],
        [
            -0.69314699, 0.87111659, -1.021907,   0.047435186,   -0.054716032, 0.20874065,
            -1.2149147,  -3.2647254, -0.09387844, -0.0070561101, 0.19083161,   0.3049525,
            -2.3620962
        ],
        [
            -0.59783681, 0.65669838, -0.96805488,  0.048552345,   -0.055620635, 0.21662497,
            -1.2033652,  -3.4843046, -0.095148717, -0.0055282645, 0.17866919,   0.28118815,
            -2.4539459
        ],
        [
            -0.51082543, 0.46779472, -0.92176276,  0.04954575,    -0.0562116, 0.22271059,
            -1.1942086,  -3.6936729, -0.095840558, -0.0058280548, 0.16703262, 0.26127816,
            -2.5413138
        ],
        [
            -0.43078272, 0.29961633, -0.88152545,  0.050399165,   -0.056618275, 0.22759561,
            -1.1867085,  -3.891558,  -0.096216409, -0.0074723342, 0.15615237,   0.24445011,
            -2.6241645
        ],
        [
            -0.35667475, 0.14850466, -0.84620302,  0.051121148,   -0.0569113, 0.23163564,
            -1.1804144,  -4.0779176, -0.096413858, -0.0095493761, 0.1461701,  0.2301008,
            -2.7026287
        ],
        [
            -0.28768188, 0.011624347, -0.81491951,  0.051728463, -0.057130715, 0.23505315,
            -1.1750362,  -4.2532761,  -0.096506955, -0.01201508, 0.13715408,   0.21775962,
            -2.7769253
        ],
        [
            -0.22314336, -0.11325307, -0.78699237,  0.052239031,  -0.05730031, 0.23799474,
            -1.1703768,  -4.4183866,  -0.096537313, -0.014714677, 0.12911534,  0.20705836,
            -2.8473138
        ],
        [
            -0.16251874, -0.22790097, -0.76188335,  0.052669202,  -0.057434873, 0.2405615,
            -1.166296,   -4.5740645,  -0.096529189, -0.017455372, 0.12202273,   0.19770741,
            -2.9140665
        ],
        [
            -0.10536032, -0.33374803, -0.73916295,  0.053032957,  -0.057543979, 0.24282582,
            -1.1626919,  -4.7211097,  -0.096497241, -0.020225543, 0.11580408,   0.18947751,
            -2.9774518
        ],
        [
            -0.051293101, -0.43195803, -0.71848472,  0.053341923,  -0.057634088, 0.24484161,
            -1.1594842,   -4.8602713,  -0.096450558, -0.023074181, 0.11033938,   0.18218579,
            -3.0377253
        ],
        [
            -0.023973565, -0.4809029, -0.70832572,  0.053485796,  -0.057675447, 0.24581681,
            -1.1579355,   -4.9305705, -0.096422185, -0.024400994, 0.10772633,   0.17867453,
            -3.0682749
        ],
        [
            0.095310373, -0.68957117, -0.66611992,  0.054026578,  -0.057829191, 0.24976451,
            -1.1516962,  -5.2369822,  -0.096271226, -0.030939388, 0.09739439,   0.16460352,
            -3.2021744
        ],
        [
            0.18232175, -0.83696713, -0.63739475,  0.054341974,  -0.057919322, 0.25235659,
            -1.1476473, -5.459581,   -0.096143085, -0.036087735, 0.09085383,   0.15551591,
            -3.3001421
        ],
        [
            0.26236446, -0.96927758, -0.61237844,  0.054582451,  -0.057989678, 0.25455028,
            -1.1442679, -5.6634208,  -0.096018325, -0.041025039, 0.085482551,  0.14793218,
            -3.3903096
        ],
        [
            0.4054653,  -1.1987183, -0.57071196,  0.054915121,  -0.058092498, 0.25806597,
            -1.1389811, -6.0251996, -0.095789738, -0.050036327, 0.077197722,  0.13599382,
            -3.5512832
        ],
        [
            0.53062844, -1.3927117, -0.53715408,  0.055125331,  -0.058164236, 0.26076269,
            -1.1350753, -6.3385407, -0.095592544, -0.058148179, 0.071110003,  0.12701229,
            -3.6915493
        ],
        [
            0.69314737, -1.6364856, -0.4970947,   0.05531563,   -0.058239154, 0.26380801,
            -1.13088,   -6.7408533, -0.095349193, -0.069179708, 0.064497104,  0.11702828,
            -3.8726142
        ],
        [
            0.99325197, -2.0665371, -0.43185523,  0.055507214,  -0.058341071, 0.26830924,
            -1.1252671, -7.4702585, -0.094953707, -0.090486595, 0.055152191,  0.10245277,
            -4.2031693
        ],
        [
            1.0986125,  -2.2123265, -0.41121676,  0.055543104,  -0.058369023, 0.26960204,
            -1.1238261, -7.7223516, -0.094832985, -0.098448008, 0.052537712,  0.098263262,
            -4.3179911
        ],
        [
            1.2089605, -2.3625411, -0.3906861,   0.055569347, -0.058394997, 0.27082036,
            -1.122553, -7.9842878, -0.094716578, -0.10656697, 0.050088881,  0.094289678,
            -4.4375729
        ],
        [
            1.3862946,  -2.5991976, -0.35977791,  0.055594096, -0.058430971, 0.27252219,
            -1.1209326, -8.4009997, -0.094550095, -0.11909134, 0.046673572,  0.088661774,
            -4.6283545
        ],
        [
            1.6094381,  -2.8897704, -0.32407266,  0.055605549, -0.058468078, 0.2742823,
            -1.1194806, -8.9185357, -0.094373187, -0.13538803, 0.043107579,  0.082668784,
            -4.8661496
        ],
        [
            1.7917597,  -3.1220773, -0.29716867,  0.055605431, -0.058492945, 0.27545897,
            -1.1186561, -9.3362992, -0.094252802, -0.14993212, 0.040669005,  0.078495151,
            -5.0587453
        ],
        [
            1.9459103,  -3.3153875, -0.27580867,  0.055601754, -0.058510939, 0.27630142,
            -1.1181455, -9.6862923, -0.094165556, -0.16185498, 0.038873871,  0.075380325,
            -5.2205166
        ],
        [
            2.0794417,  -3.4807963, -0.25822501,  0.055596915, -0.058524426, 0.27693378,
            -1.1178088, -9.9872898, -0.094099154, -0.17230094, 0.037484033,  0.072942487,
            -5.3599387
        ],
        [
            2.302787,   -3.7537,    -0.23051656,  0.055587624, -0.058544612, 0.27782275,
            -1.1174795, -10.486642, -0.094006439, -0.18863437, 0.035444573,  0.068559879,
            -5.6012587
        ],
        [
            2.9959114,  -4.5767504, -0.15548987, 0.055559649, -0.058581495, 0.27960132,
            -1.1167158, -12.009611, -0.09381557, -0.24174961, 0.030735458,  0.060445551,
            -6.3075165
        ],
        [
            3.912133,   -5.6259904, -0.074600692, 0.055537402, -0.058611811, 0.28068373,
            -1.1170503, -13.978814, -0.093704674, -0.31221103, 0.026729013,  0.053090069,
            -7.2332851
        ],
        [
            4.6054701,  -6.4005764, -0.022932864, 0.05552923,  -0.058618886, 0.2810404,
            -1.1169626, -15.447251, -0.093665421, -0.36559411, 0.02461019,   0.049046258,
            -7.930321
        ],
        [
            5.2984215,  -7.1630777, 0.022729945,  0.055525063, -0.058623012, 0.28121925,
            -1.1169694, -16.90226,  -0.093646183, -0.41896394, 0.022958373,  0.045833124,
            -8.6253095
        ],
        [
            6.2146297,  -8.158094,  0.076090623,  0.055522569, -0.05862548, 0.28132664,
            -1.1169746, -18.812308, -0.093634628, -0.48953281, 0.021239206, 0.042446255,
            -9.5428558
        ],
        [
            6.9079028,  -8.9034333, 0.11227948,   0.055521751, -0.058626301, 0.28136246,
            -1.1169764, -20.249922, -0.093630775, -0.54292982, 0.020180534,  0.040345356,
            -10.236613
        ],
        [
            7.6009443,  -9.6434279, 0.14553842,  0.055521349, -0.058626711, 0.28138036,
            -1.1169774, -21.682002, -0.09362885, -0.59630774, 0.01927276,   0.038537809,
            -10.929911
        ],
        [
            8.5173766, -10.61565,  0.18580347,   0.055521113, -0.058626963, 0.28139111,
            -1.116978, -23.569526, -0.093627695, -0.66688967, 0.0182473,    0.036491578,
            -11.846506
        ],
        [
            9.2103931,  -11.346946, 0.21388483,  0.055521035, -0.058624021, 0.28139469,
            -1.1169782, -24.993096, -0.09362731, -0.72026376, 0.017574945,  0.035148398,
            -12.53958
        ],
        [
            10.81982,   -13.035084, 0.27273557,   0.055520973, -0.058627109, 0.28139756,
            -1.1169784, -28.289275, -0.093544019, -0.83205214, 0.016267124,  0.032533957,
            -14.149056
        ],
        [
            11.51303,   -13.758622, 0.29579409,   0.055520964, -0.058627108, 0.28139791,
            -1.1169784, -29.705545, -0.093018993, -0.86172934, 0.015788702,  0.03157726,
            -14.842272
        ],
        [
            13.122449,  -15.431933, 0.34502889,  0.055520955, -0.058627117, 0.2813982,
            -1.1169784, -32.987369, -0.09077872, -0.92735144, 0.014824917,  0.029649805,
            -16.451697
        ],
        [
            13.815657,  -16.150276, 0.36463576,   0.055520953, -0.058627106, 0.28139824,
            -1.1169784, -34.398601, -0.089601366, -0.95437741, 0.014461542,  0.028923071,
            -17.144905
        ],
        [
            16.11828,   -18.528193, 0.42412699,   0.055520948, -0.058627128, 0.28139824,
            -1.1169784, -39.078291, -0.085459786, -1.0396784,  0.01342405,   0.0268481,
            -19.447529
        ],
        [
            18.420699,  -20.895875, 0.47656907,   0.055520943, -0.058627091, 0.28139824,
            -1.1192342, -43.747776, -0.081428124, -1.1192342,  0.012583695,  0.02516739,
            -21.749948
        ]
    ],

    tail_shock_table => [
        [ 7.5580352, -0.67705,    -0.041096777, -0.041096777 ],
        [ 7.6003658, -0.67926405, -0.045426413, -0.036597126 ],
        [ 8.5171358, -0.72575237, -0.05714274,  -0.020047124 ],
        [ 9.2102692, -0.75916539, -0.059383558, -0.014712894 ],
        [ 10.819794, -0.83205026, -0.060351039, -0.0078742585 ],
        [ 11.513017, -0.86172839, -0.060012291, -0.0061084398 ],
        [ 13.122447, -0.92735125, -0.058566926, -0.003376984 ],
        [ 13.815656, -0.95437731, -0.057807337, -0.0025822753 ],
        [ 16.118279, -1.0397141,  -0.055138613, -0.00087245739 ]
    ],

    blast_info => [
        0.16152246, 1.6119339,  0.37678957, -0.2398371, 1.9503888,  0.0019779009, 3.7265555,   -0.14174975,
        0.48377863, 0.38705833, 0.64146755, -0.1660755, 0.33007157, 0.25846569,   0.035819966, -0.037823937,
        1915.4,     -0.67705,   34081.991,  -0.81469587
    ],

};
1;