package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='P', gamma=6
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'P6'} = {

    table_name  => 'P6',
    symmetry    => 0,
    gamma       => 6,
    data_source => 'P4000_G6/moc3p_from_AD_r100/',

    shock_table_info => [ 0, 6, 1.52e-06, 4000, 1.23e-06, 100, 8.3e-07, 1.9e-07, 5e-07 ],

    shock_table => [
        [ -9.96706457, 12.00046263, -0.99998947, -9.96923005, 0.9989161 ],
        [ -9.00486523, 11.03827839, -0.99997245, -9.00837096, 0.99824411 ],
        [ -8.37943566, 10.41286963, -0.99995335, -8.38423147, 0.99759647 ],
        [ -7.0754681,  9.109024,    -0.99983265, -7.08469259, 0.99536741 ],
        [ -6.14888263, 8.18269372,  -0.99957792, -6.16358079, 0.99260065 ],
        [ -5.40661457, 7.44088999,  -0.99911454, -5.42798234, 0.9892134 ],
        [ -4.77618655, 6.81123901,  -0.99834022, -4.80557475, 0.98511957 ],
        [ -4.22635836, 6.26262523,  -0.99713371, -4.26519745, 0.9802712 ],
        [ -3.73901064, 5.77707291,  -0.99535757, -3.7887787,  0.97463843 ],
        [ -3.29977288, 5.34038626,  -0.99284936, -3.36204714, 0.96816831 ],
        [ -2.89677381, 4.94091828,  -0.98940744, -2.97331475, 0.96076964 ],
        [ -2.52189095, 4.57082048,  -0.9847976,  -2.61466649, 0.95234704 ],
        [ -2.16825241, 4.22357326,  -0.97873602, -2.27952196, 0.9427746 ],
        [ -1.82773381, 3.89157775,  -0.97081864, -1.96030362, 0.93181558 ],
        [ -1.49259976, 3.56787645,  -0.96048609, -1.65009008, 0.91913369 ],
        [ -1.15086474, 3.24188877,  -0.94676956, -1.33850077, 0.90406579 ],
        [ -0.77071156, 2.88555617,  -0.92708474, -0.99842341, 0.88461601 ],
        [ -0.41960447, 2.56390849,  -0.90436647, -0.69135412, 0.86414181 ],
        [ -0.09693333, 2.27598813,  -0.87964748, -0.41582451, 0.84337452 ],
        [ 0.20618523,  2.01325849,  -0.85342363, -0.16332224, 0.82245847 ],
        [ 0.48666883,  1.77753865,  -0.8271132,  0.06454267,  0.80223342 ],
        [ 0.75875134,  1.55611551,  -0.80035996, 0.2800932,   0.78217494 ],
        [ 1.03136452,  1.34164717,  -0.77303893, 0.49057338,  0.76201166 ],
        [ 1.31205279,  1.12859814,  -0.74509905, 0.7015789,   0.74155627 ],
        [ 1.59997725,  0.91809197,  -0.71733669, 0.91214753,  0.72125052 ],
        [ 1.89547672,  0.71014589,  -0.6903854,  1.12231972,  0.70143089 ],
        [ 2.20315254,  0.50178251,  -0.66443587, 1.33512704,  0.68213352 ],
        [ 2.52688966,  0.29075143,  -0.63974196, 1.55289035,  0.66346329 ],
        [ 2.8711239,   0.0746169,   -0.61651791, 1.77813183,  0.64551655 ],
        [ 3.2412802,   -0.14949113, -0.59493261, 2.01383358,  0.62837267 ],
        [ 3.64424609,  -0.38511407, -0.57511732, 2.26368662,  0.61209449 ],
        [ 4.04635569,  -0.61297527, -0.55873698, 2.50693052,  0.59809629 ],
        [ 4.45405567,  -0.83789995, -0.54510623, 2.74822811,  0.58592097 ],
        [ 4.87136772,  -1.06293069, -0.53377368, 2.99046005,  0.57527768 ],
        [ 5.30865944,  -1.29418506, -0.52425674, 3.23989559,  0.56580527 ],
        [ 5.76354834,  -1.53081477, -0.51644693, 3.49532744,  0.5574836 ],
        [ 6.24166468,  -1.77614823, -0.51008384, 3.76006047,  0.55013257 ],
        [ 6.75200174,  -2.03508659, -0.50494039, 4.03909,     0.54357864 ],
        [ 7.29012761,  -2.30567586, -0.50094606, 4.33001102,  0.53783862 ],
        [ 7.86674773,  -2.59360095, -0.49790282, 4.63862697,  0.53275631 ],
        [ 8.48755617,  -2.90196714, -0.49568754, 4.96792549,  0.52826365 ],
        [ 9.16715901,  -3.23827755, -0.4941713,  5.32552583,  0.52425643 ],
        [ 9.90310592,  -3.60159232, -0.49327422, 5.71001292,  0.52074649 ],
        [ 10.72001722, -4.00435892, -0.49288209, 6.13409133,  0.51761875 ],
        [ 11.62972237, -4.45272212, -0.49291405, 6.60366227,  0.51484707 ],
        [ 12.65969852, -4.96058396, -0.4932915,  7.13261615,  0.51237121 ],
        [ 13.84978902, -5.54802097, -0.49394638, 7.74099869,  0.51013664 ],
        [ 15.27141558, -6.25084854, -0.49482287, 8.46469529,  0.50808051 ],
        [ 15.61040736, -6.41862493, -0.49504306, 8.63685927,  0.50766488 ],
        [ 15.61520987, -6.42100244, -0.49507475, 8.63929732,  0.50765912 ]
    ],

    energy_table => [
        [ -9.96706457, 5.005e-05,  4.015e-05,  5.005e-05,  4.015e-05,  11.438871,   -0.99738195 ],
        [ -9.00486523, 0.00010796, 8.588e-05,  0.00010796, 8.588e-05,  10.478425,   -0.9990199 ],
        [ -8.37943566, 0.00017727, 0.00014008, 0.00017727, 0.00014008, 9.8532663,   -1.0001691 ],
        [ -7.0754681,  0.00049247, 0.00038225, 0.00049247, 0.00038225, 8.5474418,   -1.0027984 ],
        [ -6.14888263, 0.00100492, 0.00076644, 0.00100492, 0.00076644, 7.6173585,   -1.0048462 ],
        [ -5.40661457, 0.00176137, 0.00131903, 0.00176137, 0.00131903, 6.8708563,   -1.006626 ],
        [ -4.77618655, 0.00281162, 0.00206456, 0.00281162, 0.00206456, 6.2357576,   -1.0082401 ],
        [ -4.22635836, 0.00419284, 0.00301396, 0.00419284, 0.00301396, 5.6810001,   -1.0097293 ],
        [ -3.73901064, 0.00592878, 0.00416457, 0.00592878, 0.00416457, 5.1885807,   -1.0111169 ],
        [ -3.29977288, 0.00804156, 0.0055083,  0.00804156, 0.0055083,  4.7441788,   -1.0124254 ],
        [ -2.89677381, 0.01055921, 0.00703562, 0.01055921, 0.00703562, 4.3359252,   -1.0136768 ],
        [ -2.52189095, 0.01350569, 0.00872804, 0.01350569, 0.00872804, 3.9556925,   -1.0148873 ],
        [ -2.16825241, 0.01691007, 0.01056258, 0.01691007, 0.01056258, 3.5965836,   -1.0160723 ],
        [ -1.82773381, 0.02083587, 0.01252375, 0.02083587, 0.01252375, 3.2503942,   -1.0172684 ],
        [ -1.49259976, 0.02537673, 0.01459248, 0.02537673, 0.01459248, 2.9092701,   -1.0183012 ],
        [ -1.15086474, 0.03073429, 0.01676501, 0.03073429, 0.01676501, 2.561132,    -1.018694 ],
        [ -0.77071156, 0.03756025, 0.01912259, 0.03756025, 0.01912259, 2.1738903,   -1.0179467 ],
        [ -0.41960447, 0.0446282,  0.02109165, 0.0446282,  0.02109165, 1.8167089,   -1.0159785 ],
        [ -0.09693333, 0.05168559, 0.02259302, 0.05168559, 0.02259302, 1.489274,    -1.0129932 ],
        [ 0.20618523,  0.05870267, 0.02364193, 0.05870267, 0.02364193, 1.182721,    -1.0092884 ],
        [ 0.48666883,  0.06542755, 0.0242494,  0.06542755, 0.0242494,  0.90016216,  -1.0053195 ],
        [ 0.75875134,  0.07206547, 0.02448588, 0.07206547, 0.02448588, 0.62718116,  -1.0012654 ],
        [ 1.03136452,  0.07873381, 0.02438028, 0.07873381, 0.02438028, 0.35477935,  -0.99730868 ],
        [ 1.31205279,  0.08552285, 0.02394146, 0.08552285, 0.02394146, 0.075400221, -0.99362877 ],
        [ 1.59997725,  0.09231477, 0.02319105, 0.09231477, 0.02319105, -0.21018559, -0.99050327 ],
        [ 1.89547672,  0.09902193, 0.0221671,  0.09902193, 0.0221671,  -0.50246199, -0.98815099 ],
        [ 2.20315254,  0.10565163, 0.02090019, 0.10565163, 0.02090019, -0.80619045, -0.98675362 ],
        [ 2.52688966,  0.11218182, 0.0194242,  0.11218182, 0.0194242,  -1.1254973,  -0.98663825 ],
        [ 2.8711239,   0.11858561, 0.01777409, 0.11858561, 0.01777409, -1.4652502,  -0.98842906 ],
        [ 3.2412802,   0.12483311, 0.01598517, 0.12483311, 0.01598517, -1.8316995,  -0.99259516 ],
        [ 3.64424609,  0.13089043, 0.01409322, 0.13089043, 0.01409322, -2.2328264,  -1.0208211 ],
        [ 4.04635569,  0.13619575, 0.01231592, 0.13619575, 0.01231592, -2.6534941,  -1.056331 ],
        [ 4.45405567,  0.14087394, 0.0106594,  0.14087394, 0.0106594,  -3.0883679,  -1 ],
        [ 4.87136772,  0.14499743, 0.00913166, 0.14499743, 0.00913166, -3.50568,    -1 ],
        [ 5.30865944,  0.14867434, 0.00771645, 0.14867434, 0.00771645, -3.9429717,  -1 ],
        [ 5.76354834,  0.15188663, 0.00643903, 0.15188663, 0.00643903, -4.3978606,  -1 ],
        [ 6.24166468,  0.15468399, 0.00529485, 0.15468399, 0.00529485, -4.875977,   -1 ],
        [ 6.75200174,  0.15711754, 0.0042746,  0.15711754, 0.0042746,  -5.386314,   -1 ],
        [ 7.29012761,  0.15917261, 0.00339413, 0.15917261, 0.00339413, -5.9244399,  -1 ],
        [ 7.86674773,  0.16090334, 0.00263833, 0.16090334, 0.00263833, -6.50106,    -1 ],
        [ 8.48755617,  0.16233529, 0.00200234, 0.16233529, 0.00200234, -7.1218684,  -1 ],
        [ 9.16715901,  0.1635077,  0.00147367, 0.1635077,  0.00147367, -7.8014713,  -1 ],
        [ 9.90310592,  0.16442896, 0.00105261, 0.16442896, 0.00105261, -8.5374182,  -1 ],
        [ 10.72001722, 0.16514526, 0.00072126, 0.16514526, 0.00072126, -9.3543295,  -1 ],
        [ 11.62972237, 0.16567983, 0.00047125, 0.16567983, 0.00047125, -10.264035,  -1 ],
        [ 12.65969852, 0.16606429, 0.00028966, 0.16606429, 0.00028966, -11.294011,  -1 ],
        [ 13.84978902, 0.16632744, 0.00016421, 0.16632744, 0.00016421, -12.484101,  -1 ],
        [ 15.27141558, 0.16649661, 8.287e-05,  0.16649661, 8.287e-05,  -13.905728,  -1 ],
        [ 15.61040736, 0.16652252, 7.034e-05,  0.16652252, 7.034e-05,  -14.24472,   -1 ],
        [ 15.61520987, 0.16652285, 7.018e-05,  0.16652285, 7.018e-05,  -14.249522,  -1 ]
    ],

    impulse_table => [
        [
            -11.781497,   13.814886,  -11.782371, 4.3257353, 0, 0, 0, -7.282151,
            0.0032372486, -83.001007, 0.33605833, 1.0000577, -9.7047894
        ],
        [
            -2.302585,    4.3552224,  -2.4064279, 4.1806159,  0, 0, 0, -2.6139111,
            0.0032372525, -82.901014, 0.33095928, 0.90680299, -2.5558732
        ],
        [
            -1.6094379,   3.6803302,  -1.7577537, 4.1322016,  0, 0, 0, -2.3283704,
            0.0032372668, -82.801014, 0.32613854, 0.85775245, -2.1335171
        ],
        [
            -1.3862943,   3.4659754,  -1.5526156, 4.1152375,  0, 0, 0, -2.2449173,
            0.0032372783, -82.751014, 0.32382393, 0.83825713, -2.0052291
        ],
        [
            -1.2039728,   3.2922329,  -1.3865802, 4.1011607,  0, 0, 0, -2.1805936,
            0.0032372928, -82.701014, 0.32156903, 0.82092153, -1.9035384
        ],
        [
            -0.91629068,  3.0211204,  -1.1277769, 4.0790339,  0, 0, 0, -2.0871738,
            0.0032373313, -82.601014, 0.31722729, 0.79099521, -1.7492139
        ],
        [
            -0.69314713,  2.8138258,  -0.92997552, 4.0624276, 0, 0, 0, -2.0223851,
            0.0032373828, -82.501014, 0.31309308,  0.7656694, -1.6350141
        ],
        [
            -0.51082558,  2.6466968,  -0.77043621, 4.0495739,  0, 0, 0, -1.9749725,
            0.0032374477, -82.401014, 0.30914874,  0.74368027, -1.5454787
        ],
        [
            -0.3566749,   2.5071393,  -0.63709687, 4.0394202,  0, 0, 0, -1.9390496,
            0.0032375263, -82.301014, 0.30537884,  0.72424264, -1.4725279
        ],
        [
            -0.2231435,   2.3876511,  -0.52279533, 4.0312857,  0, 0, 0, -1.9111769,
            0.0032376186, -82.201014, 0.30176976,  0.70682871, -1.4114285
        ],
        [
            4.8035232e-08, 2.1911113,  -0.33439035,  4.0193449,  0,          0,
            0,             -1.8716305, 0.0032378454, -82.001014, 0.29498702, 0.67667038,
            -1.3137676
        ],
        [
            0.26236431,   1.965458,   -0.11722948, 4.0084317,  0, 0, 0, -1.8368332,
            0.0032382923, -81.701014, 0.28575602,  0.63982403, -1.2061994
        ],
        [
            0.5306283,    1.7412728,  0.099737587, 4.001261,   0, 0, 0, -1.8146665,
            0.0032390897, -81.301014, 0.2748918,   0.60106484, -1.1044301
        ],
        [
            0.69314723,   1.6088369,  0.22861984, 3.9990069,  0, 0, 0, -1.8078533,
            0.0032398404, -81.001014, 0.26762725, 0.57727537, -1.0468198
        ],
        [
            0.99325182,   1.3711828,  0.46147748, 3.9990186,  0, 0, 0, -1.8081361,
            0.0032421039, -80.301014, 0.25296792, 0.53321022, -0.94835588
        ],
        [
            1.0986123,    1.2898889,  0.54165057, 4.0002665,  0, 0, 0, -1.8120796,
            0.0032432944, -80.001014, 0.2474787,  0.51780841, -0.91618066
        ],
        [
            1.2089604,   1.2059373,  0.62474631, 4.0022284,  0, 0, 0, -1.8182759,
            0.003244851, -79.651014, 0.24156447, 0.50176873, -0.88378733
        ],
        [
            1.3862944,    1.0735514,  0.75643553, 4.0066966,  0, 0, 0, -1.8325135,
            0.0032482228, -79.001014, 0.23175559, 0.47626788, -0.83447053
        ],
        [
            1.609438,     0.91130967, 0.91896798, 4.0143813,  0, 0, 0, -1.8575545,
            0.0032546373, -78.001014, 0.21900289, 0.44482876, -0.77707231
        ],
        [
            1.7917595,    0.7822287,  1.0492166,  4.0221025,  0, 0, 0, -1.883534,
            0.0032625507, -77.001014, 0.20836714, 0.41981661, -0.73388276
        ],
        [
            1.9459102,    0.67543889, 1.1576129, 4.0294457,  0, 0, 0, -1.9091107,
            0.0032719789, -76.001014, 0.1993072, 0.39922233, -0.69986642
        ],
        [
            2.0794416,    0.58460477, 1.250273,   4.0362781,  0, 0, 0, -1.933778,
            0.0032829409, -75.001014, 0.19145983, 0.38183349, -0.67218817
        ],
        [
            2.3025851,    0.43610972, 1.4026578,  4.0483544, 0, 0, 0, -1.9798433,
            0.0033095581, -73.001014, 0.17845415, 0.3537832, -0.62950219
        ],
        [
            2.7080502,    0.17601065, 1.6721978,  4.0706922,  0, 0, 0, -2.0774269,
            0.0034049381, -68.001014, 0.15560503, 0.30629108, -0.562785
        ],
        [
            2.9957323, -0.0017260117, 1.8581913,    4.084825,   0,          0,
            0,         -2.1561146,    0.0035459788, -63.001013, 0.14032594, 0.27544767,
            -0.5235511
        ],
        [
            3.4011974, -0.24395988, 2.1137769,    4.0966932,  0,          0,
            0,         -2.2782673,  0.0040023803, -53.001013, 0.12045188, 0.23604289,
            -0.4797311
        ],
        [
            3.9120231, -0.53757861, 2.426291,     4.0685995,  0,           0,
            0,         -2.4479031,  0.0062828718, -33.001012, 0.098447129, 0.19294276,
            -0.4490595
        ],
        [
            4.2484953,   -0.72518881, 2.6271876,   3.9184169,  0, 0, 0, -2.5691918,
            0.017249939, -13.001011,  0.085812468, 0.16801828, -0.47133833
        ]
    ],

    blast_info => [ 0, 0, 0, 0, 0, 0, 4.3223761, 0, 0, 0, 0, 0, 0.23659643, 0.33605835, 0, 0, 0, 0, 0, 0 ],

};
1;
