package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=3.5
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S3x5'} = {

    table_name  => 'S3x5',
    symmetry    => 2,
    gamma       => 3.5,
    data_source => 'S4000_G3x5/moc2_from_moc_r200/',

    shock_table_info => [ 2, 3.5, 1.14e-06, 4000, 9.38e-07, 10, 2.04e-07, 4.36e-07, 5e-07, 152.06, -0.80648 ],

    shock_table => [
        [ -4.10855906,  12.00033666,   -2.99997134, -4.10979624, 0.9981431 ],
        [ -3.50081136,  10.17713378,   -2.99985496, -3.50389266, 0.99537114 ],
        [ -3.09639265,  8.96399061,    -2.99951961, -3.10205146, 0.99148913 ],
        [ -2.78893591,  8.04186367,    -2.99879187, -2.79792432, 0.98646217 ],
        [ -2.5197713,   7.23487089,    -2.99729834, -2.53325758, 0.9796523 ],
        [ -2.32584739,  6.65381157,    -2.9951648,  -2.34392155, 0.97268723 ],
        [ -2.15098412,  6.13033292,    -2.99184875, -2.17452972, 0.96436156 ],
        [ -1.98437974,  5.63227756,    -2.98663081, -2.01468475, 0.95405668 ],
        [ -1.83770518,  5.19470852,    -2.97938764, -1.87556255, 0.94253203 ],
        [ -1.70200323,  4.79103761,    -2.96934084, -1.74852506, 0.92931696 ],
        [ -1.57607231,  4.41790446,    -2.95587099, -1.63240582, 0.91438943 ],
        [ -1.45716515,  4.06743096,    -2.93809073, -1.52465316, 0.89750628 ],
        [ -1.34285188,  3.73283166,    -2.91484401, -1.42312323, 0.87832531 ],
        [ -1.23066264,  3.40744033,    -2.88453627, -1.32578727, 0.85633573 ],
        [ -1.11641599,  3.08009178,    -2.84430332, -1.22939908, 0.83042246 ],
        [ -0.99079091,  2.72621133,    -2.78722768, -1.12709013, 0.79761624 ],
        [ -0.87240039,  2.40008386,    -2.71987871, -1.0346948,  0.76259561 ],
        [ -0.76379404,  2.10857919,    -2.64640389, -0.95376527, 0.72725746 ],
        [ -0.66179696,  1.84257665,    -2.56810027, -0.88138161, 0.6917431 ],
        [ -0.56550679,  1.59913225,    -2.48743859, -0.81645104, 0.65670779 ],
        [ -0.47208569,  1.3705876,     -2.40479405, -0.75672225, 0.62191088 ],
        [ -0.37870662,  1.14998369,    -2.31990951, -0.7002823,  0.58694315 ],
        [ -0.28524114,  0.93714034,    -2.23468768, -0.64704572, 0.5523395 ],
        [ -0.1908976,   0.73031427,    -2.15024461, -0.59654906, 0.51834038 ],
        [ -0.09460744,  0.52729373,    -2.06724974, -0.54825366, 0.48505061 ],
        [ 0.00425462,   0.32695013,    -1.9865877,  -0.50191589, 0.45270419 ],
        [ 0.10659091,   0.1276841,     -1.90875426, -0.45720847, 0.4214168 ],
        [ 0.21340582,   -0.07215054,   -1.83407571, -0.41382809, 0.39126707 ],
        [ 0.32559506,   -0.27385258,   -1.7628943,  -0.37157962, 0.36236382 ],
        [ 0.44405436,   -0.47861517,   -1.6954881,  -0.33031558, 0.33480776 ],
        [ 0.56981605,   -0.68776781,   -1.63200984, -0.28988585, 0.30866069 ],
        [ 0.70411972,   -0.90286628,   -1.57249744, -0.25012638, 0.28394434 ],
        [ 0.84810711,   -1.12519662,   -1.51703275, -0.21095302, 0.26070255 ],
        [ 0.99929784,   -1.35066411,   -1.46678145, -0.1731819,  0.23944304 ],
        [ 1.15638633,   -1.57746443,   -1.4219088,  -0.13711011, 0.22026613 ],
        [ 1.32070166,   -1.80771455,   -1.38167176, -0.10237929, 0.20288523 ],
        [ 1.49355707,   -2.04333364,   -1.34547914, -0.06870948, 0.18707154 ],
        [ 1.67594525,   -2.28568348,   -1.31290754, -0.03593601, 0.1726641 ],
        [ 1.8692394,    -2.53654525,   -1.28353654, -0.00386497, 0.15949937 ],
        [ 2.07467019,   -2.79742677,   -1.25704036, 0.02763274,  0.14745257 ],
        [ 2.29351915,   -3.06984002,   -1.2331351,  0.0586644,   0.1364155 ],
        [ 2.5270194,    -3.35518975,   -1.21157921, 0.08930625,  0.12629763 ],
        [ 2.77741085,   -3.65604851,   -1.19208525, 0.11973435,  0.11698448 ],
        [ 3.04590652,   -3.97368831,   -1.174493,   0.14996493,  0.10841997 ],
        [ 3.33458564,   -4.31037945,   -1.1586032,  0.18009506,  0.10052812 ],
        [ 3.64551859,   -4.66833021,   -1.14425003, 0.21019143,  0.09324761 ],
        [ 3.98096893,   -5.04992985,   -1.13128318, 0.24031432,  0.08652345 ],
        [ 4.34359341,   -5.45797188,   -1.11956073, 0.27053275,  0.08030338 ],
        [ 4.73604104,   -5.89519875,   -1.10896256, 0.30088827,  0.07454502 ],
        [ 5.16195431,   -6.36541455,   -1.09936265, 0.33147023,  0.06920105 ],
        [ 5.62436778,   -6.87170389,   -1.09066943, 0.36229301,  0.06424131 ],
        [ 6.12750986,   -7.41842145,   -1.08278313, 0.3934249,   0.05962945 ],
        [ 6.67580261,   -8.0100802,    -1.0756198,  0.42491137,  0.05533614 ],
        [ 7.27406268,   -8.65157277,   -1.06910536, 0.45678868,  0.051336 ],
        [ 7.92770206,   -9.34838449,   -1.0631729,  0.48909262,  0.04760598 ],
        [ 8.64252903,   -10.1063775,   -1.05776438, 0.52184604,  0.04412653 ],
        [ 9.42514897,   -10.93221447,  -1.05282698, 0.5550758,   0.04087937 ],
        [ 10.2821652,   -11.83251359,  -1.04831763, 0.58877745,  0.03785057 ],
        [ 11.2213796,   -12.81511416,  -1.04419498, 0.62296482,  0.03502544 ],
        [ 12.25059319,  -13.88781546,  -1.0404253,  0.65762195,  0.03239262 ],
        [ 13.37820646,  -15.05901082,  -1.03697817, 0.69272871,  0.0299414 ],
        [ 14.61421962,  -16.33872146,  -1.03382383, 0.72828819,  0.02765994 ],
        [ 15.96883278,  -17.73713913,  -1.03093748, 0.76428043,  0.02553854 ],
        [ 17.45244594,  -19.26463453,  -1.02829766, 0.80066821,  0.02356868 ],
        [ 19.07745912,  -20.93361077,  -1.02588271, 0.83744095,  0.02174063 ],
        [ 20.72327136,  -22.62028183,  -1.02382584, 0.87188952,  0.02016353 ],
        [ 22.64941625,  -24.5903118,   -1.02180181, 0.90917148,  0.01859157 ],
        [ 24.65945181,  -26.64233634,  -1.02002801, 0.94510331,  0.0171979 ],
        [ 26.94111216,  -28.96770691,  -1.0183364,  0.98276723,  0.01585387 ],
        [ 29.23771107,  -31.30472387,  -1.01690097, 1.01781982,  0.0147014 ],
        [ 31.8317953,   -33.9408206,   -1.01552893, 1.05447533,  0.01358919 ],
        [ 34.77619888,  -36.92896349,  -1.01421974, 1.09286597,  0.01251771 ],
        [ 37.69068994,  -39.88326016,  -1.01312527, 1.12799823,  0.01161395 ],
        [ 40.9762406,   -43.21016591,  -1.01207785, 1.16468619,  0.01074203 ],
        [ 44.69773226,  -46.97468926,  -1.01107718, 1.20305683,  0.00990225 ],
        [ 48.29444056,  -50.60972466,  -1.01025627, 1.23739353,  0.00920823 ],
        [ 52.33065689,  -54.68570676,  -1.00946918, 1.27317324,  0.00853822 ],
        [ 56.88015869,  -59.27652682,  -1.00871557, 1.31050969,  0.00789238 ],
        [ 62.03281336,  -64.47218031,  -1.00799517, 1.34953085,  0.00727088 ],
        [ 66.86630314,  -69.34289002,  -1.00742001, 1.38343916,  0.00677168 ],
        [ 72.27232648,  -74.78749449,  -1.00686763, 1.41871146,  0.00628959 ],
        [ 78.34411971,  -80.89933502,  -1.00633783, 1.45545265,  0.0058247 ],
        [ 85.19463163,  -87.79147813,  -1.00583048, 1.49378037,  0.0053771 ],
        [ 92.96173631,  -95.60193042,  -1.00534543, 1.53382708,  0.00494687 ],
        [ 101.81512411, -104.50053232, -1.00488253, 1.57574277,  0.00453409 ],
        [ 109.81989061, -112.5429279,  -1.00452806, 1.61073487,  0.00421649 ],
        [ 118.78555882, -121.54762713, -1.00418759, 1.64713043,  0.00391015 ],
        [ 128.87116815, -131.67377966, -1.00386106, 1.68503965,  0.00361512 ],
        [ 140.2699852,  -143.114776,   -1.00354839, 1.72458638,  0.00333142 ],
        [ 153.21872109, -156.10746707, -1.00324949, 1.76591039,  0.00305909 ],
        [ 164.12078314, -167.04375558, -1.00303432, 1.79816435,  0.00286234 ],
        [ 176.20674773, -179.16510934, -1.00282682, 1.83158176,  0.00267201 ],
        [ 189.65344529, -192.64844183, -1.00262697, 1.86624564,  0.00248813 ],
        [ 204.67192135, -207.70489088, -1.00243474, 1.90224828,  0.00231072 ]
    ],

    energy_table => [
        [ -4.10855906,  0.00017084, 0.00035971, 0.00017084, 0.00035971,    5.1996244,    -1.5723208 ],
        [ -3.50081136,  0.00060964, 0.00126698, 0.00060964, 0.00126698,    4.2362146,    -1.6025964 ],
        [ -3.09639265,  0.00140498, 0.00287842, 0.00140498, 0.00287842,    3.5834166,    -1.629627 ],
        [ -2.78893591,  0.00262542, 0.00529434, 0.00262542, 0.00529434,    3.0787624,    -1.6561544 ],
        [ -2.5197713,   0.00449578, 0.00889437, 0.00449578, 0.00889437,    2.6295034,    -1.684219 ],
        [ -2.32584739,  0.00657637, 0.01277759, 0.00657637, 0.01277759,    2.3007788,    -1.7064697 ],
        [ -2.15098412,  0.0092065,  0.01752564, 0.0092065,  0.01752564,    2.0005906,    -1.7277338 ],
        [ -1.98437974,  0.01259603, 0.02340309, 0.01259603, 0.02340309,    1.7109915,    -1.7493719 ],
        [ -1.83770518,  0.01648519, 0.02983369, 0.01648519, 0.02983369,    1.4529662,    -1.7690132 ],
        [ -1.70200323,  0.02100091, 0.0368995,  0.02100091, 0.0368995,     1.2116717,    -1.7868726 ],
        [ -1.57607231,  0.02611039, 0.04439015, 0.02611039, 0.04439015,    0.98562649,   -1.8028972 ],
        [ -1.45716515,  0.03184537, 0.05217067, 0.03184537, 0.05217067,    0.77036173,   -1.8168352 ],
        [ -1.34285188,  0.03825908, 0.06009154, 0.03825908, 0.06009154,    0.56196175,   -1.828349 ],
        [ -1.23066264,  0.04544411, 0.06798584, 0.04544411, 0.06798584,    0.35625836,   -1.8372693 ],
        [ -1.11641599,  0.05365907, 0.0757375,  0.05365907, 0.0757375,     0.14592254,   -1.8429666 ],
        [ -0.99079091,  0.06366419, 0.08332928, 0.06366419, 0.08332928,    -0.085862584, -1.8450165 ],
        [ -0.87240039,  0.07388118, 0.08897585, 0.07388118, 0.08897585,    -0.3042908,   -1.8434079 ],
        [ -0.76379404,  0.08374783, 0.09241289, 0.08374783, 0.09241289,    -0.50433991,  -1.8401266 ],
        [ -0.66179696,  0.09326332, 0.09387439, 0.09326332, 0.09387439,    -0.69185133,  -1.8373649 ],
        [ -0.56550679,  0.1023028,  0.09361799, 0.1023028,  0.09361799,    -0.86867743,  -1.8368624 ],
        [ -0.47208569,  0.11097994, 0.09191993, 0.11097994, 0.09191993,    -1.0403221,   -1.8396147 ],
        [ -0.37870662,  0.11943438, 0.08896695, 0.11943438, 0.08896695,    -1.2123176,   -1.8460677 ],
        [ -0.28524114,  0.12757097, 0.08499502, 0.12757097, 0.08499502,    -1.3852506,   -1.8559986 ],
        [ -0.1908976,   0.13536996, 0.08023514, 0.13536996, 0.08023514,    -1.560901,    -1.8687374 ],
        [ -0.09460744,  0.14284138, 0.07489215, 0.14284138, 0.07489215,    -1.7415221,   -1.8834153 ],
        [ 0.00425462,   0.14996377, 0.06917794, 0.14996377, 0.06917794,    -1.9284934,   -1.8990811 ],
        [ 0.10659091,   0.15673991, 0.06327018, 0.15673991, 0.06327018,    -2.1236696,   -1.9148086 ],
        [ 0.21340582,   0.16317761, 0.0573198,  0.16317761, 0.0573198,     -2.3290475,   -1.929935 ],
        [ 0.32559506,   0.16927547, 0.05146352, 0.16927547, 0.05146352,    -2.5464125,   -1.9439688 ],
        [ 0.44405436,   0.17503161, 0.04581671, 0.17503161, 0.04581671,    -2.7775047,   -1.9564779 ],
        [ 0.56981605,   0.1804502,  0.04046668, 0.1804502,  0.04046668,    -3.0243108,   -1.9672926 ],
        [ 0.70411972,   0.18554162, 0.03547323, 0.18554162, 0.03547323,    -3.2892142,   -1.9763144 ],
        [ 0.84810711,   0.19030987, 0.03088276, 0.19030987, 0.03088276,    -3.5743795,   -1.9836776 ],
        [ 0.99929784,   0.19466237, 0.02681227, 0.19466237, 0.02681227,    -3.8747995,   -1.9893288 ],
        [ 1.15638633,   0.19858807, 0.02327604, 0.19858807, 0.02327604,    -4.1876759,   -1.9934 ],
        [ 1.32070166,   0.20215241, 0.02020489, 0.20215241, 0.02020489,    -4.5155107,   -1.9963288 ],
        [ 1.49355707,   0.20540701, 0.01753858, 0.20540701, 0.01753858,    -4.8608001,   -1.9983536 ],
        [ 1.67594525,   0.20838818, 0.01522857, 0.20838818, 0.01522857,    -5.2254285,   -1.9996794 ],
        [ 1.8692394,    0.21113179, 0.01322712, 0.21113179, 0.01322712,    -5.6120553,   -2.0004813 ],
        [ 2.07467019,   0.21366496, 0.01149461, 0.21366496, 0.01149461,    -6.0230774,   -2.0009121 ],
        [ 2.29351915,   0.21601083, 0.00999602, 0.21601083, 0.00999602,    -6.4610054,   -2.0010651 ],
        [ 2.5270194,    0.21818837, 0.00870101, 0.21818837, 0.00870101,    -6.928258,    -2.0010682 ],
        [ 2.77741085,   0.22022152, 0.00757877, 0.22022152, 0.00757877,    -7.4293053,   -2.0010124 ],
        [ 3.04590652,   0.2221214,  0.00660811, 0.2221214,  0.00660811,    -7.9665559,   -2.0008663 ],
        [ 3.33458564,   0.2239033,  0.00576747, 0.2239033,  0.00576747,    -8.5441334,   -2.000692 ],
        [ 3.64551859,   0.22557922, 0.00503892, 0.22557922, 0.00503892,    -9.1661918,   -2.0005514 ],
        [ 3.98096893,   0.22715968, 0.00440696, 0.22715968, 0.00440696,    -9.8372529,   -2.0004116 ],
        [ 4.34359341,   0.22865458, 0.00385797, 0.22865458, 0.00385797,    -10.562625,   -2.0002898 ],
        [ 4.73604104,   0.23007155, 0.00338063, 0.23007155, 0.00338063,    -11.347613,   -2.0001988 ],
        [ 5.16195431,   0.23141955, 0.00296448, 0.23141955, 0.00296448,    -12.199507,   -2.0001317 ],
        [ 5.62436778,   0.23270339, 0.00260152, 0.23270339, 0.00260152,    -13.124381,   -2.0000829 ],
        [ 6.12750986,   0.2339296,  0.00228421, 0.2339296,  0.00228421,    -14.130696,   -2.0000497 ],
        [ 6.67580261,   0.23510307, 0.00200636, 0.23526937, 0.0024584143,  -15.227302,   -2.0000289 ],
        [ 7.27406268,   0.23622787, 0.00176178, 0.23669201, 0.002316487,   -16.423834,   -2.0000158 ],
        [ 7.92770206,   0.23730771, 0.00154901, 0.23813963, 0.00212196,    -17.73112,    -2.0000082 ],
        [ 8.64252903,   0.23834528, 0.00136086, 0.23959519, 0.0019564382,  -19.160778,   -2.000004 ],
        [ 9.42514897,   0.2393432,  0.00119538, 0.24106422, 0.0018007084,  -20.72602,    -2.0000019 ],
        [ 10.2821652,   0.24030301, 0.00104982, 0.24254502, 0.0016566598,  -22.440053,   -2.0000009 ],
        [ 11.2213796,   0.24122665, 0.00092169, 0.24403164, 0.0015144578,  -24.318483,   -2.0000004 ],
        [ 12.25059319,  0.24211513, 0.00080896, 0.24551757, 0.0013776758,  -26.37691,    -2.0000002 ],
        [ 13.37820646,  0.24296938, 0.00070981, 0.24699584, 0.0012474545,  -28.632137,   -2.0000002 ],
        [ 14.61421962,  0.24379084, 0.0006226,  0.24845877, 0.0011227255,  -31.104163,   -2.0000001 ],
        [ 15.96883278,  0.24458039, 0.00054593, 0.24990024, 0.0010093767,  -33.81339,    -2.0000001 ],
        [ 17.45244594,  0.24533856, 0.00047859, 0.25131606, 0.00090212082, -36.780616,   -2.0000001 ],
        [ 19.07745912,  0.24606646, 0.00041946, 0.25273541, 0.00085358263, -40.030643,   -2 ],
        [ 20.72327136,  0.24671562, 0.00037107, 0.25405573, 0.00075429162, -43.322267,   -2 ],
        [ 22.64941625,  0.24738466, 0.00032531, 0.25541503, 0.0006605787,  -47.174557,   -2 ],
        [ 24.65945181,  0.24799847, 0.00028681, 0.25666092, 0.00058188107, -51.194628,   -2 ],
        [ 26.94111216,  0.24861116, 0.00025157, 0.2579035,  0.00050994977, -55.757949,   -2 ],
        [ 29.23771107,  0.24915477, 0.00022284, 0.25900514, 0.00045141046, -60.351147,   -2 ],
        [ 31.8317953,   0.24969738, 0.00019645, 0.26010401, 0.00039770012, -65.539315,   -2 ],
        [ 34.77619888,  0.25023891, 0.0001723,  0.26120004, 0.00034858613, -71.428122,   -2 ],
        [ 37.69068994,  0.25071183, 0.0001529,  0.26215668, 0.0003091999,  -77.257104,   -2 ],
        [ 40.9762406,   0.25118384, 0.00013506, 0.26311104, 0.00027299861, -83.828206,   -2 ],
        [ 44.69773226,  0.25165491, 0.00011871, 0.26406307, 0.00023983447, -91.271189,   -2 ],
        [ 48.29444056,  0.2520579,  0.00010581, 0.26487719, 0.00021371085, -98.464606,   -2 ],
        [ 52.33065689,  0.25246014, 9.392e-05,  0.26568954, 0.00018961795, -106.53704,   -2 ],
        [ 56.88015869,  0.2528616,  8.297e-05,  0.26650006, 0.00016746289, -115.63604,   -2 ],
        [ 62.03281336,  0.25326226, 7.293e-05,  0.26730874, 0.00014715545, -125.94135,   -2 ],
        [ 66.86630314,  0.25359551, 6.522e-05,  0.2679812,  0.00013157965, -135.60833,   -2 ],
        [ 72.27232648,  0.25392818, 5.81e-05,   0.26865234, 0.00011717516, -146.42038,   -2 ],
        [ 78.34411971,  0.25426024, 5.152e-05,  0.26932212, 0.00010389076, -158.56396,   -2 ],
        [ 85.19463163,  0.25459169, 4.547e-05,  0.26999054, 9.1677011e-05, -172.26499,   -2 ],
        [ 92.96173631,  0.2549225,  3.993e-05,  0.27065756, 8.0484668e-05, -187.7992,    -2 ],
        [ 101.81512411, 0.25525267, 3.486e-05,  0.27132318, 7.0265023e-05, -205.50597,   -2 ],
        [ 109.81989061, 0.25551633, 3.114e-05,  0.27185464, 6.2757314e-05, -221.51551,   -2 ],
        [ 118.78555882, 0.25577955, 2.77e-05,   0.27238517, 5.5816856e-05, -239.44684,   -2 ],
        [ 128.87116815, 0.25604234, 2.453e-05,  0.27291477, 4.9419697e-05, -259.61806,   -2 ],
        [ 140.2699852,  0.25630469, 2.162e-05,  0.27344342, 4.3541826e-05, -282.41569,   -2 ],
        [ 153.21872109, 0.25656658, 1.895e-05,  0.2739711,  3.8159353e-05, -308.31317,   -2 ],
        [ 164.12078314, 0.25676269, 1.71e-05,   0.27436623, 3.4433631e-05, -330.11729,   -2 ],
        [ 176.20674773, 0.25695855, 1.538e-05,  0.27476081, 3.0963516e-05, -354.28922,   -2 ],
        [ 189.65344529, 0.25715414, 1.378e-05,  0.27515482, 2.7739389e-05, -381.18262,   -2 ],
        [ 204.67192135, 0.25734946, 1.229e-05,  0.27554828, 2.4751579e-05, -411.21957,   -2 ]
    ],

    impulse_table => [
        [
            -4.7133956, 13.814838,   -4.7138948,    0.037498554, -0.0031767721, -0.083799349,
            -2.542722,  -0.26516414, -0.0041885104, -0.31013512, 0.43092908,    0.99983223,
            -3.7428665
        ],
        [
            -4.6051695, 13.490161,   -4.6057567,    0.039360313, -0.0035398739, -0.082773995,
            -2.5416964, -0.31927979, -0.0046672536, -0.30910937, 0.43092892,    0.99978932,
            -3.6775331
        ],
        [
            -4.1997044, 12.27377,    -4.2007834,    0.047039706, -0.0053097814, -0.077777625,
            -2.5366974, -0.52203176, -0.0070008513, -0.30410936, 0.4309278,     0.99950531,
            -3.4377165
        ],
        [
            -3.9120223, 11.410733,   -3.9136814,    0.053181796, -0.007079627, -0.072785176,
            -2.5316995, -0.66590919, -0.0093343881, -0.29910936, 0.4309257,    0.9990949,
            -3.2730448
        ],
        [
            -3.5065572, 10.194371,   -3.509612,   0.062805763, -0.0106189, -0.062818539,
            -2.521709,  -0.86878755, -0.01400105, -0.28924593, 0.43091716, 0.99788982,
            -3.0504861
        ],
        [
            -2.9957316, 8.6620637,  -3.0023149,   0.076340062, -0.017693052, -0.043007771,
            -2.5017629, -1.1249223, -0.023330045, -0.26938254, 0.43087336,   0.99395258,
            -2.7906253
        ],
        [
            -2.3025844, 6.5841391, -2.3213053,   0.095931432, -0.035283169, 0.0050938526,
            -2.4523107, -1.477408, -0.046558087, -0.22116298, 0.43048342,   0.97581525,
            -2.4924982
        ],
        [
            -1.8971193, 5.3718269,  -1.9317123,   0.10684509,  -0.052480905, 0.049691129,
            -2.4039421, -1.6946972, -0.069388666, -0.17630948, 0.42943303,   0.94769437,
            -2.3639836
        ],
        [
            -1.6094372, 4.516595,   -1.6629854,  0.11381975,  -0.068837444, 0.089467979,
            -2.3573109, -1.8643304, -0.09132509, -0.13658932, 0.4274202,    0.91219644,
            -2.3047628
        ],
        [
            -1.3862937, 3.8596708,  -1.4614471,  0.11868843,  -0.083871314, 0.1239946,
            -2.3130695, -2.0146125, -0.11175141, -0.10316632, 0.42419186,   0.87193031,
            -2.2840195
        ],
        [
            -1.2039721, 3.3305625,  -1.3030072,  0.1223641,    -0.097203719, 0.15349807,
            -2.2718053, -2.1582111, -0.13007431, -0.076826391, 0.41957009,   0.82927171,
            -2.2878547
        ],
        [
            -1.0498215, 2.8916026,  -1.1746505,  0.12536075,   -0.10863157, 0.17845765,
            -2.233963,  -2.3013714, -0.14586884, -0.056971328, 0.41347214,  0.78619116,
            -2.3084859
        ],
        [
            -0.91629006, 2.5200476,  -1.0684616, 0.12797374,   -0.11813379, 0.19939126,
            -2.1997843,  -2.4469073, -0.1589659, -0.043439327, 0.40591758,  0.74417183,
            -2.3408975
        ],
        [
            -0.79850703, 2.2008751,  -0.97921281, 0.13035889,   -0.12583711, 0.21680851,
            -2.1692803,  -2.5955191, -0.16945597, -0.034756975, 0.39701997,  0.70421576,
            -2.3815623
        ],
        [
            -0.69314651, 1.9234786,  -0.90324231, 0.13258278, -0.13196372, 0.23121087, -2.1422547, -2.7465959,
            -0.17762485, -0.0303083, 0.38696647,  0.66691075, -2.427881
        ],
        [
            -0.59783633, 1.6799996,  -0.83787469, 0.13466183,   -0.13677693, 0.2430837,
            -2.118378,   -2.8988304, -0.18385917, -0.029009987, 0.37599126,  0.63252438,
            -2.4779021
        ],
        [
            -0.51082496, 1.464421,   -0.78109553, 0.13659091,   -0.14053742, 0.25287742,
            -2.0972562,  -3.0506968, -0.18856032, -0.029832974, 0.3643494,   0.60109905,
            -2.5301581
        ],
        [
            -0.43078225, 1.2720327,  -0.73135497, 0.13836092,   -0.1434762, 0.26098902,
            -2.0784993,  -3.2007718, -0.19208802, -0.032251975, 0.35229506, 0.57253293,
            -2.5835541
        ],
        [
            -0.35667428, 1.0990927, -0.68744116, 0.13996702,   -0.14578363, 0.26775294,
            -2.0617539,  -3.3479,   -0.19473672, -0.035929072, 0.34006588,  0.54664165,
            -2.6372836
        ],
        [
            -0.2876814, 0.94259627, -0.64839467, 0.14141087,   -0.14760932, 0.2734412,
            -2.0467155, -3.4912407, -0.19673396, -0.040072945, 0.32787322,  0.52320132,
            -2.6907627
        ],
        [
            -0.22314288, 0.80010971, -0.61344795, 0.14269967,   -0.14906754, 0.27826983,
            -2.0331317,  -3.6302449, -0.19824948, -0.044659436, 0.31589715,  0.50197651,
            -2.743577
        ],
        [
            -0.16251826, 0.66964535, -0.58198083, 0.14384443, -0.15024429, 0.28240812, -2.0207953, -3.7646019,
            -0.19940787, -0.04968,   0.30428476,  0.48273712, -2.7954401
        ],
        [
            -0.10535985, 0.54957039, -0.55348866, 0.14485812,   -0.1512039, 0.28598821,
            -2.0095365,  -3.8941822, -0.20029988, -0.054886945, 0.29315091, 0.46526754,
            -2.8461614
        ],
        [
            -0.051292626, 0.43853361, -0.52755741, 0.1457543,    -0.15199453, 0.28911295,
            -1.9992155,   -4.0189851, -0.20099145, -0.059992629, 0.28257905,  0.44937109,
            -2.8956205
        ],
        [
            0.0098536153, 0.3158396, -0.49938616, 0.14669177,   -0.15276942, 0.29236881,
            -1.9879258,   -4.16234,  -0.20162311, -0.066314935, 0.27071898,  0.43215162,
            -2.9531188
        ],
        [
            0.18232222, -0.014815486, -0.42612173, 0.14890178,  -0.15441697, 0.30020385,
            -1.9582547, -4.5759685,   -0.20275806, -0.08504973, 0.23890142,  0.38797527,
            -3.1226913
        ],
        [
            0.26236493, -0.16115704, -0.39499137, 0.14972131,   -0.15497554, 0.30327663,
            -1.9455724, -4.7710365,  -0.20305187, -0.094201861, 0.22535511,  0.36963516,
            -3.2043219
        ],
        [
            0.40546578,  -0.41278401, -0.34340126, 0.15090077, -0.15574447, 0.30805937, -1.9246265, -5.1222609,
            -0.20333985, -0.11101491, 0.20342799,  0.34005621, -3.3535568
        ],
        [
            0.53062892,  -0.62344568, -0.3021336, 0.15167631, -0.15623852, 0.31162166, -1.9081221, -5.4304403,
            -0.20342351, -0.126273,   0.18667312, 0.3172935,  -3.4865459
        ],
        [
            0.69314785,  -0.88558807, -0.25325219, 0.15240177, -0.15670901, 0.3155497, -1.8891915, -5.8299742,
            -0.20339259, -0.1468648,  0.16803644,  0.29156019, -3.6612917
        ],
        [
            0.99325244,  -1.3417912,  -0.17463182, 0.15316069, -0.15726463, 0.32120811, -1.8612995, -6.5611875,
            -0.20314085, -0.18510894, 0.14129749,  0.25339772, -3.9863601
        ],
        [
            1.098613,    -1.4948667,  -0.15002842, 0.15330755, -0.15739978, 0.32280789, -1.8535052, -6.8151202,
            -0.20303029, -0.19835206, 0.13380192,  0.24234671, -4.1004869
        ],
        [
            1.208961,    -1.6518617,  -0.12568417, 0.1534168,  -0.15751949, 0.32430772, -1.8463409, -7.0793465,
            -0.20291257, -0.21265036, 0.12679785,  0.23184979, -4.219786
        ],
        [
            1.386295,   -1.8978661, -0.089278904, 0.15352244,  -0.15767679, 0.32639203,
            -1.8367531, -7.500196,  -0.20272854,  -0.23631478, 0.11707832,  0.21697762,
            -4.4107949
        ],
        [
            1.6094386,  -2.1979954, -0.047584273, 0.15357489,  -0.15783093, 0.3285376,
            -1.8275716, -8.0232445, -0.20251805,  -0.26475916, 0.10702181,  0.20116662,
            -4.6496239
        ],
        [
            1.7917601,  -2.4366689, -0.016416132, 0.15357837,  -0.15793134, 0.32996806,
            -1.8219996, -8.4454711, -0.20236686,  -0.28783139, 0.10021753,  0.19019294,
            -4.8433935
        ],
        [
            1.9459108,  -2.6345532, 0.0081819852, 0.15356578,  -0.158003,   0.33099077,
            -1.8183683, -8.7990932, -0.20225572,  -0.30889592, 0.095254659, 0.18203358,
            -5.0062615
        ],
        [
            2.0794422,   -2.803424,   0.02833577,  0.15354807, -0.15805694, 0.33175838, -1.8158726, -9.1030872,
            -0.20216963, -0.32621265, 0.091442644, 0.17567049, -5.1466575
        ],
        [
            2.302788,    -3.0812655,  0.059926817, 0.15351183, -0.15813892, 0.33283812, -1.8135048, -9.6070998,
            -0.20204873, -0.35520622, 0.085902193, 0.16469116, -5.3881827
        ],
        [
            2.995796,    -3.9147579,  0.1444946,   0.15340089, -0.15829641, 0.334997, -1.8085609, -11.141392,
            -0.20179496, -0.44551237, 0.073390076, 0.14364814, -6.0993127
        ],
        [
            3.9122197,   -4.9720706,  0.23432171,  0.15331103, -0.15839684, 0.33630077, -1.8068054, -13.12085,
            -0.20164008, -0.56464887, 0.063097162, 0.12506519, -7.0291533
        ],
        [
            4.6053126,   -5.7500112,  0.29102535,  0.15327781, -0.15843729, 0.33673814, -1.8071351, -14.593459,
            -0.20159113, -0.65461217, 0.057800084, 0.11506281, -7.7274951
        ],
        [
            5.2986285,   -6.5154806,  0.34082174,  0.15326082, -0.15845357, 0.33695616, -1.8069151, -16.052481,
            -0.20156482, -0.74453449, 0.053733512, 0.10720755, -8.423722
        ],
        [
            6.2147281,   -7.5128062,  0.39859379,  0.15325064,  -0.15848994, 0.3370873, -1.806912, -17.965238,
            -0.20154957, -0.86329174, 0.049559468, 0.099018402, -9.3417388
        ],
        [
            6.9078297,   -8.2593403,  0.43756032,  0.1532473,   -0.15846826, 0.33713106, -1.8069124, -19.404089,
            -0.20154449, -0.95311596, 0.047013414, 0.093977443, -10.035537
        ],
        [
            7.6010849,  -9.0006792, 0.47325541,  0.15324565,  -0.15847007, 0.33715294, -1.8069129, -20.837858,
            -0.2012074, -1.0221407, 0.044841371, 0.089658398, -10.729162
        ],
        [
            8.5173251,   -9.9738867, 0.51628576,  0.15324468,  -0.15847112, 0.33716607, -1.8069132, -22.726281,
            -0.19864494, -1.0899972, 0.042400642, 0.084791689, -11.645637
        ],
        [
            9.2103442, -10.705927, 0.54620525,  0.15324436,  -0.15847163, 0.33717045, -1.8069134, -24.150652,
            -0.195789, -1.1390458, 0.040805764, 0.081606775, -12.33874
        ],
        [
            10.819974,   -12.395638, 0.60867707,  0.15324409,  -0.15847188, 0.33717395, -1.8069135, -27.448686,
            -0.18801413, -1.2466271, 0.037714413, 0.075427892, -13.94844
        ],
        [
            11.512984,   -13.119439, 0.63306246, 0.15324405,  -0.15847169, 0.33717439, -1.8069135, -28.865046,
            -0.18453917, -1.2905921, 0.03658733, 0.073174195, -14.641459
        ],
        [
            13.122404,   -14.793657, 0.68500342,  0.153244,    -0.15847179, 0.33717474, -1.8069135, -32.14783,
            -0.17666597, -1.3881063, 0.034320992, 0.068641893, -16.250887
        ],
        [
            13.815611,   -15.512329, 0.70563759,  0.15324398,  -0.158472, 0.33717479, -1.8069135, -33.559409,
            -0.17342716, -1.4283538, 0.033468117, 0.066936189, -16.944096
        ],
        [
            16.118234,   -17.891141, 0.76807996,  0.15324394,  -0.15847164, 0.33717478, -1.8069136, -38.24005,
            -0.16347873, -1.5556091, 0.031037549, 0.062075093, -19.24672
        ],
        [
            18.420854,   -20.259712, 0.82293825,  0.1532439,   -0.15847217, 0.33717506, -1.8069133, -42.910667,
            -0.15476791, -1.6745133, 0.029073314, 0.058146628, -21.54934
        ]
    ],

    tail_shock_table => [
        [ 5.0295649, -0.80648,    -0.0409965,   -0.0409965 ],
        [ 5.2915753, -0.83067834, -0.050773541, -0.029925462 ],
        [ 6.2117446, -0.91112148, -0.056578589, -0.017983767 ],
        [ 6.9062797, -0.96796161, -0.057517266, -0.013322981 ],
        [ 7.6002821, -1.02204,    -0.057490301, -0.010177573 ],
        [ 8.51699,   -1.0899568,  -0.056757001, -0.0073436382 ],
        [ 9.2101715, -1.1390256,  -0.055940428, -0.0058263853 ],
        [ 10.819937, -1.2466231,  -0.053718474, -0.0035152702 ],
        [ 11.512965, -1.2905901,  -0.052725553, -0.0028523835 ],
        [ 13.1224,   -1.3881059,  -0.050476006, -0.0017672595 ],
        [ 13.815609, -1.4283536,  -0.049550625, -0.0014360169 ],
        [ 16.118234, -1.5557033,  -0.046711278, -0.00069101614 ]
    ],

    blast_info => [
        0.092772635, 2.551696,    0.31914591,  -0.46672578,  3.148402,   0.00064461156,
        2.146195,    -0.35398813, 0.40295223,  0.56748616,   0.65003439, -0.32995987,
        0.11485473,  0.43092995,  0.043783984, -0.045277607, 152.06,     -0.80648,
        1228.3611,   -0.97975694
    ],

};
1;