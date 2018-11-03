package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=2
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S2'} = {

    table_name => 'S2',
    symmetry   => 2,
    gamma      => 2,

    shock_table_info => [ 2, 2, 8.27e-07, 8000, 2.2e-07, 10, 9.7e-08, 2.3e-07, 5e-07, 1148.5, -0.75244 ],

    shock_table => [
        [ -4.35780834,  12.00017668,   -2.99997543, -4.35895378, 0.99828087 ],
        [ -4.01113331,  10.96016446,   -2.99993048, -4.01306074, 0.99710613 ],
        [ -3.62605081,  9.8049597,     -2.99980493, -3.62948762, 0.99483623 ],
        [ -3.19953751,  8.52558181,    -2.9993263,  -3.2060634,  0.99018134 ],
        [ -2.9122093,   7.66390422,    -2.99840635, -2.92226743, 0.98484465 ],
        [ -2.67526991,  6.95363573,    -2.99676138, -2.68964754, 0.978301 ],
        [ -2.47326676,  6.34852441,    -2.9940796,  -2.4927741,  0.97050953 ],
        [ -2.2967624,   5.82038611,    -2.98998867, -2.32224102, 0.96141921 ],
        [ -2.13932535,  5.35008341,    -2.98404295, -2.17166888, 0.95095229 ],
        [ -1.99655705,  4.9246101,     -2.97572164, -2.03672405, 0.93901876 ],
        [ -1.86498626,  4.53379118,    -2.96439525, -1.91403765, 0.925483 ],
        [ -1.74181361,  4.16953664,    -2.94929207, -1.80095868, 0.91015671 ],
        [ -1.62457564,  3.82487188,    -2.92943086, -1.69524431, 0.89276918 ],
        [ -1.51067585,  3.49261958,    -2.90344907, -1.594661,   0.87288058 ],
        [ -1.39702109,  3.16448443,    -2.86928041, -1.49673587, 0.84975059 ],
        [ -1.27790244,  2.8253497,     -2.82282818, -1.39714562, 0.82170117 ],
        [ -1.15123253,  2.47162002,    -2.75978328, -1.29518207, 0.78744932 ],
        [ -1.03679393,  2.15966972,    -2.69003396, -1.20702281, 0.75269976 ],
        [ -0.93052953,  1.87772691,    -2.61481863, -1.12888296, 0.71754578 ],
        [ -0.82984594,  1.61839059,    -2.53553458, -1.05840433, 0.68216971 ],
        [ -0.73421917,  1.37975418,    -2.4547302,  -0.99483161, 0.64727029 ],
        [ -0.63984703,  1.15199677,    -2.37168539, -0.93540027, 0.61218393 ],
        [ -0.54534001,  0.93184338,    -2.28725215, -0.87920875, 0.57700892 ],
        [ -0.45047104,  0.71885667,    -2.20312162, -0.82612437, 0.54223453 ],
        [ -0.3542115,   0.51080716,    -2.12006405, -0.77558659, 0.5080134 ],
        [ -0.25578536,  0.30616814,    -2.03888193, -0.72724534, 0.47456042 ],
        [ -0.15427885,  0.10325196,    -1.96010733, -0.68074375, 0.4420212 ],
        [ -0.04891025,  -0.09923121,   -1.88425778, -0.63584702, 0.41056931 ],
        [ 0.06125039,   -0.30274013,   -1.81164981, -0.59230937, 0.38032043 ],
        [ 0.17703268,   -0.50842947,   -1.74260272, -0.54997678, 0.35140863 ],
        [ 0.29942151,   -0.71763021,   -1.67728838, -0.50868178, 0.32391676 ],
        [ 0.42942464,   -0.9316046,    -1.61583925, -0.46829578, 0.29791588 ],
        [ 0.56825189,   -1.15184134,   -1.55828047, -0.42867337, 0.27343257 ],
        [ 0.71722182,   -1.37988374,   -1.50460517, -0.38968963, 0.25047805 ],
        [ 0.87320319,   -1.61069662,   -1.45609485, -0.35228566, 0.22961661 ],
        [ 1.03557836,   -1.84351443,   -1.41265647, -0.31656415, 0.21082719 ],
        [ 1.20566051,   -2.08037645,   -1.37361306, -0.28218698, 0.19383311 ],
        [ 1.38464933,   -2.32300797,   -1.33844094, -0.24890684, 0.17841908 ],
        [ 1.57374073,   -2.57301604,   -1.30670887, -0.21652746, 0.16440667 ],
        [ 1.7742653,    -2.8320915,    -1.27803916, -0.18487255, 0.15163849 ],
        [ 1.98745291,   -3.10171665,   -1.25213092, -0.15381859, 0.13998945 ],
        [ 2.21459552,   -3.38339658,   -1.22871915, -0.12326051, 0.12934947 ],
        [ 2.45717424,   -3.67881936,   -1.20755879, -0.0930941,  0.11961666 ],
        [ 2.71687029,   -3.98986301,   -1.18842688, -0.06321823, 0.11069843 ],
        [ 2.99531298,   -4.31829579,   -1.17113758, -0.0335635,  0.10251887 ],
        [ 3.29437885,   -4.66613941,   -1.15551537, -0.00405626, 0.09500661 ],
        [ 3.61633693,   -5.03582692,   -1.14139241, 0.02538989,  0.0880934 ],
        [ 3.96344727,   -5.42973528,   -1.12862768, 0.05483355,  0.08172354 ],
        [ 4.33816235,   -5.85042362,   -1.11709248, 0.08432678,  0.07584716 ],
        [ 4.74332726,   -6.30085457,   -1.10666457, 0.11392859,  0.07041728 ],
        [ 5.18217971,   -6.7843848,    -1.09723106, 0.14369963,  0.06539104 ],
        [ 5.65835423,   -7.30476315,   -1.0886892,  0.17369844,  0.06073017 ],
        [ 6.17586712,   -7.8661098,    -1.08094684, 0.2039776,   0.05640118 ],
        [ 6.73914087,   -8.47294015,   -1.0739214,  0.23458322,  0.05237482 ],
        [ 7.35318913,   -9.13036026,   -1.06753756, 0.26556229,  0.0486247 ],
        [ 8.02342155,   -9.84385315,   -1.06172951, 0.29694971,  0.04512854 ],
        [ 8.75584465,   -10.61949197,  -1.05643835, 0.32877682,  0.04186659 ],
        [ 9.55726257,   -11.46414766,  -1.0516109,  0.36107669,  0.03882088 ],
        [ 10.43467781,  -12.38485486,  -1.04720298, 0.39385805,  0.03597731 ],
        [ 11.39589177,  -13.38944649,  -1.04317462, 0.42713005,  0.03332275 ],
        [ 12.44930517,  -14.4863419,   -1.03949092, 0.46089294,  0.03084566 ],
        [ 13.60331836,  -15.68392626,  -1.03612287, 0.4951211,   0.02853687 ],
        [ 14.8681315,   -16.99242068,  -1.03304132, 0.52981755,  0.02638562 ],
        [ 16.25414464,  -18.42221695,  -1.03022193, 0.56496283,  0.02438311 ],
        [ 17.77235779,  -19.984297,    -1.02764302, 0.60052951,  0.02252121 ],
        [ 19.36628238,  -21.62042176,  -1.02537363, 0.63506581,  0.02085765 ],
        [ 21.15253047,  -23.45003047,  -1.02323879, 0.67086642,  0.01927019 ],
        [ 23.05956571,  -25.3995062,   -1.02132576, 0.70620375,  0.01782821 ],
        [ 25.22524864,  -27.60934292,  -1.01950462, 0.74326622,  0.0164375 ],
        [ 27.53291554,  -29.96009015,  -1.01787946, 0.77970926,  0.01518096 ],
        [ 29.97241574,  -32.44139368,  -1.01643342, 0.81532726,  0.01404998 ],
        [ 32.73708992,  -35.24953923,  -1.01505475, 0.85262242,  0.01295983 ],
        [ 35.64709726,  -38.20152971,  -1.01383411, 0.88887872,  0.01198451 ],
        [ 38.94675308,  -41.54485648,  -1.01267012, 0.92683068,  0.01104516 ],
        [ 42.40004628,  -45.04008392,  -1.01164532, 0.96349282,  0.01021029 ],
        [ 46.3152187,   -48.99887979,  -1.01066766, 1.00185083,  0.00940662 ],
        [ 50.38201442,  -53.10727106,  -1.00981241, 1.03862555,  0.00869757 ],
        [ 54.98818325,  -57.75670251,  -1.00899593, 1.07707285,  0.00801515 ],
        [ 59.72602317,  -62.53543816,  -1.00828696, 1.11359578,  0.00741804 ],
        [ 65.08111849,  -67.93304503,  -1.00760939, 1.15173903,  0.00684323 ],
        [ 71.16411895,  -74.06030888,  -1.00696292, 1.19163858,  0.00629086 ],
        [ 77.37332821,  -80.31097995,  -1.00640744, 1.22917669,  0.005813 ],
        [ 84.40455099,  -87.38533511,  -1.00587669, 1.26838892,  0.00535351 ],
        [ 91.46519704,  -94.48583913,  -1.00542551, 1.30476524,  0.00496056 ],
        [ 99.4221953,   -102.48424206, -1.00499356, 1.342691,    0.00458224 ],
        [ 108.43236162, -111.53748753, -1.00458067, 1.38229393,  0.00421861 ],
        [ 117.33047972, -120.47478649, -1.00423493, 1.4184353,   0.00391253 ],
        [ 127.33897586, -130.52396469, -1.00390358, 1.45607987,  0.00361777 ],
        [ 138.64890262, -141.87619698, -1.00358651, 1.49534991,  0.00333436 ],
        [ 151.49409615, -154.76545752, -1.00328364, 1.53638338,  0.00306235 ],
        [ 163.9431978,  -167.25385087, -1.00303525, 1.57307646,  0.00283828 ],
        [ 177.95355597, -181.30502062, -1.00279715, 1.61129081,  0.00262261 ],
        [ 193.79558501, -197.18950605, -1.00256927, 1.65115135,  0.00241537 ],
        [ 211.80070675, -215.23887004, -1.00235156, 1.69279912,  0.00221658 ]
    ],

    energy_table => [
        [ -4.35780834,  0.00235759, 0.00352326, 0.00235759, 0.00352326,    5.4571564,   -1.5959715 ],
        [ -4.01113331,  0.00395555, 0.00589627, 0.00395555, 0.00589627,    4.9004119,   -1.6186003 ],
        [ -3.62605081,  0.00701367, 0.01040344, 0.00701367, 0.01040344,    4.2717083,   -1.6513735 ],
        [ -3.19953751,  0.01316638, 0.01933339, 0.01316638, 0.01933339,    3.5585289,   -1.7011844 ],
        [ -2.9122093,   0.02003166, 0.02907472, 0.02003166, 0.02907472,    3.0641039,   -1.744614 ],
        [ -2.67526991,  0.02818871, 0.04033423, 0.02818871, 0.04033423,    2.6460768,   -1.7903935 ],
        [ -2.47326676,  0.03754864, 0.05281848, 0.03754864, 0.05281848,    2.2799142,   -1.8382202 ],
        [ -2.2967624,   0.04801766, 0.06620862, 0.04801766, 0.06620862,    1.951518,    -1.8868463 ],
        [ -2.13932535,  0.05951688, 0.080188,   0.05951688, 0.080188,      1.6507665,   -1.9331082 ],
        [ -1.99655705,  0.07196488, 0.09442013, 0.07196488, 0.09442013,    1.3718266,   -1.9772121 ],
        [ -1.86498626,  0.0853107,  0.10858566, 0.0853107,  0.10858566,    1.1088438,   -2.0224922 ],
        [ -1.74181361,  0.09953167, 0.12236867, 0.09953167, 0.12236867,    0.85699537,  -2.0603435 ],
        [ -1.62457564,  0.11464817, 0.13545402, 0.11464817, 0.13545402,    0.61369683,  -2.0908061 ],
        [ -1.51067585,  0.1307734,  0.14753378, 0.1307734,  0.14753378,    0.37383348,  -2.1111874 ],
        [ -1.39702109,  0.14816726, 0.1582662,  0.14816726, 0.1582662,     0.13328911,  -2.114018 ],
        [ -1.27790244,  0.16758303, 0.16727044, 0.16758303, 0.16727044,    -0.11822651, -2.0958541 ],
        [ -1.15123253,  0.18920151, 0.17340169, 0.18920151, 0.17340169,    -0.38160453, -2.0506638 ],
        [ -1.03679393,  0.2091867,  0.17525266, 0.2091867,  0.17525266,    -0.61332447, -1.9915513 ],
        [ -0.93052953,  0.22775021, 0.17358414, 0.22775021, 0.17358414,    -0.82167066, -1.9297109 ],
        [ -0.82984594,  0.24502266, 0.16906177, 0.24502266, 0.16906177,    -1.0130104,  -1.8760711 ],
        [ -0.73421917,  0.26088764, 0.1623959,  0.26088764, 0.1623959,     -1.1902028,  -1.8372129 ],
        [ -0.63984703,  0.27582819, 0.15396794, 0.27582819, 0.15396794,    -1.3621181,  -1.8142292 ],
        [ -0.54534001,  0.28992537, 0.14418198, 0.28992537, 0.14418198,    -1.5328709,  -1.8065399 ],
        [ -0.45047104,  0.30310287, 0.13352574, 0.30310287, 0.13352574,    -1.7042332,  -1.8114755 ],
        [ -0.3542115,   0.31541974, 0.12236386, 0.31541974, 0.12236386,    -1.87911,    -1.8253833 ],
        [ -0.25578536,  0.32690376, 0.111036,   0.32690376, 0.111036,      -2.0596471,  -1.8446906 ],
        [ -0.15427885,  0.33759989, 0.09981494, 0.33759989, 0.09981494,    -2.2479893,  -1.8665048 ],
        [ -0.04891025,  0.34753648, 0.0889393,  0.34753648, 0.0889393,     -2.4458674,  -1.8886933 ],
        [ 0.06125039,   0.3567538,  0.07858582, 0.3567538,  0.07858582,    -2.6551631,  -1.9097963 ],
        [ 0.17703268,   0.36527995, 0.068898,   0.36527995, 0.068898,      -2.8774863,  -1.9289039 ],
        [ 0.29942151,   0.37315236, 0.05996582, 0.37315236, 0.05996582,    -3.1146908,  -1.9455568 ],
        [ 0.42942464,   0.3804057,  0.05184368, 0.3804057,  0.05184368,    -3.3686459,  -1.9595805 ],
        [ 0.56825189,   0.38708115, 0.04454541, 0.38708115, 0.04454541,    -3.6415974,  -1.9710315 ],
        [ 0.71722182,   0.39321819, 0.03805953, 0.39321819, 0.03805953,    -3.9360065,  -1.9800831 ],
        [ 0.87320319,   0.39870651, 0.03250108, 0.39870651, 0.03250108,    -4.2454806,  -1.9869056 ],
        [ 1.03557836,   0.40358852, 0.02779511, 0.40358852, 0.02779511,    -4.5685882,  -1.9918325 ],
        [ 1.20566051,   0.40796443, 0.02380287, 0.40796443, 0.02380287,    -4.9077107,  -1.9952414 ],
        [ 1.38464933,   0.41191053, 0.02041227, 0.41191053, 0.02041227,    -5.265094,   -1.997582 ],
        [ 1.57374073,   0.41548788, 0.01752992, 0.41548788, 0.01752992,    -5.6429998,  -1.9991154 ],
        [ 1.7742653,    0.418748,   0.01507628, 0.418748,   0.01507628,    -6.043995,   -2.0000571 ],
        [ 1.98745291,   0.42173106, 0.01298644, 0.42173106, 0.01298644,    -6.4704563,  -2.00057 ],
        [ 2.21459552,   0.42447099, 0.01120521, 0.42447099, 0.01120521,    -6.9249112,  -2.0007994 ],
        [ 2.45717424,   0.42699786, 0.00968505, 0.42699786, 0.00968505,    -7.4102761,  -2.0008467 ],
        [ 2.71687029,   0.42933795, 0.00838549, 0.42933795, 0.00838549,    -7.9298857,  -2.0008147 ],
        [ 2.99531298,   0.43151217, 0.0072733,  0.43151217, 0.0072733,     -8.4869912,  -2.0007063 ],
        [ 3.29437885,   0.43353945, 0.00631987, 0.43353945, 0.00631987,    -9.0853071,  -2.0005605 ],
        [ 3.61633693,   0.43543734, 0.00550054, 0.43543734, 0.00550054,    -9.7293845,  -2.0004453 ],
        [ 3.96344727,   0.43721963, 0.00479519, 0.43721963, 0.00479519,    -10.423739,  -2.0003303 ],
        [ 4.33816235,   0.43889823, 0.00418684, 0.43889823, 0.00418684,    -11.173271,  -2.0002301 ],
        [ 4.74332726,   0.44048408, 0.00366091, 0.44048408, 0.00366091,    -11.983676,  -2.0001559 ],
        [ 5.18217971,   0.44198695, 0.00320507, 0.44198695, 0.00320507,    -12.861435,  -2.0001019 ],
        [ 5.65835423,   0.44341533, 0.00280892, 0.44341533, 0.00280892,    -13.813821,  -2.0000633 ],
        [ 6.17586712,   0.44477642, 0.00246403, 0.44477642, 0.00246403,    -14.848871,  -2.0000373 ],
        [ 6.73914087,   0.44607635, 0.00216332, 0.44607635, 0.00216332,    -15.975434,  -2.0000215 ],
        [ 7.35318913,   0.44732059, 0.00190043, 0.44732059, 0.00190043,    -17.20354,   -2.0000115 ],
        [ 8.02342155,   0.4485133,  0.00166821, 0.44856554, 0.0018434468,  -18.54401,   -2.0000059 ],
        [ 8.75584465,   0.44965814, 0.00146536, 0.44988786, 0.0017766272,  -20.008859,  -2.0000028 ],
        [ 9.55726257,   0.45075846, 0.0012871,  0.45127321, 0.001685257,   -21.611696,  -2.0000013 ],
        [ 10.43467781,  0.4518165,  0.00113033, 0.45271485, 0.0015965979,  -23.366527,  -2.0000005 ],
        [ 11.39589177,  0.45283428, 0.00099241, 0.45419985, 0.0014962999,  -25.288955,  -2.0000002 ],
        [ 12.44930517,  0.45381343, 0.00087104, 0.45571994, 0.0013931231,  -27.395782,  -2.0000001 ],
        [ 13.60331836,  0.45475478, 0.00076431, 0.45726258, 0.0012835384,  -29.703809,  -2.0000001 ],
        [ 14.8681315,   0.45565995, 0.00067044, 0.45881521, 0.0011749562,  -32.233435,  -2.0000001 ],
        [ 16.25414464,  0.4565299,  0.00058791, 0.46038264, 0.001197313,   -35.005462,  -2 ],
        [ 17.77235779,  0.45736542, 0.00051542, 0.46208297, 0.0010481413,  -38.041888,  -2 ],
        [ 19.36628238,  0.45813639, 0.00045414, 0.46364978, 0.00092233707, -41.229737,  -2 ],
        [ 21.15253047,  0.45889631, 0.00039874, 0.46519226, 0.00080888968, -44.802233,  -2 ],
        [ 23.05956571,  0.45960965, 0.00035107, 0.46663861, 0.00071145696, -48.616304,  -2 ],
        [ 25.22524864,  0.46032104, 0.00030752, 0.4680796,  0.00062259533, -52.94767,   -2 ],
        [ 27.53291554,  0.46098612, 0.00027024, 0.46942554, 0.00054665486, -57.563003,  -2 ],
        [ 29.97241574,  0.46160514, 0.00023839, 0.47067731, 0.00048187214, -62.442004,  -2 ],
        [ 32.73708992,  0.46222243, 0.00020924, 0.47192469, 0.0004226667,  -67.971352,  -2 ],
        [ 35.64709726,  0.46279401, 0.00018448, 0.47307898, 0.00037242829, -73.791367,  -2 ],
        [ 38.94675308,  0.46336398, 0.00016183, 0.47422935, 0.00032651724, -80.390679,  -2 ],
        [ 42.40004628,  0.46388861, 0.0001427,  0.47528768, 0.00028778291, -87.297265,  -2 ],
        [ 46.3152187,   0.46441176, 0.00012519, 0.47634254, 0.00025236976, -95.12761,   -2 ],
        [ 50.38201442,  0.46488998, 0.00011051, 0.4773064,  0.00022267893, -103.2612,   -2 ],
        [ 54.98818325,  0.46536688, 9.706e-05,  0.47826725, 0.00019550631, -112.47354,  -2 ],
        [ 59.72602317,  0.46579925, 8.585e-05,  0.4791381,  0.00017288176, -121.94922,  -2 ],
        [ 65.08111849,  0.46623046, 7.557e-05,  0.48000638, 0.00015213924, -132.65941,  -2 ],
        [ 71.16411895,  0.46666051, 6.617e-05,  0.48087207, 0.00013318627, -144.82541,  -2 ],
        [ 77.37332821,  0.46704652, 5.843e-05,  0.48164893, 0.00011758339, -157.24383,  -2 ],
        [ 84.40455099,  0.46743153, 5.134e-05,  0.48242363, 0.0001032909,  -171.30627,  -2 ],
        [ 91.46519704,  0.46777292, 4.556e-05,  0.48311041, 9.1634652e-05, -185.42757,  -2 ],
        [ 99.4221953,   0.46811348, 4.024e-05,  0.48379544, 8.0918462e-05, -201.34156,  -2 ],
        [ 108.43236162, 0.46845322, 3.536e-05,  0.48447869, 7.1097988e-05, -219.3619,   -2 ],
        [ 117.33047972, 0.4687498,  3.144e-05,  0.48507507, 6.320603e-05,  -237.15813,  -2 ],
        [ 127.33897586, 0.46904572, 2.783e-05,  0.48567006, 5.5938056e-05, -257.17512,  -2 ],
        [ 138.64890262, 0.46934098, 2.451e-05,  0.48626365, 4.9265543e-05, -279.79498,  -2 ],
        [ 151.49409615, 0.46963556, 2.147e-05,  0.48685583, 4.3160722e-05, -305.48536,  -2 ],
        [ 163.9431978,  0.46988752, 1.909e-05,  0.48736227, 3.8358663e-05, -330.38357,  -2 ],
        [ 177.95355597, 0.47013897, 1.689e-05,  0.48786765, 3.3935989e-05, -358.40428,  -2 ],
        [ 193.79558501, 0.4703899,  1.487e-05,  0.48837197, 2.9875551e-05, -390.08834,  -2 ],
        [ 211.80070675, 0.47064031, 1.302e-05,  0.48887521, 2.6160623e-05, -426.09859,  -2 ]
    ],

    impulse_table => [
        [
            -4.962698,  13.814839,  -4.9631602,    0.02982467, -0.0013594894, -0.14225902,
            -1.9401831, -0.2155842, -0.0020310088, -0.3728201, 0.33291368,    0.99809203,
            -2.8526329
        ],
        [
            -4.6051698, 12.742257,   -4.60596,     0.034726066, -0.0019437856, -0.13925305,
            -1.9371772, -0.39435313, -0.002903918, -0.36981413, 0.33291319,    0.99674207,
            -2.6978838
        ],
        [
            -4.1997047, 11.52587,    -4.2011566,   0.040956213, -0.0029156784, -0.13425306,
            -1.9321772, -0.59710092, -0.004355877, -0.36481413, 0.3329114,     0.99402816,
            -2.5290898
        ],
        [
            -3.9120226, 10.66284,   -3.9142589,    0.045760155, -0.0038875711, -0.12925308,
            -1.9271772, -0.7409712, -0.0058078358, -0.35981413, 0.33290796,    0.99082932,
            -2.4147808
        ],
        [
            -3.5065575, 9.4465067,   -3.510669,     0.052899933, -0.0058313546, -0.11925328,
            -1.9171773, -0.94382353, -0.0087117514, -0.34981412, 0.33289406,    0.98325425,
            -2.2636238
        ],
        [
            -2.9957319, 7.9143538,  -3.0046004,   0.061879021, -0.0097188791, -0.099256116,
            -1.8971786, -1.1998541, -0.014519534, -0.32981412, 0.33282259,    0.96451266,
            -2.095791
        ],
        [
            -2.3025847, 5.8377952,  -2.3278397,   0.071171613, -0.01943472, -0.049351462,
            -1.8472234, -1.5519373, -0.029035622, -0.27996523, 0.33218694,  0.9046402,
            -1.9322344
        ],
        [
            -1.8971196, 4.6290984,  -1.943834,   0.072773557, -0.029121182, 2.703517e-06,
            -1.7975285, -1.7695647, -0.04351831, -0.23041888, 0.33048633,   0.83582,
            -1.8934313
        ],
        [
            -1.6094375, 3.7805485,  -1.681748,   0.071588488, -0.038678285, 0.047661516,
            -1.7486526, -1.9421528, -0.05785254, -0.18223788, 0.32727359,   0.76542122,
            -1.9061883
        ],
        [
            -1.386294,  3.1337251,  -1.4876331,   0.069808845, -0.047867235, 0.091719224,
            -1.7016331, -2.1003168, -0.071752545, -0.13716096, 0.32223954,   0.69790882,
            -1.9472423
        ],
        [
            -1.2039724, 2.6179212,  -1.3371062,   0.068517644, -0.056290878, 0.13005005,
            -1.6579858, -2.2589276, -0.084697128, -0.09723682, 0.31526526,   0.63582802,
            -2.0053979
        ],
        [
            -1.0498217, 2.1947708,  -1.2168558,   0.068086017,  -0.063503894, 0.16135061,
            -1.6193745, -2.4255346, -0.095994449, -0.064925089, 0.30643614,   0.58033737,
            -2.0739321
        ],
        [
            -0.91629034, 1.8405707, -1.1187005,  0.068448138,  -0.069225576, 0.18569166,
            -1.586918,   -2.602482, -0.10505808, -0.041252944, 0.29601087,   0.53164647,
            -2.1483442
        ],
        [
            -0.7985073, 1.5393365,  -1.0372033,  0.069353468,  -0.07347155, 0.20413064,
            -1.560647,  -2.7879743, -0.11173039, -0.025827835, 0.28436064,  0.4893905,
            -2.2255157
        ],
        [
            -0.69314679, 1.2796676,  -0.9685592,  0.07053896,   -0.076490836, 0.21803372,
            -1.539681,   -2.9778054, -0.11632389, -0.016927299, 0.2719026,    0.45291732,
            -2.3033018
        ],
        [
            -0.59783661, 1.0531487,  -0.91001144, 0.071806708,  -0.078605381, 0.22863625,
            -1.5228397,  -3.1674206, -0.11936865, -0.012641221, 0.25904547,   0.42147885,
            -2.3802687
        ],
        [
            -0.51082523, 0.85343116, -0.85951366, 0.073035884,  -0.080095324, 0.2368885,
            -1.5090763,  -3.353234,  -0.12136138, -0.011494549, 0.24615402,   0.39434168,
            -2.455498
        ],
        [
            -0.43078252, 0.67565036, -0.81551861, 0.074165564,  -0.08116456, 0.2434645,
            -1.4976007,  -3.5329516, -0.12266858, -0.011808454, 0.23353095,  0.37084104,
            -2.528436
        ],
        [
            -0.35667455, 0.51603156, -0.77683891, 0.075172835,  -0.081949988, 0.24882638,
            -1.487853,   -3.7053538, -0.12353296, -0.013414023, 0.22141144,   0.35040023,
            -2.5987807
        ],
        [
            -0.28768168, 0.37161178, -0.74255145, 0.076055818,  -0.082541006, 0.25328937,
            -1.4794419,  -3.8699566, -0.12410896, -0.015931241, 0.20996566,   0.33253135,
            -2.666399
        ],
        [
            -0.22314316, 0.24003909, -0.71192979, 0.07682286,   -0.082995994, 0.25707082,
            -1.4720878,  -4.02673,   -0.12449411, -0.018690647, 0.19930553,   0.31682747,
            -2.7312685
        ],
        [
            -0.16251854, 0.11942812, -0.68439641, 0.077486369,  -0.083353584, 0.26032341,
            -1.4655895,  -4.1759027, -0.12475067, -0.021670016, 0.18949336,   0.30295165,
            -2.7934369
        ],
        [
            -0.10536012, 0.0082516612, -0.65948735, 0.078059643,  -0.083639885, 0.26315694,
            -1.4597968,  -4.3178401,   -0.124919,   -0.024602289, 0.18055095,   0.2906257,
            -2.852995
        ],
        [
            -0.051292902, -0.094739718, -0.63682608, 0.078555273,  -0.083872882, 0.26565199,
            -1.4545957,   -4.4529688,   -0.12502589, -0.027816838, 0.17246392,   0.27961988,
            -2.9100574
        ],
        [
            3.9198331e-07, -0.19057692, -0.61610458, 0.078984579,  -0.084065281, 0.26786924,
            -1.4498984,    -4.5817343,  -0.12508945, -0.031117444, 0.16517281,   0.26974415,
            -2.964751
        ],
        [
            0.050384782, -0.28301826, -0.59645732, 0.079369005,  -0.08423115, 0.26991807,
            -1.4454982,  -4.7085933,  -0.12512271, -0.034333427, 0.15838751,  0.26055766,
            -3.0189273
        ],
        [
            0.18232195, -0.51763868, -0.5481214,  0.080214573,  -0.084579273, 0.27475148,
            -1.434917,  -5.041608,   -0.12510954, -0.043544336, 0.14237348,   0.23875309,
            -3.1624297
        ],
        [
            0.26236466, -0.65512578, -0.52083251, 0.080626394,  -0.084743312, 0.27735658,
            -1.4291233, -5.2436271,  -0.1250525,  -0.049101088, 0.13382955,   0.22698019,
            -3.2503101
        ],
        [
            0.4054655,  -0.89276233, -0.47548795, 0.081204896,  -0.084973123, 0.28149402,
            -1.4198583, -5.6036799,  -0.12489608, -0.060354057, 0.12051535,   0.20833943,
            -3.4082825
        ],
        [
            0.53062864, -1.0929367, -0.43907894, 0.081576179,  -0.085125268, 0.28464135,
            -1.4128276, -5.9166728, -0.12472778, -0.069933513, 0.11065964,   0.19424605,
            -3.5468202
        ],
        [
            0.69314757, -1.3435638, -0.39576158, 0.081916926,  -0.085276023, 0.28817444,
            -1.4050642, -6.3195918, -0.12449312, -0.083310727, 0.099917547,  0.17853129,
            -3.7265483
        ],
        [
            0.99325216, -1.7834974, -0.32558507, 0.08226656,  -0.085466968, 0.29336647,
            -1.3942379, -7.0518296, -0.12407166, -0.10921959, 0.084737684,  0.15554992,
            -4.0563502
        ],
        [
            1.0986127, -1.9320777, -0.30348434, 0.082333425, -0.085516559, 0.2948528,
            -1.391353, -7.305173,  -0.12393637, -0.11788486, 0.080502722,  0.14894656,
            -4.1712325
        ],
        [
            1.2089607,  -2.0849085, -0.28154779, 0.082383119, -0.085561835, 0.29625254,
            -1.3887605, -7.5684752, -0.12380361, -0.12769681, 0.076544742,  0.14268847,
            -4.2909872
        ],
        [
            1.3862948,  -2.32521,   -0.24861337, 0.082431163, -0.085623169, 0.2982056,
            -1.3853877, -7.9874166, -0.12361037, -0.14389152, 0.071042637,  0.13383757,
            -4.4822044
        ],
        [
            1.6094383,  -2.6195648, -0.21070183, 0.0824556,   -0.085685413, 0.3002246,
            -1.3822722, -8.5076915, -0.1234022,  -0.16371248, 0.065325916,  0.12443566,
            -4.7207003
        ],
        [
            1.7917599,  -2.8544302, -0.18222867, 0.082457895, -0.085726697, 0.30157355,
            -1.3804484, -8.9275674, -0.1232585,  -0.18000819, 0.061436769,  0.11790666,
            -4.9139152
        ],
        [
            1.9459105,  -3.0496028, -0.15967808, 0.082453032, -0.085756485, 0.3025394,
            -1.3792929, -9.2792343, -0.12315377, -0.19366752, 0.058585952,  0.11304618,
            -5.0762121
        ],
        [
            2.0794419,  -3.2164365, -0.14115046, 0.082445891, -0.085779069, 0.30326514,
            -1.3785142, -9.5815922, -0.12307438, -0.2059823,  0.056386595,  0.10925036,
            -5.2160741
        ],
        [
            2.3027874,  -3.4913988, -0.11201814, 0.0824307,   -0.085812831, 0.30428409,
            -1.3777816, -10.083031, -0.12296212, -0.22653573, 0.053172754,  0.10248359,
            -5.4580995
        ],
        [
            2.995916,   -4.319002,  -0.033501683, 0.082383866, -0.08587761, 0.30632395,
            -1.3763405, -11.611104, -0.12273171,  -0.29002714, 0.045821657, 0.089946417,
            -6.1660669
        ],
        [
            3.9121588,  -5.3718052, 0.050619646, 0.082345822, -0.085920918, 0.30755591,
            -1.3760986, -13.584649, -0.122593,   -0.37416404, 0.039655883,  0.07870278,
            -7.0931862
        ],
        [
            4.6052203,  -6.147788,  0.10408342,  0.082331752, -0.085941688, 0.30797245,
            -1.3766306, -15.054581, -0.12255001, -0.43784777, 0.036432203,  0.072575829,
            -7.7904805
        ],
        [
            5.2984333,  -6.9118114, 0.15123108,  0.082324556, -0.085942739, 0.30817317,
            -1.3761139, -16.511603, -0.12252309, -0.50155172, 0.033934061,  0.067729371,
            -8.4860333
        ],
        [
            6.2146735,  -7.9080473, 0.20616054,  0.082320247, -0.085947007, 0.30829655,
            -1.3761172, -18.423065, -0.12250902, -0.58574681, 0.031349534,  0.062645606,
            -9.4038123
        ],
        [
            6.9077567,  -8.6538615, 0.24332204,  0.082318833, -0.085948427, 0.30833769,
            -1.3761188, -19.861029, -0.12250433, -0.64943124, 0.029764887,  0.059503607,
            -10.097453
        ],
        [
            7.6010029,  -9.3946265, 0.27744247,  0.082318138, -0.085949136, 0.30835827,
            -1.3761196, -21.29412,  -0.12250199, -0.71312784, 0.028408729,  0.056804601,
            -10.790995
        ],
        [
            8.5172381,  -10.367228, 0.31866825,  0.082317729, -0.085949561, 0.30837062,
            -1.3761201, -23.181854, -0.12250058, -0.79731022, 0.026880838,  0.053756628,
            -11.707419
        ],
        [
            9.2104555,  -11.099103, 0.34739771,  0.082317594, -0.085949708, 0.30837473,
            -1.3761203, -24.606217, -0.12250011, -0.86100079, 0.025880234,  0.051757974,
            -12.400703
        ],
        [
            10.819884,  -12.787912,  0.40749979,  0.082317486, -0.08594982, 0.30837802, -1.3761204, -27.903086,
            -0.1215948, -0.96632357, 0.023937869, 0.04787525,  -14.010187
        ],
        [
            11.513093,  -13.511682, 0.43101817,  0.082317471, -0.085949823, 0.30837843,
            -1.3761205, -29.319594, -0.12038764, -1.0011653,  0.023228325,  0.046456407,
            -14.703404
        ],
        [
            13.122513,  -15.185438, 0.48118269,  0.082317456, -0.085949855, 0.30837876,
            -1.3761205, -32.601876, -0.11677215, -1.0782228,  0.021800432,  0.043600817,
            -16.31283
        ],
        [
            13.815521,  -15.903735, 0.50113563,  0.082317452, -0.085949821, 0.3083788,
            -1.3761205, -34.012866, -0.11506674, -1.1099528,  0.021262716,  0.042525409,
            -17.005839
        ],
        [
            16.118143,  -18.282088, 0.56163436,  0.082317442, -0.085949897, 0.30837888,
            -1.3761204, -38.693008, -0.10937045, -1.2101334,  0.019728386,  0.039456771,
            -19.308462
        ],
        [
            18.420763,  -20.650308, 0.61490031,  0.082317433, -0.085949878, 0.30837908,
            -1.3761202, -43.363242, -0.10402521, -1.3035626,  0.01848693,   0.036973861,
            -21.611082
        ]
    ],

    tail_shock_table => [
        [ 7.046867,  -0.75244,    -0.041766853, -0.041766853 ],
        [ 7.6003429, -0.78707213, -0.055839367, -0.025878041 ],
        [ 8.516963,  -0.84157238, -0.060112633, -0.016839462 ],
        [ 9.2103139, -0.88078081, -0.061080159, -0.012813714 ],
        [ 10.819854, -0.96632103, -0.060797493, -0.0072803827 ],
        [ 11.513078, -1.0011641,  -0.060193871, -0.0057943111 ],
        [ 13.12251,  -1.0782225,  -0.058386088, -0.0034471643 ],
        [ 13.815519, -1.1099526,  -0.057533375, -0.0027513011 ],
        [ 16.118143, -1.2101333,  -0.054685228, -0.001224355 ]
    ],

    blast_info => [
        0.14925306, 1.9471776,  0.37986149, -0.29039171, 2.3512831,  0.0012567336, 3.4779416,   -0.19437851,
        0.50086761, 0.41006533, 0.67425467, -0.20311327, 0.25081663, 0.33291408,   0.041158721, -0.042974906,
        1148.5,     -0.75244,   16671.565,  -0.90794957
    ],

};
1;
