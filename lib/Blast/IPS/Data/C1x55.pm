package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='C', gamma=1.55
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'C1x55'} = {

    table_name  => 'C1x55',
    symmetry    => 1,
    gamma       => 1.55,
    data_source => 'C4000_G1x55/moc2_from_moc_r1000/',

    shock_table_info => [ 1, 1.55, 1.14e-06, 4000, 8.44e-07, 50, 5.47e-07, 9.15e-08, 5e-07, 81.242, -1.2002 ],

    shock_table => [
        [ -6.63557243, 12.00033657,  -1.99998507, -6.63693965, 0.99863186 ],
        [ -6.05492033, 10.83904772,  -1.9999523,  -6.05736512, 0.99755227 ],
        [ -5.2919887,  9.31326651,   -1.99979159, -5.29723883, 0.99473663 ],
        [ -4.79911741, 8.32769972,   -1.99943914, -4.80772556, 0.99135717 ],
        [ -4.38753292, 7.5048888,    -1.99872422, -4.40055065, 0.98690573 ],
        [ -4.04626266, 6.82297206,   -1.99748036, -4.06461807, 0.98149901 ],
        [ -3.75358534, 6.2386182,    -1.99548974, -3.77824652, 0.97508997 ],
        [ -3.49651298, 5.72598587,   -1.9924931,  -3.52849487, 0.96762649 ],
        [ -3.26608837, 5.26732657,   -1.98817861, -3.30648105, 0.95903118 ],
        [ -3.05648931, 4.85119191,   -1.98219136, -3.10645926, 0.94923129 ],
        [ -2.86234146, 4.46709385,   -1.97407759, -2.92321518, 0.93807774 ],
        [ -2.67983801, 4.10774744,   -1.96329773, -2.75313421, 0.92540289 ],
        [ -2.50539605, 3.76643707,   -1.94915673, -2.59292955, 0.9109561 ],
        [ -2.33498485, 3.43578112,   -1.9306633,  -2.4390725,  0.89431696 ],
        [ -2.16330031, 3.10632255,   -1.90622333, -2.28717111, 0.87472889 ],
        [ -1.97945327, 2.75885152,   -1.87235717, -2.12854093, 0.85034198 ],
        [ -1.79426265, 2.41595877,   -1.82921222, -1.97362563, 0.8220747 ],
        [ -1.62586682, 2.11180863,   -1.78183699, -1.83758016, 0.79322696 ],
        [ -1.46858562, 1.83547822,   -1.73101208, -1.71510265, 0.76384551 ],
        [ -1.31911524, 1.58067117,   -1.67774533, -1.60313223, 0.73413065 ],
        [ -1.17575761, 1.34402497,   -1.62329555, -1.50000587, 0.70444568 ],
        [ -1.03333774, 1.11681066,   -1.56726118, -1.40182054, 0.67430405 ],
        [ -0.89165843, 0.89875404,   -1.510916,   -1.30842068, 0.64418387 ],
        [ -0.74923855, 0.68757319,   -1.45488716, -1.21881385, 0.61426117 ],
        [ -0.60469465, 0.48129276,   -1.39969305, -1.13217451, 0.5847051 ],
        [ -0.45663448, 0.27808487,   -1.34573376, -1.04777012, 0.55566777 ],
        [ -0.30405131, 0.07678264,   -1.29345588, -0.96516629, 0.52736569 ],
        [ -0.14543969, -0.12432331,  -1.24306094, -0.88372458, 0.49991404 ],
        [ 0.02037762,  -0.3263818,   -1.1948169,  -0.80305375, 0.47348603 ],
        [ 0.19459564,  -0.53047402,  -1.14894044, -0.72280019, 0.44824151 ],
        [ 0.37866662,  -0.73789,     -1.1055513,  -0.6425378,  0.42429424 ],
        [ 0.57429466,  -0.95008759,  -1.06470913, -0.56178795, 0.40172599 ],
        [ 0.78339116,  -1.16862305,  -1.02644835, -0.48004809, 0.38060216 ],
        [ 1.0078072,   -1.39487981,  -0.99083188, -0.39689017, 0.36099879 ],
        [ 1.24389346,  -1.62490069,  -0.95861425, -0.31380101, 0.34336265 ],
        [ 1.48930238,  -1.85654647,  -0.92997804, -0.23149854, 0.32780951 ],
        [ 1.74558766,  -2.09153502,  -0.90451898, -0.14929062, 0.31412358 ],
        [ 2.01443285,  -2.33158358,  -0.88189067, -0.06650322, 0.30211599 ],
        [ 2.29799776,  -2.57872244,  -0.86177783, 0.01762947,  0.29161248 ],
        [ 2.59750166,  -2.83408361,  -0.84398487, 0.10355806,  0.28250044 ],
        [ 2.91518596,  -3.09963262,  -0.82828603, 0.19201257,  0.27465039 ],
        [ 3.25213277,  -3.37633003,  -0.81454339, 0.28338795,  0.26797634 ],
        [ 3.61012459,  -3.66571836,  -0.80259715, 0.37827853,  0.26237965 ],
        [ 3.99125222,  -3.96957647,  -0.79229571, 0.47736024,  0.25776501 ],
        [ 4.39753496,  -4.28961861,  -0.7835048,  0.58129294,  0.2540441 ],
        [ 4.83100572,  -4.6275744,   -0.77609705, 0.69074769,  0.25113043 ],
        [ 5.29497875,  -4.9861704,   -0.76993299, 0.8067236,   0.24893209 ],
        [ 5.79114803,  -5.36688097,  -0.76490129, 0.92981904,  0.24736678 ],
        [ 6.32480258,  -5.77394059,  -0.76085806, 1.06152832,  0.24634134 ],
        [ 6.90962108,  -6.21791055,  -0.75764007, 1.20540206,  0.24576577 ],
        [ 7.522239,    -6.68128098,  -0.75525752, 1.35588668,  0.24557375 ],
        [ 8.20273119,  -7.194571,    -0.75344896, 1.52301578,  0.24566794 ],
        [ 8.95252733,  -7.75898065,  -0.75215086, 1.7073248,   0.24598189 ],
        [ 9.78905764,  -8.38777239,  -0.75125874, 1.91328694,  0.24644953 ],
        [ 10.7401453,  -9.1019813,   -0.75067747, 2.14795097,  0.24701365 ],
        [ 11.84780762, -9.9332603,   -0.75032529, 2.4219041,   0.24762429 ],
        [ 13.18185788, -10.93408469, -0.75013249, 2.75267276,  0.2482382 ],
        [ 14.87670733, -12.20535491, -0.75004156, 3.1739211,   0.24881946 ],
        [ 17.22996845, -13.970349,   -0.7500081,  3.76012035,  0.24933513 ],
        [ 21.16826735, -16.92408397, -0.7500005,  4.74302276,  0.24974988 ]
    ],

    energy_table => [
        [ -6.63557243, 0.01377392, 0.00976264, 0.01377392, 0.00976264,    6.7287564,   -1.182424 ],
        [ -6.05492033, 0.02078179, 0.01470908, 0.02078179, 0.01470908,    6.0364816,   -1.204668 ],
        [ -5.2919887,  0.03562624, 0.0251064,  0.03562624, 0.0251064,     5.1049409,   -1.2420372 ],
        [ -4.79911741, 0.05037426, 0.03528682, 0.05037426, 0.03528682,    4.4860785,   -1.2738937 ],
        [ -4.38753292, 0.06713082, 0.04661619, 0.06713082, 0.04661619,    3.9554852,   -1.3082437 ],
        [ -4.04626266, 0.08497608, 0.05833873, 0.08497608, 0.05833873,    3.503617,    -1.3433008 ],
        [ -3.75358534, 0.1037469,  0.07021336, 0.1037469,  0.07021336,    3.1056381,   -1.3795286 ],
        [ -3.49651298, 0.12328559, 0.08199514, 0.12328559, 0.08199514,    2.7465419,   -1.4174023 ],
        [ -3.26608837, 0.14348674, 0.09346545, 0.14348674, 0.09346545,    2.4156964,   -1.4567812 ],
        [ -3.05648931, 0.16421629, 0.1043858,  0.16421629, 0.1043858,     2.1063578,   -1.4980902 ],
        [ -2.86234146, 0.18547435, 0.11458274, 0.18547435, 0.11458274,    1.8115088,   -1.5401499 ],
        [ -2.67983801, 0.20724009, 0.12385452, 0.20724009, 0.12385452,    1.5267433,   -1.5834127 ],
        [ -2.50539605, 0.22957014, 0.13200891, 0.22957014, 0.13200891,    1.2466809,   -1.6245791 ],
        [ -2.33498485, 0.25266824, 0.13885374, 0.25266824, 0.13885374,    0.96665425,  -1.6613666 ],
        [ -2.16330031, 0.27698829, 0.14414949, 0.27698829, 0.14414949,    0.67828809,  -1.6890268 ],
        [ -1.97945327, 0.30383946, 0.14751547, 0.30383946, 0.14751547,    0.36591517,  -1.6941904 ],
        [ -1.79426265, 0.33125398, 0.14804574, 0.33125398, 0.14804574,    0.05308126,  -1.6622875 ],
        [ -1.62586682, 0.35603404, 0.14582882, 0.35603404, 0.14582882,    -0.22271187, -1.5912984 ],
        [ -1.46858562, 0.37865509, 0.14146643, 0.37865509, 0.14146643,    -0.46616694, -1.4989639 ],
        [ -1.31911524, 0.39937246, 0.13546621, 0.39937246, 0.13546621,    -0.68326685, -1.4090895 ],
        [ -1.17575761, 0.41829381, 0.12830711, 0.41829381, 0.12830711,    -0.8793072,  -1.3386328 ],
        [ -1.03333774, 0.4359986,  0.120184,   0.4359986,  0.120184,      -1.0658719,  -1.2954417 ],
        [ -0.89165843, 0.45241437, 0.11147321, 0.45241437, 0.11147321,    -1.2473626,  -1.2776738 ],
        [ -0.74923855, 0.46764858, 0.10244286, 0.46764858, 0.10244286,    -1.4288531,  -1.2784042 ],
        [ -0.60469465, 0.48179461, 0.09332187, 0.48179461, 0.09332187,    -1.6142357,  -1.2909454 ],
        [ -0.45663448, 0.49493852, 0.08429957, 0.49493852, 0.08429957,    -1.8066489,  -1.3100121 ],
        [ -0.30405131, 0.50712581, 0.07555352, 0.50712581, 0.07555352,    -2.0081781,  -1.3319125 ],
        [ -0.14543969, 0.51843643, 0.0672001,  0.51843643, 0.0672001,     -2.2212699,  -1.3544118 ],
        [ 0.02037762,  0.52891559, 0.05934607, 0.52891559, 0.05934607,    -2.4477496,  -1.3760502 ],
        [ 0.19459564,  0.53860653, 0.05206826, 0.53860653, 0.05206826,    -2.6893515,  -1.396003 ],
        [ 0.37866662,  0.54756219, 0.04540726, 0.54756219, 0.04540726,    -2.948109,   -1.4138821 ],
        [ 0.57429466,  0.55583865, 0.03937599, 0.55583865, 0.03937599,    -3.2263951,  -1.429567 ],
        [ 0.78339116,  0.5634894,  0.0339681,  0.5634894,  0.0339681,     -3.526886,   -1.443108 ],
        [ 1.0078072,   0.57055615, 0.02916824, 0.57055615, 0.02916824,    -3.8521916,  -1.4546127 ],
        [ 1.24389346,  0.57693805, 0.02503653, 0.57693805, 0.02503653,    -4.1968594,  -1.4640419 ],
        [ 1.48930238,  0.58263844, 0.02154018, 0.58263844, 0.02154018,    -4.5571991,  -1.471621 ],
        [ 1.74558766,  0.58776559, 0.01857423, 0.58776559, 0.01857423,    -4.9352364,  -1.4776961 ],
        [ 2.01443285,  0.59240806, 0.01605021, 0.59240806, 0.01605021,    -5.33325,    -1.4825585 ],
        [ 2.29799776,  0.59664264, 0.01389214, 0.59664264, 0.01389214,    -5.7542803,  -1.4864483 ],
        [ 2.59750166,  0.60051705, 0.01204457, 0.60051705, 0.01204457,    -6.2000065,  -1.4895544 ],
        [ 2.91518596,  0.60408217, 0.01045514, 0.60408217, 0.01045514,    -6.6736659,  -1.4920298 ],
        [ 3.25213277,  0.60736639, 0.00908609, 0.60736639, 0.00908609,    -7.1767774,  -1.493971 ],
        [ 3.61012459,  0.61040001, 0.00790226, 0.61040001, 0.00790226,    -7.7119198,  -1.4954903 ],
        [ 3.99125222,  0.6132092,  0.00687378, 0.6132092,  0.00687378,    -8.2821543,  -1.4966732 ],
        [ 4.39753496,  0.61581356, 0.00597637, 0.61581356, 0.00597637,    -8.8904415,  -1.4975879 ],
        [ 4.83100572,  0.61822817, 0.00518998, 0.61822817, 0.00518998,    -9.5397806,  -1.4982896 ],
        [ 5.29497875,  0.62047013, 0.00449641, 0.62052814, 0.0049050881,  -10.235091,  -1.498802 ],
        [ 5.79114803,  0.62254425, 0.00388342, 0.62296624, 0.0048357959,  -10.978857,  -1.4991794 ],
        [ 6.32480258,  0.62446639, 0.00333737, 0.6254027,  0.0042825973,  -11.778994,  -1.4994595 ],
        [ 6.90962108,  0.62626874, 0.0028422,  0.62777615, 0.0038139333,  -12.655974,  -1.4996528 ],
        [ 7.522239,    0.62787431, 0.00241285, 0.62997982, 0.0033745363,  -13.574737,  -1.4997843 ],
        [ 8.20273119,  0.62937789, 0.00201916, 0.63212248, 0.002933595,   -14.595367,  -1.4998768 ],
        [ 8.95252733,  0.63075431, 0.00166442, 0.63415413, 0.0024993571,  -15.719999,  -1.4999396 ],
        [ 9.78905764,  0.63200819, 0.00134543, 0.63606716, 0.0020872169,  -16.974765,  -1.4999845 ],
        [ 10.7401453,  0.63314531, 0.00105454, 0.63785918, 0.0016892963,  -18.401402,  -1.5000223 ],
        [ 11.84780762, 0.63416804, 0.00080083, 0.63951781, 0.0013178047,  -20.062942,  -1.5000665 ],
        [ 13.18185788, 0.63507606, 0.00057322, 0.64103704, 0.00098066955, -22.064146,  -1.5001451 ],
        [ 14.87670733, 0.63586782, 0.00037509, 0.64244229, 0.00067986227, -24.606771,  -1.5003923 ],
        [ 17.22996845, 0.63653503, 0.00020824, 0.64365161, 0.00037743293, -28.138192,  -1.5005652 ],
        [ 21.16826735, 0.63705679, 7.78e-05,   0.64459728, 0.00014100241, -34.047307,  -1.5 ]
    ],

    impulse_table => [
        [
            -7.5429392, 13.815064,  -7.5434908,    0.029436389, 0,          -0.19742984,
            0,          -2.3041717, -0.0034405089, -0.44910601, 0.24481573, 0.98878141,
            -3.29739
        ],
        [
            -6.9077553, 12.544699,  -6.9087965,    0.036429327, 0,          -0.19695968,
            0,          -2.3041768, -0.0047266247, -0.44863585, 0.24481515, 0.98239701,
            -3.079283
        ],
        [
            -6.5022902, 11.733774, -6.5038522,    0.041480974, 0,          -0.19645968,
            0,          -2.304184, -0.0057889094, -0.44813585, 0.24481414, 0.97653568,
            -2.9419266
        ],
        [
            -6.2146081, 11.158417, -6.2166913,    0.045329209, 0,          -0.19595968,
            0,          -2.304194, -0.0066844568, -0.44763585, 0.24481279, 0.97123177,
            -2.8456555
        ],
        [
            -5.809143, 10.347507,  -5.8122676,    0.051082574, 0,          -0.19495968,
            0,         -2.3042223, -0.0081867542, -0.44663585, 0.24480903, 0.96167111,
            -2.712113
        ],
        [
            -5.2983174, 9.3259225,  -5.3035343,   0.058736523, 0,        -0.19295968,
            0,          -2.3043127, -0.010569054, -0.44463585, 0.244797, 0.94502525,
            -2.5486021
        ],
        [
            -4.6051702, 7.9399376,   -4.6156289, 0.069237725, 0, -0.18795968, 0, -2.304737,
            -0.0149469, -0.43963585, 0.2447406,  0.91057805,  -2.3391526
        ],
        [
            -3.912023, 6.55488,   -3.9330387,   0.07841468,  0,          -0.17795971,
            0,         -2.306435, -0.021138108, -0.42963585, 0.24451543, 0.85560535,
            -2.1525357
        ],
        [
            -3.5065579, 5.746001,   -3.5382163,   0.082161012, 0,          -0.16795995,
            0,          -2.3092684, -0.025888783, -0.41963585, 0.24414178, 0.81016585,
            -2.059866
        ],
        [
            -2.9957323, 4.7308259,  -3.0488845,   0.083864404, 0,          -0.14796454,
            0,          -2.3183633, -0.033422116, -0.39963585, 0.24295959, 0.73542698,
            -1.9692821
        ],
        [
            -2.65926, 4.0673615,  -2.7341075,   0.082518766, 0,          -0.12799205,
            0,        -2.3320826, -0.039544459, -0.37978273, 0.24122371, 0.6746759,
            -1.9316038
        ],
        [
            -2.3025851, 3.3732943,  -2.4101528,   0.078778252, 0,          -0.098200762,
            0,          -2.3615226, -0.047254481, -0.34992965, 0.23767484, 0.60103603,
            -1.9172102
        ],
        [
            -2.0402208, 2.8730036,  -2.1804733,   0.074950617, 0,          -0.069007142,
            0,          -2.4018914, -0.053833287, -0.32051776, 0.23314853, 0.54231094,
            -1.9272113
        ],
        [
            -1.89712, 2.6054235,  -2.0590229,   0.0728342,   0,          -0.050270623,
            0,        -2.4350079, -0.057749515, -0.30140116, 0.22968317, 0.50936692,
            -1.9410438
        ],
        [
            -1.6094379, 2.0825762, -1.8245727,   0.069733291, 0,         -0.0085520603,
            0,          -2.539294, -0.066077982, -0.25822389, 0.2198688, 0.44301064,
            -1.9881418
        ],
        [
            -1.3862944,  1.6942073,   -1.6529071, 0.069453126, 0, 0.022216566, 0, -2.6705859,
            -0.07216759, -0.22536415, 0.20905199, 0.39314094,  -2.0427724
        ],
        [
            -1.2039728, 1.3899809,  -1.5199655,   0.070829236, 0,          0.042110314,
            0,          -2.8190086, -0.075937918, -0.20575793, 0.19788393, 0.35454099,
            -2.0986924
        ],
        [
            -1.0498221,  1.1427,     -1.412965,  0.072854629, 0, 0.054607612, 0, -2.9728318,
            -0.07789083, -0.1961082, 0.18683195, 0.3239174,   -2.1532263
        ],
        [
            -0.91629074, 0.93609201, -1.3243528,   0.074991426, 0,          0.062770175,
            0,           -3.1237279, -0.078728357, -0.19247207, 0.17620584, 0.29910357,
            -2.2052631
        ],
        [
            -0.7985077, 0.75972792, -1.2493311,   0.077013228, 0,          0.068420919,
            0,          -3.2675502, -0.078961344, -0.19210621, 0.16619439, 0.27862787,
            -2.2544212
        ],
        [
            -0.69314719, 0.60657587, -1.1846849,   0.078844152, 0,          0.072555568,
            0,           -3.4028186, -0.078882459, -0.19359351, 0.15689876, 0.26146303,
            -2.3006632
        ],
        [
            -0.59783701, 0.47170296, -1.1281695,   0.080472518, 0,          0.075725143,
            0,           -3.5293799, -0.078647408, -0.19570465, 0.14835911, 0.2468744,
            -2.3441097
        ],
        [
            -0.51082563, 0.35153583, -1.078165,    0.081911494, 0,          0.07824771,
            0,           -3.6476599, -0.078338235, -0.19847776, 0.14057458, 0.23432512,
            -2.3849471
        ],
        [
            -0.35667495, 0.14531083, -0.99316884,  0.084306611, 0,          0.082053595,
            0,           -3.8619834, -0.077651565, -0.20422109, 0.12714519, 0.21384095,
            -2.45962
        ],
        [
            -0.22314356, -0.026799928, -0.92307898,  0.086195466, 0,          0.084834346,
            0,           -4.051087,    -0.076979874, -0.21026906, 0.11623895, 0.19781597,
            -2.5262612
        ],
        [
            -0.10536052, -0.17390126, -0.86382122,  0.087710189, 0,          0.086988275,
            0,           -4.2195838,  -0.076365731, -0.21595985, 0.10740919, 0.18491762,
            -2.5862211
        ],
        [
            -8.1811448e-09, -0.30197647, -0.81273401,  0.088946077, 0,          0.088725352,
            0,              -4.3711428,  -0.075816245, -0.22171528, 0.10020698, 0.17429384,
            -2.6405937
        ],
        [
            0.080838946, -0.39812122, -0.77470134,  0.089820724, 0,           0.089955926,
            0,           -4.4877396,  -0.075400969, -0.22607831, 0.095109725, 0.16668948,
            -2.6826888
        ],
        [
            0.18232155, -0.51635298, -0.72831233,  0.090832688, 0,           0.09139233,
            0,          -4.6342864,  -0.074892614, -0.23192938, 0.089202332, 0.15777062,
            -2.7359013
        ],
        [
            0.26236426,  -0.60777113, -0.69273507, 0.09156696, 0, 0.0924498, 0, -4.7498879,
            -0.07450444, -0.23676633, 0.084901543, 0.15119716, -2.7781004
        ],
        [
            0.33647223, -0.69104242, -0.66055086,  0.092199569, 0,           0.093375787,
            0,          -4.8568524,  -0.074156499, -0.24151826, 0.081180505, 0.14545018,
            -2.8173103
        ],
        [
            0.4054651, -0.76743781, -0.63121116,  0.092750024, 0,           0.094195824,
            0,         -4.9563283,  -0.073843279, -0.24631194, 0.077926333, 0.14037584,
            -2.8539074
        ],
        [
            0.53062824, -0.90340724, -0.57943395,  0.093660764, 0,           0.095589615,
            0,          -5.1364142,  -0.073302549, -0.25474653, 0.072496473, 0.1318013,
            -2.9204622
        ],
        [
            0.69314717, -1.0752835, -0.51478615,  0.094691223, 0,           0.097239001,
            0,          -5.3692882, -0.072655194, -0.26715098, 0.066255507, 0.12176664,
            -3.0070478
        ],
        [
            0.99325176, -1.3804422, -0.40215325,  0.096217876, 0,           0.099883478,
            0,          -5.7957169, -0.071620487, -0.29288145, 0.056705505, 0.1060005,
            -3.1669051
        ],
        [
            1.0986123, -1.4842575, -0.36443606,  0.096657836, 0,           0.10070592,
            0,         -5.9441945, -0.071304313, -0.30265871, 0.053853187, 0.10118666,
            -3.2229093
        ],
        [
            1.2089603, -1.5913358, -0.32583801,  0.097073735, 0,          0.10151533,
            0,         -6.0989835, -0.070997671, -0.31315328, 0.05109964, 0.096490592,
            -3.2814616
        ],
        [
            1.3862944, -1.7601676, -0.26558148,  0.097657713, 0,           0.10271431,
            0,         -6.3461891, -0.070554151, -0.33261651, 0.047114991, 0.089606331,
            -3.3752989
        ],
        [
            1.6094379, -1.9675128, -0.19252548,  0.098267838, 0,           0.10406322,
            0,         -6.6545783, -0.070074376, -0.35923229, 0.042755952, 0.081949767,
            -3.4928723
        ],
        [
            1.7917595, -2.1332019, -0.13483843,  0.098681463, 0,           0.10504864,
            0,         -6.9044116, -0.069740105, -0.38221468, 0.039645312, 0.076402384,
            -3.5884982
        ],
        [
            1.9459101, -2.2709702, -0.087301935, 0.098981698, 0,           0.10580862,
            0,         -7.1142116, -0.069492925, -0.40487674, 0.037281782, 0.072139297,
            -3.66904
        ],
        [
            2.0794415, -2.3887532, -0.046947684, 0.099210456, 0,           0.10641819,
            0,         -7.2949394, -0.069302508, -0.4253325,  0.035406075, 0.068725967,
            -3.7385852
        ],
        [
            2.3025851, -2.582675,  0.018966839,  0.099537334, 0,           0.10734357,
            0,         -7.5949837, -0.069027392, -0.46358743, 0.032579161, 0.063530575,
            -3.8543593
        ],
        [
            2.7080502, -2.9270612, 0.13462465,   0.099991624, 0,           0.10876977,
            0,         -8.1344785, -0.068641696, -0.54164559, 0.028245286, 0.055445311,
            -4.0634307
        ],
        [
            2.9957323, -3.1662048, 0.21406422,   0.10023119,  0,           0.10961223,
            0,         -8.5133474, -0.068438598, -0.60969754, 0.025666228, 0.050565066,
            -4.2108874
        ],
        [
            3.4011974, -3.4973514, 0.32314465,   0.10048487,  0,           0.1106047,
            0,         -9.0426925, -0.068224615, -0.72097659, 0.022567205, 0.044635019,
            -4.4176853
        ],
        [
            3.9122247, -3.9068867, 0.45695623,   0.1007065,   0, 0.11159753,
            0,         -9.7035066, -0.068055587, -0.90570727, 0, 0.038037392,
            -4.6853565
        ],
        [
            4.6052133, -4.4519365, 0.63389129,   0.10089331, 0,           0.11258777,
            0,         -10.591034, -0.067904426, -1.2349338, 0.015872534, 0.031425895,
            -5.0317345
        ],
        [
            5.2983884, -4.9887956, 0.80757236,   0.10100164, 0,           0.11328768,
            0,         -11.471754, -0.067748706, -1.6876623, 0.013132148, 0.026107944,
            -5.3783443
        ],
        [
            6.214653, -5.6900925, 1.0343851,    0.10107914, 0,           0.1139198,
            0,        -12.628837, -0.067759143, -2.6227614, 0.010307139, 0.020552775,
            -5.8365486
        ],
        [
            6.9080198, -6.2166974, 1.2050085,    0.10111037, 0,            0.11423204,
            0,         -13.50098,  -0.066194872, -3.5119579, 0.0086148899, 0.017199278,
            -6.1832757
        ],
        [
            7.6010806, -6.7408169, 1.375248,     0.10112853, 0,            0.11445205,
            0,         -14.370848, -0.065695145, -4.8975307, 0.0072163604, 0.014417543,
            -6.5298362
        ],
        [
            8.5172258,   -7.431427,  1.6002941,    0.10114107, 0, 0.11464696, 0, -15.518872,
            -0.06120355, -7.1401059, 0.0057214932, 0.01143694, -6.9879367
        ],
        [
            9.2103751, -7.9528784, 1.7707682,   0.10114548, 0,            0.11474509,
            0,         -16.386587, -0.38881623, -10.216024, 0.0048045563, 0.0096060964,
            -7.3345276
        ],
        [
            10.819951, -9.1618882, 2.1676659,    0.10114639, 0,            0.1148759,
            0,         -18.400014, -0.060958705, -21.381422, 0.0032078838, 0.0064151702,
            -8.1393535
        ],
        [
            11.512992, -9.6820267, 2.3390249,    0.10114379, 0,            0.11490685,
            0,         -19.266597, -0.057482167, -29.914935, 0.0026967605, 0.0053932269,
            -8.4858964
        ],
        [
            13.122456, -10.889525, 2.7379277,    0.10112969, 0,            0.1149481,
            0,         -21.278728, -0.012218097, -14.1034,   0.0018028401, 0.0036056358,
            -9.2907183
        ],
        [
            13.815678, -11.409518, 2.9100892,    0.10111854, 0,            0.11495786,
            0,         -22.145324, -0.010343247, -16.905949, 0.0015158595, 0.0030317112,
            -9.6373953
        ],
        [
            16.11814, -13.136465, 3.4830125,     0.10103516, 0,             0.1149739,
            0,        -25.023439, -0.0059298431, -30.736588, 0.00085239301, 0.0017048538,
            -10.789114
        ],
        [
            18.420799, -14.863479, 4.0571391,    0.10077084, 0,             0.11497972,
            0,         -27.901608, -0.067688724, -1113.8633, 0.00047930863, 0.00095884071,
            -11.941983
        ]
    ],

    tail_shock_table => [
        [ 4.4120974, -1.2002,    -0.034278058,  -0.034278058 ],
        [ 4.5861849, -1.2715292, -0.042355925,  -0.023965451 ],
        [ 5.2871136, -1.5909012, -0.042929451,  -0.012175085 ],
        [ 6.2090107, -2.1046255, -0.037377142,  -0.0065564935 ],
        [ 6.9046783, -2.5812995, -0.032674518,  -0.0044943493 ],
        [ 7.599101,  -3.1512015, -0.028255891,  -0.0032082973 ],
        [ 8.5162345, -4.0793934, -0.023102674,  -0.0021233301 ],
        [ 9.2097874, -4.9430719, -0.019752744,  -0.0015758551 ],
        [ 10.819776, -7.6538548, -0.013605375,  -0.0008107597 ],
        [ 11.512888, -9.2147337, -0.011559911,  -0.00061359698 ],
        [ 13.122425, -14.107988, -0.0078854757, -0.00032400111 ],
        [ 13.815659, -16.901805, -0.0066715834, -0.00024854072 ],
        [ 16.118137, -30.736476, -0.0038257064, -0.00010310469 ]
    ],

    blast_info => [
        0.19795968, 0,          0.44963585, -0.149469,   0,          0,          1.2302197,   -0.16175777,
        0.62326482, 0.16723098, 0.72264503, -0.11624147, 0.34188816, 0.24481572, 0.065183974, 0,
        81.242,     -1.2002,    129.51508,  -1.3895564
    ],

};
1;