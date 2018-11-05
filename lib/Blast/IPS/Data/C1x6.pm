package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='C', gamma=1.6
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'C1x6'} = {

    table_name  => 'C1x6',
    symmetry    => 1,
    gamma       => 1.6,
    data_source => 'C8000_G1x6/moc2_from_moc_r5000/',

    shock_table_info => [ 1, 1.6, 7.63e-07, 8000, 2.16e-07, 50, 1.33e-07, 1.3e-07, 5e-07, 84.682, -1.239 ],

    shock_table => [
        [ -6.59721775, 12.00017708,  -1.99998488, -6.59859354, 0.99862328 ],
        [ -5.99411966, 10.79399743,  -1.99994949, -5.9966357,  0.99748085 ],
        [ -5.82022227, 10.44621259,  -1.9999323,  -5.82321689, 0.99700099 ],
        [ -5.18151608, 9.16888781,   -1.99975682, -5.18719534, 0.9943053 ],
        [ -4.68221894, 8.17050163,   -1.99934135, -4.69159214, 0.99058592 ],
        [ -4.28448325, 7.37543038,   -1.99854254, -4.29846404, 0.98593158 ],
        [ -3.95494042, 6.71702405,   -1.99718808, -3.97442497, 0.98035286 ],
        [ -3.66917986, 6.14658507,   -1.99503682, -3.69517918, 0.9737271 ],
        [ -3.41757173, 5.64498782,   -1.99183139, -3.45110684, 0.96603996 ],
        [ -3.19138195, 5.19493635,   -1.98724983, -3.23355891, 0.95720486 ],
        [ -2.98467584, 4.78477179,   -1.98091263, -3.03670406, 0.94712289 ],
        [ -2.79286243, 4.40557748,   -1.97236313, -2.85609675, 0.93566212 ],
        [ -2.61213456, 4.05008426,   -1.96104022, -2.68813779, 0.92264236 ],
        [ -2.43879704, 3.71138417,   -1.94620518, -2.52946092, 0.90778537 ],
        [ -2.26865664, 3.38183189,   -1.9267853,  -2.37643146, 0.89062432 ],
        [ -2.09610295, 3.05148448,   -1.90102501, -2.22445912, 0.87031679 ],
        [ -1.90820867, 2.69755605,   -1.86476944, -2.06328573, 0.84462112 ],
        [ -1.72669697, 2.36293715,   -1.82073787, -1.91250667, 0.81615655 ],
        [ -1.56072874, 2.06464647,   -1.77258597, -1.77942218, 0.78712117 ],
        [ -1.40530725, 1.79306931,   -1.72119208, -1.6593526,  0.75762349 ],
        [ -1.25763494, 1.54279849,   -1.66769411, -1.54964757, 0.72793635 ],
        [ -1.11495258, 1.30872849,   -1.61290259, -1.44789703, 0.69818268 ],
        [ -0.97285613, 1.08352616,   -1.55663819, -1.3508281,  0.66801031 ],
        [ -0.83130161, 0.86717338,   -1.50021593, -1.25840022, 0.63792425 ],
        [ -0.68866819, 0.65720455,   -1.44418537, -1.16954849, 0.60806488 ],
        [ -0.54366038, 0.45180939,   -1.3890831,  -1.08352291, 0.57861605 ],
        [ -0.39493137, 0.24924938,   -1.33531288, -0.99963216, 0.54973578 ],
        [ -0.24134056, 0.04820244,   -1.28326437, -0.91738111, 0.52161256 ],
        [ -0.08160699, -0.15272368,  -1.23320725, -0.83626395, 0.49440017 ],
        [ 0.08554899,  -0.35479868,  -1.18536454, -0.75584109, 0.46824824 ],
        [ 0.26146194,  -0.55924976,  -1.13990992, -0.67570326, 0.44329281 ],
        [ 0.44767659,  -0.76743839,  -1.09694652, -0.59540093, 0.41963763 ],
        [ 0.64573868,  -0.98061781,  -1.05657451, -0.51453763, 0.39738589 ],
        [ 0.85763865,  -1.20041186,  -1.01881155, -0.43258662, 0.37659262 ],
        [ 1.0854073,   -1.42836578,  -0.98369372, -0.34906222, 0.35731813 ],
        [ 1.32283905,  -1.65809918,  -0.95227306, -0.26631445, 0.34016895 ],
        [ 1.5698216,   -1.88975438,  -0.92435287, -0.18421812, 0.32505094 ],
        [ 1.82826335,  -2.12534728,  -0.89950771, -0.10198227, 0.31173736 ],
        [ 2.09982101,  -2.36653126,  -0.8774172,  -0.01896263, 0.30005429 ],
        [ 2.38623051,  -2.61494211,  -0.85781133, 0.06546768,  0.28985222 ],
        [ 2.68901128,  -2.87196712,  -0.8404777,  0.15184524,  0.28101013 ],
        [ 3.00995311,  -3.13918565,  -0.8252191,  0.24077015,  0.27341358 ],
        [ 3.35080543,  -3.41811304,  -0.81186816, 0.33282197,  0.2669622 ],
        [ 3.71316011,  -3.71012516,  -0.80028052, 0.42853869,  0.26156559 ],
        [ 4.09948777,  -4.01729733,  -0.7902998,  0.52869261,  0.2571269 ],
        [ 4.511045,    -4.34073708,  -0.78181325, 0.63374621,  0.25356759 ],
        [ 4.95041288,  -4.68260942,  -0.7746811,  0.74451282,  0.2507958 ],
        [ 5.42030073,  -5.04517345,  -0.76877308, 0.86184058,  0.2487231 ],
        [ 5.92400454,  -5.43113758,  -0.76395977, 0.98672703,  0.24726116 ],
        [ 6.46560698,  -5.84380474,  -0.76011317, 1.12036528,  0.24632252 ],
        [ 7.04957736,  -6.28676223,  -0.7571103,  1.26404307,  0.24582199 ],
        [ 7.68357462,  -6.76599802,  -0.75482282, 1.41982997,  0.24567597 ],
        [ 8.37684531,  -7.28866853,  -0.75313333, 1.59018221,  0.24580709 ],
        [ 9.14302655,  -7.86520839,  -0.75193002, 1.77863541,  0.24614421 ],
        [ 10.00154701, -8.51037294,  -0.75111105, 1.99015714,  0.24662352 ],
        [ 10.98142848, -9.24608849,  -0.7505852,  2.2320976,   0.24718878 ],
        [ 12.12928747, -10.10745296, -0.75027265, 2.51619001,  0.24779168 ],
        [ 13.52553678, -11.15488656, -0.75010602, 2.86260535,  0.24839119 ],
        [ 15.32598756, -12.50533431, -0.75003074, 3.31035922,  0.2489521 ],
        [ 17.89685356, -14.43352058, -0.75000509, 3.95107495,  0.249442 ],
        [ 22.51756936, -17.89906433, -0.75000019, 5.10472023,  0.24982321 ]
    ],

    energy_table => [
        [ -6.59721775, 0.01083281, 0.00811266, 0.01083281, 0.00811266,    6.7097422,    -1.1815096 ],
        [ -5.99411966, 0.0170119,  0.01271931, 0.0170119,  0.01271931,    5.9906788,    -1.204261 ],
        [ -5.82022227, 0.01937308, 0.01447391, 0.01937308, 0.01447391,    5.7806603,    -1.2118018 ],
        [ -5.18151608, 0.03119058, 0.02319602, 0.03119058, 0.02319602,    4.9970874,    -1.2457089 ],
        [ -4.68221894, 0.04516188, 0.0333488,  0.04516188, 0.0333488,     4.367733,     -1.2790618 ],
        [ -4.28448325, 0.06050063, 0.04424709, 0.06050063, 0.04424709,    3.8531174,    -1.3132945 ],
        [ -3.95494042, 0.07688325, 0.05554309, 0.07688325, 0.05554309,    3.4150251,    -1.3480039 ],
        [ -3.66917986, 0.09437263, 0.06714538, 0.09437263, 0.06714538,    3.0252073,    -1.3839729 ],
        [ -3.41757173, 0.11269814, 0.07872387, 0.11269814, 0.07872387,    2.672595,     -1.4216122 ],
        [ -3.19138195, 0.13177312, 0.09006663, 0.13177312, 0.09006663,    2.3469379,    -1.4603984 ],
        [ -2.98467584, 0.1515098,  0.1009523,  0.1515098,  0.1009523,     2.0411644,    -1.5010286 ],
        [ -2.79286243, 0.17185508, 0.11116985, 0.17185508, 0.11116985,    1.7493736,    -1.5421068 ],
        [ -2.61213456, 0.19279842, 0.12051491, 0.19279842, 0.12051491,    1.4671157,    -1.5829665 ],
        [ -2.43879704, 0.21441871, 0.12879375, 0.21441871, 0.12879375,    1.1892066,    -1.6208611 ],
        [ -2.26865664, 0.23694736, 0.13580652, 0.23694736, 0.13580652,    0.91049704,   -1.6532113 ],
        [ -2.09610295, 0.26088194, 0.14129565, 0.26088194, 0.14129565,    0.62258748,   -1.6740142 ],
        [ -1.90820867, 0.28780852, 0.1448542,  0.28780852, 0.1448542,     0.30692601,   -1.6692049 ],
        [ -1.72669697, 0.31420018, 0.14545214, 0.31420018, 0.14545214,    0.0058377995, -1.6274942 ],
        [ -1.56072874, 0.33820241, 0.14336436, 0.33820241, 0.14336436,    -0.25952652,  -1.5551588 ],
        [ -1.40530725, 0.36018527, 0.1391686,  0.36018527, 0.1391686,     -0.49486732,  -1.4681958 ],
        [ -1.25763494, 0.38032876, 0.13337462, 0.38032876, 0.13337462,    -0.70522292,  -1.3881927 ],
        [ -1.11495258, 0.39887533, 0.12639789, 0.39887533, 0.12639789,    -0.89829153,  -1.3295803 ],
        [ -0.97285613, 0.41628118, 0.11845377, 0.41628118, 0.11845377,    -1.0838868,   -1.2949341 ],
        [ -0.83130161, 0.43245001, 0.10992027, 0.43245001, 0.10992027,    -1.2656129,   -1.2822543 ],
        [ -0.68866819, 0.44749719, 0.10105362, 0.44749719, 0.10105362,    -1.4482847,   -1.2855411 ],
        [ -0.54366038, 0.46149829, 0.09208615, 0.46149829, 0.09208615,    -1.6354107,   -1.2989368 ],
        [ -0.39493137, 0.47452862, 0.08320932, 0.47452862, 0.08320932,    -1.8298941,   -1.3178651 ],
        [ -0.24134056, 0.48663891, 0.07459268, 0.48663891, 0.07459268,    -2.0339285,   -1.3391983 ],
        [ -0.08160699, 0.49788622, 0.06636565, 0.49788622, 0.06636565,    -2.2496326,   -1.3608986 ],
        [ 0.08554899,  0.50832041, 0.05862887, 0.50832041, 0.05862887,    -2.4789516,   -1.3816682 ],
        [ 0.26146194,  0.51798861, 0.05145368, 0.51798861, 0.05145368,    -2.7238158,   -1.4007915 ],
        [ 0.44767659,  0.52694235, 0.04487987, 0.52694235, 0.04487987,    -2.9864041,   -1.4179226 ],
        [ 0.64573868,  0.53522535, 0.03892811, 0.53522535, 0.03892811,    -3.2688808,   -1.4329349 ],
        [ 0.85763865,  0.54289139, 0.03359077, 0.54289139, 0.03359077,    -3.5740458,   -1.4458864 ],
        [ 1.0854073,   0.54998466, 0.02885012, 0.54998466, 0.02885012,    -3.904781,    -1.4568919 ],
        [ 1.32283905,  0.55633822, 0.02480563, 0.55633822, 0.02480563,    -4.2518877,   -1.4658266 ],
        [ 1.5698216,   0.56202674, 0.02137595, 0.56202674, 0.02137595,    -4.6149245,   -1.4730144 ],
        [ 1.82826335,  0.56716095, 0.01845692, 0.56716095, 0.01845692,    -4.9964583,   -1.4787893 ],
        [ 2.09982101,  0.57182309, 0.01596588, 0.57182309, 0.01596588,    -5.3987492,   -1.4834183 ],
        [ 2.38623051,  0.5760799,  0.01383371, 0.5760799,  0.01383371,    -5.824218,    -1.4871311 ],
        [ 2.68901128,  0.57998193, 0.01200447, 0.57998193, 0.01200447,    -6.2750065,   -1.4900945 ],
        [ 3.00995311,  0.58357325, 0.01042981, 0.58357325, 0.01042981,    -6.7536709,   -1.4924385 ],
        [ 3.35080543,  0.58688855, 0.00906982, 0.58688855, 0.00906982,    -7.2627347,   -1.4942889 ],
        [ 3.71316011,  0.58995435, 0.00789164, 0.58995435, 0.00789164,    -7.8045002,   -1.495742 ],
        [ 4.09948777,  0.59279825, 0.0068654,  0.59279825, 0.0068654,     -8.3826008,   -1.4968761 ],
        [ 4.511045,    0.59543333, 0.0059695,  0.59543333, 0.0059695,     -8.9988599,   -1.4977391 ],
        [ 4.95041288,  0.59787778, 0.00518304, 0.59787778, 0.00518304,    -9.6570851,   -1.4983945 ],
        [ 5.42030073,  0.60014503, 0.00448918, 0.60023296, 0.0050759759,  -10.361301,   -1.498883 ],
        [ 5.92400454,  0.60224642, 0.00387392, 0.60276038, 0.0048655028,  -11.116399,   -1.4992385 ],
        [ 6.46560698,  0.60419146, 0.00332571, 0.60522724, 0.004304928,   -11.928476,   -1.4994979 ],
        [ 7.04957736,  0.60598603, 0.00283575, 0.60760037, 0.0038270486,  -12.804202,   -1.4996759 ],
        [ 7.68357462,  0.60763993, 0.00239565, 0.60987821, 0.0033600329,  -13.75504,    -1.4997988 ],
        [ 8.37684531,  0.6091589,  0.00199948, 0.61205145, 0.0029126616,  -14.794842,   -1.4998816 ],
        [ 9.14302655,  0.61054924, 0.0016423,  0.6141107,  0.0024743005,  -15.94405,    -1.4999352 ],
        [ 10.00154701, 0.61181576, 0.00132061, 0.61604971, 0.002056063,   -17.231793,   -1.4999702 ],
        [ 10.98142848, 0.61296181, 0.0010313,  0.6178628,  0.0016564863,  -18.701602,   -1.4999941 ],
        [ 12.12928747, 0.61399011, 0.00077294, 0.61953841, 0.0012759168,  -20.423395,   -1.5000135 ],
        [ 13.52553678, 0.61490072, 0.00054479, 0.62106771, 0.00093152401, -22.517803,   -1.5000408 ],
        [ 15.32598756, 0.61569039, 0.00034722, 0.62248407, 0.00062905982, -25.21859,    -1.5001221 ],
        [ 17.89685356, 0.61634887, 0.00018257, 0.62367702, 0.00033075254, -29.075426,   -1.5001496 ],
        [ 22.51756936, 0.6168491,  5.751e-05,  0.62458326, 0.00010418468, -36.0067,     -1.5 ]
    ],

    impulse_table => [
        [
            -7.5046643, 13.815064,  -7.5052192,    0.03058362,  0,          -0.20030092,
            0,          -2.2866741, -0.0035429927, -0.45829177, 0.25467471, 0.99121891,
            -3.3953418
        ],
        [
            -6.9077552, 12.621249, -6.9087636,    0.037384125, 0,          -0.19985144,
            0,          -2.286679, -0.0047751543, -0.45784228, 0.25467417, 0.98626483,
            -3.1789503
        ],
        [
            -6.5022901, 11.810323,  -6.5038028,    0.042595008, 0,         -0.19935144,
            0,          -2.2866861, -0.0058483458, -0.45734228, 0.2546732, 0.98139069,
            -3.0338496
        ],
        [
            -6.214608, 11.234965,  -6.2166252,   0.046570236, 0,          -0.19885144,
            0,         -2.2866958, -0.006753088, -0.45684228, 0.25467189, 0.97691935,
            -2.9321221
        ],
        [
            -5.8091429, 10.424055,  -5.8121709,    0.052524233, 0,          -0.19785144,
            0,          -2.2867233, -0.0082708099, -0.45584228, 0.25466826, 0.96874558,
            -2.7909446
        ],
        [
            -5.2983173,  9.4024648,   -5.3033683, 0.060470363, 0, -0.19585143, 0, -2.2868114,
            -0.01067757, -0.45384228, 0.25465669, 0.95424946,  -2.6179029
        ],
        [
            -4.6051701, 8.0164587,  -4.6152973,   0.071445314, 0,          -0.19085144,
            0,          -2.2872245, -0.015100364, -0.44884228, 0.25460241, 0.9235135,
            -2.3956407
        ],
        [
            -3.9120229, 6.631315,   -3.9323692,   0.081190517, 0,          -0.18085149,
            0,          -2.2888779, -0.021355138, -0.43884228, 0.25438568, 0.87311287,
            -2.1962141
        ],
        [
            -3.5065578, 5.822293,   -3.5372039,   0.085311841, 0,          -0.17085194,
            0,          -2.2916368, -0.026154584, -0.42884228, 0.25402592, 0.83057573,
            -2.0959336
        ],
        [
            -2.9957322, 4.8066757,  -3.0471791,   0.087523978, 0,          -0.15085916,
            0,          -2.3004916, -0.033765176, -0.40884228, 0.25288681, 0.75931363,
            -1.9954448
        ],
        [
            -2.65926, 4.142578,   -2.7317036,   0.086524858, 0,          -0.13089787,
            0,        -2.3138426, -0.039949925, -0.38884228, 0.25121184, 0.70039533,
            -1.9509314
        ],
        [
            -2.302585,   3.4472782,   -2.4067109, 0.083139906, 0, -0.1011617, 0, -2.3424534,
            -0.04773649, -0.35914257, 0.24777896, 0.62791715,  -1.9289443
        ],
        [
            -2.0402208, 2.9455219,  -2.1760247,   0.079536716, 0,         -0.072095885,
            0,          -2.3815736, -0.054376012, -0.32974366, 0.2433847, 0.56932944,
            -1.9331459
        ],
        [
            -1.8971199, 2.6768915, -2.053929,    0.077511518, 0,          -0.053480754,
            0,          -2.413559, -0.058325625, -0.31094789, 0.24000893, 0.53616201,
            -1.9437697
        ],
        [
            -1.6094378, 2.1513529,  -1.8179776,   0.07447119,  0,          -0.011974844,
            0,          -2.5137017, -0.066733749, -0.26792531, 0.23039727, 0.46870357,
            -1.9844547
        ],
        [
            -1.3862943, 1.7604076,   -1.6449835, 0.074112125, 0, 0.019148306, 0, -2.6390235,
            -0.0729649, -0.23497899, 0.21972205, 0.41743317,  -2.0343202
        ],
        [
            -1.2039727, 1.4538502,  -1.5108821,   0.075406135, 0,         0.039862711,
            0,          -2.7807482, -0.076972073, -0.21421675, 0.2086155, 0.37741444,
            -2.0866212
        ],
        [
            -1.0498221, 1.2045138,  -1.4028729,   0.077405518, 0,          0.053203311,
            0,          -2.9284856, -0.079182989, -0.20347344, 0.19754273, 0.34546083,
            -2.13836
        ],
        [
            -0.91629067, 0.99611292, -1.3133827,   0.07956841, 0,          0.062042199,
            0,           -3.0744657, -0.080232857, -0.1991621, 0.18682086, 0.31943976,
            -2.1882043
        ],
        [
            -0.79850763, 0.81818874, -1.2375938,   0.081648814, 0,          0.068195744,
            0,           -3.214491,  -0.080619368, -0.19820721, 0.17665038, 0.29788304,
            -2.2356167
        ],
        [
            -0.69314711, 0.66367684, -1.172274,    0.083554611, 0,          0.072699421,
            0,           -3.3468487, -0.080645968, -0.1992607,  0.16714552, 0.27975481,
            -2.280449
        ],
        [
            -0.59783693, 0.52761403, -1.1151642,   0.085263832, 0,          0.076142524,
            0,           -3.4711612, -0.080482145, -0.20100188, 0.15835917, 0.26430783,
            -2.3227418
        ],
        [
            -0.51082556, 0.40639876, -1.0646311,   0.086783803, 0,          0.078871555,
            0,           -3.5876807, -0.080220683, -0.20359533, 0.15030201, 0.25099221,
            -2.3626235
        ],
        [
            -0.35667488, 0.19842039, -0.97873885,  0.089330632, 0,          0.082960843,
            0,           -3.7994805, -0.079586502, -0.20917306, 0.13628788, 0.22920432,
            -2.4358186
        ],
        [
            -0.22314349, 0.024904696, -0.90791852,  0.091351812, 0,          0.085921507,
            0,           -3.9869098,  -0.078935839, -0.21527265, 0.12479207, 0.21211669,
            -2.5013834
        ],
        [
            -0.10536045, -0.1233457, -0.84805402,  0.092979522, 0,          0.088196633,
            0,           -4.1542466, -0.078327185, -0.22087256, 0.11540591, 0.198338,
            -2.5605337
        ],
        [
            6.6000847e-08, -0.25237799, -0.79645446,  0.094311593, 0,          0.090019499,
            0,             -4.3049738,  -0.077775332, -0.22639687, 0.10771023, 0.18697375,
            -2.6142811
        ],
        [
            0.095310246, -0.36635638, -0.7512775,   0.095418477, 0,          0.091524862,
            0,           -4.441847,   -0.077280162, -0.23182772, 0.10131153, 0.17742452,
            -2.6634474
        ],
        [
            0.10785036, -0.38116651, -0.74543547,  0.095557386, 0,          0.091713939,
            0,          -4.4598799,  -0.077215668, -0.23282046, 0.10050784, 0.17621734,
            -2.6699492
        ],
        [
            0.26236433, -0.5602783, -0.67530329,  0.097146622, 0,           0.093891729,
            0,          -4.6822895, -0.076439389, -0.24204857, 0.091288229, 0.16222539,
            -2.7505718
        ],
        [
            0.3364723, -0.64407758, -0.64282433,  0.097832721, 0,           0.094847187,
            0,         -4.7889784,  -0.076081852, -0.24704038, 0.087285717, 0.15605976,
            -2.7895121
        ],
        [
            0.40546517,  -0.72093953, -0.613222,  0.09843022, 0, 0.095690975, 0, -4.8882374,
            -0.07575911, -0.25155041, 0.08378303, 0.1506136,  -2.825883
        ],
        [
            0.53062832, -0.85769808, -0.56099577,  0.09941974,  0,           0.097120588,
            0,          -5.0680131,  -0.075200231, -0.26020398, 0.077934492, 0.14140674,
            -2.8920817
        ],
        [
            0.69314725, -1.0304963, -0.49581505,  0.10054061,  0,           0.098805222,
            0,          -5.3006141, -0.074528667, -0.27262026, 0.071207708, 0.13062671,
            -2.9782936
        ],
        [
            0.99325184, -1.3370942,  -0.3823301,  0.10220331, 0, 0.10149175, 0, -5.7268096,
            -0.0734509, -0.29840889, 0.060909403, 0.11368161, -3.1376638
        ],
        [
            1.0986124,    -1.4413431, -0.3443506,  0.10268287, 0, 0.10232403, 0, -5.8752635,
            -0.073120675, -0.3087307, 0.057833395, 0.10850679, -3.1935444
        ],
        [
            1.2089604,    -1.5488429, -0.30549494, 0.10313631, 0, 0.10314179, 0, -6.0300504,
            -0.072800146, -0.3198716, 0.054864166, 0.10345846, -3.2519878
        ],
        [
            1.3862944, -1.7182834, -0.24486076,  0.10377318,  0,           0.1043507,
            0,         -6.2772901, -0.072335912, -0.33892453, 0.050568244, 0.096058064,
            -3.3456864
        ],
        [
            1.609438, -1.9262929, -0.1713844,   0.10443877,  0,           0.10570787,
            0,        -6.5857677, -0.071833377, -0.36542164, 0.045870446, 0.087828557,
            -3.4631348
        ],
        [
            1.7917595, -2.0924525, -0.11339346,  0.10489005,  0,           0.10669742,
            0,         -6.8356964, -0.071482643, -0.39016934, 0.042519601, 0.081867294,
            -3.5586888
        ],
        [
            1.9459102, -2.2305746, -0.065624586, 0.10521776,  0,           0.10746014,
            0,         -7.0455866, -0.071223473, -0.41243057, 0.039974609, 0.077287114,
            -3.6391865
        ],
        [
            2.0794416, -2.3486345, -0.025085709, 0.10546723,  0,           0.10807018,
            0,         -7.2263962, -0.071023597, -0.43163238, 0.037955626, 0.073620622,
            -3.7087014
        ],
        [
            2.3025852, -2.5429676, 0.041107985,  0.10582404,  0,           0.10899704,
            0,         -7.5265805, -0.070735041, -0.46981303, 0.034914146, 0.068041331,
            -3.8244405
        ],
        [
            2.7080503, -2.8879596, 0.15719065,   0.10631963,  0,           0.11042199,
            0,         -8.0663232, -0.070330923, -0.55134829, 0.030254893, 0.059362643,
            -4.0334814
        ],
        [
            2.9957323, -3.127446,  0.23687986,   0.10658119,  0,           0.11126377,
            0,         -8.4453543, -0.070117583, -0.61891009, 0.027484461, 0.054126912,
            -4.1809288
        ],
        [
            3.4011974, -3.4589803, 0.34625369,   0.10685847,  0,           0.112255,
            0,         -8.9749027, -0.069893504, -0.73782743, 0.024157862, 0.047767918,
            -4.3877209
        ],
        [
            3.9122247, -3.868884,  0.48035925,   0.10709974,  0,           0.11324205,
            0,         -9.6359286, -0.069704781, -0.91680725, 0.020704842, 0.040695407,
            -4.6554497
        ],
        [
            4.605325,    -4.4143666, 0.6576205,   0.10730353,  0, 0.11422826, 0, -10.523815,
            -0.06954743, -1.2499155, 0.016980105, 0.033613941, -5.0018749
        ],
        [
            5.2983642, -4.9513487, 0.83148456,   0.10742183, 0,           0.11492569,
            0,         -11.404513, -0.069420298, -1.7102303, 0.014046182, 0.027922741,
            -5.3484155
        ],
        [
            6.2146752,   -5.6528685, 1.0585116,   0.10750665,  0, 0.11556301, 0, -12.561784,
            -0.06939744, -2.6512794, 0.011022592, 0.021978473, -5.8066357
        ],
        [
            6.9078606,    -6.1794227, 1.2292005,    0.10754101, 0, 0.11591403, 0, -13.43376,
            -0.068610177, -3.6169439, 0.0092126174, 0.01839214, -6.1532703
        ],
        [
            7.6009326, -6.7036077, 1.3995268,    0.10756125, 0,            0.11609718,
            0,         -14.303682, -0.069333959, -5.1703899, 0.0077166516, 0.015416838,
            -6.4998363
        ],
        [
            8.517486,    -7.394571,  1.6247563,    0.10757578, 0, 0.11629176, 0, -15.45225,
            -0.12350245, -7.9589135, 0.0061172834, 0.01222801, -6.9581384
        ],
        [
            9.2104385, -7.9158948, 1.7952297,    0.10758154, 0,            0.11638956,
            0,         -16.319733, -0.061788257, -10.059227, 0.0051370836, 0.010270894,
            -7.3046269
        ],
        [
            10.819818, -9.1247808, 2.1921566,    0.10758651, 0,            0.11651991,
            0,         -18.332931, -0.085521777, -21.129081, 0.0034300104, 0.0068593685,
            -8.1093349
        ],
        [
            11.513059, -9.6450742, 2.3635887,    0.10758666, 0,            0.11655076,
            0,         -19.199768, -0.057531016, -29.567066, 0.0028833409, 0.0057663568,
            -8.455963
        ],
        [
            13.122524, -10.852578, 2.7625318,    0.10758384, 0,            0.1165919,
            0,         -21.211916, -0.012658072, -14.43387,  0.0019275446, 0.0038550272,
            -9.2607202
        ],
        [
            13.815546, -11.372421, 2.9346564,    0.1075811,  0,            0.11660164,
            0,         -22.078253, -0.010718209, -17.304151, 0.0016208142, 0.0032416004,
            -9.6072483
        ],
        [
            16.118208, -13.099519, 3.5076584,     0.10755967, 0,             0.11661771,
            0,         -24.956654, -0.0061432413, -31.454649, 0.00091135726, 0.0018227256,
            -10.758705
        ],
        [
            18.420867, -14.826533, 4.081803,     0.10749119, 0,             0.11662242,
            0,         -27.835007, -0.069326521, -1126.8729, 0.00051245015, 0.0010249454,
            -11.910431
        ]
    ],

    tail_shock_table => [
        [ 4.4534283, -1.239,     -0.034155594,  -0.034155594 ],
        [ 4.5858373, -1.2949912, -0.041325127,  -0.025269684 ],
        [ 5.2868146, -1.6223612, -0.042766228,  -0.012638045 ],
        [ 6.2088948, -2.1488997, -0.037374922,  -0.0067718585 ],
        [ 6.9044366, -2.6370725, -0.032711789,  -0.0046282568 ],
        [ 7.5989039, -3.2207205, -0.028309545,  -0.0032975526 ],
        [ 8.5164703, -4.1719645, -0.02316082,   -0.0021776605 ],
        [ 9.2098363, -5.0562552, -0.019810205,  -0.0016147546 ],
        [ 10.819639, -7.8325816, -0.013654004,  -0.0008294511 ],
        [ 11.512953, -9.4295553, -0.01160002,   -0.00062762149 ],
        [ 13.122492, -14.444386, -0.0079173114, -0.00033102989 ],
        [ 13.815527, -17.318679, -0.0067046181, -0.00025285575 ],
        [ 16.118205, -31.51183,  -0.0038464957, -0.00010535648 ]
    ],

    blast_info => [
        0.20085138, 0,          0.45884809, -0.15100366, 0,          0,          1.2489497,   -0.16929751,
        0.63915761, 0.16334222, 0.73871864, -0.11770941, 0.33279978, 0.25467474, 0.067224792, 0,
        84.682,     -1.239,     136.10592,  -1.4386906
    ],

};
1;
