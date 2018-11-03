package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='P', gamma=4.5
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'P4x5'} = {

    table_name => 'P4x5',
    symmetry   => 0,
    gamma      => 4.5,

    shock_table_info => [ 0, 4.5, 1.73e-06, 4000, 1.36e-06, 100, 9.22e-07, 3.05e-07, 5e-07 ],

    shock_table => [
        [ -10.31633794, 12.00046241, -0.99998995, -10.31845357, 0.99894108 ],
        [ -9.33039601,  11.01453573, -0.99997307, -9.33386197,  0.99826407 ],
        [ -8.76744692,  10.45160471, -0.99995489, -8.77204214,  0.99769723 ],
        [ -7.37163306,  9.05592024,  -0.99982807, -7.38088811,  0.99535201 ],
        [ -6.46441654,  8.1489576,   -0.99957447, -6.4790204,   0.99264846 ],
        [ -5.72339808,  7.40840599,  -0.9991085,  -5.74461494,  0.98929034 ],
        [ -5.08859644,  6.77439377,  -0.99832215, -5.11784099,  0.98519337 ],
        [ -4.54001785,  6.22703919,  -0.99710618, -4.57864169,  0.98038257 ],
        [ -4.05391974,  5.7427478,   -0.99531905, -4.10337902,  0.97479932 ],
        [ -3.61530679,  5.30670253,  -0.99279499, -3.67717275,  0.9683826 ],
        [ -3.21293252,  4.90787936,  -0.98933437, -3.28894426,  0.96104928 ],
        [ -2.83867446,  4.53842953,  -0.98470369, -2.93077486,  0.95270584 ],
        [ -2.48503496,  4.19122039,  -0.97860773, -2.5954894,   0.94321022 ],
        [ -2.14450987,  3.8592692,   -0.97064783, -2.27610191,  0.9323403 ],
        [ -1.80898971,  3.53526267,  -0.96025066, -1.96534034,  0.91974748 ],
        [ -1.466358,    3.20851755,  -0.94642893, -1.65270947,  0.90476066 ],
        [ -1.08338587,  2.84971931,  -0.92648952, -1.30984368,  0.8853107 ],
        [ -0.73322332,  2.52915473,  -0.90372297, -1.00332291,  0.86503771 ],
        [ -0.41085724,  2.24172382,  -0.87894757, -0.72773754,  0.84444323 ],
        [ -0.10880466,  1.98012365,  -0.85277015, -0.47576635,  0.82375489 ],
        [ 0.17182093,   1.74447221,  -0.82643289, -0.24740332,  0.80366423 ],
        [ 0.44418207,   1.52301019,  -0.79966483, -0.03122676,  0.7837185 ],
        [ 0.71688938,   1.30865394,  -0.7723733,  0.17976287,   0.76366872 ],
        [ 0.99824983,   1.0952827,   -0.74443283, 0.39174832,   0.74326959 ],
        [ 1.28643575,   0.8847683,   -0.71673697, 0.60301191,   0.72302882 ],
        [ 1.58310748,   0.67617464,  -0.68979427, 0.81454296,   0.70318875 ],
        [ 1.89172785,   0.4673455,   -0.66390363, 1.02854244,   0.68386347 ],
        [ 2.2162734,    0.2559451,   -0.63930493, 1.24740434,   0.66514756 ],
        [ 2.56187869,   0.03908731,  -0.61616167, 1.47410794,   0.64709833 ],
        [ 2.93360434,   -0.18585609, -0.59467361, 1.71137231,   0.62981862 ],
        [ 3.33874971,   -0.42266748, -0.57495477, 1.96312572,   0.61335719 ],
        [ 3.74109079,   -0.65063117, -0.55875339, 2.20699307,   0.59923264 ],
        [ 4.14882363,   -0.87561531, -0.54528857, 2.44874503,   0.58691509 ],
        [ 4.5683949,    -1.10196016, -0.53404869, 2.6926606,    0.57605918 ],
        [ 5.0066244,    -1.33385867, -0.5246555,  2.94293311,   0.56640255 ],
        [ 5.46036439,   -1.57010537, -0.51698857, 3.19795681,   0.55793002 ],
        [ 5.94068391,   -1.8168486,  -0.5107035,  3.464073,     0.55036914 ],
        [ 6.4515202,    -2.07638022, -0.50564797, 3.74344984,   0.54363572 ],
        [ 6.99098665,   -2.34804363, -0.50171838, 4.03507702,   0.53771654 ],
        [ 7.5697284,    -2.63748936, -0.49872145, 4.34470688,   0.53246311 ],
        [ 8.19992407,   -2.95103764, -0.4965172,  4.67873443,   0.5277715 ],
        [ 8.87607306,   -3.28621208, -0.49503149, 5.0341541,    0.52367746 ],
        [ 9.61523889,   -3.65175053, -0.49412787, 5.41985547,   0.52006734 ],
        [ 10.44158224,  -4.05986208, -0.49370764, 5.84822883,   0.51685068 ],
        [ 11.36510001,  -4.51577645, -0.49369696, 6.3241904,    0.51402306 ],
        [ 12.41467149,  -5.0340937,  -0.49401943, 6.86232824,   0.51152778 ],
        [ 13.62894079,  -5.63430887, -0.49460709, 7.4820572,    0.50931683 ],
        [ 15.08254084,  -6.35384608, -0.4954041,  8.22088058,   0.5073234 ],
        [ 16.43092891,  -7.02233357, -0.49612595, 8.90395812,   0.50591072 ]
    ],

    energy_table => [
        [ -10.31633794, 0.00010188, 7.714e-05,  0.00010188, 7.714e-05,  11.481895,   -0.99281538 ],
        [ -9.33039601,  0.00021435, 0.00016109, 0.00021435, 0.00016109, 10.500741,   -0.99739088 ],
        [ -8.76744692,  0.00032689, 0.00024436, 0.00032689, 0.00024436, 9.9385386,   -0.99986752 ],
        [ -7.37163306,  0.00092016, 0.00067584, 0.00092016, 0.00067584, 8.5387781,   -1.0056063 ],
        [ -6.46441654,  0.00178238, 0.0012878,  0.00178238, 0.0012878,  7.6248369,   -1.009111 ],
        [ -5.72339808,  0.00303082, 0.00215134, 0.00303082, 0.00215134, 6.8760389,   -1.0118288 ],
        [ -5.08859644,  0.00473622, 0.00329701, 0.00473622, 0.00329701, 6.2330048,   -1.014071 ],
        [ -4.54001785,  0.006912,   0.00471106, 0.006912,   0.00471106, 5.6761846,   -1.0159506 ],
        [ -4.05391974,  0.0095918,  0.00638853, 0.0095918,  0.00638853, 5.1819337,   -1.0175743 ],
        [ -3.61530679,  0.01280064, 0.00831287, 0.01280064, 0.00831287, 4.7352948,   -1.0190076 ],
        [ -3.21293252,  0.01656512, 0.01046177, 0.01656512, 0.01046177, 4.3250105,   -1.0202972 ],
        [ -2.83867446,  0.02090843, 0.01280315, 0.02090843, 0.01280315, 3.9429337,   -1.0214756 ],
        [ -2.48503496,  0.02587063, 0.0153044,  0.02587063, 0.0153044,  3.5815044,   -1.022571 ],
        [ -2.14450987,  0.03152477, 0.01793501, 0.03152477, 0.01793501, 3.2331151,   -1.0236095 ],
        [ -1.80898971,  0.03799792, 0.02066564, 0.03799792, 0.02066564, 2.8895032,   -1.0246397 ],
        [ -1.466358,    0.04556236, 0.0234826,  0.04556236, 0.0234826,  2.5382463,   -1.0252323 ],
        [ -1.08338587,  0.05513797, 0.02647902, 0.05513797, 0.02647902, 2.1455857,   -1.0243521 ],
        [ -0.73322332,  0.06484236, 0.02887757, 0.06484236, 0.02887757, 1.7871988,   -1.0216542 ],
        [ -0.41085724,  0.07444802, 0.03063332, 0.07444802, 0.03063332, 1.4583961,   -1.0174789 ],
        [ -0.10880466,  0.0838868,  0.0317769,  0.0838868,  0.0317769,  1.1517681,   -1.0122752 ],
        [ 0.17182093,   0.09289603, 0.03235049, 0.09289603, 0.03235049, 0.86844672,  -1.0066804 ],
        [ 0.44418207,   0.10172999, 0.03244346, 0.10172999, 0.03244346, 0.59503959,  -1.0009849 ],
        [ 0.71688938,   0.11053991, 0.03209701, 0.11053991, 0.03209701, 0.32284343,  -0.9954628 ],
        [ 0.99824983,   0.11947123, 0.03132443, 0.11947123, 0.03132443, 0.043532557, -0.99035399 ],
        [ 1.28643575,   0.12833945, 0.0301652,  0.12833945, 0.0301652,  -0.24117646, -0.98603545 ],
        [ 1.58310748,   0.13707315, 0.02866902, 0.13707315, 0.02866902, -0.5331251,  -0.9827208 ],
        [ 1.89172785,   0.14565019, 0.02688301, 0.14565019, 0.02688301, -0.83597647, -0.98053283 ],
        [ 2.2162734,    0.15404851, 0.02485346, 0.15404851, 0.02485346, -1.1539398,  -0.97958174 ],
        [ 2.56187869,   0.16225327, 0.02262271, 0.16225327, 0.02262271, -1.492438,   -0.9801019 ],
        [ 2.93360434,   0.17021798, 0.02023968, 0.17021798, 0.02023968, -1.8570333,  -0.98276008 ],
        [ 3.33874971,   0.17790842, 0.01774857, 0.17790842, 0.01774857, -2.2560515,  -0.98864764 ],
        [ 3.74109079,   0.1845795,  0.01544446, 0.1845795,  0.01544446, -2.6553318,  -1.0186672 ],
        [ 4.14882363,   0.19043546, 0.01331644, 0.19043546, 0.01331644, -3.0815324,  -1.0353077 ],
        [ 4.5683949,    0.19560358, 0.01135828, 0.19560358, 0.01135828, -3.5116063,  -1 ],
        [ 5.0066244,    0.20017852, 0.00956262, 0.20017852, 0.00956262, -3.9498358,  -1 ],
        [ 5.46036439,   0.20414407, 0.00795841, 0.20414407, 0.00795841, -4.4035757,  -1 ],
        [ 5.94068391,   0.20761084, 0.00651882, 0.20761084, 0.00651882, -4.8838953,  -1 ],
        [ 6.4515202,    0.2106053,  0.00524638, 0.2106053,  0.00524638, -5.3947316,  -1 ],
        [ 6.99098665,   0.21312971, 0.00415172, 0.21312971, 0.00415172, -5.934198,   -1 ],
        [ 7.5697284,    0.21525076, 0.00321537, 0.21525076, 0.00321537, -6.5129398,  -1 ],
        [ 8.19992407,   0.2170164,  0.00242338, 0.2170164,  0.00242338, -7.1431354,  -1 ],
        [ 8.87607306,   0.21842732, 0.00178147, 0.21842732, 0.00178147, -7.8192844,  -1 ],
        [ 9.61523889,   0.21954365, 0.00126717, 0.21954365, 0.00126717, -8.5584502,  -1 ],
        [ 10.44158224,  0.220413,   0.00086212, 0.220413,   0.00086212, -9.3847936,  -1 ],
        [ 11.36510001,  0.22105888, 0.00055811, 0.22105888, 0.00055811, -10.308311,  -1 ],
        [ 12.41467149,  0.22152028, 0.00033892, 0.22152028, 0.00033892, -11.357883,  -1 ],
        [ 13.62894079,  0.2218324,  0.00018938, 0.2218324,  0.00018938, -12.572152,  -1 ],
        [ 15.08254084,  0.22203026, 9.383e-05,  0.22203026, 9.383e-05,  -14.025752,  -1 ],
        [ 16.43092891,  0.22212309, 4.87e-05,   0.22212309, 4.87e-05,   -15.37414,   -1 ]
    ],

    impulse_table => [
        [
            -12.130771,   13.814887,  -12.131624, 3.1393491,  0, 0, 0, -7.35396,
            0.0017592223, -84.970134, 0.32513053, 0.99995392, -9.3075002
        ],
        [
            -2.3025851,   4.0130261,  -2.4239036, 2.9811223, 0, 0, 0, -2.5288146,
            0.0017592241, -84.870139, 0.31824651, 0.870468,  -2.295701
        ],
        [
            -1.6094379,   3.3443834,  -1.7826343, 2.9315795,  0, 0, 0, -2.2583307,
            0.0017592311, -84.770139, 0.31186221, 0.80972788, -1.9113955
        ],
        [
            -1.3862943,   3.1328918,  -1.5804234, 2.9149193,  0, 0, 0, -2.1816531,
            0.0017592368, -84.720139, 0.3088371,  0.78640703, -1.795856
        ],
        [
            -1.2039728,   2.9618581,  -1.4169915, 2.9014477,  0, 0, 0, -2.1236816,
            0.0017592442, -84.670139, 0.30591404, 0.76601638, -1.7047453
        ],
        [
            -0.91629072,  2.6957551,  -1.1626893, 2.8810559,  0, 0, 0, -2.0418678,
            0.0017592641, -84.570139, 0.30034985, 0.73152746, -1.5673676
        ],
        [
            -0.69314717, 2.4929947,  -0.96870487, 2.8665357,  0, 0, 0, -1.9874409,
            0.001759291, -84.470139, 0.29512647,  0.70299355, -1.4664693
        ],
        [
            -0.51082561,  2.3299965,  -0.81248491, 2.8558932,  0, 0, 0, -1.9493516,
            0.0017593253, -84.370139, 0.29020769,  0.67867041, -1.3878546
        ],
        [
            -0.35667493, 2.1942221,  -0.68208169, 2.847959,   0, 0, 0, -1.9218724,
            0.001759367, -84.270139, 0.28556283,  0.65749891, -1.3241406
        ],
        [
            -0.22314354,  2.0782159,  -0.57041073, 2.8419907,  0, 0, 0, -1.9016905,
            0.0017594163, -84.170139, 0.28116562,  0.63878126, -1.2710213
        ],
        [
            1.2663469e-08, 1.8878812,  -0.38655662,  2.8341097,  0,          0,
            0,             -1.8756566, 0.0017595381, -83.970139, 0.27302643, 0.60689318,
            -1.1865941
        ],
        [
            0.26236428,   1.6700421,  -0.17493554, 2.8284911,  0, 0, 0, -1.8575194,
            0.0017597791, -83.670139, 0.26219472,  0.56878242, -1.0943041
        ],
        [
            0.53062826,   1.4542554,  0.036247471, 2.8270017,  0, 0, 0, -1.8528085,
            0.0017602108, -83.270139, 0.24978069,  0.52961921, -1.0076794
        ],
        [
            0.69314719,   1.32702,    0.16161103, 2.8281646,  0, 0, 0, -1.8565705,
            0.0017606178, -82.970139, 0.24166649, 0.50601613, -0.95893949
        ],
        [
            0.99325179,   1.0990046,  0.38803253, 2.834126,   0, 0, 0, -1.8759807,
            0.0017618462, -82.270139, 0.22571109, 0.46309829, -0.87612325
        ],
        [
            1.0986123,    1.0210619,  0.46598547, 2.837289,   0, 0, 0, -1.8864338,
            0.0017624924, -81.970139, 0.21987118, 0.44832743, -0.84919048
        ],
        [
            1.2089604,    0.94058089, 0.54678768, 2.8411345,  0, 0, 0, -1.8992987,
            0.0017633373, -81.620139, 0.21365627, 0.43306407, -0.82213699
        ],
        [
            1.3862944,    0.81366069, 0.67487171, 2.8483402,  0, 0, 0, -1.92389,
            0.0017651672, -80.970139, 0.20351661, 0.40903693, -0.78106468
        ],
        [
            1.6094379,   0.65804237, 0.83303577, 2.8589258,  0, 0, 0, -1.9612507,
            0.001768647, -79.970139, 0.19062728, 0.37979701, -0.73342674
        ],
        [
            1.7917595,   0.53411919, 0.95987391, 2.8685623,  0, 0, 0, -1.9966675,
            0.001772938, -78.970139, 0.18011245, 0.35681755, -0.69768892
        ],
        [
            1.9459102,    0.43149049, 1.0655077,  2.8772216,  0, 0, 0, -2.0297696,
            0.0017780484, -77.970139, 0.17131153, 0.33807452, -0.66960024
        ],
        [
            2.0794416,    0.34410322, 1.1558687,  2.8849905,  0, 0, 0, -2.0606196,
            0.0017839879, -76.970139, 0.16379705, 0.32236799, -0.64677921
        ],
        [
            2.3025851,    0.20102909, 1.3046111,  2.8982816,  0, 0, 0, -2.1163339,
            0.0017984031, -74.970139, 0.15154903, 0.29725087, -0.61163461
        ],
        [
            2.7080502, -0.050324198, 1.5681757,    2.9222198,  0,          0,
            0,         -2.2293854,   0.0018500105, -69.970139, 0.13058746, 0.25530495,
            -0.55677944
        ],
        [
            2.9957323, -0.22270019, 1.7504182,    2.9377627,  0,          0,
            0,         -2.3175818,  0.0019262341, -64.970139, 0.11691583, 0.22842904,
            -0.52450483
        ],
        [
            3.4011974, -0.45848675, 2.0013555,    2.9544139,  0,          0,
            0,         -2.4513448,  0.0021723362, -54.970139, 0.09948575, 0.194487,
            -0.48826793
        ],
        [
            3.912023, -0.74562238, 2.3089557,    2.9505499,  0,           0,
            0,        -2.6330318,  0.0033912985, -34.970138, 0.080585773, 0.15785881,
            -0.46175795
        ],
        [
            4.2484953, -0.92981999, 2.5071068,    2.8771167,  0,           0,
            0,         -2.7600909,  0.0091041253, -14.970138, 0.069897945, 0.13703991,
            -0.47540088
        ]
    ],

    blast_info => [ 0, 0, 0, 0, 0, 0, 3.1371825, 0, 0, 0, 0, 0, 0.24998798, 0.32513054, 0, 0, 0, 0, 0, 0 ],

};
1;
