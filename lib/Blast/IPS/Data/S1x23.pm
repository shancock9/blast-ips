package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='S', gamma=1.23
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'S1x23'} = {

    table_name => 'S1x23',
    symmetry   => 2,
    gamma      => 1.23,

    shock_table_info => [ 2, 1.23, 1.15e-06, 8000, 2e-07, 10, 3.17e-07, 3.3e-07, 5e-07, 2027.8, -0.5482 ],

    shock_table => [
        [ -4.78132852,  12.00017659,   -2.99997967, -4.78237035, 0.99843645 ],
        [ -4.66241093,  11.64342599,   -2.99997096, -4.66365632, 0.99813076 ],
        [ -4.26064413,  10.43814761,   -2.99990542, -4.26292056, 0.99658156 ],
        [ -3.87600953,  9.28431093,    -2.99970067, -3.88006643, 0.99390284 ],
        [ -3.48532827,  8.11248567,    -2.99905351, -3.4926289,  0.98901221 ],
        [ -3.22050707,  7.31840582,    -2.99790474, -3.2313857,  0.98360364 ],
        [ -3.00775431,  6.68077171,    -2.99603542, -3.02274895, 0.97736643 ],
        [ -2.8150451,   6.1036767,     -2.99295661, -2.83510659, 0.96967037 ],
        [ -2.64540276,  5.59630342,    -2.98834427, -2.67133444, 0.96073663 ],
        [ -2.49490817,  5.14702847,    -2.98182755, -2.52748104, 0.95061763 ],
        [ -2.3538195,   4.72693714,    -2.9725478,  -2.39416533, 0.93877426 ],
        [ -2.2242025,   4.342404,      -2.9600803,  -2.27332102, 0.92543035 ],
        [ -2.10192565,  3.98141326,    -2.94349534, -2.16106245, 0.91025248 ],
        [ -1.98446373,  3.63688192,    -2.92167445, -2.05513285, 0.89290389 ],
        [ -1.86890539,  3.30083881,    -2.89298443, -1.95308058, 0.8728133 ],
        [ -1.75099557,  2.9618908,     -2.85462077, -1.85154299, 0.84887436 ],
        [ -1.62099243,  2.59419332,    -2.79979982, -1.74313403, 0.81813823 ],
        [ -1.49813463,  2.2540795,     -2.73470729, -1.64462527, 0.78479269 ],
        [ -1.38562252,  1.95028822,    -2.66367638, -1.55821247, 0.75073385 ],
        [ -1.28014194,  1.67324064,    -2.58803257, -1.48083064, 0.71609154 ],
        [ -1.17945025,  1.41657251,    -2.50912631, -1.4104765,  0.68105098 ],
        [ -1.08169596,  1.17522891,    -2.42806207, -1.34561851, 0.64574737 ],
        [ -0.98494944,  0.94430822,    -2.34542475, -1.28486342, 0.61016076 ],
        [ -0.8885209,   0.72214151,    -2.26253975, -1.22774044, 0.57465594 ],
        [ -0.79143846,  0.50649821,    -2.18025973, -1.17366683, 0.53945269 ],
        [ -0.69281126,  0.2954843,     -2.09930521, -1.12218339, 0.50476815 ],
        [ -0.59180969,  0.08747847,    -2.02028855, -1.07293047, 0.47081724 ],
        [ -0.4875348,   -0.11914506,   -1.94365097, -1.0255777,  0.43777441 ],
        [ -0.37925484,  -0.3255544,    -1.86986817, -0.97992642, 0.40585161 ],
        [ -0.26608706,  -0.53310425,   -1.79923604, -0.93575791, 0.37519645 ],
        [ -0.14722441,  -0.74290356,   -1.73204684, -0.89292819, 0.34596045 ],
        [ -0.02158831,  -0.95643681,   -1.66840432, -0.85123895, 0.31821323 ],
        [ 0.11177408,   -1.17485957,   -1.60846653, -0.8105821,  0.29204436 ],
        [ 0.25409448,   -1.39968777,   -1.55225194, -0.77080518, 0.26747694 ],
        [ 0.40658824,   -1.63230169,   -1.49979355, -0.7318073,  0.24453415 ],
        [ 0.56746807,   -1.86964382,   -1.45193923, -0.69419278, 0.22359023 ],
        [ 0.73480193,   -2.10892575,   -1.40906992, -0.65838885, 0.20481125 ],
        [ 0.90970709,   -2.35192615,   -1.37056239, -0.62408056, 0.18792111 ],
        [ 1.09380544,   -2.6009629,    -1.33581013, -0.59092616, 0.17264897 ],
        [ 1.28879861,   -2.85828703,   -1.30433162, -0.55864802, 0.15877774 ],
        [ 1.49524617,   -3.12455003,   -1.27589974, -0.52720049, 0.14620242 ],
        [ 1.71424664,   -3.40108515,   -1.25022252, -0.49646433, 0.13478987 ],
        [ 1.94781372,   -3.69030486,   -1.22695969, -0.46622875, 0.12438564 ],
        [ 2.19679148,   -3.99310138,   -1.20594061, -0.43646988, 0.11491179 ],
        [ 2.46345995,   -4.31207755,   -1.18689967, -0.40701214, 0.10624832 ],
        [ 2.74913587,   -4.64861739,   -1.16968517, -0.37782084, 0.09832749 ],
        [ 3.05566257,   -5.0047042,    -1.15412796, -0.34882177, 0.09107485 ],
        [ 3.38510706,   -5.38254261,   -1.14007173, -0.31994235, 0.0844227 ],
        [ 3.73992776,   -5.78474369,   -1.1273679,  -0.29110091, 0.0783077 ],
        [ 4.1263777,    -6.21811017,   -1.11578718, -0.26196587, 0.07262635 ],
        [ 4.5284835,    -6.66468543,   -1.10568985, -0.2338047,  0.06757139 ],
        [ 4.97490057,   -7.15612576,   -1.09630134, -0.20474012, 0.06276723 ],
        [ 5.45864664,   -7.68434114,   -1.0878087,  -0.17548087, 0.05831866 ],
        [ 5.98374325,   -8.25346574,   -1.08011879, -0.14596991, 0.05419121 ],
        [ 6.55460664,   -8.86801608,   -1.07314816, -0.11615766, 0.05035515 ],
        [ 7.17644864,   -9.533317,     -1.06681868, -0.08598445, 0.04678278 ],
        [ 7.85467748,   -10.25485366,  -1.06106445, -0.055414,   0.04345223 ],
        [ 8.59569863,   -11.03912424,  -1.05582365, -0.02439733, 0.04034265 ],
        [ 9.40571556,   -11.8923662,   -1.05104717, 0.00707405,  0.03743923 ],
        [ 10.29253033,  -12.82245811,  -1.04668598, 0.0390407,   0.03472563 ],
        [ 11.26414408,  -13.8374378,   -1.04270018, 0.07151578,  0.03218943 ],
        [ 12.32835739,  -14.94509654,  -1.0390576,  0.10447862,  0.02982152 ],
        [ 13.49497055,  -16.15527241,  -1.03572513, 0.13794461,  0.02761066 ],
        [ 14.77298367,  -17.47693687,  -1.03267779, 0.17187879,  0.02554946 ],
        [ 16.17359679,  -18.92130759,  -1.0298895,  0.20628203,  0.02362854 ],
        [ 17.70840994,  -20.49997885,  -1.0273382,  0.24113748,  0.02183997 ],
        [ 19.37932623,  -22.2145782,   -1.02501743, 0.27621092,  0.02018592 ],
        [ 21.15679684,  -24.03462421,  -1.02294822, 0.31072447,  0.01868828 ],
        [ 23.01002516,  -25.928648,    -1.02112863, 0.34408882,  0.0173525 ],
        [ 25.10255178,  -28.0635202,   -1.01939438, 0.37901355,  0.01606216 ],
        [ 27.47637103,  -30.48136095,  -1.01774401, 0.41562456,  0.01481787 ],
        [ 29.94329613,  -32.99023224,  -1.01630371, 0.45078699,  0.01371822 ],
        [ 32.73841169,  -35.82894465,  -1.01493169, 0.48760888,  0.01265822 ],
        [ 35.61414061,  -38.74585617,  -1.01374283, 0.52264013,  0.01172933 ],
        [ 38.86497971,  -42.03947457,  -1.01260879, 0.55927541,  0.01083388 ],
        [ 42.16644647,  -45.38091712,  -1.01163441, 0.59372786,  0.01005681 ],
        [ 45.88510515,  -49.14106069,  -1.01070331, 0.62969538,  0.00930728 ],
        [ 50.0930693,   -53.39214075,  -1.00981495, 0.66729987,  0.00858552 ],
        [ 54.31519486,  -57.65407155,  -1.00906075, 0.70221106,  0.00796742 ],
        [ 59.06895491,  -62.44914253,  -1.00833954, 0.73863261,  0.00737155 ],
        [ 64.44611425,  -67.86924036,  -1.00765097, 0.77668638,  0.00679805 ],
        [ 69.75021324,  -73.21235327,  -1.00707495, 0.81142999,  0.00631469 ],
        [ 75.70873966,  -79.21134817,  -1.00652339, 0.8476317,   0.00584863 ],
        [ 82.43249935,  -85.97714686,  -1.00599606, 0.88540745,  0.00539997 ],
        [ 90.05666036,  -93.64504788,  -1.00549272, 0.9248877,   0.00496877 ],
        [ 97.43329385,  -101.06063806, -1.00508021, 0.96019592,  0.00461314 ],
        [ 105.72055377, -109.38831713, -1.004685,   0.99696917,  0.00427044 ],
        [ 115.07304204, -118.78280434, -1.00430695, 1.0353247,   0.00394071 ],
        [ 125.67947384, -129.43294696, -1.00394591, 1.07539459,  0.00362401 ],
        [ 135.64202024, -139.43333397, -1.00365794, 1.110203,    0.00337006 ],
        [ 146.79916044, -150.62970431, -1.00338159, 1.14640264,  0.00312521 ],
        [ 159.3480032,  -163.2192755,  -1.00311679, 1.18410224,  0.00288947 ],
        [ 173.52790204, -177.44152313, -1.00286345, 1.22342388,  0.00266287 ],
        [ 189.63178402, -193.58951201, -1.00262148, 1.2645053,   0.00244542 ],
        [ 204.14032109, -208.1347047,  -1.00243603, 1.29874169,  0.00227806 ]
    ],

    energy_table => [
        [ -4.78132852,  0.10220763, 0.05731347, 0.10220763, 0.05731347,    5.7087174,   -1.6006957 ],
        [ -4.66241093,  0.10925514, 0.06125746, 0.10925514, 0.06125746,    5.5179548,   -1.6082771 ],
        [ -4.26064413,  0.13684803, 0.07666519, 0.13684803, 0.07666519,    4.8662112,   -1.639969 ],
        [ -3.87600953,  0.16972396, 0.09489939, 0.16972396, 0.09489939,    4.228877,    -1.6803426 ],
        [ -3.48532827,  0.2110679,  0.11746359, 0.2110679,  0.11746359,    3.5631303,   -1.7350702 ],
        [ -3.22050707,  0.24447916, 0.1351847,  0.24447916, 0.1351847,     3.0980814,   -1.787313 ],
        [ -3.00775431,  0.27487188, 0.15068877, 0.27487188, 0.15068877,    2.7124883,   -1.8394949 ],
        [ -2.8150451,   0.30532164, 0.16539865, 0.30532164, 0.16539865,    2.3532712,   -1.8969463 ],
        [ -2.64540276,  0.33449055, 0.17846395, 0.33449055, 0.17846395,    2.0265543,   -1.9604925 ],
        [ -2.49490817,  0.36220046, 0.18968525, 0.36220046, 0.18968525,    1.7268946,   -2.0269408 ],
        [ -2.3538195,   0.38965962, 0.1993813,  0.38965962, 0.1993813,     1.4361842,   -2.1046419 ],
        [ -2.2242025,   0.41601656, 0.20706687, 0.41601656, 0.20706687,    1.1581281,   -2.1867551 ],
        [ -2.10192565,  0.44170067, 0.21272978, 0.44170067, 0.21272978,    0.88594662,  -2.265091 ],
        [ -1.98446373,  0.4669146,  0.21622398, 0.4669146,  0.21622398,    0.61546753,  -2.3501968 ],
        [ -1.86890539,  0.4919897,  0.21733705, 0.4919897,  0.21733705,    0.33848252,  -2.4537746 ],
        [ -1.75099557,  0.51754744, 0.21567512, 0.51754744, 0.21567512,    0.042319224, -2.5426787 ],
        [ -1.62099243,  0.54527302, 0.2102071,  0.54527302, 0.2102071,     -0.29266619, -2.6118927 ],
        [ -1.49813463,  0.57059485, 0.20142994, 0.57059485, 0.20142994,    -0.617638,   -2.6077053 ],
        [ -1.38562252,  0.59266749, 0.1904986,  0.59266749, 0.1904986,     -0.90718283, -2.4647659 ],
        [ -1.28014194,  0.61212443, 0.17811732, 0.61212443, 0.17811732,    -1.1564197,  -2.2304979 ],
        [ -1.17945025,  0.62940105, 0.16484943, 0.62940105, 0.16484943,    -1.3682887,  -1.9621818 ],
        [ -1.08169596,  0.6448496,  0.15112406, 0.6448496,  0.15112406,    -1.5466286,  -1.757331 ],
        [ -0.98494944,  0.65879843, 0.13722134, 0.65879843, 0.13722134,    -1.7102254,  -1.6681302 ],
        [ -0.8885209,   0.67136636, 0.12350479, 0.67136636, 0.12350479,    -1.8688842,  -1.6531933 ],
        [ -0.79143846,  0.68270619, 0.11022358, 0.68270619, 0.11022358,    -2.0301468,  -1.6798374 ],
        [ -0.69281126,  0.69294552, 0.0975747,  0.69294552, 0.0975747,     -2.197703,   -1.7209842 ],
        [ -0.59180969,  0.70219163, 0.08570894, 0.70219163, 0.08570894,    -2.3738109,  -1.7655641 ],
        [ -0.4875348,   0.71054479, 0.07472347, 0.71054479, 0.07472347,    -2.5602777,  -1.8083799 ],
        [ -0.37925484,  0.71808033, 0.06469347, 0.71808033, 0.06469347,    -2.7583554,  -1.8469743 ],
        [ -0.26608706,  0.72487617, 0.05564352, 0.72487617, 0.05564352,    -2.9694626,  -1.8804305 ],
        [ -0.14722441,  0.73099687, 0.04757555, 0.73099687, 0.04757555,    -3.1948476,  -1.9085908 ],
        [ -0.02158831,  0.73651255, 0.04045088, 0.73651255, 0.04045088,    -3.4362836,  -1.9317486 ],
        [ 0.11177408,   0.74147805, 0.03422446, 0.74147805, 0.03422446,    -3.6953275,  -1.9503498 ],
        [ 0.25409448,   0.74595113, 0.0288277,  0.74595113, 0.0288277,     -3.9741086,  -1.9649656 ],
        [ 0.40658824,   0.7499802,  0.02418962, 0.7499802,  0.02418962,    -4.2747565,  -1.9761641 ],
        [ 0.56746807,   0.75354688, 0.02030094, 0.75354688, 0.02030094,    -4.5934658,  -1.9843872 ],
        [ 0.73480193,   0.75666555, 0.01709869, 0.75666555, 0.01709869,    -4.9261036,  -1.9901908 ],
        [ 0.90970709,   0.75941559, 0.01445091, 0.75941559, 0.01445091,    -5.2746227,  -1.9942192 ],
        [ 1.09380544,   0.76186524, 0.01224787, 0.76186524, 0.01224787,    -5.6420656,  -1.9969205 ],
        [ 1.28879861,   0.76406676, 0.01040519, 0.76406676, 0.01040519,    -6.0316612,  -1.9986469 ],
        [ 1.49524617,   0.76604968, 0.00886545, 0.76604968, 0.00886545,    -6.4444192,  -1.999721 ],
        [ 1.71424664,   0.76784447, 0.00757605, 0.76784447, 0.00757605,    -6.8824485,  -2.0003293 ],
        [ 1.94781372,   0.76948219, 0.00649021, 0.76948219, 0.00649021,    -7.3497093,  -2.0006111 ],
        [ 2.19679148,   0.77097982, 0.00557587, 0.77097982, 0.00557587,    -7.8478352,  -2.0006933 ],
        [ 2.46345995,   0.77235948, 0.00480176, 0.77235948, 0.00480176,    -8.3813596,  -2.0006709 ],
        [ 2.74913587,   0.77363387, 0.00414574, 0.77363387, 0.00414574,    -8.9528934,  -2.0005951 ],
        [ 3.05566257,   0.7748159,  0.00358829, 0.7748159,  0.00358829,    -9.5661154,  -2.0004876 ],
        [ 3.38510706,   0.77591677, 0.00311318, 0.77591677, 0.00311318,    -10.225143,  -2.0003725 ],
        [ 3.73992776,   0.77694654, 0.00270683, 0.77694654, 0.00270683,    -10.934898,  -2.0002777 ],
        [ 4.1263777,    0.77792205, 0.00235525, 0.77792205, 0.00235525,    -11.707887,  -2.0001939 ],
        [ 4.5284835,    0.77880826, 0.00206331, 0.77880826, 0.00206331,    -12.512161,  -2.0001313 ],
        [ 4.97490057,   0.7796693,  0.00180397, 0.7796693,  0.00180397,    -13.405042,  -2.0000852 ],
        [ 5.45864664,   0.78048563, 0.00157939, 0.78048563, 0.00157939,    -14.372565,  -2.0000525 ],
        [ 5.98374325,   0.78126188, 0.00138454, 0.78126188, 0.00138454,    -15.422779,  -2.0000311 ],
        [ 6.55460664,   0.7820019,  0.00121537, 0.7820019,  0.00121537,    -16.564519,  -2.0000174 ],
        [ 7.17644864,   0.78270948, 0.00106633, 0.78270948, 0.00106633,    -17.80821,   -2.0000094 ],
        [ 7.85467748,   0.78338703, 0.00093641, 0.78338973, 0.00095926135, -19.164672,  -2.0000047 ],
        [ 8.59569863,   0.78403717, 0.00082247, 0.78408567, 0.00092666377, -20.646717,  -2.0000023 ],
        [ 9.40571556,   0.78466139, 0.00072244, 0.78482595, 0.0009007123,  -22.266752,  -2.0000011 ],
        [ 10.29253033,  0.78526164, 0.00063448, 0.78561604, 0.0008771719,  -24.040382,  -2.0000005 ],
        [ 11.26414408,  0.78583913, 0.00055708, 0.78644849, 0.00083985609, -25.98361,   -2.0000002 ],
        [ 12.32835739,  0.78639443, 0.000489,   0.7873163,  0.00079385192, -28.112037,  -2.0000001 ],
        [ 13.49497055,  0.78692868, 0.0004291,  0.78826687, 0.00088388458, -30.445263,  -2.0000001 ],
        [ 14.77298367,  0.78744219, 0.00037643, 0.78932363, 0.00077395708, -33.001289,  -2.0000001 ],
        [ 16.17359679,  0.78793579, 0.00033012, 0.79033766, 0.000677609,   -35.802516,  -2.0000001 ],
        [ 17.70840994,  0.78841008, 0.00028942, 0.79131046, 0.00059317535, -38.872142,  -2.0000001 ],
        [ 19.37932623,  0.7888629,  0.00025388, 0.79223787, 0.00051961845, -42.213975,  -2 ],
        [ 21.15679684,  0.78928613, 0.00022344, 0.79310361, 0.0004567732,  -45.768916,  -2 ],
        [ 23.01002516,  0.78967554, 0.0001977,  0.79389929, 0.00040374365, -49.475373,  -2 ],
        [ 25.10255178,  0.79006367, 0.00017412, 0.79469157, 0.00035523832, -53.660426,  -2 ],
        [ 27.47637103,  0.79045047, 0.00015258, 0.7954804,  0.00031101217, -58.408064,  -2 ],
        [ 29.94329613,  0.79080382, 0.00013454, 0.79620042, 0.00027402807, -63.341915,  -2 ],
        [ 32.73841169,  0.79115595, 0.00011804, 0.79691744, 0.00024026124, -68.932146,  -2 ],
        [ 35.61414061,  0.79147499, 0.00010432, 0.79756663, 0.00021220985, -74.683604,  -2 ],
        [ 38.86497971,  0.79179296, 9.176e-05,  0.79821328, 0.00018654515, -81.185282,  -2 ],
        [ 42.16644647,  0.7920782,  8.139e-05,  0.79879305, 0.00016538406, -87.788215,  -2 ],
        [ 45.88510515,  0.79236252, 7.186e-05,  0.7993707,  0.00014596567, -95.225533,  -2 ],
        [ 50.0930693,   0.7926459,  6.315e-05,  0.79994619, 0.00012820227, -103.64146,  -2 ],
        [ 54.31519486,  0.79289699, 5.604e-05,  0.80045591, 0.0001137325,  -112.08571,  -2 ],
        [ 59.06895491,  0.79314731, 4.951e-05,  0.80096389, 0.0001004434,  -121.59323,  -2 ],
        [ 64.44611425,  0.79339683, 4.352e-05,  0.80147011, 8.8277063e-05, -132.34755,  -2 ],
        [ 69.75021324,  0.7936145,  3.872e-05,  0.80191158, 7.8507513e-05, -142.95575,  -2 ],
        [ 75.70873966,  0.79383155, 3.429e-05,  0.80235168, 6.9516188e-05, -154.8728,   -2 ],
        [ 82.43249935,  0.79404795, 3.023e-05,  0.80279038, 6.1266246e-05, -168.32032,  -2 ],
        [ 90.05666036,  0.7942637,  2.651e-05,  0.80322767, 5.3721325e-05, -183.56864,  -2 ],
        [ 97.43329385,  0.79444811, 2.359e-05,  0.80360136, 4.7788373e-05, -198.32191,  -2 ],
        [ 105.72055377, 0.79463202, 2.089e-05,  0.80397399, 4.2324941e-05, -214.89643,  -2 ],
        [ 115.07304204, 0.79481544, 1.842e-05,  0.80434556, 3.7309265e-05, -233.60141,  -2 ],
        [ 125.67947384, 0.79499834, 1.616e-05,  0.80471605, 3.2719802e-05, -254.81427,  -2 ],
        [ 135.64202024, 0.79515038, 1.442e-05,  0.80502396, 2.9205598e-05, -274.73936,  -2 ],
        [ 146.79916044, 0.79530205, 1.282e-05,  0.80533112, 2.5960499e-05, -297.05364,  -2 ],
        [ 159.3480032,  0.79545336, 1.135e-05,  0.80563751, 2.2972671e-05, -322.15133,  -2 ],
        [ 173.52790204, 0.7956043,  9.99e-06,   0.80594313, 2.0230272e-05, -350.51113,  -2 ],
        [ 189.63178402, 0.79575486, 8.76e-06,   0.80624797, 1.7721614e-05, -382.71889,  -2 ],
        [ 204.14032109, 0.79587505, 7.84e-06,   0.80649129, 1.587523e-05,  -411.73596,  -2 ]
    ],

    impulse_table => [
        [
            -5.3862182, 13.81484,    -5.3866386,     0.016482735, -0.0004434523, -0.14738127,
            -1.2153852, -0.37925003, -0.00093968831, -0.32786389, 0.15506093,    0.91045348,
            -1.805932
        ],
        [
            -4.9618452, 12.541724,   -4.9626398,    0.019700671, -0.00067787526, -0.14496053,
            -1.2129644, -0.59143997, -0.0014364374, -0.32544315, 0.1550606,      0.88638865,
            -1.7381246
        ],
        [
            -4.6051702, 11.471706,   -4.6065272,    0.022719096, -0.00096839323, -0.14196053,
            -1.2099644, -0.76978467, -0.0020520534, -0.32244315, 0.1550597,      0.86123386,
            -1.6849067
        ],
        [
            -4.1997051, 10.255337,   -4.2021995,    0.026430575, -0.0014525898, -0.13696053,
            -1.2049644, -0.97254197, -0.0030780801, -0.31744315, 0.15505668,    0.82583194,
            -1.630352
        ],
        [
            -3.9120231, 9.3923413,  -3.9158661,    0.029162531, -0.0019367865, -0.13196053,
            -1.1999644, -1.1164309, -0.0041041067, -0.31244315, 0.15505085,    0.79540127,
            -1.5969144
        ],
        [
            -3.506558,  8.1761552,  -3.5136289,    0.032932296, -0.0029051797, -0.12196053,
            -1.1899644, -1.3193604, -0.0061561601, -0.30244315, 0.15502695,    0.74343334,
            -1.5602909
        ],
        [
            -2.9957323, 6.6447543,  -3.0110015,   0.036841408, -0.0048419662, -0.10196053,
            -1.1699644, -1.5757891, -0.010260267, -0.28244315, 0.15490389,    0.65965508,
            -1.5404116
        ],
        [
            -2.3025851, 4.5747501,  -2.3461933,   0.037734023, -0.0096839323, -0.051960537,
            -1.1199644, -1.9314457, -0.020520534, -0.23244316, 0.15382079,    0.50787647,
            -1.5973513
        ],
        [
            -1.89712,   3.3825737,  -1.9777805,   0.033834977, -0.01452589, -0.0019609241,
            -1.0699647, -2.1589716, -0.030780789, -0.18244316, 0.15101445,  0.40177481,
            -1.7109299
        ],
        [
            -1.609438,  2.5618751,  -1.7336979,   0.028966455, -0.019366967, 0.048010267,
            -1.0199813, -2.3519346, -0.041039915, -0.13256944, 0.14603295,   0.32494177,
            -1.8480149
        ],
        [
            -1.3862944,  1.9520781,  -1.558717,    0.025270589,  -0.024178107, 0.097223198,
            -0.97044084, -2.5481861, -0.051260218, -0.082822117, 0.13891853,   0.2687805,
            -1.994057
        ],
        [
            -1.2039729,  1.4783455,  -1.4272843,   0.0240345,    -0.02856695, 0.13869531,
            -0.92606869, -2.7772256, -0.060865721, -0.037013294, 0.13014457,  0.2273547,
            -2.1400698
        ],
        [
            -1.0498222,  1.0982683,  -1.3252223,   0.024696895,   -0.031270749, 0.1617946,
            -0.89801154, -3.0576357, -0.067171033, -0.0090046068, 0.12039944,   0.19639252,
            -2.280657
        ],
        [
            -0.91629079, 0.78530282, -1.2438401,   0.025844908,    -0.032345339, 0.17363415,
            -0.88266186, -3.3612127, -0.069432738, -8.5261029e-05, 0.11037277,   0.17286488,
            -2.4131084
        ],
        [
            -0.79850775, 0.52193204, -1.1774893,  0.026908861,  -0.032788264, 0.18111612,
            -0.87276404, -3.6541279, -0.07006203, 0.0011857099, 0.10062858,   0.15465863,
            -2.5364851
        ],
        [
            -0.69314724, 0.29618966, -1.122353,    0.027780549,   -0.033015338, 0.18658334,
            -0.86554747, -3.9259279, -0.070192544, 0.00017736968, 0.091562372,  0.14030903,
            -2.6508484
        ],
        [
            -0.59783706, 0.099669321, -1.0757742,   0.028471646,   -0.03315434, 0.19089241,
            -0.85993264, -4.1757961,  -0.070157053, -0.0018228397, 0.083408659, 0.1287979,
            -2.7567424
        ],
        [
            -0.51082568, -0.073682254, -1.0358574,   0.029016464,   -0.033250127, 0.19443326,
            -0.85539773, -4.4054616,   -0.070062051, -0.0040964582, 0.076270276,  0.11941085,
            -2.854898
        ],
        [
            -0.43078297, -0.22832425, -1.0012201,   0.029447761,   -0.033321385, 0.1974186,
            -0.85164631, -4.6171531,  -0.069946545, -0.0065413209, 0.070152887,  0.11164013,
            -2.9460807
        ],
        [
            -0.356675,   -0.36761016, -0.97083413,  0.029791927,   -0.033377152, 0.19998051,
            -0.84848964, -4.813015,   -0.069826336, -0.0092473887, 0.064996911,  0.10511939,
            -3.0310191
        ],
        [
            -0.28768213, -0.49411004, -0.94392096,  0.030069112,  -0.033422339, 0.20220809,
            -0.8457994,  -4.9949525,  -0.069708209, -0.011858505, 0.060697929,  0.099580221,
            -3.1103765
        ],
        [
            -0.22000858, -0.61538698, -0.91874044,  0.030304661,  -0.03346162, 0.20425754,
            -0.843375,   -5.1728374,  -0.069589527, -0.014396757, 0.056923436, 0.094601591,
            -3.1883564
        ],
        [
            -0.10536057, -0.81494998, -0.878647,    0.030632515,  -0.033519081, 0.20744881,
            -0.83970995, -5.4725215,  -0.069387543, -0.019326506, 0.051376874,  0.08708538,
            -3.3205097
        ],
        [
            -0.05129335, -0.9066641, -0.86078428,  0.030760278,  -0.033542911, 0.20884027,
            -0.83816008, -5.6130048, -0.069293371, -0.021720387, 0.04907973,   0.083899948,
            -3.3827609
        ],
        [
            -5.5584269e-08, -0.9923442, -0.84441752,  0.030867716,  -0.033563877, 0.21009752,
            -0.83678763,    -5.7457337, -0.069205297, -0.024055852, 0.047064978,  0.081069621,
            -3.4417353
        ],
        [
            0.095310124, -1.1483204, -0.8154154,   0.031036245,  -0.033599126, 0.21228055,
            -0.83447406, -5.9908768, -0.069045826, -0.028572392, 0.043698298,  0.076260959,
            -3.5510288
        ],
        [
            0.1823215,   -1.2873058, -0.7904278,  0.031160074,  -0.033627683, 0.21411116,
            -0.83260976, -6.2129252, -0.06890598, -0.032646722, 0.040997484,  0.072328382,
            -3.6504001
        ],
        [
            0.26236421,  -1.4125119, -0.76859871,  0.03125303,   -0.033651344, 0.21566843,
            -0.83108606, -6.4156835, -0.068782754, -0.036837179, 0.03878238,   0.069050478,
            -3.7414163
        ],
        [
            0.40546505,  -1.6306169, -0.73208205, 0.031379318,  -0.033688408, 0.21817595,
            -0.82876662, -6.7745656, -0.06857642, -0.044037651, 0.035362058,  0.063890133,
            -3.9030897
        ],
        [
            0.5306282,   -1.8159654, -0.70251264,  0.031457453,  -0.033716236, 0.22010668,
            -0.82711099, -7.0847191, -0.068411141, -0.050953806, 0.032839077,  0.060001101,
            -4.0433348
        ],
        [
            0.69314712,  -2.0500229, -0.66701152,  0.031526611,  -0.033747194, 0.22229223,
            -0.82539514, -7.4824152, -0.068217556, -0.059870864, 0.030081513,  0.055663168,
            -4.2237871
        ],
        [
            0.99325172,  -2.4657349, -0.60868568,  0.031593318,  -0.033792459, 0.22552851,
            -0.82322087, -8.2029341, -0.067919103, -0.076764878, 0.026138284,  0.049283324,
            -4.5522756
        ],
        [
            1.0986122,   -2.6073819, -0.59009716,  0.031604971, -0.033805405, 0.22645828,
            -0.82269043, -8.4519834, -0.067830721, -0.08305888, 0.025021965,  0.047435159,
            -4.6662401
        ],
        [
            1.2089603,   -2.7536664, -0.57153855,  0.031613053,  -0.033817649, 0.22733481,
            -0.82223383, -8.7108285, -0.067746599, -0.089470686, 0.023969965,  0.045674817,
            -4.7849027
        ],
        [
            1.3862943,   -2.9847633, -0.54347281,  0.031619727,  -0.033834805, 0.22855898,
            -0.82167138, -9.1228255, -0.067627535, -0.099918453, 0.022491079,  0.043167502,
            -4.9742123
        ],
        [
            1.6094379,   -3.2694424, -0.51086248,  0.031621235, -0.033852713, 0.22982512,
            -0.82119149, -9.6349236, -0.067502852, -0.11296973, 0.020930302,  0.04047664,
            -5.2102395
        ],
        [
            1.7917594,   -3.497675,  -0.48615841,  0.03161921, -0.033864709, 0.23067084,
            -0.82093027, -10.048667, -0.067418413, -0.1241706, 0.01985178,   0.038588091,
            -5.401501
        ],
        [
            1.9459101,   -3.687969,  -0.46646561,  0.031616357, -0.033873347, 0.23127612,
            -0.82077693, -10.395549, -0.067358047, -0.13293724, 0.019051333,  0.037169922,
            -5.5622373
        ],
        [
            2.0794415,   -3.8510361, -0.45020261,  0.031613559, -0.03387995, 0.23173112,
            -0.82068368, -10.694058, -0.067312139, -0.14044324, 0.018427495, 0.036054353,
            -5.700832
        ],
        [
            2.3027866,  -4.1204992, -0.42448291,  0.031608489, -0.033889409, 0.23236912,
            -0.8205848, -11.189643, -0.067247692, -0.15443794, 0.017505039,  0.034056244,
            -5.9397848
        ],
        [
            2.9958726,   -4.9356151, -0.3543065,   0.031594842, -0.0339083, 0.23364861,
            -0.82045307, -12.703498, -0.067117809, -0.19632687, 0.01533919, 0.030255061,
            -6.643809
        ],
        [
            3.9121518,   -5.9784311, -0.27784469,  0.031584572, -0.033924986, 0.23443007,
            -0.82080545, -14.665243, -0.067044023, -0.25188798, 0.013450419,  0.026748182,
            -7.5679257
        ],
        [
            4.6055286,   -6.749806,  -0.22863282,  0.031580871, -0.033925684, 0.23467972,
            -0.82061037, -16.130052, -0.067014858, -0.29398137, 0.012432586,  0.024792839,
            -8.2643123
        ],
        [
            5.2984059,   -7.5098204, -0.18493674,  0.031578994, -0.03392759, 0.23480792,
            -0.82061893, -17.582252, -0.067001716, -0.33605671, 0.011630005, 0.023225261,
            -8.9588394
        ],
        [
            6.2147728,   -8.5026576, -0.13363973,  0.031577874, -0.033928729, 0.23488491,
            -0.82062439, -19.490133, -0.066993822, -0.39170738, 0.010786374,  0.021559394,
            -9.8762856
        ],
        [
            6.907833,    -9.2464094, -0.098746731, 0.031577506, -0.033929107, 0.23491057,
            -0.82062627, -20.925908, -0.066991192, -0.43379642, 0.010263366,  0.020520209,
            -10.569735
        ],
        [
            7.6010681,   -9.9855036, -0.06658251,  0.031577325, -0.033929295, 0.2349234,
            -0.82062722, -22.357269, -0.066989877, -0.47589532, 0.0098125988, 0.019622,
            -11.263177
        ],
        [
            8.5172967,   -10.956326, -0.027572224, 0.031577218, -0.033929406, 0.2349311,
            -0.82062779, -24.243194, -0.066989088, -0.53153521, 0.0093016065, 0.018601961,
            -12.179536
        ],
        [
            9.2105119,   -11.687094, -0.00029799171, 0.031577182, -0.03392944,  0.23493367,
            -0.82062798, -25.666443, -0.066988825,   -0.57363163, 0.0089653252, 0.017930033,
            -12.872797
        ],
        [
            10.819938,   -13.373885, 0.056973772,  0.031577151, -0.03392947,  0.23493572,
            -0.82062813, -28.961291, -0.066971779, -0.66591352, 0.0083092187, 0.016618317,
            -14.482263
        ],
        [
            11.512947,   -14.096751, 0.079451096,  0.031577146, -0.033929477, 0.23493598,
            -0.82062815, -30.376694, -0.066679075, -0.68893235, 0.0080686017, 0.016137144,
            -15.175277
        ],
        [
            13.122366,   -15.769173, 0.12753354,   0.03157714,  -0.033929471, 0.23493619,
            -0.82062817, -33.65764,  -0.065185992, -0.73989028, 0.0075827264, 0.015165441,
            -16.784701
        ],
        [
            13.815574,   -16.487198, 0.14670822,   0.031577138, -0.033929488, 0.2349362,
            -0.82062818, -35.068555, -0.064369977, -0.76089412, 0.0073992057, 0.014798406,
            -17.477909
        ],
        [
            16.118196,   -18.864248, 0.20497106,   0.031577131, -0.033929033, 0.23493628,
            -0.82724611, -39.747385, -0.061448402, -0.82724611, 0.0068742834, 0.013748566,
            -19.780532
        ],
        [
            18.420816,  -21.231485, 0.25643059,   0.031577125, -0.03389039,  0.23493651,
            -0.8892132, -44.416625, -0.058570572, -0.8892132,  0.0064481029, 0.012896206,
            -22.083152
        ]
    ],

    tail_shock_table => [
        [ 7.614977,  -0.5482,     -0.035102082, -0.035102082 ],
        [ 8.5171022, -0.58355341, -0.050927729, -0.017370808 ],
        [ 9.2104119, -0.60942456, -0.053260975, -0.012354726 ],
        [ 10.819917, -0.66591234, -0.054448639, -0.0060627521 ],
        [ 11.512937, -0.68893176, -0.054210655, -0.0044678173 ],
        [ 13.122364, -0.73989016, -0.052996748, -0.0020389706 ],
        [ 13.815573, -0.76089406, -0.052333317, -0.0013457077 ]
    ],

    blast_info => [
        0.15196047, 1.2199638, 0.33250735, -0.20520545, 1.4860429,  0.0031118032, 3.4734582,   -0.096839318,
        0.40336997, 0.4179451, 0.5496372,  -0.14009504, 0.41389801, 0.15506091,   0.025672464, -0.027584928,
        2027.8,     -0.5482,   40426.068,  -0.6584485
    ],

};
1;
