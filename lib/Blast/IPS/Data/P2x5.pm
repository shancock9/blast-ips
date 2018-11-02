package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='P', gamma=2.5
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'P2.5'} = {

    shock_table_info => [ 0, 2.5, 2.07e-06, 4000, 1.72e-06, 100, 1.2e-06, 3.7e-07, 5e-07 ],

    shock_table => [
        [ -11.13471153, 12.00046197, -0.99999123, -11.13668814, 0.99901073 ],
        [ -9.45523631,  10.32102358, -0.99995471, -9.45981953,  0.99770326 ],
        [ -8.26435842,  9.13025006,  -0.99985003, -8.27268662,  0.99581928 ],
        [ -7.31652958,  8.18265795,  -0.99961379, -7.32993895,  0.99325334 ],
        [ -6.55551737,  7.42208596,  -0.99917434, -6.5751917,   0.99007544 ],
        [ -5.91196846,  6.77928197,  -0.99843207, -5.93920062,  0.98622325 ],
        [ -5.35526741,  6.22374689,  -0.99727341, -5.3913731,   0.98167854 ],
        [ -4.86229645,  5.73250887,  -0.99555849, -4.90868167,  0.97639007 ],
        [ -4.41868506,  5.29136843,  -0.99312802, -4.4768398,   0.97031262 ],
        [ -4.01256196,  4.8886711,   -0.98978676, -4.08413394,  0.96336759 ],
        [ -3.63517988,  4.51594009,  -0.98530254, -3.7220199,   0.95546018 ],
        [ -3.27904117,  4.16603295,  -0.97938675, -3.38329837,  0.94645776 ],
        [ -2.93664815,  3.83195695,  -0.97165157, -3.06095257,  0.93615223 ],
        [ -2.60050319,  3.50695985,  -0.96156139, -2.74821937,  0.92423798 ],
        [ -2.25822434,  3.1800325,   -0.94816377, -2.43423497,  0.91007661 ],
        [ -1.88021745,  2.82507943,  -0.92905916, -2.09357138,  0.89188935 ],
        [ -1.52508965,  2.49900388,  -0.90657713, -1.78023929,  0.87233477 ],
        [ -1.19959956,  2.20780906,  -0.88209469, -1.49949258,  0.85244073 ],
        [ -0.89379815,  1.94197731,  -0.85604896, -1.24186515,  0.83228176 ],
        [ -0.61129814,  1.70379672,  -0.82990592, -1.00949393,  0.81269858 ],
        [ -0.33793119,  1.48053767,  -0.80334642, -0.78998598,  0.79319746 ],
        [ -0.06444353,  1.26453608,  -0.77622864, -0.57575101,  0.77348565 ],
        [ 0.21680687,   1.05013322,  -0.74849617, -0.36104595,  0.75335626 ],
        [ 0.50541293,   0.83812374,  -0.72089749, -0.14654802,  0.73319453 ],
        [ 0.80248763,   0.62800281,  -0.69399937, 0.06828242,   0.71327458 ],
        [ 1.111108,     0.41786996,  -0.66813875, 0.28536503,   0.69373163 ],
        [ 1.43565355,   0.20509639,  -0.64352535, 0.50737128,   0.67463244 ],
        [ 1.78063403,   -0.01282444, -0.62037056, 0.73685203,   0.656068 ],
        [ 2.15191851,   -0.23905328, -0.59882327, 0.97704048,   0.6380997 ],
        [ 2.55625544,   -0.47705687, -0.57903672, 1.23147463,   0.62080828 ],
        [ 2.95984614,   -0.70734431, -0.56268262, 1.47891559,   0.60573544 ],
        [ 3.36811553,   -0.93420647, -0.54911178, 1.7234494,    0.59248342 ],
        [ 3.78795115,   -1.16228079, -0.53778724, 1.96966027,   0.5806989 ],
        [ 4.22317539,   -1.39421418, -0.528382,   2.22004973,   0.57019277 ],
        [ 4.6785161,    -1.63296767, -0.52061564, 2.4774853,    0.56079578 ],
        [ 5.15710146,   -1.88054519, -0.51428472, 2.74381064,   0.55240594 ],
        [ 5.66370398,   -2.13973205, -0.50919383, 3.0217083,    0.5449145 ],
        [ 6.20716835,   -2.41530483, -0.50515339, 3.31596773,   0.53819374 ],
        [ 6.78345415,   -2.70548054, -0.50208034, 3.62436476,   0.53228697 ],
        [ 7.40409313,   -3.01633572, -0.49979622, 3.95304747,   0.52706501 ],
        [ 8.0894025,    -3.35824425, -0.49815455, 4.31259247,   0.52239902 ],
        [ 8.83749753,   -3.73047275, -0.49708519, 4.70182074,   0.51834589 ],
        [ 9.64021832,   -4.12922353, -0.49648847, 5.1164817,    0.5149344 ],
        [ 10.58748417,  -4.59938285, -0.49624033, 5.60273255,   0.51185348 ],
        [ 11.65724245,  -5.13024871, -0.49629981, 6.14884075,   0.50927421 ],
        [ 12.91603584,  -5.75516081, -0.49660205, 6.78845806,   0.50709048 ],
        [ 14.41868094,  -6.50173444, -0.49708715, 7.54898636,   0.50527647 ],
        [ 16.34782696,  -7.46130512, -0.4977176,  8.52213089,   0.5037209 ],
        [ 18.94944892,  -8.75714151, -0.49842674, 9.83076142,   0.50239997 ],
        [ 18.96705075,  -8.7659149,  -0.49847801, 9.83960452,   0.50239296 ]
    ],

    energy_table => [
        [ -11.13471153, 0.00083502, 0.00049713, 0.00083502, 0.00049713, 11.628379,    -0.96013024 ],
        [ -9.45523631,  0.00226188, 0.00133634, 0.00226188, 0.00133634, 9.9824963,    -0.9969959 ],
        [ -8.26435842,  0.00455684, 0.00266598, 0.00455684, 0.00266598, 8.7808427,    -1.0183936 ],
        [ -7.31652958,  0.00790898, 0.0045707,  0.00790898, 0.0045707,  7.8085301,    -1.0319662 ],
        [ -6.55551737,  0.01224165, 0.00697417, 0.01224165, 0.00697417, 7.0194428,    -1.0410065 ],
        [ -5.91196846,  0.01761249, 0.00986928, 0.01761249, 0.00986928, 6.3472649,    -1.0473679 ],
        [ -5.35526741,  0.02399306, 0.01319376, 0.02399306, 0.01319376, 5.7628062,    -1.0518681 ],
        [ -4.86229645,  0.03137748, 0.01689051, 0.03137748, 0.01689051, 5.2433896,    -1.0550033 ],
        [ -4.41868506,  0.03973165, 0.02088104, 0.03973165, 0.02088104, 4.7748361,    -1.057071 ],
        [ -4.01256196,  0.04904957, 0.02509309, 0.04904957, 0.02509309, 4.3452207,    -1.0582718 ],
        [ -3.63517988,  0.05932847, 0.02944623, 0.05932847, 0.02944623, 3.9456982,    -1.0587352 ],
        [ -3.27904117,  0.07059424, 0.03385999, 0.07059424, 0.03385999, 3.5686192,    -1.0585416 ],
        [ -2.93664815,  0.08293879, 0.03826048, 0.08293879, 0.03826048, 3.2062664,    -1.0577276 ],
        [ -2.60050319,  0.09652673, 0.0425671,  0.09652673, 0.0425671,  2.8509038,    -1.0562805 ],
        [ -2.25822434,  0.1118183,  0.04672673, 0.1118183,  0.04672673, 2.4896716,    -1.0540915 ],
        [ -1.88021745,  0.13026695, 0.05075745, 0.13026695, 0.05075745, 2.0917531,    -1.0508818 ],
        [ -1.52508965,  0.14884439, 0.05370806, 0.14884439, 0.05370806, 1.7191541,    -1.0461239 ],
        [ -1.19959956,  0.1666413,  0.05548511, 0.1666413,  0.05548511, 1.379568,     -1.0379367 ],
        [ -0.89379815,  0.18374629, 0.05623023, 0.18374629, 0.05623023, 1.0637084,    -1.0263005 ],
        [ -0.61129814,  0.19963028, 0.05609112, 0.19963028, 0.05609112, 0.77549849,   -1.0134044 ],
        [ -0.33793119,  0.21486127, 0.05522636, 0.21486127, 0.05522636, 0.50026805,   -1.0002466 ],
        [ -0.06444353,  0.2297716,  0.05371126, 0.2297716,  0.05371126, 0.22851059,   -0.98774866 ],
        [ 0.21680687,   0.24459037, 0.05157958, 0.24459037, 0.05157958, -0.047580432, -0.97666455 ],
        [ 0.50541293,   0.25910251, 0.04891941, 0.25910251, 0.04891941, -0.32797363,  -0.96777625 ],
        [ 0.80248763,   0.27318323, 0.04582923, 0.27318323, 0.04582923, -0.61432363,  -0.96142185 ],
        [ 1.111108,     0.28680106, 0.04239457, 0.28680106, 0.04239457, -0.91024379,  -0.95758535 ],
        [ 1.43565355,   0.29995922, 0.03868721, 0.29995922, 0.03868721, -1.2205928,   -0.95605381 ],
        [ 1.78063403,   0.31262895, 0.03478079, 0.31262895, 0.03478079, -1.5503381,   -0.95653822 ],
        [ 2.15191851,   0.32478584, 0.03074197, 0.32478584, 0.03074197, -1.9057658,   -0.95875258 ],
        [ 2.55625544,   0.33637531, 0.02664215, 0.33637531, 0.02664215, -2.294068,    -0.96248642 ],
        [ 2.95984614,   0.34636407, 0.02292235, 0.34636407, 0.02292235, -2.683382,    -0.96732189 ],
        [ 3.36811553,   0.35502286, 0.01956276, 0.35502286, 0.01956276, -3.0794234,   -0.9735988 ],
        [ 3.78795115,   0.36258432, 0.01652789, 0.36258432, 0.01652789, -3.4897073,   -0.9905111 ],
        [ 4.22317539,   0.36917056, 0.01380758, 0.36917056, 0.01380758, -3.9267854,   -1.0029089 ],
        [ 4.6785161,    0.37489078, 0.01138593, 0.37489078, 0.01138593, -4.3828072,   -1 ],
        [ 5.15710146,   0.37981441, 0.00925615, 0.37981441, 0.00925615, -4.8613926,   -1 ],
        [ 5.66370398,   0.38401812, 0.00740302, 0.38401812, 0.00740302, -5.3679951,   -1 ],
        [ 6.20716835,   0.3875898,  0.00580182, 0.3875898,  0.00580182, -5.9114595,   -1 ],
        [ 6.78345415,   0.39053164, 0.0044632,  0.39053164, 0.0044632,  -6.4877453,   -1 ],
        [ 7.40409313,   0.39294122, 0.00335228, 0.39294122, 0.00335228, -7.1083843,   -1 ],
        [ 8.0894025,    0.39490804, 0.00243481, 0.39490804, 0.00243481, -7.7936936,   -1 ],
        [ 8.83749753,   0.39644338, 0.00171114, 0.39644338, 0.00171114, -8.5417887,   -1 ],
        [ 9.64021832,   0.39758541, 0.00116799, 0.39758541, 0.00116799, -9.3445094,   -1 ],
        [ 10.58748417,  0.39847487, 0.00074159, 0.39847487, 0.00074159, -10.291775,   -1 ],
        [ 11.65724245,  0.3990946,  0.00044237, 0.3990946,  0.00044237, -11.361534,   -1 ],
        [ 12.91603584,  0.39951124, 0.00023992, 0.39951124, 0.00023992, -12.620327,   -1 ],
        [ 14.41868094,  0.39976668, 0.00011509, 0.39976668, 0.00011509, -14.122972,   -1 ],
        [ 16.34782696,  0.39991017, 4.459e-05,  0.39991017, 4.459e-05,  -16.052118,   -1 ],
        [ 18.94944892,  0.3999755,  1.233e-05,  0.3999755,  1.233e-05,  -18.65374,    -1 ],
        [ 18.96705075,  0.39997572, 1.222e-05,  0.39997572, 1.222e-05,  -18.671342,   -1 ]
    ],

    impulse_table => [
        [
            -12.949144, 13.814887, -12.949942,    1.4390335,  0,          0,
            0,          -7.611179, 0.00031032806, -88.627795, 0.28525783, 0.99929513,
            -7.7838932
        ],
        [
            -2.3025851, 3.2221368, -2.4746501,    1.2634678,  0,          0,
            0,          -2.427266, 0.00031032824, -88.527798, 0.27224838, 0.72560892,
            -1.8042201
        ],
        [
            -1.6094379, 2.5757162,  -1.854025,     1.2191109,  0,          0,
            0,          -2.2031354, 0.00031032917, -88.427798, 0.26098807, 0.63915218,
            -1.5303587
        ],
        [
            -1.3862943, 2.3738617,  -1.6597333,    1.2064848,  0,         0,
            0,          -2.1475547, 0.00031033004, -88.377798, 0.2558881, 0.60909739,
            -1.4503745
        ],
        [
            -1.2039728,   2.2116674,  -1.5032211, 1.1973817,  0, 0, 0, -2.1093165,
            0.0003103312, -88.327798, 0.25109118, 0.58400455, -1.3881757
        ],
        [
            -0.9162907, 1.9612546,  -1.2606025,    1.1859216,  0,          0,
            0,          -2.0631676, 0.00031033448, -88.227798, 0.24228883, 0.54379684,
            -1.2959443
        ],
        [
            -0.69314715, 1.7720417,  -1.0762483,    1.1799437,  0,          0,
            0,           -2.0399152, 0.00031033912, -88.127798, 0.23438509, 0.51241663,
            -1.2294521
        ],
        [
            -0.5108256, 1.6208981,  -0.9281975,    1.1771404,  0,          0,
            0,          -2.0291952, 0.00031034517, -88.027798, 0.22723225, 0.48686433,
            -1.178403
        ],
        [
            -0.35667492, 1.4956127,  -0.8048661,    1.1762689,  0,          0,
            0,           -2.0258862, 0.00031035263, -87.927798, 0.22071588, 0.46543791,
            -1.1375245
        ],
        [
            -0.22314352, 1.388975,   -0.69941126,   1.1766108,  0,          0,
            0,           -2.0271847, 0.00031036153, -87.827798, 0.21474542, 0.44707782,
            -1.1037836
        ],
        [
            2.8108213e-08, 1.2147195,  -0.52605444,   1.1793508,  0,          0,
            0,             -2.0376439, 0.00031038366, -87.627798, 0.20416367, 0.41695408,
            -1.0507906
        ],
        [
            0.26236429, 1.0161348,  -0.32679854,   1.1857132,  0,          0,
            0,          -2.0623608, 0.00031042766, -87.327798, 0.19093466, 0.38267298,
            -0.9937386
        ],
        [
            0.53062828, 0.81997574, -0.12808208,   1.1951503,  0,          0,
            0,          -2.1001877, 0.00031050644, -86.927798, 0.17682735, 0.34916911,
            -0.94099255
        ],
        [
            0.69314721, 0.70441271, -0.010100711,  1.2020191,  0,         0,
            0,          -2.1286508, 0.00031058058, -86.627798, 0.1681455, 0.32973127,
            -0.91164074
        ],
        [
            0.9932518, 0.49717699, 0.20317536,    1.2163819,  0,          0,
            0,         -2.1909336, 0.00031080357, -85.927798, 0.15216129, 0.29568603,
            -0.86227144
        ],
        [
            1.0986123, 0.42622507, 0.27669159,    1.2218056,  0,         0,
            0,         -2.2155054, 0.00031092055, -85.627798, 0.1466343, 0.28432522,
            -0.84634419
        ],
        [
            1.2089604, 0.35287016, 0.35295743,    1.2276355,  0,          0,
            0,         -2.2426131, 0.00031107327, -85.277797, 0.14092438, 0.2727637,
            -0.83040469
        ],
        [
            1.3862944, 0.23694804, 0.47400292,    1.2372325,  0,          0,
            0,         -2.2889166, 0.00031140335, -84.627797, 0.13196015, 0.254912,
            -0.80631223
        ],
        [
            1.6094379, 0.094320202, 0.62377293,    1.2495213,  0,          0,
            0,         -2.3515595,  0.00031202947, -83.627797, 0.12113097, 0.2337259,
            -0.77851553
        ],
        [
            1.7917595, -0.019722483, 0.74414793,   1.2595849,  0,          0,
            0,         -2.4060046,   0.0003127999, -82.627797, 0.11271148, 0.21746325,
            -0.75775321
        ],
        [
            1.9459102, -0.11452035, 0.84459679,    1.2680224,  0,          0,
            0,         -2.4541098,  0.00031371601, -81.627797, 0.10591836, 0.20443624,
            -0.74148024
        ],
        [
            2.0794416, -0.19551042, 0.9306722,     1.2752324,  0,          0,
            0,         -2.4971949,  0.00031477952, -80.627797, 0.10028381, 0.19367649,
            -0.72828245
        ],
        [
            2.3025851, -0.32868215, 1.0726712,     1.286988,   0,           0,
            0,         -2.5718811,  0.00031735731, -78.627797, 0.091393797, 0.17675396,
            -0.70798377
        ],
        [
            2.7080502, -0.56445101, 1.3252576,     1.3070611,  0,          0,
            0,         -2.7151883,  0.00032656714, -73.627797, 0.07690721, 0.14922836,
            -0.67628043
        ],
        [
            2.9957323, -0.72751338, 1.5006308,     1.320004,   0,           0,
            0,         -2.8220041,  0.00034014267, -68.627797, 0.067876875, 0.13204839,
            -0.65750392
        ],
        [
            3.4011974, -0.95235586, 1.7430334,     1.335826,   0,           0,
            0,         -2.9786754,  0.00038383269, -58.627797, 0.056764766, 0.11083621,
            -0.6359729
        ],
        [
            3.912023, -1.2288233, 2.0415113,     1.3489865,  0,           0,
            0,        -3.1846361, 0.00059748639, -38.627796, 0.045151655, 0.088539807,
            -0.61806205
        ],
        [
            4.2484953, -1.4075866, 2.2344798,    1.3440788,  0,           0,
            0,         -3.3247538, 0.0015499806, -18.627796, 0.038762162, 0.076196689,
            -0.61857429
        ]
    ],

    blast_info => [ 0, 0, 0, 0, 0, 0, 1.4383435, 0, 0, 0, 0, 0, 0.29559158, 0.28525783, 0, 0, 0, 0, 0, 0 ],

};
1;
