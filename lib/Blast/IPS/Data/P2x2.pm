package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='P', gamma=2.2
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'P2.2'} = {

    shock_table_info => [ 0, 2.2, 3.67e-06, 4000, 1.8e-06, 100, 1.28e-06, 1.89e-06, 5e-07 ],

    shock_table => [
        [ -11.34703608, 12.00046192, -0.99999156, -11.34897524, 0.99902949 ],
        [ -9.60570518,  10.25917021, -0.99995224, -9.61034301,  0.99767584 ],
        [ -8.40545522,  9.05903129,  -0.99984117, -8.41392256,  0.99574917 ],
        [ -7.47949456,  8.13331257,  -0.99959953, -7.49297924,  0.99321527 ],
        [ -6.72098157,  7.37525376,  -0.99914631, -6.7407417,   0.99003194 ],
        [ -6.08243109,  6.73746539,  -0.99838701, -6.10971313,  0.9861981 ],
        [ -5.52885407,  6.18507902,  -0.99720407, -5.56496846,  0.98167482 ],
        [ -5.03900713,  5.69699475,  -0.99546029, -5.08532924,  0.97642417 ],
        [ -4.59727016,  5.25776994,  -0.99299037, -4.65528891,  0.97038597 ],
        [ -4.19239667,  4.85637738,  -0.98959749, -4.26375335,  0.96348464 ],
        [ -3.81626419,  4.48496272,  -0.98505292, -3.90278393,  0.9556356 ],
        [ -3.46075028,  4.13577184,  -0.97905743, -3.56458426,  0.94669218 ],
        [ -3.11898207,  3.80243338,  -0.97122978, -3.24273439,  0.93646109 ],
        [ -2.7822123,   3.47699811,  -0.96099246, -2.92930832,  0.92459312 ],
        [ -2.43856968,  3.14899928,  -0.94737978, -2.6139467,   0.91045607 ],
        [ -2.0551628,   2.78936945,  -0.92776377, -2.26830215,  0.89209815 ],
        [ -1.7021481,   2.46571123,  -0.90518201, -1.95672479,  0.87275637 ],
        [ -1.37754749,  2.17578259,  -0.88060236, -1.67658091,  0.85303709 ],
        [ -1.0740138,   1.9123615,   -0.85466278, -1.42064034,  0.8331642 ],
        [ -0.79195497,  1.67494165,  -0.82854574, -1.18835924,  0.81375132 ],
        [ -0.51820224,  1.45174056,  -0.80198988, -0.96823839,  0.79435852 ],
        [ -0.24358567,  1.23522361,  -0.77485326, -0.75279351,  0.774694 ],
        [ 0.03936231,   1.01992,     -0.74709979, -0.53645064,  0.75456211 ],
        [ 0.32921798,   0.80738739,  -0.71957703, -0.32067172,  0.73441619 ],
        [ 0.6274726,    0.59680743,  -0.69280738, -0.10462271,  0.7144988 ],
        [ 0.93744545,   0.3861002,   -0.66709923, 0.11378675,   0.69492627 ],
        [ 1.26394277,   0.17236382,  -0.64262801, 0.33750459,   0.67574062 ],
        [ 1.61131193,   -0.04677961, -0.61962236, 0.5689374,    0.65704778 ],
        [ 1.98527246,   -0.27439028, -0.59824276, 0.81119424,   0.63892146 ],
        [ 2.39255659,   -0.51393002, -0.57864169, 1.06778052,   0.62144666 ],
        [ 2.79739689,   -0.74482057, -0.56252904, 1.31622089,   0.6062466 ],
        [ 3.2075407,    -0.97270166, -0.54915078, 1.56205424,   0.59283402 ],
        [ 3.62925073,   -1.20184812, -0.53799924, 1.80947745,   0.58088505 ],
        [ 4.0669742,    -1.43523792, -0.52873555, 2.06134681,   0.5702007 ],
        [ 4.52465856,   -1.67541266, -0.52109861, 2.32007225,   0.56063662 ],
        [ 5.0071372,    -1.92526066, -0.51486242, 2.58844237,   0.55206473 ],
        [ 5.51546999,   -2.18565396, -0.50987501, 2.86708166,   0.54444218 ],
        [ 6.05865537,   -2.46148249, -0.50592832, 3.16090682,   0.53762526 ],
        [ 6.6430263,    -2.7561904,  -0.50288272, 3.47324742,   0.53155216 ],
        [ 7.26596921,   -3.06871363, -0.50064011, 3.80266603,   0.5262494 ],
        [ 7.95620904,   -3.41367124, -0.49901452, 4.16420657,   0.52150729 ],
        [ 8.70861634,   -3.78869367, -0.49794482, 4.5549896,    0.51741153 ],
        [ 9.53475879,   -4.19977773, -0.49732022, 4.98093772,   0.5139146 ],
        [ 10.50155799,  -4.68042046, -0.49703793, 5.47622145,   0.5108249 ],
        [ 11.58983197,  -5.22131606, -0.49704669, 6.03068246,   0.50828765 ],
        [ 12.87195134,  -5.85872097, -0.49727803, 6.68093425,   0.50618122 ],
        [ 14.43395407,  -6.63577796, -0.49768154, 7.47014514,   0.50445049 ],
        [ 16.45626197,  -7.64278981, -0.49821349, 8.48873665,   0.50301922 ],
        [ 19.84771533,  -9.33372507, -0.49891544, 10.19215542,  0.50166923 ],
        [ 19.8877186,   -9.35369198, -0.49977018, 10.21222358,  0.50165375 ],
        [ 19.9445232,   -9.38209107, -0.49999229, 10.24071858,  0.50160953 ],
        [ 20.26514831,  -9.54240782, -0.5000169,  10.4015081,   0.50137164 ]
    ],

    energy_table => [
        [ -11.34703608, 0.0015718,  0.00085289, 0.0015718,  0.00085289, 11.688982,   -0.94401349 ],
        [ -9.60570518,  0.00403293, 0.00217444, 0.00403293, 0.00217444, 9.996218,    -0.99591205 ],
        [ -8.40545522,  0.0076839,  0.00410759, 0.0076839,  0.00410759, 8.7811824,   -1.0248355 ],
        [ -7.47949456,  0.01257412, 0.00664932, 0.01257412, 0.00664932, 7.8232837,   -1.0424932 ],
        [ -6.72098157,  0.01873268, 0.00977593, 0.01873268, 0.00977593, 7.0275674,   -1.0545708 ],
        [ -6.08243109,  0.02607717, 0.01339862, 0.02607717, 0.01339862, 6.3512018,   -1.0630915 ],
        [ -5.52885407,  0.03457246, 0.0174461,  0.03457246, 0.0174461,  5.7608409,   -1.0691469 ],
        [ -5.03900713,  0.04415984, 0.02182905, 0.04415984, 0.02182905, 5.2359525,   -1.0733343 ],
        [ -4.59727016,  0.05480199, 0.02646156, 0.05480199, 0.02646156, 4.7611045,   -1.0760275 ],
        [ -4.19239667,  0.06646892, 0.03125318, 0.06646892, 0.03125318, 4.3250521,   -1.0774628 ],
        [ -3.81626419,  0.07912597, 0.03610327, 0.07912597, 0.03610327, 3.9196248,   -1.0777897 ],
        [ -3.46075028,  0.0928141,  0.04092855, 0.0928141,  0.04092855, 3.5364874,   -1.0770959 ],
        [ -3.11898207,  0.1076074,  0.04563705, 0.1076074,  0.04563705, 3.1685688,   -1.0754063 ],
        [ -2.7822123,   0.1237435,  0.05015316, 0.1237435,  0.05015316, 2.8067723,   -1.0726597 ],
        [ -2.43856968,  0.14172183, 0.05439787, 0.14172183, 0.05439787, 2.4387408,   -1.0686373 ],
        [ -2.0551628,   0.1633717,  0.05837338, 0.1633717,  0.05837338, 2.0300165,   -1.0626467 ],
        [ -1.7021481,   0.1844786,  0.06102342, 0.1844786,  0.06102342, 1.6559872,   -1.0556493 ],
        [ -1.37754749,  0.20454011, 0.0624036,  0.20454011, 0.0624036,  1.3144807,   -1.0455364 ],
        [ -1.0740138,   0.22354965, 0.06268581, 0.22354965, 0.06268581, 0.99898305,  -1.0311314 ],
        [ -0.79195497,  0.24116296, 0.06206634, 0.24116296, 0.06206634, 0.71031468,  -1.0148335 ],
        [ -0.51820224,  0.25798322, 0.0607012,  0.25798322, 0.0607012,  0.4347861,   -0.99826971 ],
        [ -0.24358567,  0.27438723, 0.05866512, 0.27438723, 0.05866512, 0.16290844,  -0.98276513 ],
        [ 0.03936231,   0.29062052, 0.05599333, 0.29062052, 0.05599333, -0.11304606, -0.96935734 ],
        [ 0.32921798,   0.30639717, 0.05280158, 0.30639717, 0.05280158, -0.39226011, -0.95901906 ],
        [ 0.6274726,    0.32161426, 0.04919867, 0.32161426, 0.04919867, -0.67698232, -0.95201771 ],
        [ 0.93744545,   0.33625913, 0.04527484, 0.33625913, 0.04527484, -0.97124106, -0.94817258 ],
        [ 1.26394277,   0.35035911, 0.04110104, 0.35035911, 0.04110104, -1.2804275,  -0.94709295 ],
        [ 1.61131193,   0.36387734, 0.03675799, 0.36387734, 0.03675799, -1.6094601,  -0.94833192 ],
        [ 1.98527246,   0.37678412, 0.03231851, 0.37678412, 0.03231851, -1.9645492,  -0.95144399 ],
        [ 2.39255659,   0.38902425, 0.02785845, 0.38902425, 0.02785845, -2.3529035,  -0.95604307 ],
        [ 2.79739689,   0.39947885, 0.02386474, 0.39947885, 0.02386474, -2.7409629,  -0.96132784 ],
        [ 3.2075407,    0.40851585, 0.02027961, 0.40851585, 0.02027961, -3.1364,     -0.96710584 ],
        [ 3.62925073,   0.41637352, 0.01706331, 0.41637352, 0.01706331, -3.5455231,  -0.97730937 ],
        [ 4.0669742,    0.42319804, 0.01419479, 0.42319804, 0.01419479, -3.976566,   -0.99302383 ],
        [ 4.52465856,   0.42909723, 0.01165781, 0.42909723, 0.01165781, -4.4350227,  -1 ],
        [ 5.0071372,    0.43416792, 0.00943307, 0.43416792, 0.00943307, -4.9175014,  -1 ],
        [ 5.51546999,   0.43845898, 0.007517,   0.43845898, 0.007517,   -5.4258342,  -1 ],
        [ 6.05865537,   0.44207901, 0.00587524, 0.44207901, 0.00587524, -5.9690195,  -1 ],
        [ 6.6430263,    0.44509049, 0.00449046, 0.44509049, 0.00449046, -6.5533905,  -1 ],
        [ 7.26596921,   0.44751932, 0.00335989, 0.44751932, 0.00335989, -7.1763334,  -1 ],
        [ 7.95620904,   0.44949997, 0.00242792, 0.44949997, 0.00242792, -7.8665732,  -1 ],
        [ 8.70861634,   0.45103627, 0.00169809, 0.45103627, 0.00169809, -8.6189805,  -1 ],
        [ 9.53475879,   0.45219506, 0.00114303, 0.45219506, 0.00114303, -9.445123,   -1 ],
        [ 10.50155799,  0.45307838, 0.00071681, 0.45307838, 0.00071681, -10.411922,  -1 ],
        [ 11.58983197,  0.45368443, 0.0004225,  0.45368443, 0.0004225,  -11.500196,  -1 ],
        [ 12.87195134,  0.4540871,  0.00022586, 0.4540871,  0.00022586, -12.782316,  -1 ],
        [ 14.43395407,  0.45433354, 0.0001049,  0.45433354, 0.0001049,  -14.344318,  -1 ],
        [ 16.45626197,  0.45446783, 3.869e-05,  0.45446783, 3.869e-05,  -16.366626,  -1 ],
        [ 19.84771533,  0.4545314,  7.21e-06,   0.4545314,  7.21e-06,   -19.758079,  -1 ],
        [ 19.8877186,   0.45453169, 7.06e-06,   0.45453169, 7.06e-06,   -19.798083,  -1 ],
        [ 19.9445232,   0.45453209, 6.87e-06,   0.45453209, 6.87e-06,   -19.854887,  -1 ],
        [ 20.26514831,  0.45453412, 5.85e-06,   0.45453412, 5.85e-06,   -20.175512,  -1 ]
    ],

    impulse_table => [
        [
            -13.161469, 13.814887,  -13.162251,    1.1665846,  0,          0,
            0,          -7.7022099, 0.00017229856, -89.413164, 0.27111323, 0.99870908,
            -7.2716447
        ],
        [
            -2.3025851, 3.0205976,  -2.4905558,    0.99097439, 0,          0,
            0,          -2.4275023, 0.00017229864, -89.313166, 0.25607665, 0.67170911,
            -1.7201092
        ],
        [
            -1.6094379, 2.3820982,  -1.8760631,    0.95004344, 0,         0,
            0,          -2.2189567, 0.00017229913, -89.213166, 0.2433861, 0.58164896,
            -1.4777138
        ],
        [
            -1.3862943,   2.1834882,  -1.6840447, 0.93926454, 0, 0, 0, -2.1705755,
            0.0001722996, -89.163166, 0.23772804, 0.55121358, -1.4074991
        ],
        [
            -1.2039728, 2.02418,    -1.5294826,    0.93195281, 0,          0,
            0,          -2.1390406, 0.00017230024, -89.113166, 0.23245432, 0.5261189,
            -1.3531007
        ],
        [
            -0.91629071,  1.7786908,  -1.2900768, 0.9237482,  0, 0, 0, -2.1047999,
            0.0001723021, -89.013166, 0.22289319, 0.48647524, -1.27278
        ],
        [
            -0.69314715, 1.5935427,  -1.1082975,    0.92051786, 0,          0,
            0,           -2.0916341, 0.00017230474, -88.913166, 0.21443085, 0.45599369,
            -1.2151377
        ],
        [
            -0.5108256,   1.4458272,  -0.96238064, 0.91998267, 0, 0, 0, -2.0894699,
            0.0001723082, -88.813166, 0.20686853,  0.43145339, -1.1710348
        ],
        [
            -0.35667492, 1.3234843,  -0.84086133,   0.92097007, 0,         0,
            0,           -2.0934679, 0.00017231247, -88.713166, 0.2000561, 0.41106047,
            -1.1358138
        ],
        [
            -0.22314353, 1.2194046,  -0.73697208,   0.92283475, 0,          0,
            0,           -2.1010614, 0.00017231757, -88.613166, 0.19387723, 0.39371488,
            -1.1068058
        ],
        [
            2.6144161e-08, 1.0494028,  -0.56620667,   0.92785467, 0,          0,
            0,             -2.1217952, 0.00017233021, -88.413166, 0.18306992, 0.36550234,
            -1.0613595
        ],
        [
            0.26236429,   0.85570195, -0.36992338, 0.93636086, 0, 0, 0, -2.1579412,
            0.0001723553, -88.113166, 0.16980971,  0.33375011, -1.0125823
        ],
        [
            0.53062828,   0.6643116,  -0.17412513, 0.94731523, 0, 0, 0, -2.2065041,
            0.0001724001, -87.713166, 0.15596434,  0.30305966, -0.96761752
        ],
        [
            0.69314721, 0.55149336, -0.057838189,  0.9547702,  0,          0,
            0,          -2.2409578, 0.00017244217, -87.413166, 0.14758895, 0.28540059,
            -0.94264695
        ],
        [
            0.9932518, 0.34899435, 0.15247331,    0.96960639, 0,          0,
            0,         -2.3132798, 0.00017256845, -86.713166, 0.13244264, 0.25471632,
            -0.90072281
        ],
        [
            1.0986123, 0.27959603, 0.22500223,    0.97502772, 0,          0,
            0,         -2.3410725, 0.00017263461, -86.413166, 0.12728351, 0.24454391,
            -0.88721627
        ],
        [
            1.2089604, 0.20780469, 0.30026508,    0.98077224, 0,         0,
            0,         -2.3713932, 0.00017272092, -86.063166, 0.1219934, 0.23422502,
            -0.87370775
        ],
        [
            1.3862944, 0.094258589, 0.41976496,    0.99007358, 0,          0,
            0,         -2.4225341,  0.00017290734, -85.413166, 0.11376569, 0.21835637,
            -0.85330442
        ],
        [
            1.6094379, -0.045618335, 0.56770601,    1.001756,   0,          0,
            0,         -2.4907356,   0.00017326067, -84.413166, 0.10394365, 0.19962251,
            -0.82978323
        ],
        [
            1.7917595, -0.15760792, 0.68668068,    1.0111699,  0,           0,
            0,         -2.549326,   0.00017369515, -83.413166, 0.096388366, 0.18531311,
            -0.81222496
        ],
        [
            1.9459102, -0.25080107, 0.78600925,    1.0189745,  0,           0,
            0,         -2.6006848,  0.00017421158, -82.413166, 0.090339824, 0.17389421,
            -0.79846731
        ],
        [
            2.0794416, -0.33049613, 0.87116013,    1.0255891,  0,           0,
            0,         -2.6464177,  0.00017481093, -81.413166, 0.085352604, 0.16449149,
            -0.7873105
        ],
        [
            2.3025851, -0.4616894, 1.0117041,     1.0362819,  0,           0,
            0,         -2.7252088, 0.00017626327, -79.413166, 0.077534696, 0.1497553,
            -0.77014829
        ],
        [
            2.7080502, -0.69441475, 1.2619136,     1.0543421,  0,           0,
            0,         -2.8750885,  0.00018145051, -74.413165, 0.064915273, 0.12592071,
            -0.74331129
        ],
        [
            2.9957323, -0.85570413, 1.4357872,     1.0659221,  0,           0,
            0,         -2.9859862,  0.00018909574, -69.413165, 0.057116475, 0.1111289,
            -0.72736242
        ],
        [
            3.4011974, -1.0785154, 1.6763023,     1.0802011,  0,           0,
            0,         -3.1477439, 0.00021370246, -59.413165, 0.047584349, 0.092955976,
            -0.70891656
        ],
        [
            3.912023, -1.353078, 1.9727175,     1.093332,   0,           0,
            0,        -3.359211, 0.00033404141, -39.413165, 0.037693718, 0.073965564,
            -0.69290685
        ],
        [
            4.2484953, -1.5309148, 2.1644859,     1.0937759,  0,           0,
            0,         -3.5024708, 0.00086785938, -19.413164, 0.032281986, 0.063506638,
            -0.69126068
        ]
    ],

    blast_info => [ 0, 0, 0, 0, 0, 0, 1.166088, 0, 0, 0, 0, 0, 0.31055756, 0.27111323, 0.5370386, 0, 0, 0, 0, 0 ],

};
1;
