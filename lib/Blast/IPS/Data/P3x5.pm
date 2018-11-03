package Blast::IPS::Data;
use strict;
use warnings;

# Blast data for symmetry='P', gamma=3.5
# For data definitions see the README file in Blast::IPS::Data.

our $rBlastData;
$rBlastData->{'P3x5'} = {

    table_name => 'P3x5',
    symmetry   => 0,
    gamma      => 3.5,

    shock_table_info => [ 0, 3.5, 1.83e-06, 4000, 1.49e-06, 100, 1.02e-06, 3.07e-07, 5e-07 ],

    shock_table => [
        [ -10.6434665, 12.00046224, -0.99999045, -10.64552918, 0.99896761 ],
        [ -9.62565953, 10.98267073, -0.99997357, -9.62909304,  0.99828034 ],
        [ -9.08582821, 10.44285674, -0.99995739, -9.09032796,  0.99774517 ],
        [ -7.77873663, 9.13587802,  -0.99984545, -7.78740412,  0.99564826 ],
        [ -6.85027674, 8.20765475,  -0.99960901, -6.86409835,  0.99304462 ],
        [ -6.09051415, 7.44833666,  -0.99916537, -6.11078229,  0.98977322 ],
        [ -5.44883965, 6.80741222,  -0.99841797, -5.47686989,  0.98581475 ],
        [ -4.89338821, 6.25313381,  -0.99725227, -4.9305324,   0.98114378 ],
        [ -4.40166685, 5.76315325,  -0.99552955, -4.44936132,  0.97571183 ],
        [ -3.95868027, 5.32264942,  -0.99308731, -4.01846449,  0.96946367 ],
        [ -3.55318198, 4.92059083,  -0.98973217, -3.62674494,  0.96232502 ],
        [ -3.1764247,  4.54850037,  -0.9852321,  -3.26566307,  0.95419916 ],
        [ -2.82091079, 4.1992352,   -0.97929863, -2.92802656,  0.94495087 ],
        [ -2.47914258, 3.86580221,  -0.97154375, -2.60683019,  0.93436785 ],
        [ -2.14326026, 3.54110087,  -0.96141849, -2.29499465,  0.92212433 ],
        [ -1.80211686, 3.21530852,  -0.94801133, -1.9828306,   0.90761826 ],
        [ -1.4263705,  2.86252679,  -0.92895726, -1.64519627,  0.8890635 ],
        [ -1.07213939, 2.53731283,  -0.90646126, -1.33373225,  0.86907664 ],
        [ -0.74696903, 2.24645183,  -0.8819121,  -1.05439201,  0.84874368 ],
        [ -0.44147268, 1.98095398,  -0.85577712, -0.79821346,  0.82817874 ],
        [ -0.15941385, 1.74323446,  -0.82954242, -0.56741101,  0.80825673 ],
        [ 0.11371408,  1.52028677,  -0.80285622, -0.34935036,  0.78845474 ],
        [ 0.38670003,  1.30483086,  -0.77562201, -0.13683548,  0.76851594 ],
        [ 0.66746011,  1.09099056,  -0.74775673, 0.07607492,   0.74822044 ],
        [ 0.9555648,   0.87958272,  -0.72001692, 0.28870558,   0.72796864 ],
        [ 1.25144668,  0.67057834,  -0.69303483, 0.50113089,   0.70808747 ],
        [ 1.55874009,  0.46165805,  -0.66709138, 0.71570398,   0.68867391 ],
        [ 1.88173097,  0.2502544,   -0.6424008,  0.93504751,   0.66979974 ],
        [ 2.22502066,  0.03380166,  -0.61916417, 1.1617964,    0.65154821 ],
        [ 2.59436853,  -0.19078639, -0.59753419, 1.39913552,   0.63398292 ],
        [ 2.99635669,  -0.42687323, -0.5776674,  1.65053489,   0.61718495 ],
        [ 3.39873096,  -0.65589004, -0.56118894, 1.89586703,   0.60258745 ],
        [ 3.80607048,  -0.88160219, -0.54749915, 2.13866004,   0.58981973 ],
        [ 4.224425,    -1.10817793, -0.53608433, 2.38299485,   0.57854386 ],
        [ 4.65918102,  -1.33909862, -0.52657816, 2.63228377,   0.56852085 ],
        [ 5.11239747,  -1.57590733, -0.5187543,  2.88787731,   0.55963174 ],
        [ 5.58927564,  -1.82170024, -0.51236801, 3.15281322,   0.55171891 ],
        [ 6.09712076,  -2.08053062, -0.5072088,  3.43114972,   0.54463805 ],
        [ 6.63112908,  -2.35025695, -0.50319703, 3.72028039,   0.53842086 ],
        [ 7.20565562,  -2.63841689, -0.50010671, 4.02797318,   0.53287414 ],
        [ 7.82870444,  -2.94924563, -0.49781622, 4.3583891,    0.52793459 ],
        [ 8.50452231,  -3.28509748, -0.49622539, 4.71365355,   0.52358154 ],
        [ 9.23820639,  -3.6487657,  -0.49522608, 5.09635173,   0.51978361 ],
        [ 10.06001356, -4.05549656, -0.49470227, 5.52207299,   0.51641343 ],
        [ 10.98214342, -4.51159192, -0.49458238, 5.99685758,   0.51347017 ],
        [ 12.03086526, -5.0303578,  -0.49479161, 6.53394116,   0.5109114 ],
        [ 13.24723766, -5.63248006, -0.49526395, 7.15397891,   0.50868549 ],
        [ 14.69629556, -6.35063039, -0.49593976, 7.88960701,   0.50673708 ],
        [ 16.55266639, -7.27207977, -0.49678894, 8.82854968,   0.50495909 ],
        [ 17.31664576, -7.65173924, -0.49711183, 9.21410526,   0.50438886 ]
    ],

    energy_table => [
        [ -10.6434665, 0.00021932, 0.00015394, 0.00021932, 0.00015394, 11.528889,   -0.98443573 ],
        [ -9.62565953, 0.0004471,  0.00031182, 0.0004471,  0.00031182, 10.521688,   -0.99438011 ],
        [ -9.08582821, 0.00065098, 0.00045206, 0.00065098, 0.00045206, 9.983516,    -0.99913851 ],
        [ -7.77873663, 0.00160387, 0.00109794, 0.00160387, 0.00109794, 8.6705496,   -1.0091786 ],
        [ -6.85027674, 0.0030148,  0.00203314, 0.0030148,  0.00203314, 7.7304786,   -1.0153652 ],
        [ -6.09051415, 0.0050121,  0.00332408, 0.0050121,  0.00332408, 6.9572643,   -1.0197628 ],
        [ -5.44883965, 0.00764211, 0.00497542, 0.00764211, 0.00497542, 6.3017934,   -1.0230392 ],
        [ -4.89338821, 0.01093291, 0.00697413, 0.01093291, 0.00697413, 5.732805,    -1.0255512 ],
        [ -4.40166685, 0.01491048, 0.00929956, 0.01491048, 0.00929956, 5.2280059,   -1.0275182 ],
        [ -3.95868027, 0.01959066, 0.01191818, 0.01959066, 0.01191818, 4.7724611,   -1.0290767 ],
        [ -3.55318198, 0.02498999, 0.01478945, 0.02498999, 0.01478945, 4.3549021,   -1.0303188 ],
        [ -3.1764247,  0.03112961, 0.01786647, 0.03112961, 0.01786647, 3.9665204,   -1.0313078 ],
        [ -2.82091079, 0.03804745, 0.02109961, 0.03804745, 0.02109961, 3.5997238,   -1.0320887 ],
        [ -2.47914258, 0.04582466, 0.02444269, 0.04582466, 0.02444269, 3.2468729,   -1.0326934 ],
        [ -2.14326026, 0.05460502, 0.02784888, 0.05460502, 0.02784888, 2.8999217,   -1.0331414 ],
        [ -1.80211686, 0.06469443, 0.03128275, 0.06469443, 0.03128275, 2.5474076,   -1.0334616 ],
        [ -1.4263705,  0.07712383, 0.03480932, 0.07712383, 0.03480932, 2.1590342,   -1.0327257 ],
        [ -1.07213939, 0.08997399, 0.03764319, 0.08997399, 0.03764319, 1.7935044,   -1.0295707 ],
        [ -0.74696903, 0.10255531, 0.0396276,  0.10255531, 0.0396276,  1.4594127,   -1.024109 ],
        [ -0.44147268, 0.11486222, 0.04082786, 0.11486222, 0.04082786, 1.147506,    -1.0170226 ],
        [ -0.15941385, 0.12646025, 0.04130809, 0.12646025, 0.04130809, 0.8616773,   -1.0093001 ],
        [ 0.11371408,  0.1377396,  0.04119221, 0.1377396,  0.04119221, 0.58708447,  -1.0014004 ],
        [ 0.38670003,  0.14890675, 0.04053705, 0.14890675, 0.04053705, 0.31479752,  -0.993763 ],
        [ 0.66746011,  0.16013431, 0.03936529, 0.16013431, 0.03936529, 0.036850123, -0.98675083 ],
        [ 0.9555648,   0.17124931, 0.03772953, 0.17124931, 0.03772953, -0.24648218, -0.98085642 ],
        [ 1.25144668,  0.18212021, 0.03570267, 0.18212021, 0.03570267, -0.5359162,  -0.97637329 ],
        [ 1.55874009,  0.1927343,  0.03334543, 0.1927343,  0.03334543, -0.83536374, -0.97338292 ],
        [ 1.88173097,  0.20308229, 0.03071392, 0.20308229, 0.03071392, -1.1493892,  -0.97188638 ],
        [ 2.22502066,  0.21313702, 0.02786475, 0.21313702, 0.02786475, -1.4828978,  -0.97187102 ],
        [ 2.59436853,  0.22286924, 0.02485182, 0.22286924, 0.02485182, -1.8420007,  -0.97344868 ],
        [ 2.99635669,  0.23222569, 0.0217336,  0.23222569, 0.0217336,  -2.2338348,  -0.97719181 ],
        [ 3.39873096,  0.24038144, 0.01884758, 0.24038144, 0.01884758, -2.628019,   -0.98396664 ],
        [ 3.80607048,  0.24751017, 0.01620163, 0.24751017, 0.01620163, -3.0306093,  -1.0115547 ],
        [ 4.224425,    0.25377199, 0.0137845,  0.25377199, 0.0137845,  -3.4637719,  -1.0200892 ],
        [ 4.65918102,  0.25927612, 0.01158837, 0.25927612, 0.01158837, -3.9003462,  -1 ],
        [ 5.11239747,  0.26407034, 0.00962043, 0.26407034, 0.00962043, -4.3535626,  -1 ],
        [ 5.58927564,  0.26822862, 0.00787093, 0.26822862, 0.00787093, -4.8304408,  -1 ],
        [ 6.09712076,  0.27182074, 0.00632649, 0.27182074, 0.00632649, -5.3382859,  -1 ],
        [ 6.63112908,  0.27483386, 0.00500599, 0.27483386, 0.00500599, -5.8722943,  -1 ],
        [ 7.20565562,  0.27737197, 0.00387474, 0.27737197, 0.00387474, -6.4468208,  -1 ],
        [ 7.82870444,  0.27947629, 0.00292259, 0.27947629, 0.00292259, -7.0698696,  -1 ],
        [ 8.50452231,  0.28117515, 0.00214352, 0.28117515, 0.00214352, -7.7456875,  -1 ],
        [ 9.23820639,  0.28250841, 0.00152481, 0.28250841, 0.00152481, -8.4793716,  -1 ],
        [ 10.06001356, 0.28354855, 0.001037,   0.28354855, 0.001037,   -9.3011787,  -1 ],
        [ 10.98214342, 0.28432359, 0.00067005, 0.28432359, 0.00067005, -10.223309,  -1 ],
        [ 12.03086526, 0.28487651, 0.00040602, 0.28487651, 0.00040602, -11.27203,   -1 ],
        [ 13.24723766, 0.28525045, 0.00022607, 0.28525045, 0.00022607, -12.488403,  -1 ],
        [ 14.69629556, 0.28548586, 0.00011196, 0.28548586, 0.00011196, -13.937461,  -1 ],
        [ 16.55266639, 0.2856226,  4.523e-05,  0.2856226,  4.523e-05,  -15.793832,  -1 ],
        [ 17.31664576, 0.28565143, 3.11e-05,   0.28565143, 3.11e-05,   -16.557811,  -1 ]
    ],

    impulse_table => [
        [
            -12.457899, 13.814887,  -12.458731,    2.3103609,  0,          0,
            0,          -7.4402887, 0.00094001591, -86.573122, 0.31183246, 0.99981817,
            -8.8003289
        ],
        [
            -2.3025851,   3.694697,   -2.4423996, 2.1424086,  0, 0, 0, -2.470254,
            0.0009400167, -86.473126, 0.30283912, 0.82396456, -2.0739687
        ],
        [
            -1.6094379, 3.0335103,  -1.8088236,    2.0934879,  0,          0,
            0,          -2.2161902, 0.00094002004, -86.373126, 0.29468816, 0.75184656,
            -1.7307153
        ],
        [
            -1.3862943, 2.825344,   -1.6096088,    2.0778508,  0,          0,
            0,          -2.1469311, 0.00094002289, -86.323126, 0.29088501, 0.72516425,
            -1.628632
        ],
        [
            -1.2039728, 2.6574081,  -1.4488178,    2.0656173,  0,          0,
            0,          -2.0959052, 0.00094002662, -86.273126, 0.28724439, 0.70224166,
            -1.548563
        ],
        [
            -0.91629072, 2.3969129,  -1.1990246,    2.0479997,  0,          0,
            0,           -2.0267159, 0.00094003682, -86.173126, 0.28040317, 0.66427357,
            -1.4286214
        ],
        [
            -0.69314717, 2.1991039,  -1.0088059,    2.036347,   0,          0,
            0,           -1.9834506, 0.00094005087, -86.073126, 0.27408177, 0.63357131,
            -1.3411813
        ],
        [
            -0.51082561,  2.0405194,  -0.85581613, 2.0284869,  0, 0, 0, -1.9552903,
            0.0009400689, -85.973126, 0.26821348,  0.60787232, -1.2734615
        ],
        [
            -0.35667493, 1.9087118,  -0.72823604, 2.0231741,  0, 0, 0, -1.9366972,
            0.000940091, -85.873126, 0.26274372,  0.58583716, -1.218851
        ],
        [
            -0.22314354, 1.796295,   -0.61906659,   2.0196384,  0,          0,
            0,           -1.9245133, 0.00094011724, -85.773126, 0.25762712, 0.56660223,
            -1.1735136
        ],
        [
            1.5025249e-08, 1.6122224,  -0.43948028,   2.0160492,  0,          0,
            0,             -1.9122998, 0.00094018228, -85.573126, 0.24830683, 0.53433691,
            -1.1018216
        ],
        [
            0.26236428, 1.4020416,  -0.23295377,   2.0155954,  0,          0,
            0,          -1.9107809, 0.00094031155, -85.273126, 0.23618846, 0.49655434,
            -1.0239727
        ],
        [
            0.53062827, 1.194231,  -0.026977798,  2.0191763,  0,          0,
            0,          -1.922984, 0.00094054348, -84.873126, 0.22267138, 0.45854205,
            -0.95139352
        ],
        [
            0.6931472, 1.0718153,  0.095270932,   2.0231661,  0,          0,
            0,         -1.9367531, 0.00094076223, -84.573126, 0.21403626, 0.43600147,
            -0.91076055
        ],
        [
            0.99325179,   0.85251416, 0.31609164, 2.0336384,  0, 0, 0, -1.9738336,
            0.0009414222, -83.873126, 0.19748407, 0.39567156, -0.8420419
        ],
        [
            1.0986123, 0.77753994, 0.39213976,    2.0381315,  0,          0,
            0,         -1.9901829, 0.00094176919, -83.573126, 0.19155904, 0.38197547,
            -0.81977734
        ],
        [
            1.2089604, 0.70010269, 0.47098768,    2.0432188,  0,          0,
            0,         -2.0090309, 0.00094222272, -83.223126, 0.18532761, 0.36791642,
            -0.79745209
        ],
        [
            1.3862944, 0.57791574, 0.59602628,    2.0520916,  0,         0,
            0,         -2.0428051, 0.00094320437, -82.573126, 0.1753176, 0.34597048,
            -0.76362952
        ],
        [
            1.6094379, 0.42794106, 0.75054026,    2.0642116,  0,          0,
            0,         -2.0909405, 0.00094506969, -81.573126, 0.16285629, 0.31955423,
            -0.72450034
        ],
        [
            1.7917595, 0.30834765, 0.87455641,    2.0746663,  0,          0,
            0,         -2.1345141, 0.00094736822, -80.573126, 0.15289298, 0.29900613,
            -0.69520919
        ],
        [
            1.9459102, 0.20917315, 0.9779203,     2.0837442,  0,         0,
            0,         -2.1740695, 0.00095010409, -79.573126, 0.1446831, 0.28237792,
            -0.67222063
        ],
        [
            2.0794416, 0.12462209, 1.066402,      2.0916989,  0,          0,
            0,         -2.2101927, 0.00095328247, -78.573126, 0.13776051, 0.26853128,
            -0.65356163
        ],
        [
            2.3025851, -0.014035422, 1.2121829,     2.1050023,  0,          0,
            0,         -2.2740955,   0.00096099241, -76.573126, 0.12663706, 0.24654845,
            -0.62485117
        ],
        [
            2.7080502, -0.25837204, 1.4709226,     2.1284094,  0,          0,
            0,         -2.4002079,  0.00098857013, -71.573126, 0.10801143, 0.2102552,
            -0.58005761
        ],
        [
            2.9957323, -0.42651252, 1.6501495,    2.1436438,  0,           0,
            0,         -2.4964236,  0.0010292627, -66.573126, 0.096108233, 0.18726325,
            -0.55366066
        ],
        [
            3.4011974, -0.65727406, 1.8973532,    2.1613889,  0,           0,
            0,         -2.6400088,  0.0011604234, -56.573126, 0.081174237, 0.15850717,
            -0.52382606
        ],
        [
            3.912023, -0.93944433, 2.200992,     2.168972,   0,           0,
            0,        -2.8320085,  0.0018057829, -36.573125, 0.065247733, 0.12782399,
            -0.5009753
        ],
        [
            4.2484953,    -1.1210746, 2.3969134,   2.1357862,  0, 0, 0, -2.9644252,
            0.0047615059, -16.573125, 0.056350793, 0.11059546, -0.50834228
        ]
    ],

    blast_info => [ 0, 0, 0, 0, 0, 0, 2.3089594, 0, 0, 0, 0, 0, 0.26576444, 0.31183247, 0, 0, 0, 0, 0, 0 ],

};
1;
