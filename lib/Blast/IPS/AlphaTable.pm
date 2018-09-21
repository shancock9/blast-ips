package Blast::IPS::AlphaTable;

# MIT License
# Copyright (c) 2018 Steven Hancock
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# This module provides a routine 'alpha_interpolate' which returns an accurate
# value of the similarity parameter alpha for the three symmetries and
# gamma between 1.1 and 7.

use strict;
use warnings;
use 5.006;

our @EXPORT_OK = qw(
  alpha_interpolate
);

use Exporter;
our @ISA = qw(Exporter);

use Carp;
our $ralpha_table;

use Blast::IPS::MathUtils qw(
  polint
  set_interpolation_points
  locate_2d
);

##################################################################
# ALPHA TABLE
##################################################################
BEGIN {

    # This table of values was created by integrating the Von Neumann
    # analytical solution using the script 'alpha_integrate.pl'. The values
    # are in good agreement with values produced by other software and
    # with the values extracted from the finite difference solutions.

    # The table format is:
    # $ralpha_table->[$symmetry]= [gamma, alpha, err];

    # where:
    #   symmetry = 0,1,2 for plane, cylindrical, spherical
    #   err = the estimated error in the value of alpha

    # The errors for some of the very low gamma values got truncated to zero.
    # I don't know what the actual errors are for them.

    $ralpha_table = [
        [
            [ 1.01, 55.18267536362, 0 ],
            [ 1.02, 11.1702575238,  5.14e-11 ],
            [ 1.03, 7.4623970603,   5.84e-11 ],
            [ 1.04, 5.4548897110,   2.14e-10 ],
            [ 1.05, 4.4917631083,   3.44e-10 ],
            [ 1.06, 3.7474782153,   2.05e-10 ],
            [ 1.07, 3.2150100919,   1.27e-10 ],
            [ 1.08, 2.8149818824,   8.44e-11 ],
            [ 1.09, 2.5032911244,   5.7e-11 ],
            [ 1.1,  2.25347295355,  4e-11 ],
            [ 1.11, 2.0486833331,   4.76e-10 ],
            [ 1.12, 1.8776906316,   3.55e-10 ],
            [ 1.13, 1.7327172142,   3.11e-10 ],
            [ 1.14, 1.6082060428,   2.41e-10 ],
            [ 1.15, 1.5000806175,   2.25e-10 ],
            [ 1.16, 1.4052824260,   2.24e-10 ],
            [ 1.17, 1.3214716330,   1.77e-10 ],
            [ 1.18, 1.2468275306,   1.43e-10 ],
            [ 1.19, 1.1799119950,   1.17e-10 ],
            [ 1.2,  1.1195739005,   9.65e-11 ],
            [ 1.21, 1.0648808370,   8.04e-11 ],
            [ 1.22, 1.0150694456,   6.76e-11 ],
            [ 1.23, 0.9695087058,   5.73e-11 ],
            [ 1.24, 0.89629622369,  3.77e-11 ],
            [ 1.25, 0.88911816980,  4.19e-11 ],
            [ 1.26, 0.85347143540,  3.62e-11 ],
            [ 1.27, 0.82041284472,  3.15e-11 ],
            [ 1.28, 0.7896684402,   4.37e-10 ],
            [ 1.29, 0.7610018468,   3.83e-10 ],
            [ 1.3,  0.7342080226,   3.38e-10 ],
            [ 1.31, 0.7091082184,   2.99e-10 ],
            [ 1.32, 0.6855458807,   2.66e-10 ],
            [ 1.33, 0.6633832988,   2.37e-10 ],
            [ 1.34, 0.6424988427,   2.12e-10 ],
            [ 1.35, 0.6227846745,   1.91e-10 ],
            [ 1.36, 0.6041448395,   1.72e-10 ],
            [ 1.37, 0.5864936668,   1.55e-10 ],
            [ 1.38, 0.5697544223,   1.41e-10 ],
            [ 1.39, 0.5538581681,   1.28e-10 ],
            [ 1.4,  0.5387427924,   1.16e-10 ],
            [ 1.41, 0.5243521814,   1.06e-10 ],
            [ 1.42, 0.5106355091,   9.71e-11 ],
            [ 1.43, 0.4975466262,   8.9e-11 ],
            [ 1.44, 0.4850435310,   8.18e-11 ],
            [ 1.45, 0.4730879124,   7.53e-11 ],
            [ 1.46, 0.4616447504,   6.94e-11 ],
            [ 1.47, 0.4506819689,   6.41e-11 ],
            [ 1.48, 0.4257185242,   5.25e-11 ],
            [ 1.49, 0.41597659404,  4.88e-11 ],
            [ 1.5,  0.4203931636,   5.11e-11 ],
            [ 1.51, 0.41108010536,  4.75e-11 ],
            [ 1.52, 0.40212174244,  4.42e-11 ],
            [ 1.53, 0.39349840117,  4.13e-11 ],
            [ 1.54, 0.38519184421,  3.85e-11 ],
            [ 1.55, 0.37718514107,  3.6e-11 ],
            [ 1.56, 0.36946255256,  3.37e-11 ],
            [ 1.57, 0.36200942720,  3.16e-11 ],
            [ 1.58, 0.3548121084,   4.71e-10 ],
            [ 1.59, 0.3478578510,   4.43e-10 ],
            [ 1.6,  0.3411347458,   4.16e-10 ],
            [ 1.61, 0.3346316524,   3.92e-10 ],
            [ 1.62, 0.3283381372,   3.69e-10 ],
            [ 1.63, 0.3222444181,   3.49e-10 ],
            [ 1.64, 0.3163413144,   3.29e-10 ],
            [ 1.65, 0.3106202007,   3.11e-10 ],
            [ 1.66, 0.3050729656,   2.94e-10 ],
            [ 1.67, 0.2996919737,   2.79e-10 ],
            [ 1.68, 0.2944700309,   2.64e-10 ],
            [ 1.69, 0.2894003530,   2.51e-10 ],
            [ 1.7,  0.2844765367,   2.38e-10 ],
            [ 1.71, 0.2796925330,   2.26e-10 ],
            [ 1.72, 0.2750426231,   2.15e-10 ],
            [ 1.73, 0.2705213959,   2.04e-10 ],
            [ 1.74, 0.2661237275,   1.95e-10 ],
            [ 1.75, 0.2618447627,   1.85e-10 ],
            [ 1.76, 0.2576798968,   1.77e-10 ],
            [ 1.77, 0.2536247603,   1.69e-10 ],
            [ 1.78, 0.2496752035,   1.61e-10 ],
            [ 1.79, 0.2458272831,   1.54e-10 ],
            [ 1.8,  0.2420772493,   1.47e-10 ],
            [ 1.81, 0.2384215343,   1.41e-10 ],
            [ 1.82, 0.2348567409,   1.34e-10 ],
            [ 1.83, 0.2313796328,   1.29e-10 ],
            [ 1.84, 0.2279871249,   1.23e-10 ],
            [ 1.85, 0.2246762746,   1.18e-10 ],
            [ 1.86, 0.2214442738,   1.13e-10 ],
            [ 1.87, 0.2182884411,   1.09e-10 ],
            [ 1.88, 0.2152062147,   1.04e-10 ],
            [ 1.89, 0.2121951459,   1e-10 ],
            [ 1.9,  0.2092528928,   9.64e-11 ],
            [ 1.91, 0.2063772147,   9.27e-11 ],
            [ 1.92, 0.2035659664,   8.92e-11 ],
            [ 1.93, 0.2008170935,   8.58e-11 ],
            [ 1.94, 0.1981286272,   8.27e-11 ],
            [ 1.95, 0.1954986804,   7.96e-11 ],
            [ 1.96, 0.1929254430,   7.67e-11 ],
            [ 1.97, 0.1904071784,   7.4e-11 ],
            [ 1.98, 0.1879422195,   7.14e-11 ],
            [ 1.99, 0.1855289655,   6.89e-11 ],
            [ 2,    0.183165869,    4.65e-09 ],
            [ 2.02, 0.1785843495,   6.21e-11 ],
            [ 2.04, 0.1741864757,   5.81e-11 ],
            [ 2.06, 0.1699619278,   5.44e-11 ],
            [ 2.08, 0.1659011292,   5.1e-11 ],
            [ 2.1,  0.16199518876,  4.8e-11 ],
            [ 2.12, 0.15823584020,  4.52e-11 ],
            [ 2.14, 0.15461538845,  4.26e-11 ],
            [ 2.16, 0.15112666112,  4.02e-11 ],
            [ 2.18, 0.14776296484,  3.81e-11 ],
            [ 2.2,  0.14451804607,  3.61e-11 ],
            [ 2.22, 0.14138605553,  3.42e-11 ],
            [ 2.24, 0.13836151621,  3.26e-11 ],
            [ 2.26, 0.1354392942,   4.9e-10 ],
            [ 2.28, 0.1326145723,   4.67e-10 ],
            [ 2.3,  0.1298828263,   4.46e-10 ],
            [ 2.32, 0.1272398025,   4.26e-10 ],
            [ 2.34, 0.1246814986,   4.08e-10 ],
            [ 2.36, 0.1222041450,   3.91e-10 ],
            [ 2.38, 0.1198041883,   3.75e-10 ],
            [ 2.4,  0.1174782761,   3.61e-10 ],
            [ 2.42, 0.1152232432,   3.47e-10 ],
            [ 2.44, 0.1130360985,   3.35e-10 ],
            [ 2.46, 0.1109140134,   3.23e-10 ],
            [ 2.48, 0.1088543111,   3.12e-10 ],
            [ 2.5,  0.1068544564,   3.02e-10 ],
            [ 2.52, 0.1049120466,   2.93e-10 ],
            [ 2.54, 0.1030248030,   2.84e-10 ],
            [ 2.56, 0.0980492077,   2.71e-10 ],
            [ 2.58, 0.0994072737,   2.68e-10 ],
            [ 2.6,  0.0946437156,   2.57e-10 ],
            [ 2.62, 0.0959858357,   2.55e-10 ],
            [ 2.64, 0.0943440666,   2.49e-10 ],
            [ 2.66, 0.0927459954,   2.43e-10 ],
            [ 2.68, 0.0911900214,   2.37e-10 ],
            [ 2.7,  0.0896746193,   2.32e-10 ],
            [ 2.72, 0.0881983346,   2.28e-10 ],
            [ 2.74, 0.0867597795,   2.23e-10 ],
            [ 2.76, 0.0853576296,   2.19e-10 ],
            [ 2.78, 0.0839906197,   2.15e-10 ],
            [ 2.8,  0.0826575411,   2.11e-10 ],
            [ 2.82, 0.0813572383,   2.08e-10 ],
            [ 2.84, 0.0800886062,   2.05e-10 ],
            [ 2.86, 0.0788505873,   2.02e-10 ],
            [ 2.88, 0.0776421693,   1.99e-10 ],
            [ 2.9,  0.0764623828,   1.96e-10 ],
            [ 2.92, 0.0753102987,   1.94e-10 ],
            [ 2.94, 0.0741850268,   1.91e-10 ],
            [ 2.96, 0.0730857133,   1.89e-10 ],
            [ 2.98, 0.0720115391,   1.87e-10 ],
            [ 3,    0.0709617180,   1.85e-10 ],
            [ 3.05, 0.0684388313,   1.81e-10 ],
            [ 3.1,  0.0660525126,   1.77e-10 ],
            [ 3.15, 0.0637927913,   1.74e-10 ],
            [ 3.2,  0.0616506137,   1.71e-10 ],
            [ 3.25, 0.0596177415,   4.59e-10 ],
            [ 3.3,  0.0576866629,   3.92e-10 ],
            [ 3.35, 0.0558505163,   3.33e-10 ],
            [ 3.4,  0.0541030218,   2.79e-10 ],
            [ 3.45, 0.0524384220,   2.31e-10 ],
            [ 3.5,  0.0508514298,   1.88e-10 ],
            [ 3.55, 0.0493371822,   1.49e-10 ],
            [ 3.6,  0.0478911989,   1.13e-10 ],
            [ 3.65, 0.0465093466,   8.15e-11 ],
            [ 3.7,  0.0451878062,   2.15e-10 ],
            [ 3.75, 0.0439230448,   1.9e-10 ],
            [ 3.8,  0.0427117888,   1.67e-10 ],
            [ 3.85, 0.0415510018,   1.46e-10 ],
            [ 3.9,  0.0404378639,   1.28e-10 ],
            [ 3.95, 0.0381714768,   1.1e-10 ],
            [ 4,    0.0383442266,   9.49e-11 ],
            [ 4.1,  0.0364119827,   6.79e-11 ],
            [ 4.2,  0.03462468813,  4.56e-11 ],
            [ 4.3,  0.03296797170,  2.71e-11 ],
            [ 4.4,  0.03142922556,  1.19e-11 ],
            [ 4.5,  0.02999734835,  7.53e-13 ],
            [ 4.6,  0.02866253158,  1.12e-11 ],
            [ 4.7,  0.02741608099,  1.98e-11 ],
            [ 4.8,  0.02625026624,  2.68e-11 ],
            [ 4.9,  0.02515819410,  3.26e-11 ],
            [ 5,    0.02413370080,  3.72e-11 ],
            [ 5.1,  0.02317126035,  4.1e-11 ],
            [ 5.2,  0.02226590620,  4.4e-11 ],
            [ 5.3,  0.02141316395,  4.63e-11 ],
            [ 5.4,  0.02060899342,  4.81e-11 ],
            [ 5.5,  0.01984973861,  4.95e-11 ],
            [ 5.6,  0.0191320843,   5.05e-11 ],
            [ 5.7,  0.0184530182,   5.11e-11 ],
            [ 5.8,  0.0178097983,   5.14e-11 ],
            [ 5.9,  0.0171999238,   5.16e-11 ],
            [ 6,    0.0166211100,   5.15e-11 ],
            [ 6.1,  0.0160712661,   5.13e-11 ],
            [ 6.2,  0.0155484757,   5.09e-11 ],
            [ 6.3,  0.01459804973,  4.99e-11 ],
            [ 6.4,  0.01457715945,  4.99e-11 ],
            [ 6.5,  0.01412552647,  4.92e-11 ],
            [ 6.6,  0.01369470705,  4.85e-11 ],
            [ 6.7,  0.01328343314,  4.78e-11 ],
            [ 6.8,  0.01289053243,  4.7e-11 ],
            [ 6.9,  0.01251491973,  4.62e-11 ],
            [ 7,    0.01215558929,  4.53e-11 ],
        ],
        [
            [ 1.01, 97.51488588255, 0 ],
            [ 1.02, 48.43243012129, 0 ],
            [ 1.03, 13.1943384920,  2.62e-10 ],
            [ 1.04, 9.9169080078,   4.98e-10 ],
            [ 1.05, 7.9486739753,   4.15e-10 ],
            [ 1.06, 6.6351271682,   2.83e-10 ],
            [ 1.07, 5.6957582344,   1.85e-10 ],
            [ 1.08, 4.9903081760,   1.22e-10 ],
            [ 1.09, 4.4408518079,   8.57e-11 ],
            [ 1.1,  4.000631113,    9e-10 ],
            [ 1.11, 3.639888560,    6.5e-10 ],
            [ 1.12, 3.338783748,    5.28e-10 ],
            [ 1.13, 3.0835792739,   4.01e-10 ],
            [ 1.14, 2.8644612137,   3.54e-10 ],
            [ 1.15, 2.6742314041,   3.35e-10 ],
            [ 1.16, 2.5074901076,   2.61e-10 ],
            [ 1.17, 2.3601071411,   2.56e-10 ],
            [ 1.18, 2.2288692896,   2.04e-10 ],
            [ 1.19, 2.1112390567,   1.65e-10 ],
            [ 1.2,  2.0051857851,   1.35e-10 ],
            [ 1.21, 1.9090650254,   1.11e-10 ],
            [ 1.22, 1.8215307996,   9.23e-11 ],
            [ 1.23, 1.7414707503,   7.73e-11 ],
            [ 1.24, 1.6679574997,   6.52e-11 ],
            [ 1.25, 1.600211680,    8.81e-10 ],
            [ 1.26, 1.537573493,    7.53e-10 ],
            [ 1.27, 1.479480589,    6.47e-10 ],
            [ 1.28, 1.425450685,    5.59e-10 ],
            [ 1.29, 1.3750677730,   4.85e-10 ],
            [ 1.3,  1.3279710933,   4.23e-10 ],
            [ 1.31, 1.2838462342,   3.7e-10 ],
            [ 1.32, 1.2424179027,   3.25e-10 ],
            [ 1.33, 1.2034440083,   2.87e-10 ],
            [ 1.34, 1.1667107900,   2.54e-10 ],
            [ 1.35, 1.1320287779,   2.25e-10 ],
            [ 1.36, 1.0992294273,   2.01e-10 ],
            [ 1.37, 1.0681622984,   1.8e-10 ],
            [ 1.38, 1.0386926805,   1.61e-10 ],
            [ 1.39, 1.0106995810,   1.44e-10 ],
            [ 1.4,  0.9840740169,   1.3e-10 ],
            [ 1.41, 0.9587175542,   1.17e-10 ],
            [ 1.42, 0.9345410571,   1.06e-10 ],
            [ 1.43, 0.9114636109,   9.64e-11 ],
            [ 1.44, 0.8894115913,   8.77e-11 ],
            [ 1.45, 0.8683178582,   7.99e-11 ],
            [ 1.46, 0.8481210535,   7.3e-11 ],
            [ 1.47, 0.8287649894,   6.68e-11 ],
            [ 1.48, 0.810198112,    9.8e-10 ],
            [ 1.49, 0.792373032,    9.02e-10 ],
            [ 1.5,  0.775246109,    8.34e-10 ],
            [ 1.51, 0.758777086,    7.72e-10 ],
            [ 1.52, 0.742928768,    7.18e-10 ],
            [ 1.53, 0.727666734,    6.69e-10 ],
            [ 1.54, 0.712959080,    6.26e-10 ],
            [ 1.55, 0.698776193,    5.89e-10 ],
            [ 1.56, 0.685090551,    5.56e-10 ],
            [ 1.57, 0.671876534,    5.27e-10 ],
            [ 1.58, 0.659110267,    5.02e-10 ],
            [ 1.59, 0.6467694722,   4.81e-10 ],
            [ 1.6,  0.6348333341,   4.63e-10 ],
            [ 1.61, 0.6232823840,   4.49e-10 ],
            [ 1.62, 0.6120983905,   4.38e-10 ],
            [ 1.63, 0.6012642627,   4.29e-10 ],
            [ 1.64, 0.5907639617,   4.23e-10 ],
            [ 1.65, 0.5805824203,   4.2e-10 ],
            [ 1.66, 0.5707054705,   4.19e-10 ],
            [ 1.67, 0.5611197766,   4.21e-10 ],
            [ 1.68, 0.5518127747,   4.25e-10 ],
            [ 1.69, 0.5427726175,   4.31e-10 ],
            [ 1.7,  0.5339881234,   4.38e-10 ],
            [ 1.71, 0.5254487305,   4.48e-10 ],
            [ 1.72, 0.5171444535,   4.6e-10 ],
            [ 1.73, 0.5090658448,   4.74e-10 ],
            [ 1.74, 0.5012039590,   4.84e-10 ],
            [ 1.75, 0.4935503193,   2.31e-10 ],
            [ 1.76, 0.4860968870,   2.32e-10 ],
            [ 1.77, 0.4788360337,   2.34e-10 ],
            [ 1.78, 0.4717605154,   2.67e-10 ],
            [ 1.79, 0.4648634479,   2.73e-10 ],
            [ 1.8,  0.4581382851,   2.81e-10 ],
            [ 1.81, 0.4515787982,   2.89e-10 ],
            [ 1.82, 0.4451790567,   2.99e-10 ],
            [ 1.83, 0.4389334103,   3.1e-10 ],
            [ 1.84, 0.4328364730,   3.21e-10 ],
            [ 1.85, 0.4268831074,   3.34e-10 ],
            [ 1.86, 0.4210684105,   3.47e-10 ],
            [ 1.87, 0.4153877004,   3.61e-10 ],
            [ 1.88, 0.4098365040,   3.76e-10 ],
            [ 1.89, 0.4044105452,   3.92e-10 ],
            [ 1.9,  0.3991057345,   4.09e-10 ],
            [ 1.91, 0.3939181586,   4.26e-10 ],
            [ 1.92, 0.3888440710,   4.44e-10 ],
            [ 1.93, 0.3838798831,   4.63e-10 ],
            [ 1.94, 0.3790221562,   4.83e-10 ],
            [ 1.95, 0.374267593,    5.03e-10 ],
            [ 1.96, 0.369613033,    5.24e-10 ],
            [ 1.97, 0.365055439,    5.45e-10 ],
            [ 1.98, 0.360591901,    5.67e-10 ],
            [ 1.99, 0.356219619,    5.9e-10 ],
            [ 2,    0.35193590,     6.51e-09 ],
            [ 2.02, 0.343623955,    6.6e-10 ],
            [ 2.04, 0.335636540,    7.1e-10 ],
            [ 2.06, 0.327955601,    7.61e-10 ],
            [ 2.08, 0.320564388,    8.14e-10 ],
            [ 2.1,  0.313447345,    8.68e-10 ],
            [ 2.12, 0.306590017,    9.61e-10 ],
            [ 2.14, 0.299978916,    7.14e-10 ],
            [ 2.16, 0.293601491,    5.23e-10 ],
            [ 2.18, 0.2874460198,   3.83e-10 ],
            [ 2.2,  0.2815015440,   2.89e-10 ],
            [ 2.22, 0.2757578120,   2.35e-10 ],
            [ 2.24, 0.2702052207,   2.18e-10 ],
            [ 2.26, 0.2648347655,   2.33e-10 ],
            [ 2.28, 0.2596379942,   2.78e-10 ],
            [ 2.3,  0.2546069652,   3.48e-10 ],
            [ 2.32, 0.2497342096,   4.41e-10 ],
            [ 2.34, 0.245012697,    5.55e-10 ],
            [ 2.36, 0.240435801,    6.86e-10 ],
            [ 2.38, 0.235997277,    8.33e-10 ],
            [ 2.4,  0.231691227,    9.94e-10 ],
            [ 2.42, 0.2275120666,   2.52e-10 ],
            [ 2.44, 0.2234545599,   2.61e-10 ],
            [ 2.46, 0.2195137108,   2.71e-10 ],
            [ 2.48, 0.2156848027,   2.81e-10 ],
            [ 2.5,  0.2119633666,   2.9e-10 ],
            [ 2.52, 0.2083451651,   3e-10 ],
            [ 2.54, 0.2048261776,   3.09e-10 ],
            [ 2.56, 0.2014025868,   3.18e-10 ],
            [ 2.58, 0.1980707656,   3.27e-10 ],
            [ 2.6,  0.1948272661,   3.36e-10 ],
            [ 2.62, 0.1916688079,   3.45e-10 ],
            [ 2.64, 0.1885922689,   3.54e-10 ],
            [ 2.66, 0.1855946754,   3.62e-10 ],
            [ 2.68, 0.1826731932,   3.71e-10 ],
            [ 2.7,  0.1798251204,   3.79e-10 ],
            [ 2.72, 0.1770478788,   3.87e-10 ],
            [ 2.74, 0.1743390076,   3.95e-10 ],
            [ 2.76, 0.1716961568,   4.03e-10 ],
            [ 2.78, 0.1691170804,   4.1e-10 ],
            [ 2.8,  0.1665996318,   4.18e-10 ],
            [ 2.82, 0.1641417573,   4.25e-10 ],
            [ 2.84, 0.1617414920,   4.32e-10 ],
            [ 2.86, 0.1593969546,   4.39e-10 ],
            [ 2.88, 0.1571063432,   4.46e-10 ],
            [ 2.9,  0.1548679312,   4.52e-10 ],
            [ 2.92, 0.1526800632,   4.59e-10 ],
            [ 2.94, 0.1505411519,   4.65e-10 ],
            [ 2.96, 0.1484496740,   4.71e-10 ],
            [ 2.98, 0.1464041676,   4.76e-10 ],
            [ 3,    0.1444032288,   4.82e-10 ],
            [ 3.05, 0.1395871438,   4.95e-10 ],
            [ 3.1,  0.135021638,    5.07e-10 ],
            [ 3.15, 0.130688852,    5.18e-10 ],
            [ 3.2,  0.126572545,    5.28e-10 ],
            [ 3.25, 0.122657919,    5.37e-10 ],
            [ 3.3,  0.118931459,    5.45e-10 ],
            [ 3.35, 0.115380804,    5.52e-10 ],
            [ 3.4,  0.111994623,    5.59e-10 ],
            [ 3.45, 0.108762515,    5.64e-10 ],
            [ 3.5,  0.105674910,    5.68e-10 ],
            [ 3.55, 0.102722998,    5.72e-10 ],
            [ 3.6,  0.099898647,    5.75e-10 ],
            [ 3.65, 0.097194345,    5.77e-10 ],
            [ 3.7,  0.094603141,    5.79e-10 ],
            [ 3.75, 0.092118594,    5.8e-10 ],
            [ 3.8,  0.089734729,    5.81e-10 ],
            [ 3.85, 0.087445995,    5.81e-10 ],
            [ 3.9,  0.085247229,    5.8e-10 ],
            [ 3.95, 0.083133623,    5.79e-10 ],
            [ 4,    0.081100695,    5.78e-10 ],
            [ 4.1,  0.077260411,    5.74e-10 ],
            [ 4.2,  0.073696077,    5.69e-10 ],
            [ 4.3,  0.070381137,    5.62e-10 ],
            [ 4.4,  0.067292225,    5.55e-10 ],
            [ 4.5,  0.0644087118,   3.58e-10 ],
            [ 4.6,  0.0617123166,   3.57e-10 ],
            [ 4.7,  0.0591867926,   3.55e-10 ],
            [ 4.8,  0.056817654,    8.01e-10 ],
            [ 4.9,  0.054591948,    7.88e-10 ],
            [ 5,    0.052498063,    7.74e-10 ],
            [ 5.1,  0.050525561,    7.6e-10 ],
            [ 5.2,  0.048665034,    7.45e-10 ],
            [ 5.3,  0.046907986,    7.3e-10 ],
            [ 5.4,  0.045246727,    7.16e-10 ],
            [ 5.5,  0.043674279,    7e-10 ],
            [ 5.6,  0.042184298,    6.85e-10 ],
            [ 5.7,  0.040771007,    6.7e-10 ],
            [ 5.8,  0.039429134,    6.55e-10 ],
            [ 5.9,  0.038153860,    6.41e-10 ],
            [ 6,    0.036940771,    6.26e-10 ],
            [ 6.1,  0.035785819,    6.11e-10 ],
            [ 6.2,  0.034685283,    5.97e-10 ],
            [ 6.3,  0.033635742,    5.83e-10 ],
            [ 6.4,  0.032634042,    5.69e-10 ],
            [ 6.5,  0.031677272,    5.56e-10 ],
            [ 6.6,  0.030762743,    5.42e-10 ],
            [ 6.7,  0.029887967,    5.29e-10 ],
            [ 6.8,  0.029050638,    5.17e-10 ],
            [ 6.9,  0.028248617,    5.04e-10 ],
            [ 7,    0.0274799189,   4.92e-10 ],
        ],
        [
            [ 1.01, 83.21242947121, 0 ],
            [ 1.02, 41.32847258068, 0 ],
            [ 1.03, 27.36972580397, 0 ],
            [ 1.04, 8.465143946,    9.94e-10 ],
            [ 1.05, 6.7861575847,   1.72e-10 ],
            [ 1.06, 5.6658024124,   1.7e-10 ],
            [ 1.07, 4.8647123879,   1.3e-10 ],
            [ 1.08, 4.2632021819,   9.15e-11 ],
            [ 1.09, 3.7947777877,   6.37e-11 ],
            [ 1.1,  3.419541003,    7.12e-10 ],
            [ 1.11, 3.112100548,    5.16e-10 ],
            [ 1.12, 2.8555276124,   4.22e-10 ],
            [ 1.13, 2.6381011464,   3.2e-10 ],
            [ 1.14, 2.4514480163,   2.83e-10 ],
            [ 1.15, 2.2894270949,   2.68e-10 ],
            [ 1.16, 2.1474318125,   2.08e-10 ],
            [ 1.17, 2.0219388622,   2.04e-10 ],
            [ 1.18, 1.9102073341,   1.62e-10 ],
            [ 1.19, 1.8100728561,   1.3e-10 ],
            [ 1.2,  1.7198034900,   1.06e-10 ],
            [ 1.21, 1.6379967981,   8.68e-11 ],
            [ 1.22, 1.5635049809,   7.18e-11 ],
            [ 1.23, 1.495379542,    9.54e-10 ],
            [ 1.24, 1.432829784,    8.02e-10 ],
            [ 1.25, 1.375191268,    6.78e-10 ],
            [ 1.26, 1.321901545,    5.77e-10 ],
            [ 1.27, 1.2724812870,   4.93e-10 ],
            [ 1.28, 1.2265194490,   4.24e-10 ],
            [ 1.29, 1.1836615126,   3.66e-10 ],
            [ 1.3,  1.1436000722,   3.18e-10 ],
            [ 1.31, 1.1060672451,   2.77e-10 ],
            [ 1.32, 1.0708285042,   2.42e-10 ],
            [ 1.33, 1.0376776309,   2.13e-10 ],
            [ 1.34, 1.0064325591,   1.88e-10 ],
            [ 1.35, 0.9769319301,   1.66e-10 ],
            [ 1.36, 0.9490322225,   1.47e-10 ],
            [ 1.37, 0.9226053460,   1.31e-10 ],
            [ 1.38, 0.8975366162,   1.17e-10 ],
            [ 1.39, 0.8737230393,   1.04e-10 ],
            [ 1.4,  0.8510718548,   9.34e-11 ],
            [ 1.41, 0.8294992907,   8.39e-11 ],
            [ 1.42, 0.8089294969,   7.55e-11 ],
            [ 1.43, 0.7892936267,   6.81e-11 ],
            [ 1.44, 0.770529044,    9.79e-10 ],
            [ 1.45, 0.752578634,    8.86e-10 ],
            [ 1.46, 0.735390207,    8.03e-10 ],
            [ 1.47, 0.718915978,    7.28e-10 ],
            [ 1.48, 0.703112106,    6.61e-10 ],
            [ 1.49, 0.687938294,    6.37e-10 ],
            [ 1.5,  0.673357440,    5.94e-10 ],
            [ 1.51, 0.659335321,    5.58e-10 ],
            [ 1.52, 0.645840319,    5.88e-10 ],
            [ 1.53, 0.632843178,    5.7e-10 ],
            [ 1.54, 0.620316783,    5.59e-10 ],
            [ 1.55, 0.608235971,    5.54e-10 ],
            [ 1.56, 0.596577355,    5.55e-10 ],
            [ 1.57, 0.585319167,    5.61e-10 ],
            [ 1.58, 0.574441126,    5.72e-10 ],
            [ 1.59, 0.563924304,    5.89e-10 ],
            [ 1.6,  0.553751023,    6.12e-10 ],
            [ 1.61, 0.543904746,    6.39e-10 ],
            [ 1.62, 0.534369991,    6.71e-10 ],
            [ 1.63, 0.525132243,    7.08e-10 ],
            [ 1.64, 0.516177883,    7.5e-10 ],
            [ 1.65, 0.507494118,    7.97e-10 ],
            [ 1.66, 0.499068919,    8.48e-10 ],
            [ 1.67, 0.490890965,    9.04e-10 ],
            [ 1.68, 0.482949591,    9.64e-10 ],
            [ 1.69, 0.4752347392,   9.47e-11 ],
            [ 1.7,  0.4677369197,   1.03e-10 ],
            [ 1.71, 0.4604471680,   1.13e-10 ],
            [ 1.72, 0.4533570098,   1.23e-10 ],
            [ 1.73, 0.4464584275,   1.33e-10 ],
            [ 1.74, 0.4397438302,   1.44e-10 ],
            [ 1.75, 0.4332060247,   1.56e-10 ],
            [ 1.76, 0.4268381904,   1.69e-10 ],
            [ 1.77, 0.4206338550,   1.82e-10 ],
            [ 1.78, 0.4145868722,   1.96e-10 ],
            [ 1.79, 0.4086914017,   2.1e-10 ],
            [ 1.8,  0.4029419103,   1.04e-10 ],
            [ 1.81, 0.3973330536,   2.41e-10 ],
            [ 1.82, 0.3918598612,   2.57e-10 ],
            [ 1.83, 0.3865175203,   2.74e-10 ],
            [ 1.84, 0.3813014621,   2.91e-10 ],
            [ 1.85, 0.3762073293,   3.09e-10 ],
            [ 1.86, 0.3712309631,   3.28e-10 ],
            [ 1.87, 0.3663683928,   3.47e-10 ],
            [ 1.88, 0.3616158244,   3.66e-10 ],
            [ 1.89, 0.3569696317,   3.86e-10 ],
            [ 1.9,  0.3524263461,   4.07e-10 ],
            [ 1.91, 0.3479826491,   4.28e-10 ],
            [ 1.92, 0.3436353633,   4.49e-10 ],
            [ 1.93, 0.3393814455,   4.71e-10 ],
            [ 1.94, 0.3352179797,   4.93e-10 ],
            [ 1.95, 0.331142170,    5.15e-10 ],
            [ 1.96, 0.327151335,    5.38e-10 ],
            [ 1.97, 0.323242901,    5.61e-10 ],
            [ 1.98, 0.319414399,    5.85e-10 ],
            [ 1.99, 0.315663456,    6.1e-10 ],
            [ 2,    0.31198787,     1.08e-08 ],
            [ 2.02, 0.304853627,    6.82e-10 ],
            [ 2.04, 0.297995363,    7.35e-10 ],
            [ 2.06, 0.2913976834,   4.05e-10 ],
            [ 2.08, 0.2850463862,   4.37e-10 ],
            [ 2.1,  0.278928280,    5.33e-10 ],
            [ 2.12, 0.273031100,    5.69e-10 ],
            [ 2.14, 0.267343421,    6.06e-10 ],
            [ 2.16, 0.261854593,    6.43e-10 ],
            [ 2.18, 0.256554673,    6.8e-10 ],
            [ 2.2,  0.251434370,    7.17e-10 ],
            [ 2.22, 0.246484987,    7.53e-10 ],
            [ 2.24, 0.241698379,    7.89e-10 ],
            [ 2.26, 0.237066908,    8.25e-10 ],
            [ 2.28, 0.232583404,    8.59e-10 ],
            [ 2.3,  0.228241130,    8.94e-10 ],
            [ 2.32, 0.224033749,    9.27e-10 ],
            [ 2.34, 0.219955295,    9.6e-10 ],
            [ 2.36, 0.216000148,    9.91e-10 ],
            [ 2.38, 0.2121630066,   1.79e-10 ],
            [ 2.4,  0.2084388659,   1.86e-10 ],
            [ 2.42, 0.2048229991,   1.92e-10 ],
            [ 2.44, 0.2013109367,   1.99e-10 ],
            [ 2.46, 0.1978984499,   2.05e-10 ],
            [ 2.48, 0.1945815344,   2.11e-10 ],
            [ 2.5,  0.1913563957,   2.16e-10 ],
            [ 2.52, 0.1882194356,   2.22e-10 ],
            [ 2.54, 0.1851672398,   2.27e-10 ],
            [ 2.56, 0.1821965658,   2.32e-10 ],
            [ 2.58, 0.1793043330,   2.37e-10 ],
            [ 2.6,  0.1764876120,   2.42e-10 ],
            [ 2.62, 0.1737436160,   2.46e-10 ],
            [ 2.64, 0.1710696919,   2.5e-10 ],
            [ 2.66, 0.1684633125,   2.53e-10 ],
            [ 2.68, 0.1659220691,   2.57e-10 ],
            [ 2.7,  0.1634436647,   2.6e-10 ],
            [ 2.72, 0.1610259077,   2.63e-10 ],
            [ 2.74, 0.1586667056,   2.65e-10 ],
            [ 2.76, 0.1563640601,   2.67e-10 ],
            [ 2.78, 0.1541160610,   2.69e-10 ],
            [ 2.8,  0.1519208822,   2.71e-10 ],
            [ 2.82, 0.1497767767,   2.72e-10 ],
            [ 2.84, 0.1476820727,   2.73e-10 ],
            [ 2.86, 0.1456351692,   2.74e-10 ],
            [ 2.88, 0.1436345326,   2.75e-10 ],
            [ 2.9,  0.1416786931,   2.75e-10 ],
            [ 2.92, 0.1397662417,   2.75e-10 ],
            [ 2.94, 0.1378958264,   2.74e-10 ],
            [ 2.96, 0.1360661503,   2.74e-10 ],
            [ 2.98, 0.1342759680,   2.73e-10 ],
            [ 3,    0.1325240837,   2.71e-10 ],
            [ 3.05, 0.1283044957,   2.67e-10 ],
            [ 3.1,  0.1243004359,   2.62e-10 ],
            [ 3.15, 0.1204966777,   2.55e-10 ],
            [ 3.2,  0.1168793690,   2.46e-10 ],
            [ 3.25, 0.1134358819,   2.36e-10 ],
            [ 3.3,  0.1101546808,   2.26e-10 ],
            [ 3.35, 0.107025208,    9.68e-10 ],
            [ 3.4,  0.104037783,    9.05e-10 ],
            [ 3.45, 0.101183514,    8.38e-10 ],
            [ 3.5,  0.098454220,    7.7e-10 ],
            [ 3.55, 0.095842359,    6.99e-10 ],
            [ 3.6,  0.093340974,    6.28e-10 ],
            [ 3.65, 0.090943630,    5.54e-10 ],
            [ 3.7,  0.0886443733,   4.81e-10 ],
            [ 3.75, 0.0864376823,   4.06e-10 ],
            [ 3.8,  0.0843184331,   3.31e-10 ],
            [ 3.85, 0.082281871,    9.33e-10 ],
            [ 3.9,  0.0803235458,   4.09e-10 ],
            [ 3.95, 0.078439335,    5.94e-10 ],
            [ 4,    0.0766253961,   2.92e-10 ],
            [ 4.1,  0.0731941674,   1.72e-10 ],
            [ 4.2,  0.070003788,    8.54e-10 ],
            [ 4.3,  0.0670313894,   3.87e-10 ],
            [ 4.4,  0.064256832,    5.16e-10 ],
            [ 4.5,  0.061662314,    6.39e-10 ],
            [ 4.6,  0.059232050,    7.54e-10 ],
            [ 4.7,  0.056951992,    8.62e-10 ],
            [ 4.8,  0.054809604,    9.63e-10 ],
            [ 4.9,  0.0527936671,   2.78e-10 ],
            [ 5,    0.0508941115,   3.01e-10 ],
            [ 5.1,  0.0491018776,   3.23e-10 ],
            [ 5.2,  0.0474087939,   3.43e-10 ],
            [ 5.3,  0.0458074727,   3.6e-10 ],
            [ 5.4,  0.0442912202,   3.76e-10 ],
            [ 5.5,  0.0428539588,   3.89e-10 ],
            [ 5.6,  0.0414901592,   4e-10 ],
            [ 5.7,  0.0401947812,   4.08e-10 ],
            [ 5.8,  0.0389632226,   4.13e-10 ],
            [ 5.9,  0.0377912739,   4.15e-10 ],
            [ 6,    0.0366750783,   4.14e-10 ],
            [ 6.1,  0.0356110976,   4.08e-10 ],
            [ 6.2,  0.0345960808,   3.98e-10 ],
            [ 6.3,  0.0336270381,   3.81e-10 ],
            [ 6.4,  0.0327012166,   3.55e-10 ],
            [ 6.5,  0.0318160801,   3.21e-10 ],
            [ 6.6,  0.0309692916,   1.21e-10 ],
            [ 6.7,  0.0301587005,   2.25e-10 ],
            [ 6.8,  0.0293823307,   2.48e-10 ],
            [ 6.9,  0.0286383831,   3.49e-10 ],
            [ 7,    0.027925269,    8.44e-10 ],
        ],
    ];
}

sub alpha_interpolate {
    my ( $sym, $gamma ) = @_;

    # Given: a 1d symmetry (0, 1, or 2) and an ideal gas gamma
    # return: parameter alpha for the similarity solution
    # returns undef if out of bounds of table
    # (currently gamma<1.1 || $gamma>7 )

    # alpha is obtained by interpolating a table of pre-computed values.
    # The table is spaced closely enough that cubic interpolation of the
    # interpolated values gives comparable accuracy to the tabulated values
    # (error below about 1.e-7 over most of the range).

    return if ( $sym != 0 && $sym != 1 && $sym != 2 );

    my $rtab = $ralpha_table->[$sym];
    my $ntab = @{$rtab};
    my ( $jl, $ju );
    my $icol = 0;

    # A small tolerance to avoid interpolations
    my $eps = 1.e-6;

    ( $jl, $ju ) = locate_2d( $gamma, $icol, $rtab, $jl, $ju );
    my ( $gamma_min, $alpha_min ) = @{ $rtab->[0] };
    my ( $gamma_max, $alpha_max ) = @{ $rtab->[1] };
    if ( $jl < 0 ) {
        return if ( $gamma + $eps < $gamma_min );
        return $alpha_min;
    }
    if ( $ju >= $ntab ) {
        return if ( $gamma - $eps > $gamma_max );
        return $alpha_max;
    }

    # Define N consecutive lagrange interpolation points;
    # Using 4 points gives sufficient accuracy
    my $NLAG = 4;
    my $rj_interp = set_interpolation_points( $jl, $ntab, $NLAG );

    my ( $rx, $ry );

    # alpha varies approximately as 1/{ (gamma-1)*sqrt(gamma+1) },
    # so we can improve accuracy by interpolating the function
    # alpha*(gamma-1)*sqrt(gamma+1)
    foreach my $jj ( @{$rj_interp} ) {
        my ( $xx, $yy ) = @{ $rtab->[$jj] };
        push @{$rx}, $xx;
        push @{$ry}, $yy * ( $xx - 1 ) * sqrt( $xx + 1 );
    }

    my $ff = polint( $gamma, $rx, $ry );
    my $alpha = $ff / ( ( $gamma - 1 ) * sqrt( $gamma + 1 ) );

    return ($alpha);
}
1;
