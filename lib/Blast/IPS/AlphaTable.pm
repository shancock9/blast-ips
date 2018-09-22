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
    #   err = the estimated integration error in the value of alpha
    #         the table is printed with up to 12 decimal digits

    $ralpha_table = [
        [
            [ 1.01, 22.28576216759, 1.24e-11 ],
            [ 1.02, 11.1702575237,  1.78e-10 ],
            [ 1.03, 7.4623970604,   9.98e-11 ],
            [ 1.04, 5.606607857586, 1.41e-13 ],
            [ 1.05, 4.491763108289, 4e-14 ],
            [ 1.06, 3.747478215300, 1.04e-13 ],
            [ 1.07, 3.215010091881, 8.79e-14 ],
            [ 1.08, 2.814981882407, 2.44e-14 ],
            [ 1.09, 2.503291124368, 2.11e-13 ],
            [ 1.1,  2.253472953549, 1.39e-13 ],
            [ 1.11, 2.048683333022, 1.19e-13 ],
            [ 1.12, 1.877690631556, 9.64e-14 ],
            [ 1.13, 1.732717214216, 6.13e-14 ],
            [ 1.14, 1.608206042757, 5.86e-14 ],
            [ 1.15, 1.500080617477, 3.51e-14 ],
            [ 1.16, 1.405282425967, 4.49e-14 ],
            [ 1.17, 1.321471633021, 4.04e-14 ],
            [ 1.18, 1.246827530548, 3.86e-14 ],
            [ 1.19, 1.179911994969, 3.38e-14 ],
            [ 1.2,  1.119573900472, 2.26e-14 ],
            [ 1.21, 1.064880836969, 1.47e-14 ],
            [ 1.22, 1.015069445546, 1.29e-14 ],
            [ 1.23, 0.969508705766, 2.22e-13 ],
            [ 1.24, 0.927672397650, 1.87e-13 ],
            [ 1.25, 0.889118169795, 1.65e-13 ],
            [ 1.26, 0.853471435395, 1.43e-13 ],
            [ 1.27, 0.820412844721, 1.25e-13 ],
            [ 1.28, 0.789668440199, 1.06e-13 ],
            [ 1.29, 0.761001846735, 9.48e-14 ],
            [ 1.3,  0.734208022570, 8.38e-14 ],
            [ 1.31, 0.709108218419, 7.64e-14 ],
            [ 1.32, 0.685545880708, 6.94e-14 ],
            [ 1.33, 0.663383298737, 5.77e-14 ],
            [ 1.34, 0.642498842709, 5.34e-14 ],
            [ 1.35, 0.622784674521, 4.96e-14 ],
            [ 1.36, 0.604144839455, 4.67e-14 ],
            [ 1.37, 0.586493666769, 4.13e-14 ],
            [ 1.38, 0.569754422335, 3.34e-14 ],
            [ 1.39, 0.553858168121, 3.2e-14 ],
            [ 1.4,  0.538742792368, 2.99e-14 ],
            [ 1.41, 0.524352181345, 2.64e-14 ],
            [ 1.42, 0.510635509121, 2.38e-14 ],
            [ 1.43, 0.497546626155, 2.15e-14 ],
            [ 1.44, 0.485043531025, 2.09e-14 ],
            [ 1.45, 0.473087912366, 1.71e-14 ],
            [ 1.46, 0.461644750361, 1.54e-14 ],
            [ 1.47, 0.450681968926, 1.41e-14 ],
            [ 1.48, 0.440170131216, 2.34e-13 ],
            [ 1.49, 0.430082172265, 2.18e-13 ],
            [ 1.5,  0.420393163590, 2.02e-13 ],
            [ 1.51, 0.411080105357, 1.88e-13 ],
            [ 1.52, 0.402121742435, 1.75e-13 ],
            [ 1.53, 0.393498401168, 1.62e-13 ],
            [ 1.54, 0.385191844203, 1.51e-13 ],
            [ 1.55, 0.377185141069, 1.43e-13 ],
            [ 1.56, 0.369462552556, 1.32e-13 ],
            [ 1.57, 0.362009427195, 1.24e-13 ],
            [ 1.58, 0.354812108387, 1.17e-13 ],
            [ 1.59, 0.347857850924, 1.09e-13 ],
            [ 1.6,  0.341134745805, 1.02e-13 ],
            [ 1.61, 0.334631652405, 9.7e-14 ],
            [ 1.62, 0.328338137163, 9.2e-14 ],
            [ 1.63, 0.322244418070, 8.64e-14 ],
            [ 1.64, 0.316341314331, 8.08e-14 ],
            [ 1.65, 0.310620200635, 7.63e-14 ],
            [ 1.66, 0.305072965558, 7.19e-14 ],
            [ 1.67, 0.299691973657, 6.87e-14 ],
            [ 1.68, 0.294470030894, 6.52e-14 ],
            [ 1.69, 0.289400353031, 6.17e-14 ],
            [ 1.7,  0.284476536720, 5.92e-14 ],
            [ 1.71, 0.279692533017, 5.57e-14 ],
            [ 1.72, 0.275042623079, 5.42e-14 ],
            [ 1.73, 0.270521395847, 5.11e-14 ],
            [ 1.74, 0.266123727523, 4.81e-14 ],
            [ 1.75, 0.261844762667, 4.57e-14 ],
            [ 1.76, 0.257679896783, 4.42e-14 ],
            [ 1.77, 0.253624760244, 4.22e-14 ],
            [ 1.78, 0.249675203440, 4e-14 ],
            [ 1.79, 0.245827283048, 3.82e-14 ],
            [ 1.8,  0.242077249322, 3.65e-14 ],
            [ 1.81, 0.238421534307, 3.47e-14 ],
            [ 1.82, 0.234856740918, 3.21e-14 ],
            [ 1.83, 0.231379632794, 3.16e-14 ],
            [ 1.84, 0.227987124870, 3.08e-14 ],
            [ 1.85, 0.224676274605, 2.95e-14 ],
            [ 1.86, 0.221444273825, 2.7e-14 ],
            [ 1.87, 0.218288441110, 2.56e-14 ],
            [ 1.88, 0.215206214705, 2.53e-14 ],
            [ 1.89, 0.212195145893, 2.44e-14 ],
            [ 1.9,  0.209252892818, 2.36e-14 ],
            [ 1.91, 0.206377214699, 2.24e-14 ],
            [ 1.92, 0.203565966421, 2.12e-14 ],
            [ 1.93, 0.200817093476, 2.13e-14 ],
            [ 1.94, 0.198128627213, 2.04e-14 ],
            [ 1.95, 0.195498680395, 1.88e-14 ],
            [ 1.96, 0.192925443022, 1.87e-14 ],
            [ 1.97, 0.190407178413, 1.8e-14 ],
            [ 1.98, 0.187942219526, 1.76e-14 ],
            [ 1.99, 0.185528965498, 1.73e-14 ],
            [ 2,    0.1831658778,   1.53e-10 ],
            [ 2.02, 0.178584349525, 2.4e-13 ],
            [ 2.04, 0.174186475715, 2.22e-13 ],
            [ 2.06, 0.169961927797, 2.07e-13 ],
            [ 2.08, 0.165901129231, 1.93e-13 ],
            [ 2.1,  0.161995188760, 1.8e-13 ],
            [ 2.12, 0.158235840193, 1.68e-13 ],
            [ 2.14, 0.154615388450, 1.57e-13 ],
            [ 2.16, 0.151126661115, 1.47e-13 ],
            [ 2.18, 0.147762964840, 1.38e-13 ],
            [ 2.2,  0.144518046064, 1.28e-13 ],
            [ 2.22, 0.141386055531, 1.2e-13 ],
            [ 2.24, 0.138361516203, 1.13e-13 ],
            [ 2.26, 0.135439294180, 1.05e-13 ],
            [ 2.28, 0.132614572315, 9.88e-14 ],
            [ 2.3,  0.129882826225, 9.22e-14 ],
            [ 2.32, 0.127239802452, 8.62e-14 ],
            [ 2.34, 0.124681498559, 8.11e-14 ],
            [ 2.36, 0.122204144951, 7.55e-14 ],
            [ 2.38, 0.119804188254, 7.07e-14 ],
            [ 2.4,  0.117478276102, 6.59e-14 ],
            [ 2.42, 0.115223243191, 6.17e-14 ],
            [ 2.44, 0.113036098474, 5.75e-14 ],
            [ 2.46, 0.110914013402, 5.4e-14 ],
            [ 2.48, 0.108854311094, 5.04e-14 ],
            [ 2.5,  0.106854456372, 4.66e-14 ],
            [ 2.52, 0.104912046561, 4.35e-14 ],
            [ 2.54, 0.103024803007, 4.06e-14 ],
            [ 2.56, 0.101190563227, 3.75e-14 ],
            [ 2.58, 0.099407273658, 3.47e-14 ],
            [ 2.6,  0.097672982930, 3.21e-14 ],
            [ 2.62, 0.095985835641, 2.96e-14 ],
            [ 2.64, 0.094344066576, 2.73e-14 ],
            [ 2.66, 0.092745995339, 2.5e-14 ],
            [ 2.68, 0.091190021363, 2.29e-14 ],
            [ 2.7,  0.089674619270, 2.05e-14 ],
            [ 2.72, 0.088198334542, 1.87e-14 ],
            [ 2.74, 0.086759779502, 1.67e-14 ],
            [ 2.76, 0.085357629548, 1.5e-14 ],
            [ 2.78, 0.083990619649, 1.32e-14 ],
            [ 2.8,  0.082657541067, 1.17e-14 ],
            [ 2.82, 0.081357238296, 1.02e-14 ],
            [ 2.84, 0.080088606194, 2.33e-13 ],
            [ 2.86, 0.078850587297, 2.13e-13 ],
            [ 2.88, 0.077642169310, 1.94e-13 ],
            [ 2.9,  0.076462382744, 1.75e-13 ],
            [ 2.92, 0.075310298707, 1.58e-13 ],
            [ 2.94, 0.074185026827, 1.41e-13 ],
            [ 2.96, 0.073085713298, 1.25e-13 ],
            [ 2.98, 0.072011539050, 1.1e-13 ],
            [ 3,    0.070961718021, 9.52e-14 ],
            [ 3.05, 0.068438831313, 6.18e-14 ],
            [ 3.1,  0.066052512606, 3.2e-14 ],
            [ 3.15, 0.063792791291, 5.7e-15 ],
            [ 3.2,  0.061650613654, 1.78e-14 ],
            [ 3.25, 0.059617741307, 3.88e-14 ],
            [ 3.3,  0.057686662751, 5.73e-14 ],
            [ 3.35, 0.055850516133, 7.37e-14 ],
            [ 3.4,  0.054103021589, 9.32e-14 ],
            [ 3.45, 0.052438421800, 9.81e-14 ],
            [ 3.5,  0.050851429648, 1.13e-13 ],
            [ 3.55, 0.049337181981, 1.23e-13 ],
            [ 3.6,  0.047891198693, 1.32e-13 ],
            [ 3.65, 0.046509346409, 1.4e-13 ],
            [ 3.7,  0.045187806193, 1.47e-13 ],
            [ 3.75, 0.043923044760, 1.53e-13 ],
            [ 3.8,  0.042711788774, 1.58e-13 ],
            [ 3.85, 0.041551001836, 1.63e-13 ],
            [ 3.9,  0.040437863853, 1.67e-13 ],
            [ 3.95, 0.039369752510, 1.7e-13 ],
            [ 4,    0.038344226579, 1.72e-13 ],
            [ 4.1,  0.036411982753, 1.77e-13 ],
            [ 4.2,  0.034624688139, 1.79e-13 ],
            [ 4.3,  0.032967971705, 1.8e-13 ],
            [ 4.4,  0.031429225568, 1.81e-13 ],
            [ 4.5,  0.029997348356, 1.8e-13 ],
            [ 4.6,  0.028662531591, 1.78e-13 ],
            [ 4.7,  0.027416080995, 1.76e-13 ],
            [ 4.8,  0.026250266248, 1.74e-13 ],
            [ 4.9,  0.025158194112, 1.71e-13 ],
            [ 5,    0.024133700811, 1.68e-13 ],
            [ 5.1,  0.023171260364, 1.65e-13 ],
            [ 5.2,  0.022265906215, 1.62e-13 ],
            [ 5.3,  0.021413163961, 1.58e-13 ],
            [ 5.4,  0.020608993431, 1.55e-13 ],
            [ 5.5,  0.019849738617, 1.51e-13 ],
            [ 5.6,  0.019132084278, 1.47e-13 ],
            [ 5.7,  0.018453018203, 1.44e-13 ],
            [ 5.8,  0.017809798293, 1.4e-13 ],
            [ 5.9,  0.017199923780, 1.36e-13 ],
            [ 6,    0.016621109981, 1.33e-13 ],
            [ 6.1,  0.016071266108, 1.29e-13 ],
            [ 6.2,  0.015548475706, 1.26e-13 ],
            [ 6.3,  0.015050979371, 1.23e-13 ],
            [ 6.4,  0.014577159453, 1.19e-13 ],
            [ 6.5,  0.014125526475, 1.16e-13 ],
            [ 6.6,  0.013694707054, 1.13e-13 ],
            [ 6.7,  0.013283433146, 1.1e-13 ],
            [ 6.8,  0.012890532434, 1.07e-13 ],
            [ 6.9,  0.012514919736, 1.04e-13 ],
            [ 7,    0.012155589298, 1.01e-13 ],
        ],
        [
            [ 1.01, 39.3847992396,  2.92e-10 ],
            [ 1.02, 19.7444723969,  9.19e-11 ],
            [ 1.03, 13.19433848888, 1.3e-11 ],
            [ 1.04, 9.916908007790, 0 ],
            [ 1.05, 7.948673975244, 1.85e-13 ],
            [ 1.06, 6.635127168142, 1.59e-13 ],
            [ 1.07, 5.695758234365, 8.44e-14 ],
            [ 1.08, 4.990308176019, 4.57e-13 ],
            [ 1.09, 4.440851807929, 2.76e-13 ],
            [ 1.1,  4.000631112457, 2.44e-13 ],
            [ 1.11, 3.639888559552, 1.5e-13 ],
            [ 1.12, 3.338783748016, 1.39e-13 ],
            [ 1.13, 3.083579273850, 1.3e-13 ],
            [ 1.14, 2.864461213671, 9.55e-14 ],
            [ 1.15, 2.674231404118, 5.95e-14 ],
            [ 1.16, 2.507490107540, 4.75e-14 ],
            [ 1.17, 2.360107141074, 6.26e-14 ],
            [ 1.18, 2.228869289627, 5.28e-14 ],
            [ 1.19, 2.111239056670, 3.73e-14 ],
            [ 1.2,  2.005185785122, 3.29e-14 ],
            [ 1.21, 1.909065025395, 4.23e-13 ],
            [ 1.22, 1.821530799595, 3.61e-13 ],
            [ 1.23, 1.741470750322, 3.07e-13 ],
            [ 1.24, 1.667957499666, 2.51e-13 ],
            [ 1.25, 1.600211679762, 2.21e-13 ],
            [ 1.26, 1.537573492782, 1.93e-13 ],
            [ 1.27, 1.479480589175, 1.61e-13 ],
            [ 1.28, 1.425450684766, 1.41e-13 ],
            [ 1.29, 1.375067772958, 1.19e-13 ],
            [ 1.3,  1.327971093308, 9.46e-14 ],
            [ 1.31, 1.283846234169, 9.3e-14 ],
            [ 1.32, 1.242417902641, 7.7e-14 ],
            [ 1.33, 1.203444008251, 7.04e-14 ],
            [ 1.34, 1.166710789955, 6.59e-14 ],
            [ 1.35, 1.132028777838, 6.13e-14 ],
            [ 1.36, 1.099229427283, 5.35e-14 ],
            [ 1.37, 1.068162298422, 4.82e-14 ],
            [ 1.38, 1.038692680454, 3.86e-14 ],
            [ 1.39, 1.010699581038, 3.33e-14 ],
            [ 1.4,  0.984074016880, 3.66e-14 ],
            [ 1.41, 0.958717554147, 4.62e-13 ],
            [ 1.42, 0.934541057066, 4.17e-13 ],
            [ 1.43, 0.911463610869, 3.79e-13 ],
            [ 1.44, 0.889411591342, 3.43e-13 ],
            [ 1.45, 0.868317858211, 3.17e-13 ],
            [ 1.46, 0.848121053520, 2.86e-13 ],
            [ 1.47, 0.828764989399, 2.61e-13 ],
            [ 1.48, 0.810198112172, 2.4e-13 ],
            [ 1.49, 0.792373031939, 2.2e-13 ],
            [ 1.5,  0.775246108457, 1.98e-13 ],
            [ 1.51, 0.758777085611, 1.82e-13 ],
            [ 1.52, 0.742928767940, 1.7e-13 ],
            [ 1.53, 0.727666733674, 1.57e-13 ],
            [ 1.54, 0.712959079551, 1.41e-13 ],
            [ 1.55, 0.698776193389, 1.34e-13 ],
            [ 1.56, 0.685090550924, 1.2e-13 ],
            [ 1.57, 0.671876533976, 1.13e-13 ],
            [ 1.58, 0.659110267333, 1.01e-13 ],
            [ 1.59, 0.646769472174, 9.4e-14 ],
            [ 1.6,  0.634833334077, 8.68e-14 ],
            [ 1.61, 0.623282383959, 8.09e-14 ],
            [ 1.62, 0.612098390491, 7.02e-14 ],
            [ 1.63, 0.601264262697, 6.66e-14 ],
            [ 1.64, 0.590763961657, 5.88e-14 ],
            [ 1.65, 0.580582420305, 5.17e-14 ],
            [ 1.66, 0.570705470495, 4.8e-14 ],
            [ 1.67, 0.561119776548, 4.45e-14 ],
            [ 1.68, 0.551812774642, 3.9e-14 ],
            [ 1.69, 0.542772617435, 3.23e-14 ],
            [ 1.7,  0.533988123411, 2.76e-14 ],
            [ 1.71, 0.525448730481, 4.81e-13 ],
            [ 1.72, 0.517144453420, 4.23e-13 ],
            [ 1.73, 0.509065844788, 3.68e-13 ],
            [ 1.74, 0.501203958992, 3.15e-13 ],
            [ 1.75, 0.493550319212, 2.66e-13 ],
            [ 1.76, 0.486096886909, 2.22e-13 ],
            [ 1.77, 0.478836033699, 1.77e-13 ],
            [ 1.78, 0.471760515376, 1.37e-13 ],
            [ 1.79, 0.464863447888, 1e-13 ],
            [ 1.8,  0.458138285105, 6.46e-14 ],
            [ 1.81, 0.451578798220, 3.51e-14 ],
            [ 1.82, 0.445179056642, 6.55e-15 ],
            [ 1.83, 0.438933410267, 1.83e-14 ],
            [ 1.84, 0.432836472994, 3.9e-14 ],
            [ 1.85, 0.426883107403, 5.77e-14 ],
            [ 1.86, 0.421068410492, 7.11e-14 ],
            [ 1.87, 0.415387700384, 8.29e-14 ],
            [ 1.88, 0.409836503937, 8.96e-14 ],
            [ 1.89, 0.404410545177, 9.33e-14 ],
            [ 1.9,  0.399105734502, 9.28e-14 ],
            [ 1.91, 0.393918158586, 8.8e-14 ],
            [ 1.92, 0.388844070933, 7.87e-14 ],
            [ 1.93, 0.383879883040, 6.65e-14 ],
            [ 1.94, 0.379022156118, 4.81e-14 ],
            [ 1.95, 0.374267593323, 2.64e-14 ],
            [ 1.96, 0.369613032479, 2.78e-16 ],
            [ 1.97, 0.365055439238, 3.15e-14 ],
            [ 1.98, 0.360591900655, 6.79e-14 ],
            [ 1.99, 0.356219619154, 1.1e-13 ],
            [ 2,    0.3519359069,   3.65e-10 ],
            [ 2.02, 0.343623954944, 2.61e-13 ],
            [ 2.04, 0.335636540064, 3.88e-13 ],
            [ 2.06, 0.327955600766, 6.13e-14 ],
            [ 2.08, 0.320564387483, 8.54e-14 ],
            [ 2.1,  0.313447345146, 1.12e-13 ],
            [ 2.12, 0.306590008240, 1.42e-13 ],
            [ 2.14, 0.299978906779, 1.76e-13 ],
            [ 2.16, 0.293601481877, 2.13e-13 ],
            [ 2.18, 0.287446009817, 2.55e-13 ],
            [ 2.2,  0.281501533602, 2.99e-13 ],
            [ 2.22, 0.275757801162, 3.45e-13 ],
            [ 2.24, 0.270205209459, 3.98e-13 ],
            [ 2.26, 0.264834753847, 4.52e-13 ],
            [ 2.28, 0.259637982120, 7.55e-14 ],
            [ 2.3,  0.254606952734, 8.29e-14 ],
            [ 2.32, 0.249734196788, 9.37e-14 ],
            [ 2.34, 0.245012683350, 1.05e-13 ],
            [ 2.36, 0.240435787814, 1.16e-13 ],
            [ 2.38, 0.235997262953, 1.28e-13 ],
            [ 2.4,  0.231691212435, 1.41e-13 ],
            [ 2.42, 0.227512066530, 1.55e-13 ],
            [ 2.44, 0.223454559822, 1.7e-13 ],
            [ 2.46, 0.219513710723, 1.83e-13 ],
            [ 2.48, 0.215684802623, 1.97e-13 ],
            [ 2.5,  0.211963366526, 2.14e-13 ],
            [ 2.52, 0.208345165038, 2.29e-13 ],
            [ 2.54, 0.204826177581, 2.45e-13 ],
            [ 2.56, 0.201402586732, 2.62e-13 ],
            [ 2.58, 0.198070765578, 2.8e-13 ],
            [ 2.6,  0.194827266011, 2.97e-13 ],
            [ 2.62, 0.191668807872, 3.16e-13 ],
            [ 2.64, 0.188592268882, 3.34e-13 ],
            [ 2.66, 0.185594675285, 3.54e-13 ],
            [ 2.68, 0.182673193151, 3.72e-13 ],
            [ 2.7,  0.179825120287, 3.91e-13 ],
            [ 2.72, 0.177047878698, 4.12e-13 ],
            [ 2.74, 0.174339007563, 4.3e-13 ],
            [ 2.76, 0.171696156685, 4.51e-13 ],
            [ 2.78, 0.169117080367, 4.71e-13 ],
            [ 2.8,  0.166599631703, 4.91e-13 ],
            [ 2.82, 0.164141757230, 8.7e-14 ],
            [ 2.84, 0.161741491924, 9.15e-14 ],
            [ 2.86, 0.159396954521, 9.61e-14 ],
            [ 2.88, 0.157106343123, 9.98e-14 ],
            [ 2.9,  0.154867931084, 1.03e-13 ],
            [ 2.92, 0.152680063147, 1.06e-13 ],
            [ 2.94, 0.150541151814, 1.12e-13 ],
            [ 2.96, 0.148449673943, 1.16e-13 ],
            [ 2.98, 0.146404167534, 1.19e-13 ],
            [ 3,    0.144403228720, 1.25e-13 ],
            [ 3.05, 0.139587143750, 1.35e-13 ],
            [ 3.1,  0.135021638057, 1.45e-13 ],
            [ 3.15, 0.130688852079, 1.55e-13 ],
            [ 3.2,  0.126572545356, 1.65e-13 ],
            [ 3.25, 0.122657918644, 1.75e-13 ],
            [ 3.3,  0.118931458916, 1.85e-13 ],
            [ 3.35, 0.115380803894, 1.94e-13 ],
            [ 3.4,  0.111994623297, 2.05e-13 ],
            [ 3.45, 0.108762514453, 2.13e-13 ],
            [ 3.5,  0.105674910278, 2.23e-13 ],
            [ 3.55, 0.102722997960, 2.31e-13 ],
            [ 3.6,  0.099898646915, 2.4e-13 ],
            [ 3.65, 0.097194344810, 2.48e-13 ],
            [ 3.7,  0.094603140617, 2.55e-13 ],
            [ 3.75, 0.092118593816, 2.62e-13 ],
            [ 3.8,  0.089734728990, 2.71e-13 ],
            [ 3.85, 0.087445995150, 2.77e-13 ],
            [ 3.9,  0.085247229231, 2.83e-13 ],
            [ 3.95, 0.083133623274, 2.9e-13 ],
            [ 4,    0.081100694851, 2.95e-13 ],
            [ 4.1,  0.077260411069, 3.06e-13 ],
            [ 4.2,  0.073696077196, 3.15e-13 ],
            [ 4.3,  0.070381136724, 3.23e-13 ],
            [ 4.4,  0.067292224992, 3.3e-13 ],
            [ 4.5,  0.064408711106, 3.36e-13 ],
            [ 4.6,  0.061712315941, 3.41e-13 ],
            [ 4.7,  0.059186791928, 3.45e-13 ],
            [ 4.8,  0.056817653328, 3.48e-13 ],
            [ 4.9,  0.054591947974, 3.5e-13 ],
            [ 5,    0.052498063280, 3.51e-13 ],
            [ 5.1,  0.050525560675, 3.52e-13 ],
            [ 5.2,  0.048665033758, 3.52e-13 ],
            [ 5.3,  0.046907986311, 3.52e-13 ],
            [ 5.4,  0.045246727032, 3.51e-13 ],
            [ 5.5,  0.043674278378, 3.49e-13 ],
            [ 5.6,  0.042184297387, 3.47e-13 ],
            [ 5.7,  0.040771006677, 3.45e-13 ],
            [ 5.8,  0.039429134158, 3.42e-13 ],
            [ 5.9,  0.038153860192, 3.39e-13 ],
            [ 6,    0.036940771166, 3.36e-13 ],
            [ 6.1,  0.035785818592, 3.33e-13 ],
            [ 6.2,  0.034685282992, 3.29e-13 ],
            [ 6.3,  0.033635741919, 3.25e-13 ],
            [ 6.4,  0.032634041599, 3.21e-13 ],
            [ 6.5,  0.031677271701, 3.17e-13 ],
            [ 6.6,  0.030762742870, 3.13e-13 ],
            [ 6.7,  0.029887966658, 3.09e-13 ],
            [ 6.8,  0.029050637575, 3.04e-13 ],
            [ 6.9,  0.028248617003, 3e-13 ],
            [ 7,    0.027479918747, 2.96e-13 ],
        ],
        [
            [ 1.01,         33.609111366107, 4.26e-14 ],
            [ 1.02,         16.8500645067,   2.66e-10 ],
            [ 1.03,         11.261271334324, 2.03e-13 ],
            [ 1.04,         8.465143317167,  9.25e-13 ],
            [ 1.05,         6.786157584719,  3.64e-14 ],
            [ 1.06,         5.665802412362,  1.25e-13 ],
            [ 1.07,         4.864712387895,  2.04e-14 ],
            [ 1.08,         4.263202181931,  3.77e-13 ],
            [ 1.09,         3.794777787671,  2.63e-13 ],
            [ 1.1,          3.419541003053,  1.72e-13 ],
            [ 1.11,         3.112100547746,  1.47e-13 ],
            [ 1.12,         2.855527612354,  9.81e-14 ],
            [ 1.13,         2.638101146388,  7.06e-14 ],
            [ 1.14,         2.451448016322,  5.46e-14 ],
            [ 1.15,         2.289427094928,  5.51e-14 ],
            [ 1.16,         2.147431812448,  5.73e-14 ],
            [ 1.17,         2.021938862183,  5.33e-14 ],
            [ 1.18,         1.910207334131,  5.11e-14 ],
            [ 1.19,         1.810072856103,  4.11e-14 ],
            [ 1.2,          1.719803489954,  4.19e-13 ],
            [ 1.21,         1.637996798096,  3.45e-13 ],
            [ 1.22,         1.563504980851,  2.86e-13 ],
            [ 1.23,         1.495379541498,  2.37e-13 ],
            [ 1.24,         1.432829783574,  2.03e-13 ],
            [ 1.25,         1.375191267517,  1.7e-13 ],
            [ 1.26,         1.321901545403,  1.45e-13 ],
            [ 1.27,         1.272481286935,  1.19e-13 ],
            [ 1.28,         1.226519448971,  1.08e-13 ],
            [ 1.29,         1.183661512623,  9.48e-14 ],
            [ 1.3,          1.143600072220,  7.95e-14 ],
            [ 1.31,         1.106067245131,  6.24e-14 ],
            [ 1.32,         1.070828504180,  6.24e-14 ],
            [ 1.33,         1.037677630924,  5.31e-14 ],
            [ 1.34,         1.006432559081,  4.09e-14 ],
            [ 1.35,         0.976931930102,  3.83e-14 ],
            [ 1.36,         0.949032222441,  3.43e-14 ],
            [ 1.37,         0.922605346025,  3.09e-14 ],
            [ 1.38,         0.897536616237,  4.61e-13 ],
            [ 1.39,         0.873723039325,  4.1e-13 ],
            [ 1.4,          0.851071854758,  3.7e-13 ],
            [ 1.41,         0.829499290670,  3.28e-13 ],
            [ 1.42,         0.808929496893,  2.98e-13 ],
            [ 1.43,         0.789293626694,  2.69e-13 ],
            [ 1.44,         0.770529043555,  2.43e-13 ],
            [ 1.45,         0.752578633564,  2.18e-13 ],
            [ 1.46,         0.735390207363,  1.97e-13 ],
            [ 1.47,         0.718915978304,  1.78e-13 ],
            [ 1.48,         0.703112105733,  1.63e-13 ],
            [ 1.49,         0.687938294082,  1.48e-13 ],
            [ 1.5,          0.673357439976,  1.34e-13 ],
            [ 1.51,         0.659335320773,  1.2e-13 ],
            [ 1.52,         0.645840318948,  1.03e-13 ],
            [ 1.53,         0.632843177620,  9.49e-14 ],
            [ 1.54,         0.620316783157,  8.44e-14 ],
            [ 1.55,         0.608235971443,  7.37e-14 ],
            [ 1.56,         0.596577354843,  6.43e-14 ],
            [ 1.57,         0.585319167327,  5.6e-14 ],
            [ 1.58,         0.574441125560,  4.79e-14 ],
            [ 1.59,         0.563924304084,  3.75e-14 ],
            [ 1.6,          0.553751022931,  3.02e-14 ],
            [ 1.61,         0.543904746248,  2.2e-14 ],
            [ 1.62,         0.534369990710,  4.72e-13 ],
            [ 1.63,         0.525132242604,  3.93e-13 ],
            [ 1.64,         0.516177882671,  3.22e-13 ],
            [ 1.65,         0.507494117842,  2.61e-13 ],
            [ 1.66,         0.499068919169,  2.07e-13 ],
            [ 1.67,         0.490890965272,  1.66e-13 ],
            [ 1.68,         0.482949590768,  1.36e-13 ],
            [ 1.69,         0.475234739152,  1.19e-13 ],
            [ 1.7,          0.467736919701,  1.15e-13 ],
            [ 1.71,         0.460447168002,  1.28e-13 ],
            [ 1.72,         0.453357009750,  1.53e-13 ],
            [ 1.73,         0.446458427507,  1.97e-13 ],
            [ 1.74,         0.439743830136,  2.61e-13 ],
            [ 1.75,         0.433206024673,  3.42e-13 ],
            [ 1.76,         0.426838190401,  4.45e-13 ],
            [ 1.77,         0.420633854935,  5.38e-14 ],
            [ 1.78,         0.414586872133,  7.47e-14 ],
            [ 1.79,         0.408691401675,  9.8e-14 ],
            [ 1.8,          0.402941890163,  1.28e-13 ],
            [ 1.81,         0.397333053610,  1.6e-13 ],
            [ 1.82,         0.391859861202,  1.96e-13 ],
            [ 1.83,         0.386517520230,  2.37e-13 ],
            [ 1.84,         0.381301462082,  2.84e-13 ],
            [ 1.85,         0.376207329224,  3.35e-13 ],
            [ 1.86,         0.371230963076,  3.92e-13 ],
            [ 1.87,         0.366368392725,  4.56e-13 ],
            [ 1.88,         0.361615824394,  7.37e-14 ],
            [ 1.89,         0.356969631618,  8.65e-14 ],
            [ 1.9,          0.352426346075,  9.87e-14 ],
            [ 1.91,         0.347982649003,  1.11e-13 ],
            [ 1.92,         0.343635363189,  1.29e-13 ],
            [ 1.93,         0.339381445457,  1.43e-13 ],
            [ 1.94,         0.335217979634,  1.6e-13 ],
            [ 1.95,         0.331142169963,  1.8e-13 ],
            [ 1.96,         0.327151334917,  1.98e-13 ],
            [ 1.97,         0.323242901393,  2.2e-13 ],
            [ 1.98,         0.319414399259,  2.44e-13 ],
            [ 1.99,         0.315663456230,  2.7e-13 ],
            [ 2,            0.3119877930,    1.57e-10 ],
            [ 2.02,         0.304853627340,  3.52e-13 ],
            [ 2.04,         0.297995362599,  4.13e-13 ],
            [ 2.06,         0.291397683242,  4.82e-13 ],
            [ 2.08,         0.285046386008,  9.28e-14 ],
            [ 2.1,          0.278928280364,  1.08e-13 ],
            [ 2.12,         0.273031099468,  1.23e-13 ],
            [ 2.14,         0.267343420407,  1.41e-13 ],
            [ 2.16,         0.261854592587,  1.58e-13 ],
            [ 2.18,         0.256554673325,  1.8e-13 ],
            [ 2.2,          0.251434369806,  2e-13 ],
            [ 2.22,         0.246484986681,  2.25e-13 ],
            [ 2.24,         0.241698378674,  2.47e-13 ],
            [ 2.26,         0.237066907648,  2.72e-13 ],
            [ 2.28,         0.232583403643,  3e-13 ],
            [ 2.3,          0.228241129467,  3.27e-13 ],
            [ 2.32,         0.224033748464,  3.56e-13 ],
            [ 2.34,         0.219955295129,  3.85e-13 ],
            [ 2.36,         0.216000148283,  4.15e-13 ],
            [ 2.38,         0.212163006546,  4.47e-13 ],
            [ 2.4,          0.208438865890,  4.79e-13 ],
            [ 2.42,         0.204822999056,  9.95e-14 ],
            [ 2.44,         0.201310936666,  1.05e-13 ],
            [ 2.46,         0.197898449863,  1.13e-13 ],
            [ 2.48,         0.194581534345,  1.22e-13 ],
            [ 2.5,          0.191356395648,  1.29e-13 ],
            [ 2.52,         0.188219435582,  1.36e-13 ],
            [ 2.54,         0.185167239709,  1.44e-13 ],
            [ 2.56,         0.182196565760,  1.49e-13 ],
            [ 2.58,         0.179304332933,  1.59e-13 ],
            [ 2.6,          0.176487611972,  1.68e-13 ],
            [ 2.62,         0.173743615970,  1.74e-13 ],
            [ 2.64,         0.171069691844,  1.82e-13 ],
            [ 2.66,         0.168463312405,  1.89e-13 ],
            [ 2.68,         0.165922068992,  1.99e-13 ],
            [ 2.7,          0.163443664623,  2.05e-13 ],
            [ 2.72,         0.161025907604,  2.12e-13 ],
            [ 2.74,         0.158666705583,  2.21e-13 ],
            [ 2.76,         0.156364059999,  2.27e-13 ],
            [ 2.78,         0.154116060902,  2.35e-13 ],
            [ 2.8,          0.151920882106,  2.4e-13 ],
            [ 2.82,         0.149776776666,  2.48e-13 ],
            [ 2.84,         0.147682072639,  2.54e-13 ],
            [ 2.86,         0.145635169120,  2.59e-13 ],
            [ 2.88,         0.143634532518,  2.65e-13 ],
            [ 2.9,          0.141678693073,  2.71e-13 ],
            [ 2.92,         0.139766241583,  2.75e-13 ],
            [ 2.94,         0.137895826331,  2.8e-13 ],
            [ 2.96,         0.136066150195,  2.86e-13 ],
            [ 2.98,         0.134275967938,  2.89e-13 ],
            [ 3,            0.132524083645,  2.93e-13 ],
            [ 3.05,         0.128304495613,  3.01e-13 ],
            [ 3.1,          0.124300435834,  3.06e-13 ],
            [ 3.15,         0.120496677615,  3.1e-13 ],
            [ 3.2,          0.116879368963,  3.11e-13 ],
            [ 3.25,         0.113435881837,  3.08e-13 ],
            [ 3.3,          0.110154680785,  3.05e-13 ],
            [ 3.35,         0.107025208115,  2.98e-13 ],
            [ 3.4,          0.104037783228,  2.88e-13 ],
            [ 3.45,         0.101183514108,  2.75e-13 ],
            [ 3.5,          0.098454219296,  2.61e-13 ],
            [ 3.55,         0.095842358920,  2.44e-13 ],
            [ 3.6,          0.093340973588,  2.25e-13 ],
            [ 3.65,         0.090943630108,  2.03e-13 ],
            [ 3.7,          0.088644373166,  1.8e-13 ],
            [ 3.75,         0.086437682222,  1.54e-13 ],
            [ 3.8,          0.084318432956,  1.26e-13 ],
            [ 3.85,         0.082281862744,  3.94e-13 ],
            [ 3.9,          0.080323539655,  2.67e-13 ],
            [ 3.95,         0.078439334573,  1.36e-13 ],
            [ 4,            0.076625396079,  4.8e-14 ],
            [ 4.1,          0.073194167803,  2.83e-13 ],
            [ 4.2,          0.070003788337,  1.47e-13 ],
            [ 4.3,          0.067031389568,  2.25e-13 ],
            [ 4.4,          0.064256831937,  3.06e-13 ],
            [ 4.5,          0.061662314608,  3.88e-13 ],
            [ 4.6,          0.059232050213,  4.71e-13 ],
            [ 4.7,          0.056951992073,  1.44e-13 ],
            [ 4.8,          0.054809604304,  1.67e-13 ],
            [ 4.9,          0.052793667174,  1.88e-13 ],
            [ 5,            0.050894111586,  2.1e-13 ],
            [ 5.1,          0.049101877754,  2.31e-13 ],
            [ 5.2,          0.047408794059,  2.51e-13 ],
            [ 5.3,          0.045807472827,  2.7e-13 ],
            [ 5.4,          0.044291220355,  2.89e-13 ],
            [ 5.5,          0.042853958967,  3.06e-13 ],
            [ 5.6,          0.041490159308,  3.21e-13 ],
            [ 5.7,          0.040194781339,  3.35e-13 ],
            [ 5.8,          0.038963222795,  3.46e-13 ],
            [ 5.9,          0.037791274043,  3.56e-13 ],
            [ 6,            0.036675078477,  3.62e-13 ],
            [ 6.1,          0.035611097715,  3.66e-13 ],
            [ 6.2,          0.034596080994,  3.65e-13 ],
            [ 6.3,          0.033627038279,  3.59e-13 ],
            [ 6.4,          0.032701216718,  3.46e-13 ],
            [ 6.5,          0.031816080193,  3.22e-13 ],
            [ 6.6,          0.030969291939,  2.85e-13 ],
            [ 6.7,          0.030158700569,  2.27e-13 ],
            [ 6.8,          0.029382330740,  1.46e-13 ],
            [ 6.9,          0.028638382976,  1.81e-13 ],
            [ 7, 0.027925268033,  2.78e-17 ],
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
