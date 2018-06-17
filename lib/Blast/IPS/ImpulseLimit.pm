package Blast::IPS::ImpulseLimit;

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

use strict;
use warnings;
use 5.006;

use Carp;
our $rimpulse_limit;

# returns a hash reference to the fits of the form:
#     $rimpulse_limit->{symmetry}->{gamma} = [ Simp+, Simp- ]
# where Sint+ = positive phase integral of Sigma times dr
# where Sint- = negative phase integral of Sigma times dr

# Zero or undefined values occur if a value is not defined:
# Only spherical waves have a complete negative phase and an Simp- value
# Plane waves with gamma > 1.6 do not have complete positive phases

# These limiting spatial integral values can be converted into the limiting
# time integrals of overpressure times r**N/2 by multiplying by gamma

# For interpolation of these values in gamma, it is possible to use
# the fact that  I*gamma*alpha**0.55  is very slowly varying
# for more accuracy than simple interpolation.


BEGIN {

    # These values were taken from "blast_run_path_best.FD.csv"
    # which was created by 'blast_scan_all.pl FD'
    # which in turn pulled the values from the moc2 PT.txt files at about X=16
    my $rtable = [

        # [symmetry,gamma,Sint_plus,Sint_minus],
        [ 0, 1.1,   0.101566411374584,  0 ],
        [ 0, 1.12,  0.118016974505561,  0 ],
        [ 0, 1.15,  0.141688442954184,  0 ],
        [ 0, 1.17,  0.156616215963808,  0 ],
        [ 0, 1.2,   0.17697248399399,   0 ],
        [ 0, 1.23,  0.197566706851827,  0 ],
        [ 0, 1.25,  0.210058630384749,  0 ],
        [ 0, 1.3,   0.239460181985128,  0 ],
        [ 0, 1.35,  0.266364705610526,  0 ],
        [ 0, 1.4,   0.29113660806483,   0 ],
        [ 0, 1.45,  0.314151232018147,  0 ],
        [ 0, 1.5,   0.335712544243396,  0 ],
        [ 0, 1.55,  0.356052553581441,  0 ],
        [ 0, 1.6,   0.375384071066153,  0 ],
        [ 0, 1.65,  0.393802559124535,  0 ],
        [ 1, 1.1,   0.0273858517659239, 0 ],
        [ 1, 1.12,  0.0306276747609202, 0 ],
        [ 1, 1.15,  0.035017616943408,  0 ],
        [ 1, 1.17,  0.0376273357039186, 0 ],
        [ 1, 1.2,   0.0411973429492618, 0 ],
        [ 1, 1.23,  0.0443786370564763, 0 ],
        [ 1, 1.25,  0.0463289999799219, 0 ],
        [ 1, 1.3,   0.050681067067733,  0 ],
        [ 1, 1.35,  0.0544252315507507, 0 ],
        [ 1, 1.4,   0.0576853229100157, 0 ],
        [ 1, 1.45,  0.0605173143774102, 0 ],
        [ 1, 1.49,  0.0625475155701852, 0 ],
        [ 1, 1.5,   0.0631539875526601, 0 ],
        [ 1, 1.51,  0.0634874738053367, 0 ],
        [ 1, 1.55,  0.0651839742301522, 0 ],
        [ 1, 1.6,   0.0672247915300456, 0 ],
        [ 1, 1.65,  0.068998976686334,  0 ],
        [ 1, 1.667, 0.0695688671321038, 0 ],
        [ 1, 1.7,   0.0705862447829164, 0 ],
        [ 1, 1.8,   0.0733349782772008, 0 ],
        [ 1, 1.9,   0.0755879231167528, 0 ],
        [ 1, 2,     0.077440670910839,  0 ],
        [ 1, 2.2,   0.080213316969746,  0 ],
        [ 1, 2.5,   0.0830347579985831, 0 ],
        [ 1, 3,     0.085519639760284,  0 ],
        [ 1, 3.5,   0.0864622073504826, 0 ],
        [ 1, 4,     0.0866329771272463, 0 ],
        [ 1, 4.5,   0.0864267479790058, 0 ],
        [ 1, 5,     0.0859665282273764, 0 ],
        [ 1, 5.5,   0.0854319825777939, 0 ],
        [ 1, 6,     0.0847776605389496, 0 ],
        [ 1, 6.5,   0.0841343551150417, 0 ],
        [ 1, 7,     0.0834406488555697, 0 ],
        [ 2, 1.1,   0.0167394912158948, -0.0182137810684807 ],
        [ 2, 1.12,  0.0185087997462891, -0.0200930860197093 ],
        [ 2, 1.15,  0.0208314434012834, -0.0225431394633464 ],
        [ 2, 1.17,  0.0222023282168947, -0.0239802518710064 ],
        [ 2, 1.2,   0.0240440876867696, -0.0258993631654344 ],
        [ 2, 1.23,  0.025672464299312,  -0.027584927841135 ],
        [ 2, 1.25,  0.0266582246742109, -0.0286000352617623 ],
        [ 2, 1.3,   0.0288337102975662, -0.0308257048214101 ],
        [ 2, 1.35,  0.0306729090972275, -0.0326910279946982 ],
        [ 2, 1.4,   0.0322462814761683, -0.0342742530460385 ],
        [ 2, 1.45,  0.0336046578486096, -0.0356314299500579 ],
        [ 2, 1.5,   0.0347861004933462, -0.0368040571930541 ],
        [ 2, 1.55,  0.0358199660967535, -0.037823936984408 ],
        [ 2, 1.6,   0.0367292028831404, -0.0387157271925653 ],
        [ 2, 1.65,  0.0375321946906636, -0.0394990271276107 ],
        [ 2, 1.667, 0.0377837908970972, -0.0397435807417583 ],
        [ 2, 1.7,   0.0382439039607811, -0.0401896761533068 ],
        [ 2, 1.8,   0.0394403157460913, -0.0413422172570564 ],
        [ 2, 1.9,   0.0403939740336807, -0.0422520537527614 ],
        [ 2, 2,     0.0411587207596495, -0.0429749057718151 ],
        [ 2, 2.2,   0.0422692693582729, -0.0440093449978854 ],
        [ 2, 2.5,   0.0432373188501933, -0.0448865975485948 ],
        [ 2, 3,     0.04382792585655,   -0.0453763685804637 ],
        [ 2, 3.5,   0.0437839839782902, -0.0452776071482813 ],
        [ 2, 4,     0.0434409426470419, -0.0449110290573198 ],
        [ 2, 4.5,   0.0429503651833,    -0.0444195411735799 ],
        [ 2, 5,     0.04238702944095,   -0.0438721343885807 ],
        [ 2, 5.5,   0.0417896149821331, -0.0432963519629968 ],
        [ 2, 6,     0.0411786186443639, -0.042691725445229 ],
        [ 2, 6.5,   0.040563722549905,  -0.0413784053699492 ],
    ];
    foreach my $item ( @{$rtable} ) {
        my ( $symmetry, $gamma, @vars ) = @{$item}[ 0 .. 9 ];
        $rimpulse_limit->{$symmetry}->{$gamma} = [@vars];
    }
}
1;
