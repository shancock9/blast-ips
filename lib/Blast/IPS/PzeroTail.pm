package Blast::IPS::PzeroTail;

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
our $rpzero_tail;

#####################################
# Tables of the second zero in spherical shocks
#####################################

# returns a hash reference of the form:

# $rpzero_tail->{symmetry}->{gamma} = [ t0, r0, zlimit ]

# (t0, r0) = (time, range) that this zero overpressure point first appears
# zlimit = limiting value of z=r-ct at long range

BEGIN {

    # Created with: ./update_PT_table.pl
    # Thu May 10 18:41:11 2018   Inspiron-3668

    my $rtable = [

        # [symmetry, gamma, r0, t0, zlimit],
        [ 2, 1.1,   0.9308, 0.2244,  -0.61434, ],
        [ 2, 1.12,  0.9879, 0.2249,  -0.65418, ],
        [ 2, 1.15,  1.062,  0.2322,  -0.70671, ],
        [ 2, 1.17,  1.106,  0.2241,  -0.73819, ],
        [ 2, 1.2,   1.166,  0.2351,  -0.78128, ],
        [ 2, 1.23,  1.22,   0.2178,  -0.82051, ],
        [ 2, 1.25,  1.253,  0.2187,  -0.84494, ],
        [ 2, 1.3,   1.329,  0.2223,  -0.90104, ],
        [ 2, 1.35,  1.396,  0.2111,  -0.95153, ],
        [ 2, 1.4,   1.457,  0.1662,  -0.99775, ],
        [ 2, 1.45,  1.513,  0.2025,  -1.0402, ],
        [ 2, 1.5,   1.564,  0.1543,  -1.0797, ],
        [ 2, 1.55,  1.612,  0.1565,  -1.1168, ],
        [ 2, 1.6,   1.657,  0.1436,  -1.1517, ],
        [ 2, 1.65,  1.7,    0.14,    -1.1847, ],
        [ 2, 1.667, 1.714,  0.1745,  -1.1955, ],
        [ 2, 1.7,   1.74,   0.1491,  -1.2159, ],
        [ 2, 1.8,   1.815,  0.148,   -1.274, ],
        [ 2, 1.9,   1.884,  0.1212,  -1.3271, ],
        [ 2, 2,     1.947,  0.1209,  -1.3759, ],
        [ 2, 2.2,   2.062,  0.107,   -1.4627, ],
        [ 2, 2.5,   2.21,   0.1326,  -1.5718, ],
        [ 2, 3,     2.406,  0.08679, -1.7099, ],
        [ 2, 3.5,   2.552,  0.08777, -1.8067, ],
        [ 2, 4,     2.656,  0.05512, -1.8718, ],
        [ 2, 4.5,   2.723,  0.08787, -1.9114, ],
        [ 2, 5,     2.757,  0.06177, -1.9301, ],
        [ 2, 5.5,   2.763,  0.04805, -1.9318, ],
        [ 2, 6,     2.743,  0.05557, -1.9258, ],
        [ 2, 6.5,   2.696,  0.03321, -1.8772, ],
    ];
    foreach my $item ( @{$rtable} ) {
        my ( $symmetry, $gamma, @vars ) = @{$item}[ 0 .. 9 ];
        $rpzero_tail->{$symmetry}->{$gamma} = [@vars];
    }
}
1;
