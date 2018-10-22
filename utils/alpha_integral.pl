#!/usr/bin/perl
use strict;
use warnings;

# This program creates a table of alpha values for the point source similarity
# solution.  Output values are in the table Blast::IPS::AlphaTable.pm.  It is
# saved for reference but should not need to be rerun.

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

use Blast::IPS::SimilaritySolution;

# Settings..

# itmax = max iterations in case the tolerance is not reached; 
my $itmax = 15;

# tol = stopping tolerance on absolute value of alpha
my $tol = 1.e-13;

# output file
my $ftable = "alpha_table_$itmax.txt";

my $fh_table;
open( $fh_table, ">", $ftable ) || die "cannot open $ftable: $!\n";

$fh_table->print("\$ralpha_table=[\n");

foreach my $symmetry ( 0, 1, 2 ) {
    $fh_table->print("    [\n");

    my $idel       = 1;
    my $igamma_max = $symmetry == 2 ? 7000 : 11000;
    for ( my $igamma = 1001 ; $igamma <= $igamma_max ; $igamma += $idel ) {
        if ( $igamma >= 1010 ) { $idel = 10 }
        if ( $igamma >= 2000 ) { $idel = 20 }
        if ( $igamma >= 3000 ) { $idel = 50 }
        if ( $igamma >= 4000 ) { $idel = 100 }
        if ( $igamma >= 7000 ) { $idel = 200 }
        my $gamma = $igamma / 1000;

        my $obj = Blast::IPS::SimilaritySolution->new(
            symmetry => $symmetry,
            gamma    => $gamma
        );

        my ( $alpha, $err, $it, $alpha_trap, $err0 ) =
          $obj->alpha_integral( $tol, $itmax );
        my $dalp = abs( $alpha - $alpha_trap );

        print "sym=$symmetry\tgamma=$gamma\talpha=$alpha\tdalp=$dalp\terr=$err\tit=$it\tI0=$err0\n";
        print STDERR "$symmetry\t$gamma\t$alpha\t$alpha_trap\t$dalp\t$err\t$it\t$err0\n";

        $err = sprintf( "%0.3g", $err );
        $dalp = sprintf( "%0.3g", $dalp );
        my $n =
            $err < 5.e-12 ? 12
          : $err < 5.e-11 ? 11
          : $err < 5.e-10 ? 10
          : $err < 5.e-9  ? 9
          :                 8;

        #$err< 5.e-8 ? 7 : 6
        my $format = "%0.$n" . "f";
        $alpha = sprintf( $format, $alpha );
        $fh_table->print("[$gamma, $alpha, $err, $dalp],\n");
    }
    $fh_table->print("   ],\n");
}
$fh_table->print("];\n");
$fh_table->close();
