#!/usr/bin/perl
use strict;
use warnings;

# This program creates a table of alpha values for the point source similarity
# solution.  Output values are in the table Blast::IPS::AlphaTable.pm.  It is
# saved for reference but should not have to be rerun.

# References consulted:

# 1. Blast Wave, H.Bethe, K.Fuchs, J.O.Hirschfelder,J. Magee,
# R.Peirels,J.vonNeumann, Aug 1947, Los Alamos LA-2000.  Has Von Neumann's
# solution in spherical symmetry.

# 2. Sadek, H. S., Gottlieb, J. J., "Initial Decay Of Flow Properties Of
# Planar, Cyl1Ndr1Cal And Spher1Cal Blast Waves", 1983.  Institute for
# Aerospace Studies, University of Toronto. Has the Von Neumann solution
# extended to cylindrical and plane symmetry - but note a typo.

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

# Settings..

# itmax = max iterations; for testing use itmax=3 then increase
my $itmax = 8;  # Use small value (such as 2 or 3) for testing

# tol = stopping tolerance on absolute value of alpha
my $tol   = 1.e-9;

# output file
my $ftable = "alpha_table_$itmax.txt";

my $fh_table;
open( $fh_table, ">", $ftable ) || die "cannot open $ftable: $!\n";

$fh_table->print("\$ralpha_table=[\n");

foreach my $symmetry ( 0, 1, 2 ) {
    $fh_table->print("    [\n");
    my $igamma=110;
    my $idel=1;
    my $igamma_max=700;
    for (my $igamma=110; $igamma<=700; $igamma+=$idel) {
	if ($igamma>=200) {$idel=2}
	if ($igamma>=300) {$idel=5}
	if ($igamma>=400) {$idel=10}
	my $gamma=$igamma/100;
        my ($alpha,$err,$it) = get_alpha($symmetry,$gamma,$tol,$itmax);

        print "sym=$symmetry\tgamma=$gamma\talpha=$alpha\terr=$err\tit=$it\n";
        print STDERR "$symmetry\t$gamma\t$alpha\t$err\t$it\n";

        $err = sprintf( "%0.3g", $err );
        my $n =
            $err < 5.e-11 ? 11
          : $err < 5.e-10 ? 10
          : $err < 5.e-9  ? 9
          :                 8;

        #$err< 5.e-8 ? 7 : 6
        my $format = "%0.$n" . "f";
        $alpha = sprintf( $format, $alpha );
        $fh_table->print("[$gamma, $alpha, $err],\n");
    }
    $fh_table->print("   ],\n");
}
$fh_table->print("];\n");
$fh_table->close();

sub get_alpha {
    my ( $N, $gamma, $tol, $itmax ) = @_;
    my $num = 10000;
    my ( $alpha, $err );
    my $it;

    # Patch for gamma=2 and gamma=7 spherical symmetry..
    # We will perturb gamma a little to avoid divdes by zero
    if ( $N == 2 && $gamma == 7 ) { $gamma -= $tol / 100 }
    if ( $gamma == 2 ) { $gamma -= $tol / 100 }

    # brute-force approach:
    # Iterate until we reach the desired tolerance or the max iteration count
    for ( $it = 0 ; $it <= $itmax ; $it++ ) {
        my $alpha_last = $alpha;
        $num *= 2;
        my $dx = 1 / $num;
        my $coef = coef( $gamma, $N );
        my $rxf;
        my $rxr;

        for ( my $i = 0 ; $i <= $num ; $i++ ) {
            my $x    = $i * $dx;
            my $f    = fofx( $x, $gamma, $N );
            my $rrat = rrat( $x, $gamma, $N );
            push @{$rxf}, [ $rrat**( $N + 1 ), $f ];
        }
        my $rxy = multiseg_integral($rxf);
        $alpha = $coef * $rxy->[-1]->[1];
        $err = ( $it > 0 ) ? abs( $alpha - $alpha_last ) : 10*$tol;
	last if ($err<$tol);
    }

    # Correct for plane symmetry: I am using the half geometry for the problem in plane symmetry
    if ( $N == 0 ) { $alpha /= 2; $err /= 2 }

    return wantarray ? ( $alpha, $err, $it ) : $alpha;
}

sub coef {

    # Leading coefficient
    my ( $gamma, $N ) = @_;
    my $pi = 4 * atan2( 1, 1 );
    my $c =
      ( 16 * ( ( 1 - $N ) * ( 2 - $N ) + 2 * $N * $pi ) ) /
      ( ( $N + 1 ) * ( $N + 3 )**2 * ( $gamma + 1 ) * ( $gamma - 1 ) );
    return $c;
}

sub rrat {
    my ( $x, $gamma, $N ) = @_;

    # Ratio of radius to shock radius
    my $t1 = $x;
    my $p1 = ( $gamma - 1 ) / ( 2 * $gamma + $N - 1 );
    my $t2 = ( $x + 1 ) / 2;
    my $p2 = (-2) / ( $N + 3 );
    my $t3 = ( $x + $gamma ) / ( 1 + $gamma );
    my $p3 = ( $gamma + 1 ) / ( ( $N + 1 ) * $gamma - $N + 1 );

    # NOTE CORRECTION based on la2000; the -N+1 should have been +N-1:
    my $t4_top = ( ( $N + 1 ) * ( 2 - $gamma ) * $x + 2 * $gamma + $N - 1 );
    my $t4_bot = ( ( 3 * $N + 1 ) - ( $N - 1 ) * $gamma );

    # avoid accidental divide by zero at symmetry=2 gamma=7
    if ($t4_bot == 0) {
	return 0;
    }
    my $t4 = $t4_top / $t4_bot;

    # avoid accidental divide by zero at gamma=2
    my $p4 = 1;
    if ( $t4 != 1 ) {
        $p4 =
          -( ( $N**2 + 2 * $N + 5 ) * $gamma**2 -
              ( 3 * $N**2 - 2 * $N - 1 ) * $gamma +
              4 * ( $N**2 - 1 ) ) /
          ( ( $N + 3 ) *
              ( 2 * $gamma + $N - 1 ) *
              ( ( $N + 1 ) * $gamma - $N + 1 ) );
    }

    my $rrat = $t1**$p1 * $t2**$p2 * $t3**$p3 * $t4**$p4;
    return $rrat;
}

sub fofx {
    my ( $x, $gamma, $N ) = @_;

    # integrate this from x=0 to x=1
    my $p1 = ( 3 * $N + 5 ) / ( $N + 3 );
    my $t1 = ( $x + 1 ) / 2;
    my $p2 = ( -2 * $N * $gamma ) / ( ( $N + 1 ) * $gamma - $N + 1 );
    my $t2 = ( $x + $gamma ) / ( 1 + $gamma );

    # NOTE CORRECTION based on la2000; the -N+1 should have been +N-1:
    my $t3_top = ( ( $N + 1 ) * ( 2 - $gamma ) * $x + 2 * $gamma + $N - 1 );
    my $t3_bot = ( ( 3 * $N + 1 ) - ( $N - 1 ) * $gamma );

    # avoid divide by zero at symmetry=2 gamma=7
    if ($t3_bot == 0) {
	return 0;
    }
    my $t3 = $t3_top / $t3_bot;

    # avoid divide by zero at gamma=2
    my $p3 = 1;
    if ( $t3 != 1 ) {
        $p3 =
          ( ( $N**2 + 2 * $N + 5 ) * $gamma**2 -
              ( 3 * $N**2 - 2 * $N - 1 ) * $gamma +
              4 * ( $N**2 - 1 ) ) /
          ( ( $N + 3 ) * ( 2 - $gamma ) * ( ( $N + 1 ) * $gamma - $N + 1 ) );
    }

    my $ff = $t1**$p1 * $t2**$p2 * $t3**$p3;
    return $ff;
}

sub multiseg_integral {

    # Given a table of irregularly spaced (x, dydx) values
    #   $rxyp = [ x, dydx ]
    # and a starting value
    #   $yy_i = a starting value of y  [default=0]
    #
    # and OPTIONALLY (see below)
    #   $rpflags = a list of same length with flags with values:
    #      0 = not smooth, use trapezoidal method
    #      1 = smooth, parabolic method okay
    #      2 = forced break point if parabolic
    # Return the integral
    #  $rxy = [x, y]
    # Using the parabolic method for regions with 1's and trapezoidal
    # otherwise;  
    # At a segment with 0 at either end, method will be 0 (trapezoidal)
    # At a segment with 2 at both ends (forced break), method will be 0

    # OPTIONS:
    # * If rpflags is not defined, flags are automatically set to have
    # break points where mesh changes significantly (2%)
    # * If rpflags is a scalar value, it is the percent change used to set
    # breakpoints (should be a small fraction say from 0 to .1)

    # This routine allows integrating functions which have some very smooth
    # sections as well as some discontinuous sections.  

    my ( $rxyp, $yy_i, $unknown ) = @_;
    return unless ($rxyp && @{$rxyp});

    my $mark_mesh_breaks = sub {
        my ($eps) = @_;

	# Routine for optional automatic setting mesh break points
        # Set flags such that
        # flag=1 where the grid spacing is almost uniform
        # flag=2 where the grid spacing suddenly changes
        $eps = 0.02 unless defined($eps);
        if ( $eps < 0 ) { $eps = 0 }

        my $num = @{$rxyp};
        my $dxl = 0;
        my $dxr = 0;
        my $rpflags;
        for ( my $i = 0 ; $i < $num ; $i++ ) {
            $dxl = $dxr;
            if ( $i < $num - 1 ) {
                $dxr = $rxyp->[ $i + 1 ]->[0] - $rxyp->[$i]->[0];
            }
            else { $dxr = 0 }
            if ( ( $dxl > 0 && $dxr > 0 || $dxl < 0 && $dxr < 0 )
                && abs( $dxl - $dxr ) < $eps * abs( $dxl + $dxr ) )
            {
                $rpflags->[$i] = 1;
            }
            else {
                $rpflags->[$i] = 2;
            }
        }
        return $rpflags;
    };

    my $rpflags = $unknown;

    # Handle input options for the flags
    my $ref_type=ref($unknown);
    if ( !$ref_type ) {
        $rpflags=$mark_mesh_breaks->($unknown); 
    }
    else {
        $rpflags = $unknown;
    }

    # Verify same array lengths
    my $num  = @{$rxyp};
    my $nump = @{$rpflags};
    if ( $nump != $num ) {
        die "Input error to 'multiseg_integral', num=$num but nump=$nump\n";
    }

    my $segment_flag = sub {

        my ($i) = @_;

        # Returns the segment flag for a segment from i to i+1
        # SO: MUST NOT BE CALLED FOR THE LAST POINT
        # Segment flag is 0 or 1:
        #  = 0 use trapezoidal method for this segment
        #  = 1 use parabolic method for this segment if possible

	# Table of input and output values:
        # Inputs: we are given a Left and Right Point flags L and R
        # Output: We want these Segment flags S
        #    L R S
        #    0 0 0
        #    0 1 0
        #    0 2 0
        #    1 0 0
        #    1 1 1
        #    1 2 1
        #    2 0 0
        #    2 1 1
        #    2 2 0

        # Cases 0-*, *-0, 2-2
        if (   $rpflags->[$i] == 0
            || $rpflags->[ $i + 1 ] == 0
            || $rpflags->[ $i + 1 ] == 2 && $rpflags->[$i] == 2 )
        {
            return 0;
        }

        # Cases 1-1, 1-2, 2-1
        else {
            return 1;
        }
    };

    my $rxy;


    # Now do the integral region by region, from i=$ib to i=$ie
    my $ie = 0;
    while ( $ie < $num - 1 ) {

        my $ib = $ie;

        # sflag indicates the type of integral for the segments we are
        # accumulating between i=$ib and i=$ie
        my $sflag = $segment_flag->($ib);
        for ( $ie = $ib + 1 ; $ie < $num ; $ie++ ) {

            # Look at next segment flag to see if we can continue
            if ( $ie < $num - 1 ) {
                my $sflag_next = $segment_flag->($ie);

                my $want_break = $rpflags->[$ie] == 2;
		if ($want_break && $ie-$ib==1) {$sflag=0}

                # Continue building the chain unless type changes or we hit
                # a forced break between parabolic secitons 
		# (a forced break in trapezoidal region means nothing)
                if ( $sflag_next == $sflag
                    && ( $sflag == 0 || $rpflags->[$ie] != 2 ) )
                {
                    next;
                }

		# If we have just 2 parabolic points, reset to trapezoid 
		# and keep going
                if ( $sflag == 1 && $ie - $ib == 1 ) { $sflag = 0; next }

            }
	    last;
        }

        # Now use method for $sflag between $ib and $ie
        my $rxyp_part;
        my $rxy_part;
        for ( my $i = $ib ; $i <= $ie ; $i++ ) {
            $rxyp_part->[ $i - $ib ] = $rxyp->[$i];
        }
        if ($sflag) {
            $rxy_part = parabolic_integral( $rxyp_part, $yy_i );
        }
        else {
            $rxy_part = trapezoidal_integral( $rxyp_part, $yy_i );
        }
        $yy_i = $rxy_part->[-1]->[1];
        for ( my $i = $ib ; $i <= $ie ; $i++ ) {
            $rxy->[$i] = $rxy_part->[ $i - $ib ];
        }
    }
    return ($rxy);
}
sub parabolic_integral {

    # symmetric version of parabolic integration for irregularly spaced points
    # given a table of irregularly spaced (x, dydx) values
    #   $rxyp = [ x, dydx ]
    # and
    #   $yy_i = a starting value of y  [default=0]
    # return the integral
    #  $rxy = [x, y]
    # using a symmetric version of simpson's method

    my ( $rxyp, $yy_i ) = @_;
    $yy_i = 0 unless defined($yy_i);

    # we are integrating dydx to get y(x)

    # loop over every two points
    my $rxy;
    my $ystore = sub {
            my ( $ix, $yy_i ) = @_;
            $rxy->[$ix] = [ $rxyp->[$ix]->[0], $yy_i ];
    };
   
    my $num = defined($rxyp) ? @{$rxyp} : 0;
    if ( $num <= 0 ) {
    }
    elsif ( $num == 1 ) {
	$ystore->(0, $yy_i);
    }
    elsif ( $num == 2 ) {
        my $ii_l = 0;
        my $ii_u = 1;
        my ( $xx_l, $dydx_l ) = @{ $rxyp->[$ii_l] };
        my ( $xx_u, $dydx_u ) = @{ $rxyp->[$ii_u] };
        my $dy = 0.5 * ( $dydx_l + $dydx_u ) * ( $xx_u - $xx_l );
        $ystore->( 0, $yy_i );
        $ystore->( 1, $yy_i + $dy );
    }
    else {

        # first pass loop
        # we loop over all interior points ii_c and compute the two parts of
        # parabolic integral, and store them.
        # notation: we are working with 3 points, l,c, u:
        # _l, _c, _u = lower, center, upper
        my $rtmp;
        for ( my $ii_c = 1 ; $ii_c < $num - 1 ; $ii_c += 1 ) {
            my $ii_l = $ii_c - 1;
            my $ii_u = $ii_c + 1;
            my ( $xx_l, $dydx_l ) = @{ $rxyp->[$ii_l] };
            my ( $xx_c, $dydx_c ) = @{ $rxyp->[$ii_c] };
            my ( $xx_u, $dydx_u ) = @{ $rxyp->[$ii_u] };

            my $dxlc = $xx_l - $xx_c;
            my $dxuc = $xx_u - $xx_c;
            my $dflc = $dydx_l - $dydx_c;
            my $dfuc = $dydx_u - $dydx_c;
            my $det  = $dxlc * $dxuc * ( $dxuc - $dxlc );

            # initialize to trapezoidal values as a backup method
            my $di_1 = -0.5 * ( $dydx_l + $dydx_c ) * $dxlc;
            my $di_2 = 0.5 *  ( $dydx_c + $dydx_u ) * $dxuc;

            # use parabolic shape if safe ...
            # avoid if sign of dx changes or if the ratio of the two
            # dx values is extreme (more than 5:1 here)
            my $trapez_flag =
                 ( $det == 0 )
              || $dxlc * $dxuc >= 0
              || -$dxlc / $dxuc < 0.2
              || -$dxlc / $dxuc > 5;

            if ( !$trapez_flag ) {
                my $aa = $dydx_c;
                my $bb =
                  ( $dxuc * $dxuc * $dflc - $dxlc * $dxlc * $dfuc ) / $det;
                my $cc = ( -$dxuc * $dflc + $dxlc * $dfuc ) / $det;
                $di_1 =
                  0 - $dxlc * ( $aa + $bb / 2 * $dxlc + $cc / 3 * $dxlc**2 );
                $di_2 = $dxuc * ( $aa + $bb / 2 * $dxuc + $cc / 3 * $dxuc**2 );
            }
 	    else {
		#print stderr "booga\n";
	    }
            $rtmp->[$ii_c] = [ $di_1, $di_2 ];
        }

        # second pass
        # we sum the pieces of the integral
        $ystore->( 0, $yy_i );

        # first segment stores the left half of the first parabola
        $yy_i += $rtmp->[1]->[0];
        $ystore->( 1, $yy_i );

        # interior points sum the average of the two overlapping halves
        for ( my $ii = 2 ; $ii < $num - 1 ; $ii++ ) {
            my $dy = 0.5 * ( $rtmp->[ $ii - 1 ]->[1] + $rtmp->[$ii]->[0] );
            $yy_i += $dy;
            $ystore->( $ii, $yy_i );
        }

        # last segment stores the right half of the last parabola
        $yy_i += $rtmp->[$num-2]->[1];
        $ystore->( $num - 1, $yy_i );
    }
    return ($rxy);
}

sub trapezoidal_integral {

    my ( $rxyp, $yy_i ) = @_;
    $yy_i = 0 unless defined($yy_i);

    # we are integrating dydx to get y(x)
    # notation: we are working with 3 points, l,c, u:
    # _l, _c, _u = lower, center, upper

    # main loop
    my $num = @{$rxyp};
    my $rxy;
    $rxy->[0]=[$rxyp->[0]->[0], $yy_i];
    for ( my $ii_u = 1 ; $ii_u < $num ; $ii_u += 1 ) {
        my $ii_l = $ii_u - 1;
        my ( $xx_l, $dydx_l ) = @{ $rxyp->[$ii_l] };
        my ( $xx_u, $dydx_u ) = @{ $rxyp->[$ii_u] };
        $yy_i += 0.5 * ( $dydx_l + $dydx_u ) * ($xx_u-$xx_l);
        $rxy->[$ii_u] = [ $xx_u, $yy_i ];
    }
    return ($rxy);
}
