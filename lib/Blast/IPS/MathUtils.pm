package Blast::IPS::MathUtils;

# Some math utilities for Blast::IPS

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
our @EXPORT_OK = qw(
  locate_2d
  macheps
  multiseg_integral
  nbrenti
  nbrentx
  parabolic_integral
  polint
  set_interpolation_points
  table_row_interpolation
  trapezoidal_integral
);
use Exporter;
our @ISA = qw(Exporter);

{
    my $macheps;

    BEGIN {
        $macheps = 1;

        # approxomate machine epsilon
        for ( my $i = 0 ; $i < 1000 ; $i++ ) {
            my $diff = 1 - $macheps;
            last if ( $diff eq 1 );
            $macheps /= 2;
        }
    }
    sub macheps { return $macheps }
}

sub multiseg_integral {

    # This is a fairly robust and accurate integrator.
    # It allows integrating functions which have some very smooth
    # sections as well as some discontinuous sections (like shock waves).  

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
    # using a symmetric version of Simpson's method

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

sub polint {

    #  Slightly modified versions of the "polint" routine from
    #  Press, William H., Brian P. Flannery, Saul A. Teukolsky and
    #  William T. Vetterling, 1986, "Numerical Recipes: The Art of
    #  Scientific Computing" (Fortran), Cambrigde University Press,
    #  pp. 80-82.

    # Given:
    # $xx = an x location where y is required
    # ($rx, $ry) = arrays of ($x,$y) lagrange interpolation points

    # Return:
    #  $yy = the interpolated value at $xx
    #  $dy = the estimated error

    # Example call:
    # my ( $yy, $dy ) = polint( $xx, $rx, $ry );

    my ( $xx, $rx, $ry ) = @_;

    my $n = @{$rx};

    # Return values
    my ( $yy, $dy );

    #..find the index ns of the closest table entry; initialize the c and d
    # tables
    my ( @c, @d );
    my $ns  = 0;
    my $dif = abs( $xx - $rx->[0] );
    for ( my $i = 0 ; $i < $n ; ++$i ) {
        my $dift = abs( $xx - $rx->[$i] );
        if ( $dift < $dif ) {
            $ns  = $i;
            $dif = $dift;
        }
        $c[$i] = $ry->[$i];
        $d[$i] = $ry->[$i];
    }

    #..first guess for y
    $yy = $ry->[$ns];

    #..for each column of the table, loop over the c's and d's and update them
    $ns = $ns - 1;
    for ( my $m = 0 ; $m < $n - 1 ; ++$m ) {
        for ( my $i = 0 ; $i < $n - $m - 1 ; ++$i ) {
            my $ho  = $rx->[$i] - $xx;
            my $hp  = $rx->[ $i + $m + 1 ] - $xx;
            my $w   = $c[ $i + 1 ] - $d[$i];
            my $den = $ho - $hp;
            if ( $den == 0.0 ) {
                print STDERR "polint: 2 rx entries are the same \n";
                return;
            }
            $den   = $w / $den;
            $d[$i] = $hp * $den;
            $c[$i] = $ho * $den;
        }

        # after each column is completed, decide which correction c or d, to
        # add to the accumulating value of y, that is, which path to take in
        # the table by forking up or down. ns is updated as we go to keep track
        # of where we are. the last dy added is the error indicator.
        if ( 2 * ( $ns + 1 ) < $n - $m - 1 ) {
            $dy = $c[ $ns + 1 ];
        }
        else {
            $dy = $d[$ns];
            $ns = $ns - 1;
        }
        $yy += $dy;
    }
    return wantarray ? ( $yy, $dy ) : $yy;
}

sub set_interpolation_points {

    my ( $jfloor, $ntab, $NLAG ) = @_;

    # Find the index range for NLAG valid lagrange interpolation points
    # Given:
    #   $jfloor is the first index before the point of interest
    #   $ntab is the number of points in the table
    #   $NLAG is the number of interpolation points desired
    # Return:
    #   a reference to a list of consecutive indexes
    #   the actual number may be less than NLAG for small tables

    return if ( $ntab <= 0 || $NLAG <= 0 );

    # First add points on both sides (will be lopsided if NLAG is odd)
    #my $j_lo = $jfloor - int( $NLAG / 2 ); # ORIGINAL, lopsided
    my $j_lo = $jfloor - int( ( $NLAG - 1 ) / 2 );    # Corrected
    my $j_hi = $j_lo + $NLAG - 1;

    #print STDERR "jfloor=$jfloor, j_lo=$j_lo, j_hi=$j_hi\n";

    # Shift points if too near an edge
    if ( $j_lo < 0 ) {
        my $nshift = 0 - $j_lo;
        $j_lo += $nshift;
        $j_hi += $nshift;
    }
    if ( $j_hi > $ntab - 1 ) {
        my $nshift = $ntab - 1 - $j_hi;
        $j_lo += $nshift;
        $j_hi += $nshift;
    }

    # Be sure points are in bounds
    if ( $j_lo < 0 )         { $j_lo = 0 }
    if ( $j_hi > $ntab - 1 ) { $j_hi = $ntab - 1 }

    #print STDERR "returning j_lo=$j_lo, j_hi=$j_hi\n";
    return [ $j_lo .. $j_hi ];
}

sub locate_2d {
    my ( $rx, $xx, $icol, $jl, $ju ) = @_;

    # Binary search for two consecutive table row indexes, jl and ju, of a
    # 2D matrix such that the value of x lies between these two table values.  
    # $icol is the column of the variable x
    # If x is out of bounds, returns either jl<0 or ju>=N

    # Feed back old $jl and $ju for improved efficiency

    $icol=0 unless defined($icol);
    my $num=@{$rx};
    if ( $icol > @{ $rx->[0] } - 1 ) {
        print STDERR
          "ERROR in MathUtils::locate_2d: column=$icol exceeds table size\n";
	return;
    }
    if (!defined($xx)) {
        print STDERR
          "ERROR in MathUtils::locate_2d: x is undefined\n";
	return;
    }
    my $dx_is_positive = $rx->[ -1 ]->[$icol] > $rx->[0]->[$icol];

    # ... but reset to beyond end of table if no longer valid
    $jl = -1
      if ( !defined($jl)
        || $jl < 0
        || $jl >= $num
        || ( $xx > $rx->[$jl]->[$icol] ne $dx_is_positive ) );
    $ju = $num
      if ( !defined($ju)
        || $ju < 0
        || $ju >= $num
        || ( $xx > $rx->[$ju]->[$icol] eq $dx_is_positive ) );
	
    # Loop until the requested point lies in a single interval
    while ( $ju - $jl > 1 ) {
        my $jm = int( ( $jl + $ju ) / 2 );
        if ( $xx > $rx->[$jm]->[$icol] eq $dx_is_positive ) {
            $jl = $jm;
        }
        else {
            $ju = $jm;
        }
    }
    return ( $jl, $ju );
}

sub OLD_locate_2d {
    my ( $xx, $icol, $rtab, $jl, $ju ) = @_;

    # Binary search for two consecutive table row indexes, jl and ju, of a
    # 2D matrix such that the value of x lies between these two table values.
    # $icol is the column of the variable x
    # If x is out of bounds, returns either jl<0 or ju>=N

    #  Based on the fortran code in
    #  Press, William H., Brian P. Flannery, Saul A. Teukolsky and
    #  William T. Vetterling, 1986, "Numerical Recipes: The Art of
    #  Scientific Computing" (Fortran), Cambrigde University Press.

    return unless $rtab;
    my $num = @{$rtab};
    return unless $num;
    my $dx_is_positive = $rtab->[-1]->[$icol] > $rtab->[0]->[$icol];

    # but reset to beyond end of table if no longer valid
    $jl = -1
      if ( !defined($jl)
        || $jl < 0
        || $jl >= $num
        || ( $xx > $rtab->[$jl]->[$icol] ne $dx_is_positive ) );
    $ju = $num
      if ( !defined($ju)
        || $ju < 0
        || $ju >= $num
        || ( $xx > $rtab->[$ju]->[$icol] eq $dx_is_positive ) );

    # Loop until the requested point lies in a single interval
    while ( $ju - $jl > 1 ) {
        my $jm = int( ( $jl + $ju ) / 2 );
        if ( $xx > $rtab->[$jm]->[$icol] eq $dx_is_positive ) {
            $jl = $jm;
        }
        else {
            $ju = $jm;
        }
    }
    return ( $jl, $ju );
}

sub table_row_interpolation {

    # interpolate all row values of a table to an arbitrary x value All values
    # in a row must be floating point, and the same NLAG is used for all

    # caller can look at $jl and $ju to see if the results are off the table

    my ( $xx, $rtable, $icx, $NLAG, $jl, $ju, $no_extrap_l, $no_extrap_u, ) =
      @_;

    # $xx = the value of the X variable
    # $rt   = table
    # $icx = column number of X variable [default is 0]
    # $NLAG = #Lagrange points [default 4]
    # $jl, $ju = optional last search location, for speed
    # $no_extrap_l = >0 to prevent extrapolation before start of OLD table
    # 		[default=0]
    # $no_extrap_u = >0 to prevent extrapolation beyond end of OLD table
    # 		[default=0]

    $icx  = 0 unless defined($icx);
    $NLAG = 4 unless defined($NLAG);

    # number of lagrange points cannot exceed number of table points
    my $npts = @{$rtable};
    if ( $NLAG > $npts ) { $NLAG = $npts }
    my $NH = int( $NLAG / 2 );

    my $jpoly_l = -1;

    ##( $jl, $ju ) = locate_2d( $xx, $icx, $rtable, $jl, $ju );
    ( $jl, $ju ) = locate_2d( $rtable, $xx, $icx, $jl, $ju );
    my $rrow;

    # loop just to simplify break
    while ( !$rrow ) {
        if ( $jl >= 0 && $ju < $npts ) {

            # We want to have about equal numbers around this point
            $jpoly_l = $jl - $NH + 1;

            # But keep the points within the old table
            my $jpoly_u = $jpoly_l + $NLAG - 1;
            my $dj = $jpoly_u - ( $npts - 1 );
            if ( $dj > 0 ) {
                $jpoly_l -= $dj;
            }
            $jpoly_l = 0 if ( $jpoly_l < 0 );
        }
        elsif ( $jl < 0 ) {
            if ($no_extrap_l) { $rrow = [ @{ $rtable->[0] } ]; last }
            $jpoly_l = 0;
        }
        else {
            #if ($no_extrap_u) { $rrow=$rtable->[-1]; last }
            if ($no_extrap_u) { $rrow = [ @{ $rtable->[-1] } ]; last }
            $jpoly_l = $npts - $NLAG;
            $jpoly_l = 0 if ( $jpoly_l < 0 );
        }

        my $nvars = @{ $rtable->[0] };

        my $rlag_x = [];
        my $jj     = $jpoly_l;
        for ( my $n = 0 ; $n < $NLAG ; $n++ ) {
            last if ( $jj > $npts - 1 );    ## For safety
            push @{$rlag_x}, $rtable->[$jj]->[$icx];
            $jj++;
        }

        # now loop to define all values
        # Note that we are also interpolating the x variable as a control
        for ( my $icy = 0 ; $icy < $nvars ; $icy++ ) {
            my $rlag_y = [];
            my $jj     = $jpoly_l;
            for ( my $n = 0 ; $n < $NLAG ; $n++ ) {
                last if ( $jj > $npts - 1 );    ## For safety
                push @{$rlag_y}, $rtable->[$jj]->[$icy];
                $jj++;
            }
            my $yy = polint( $xx, $rlag_x, $rlag_y );
            $rrow->[$icy] = $yy;
        }
        last;
    }
    return wantarray ? ( $rrow, $jl, $ju ) : $rrow;
}

sub nbrenti {

    # call arg declarations
    my ( $sa, $fa, $sb, $fb, $tol, $bvec ) = @_;

    #       initialization for Brent's great root finder
    #       (see-Algorithms for Minimization Without
    #       Derivatives, by Richard P. Brent, Prentice-Hall, 1973.)

    #       This is a direct translation of Brent's fortran code
    #       published in the above reference; translated with help of f2pl
     
    #       input parameters -
    #       bounding states:
    #         sa - first x coordinate
    #         fa = f(sa)
    #         sb - second x coordinate
    #         fb = f(sb)
    #         **  must have fa*fb <= 0. **
    #       tolx = convergence tolerance on x
    #
    #       output parameter -
    #       xnext = next value of x to try
    #
    #       there are two subroutines:
    #               brenti must be called once to start an interation
    #               brentx must be called once for each step
    #
    #****************************************************************
    #
    #       sample coding: assume we have a root trapped between sa & sb
    #
    #       dimension bvec(7)
    #       data tol/.001/,maxit/25/
    #
    #     call nbrenti (sa,fa,sb,fb,tol,xx,bvec)
    #
    #     loop over iterations
    #     do 250 iter=1,maxit
    #
    #       calculate ff(xx)
    #
    #       get next value of xx
    #     call nbrentx (ff,xx,ifconv,bvec)
    #     if (ifconv.ne.0) go to 700
    #
    #       note - it is possible to make a convergence test here on the
    #       magnitude of ff as well
    #
    # 250 continue
    #
    #     no convergence after maxit iterations -- cant happen if maxit
    #               is large enough for the given tol - see brents book
    #
    #     converged - finish up
    # 700 ...
    #****************************************************************
    #
    my $sc        = $sb;
    my $fc        = $fb;
    $bvec->[1] = $sa;
    $bvec->[2] = $fa;
    $bvec->[3] = $sb;
    $bvec->[4] = $fb;
    $bvec->[5] = $sc;
    $bvec->[6] = $fc;
    $bvec->[7] = $tol;
    my ($xnext, $ifconv) = nbrentx( $fb, $bvec );
    return ($xnext, $ifconv);
}

sub nbrentx {

    # call arg declarations
    my ( $flast, $bvec ) = @_;

    #
    #     take one iteration step with brents method
    #     flast = value of f of previous step
    #     xnext = next value of x to try
    #     ifconv = 0 if not converged
    #            = 1 if converged
    #       bvec = state vector
    #
    #
    #       get the previous state
    my $sa  = $bvec->[1];
    my $fa  = $bvec->[2];
    my $sb  = $bvec->[3];
    my $fb  = $bvec->[4];
    my $sc  = $bvec->[5];
    my $fc  = $bvec->[6];
    my $tol = $bvec->[7];

    my $slast = $sb; ## SLH

    my ($bd, $be, $bs, $bm, $bp, $bq, $br, $xnext, $ifconv);
    $be=0; # patch because abs function below 230 uses be when not defined

    $fb = $flast;
    if ( ( $fb > 0. ) && ( $fc > 0. ) ) {
        goto L210;
    }
    if ( ( $fb <= 0. ) && ( $fc <= 0. ) ) {
        goto L210;
    }
    goto L220;

    #
  L210:
    $sc = $sa;
    $fc = $fa;
    $be = $sb - $sa;
    $bd = $be;
  L220:
    if ( abs($fc) >= abs($fb) ) {
        goto L230;
    }
    $sa = $sb;
    $sb = $sc;
    $sc = $sa;
    $fa = $fb;
    $fb = $fc;
    $fc = $fa;
  L230:
    $bm = 0.5 * ( $sc - $sb );

    if ( ( abs($bm) <= $tol ) || ( $fb == 0. ) ) {
        goto L340;
    }
    if ( ( abs($be) >= $tol ) && ( abs($fa) > abs($fb) ) ) {
        goto L240;
    }

    #
    #         bisection is forced
    $be = $bm;
    $bd = $be;
    goto L300;

    #
  L240:
    $bs = $fb / $fa;
    if ( $sa != $sc ) {
        goto L250;
    }

    #
    #         linear interpolation
    $bp = 2. * $bm * $bs;
    $bq = 1. - $bs;
    goto L260;

    #
    #         inverse quadratic interpolation
  L250:
    $bq = $fa / $fc;
    $br = $fb / $fc;
    $bp = $bs *
      ( 2. * $bm * $bq * ( $bq - $br ) - ( $sb - $sa ) * ( $br - 1. ) );
    $bq = ( $bq - 1. ) * ( $br - 1. ) * ( $bs - 1. );

    #
  L260:
    if ( $bp <= 0. ) {
        goto L270;
    }
    $bq = -$bq;
    goto L280;
  L270:
    $bp = -$bp;
  L280:
    $bs = $be;
    $be = $bd;

    if (   ( 2. * $bp >= 3. * $bm * $bq - abs( $tol * $bq ) )
        || ( $bp >= abs( 0.5 * $bs * $bq ) ) )
    {
        goto L290;
    }
    $bd = $bp / $bq;
    goto L300;
  L290:
    $be = $bm;
    $bd = $be;
  L300:
    $sa = $sb;
    $fa = $fb;

    if ( abs($bd) <= $tol ) {
        goto L310;
    }
    $sb = $sb + $bd;
    goto L330;
  L310:
    if ( $bm <= 0. ) {
        goto L320;
    }
    $sb = $sb + $tol;
    goto L330;
  L320:
    $sb = $sb - $tol;
  L330:
    $xnext  = $sb;
    $ifconv = 0;
    goto L900;
  L340:
    $ifconv = 1;
    $xnext  = $slast;  # SLH: added because not defined otherwise

    #
    #       save the state
  L900:
    $bvec->[1] = $sa;
    $bvec->[2] = $fa;
    $bvec->[3] = $sb;
    $bvec->[4] = $fb;
    $bvec->[5] = $sc;
    $bvec->[6] = $fc;
    return ($xnext, $ifconv);
}


1;

