package Blast::IPS::SimilaritySolution;

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

# This is a collection of software related to the very high pressure region of
# an ideal point source explosion where it may be described with the similarity
# solution.

use strict;
use warnings;
use 5.006;

use Blast::IPS::MathUtils qw(
  _locate_2d
  multiseg_integral
  nbrenti
  nbrentx
  polint
  set_interpolation_points
);

use Blast::IPS::AlphaTable qw(alpha_interpolate);

my $VERSION = 1.00;
use Carp;

sub new {
    my $class = shift @_;
    my $self  = {};
    bless $self, $class;
    $self->_setup(@_);
    return $self;
}

sub _setup {
    my ( $self, @args)=@_;

    my @input_keys = qw(
      gamma
      symmetry
      E0
      rho_amb
    );
    my %valid_input_keys;
    @valid_input_keys{@input_keys} = (1) x scalar(@input_keys);

    my $rinput_hash;
    my $reftype;
    if (@args) {
        my $arg0    = $args[0];
        my $reftype = ref($arg0);
        if ( !$reftype ) {

            if ( defined($arg0) ) {

                # simple hash of named values
                my %input_hash = @args;
                $rinput_hash = \%input_hash;
            }

        }
        elsif ( $reftype eq 'HASH' ) {
            $rinput_hash = $arg0;
        }
        else {
            carp "Unexpected ref type: $reftype\n";
            return;
        }
    }
    else {

        # default to spherical table with gamma 1.4 if no arg given
        $rinput_hash = { symmetry => 2, gamma => 1.4, E0=>1, rho_amb=>1.4 };
    }

    # Validate input keys
    _check_keys( $rinput_hash, \%valid_input_keys,
        "Checking for valid input keys" );

    my $gamma    = $rinput_hash->{gamma};
    my $symmetry = $rinput_hash->{symmetry};
    my $E0       = $rinput_hash->{E0};
    my $rho_amb     = $rinput_hash->{rho_amb};

    my $error = "";
    if ( !defined($E0) )       { $E0       = 1 }
    if ( !defined($gamma) )    { $gamma    = 1.4 }
    if ( !defined($rho_amb) )     { $rho_amb     = 1 }
    if ( !defined($symmetry) ) { $symmetry = 2 }
    if ( !defined($rho_amb) || $rho_amb <= 0 ) {
        $error .= "Bad rho_amb='$rho_amb'\n";
    }
    if ( !defined($gamma) || $gamma <= 1.0 ) {
        $error .= "Bad gamma='$gamma'\n";
    }
    if ( !defined($symmetry) || $symmetry !~ /^[012]$/ ) {
        $error .= "Bad symmetry='$symmetry'\n";
    }
    my $alpha = alpha_interpolate( $symmetry, $gamma );

    $self->{_E0}       = $E0;
    $self->{_rho_amb}     = $rho_amb;
    $self->{_gamma}    = $gamma;
    $self->{_symmetry} = $symmetry;
    $self->{_Eoad}     = $E0 / ( $rho_amb * $alpha );
    $self->{_alpha}    = $alpha;
    if ($error) { carp "$error" }
    return; 
}

sub _check_keys {
    my ( $rtest, $rvalid, $msg, $exact_match ) = @_;

    # Check the keys of a hash for validity:
    # $rtest  = ref to hash to test
    # $rvalid = ref to has with valid keys

    # $msg = a message to write in case of error
    # $exact_match defines the type of check:
    #     = false: test hash must not have unknown key
    #     = true:  test hash must have exactly same keys as known hash
    my @unknown_keys =
      grep { $_ && !exists $rvalid->{$_} } keys %{$rtest};
    my @missing_keys =
      grep { !exists $rtest->{$_} } keys %{$rvalid};
    my $error = @unknown_keys;
    if ($exact_match) { $error ||= @missing_keys }
    if ($error) {
        local $" = ')(';
        my @expected_keys = sort keys %{$rvalid};
        @missing_keys = sort @missing_keys;
        @unknown_keys = sort @unknown_keys;
        croak(<<EOM);
------------------------------------------------------------------------
Blast::IPS::SimilaritySolution program error detected checking hash keys
Message is: '$msg'
Valid keys are: (@expected_keys)
Keys not seen : (@missing_keys)
Unknown key(s): (@unknown_keys)
------------------------------------------------------------------------
EOM
    }
    return;
}

sub get_Eoad { $_[0]->{_Eoad} }
sub get_alpha { $_[0]->{_alpha} }
sub get_error { $_[0]->{_error} }
sub get_gamma { $_[0]->{_gamma} }
sub get_rho_amb { $_[0]->{_rho_amb} }
sub get_E0 { $_[0]->{_E0} }

sub shock_front_values_as_hash {
    my ( $self, @args ) = @_;

    # given any one of these parameters:
    #  D - shock speed
    #  R - shock radius
    #  T - time
    #  P - shock pressure
    #  U - shock particle velocity

    # Returns a hash of shock front values 

    my $rinput_hash;
    my $reftype;
    if (@args) {
        my $arg0    = $args[0];
        my $reftype = ref($arg0);
        if ( !$reftype ) {

            if ( defined($arg0) ) {

                # simple hash of named values
                my %input_hash = @args;
                $rinput_hash = \%input_hash;
            }

        }
        elsif ( $reftype eq 'HASH' ) {
            $rinput_hash = $arg0;
        }
        else {
            carp "Unexpected ref type: $reftype\n";
            return;
        }
    }
    else {

        # no default 
        $rinput_hash = undef;
    }

    return unless defined($rinput_hash);

    my $Eoad     = $self->{_Eoad};
    my $rho_amb     = $self->{_rho_amb};
    my $gamma    = $self->{_gamma};
    my $symmetry = $self->{_symmetry};

    # STEP 1: determine the shock speed from the input parameter
    my $D;
    if ( defined( $rinput_hash->{D} ) ) {
        $D = $rinput_hash->{D};
    }
    elsif ( defined( $rinput_hash->{P} ) ) {
        my $P = $rinput_hash->{P};
        $D = sqrt( ( $gamma + 1 ) / 2 * $P / $rho_amb );
    }
    elsif ( defined( $rinput_hash->{U} ) ) {
        my $U = $rinput_hash->{U};
        $D = ( $gamma + 1 ) / 2 * $U;
    }
    elsif ( defined( $rinput_hash->{R} ) ) {
        my $R = $rinput_hash->{R};

        $D = 1.e99;
        if ( $R > 0 ) {

            # Reference: Korobeinikov top of p 67
            my $delta = 2 / ( $symmetry + 3 );
            my $pow = ( $symmetry + 1 );
            $D = $delta * sqrt( $Eoad / $R**$pow );
        }
    }
    elsif ( defined( $rinput_hash->{T} ) ) {
        my $T = $rinput_hash->{T};

        $D = 1.e99;
        if ( $T > 0 ) {

            # Reference: Korobeinikov top of p 67
            my $delta = 2 / ( $symmetry + 3 );
            my $pow = ( $symmetry + 1 );
            $D = $delta * ( $Eoad / $T**$pow )**( $delta / 2 );
        }
    }
    else {
        my @input_keys = qw(D P U R T );
        my %valid_input_keys;
        @valid_input_keys{@input_keys} = (1) x scalar(@input_keys);
        _check_keys( $rinput_hash, \%valid_input_keys,
            "Checking for valid input keys" );

        # shouldn't get here if keys okay
        croak "Programming error in shock_parameters\n";
    }

    # STEP 2: use the shock speed to determine the other parameters
    my $U   = 2 / ( $gamma + 1 ) * $D;
    my $P   = $rho_amb * $U * $D;
    my $rho = ( $gamma + 1 ) / ( $gamma - 1 ) * $rho_amb;
    my $R   = 1.e99;
    my $T   = 1.e99;
    if ( $D > 0 ) {
        my $delta = 2 / ( $symmetry + 3 );
        my $pow = ( $symmetry + 1 );
        $R = ( $delta**2 * $Eoad / $D**2 )**( 1 / $pow );
        $T = $R**( 1 / $delta ) / sqrt($Eoad);
    }
    my %return_hash = (
        'R'   => $R,
        'T'   => $T,
        'P'   => $P,
        'U'   => $U,
        'RHO' => $rho,
        'D'   => $D,
    );
    return \%return_hash;
}

sub shock_front_values {

    # given any one of these parameters:
    #  D - shock speed
    #  R - shock radius
    #  T - time
    #  P - shock pressure
    #  U - shock particle velocity

    # Returns [R, T, P, U, rho]

    my $rhash = shock_front_values_as_hash(@_);
    return [$rhash->{R}, $rhash->{T}, $rhash->{P}, $rhash->{U}, $rhash->{RHO}];
}


sub get_normalized_profile {
    my ( $self, $rr ) = @_;
    my $symmetry = $self->{_symmetry};
    my $gamma = $self->{_gamma};

    # Given 
    #  $symmetry=0,1,2 for plane, cylindrical, spherical
    #  $gamma=ideal gas gamma
    #  $rr = a list of r/rs radial coordinates (0 to 1)

    # Return:
    #  a list of the dimensionless profile of the point source similarity solution 
    #  with variables [ r/rs, t/ts, p/ps, u/us, rho/rhos ]
    # where
    #  rs=shock radius, 
    #  ts=time
    #  ps=shock pressure
    #  us=shock particle velocity
    #  rhos=shock density

    # Note: Whenever values of r/rs fall outside of [0,1] the profile values
    # will be set equal to zero

    # We will work with an array of increasing r; so reverse array if necessary 
    # and reverse again on exit
    return unless defined($rr) && @{$rr} > 0;
    my $is_reversed;
    if ( @{$rr} > 1 && $rr->[0] > $rr->[-1] ) {
        my @rev = reverse( @{$rr} );
        $rr          = \@rev;
        $is_reversed = 1;
    }

    # We will tweak gamma to avoid divides by zero
    if ( $gamma == 2 || $gamma == 7 && $symmetry == 2 ) { $gamma -= 1.e-10 }

    # Start by making a reference table of [theta, r/rs] values
    my $rxr;
    my $dx  = 1.e-15;
    my $geo = 1.03;
    for ( my $x = 0 ; $x < 1 ; $x += $dx ) {
        my $rrat = rrat( $x, $gamma, $symmetry );
        push @{$rxr}, [ $x, $rrat ];
        $dx *= $geo;
    }
    push @{$rxr}, [ 1, 1 ];

    #write_text_data( $rxr, "x\tr", "junk_table.txt" );

    # Now find the value of theta for each desirred radius of interest
    my $rgrid = _lookup_theta( $rr, $rxr, $symmetry, $gamma );

    # Now fill in the other variables as functions of theta
    # Special treatment is needed in the core where theta value will be undefined
    my $rprofile;
    my $gpow = ( $symmetry + 1 ) * $gamma / ( $gamma - 1 );
    my ( $x, $rrat, $prat, $urat, $rhorat, $uovr, $G );
    my $p_amb    = 0;
    my $RHOMIN   = 1.e-80;
    my $rho_to_G = sub {

        # Note: must not be called with rho=0
        my ( $rho, $rr, $pp ) = @_;
        if ( $rho <= 0 ) { $rho = $RHOMIN }
        $pp += $p_amb;
        my $rhopow = $rho**$gamma;
        if ( $rhopow <= 0 ) { $rhopow = 1.e-100; }
        my $kk  = $pp / $rhopow;
        my $Fkk = $kk * $rr**$gpow;
        return $Fkk;
    };
    my $G_to_rho = sub {

        # Invert the density interpolation function
        # Note: must not be called with r=0
        my ( $FF, $rr, $pp ) = @_;
        $pp += $p_amb;
        my $rho = 0;
        if ( $rr > 0 ) {
            my $kk = $FF / $rr**$gpow;
            $rho = ( $pp / $kk )**( 1 / $gamma );
        }
        return $rho;
    };
    my $imax = @{$rgrid} - 1;
    for ( my $i = $imax ; $i >= 0 ; $i-- ) {
        ( $x, $rrat ) = @{ $rgrid->[$i] };
        if ( defined($x) ) {
            $prat = prat( $x, $gamma, $symmetry );
            $urat = urat( $x, $gamma, $symmetry );
            $rhorat = rhorat( $x, $gamma, $symmetry );
            $uovr   = $urat / $rrat;
            $G      = $rho_to_G->( $rhorat, $rrat, $prat );
        }
	elsif ($rrat>1 || $rrat < 0) {
            $prat = 0;
            $urat = 0;
            $rhorat = 0;
	}
        else {

            # Core region..
            # prat is uniform in r at center
            # urat is linear in r
            # density goes to zero
            $urat = $rrat * $uovr;
            $rhorat = $G_to_rho->( $G, $rrat, $prat );
        }
        unshift @{$rprofile}, [ $rrat, 1, $prat, $urat, $rhorat ];
    }

    if ($is_reversed) {
        my @rev = reverse( @{$rprofile} );
        $rprofile = \@rev;
    }

    my $hdr = "r/rs\tt/ts\tp/ps\tu/us\trho/rhos";
    return wantarray ? ( $rprofile, $hdr ) : $rprofile;
}

sub _lookup_theta {
    my ( $rr, $rxr, $symmetry, $gamma ) = @_;

    # given a list of dimensionless radius ratios, 0 to 1, return a list with
    # [theta, radius] where theta is the implicit variable for that radius

    # Note: return theta undefined for r<rmin, where rmin=the first non-zero
    # radius in the reference table

    my $jmax = @{$rxr};
    my $imax = @{$rr};
    my $jhi  = $jmax - 1;
    my $jlo  = $jmax - 2;
    my $rgrid;
    for ( my $i = $imax - 1 ; $i >= 0 ; $i-- ) {
        my $r = $rr->[$i];
        if ( $r > 1 || $r < 0) {
            unshift @{$rgrid}, [ undef, $r ];
            next;
        }
        my ( $xhi, $rhi ) = @{ $rxr->[$jhi] };
        my ( $xlo, $rlo ) = @{ $rxr->[$jlo] };
        while ( $rlo > $r ) {
            if ( $jlo <= 0 ) { die "error in lookup\n" }
            ( $jhi, $xhi, $rhi ) = ( $jlo, $xlo, $rlo );
            $jlo--;
            ( $xlo, $rlo ) = @{ $rxr->[$jlo] };
        }

        if ( $jlo <= 0 ) {
            unshift @{$rgrid}, [ undef, $r ];
            next;
        }

        my $tolx  = 1.e-10 * abs($xhi-$xlo);
        my $tolr  = 1.e-10 * abs($rhi-$rlo);
        my $bvec  = [];
        my $ff_lo = $rlo - $r;
        my $ff_hi = $rhi - $r;
        my $ifconv;
        my $xx;
        ( $xx, $ifconv ) = nbrenti( $xlo, $ff_lo, $xhi, $ff_hi, $tolx, $bvec );

        # loop over iterations
        my $maxit = 100;
        my $iter;
        my $rrat;
        my $ff;
        for ( $iter = 1 ; $iter <= $maxit ; $iter += 1 ) {
            $rrat = rrat( $xx, $gamma, $symmetry );
            $ff = $rrat - $r;
            ( $xx, $ifconv ) = nbrentx( $ff, $bvec );
            if ( abs($ff) < $tolr  && $ifconv != 0) {    #$ifconv != 0 ) {
                #print STDERR "ff=$ff, xx=$xx, flo=$ff_lo, fhi=$ff_hi, xxlo=$xlo, xhi=$xhi, ifconv=$ifconv, tol=$tolr, iter=$iter\n";
                last;
            }
        }
        if ( $iter >= $maxit ) {

            # no convergence after maxit iterations -- shouldn't happen
            print STDERR "**warning-no convergence in brent iteration**; iter=$iter, ff=$ff\n";
        }
        unshift @{$rgrid}, [ $xx, $r ];
        ##print STDERR "i=$i, r=$r, jhi=$jhi, rhi=$rhi, jlo=$jlo, rlo=$rlo x=$xx, rrat=$rrat, ff=$ff, it=$iter\n";
    }
    return $rgrid;
}

sub alpha_integral {
    my ( $self, $tol, $itmax ) = @_;
    my $symmetry = $self->{_symmetry};
    my $gamma    = $self->{_gamma};
    my $num      = 10000;
    my ( $alpha, $err );
    my $it;

    # Patch for gamma=7 spherical symmetry
    if ( $symmetry == 2 && $gamma == 7 ) { $gamma -= $tol / 100 }

    # Patch for gamma=2
    if ( $gamma == 2 ) { $gamma -= $tol / 100 }

    for ( $it = 0 ; $it <= $itmax ; $it++ ) {
        my $alpha_last = $alpha;
        $num *= 2;
        my $dx = 1 / $num;
        my $coef = coef( $gamma, $symmetry );
        my $rxf;
        my $rxr;

        for ( my $i = 0 ; $i <= $num ; $i++ ) {
            my $x    = $i * $dx;
            my $f    = fofx( $x, $gamma, $symmetry );
            my $rrat = rrat( $x, $gamma, $symmetry );

            #push @{$rxr}, [$x, $rrat];
            push @{$rxf}, [ $rrat**( $symmetry + 1 ), $f ];
        }

        #write_text_data( $rxf, "rrat**N, f", "junkrxf.txt" );
        my $rxy = multiseg_integral($rxf);
        $alpha = $coef * $rxy->[-1]->[1];
        $err = ( $it > 0 ) ? abs( $alpha - $alpha_last ) : 10 * $tol;

        #print "$num\t$alpha\t$err\n";
        last if ( $err < $tol );

        #write_text_data( $rxy, "x,alpha", "junk1.txt" );
    }

    # Correct for plane symmetry
    if ( $symmetry == 0 ) { $alpha /= 2; $err /= 2 }

    return wantarray ? ( $alpha, $err, $it ) : $alpha;
}

##########################################
# Similarity Solution analytical functions 
##########################################

# The Von Neumann analytical solution is currently used to generate the
# profiles.  I switched to the Von Neumann version after encountering accuracy
# difficulties at low gamma values using with the Sedov analytical solution.

# References consulted:

# 1. Blast Wave, H.Bethe, K.Fuchs, J.O.Hirschfelder,J. Magee,
# R.Peirels,J.vonNeumann, Aug 1947, Los Alamos LA-2000.  Has Von Neumann's
# solution in spherical symmetry.

# 2. Sadek, H. S., Gottlieb, J. J., "Initial Decay Of Flow Properties Of
# Planar, Cylindrical And Spherical Blast Waves", 1983.  Institute for
# Aerospace Studies, University of Toronto. Has the Von Neumann solution
# extended to cylindrical and plane symmetry - but note a typo.

sub coef {

    # Leading coefficient
    my ( $gamma, $N ) = @_;
    my $pi = 4 * atan2( 1, 1 );
    my $c =
      ( 16 * ( ( 1 - $N ) * ( 2 - $N ) + 2 * $N * $pi ) ) /
      ( ( $N + 1 ) * ( $N + 3 )**2 * ( $gamma + 1 ) * ( $gamma - 1 ) );
    return $c;
}

sub rhorat {
    my ( $x, $gamma, $N ) = @_;
    my $t1 = $x; 
    my $p1 = ( $N + 1 ) / ( 2 * $gamma + $N - 1 );
    my $t2 = ( $x + $gamma ) / ( 1 + $gamma );
    my $p2 = ( -2 * $N ) / ( ( $N + 1 ) * $gamma - $N + 1 );

    my $t3_top = ( ( $N + 1 ) * ( 2 - $gamma ) * $x + 2 * $gamma + $N - 1 );
    my $t3_bot = ( ( 3 * $N + 1 ) - ( $N - 1 ) * $gamma );

    # avoid divide by zero at symmetry=2 gamma=7
    if ( $t3_bot == 0 ) {
        return 0;
    }
    my $t3 = $t3_top / $t3_bot;
    my $p3 = 1;
    if ( $t3 != 1 ) {
        $p3 =
          ( ( $N**2 + 2 * $N + 5 ) * $gamma**2 -
              ( 3 * $N**2 - 2 * $N - 1 ) * $gamma +
              4 * ( $N**2 - 1 ) ) /
          ( ( 2-$gamma ) *
              ( 2 * $gamma + $N - 1 ) *
              ( ( $N + 1 ) * $gamma - $N + 1 ) );
    }

    my $rhorat = $t1**$p1 * $t2**$p2 * $t3**$p3;
    return $rhorat;
}

sub prat {
    my ( $x, $gamma, $N ) = @_;
    my $t1 = ( $x + 1 ) / 2;
    my $p1 = 2 * ( $N + 1 ) / ( $N + 3 );
    my $t2 = ( $x + $gamma ) / ( 1 + $gamma );
    my $p2 = ( -2 * $N * $gamma ) / ( ( $N + 1 ) * $gamma - $N + 1 );

    my $t3_top = ( ( $N + 1 ) * ( 2 - $gamma ) * $x + 2 * $gamma + $N - 1 );
    my $t3_bot = ( ( 3 * $N + 1 ) - ( $N - 1 ) * $gamma );

    # avoid divide by zero at symmetry=2 gamma=7
    if ( $t3_bot == 0 ) {
        return 0;
    }
    my $t3 = $t3_top / $t3_bot;
    my $p3 = 1;
    if ( $t3 != 1 ) {
        $p3 =
          ( ( $N**2 + 2 * $N + 5 ) * $gamma**2 -
              ( 3 * $N**2 - 2 * $N - 1 ) * $gamma +
              4 * ( $N**2 - 1 ) ) /
          ( ( $N + 3 ) *
              ( 2 - $gamma ) *
              ( ( $N + 1 ) * $gamma - $N + 1 ) );
    }

    my $prat = $t1**$p1 * $t2**$p2 * $t3**$p3;
    return $prat;
}

sub urat {
    my ( $x, $gamma, $N ) = @_;
    my $t1 = $x;
    my $p1 = ( $gamma - 1 ) / ( 2 * $gamma + $N - 1 );
    my $t2 = ( $x + 1 ) / 2;
    my $p2 = ($N+1) / ( $N + 3 );
    my $t3 = ( $x + $gamma ) / ( 1 + $gamma );
    my $p3 = -$N * ( $gamma - 1 ) / ( ( $N + 1 ) * $gamma - $N + 1 );

    my $t4_top = ( ( $N + 1 ) * ( 2 - $gamma ) * $x + 2 * $gamma + $N - 1 );
    my $t4_bot = ( ( 3 * $N + 1 ) - ( $N - 1 ) * $gamma );

    # avoid divide by zero at symmetry=2 gamma=7
    if ($t4_bot == 0) {
	return 0;
    }
    my $t4 = $t4_top / $t4_bot;

    # avoid divide by zero at gamma=2
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

    my $urat = $t1**$p1 * $t2**$p2 * $t3**$p3 * $t4**$p4;
    return $urat;
}

sub rrat {
    my ( $x, $gamma, $N ) = @_;
    my $t1 = $x;
    my $p1 = ( $gamma - 1 ) / ( 2 * $gamma + $N - 1 );
    my $t2 = ( $x + 1 ) / 2;
    my $p2 = (-2) / ( $N + 3 );
    my $t3 = ( $x + $gamma ) / ( 1 + $gamma );
    my $p3 = ( $gamma + 1 ) / ( ( $N + 1 ) * $gamma - $N + 1 );

    # NOTE CORRECTION based on la2000; the -N+1 should have been +N-1:
    my $t4_top = ( ( $N + 1 ) * ( 2 - $gamma ) * $x + 2 * $gamma + $N - 1 );
    my $t4_bot = ( ( 3 * $N + 1 ) - ( $N - 1 ) * $gamma );

    # avoid divide by zero at symmetry=2 gamma=7
    if ($t4_bot == 0) {
	return 0;
    }
    my $t4 = $t4_top / $t4_bot;

    # avoid divide by zero at gamma=2
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

1;
