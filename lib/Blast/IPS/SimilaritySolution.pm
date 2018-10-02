package Blast::IPS::SimilaritySolution;

# This is a collection of software related to the very high pressure region of
# an ideal point source explosion where it may be described with the similarity
# solution.

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

use Blast::IPS::Utils qw(
  check_keys
);

use Blast::IPS::MathUtils qw(
  locate_2d
  multiseg_integral
  trapezoidal_integral
  parabolic_integral
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
    my ( $self, @args ) = @_;

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
        $rinput_hash = { symmetry => 2, gamma => 1.4, E0 => 1, rho_amb => 1.4 };
    }

    # Validate input keys
    check_keys( $rinput_hash, \%valid_input_keys,
        "Checking for valid input keys" );

    my $gamma    = $rinput_hash->{gamma};
    my $symmetry = $rinput_hash->{symmetry};
    my $E0       = $rinput_hash->{E0};
    my $rho_amb  = $rinput_hash->{rho_amb};

    # PATCH until coding revised to avoid the divide by gamma-2
    if ($gamma==2) {$gamma-=1.e-10}

    my $error = "";
    if ( !defined($E0) )       { $E0       = 1 }
    if ( !defined($gamma) )    { $gamma    = 1.4 }
    if ( !defined($rho_amb) )  { $rho_amb  = 1 }
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
    $alpha = 1 unless ($alpha);

    $self->{_E0}           = $E0;
    $self->{_rho_amb}      = $rho_amb;
    $self->{_gamma}        = $gamma;
    $self->{_symmetry}     = $symmetry;
    $self->{_Eoad}         = $E0 / ( $rho_amb * $alpha );
    $self->{_alpha}        = $alpha;
    $self->{_rtheta_table} = $self->make_reference_table(64);
    if ($error) { carp "$error" }
    return;
}
sub get_Eoad { $_[0]->{_Eoad} }
sub get_alpha { $_[0]->{_alpha} }
sub get_error { $_[0]->{_error} }
sub get_gamma { $_[0]->{_gamma} }
sub get_symmetry { $_[0]->{_symmetry} }
sub get_rho_amb { $_[0]->{_rho_amb} }
sub get_E0 { $_[0]->{_E0} }

sub shock_front_values_as_hash {
    my ( $self, @args ) = @_;

    # given any one of these shock front parameters:
    #  D - shock speed
    #  R - shock radius
    #  T - time
    #  P - shock pressure
    #  U - shock particle velocity

    # Return a hash of shock front values 

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
        check_keys( $rinput_hash, \%valid_input_keys,
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

    # given any one of these shock front parameters:
    #  D - shock speed
    #  R - shock radius
    #  T - time
    #  P - shock pressure
    #  U - shock particle velocity

    # Return [R, T, P, U, rho]

    my $rhash = shock_front_values_as_hash(@_);
    return [$rhash->{R}, $rhash->{T}, $rhash->{P}, $rhash->{U}, $rhash->{RHO}];
}

sub make_reference_table {

    # Make a reference table of [theta, r/rs] values approximately uniformly
    # spaced in r.  This will allow binary searches to be made over small
    # intervals.  

    my ( $obj, $max_count ) = @_;
    my $gamma    = $obj->get_gamma();
    my $symmetry = $obj->get_symmetry();
    my $count    = 0;
    my ( $jl, $ju );

    # We only search down to a very small value of theta>0.  This corresponds
    # to a point very close to the origin.  We can fill in the solution between
    # the origin and this point using analytical equations.
    my $tiny_theta = 1.e-20;

    # The point array stores each point in the order created.
    my $rpoint_table;
    my $store_point = sub {
	my ($rpoint)=@_;
	push @{$rpoint_table}, $rpoint;
        my $k = @{$rpoint_table}-1;
	return ($k);
    };
	
    # The segment table stores segments between points sorted by segment size
    my $rseg_table;
    my $insert_segment = sub {
	my ($km, $kp)=@_;

        # store the data for the segment between points k and km
        my $dr = $rpoint_table->[$km]->[1] - $rpoint_table->[$kp]->[1];
        if ( $dr < 0 ) { $dr = -$dr; ( $km, $kp ) = ( $kp, $km ) }
        my $item = [ $dr, $km, $kp ];

	if (!defined($rseg_table) || @{$rseg_table}<=0) {
	     push @{$rseg_table}, $item;
	}
        elsif ( @{$rseg_table} == 1 ) {

	    # we cannot call locate_2d with just one point
            if ( $dr < $rseg_table->[0]->[0] ) {
                unshift @{$rseg_table}, $item;
            }
            else {
                push @{$rseg_table}, $item;
            }
        }
	else { 
            ($jl, $ju) = locate_2d($dr, 0, $rseg_table, $jl, $ju);
            if ( !defined($jl) || $jl < 0 ) {
                unshift @{$rseg_table}, $item;
            }
	    elsif (!defined($ju) || $ju>@{$rseg_table}) {
		push @{$rseg_table}, $item;
	    }
	    else {
	       splice @{$rseg_table}, $ju, 0, $item;
            }
	}
    };

    my $get_values = sub {

	# get the profile values at a specific theta
	my ($theta)=@_;
        my $rrat = rrat( $theta, $gamma, $symmetry );
        my $prat = prat( $theta, $gamma, $symmetry );
        my $urat = urat( $theta, $gamma, $symmetry );
        my $rhorat = rhorat( $theta, $gamma, $symmetry );
        return [ $theta, $rrat, $prat, $urat, $rhorat ];
    };

    # Start the point table 
    my $kz = $store_point->($get_values->(0));
    my $km = $store_point->($get_values->($tiny_theta));
    my $kp = $store_point->($get_values->(1));

    # Start the segment table. Note that we are ignoring the segment between
    # theta=0 and theta = tiny_theta on purpose.  
    $insert_segment->( $km, $kp );

    # loop to keep splitting the largest segment until we have enough points.
    # This uses a very fast, clean algorithm, with one array and one stack.
    while(1) {

	# Get the largest segment. It runs from index kl to ku.
        my $big=pop @{$rseg_table};
	my ($dr, $kl, $ku)=@{$big};

	# pull out the two points
	my ($tl, $rl) = @{$rpoint_table->[$kl]};
	my ($tu, $ru) = @{$rpoint_table->[$ku]};

	# use logarithmic interpolation to find the theta which gives the mid
	# point in r. This approximation gets increasingly better as points
        # are added.
        my $rmid = 0.5 * ( $rl + $ru );
        my $slope = log( $tu / $tl ) / log( $ru / $rl );
        my $theta = $tl * exp( $slope * log( $rmid / $rl ) );

	# make and store the new point
        my $kk = $store_point->($get_values->($theta));

	# Stop if enough points
	$count++;
	last if ($count>=$max_count);

	# put the two new segments into the dr stack and keep going
	$insert_segment->($kl, $kk);
	$insert_segment->($kk, $ku);
    }
    
    # Since the points are stored in the order found, we have to sort them
    my @sorted=sort {$a->[1] <=> $b->[1]} @{$rpoint_table};

    return \@sorted;
}

sub get_normalized_point {
    my ( $self, $r, $flag ) = @_;
    my $rr;
    push @{$rr}, $r;
    my ($rarray, $hdr) = $self->get_normalized_profile($rr, $flag);
    my $num=@{$rarray};
    return wantarray ? ( $rarray->[0], $hdr ) : $rarray->[0];
}

sub get_normalized_profile {
    my ( $self, $rr, $flag ) = @_;
    my $symmetry = $self->{_symmetry};
    my $gamma = $self->{_gamma};

    # Given 
    #  $symmetry=0,1,2 for plane, cylindrical, spherical
    #  $gamma=ideal gas gamma
    #  $rr = a list of r/rs radial coordinates (0 to 1)
    #  $flag = 1 to include theta as an additional variable

    # Return:
    #  a list of the dimensionless profile of the point source similarity
    #  solution with variables [ r/rs, t/ts, p/ps, u/us, rho/rhos ]
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

    my $rprofile;

    # handle special case
    if ( $symmetry == 2 && $gamma == 7 ) {
        foreach my $rrat ( @{$rr} ) {
            my ( $trat, $prat, $urat, $rhorat ) = ( 1, 0, 0, 0 );
            if ( $rrat > 0 && $rrat <= 1 ) {
                $urat   = $rrat;
                $rhorat = $rrat;
                $prat   = $rrat**3;
            }
            unshift @{$rprofile}, [ $rrat, $trat, $prat, $urat, $rhorat ];
        }
    }
    else {

        # find the value of theta for each desired radius of interest
        my $rgrid = $self->_lookup_theta_profile( $rr); 

        $rprofile = $self->_fill_profile($rgrid, $flag);
    }

    if ($is_reversed) {
        my @rev = reverse( @{$rprofile} );
        $rprofile = \@rev;
    }

    my $hdr = "r/rs\tt/ts\tp/ps\tu/us\trho/rhos";
    return wantarray ? ( $rprofile, $hdr ) : $rprofile;
}

sub _fill_profile {
    my ( $self, $rgrid, $flag ) = @_;

    # Given an array of theta values [theta, r/rs], create a complete profile
    # of normalized values
    #     [ $rrat, 1, $prat, $urat, $rhorat ];  
    # or, if flag:
    #     [ $rrat, 1, $prat, $urat, $rhorat, $theta ];  
    # Special treatment is needed in the core where theta value will be
    # undefined

    my $symmetry       = $self->{_symmetry};
    my $gamma          = $self->{_gamma};
    my $rtheta_table   = $self->{_rtheta_table};

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

	    # Note that we are re-defining rrat to be correct for the
	    # given theta value
            $rrat = rrat( $x, $gamma, $symmetry );
            $urat = urat( $x, $gamma, $symmetry );
            $rhorat = rhorat( $x, $gamma, $symmetry );
            $uovr   = $urat / $rrat if ($rrat>0);
            $G      = $rho_to_G->( $rhorat, $rrat, $prat );
        }
        elsif ( $rrat > 1 || $rrat < 0) {
            $prat   = 0;
            $urat   = 0;
            $rhorat = 0;
        }
        else {

            # Backup logic for Core region..in case theta is not set
            # prat is uniform in r at center
            # urat is linear in r
            # density goes to zero
            $urat = $rrat * $uovr;
            $rhorat = $G_to_rho->( $G, $rrat, $prat );
        }
	my $item;
        if ($flag) { $item= [ $rrat, 1, $prat, $urat, $rhorat, $x ] }
        else { $item= [ $rrat, 1, $prat, $urat, $rhorat ] }
        unshift @{$rprofile}, $item;
    }
    return $rprofile;
}

sub _lookup_theta_profile {

    my ( $self, $rr) = @_; #, $rxr, $symmetry, $gamma ) = @_;
    my $symmetry       = $self->{_symmetry};
    my $gamma          = $self->{_gamma};
    my $rxr            = $self->{_rtheta_table};

    # given a list of dimensionless radius ratios, 0 to 1, return a list with
    # [theta, radius] where theta is the implicit variable for that radius

    # Experimental: Setup for center slopes to find center theta values
    # This will allow us to define theta in the core region
    my ($dlnt_dlnr, $lntu, $lnru);
    if ( @{$rxr} > 2 ) {
        my ( $tl, $rl) = @{ $rxr->[1] };
        my ( $tu, $ru) = @{ $rxr->[2] };
        my $lntl = log($tl);
        my $lnrl = log($rl);
        $lnru = log($ru);
        $lntu = log($tu);
        $dlnt_dlnr = ( $lntu - $lntl ) / ( $lnru - $lnrl );
    }

    # Note: return theta undefined for r<rmin, where rmin=the first non-zero
    # radius in the reference table
    my $jmax = @{$rxr} -1;
    my $imax = @{$rr} -1;
    my $jhi  = $jmax;
    my $jlo  = $jmax - 1;

    # Find the starting interval
    ($jlo, $jhi) = locate_2d($rr->[$imax], 1, $rxr, $jlo, $jhi);

    my ($rgrid, $xlast, $rlast);

    # Loop over points to fill
    my $r;
    for ( my $i = $imax ; $i >= 0 ; $i-- ) {
	my $rlast=$r;
        $r = $rr->[$i];
        if ( $r > 1 || $r < 0) {
	    $xlast = undef;
	    $rlast = $r;
            unshift @{$rgrid}, [ $xlast, $r ];
            next;
        }

        if ( $jhi >= @{$rxr} ) {
            $jhi = $jlo;
            $jlo = $jlo - 1;
        }

        my ( $xhi, $rhi ) = @{ $rxr->[$jhi] };

	# optimization: the previous computed point will form a better upper
	# bound than the previous table value if we stay in the same table
	# interval
        if (   defined($xlast)
            && defined($rlast)
            && $rlast > $r
            && $rlast < $rhi )
        {
            $xhi = $xlast;
            $rhi = $rlast;
        }

        my ( $xlo, $rlo ) = @{ $rxr->[$jlo] };
	my $count=0;
        while ( $rlo > $r ) {
	    $count++;
            if ( $jlo <= 0 ) { die "error in lookup\n" }
            ( $jhi, $xhi, $rhi ) = ( $jlo, $xlo, $rlo );
            $jlo--;
            ( $xlo, $rlo ) = @{ $rxr->[$jlo] };
        }

        if ( $jlo <= 0 ) {

	    # We are in the core, before the first table point.  Try to
	    # define theta using the slope set above.
	    my $theta;

            if ( $r <= 0 ) { $theta = 0 }
            elsif ($dlnt_dlnr) {
                my $lnr = log($r);
                my $lnt = $lntu + ( $lnr - $lnru ) * $dlnt_dlnr;
                $theta = exp($lnt);
            }
            unshift @{$rgrid}, [ $theta, $r ];
            next;
        }

        my $bvec  = [];
        my $ff_lo = $rlo - $r;
        my $ff_hi = $rhi - $r;
        my $FF_lo = log( $rlo / $r );
        my $FF_hi = log( $rhi / $r );
        my $Xlo   = log($xlo);
        my $Xhi   = log($xhi);

        my $dr = $rhi - $rlo;
	my $DR = log($rhi/$rlo);
        if ( defined($rlast) ) {
            $dr = abs( $rlast - $r);
	    $DR = log($rlast/$r);
        }

	# These tolerances should be set pretty tight. But some non-convergence
	# messages may occur even when the root is found fairly accurately.
        my $tolf = 1.e-8 * $dr;                  ##abs( $rhi - $rlo );
        my $tolF = 1.e-8 * abs($DR);             ##abs( $FF_hi - $FF_lo );

        my $tolx = 1.e-12 * abs( $xhi - $xlo );
        my $tolX = 1.e-12 * abs( $Xhi - $Xlo );

	# This version has the option of either using log or linear searches
	# log searches seem best
        my $use_log_search = 1; 
        ##my $use_log_search = ( $r < 0.95 );

	my $ifconv;
	my $XX;
	my $xx;

        if ($use_log_search) {
            ( $XX, $ifconv ) =
              nbrenti( $Xlo, $FF_lo, $Xhi, $FF_hi, $tolX, $bvec );
	    $xx=exp($XX);
        }
        else {
            ( $xx, $ifconv ) =
              nbrenti( $xlo, $ff_lo, $xhi, $ff_hi, $tolx, $bvec );
        }

        # loop over iterations
        my $maxit = 100;
        my $iter;
        my $rrat;
        my $FF;
	my $ff;
        for ( $iter = 1 ; $iter <= $maxit ; $iter += 1 ) {
            $rrat = rrat( $xx, $gamma, $symmetry );
            $ff = $rrat - $r; 

            if ($use_log_search) {
                $FF = log( $rrat / $r );
                ( $XX, $ifconv ) = nbrentx( $FF, $bvec );
                $xx = exp($XX);
            }
	    else {
                ( $xx, $ifconv ) = nbrentx( $ff, $bvec );
	    }

            ##if ( abs($FF) < $tolF  && $ifconv != 0) {    #$ifconv != 0 ) {
            ##if ( abs($ff) < $tolf  && $ifconv != 0) {    #$ifconv != 0 ) {

	    # It seems to work best to put a tight tolerance on x and
	    # a somewhat lower tolerance of f and then check both.
            #if (   $ifconv != 0 ) {
            if ( abs($ff) < $tolf  && $ifconv != 0 ) {
                #print STDERR "ff=$ff, xx=$xx, flo=$ff_lo, fhi=$ff_hi, xxlo=$xlo, xhi=$xhi, ifconv=$ifconv, tol=$tolr, iter=$iter\n";
                last;
            }
        }
        if ( $iter >= $maxit ) {

            # no convergence after maxit iterations -- shouldn't happen
            print STDERR "**warning-no convergence in brent iteration**; iter=$iter, for r=$r, last x=$xx, xlo=$xlo, xhi=$xhi, ff=$ff tolF=$tolF, tolX=$tolX, tolx=$tolx\n";
        }
        unshift @{$rgrid}, [ $xx, $r ];
	$xlast = $xx;
	$rlast = $r;

        ##print STDERR "i=$i, r=$r, jhi=$jhi, rhi=$rhi, jlo=$jlo, rlo=$rlo x=$xx, rrat=$rrat, ff=$ff, it=$iter\n";
    }
    return $rgrid;
}

sub alpha_integral {

    # Evaluate the energy integral to obtain the parameter alpha
    my ( $self, $tol, $itmax, $tiny_theta ) = @_;

    # tol = stopping tolerance on absolute value of alpha
    $tol = 1.e-13;

    # itmax = max iterations in case the tolerance is not reached; 
    $itmax = 15 unless defined($itmax);

    # We are integrating the von Neumann solution from theta=0 to 1
    # We split the integral into two parts, 
    #  Part 1: from theta=0 to tiny_theta
    #  Part 2: from theta=tiny_theta to 1

    # tiny_theta should be a very small value of theta. The recommended value
    # is 1.e-40.  A value 1.e-40 worked well for gamma down to 1.001.
    $tiny_theta = 1.e-40 unless defined($tiny_theta);

    my $symmetry     = $self->{_symmetry};
    my $gamma        = $self->{_gamma};
    my $rtheta_table = $self->{_rtheta_table};
    my ( $alpha, $err );
    my $it;

    # Patch for gamma=7 spherical symmetry
    # We cannot evaluate at gamma exactly 7, so we reduce gamma slightly
    if ( $symmetry == 2 && $gamma == 7 ) { $gamma -= $tol / 100 }

    # Patch for gamma=2:
    # Until the functions can handle gamma=2, we reduce gamma slightly
    if ( $gamma == 2 ) { $gamma -= $tol / 100 }

    my $pow        = $symmetry + 1;

    # Get the leading coefficient
    my $coef = coef( $gamma, $symmetry );

    my $get_values = sub {

        # get the volume and integrand at a specific theta
        my ($x) = @_;
        my $r = rrat( $x, $gamma, $symmetry );
        my $V = $r**$pow;
        if ( $V > 1 ) { $V = 1 }
        my $f = fofx( $x, $gamma, $symmetry );
        return [ $V, $f, $x ];
    };

    my $rVfx;

    #########################
    # Do the Part 1 integral:
    #########################
    # Use the trapezoidal rule to integrate over the first tiny segment
    # I0 = integral from theta=0 to tiny_theta
    # dI0 = estimated error of this integral
    my $rVfx0 = $get_values->(0);
    my $rVfx1 = $get_values->($tiny_theta);
    my ( $V0, $f0, $x0 ) = @{$rVfx0};
    my ( $V1, $f1, $x1 ) = @{$rVfx1};
    my $dV = $V1 - $V0;
    my $I0 = 0.5 * ( $f0 + $f1 ) * $dV;

    # The error is estimated by taking the difference of the integral made with
    # either endpoint as f.  That is:
    #   dI0 = $f1*$dV - $f2*$dV, so
    my $err0 = $coef * ( $f1 - $f0 ) * $dV;

    #########################
    # Do the Part 2 integral:
    #########################

    # Start with 2 points
    push @{$rVfx}, $rVfx1;
    push @{$rVfx}, $get_values->(1);

    # Loop to keep dividing by two until the error is small enough
    $alpha = 0;
    for ( $it = 0 ; $it <= $itmax ; $it++ ) {
        my $alpha_last = $alpha;
        my $rxy        = parabolic_integral($rVfx, $I0);
        $alpha = $coef * $rxy->[-1]->[1];
        $err = ( $it > 0 ) ? abs( $alpha - $alpha_last ) : 10 * $tol;

        print "$it\t$alpha\t$err\n";

        last if ( $err < $tol );

	# Double the number of integration points and continue...
	# Divide each segment into approximately equal energy increments 
        my $rVfx_fine;
        for ( my $i = 0 ; $i < @{$rVfx} - 1 ; $i++ ) {
            push @{$rVfx_fine}, $rVfx->[$i];
            my ( $Vl, $fl, $xl ) = @{ $rVfx->[$i] };
	    my $Il = $rxy->[$i]->[1];
            my ( $Vu, $fu, $xu ) = @{ $rVfx->[ $i + 1 ] };
	    my $Iu = $rxy->[$i+1]->[1];
            my $I     = 0.5 * ( $Il + $Iu );
            my $slope = log( $xu / $xl ) / log( $Iu / $Il );
            my $x     = $xl * exp( $slope * log( $I / $Il ) );
            push @{$rVfx_fine}, $get_values->($x);
        }
        push @{$rVfx_fine}, $rVfx->[-1];
        $rVfx = $rVfx_fine;
    }

    # include the starting error
    $err += $err0;

    # Use the trapezoidal method to get another estimate
    # this can be used to give an upper bound estimate of the error
    my $rxy_trap = trapezoidal_integral( $rVfx, $I0 );
    my $alpha_trap = $coef * $rxy_trap->[-1]->[1];

    # Correct for plane symmetry
    if ( $symmetry == 0 ) { $alpha /= 2; $err /= 2 ; $alpha_trap/=2}

    # Return values:
    #  $alpha = the best estimate of alpha (using simpsons method)
    #  $err   = the best error estimate (difference in last two iterations
    #           using simpsons method
    #  $it    = number of iterations
    #  $alpha_trap = the value of alpha computed with trapezoidal integration
    #  $err0  = the error in the first of two parts (this is included in $err
    #           but also returned separately so that it can be checked)
    return wantarray ? ( $alpha, $err, $it, $alpha_trap, $err0 ) : $alpha;
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
    if ( $gamma == 2 || $gamma == 7 && $N == 2 ) { $gamma -= 1.e-10 }

    my $t1 = $x; 
    my $p1 = ( $N + 1 ) / ( 2 * $gamma + $N - 1 );
    my $t2 = ( $x + $gamma ) / ( 1 + $gamma );
    my $p2 = ( -2 * $N ) / ( ( $N + 1 ) * $gamma - $N + 1 );

    my $t3_top = ( ( $N + 1 ) * ( 2 - $gamma ) * $x + 2 * $gamma + $N - 1 );
    my $t3_bot = ( ( 3 * $N + 1 ) - ( $N - 1 ) * $gamma );

    # shouldn't happen: avoid divide by zero at symmetry=2 gamma=7
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
    if ( $gamma == 2 || $gamma == 7 && $N == 2 ) { $gamma -= 1.e-10 }
    my $t1 = ( $x + 1 ) / 2;
    my $p1 = 2 * ( $N + 1 ) / ( $N + 3 );
    my $t2 = ( $x + $gamma ) / ( 1 + $gamma );
    my $p2 = ( -2 * $N * $gamma ) / ( ( $N + 1 ) * $gamma - $N + 1 );

    my $t3_top = ( ( $N + 1 ) * ( 2 - $gamma ) * $x + 2 * $gamma + $N - 1 );
    my $t3_bot = ( ( 3 * $N + 1 ) - ( $N - 1 ) * $gamma );

    # shouldn't happen: avoid divide by zero at symmetry=2 gamma=7
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
    if ( $gamma == 2 || $gamma == 7 && $N == 2 ) { $gamma -= 1.e-10 }
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

    # shouldn't happen: avoid divide by zero at gamma=2
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
    if ( $gamma == 2 || $gamma == 7 && $N == 2 ) { $gamma -= 1.e-10 }
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
    if ( $gamma == 2 || $gamma == 7 && $N == 2 ) { $gamma -= 1.e-10 }

    # integrate this from x=0 to x=1
    my $p1 = ( 3 * $N + 5 ) / ( $N + 3 );
    my $t1 = ( $x + 1 ) / 2;
    my $p2 = ( -2 * $N * $gamma ) / ( ( $N + 1 ) * $gamma - $N + 1 );
    my $t2 = ( $x + $gamma ) / ( 1 + $gamma );

    # NOTE CORRECTION based on la2000; the -N+1 should have been +N-1:
    my $t3_top = ( ( $N + 1 ) * ( 2 - $gamma ) * $x + 2 * $gamma + $N - 1 );
    my $t3_bot = ( ( 3 * $N + 1 ) - ( $N - 1 ) * $gamma );

    # shouldn't happen: avoid divide by zero at symmetry=2 gamma=7
    if ($t3_bot == 0) {
	return 0;
    }
    my $t3 = $t3_top / $t3_bot;


    # shouldn't happen: avoid divide by zero at gamma=2
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
