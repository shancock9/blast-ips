package Blast::IPS::Medium;
use Blast::IPS::Utils qw(check_keys);

# This is a collection of utilities related to the overall medium and geometry
# of an explosion in air

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

my $VERSION = 1.00;
use Carp;
use strict;

sub new {
    my ( $class, @args ) = @_;
    my $self = {};
    bless $self, $class;
    $self->_setup(@args);
    return $self;
}

sub _setup {
    my ( $self, @args ) = @_;

    my @input_keys = qw(
      symmetry
      ASYM
      gamma
      p_amb
      sspd_amb
      rho_amb
      E0
      ground_plane
    );
    my %valid_input_keys;
    @valid_input_keys{@input_keys} = (1) x scalar(@input_keys);

    # Convert the various input methods to a hash ref
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

        # default to spherical with gamma 1.4 if no arg given
        $rinput_hash = { symmetry => 2, gamma => 1.4 };
    }

    # Validate input keys
    check_keys( $rinput_hash, \%valid_input_keys,
        "Checking for valid input keys" );

    # The following four quantities can be specified:
    my $gamma        = $rinput_hash->{gamma};
    my $symmetry     = $rinput_hash->{symmetry};
    my $p_amb        = $rinput_hash->{p_amb};
    my $sspd_amb     = $rinput_hash->{sspd_amb};
    my $rho_amb      = $rinput_hash->{rho_amb};
    my $ground_plane = $rinput_hash->{ground_plane};
    my $E0           = $rinput_hash->{E0};

    # Allow older symbol 'ASYM' instead of 'symmetry' for compatability
    if ( !defined($symmetry) ) { $symmetry = $rinput_hash->{ASYM} }

    if ( !defined($p_amb) )        { $p_amb        = 1 }
    if ( !defined($sspd_amb) )     { $sspd_amb     = 1 }
    if ( !defined($rho_amb) )      { $rho_amb      = $gamma }
    if ( !defined($E0) )           { $E0           = 1 }
    if ( !defined($ground_plane) ) { $ground_plane = $symmetry == 0 ? 1 : 0 }

    if ( $symmetry =~ /^(SCP)/i ) {
        my $symmetry_name = uc $1;
        $symmetry =
            $symmetry_name eq 'S' ? 2
          : $symmetry_name eq 'C' ? 1
          :                         0;
    }

    if ( $symmetry !~ /^[012]$/ ) {
        confess <<EOM;
------------------------------------------------------------------------
Unknown Symmetry Flag: $symmetry; Expecting one of: (0, 1, 2)
------------------------------------------------------------------------

EOM

    }
    if ( $p_amb < 0 ) {
        confess <<EOM;
------------------------------------------------------------------------
Expecting non-negative atmospheric pressure but got '$p_amb'
------------------------------------------------------------------------
EOM
    }
    if ( $rho_amb <= 0 ) {
        confess <<EOM;
------------------------------------------------------------------------
Expecting positive atmospheric density but got '$rho_amb'
------------------------------------------------------------------------
EOM
    }
    if (!defined($gamma)) {
        if ($p_amb) {
            $gamma = $rho_amb * $sspd_amb**2 / $p_amb;
        }
	else {
	    $gamma=1.4;
	}
        $self->{_gamma} = $gamma;
    }
    if ( $gamma <= 1 ) {
        confess <<EOM;
------------------------------------------------------------------------
Expecting gamma>1 but got '$gamma'
------------------------------------------------------------------------
EOM
    }

    $self->{_gamma}        = $gamma;
    $self->{_symmetry}     = $symmetry;
    $self->{_p_amb}        = $p_amb;
    $self->{_sspd_amb}     = $sspd_amb;
    $self->{_rho_amb}      = $rho_amb;
    $self->{_ground_plane} = $ground_plane;
    $self->{_E0}           = $E0;

    ##################################################
    # FIXME: To be deleted.  Need to compute on the fly because
    # any of p, rho, gamma can change after creating the object
    my $A_amb = $p_amb / $rho_amb**$gamma;
    $self->{_A_amb} = $A_amb;

    # Define a reference value for A which is the ambient value if positive and
    # is 1 if A_amb is zero. This allows a reference value to be used
    # for all types of runs.  
    my $A_ref    = $A_amb;
    if ($A_ref<=0) { $A_ref    = 1 }
    $self->{_A_ref} = $A_ref;
    ##################################################

    return $self;
}

sub get_A_parameters {
    my ($self)  = @_;
    my $p_amb   = $self->{_p_amb};
    my $rho_amb = $self->{_rho_amb};
    my $gamma   = $self->{_gamma};
    my $A_amb   = $p_amb / $rho_amb**$gamma;

    # Define a reference value for A which is the ambient value if positive and
    # is 1 if A_amb is zero. This allows a reference value to be used
    # for all types of runs.
    my $A_ref = $A_amb;
    if ( $A_ref <= 0 ) { $A_ref = 1 }
    return ( $A_amb, $A_ref );
}

# Special get_ and set_ methods go here
sub get_symmetry {
    my ($self) = @_;
    return $self->{_symmetry};
}
sub get_gamma {
    my ($self) = @_;
    return $self->{_gamma};
}
sub get_p_amb {
    my ($self) = @_;
    return $self->{_p_amb};
}
sub get_rho_amb {
    my ($self) = @_;
    return $self->{_rho_amb};
}
sub get_sspd_amb {
    my ($self) = @_;
    return $self->{_sspd_amb};
}
sub get_E0{
    my ($self) = @_;
    return $self->{_E0};
}
sub get_ground_plane{
    my ($self) = @_;
    return $self->{_ground_plane};
}

sub set_symmetry {
    my ($self, $val) = @_;
    $self->{_symmetry} = $val;
    return
}
sub set_gamma {
    my ($self, $val) = @_;
    $self->{_gamma} = $val;
    return;
}
sub set_p_amb {
    my ($self, $val) = @_;
    $self->{_p_amb}=$val;
    return;
}
sub set_rho_amb {
    my ($self, $val) = @_;
    $self->{_rho_amb} = $val;
    return
}
sub set_sspd_amb {
    my ($self, $val) = @_;
    $self->{_sspd_amb} = $val;
    return
}
sub set_E0{
    my ($self, $val) = @_;
    $self->{_E0} = $val;
    return
}
sub set_ground_plane{
    my ($self, $val) = @_;
    $self->{_ground_plane}=$val;
    return;
}

sub get_kk {
    my ($self) = @_;
    my $gamma = $self->{_gamma};
    my $kk = ( $gamma + 1 ) / ( 2. * $gamma );
    return ($kk);
}

sub ushock_from_sigma {
    my ( $self, $sigma, $log_Arat ) = @_;
    my $pp = $self->ovprat_from_sigma( $sigma, $log_Arat );
    my ( $ushock, $upshock, $ushockx, $dup_dovp ) =
      $self->ushock_from_ovp($pp);
      #$self->ushock_from_ovprat($pp);
    return ( $ushock, $upshock, $ushockx, $dup_dovp );
}

sub ushock_from_ovp {
    my ( $self, $ovp ) = @_;
    my $p_amb = $self->{_p_amb};
    my $sspd_amb = $self->{_sspd_amb};
    if ($p_amb>0 && $sspd_amb>0) {
        return $self->ushock_from_ovprat( $ovp/$p_amb );
    }
    else {
        return $self->ushock_from_pabs( $ovp + $p_amb );
    }
}

sub ushock_from_pabs {
    my ( $self, $pabs ) = @_;

    # shock speed from absolute pressure
    if ($pabs<0) {
	croak "pabs=$pabs\n";
	$pabs=0;
    }
    my $p_amb = $self->{_p_amb};
    my $sspd_amb = $self->{_sspd_amb};
    my ( $ushock, $upshock, $ushockx, $dup_dovp, $rhoshock );
    if ($p_amb>0 && $sspd_amb>0) {
        my $ovprat = ( $pabs - $p_amb ) / $p_amb;
        ( $ushock, $upshock, $ushockx, $dup_dovp, $rhoshock ) =
          $self->ushock_from_ovprat($ovprat);
    }
    else {
        my $gamma   = $self->{_gamma};
        my $rho_amb = $self->{_rho_amb};

        # FIXME!!
        $ushock  = 0;
        $upshock = 0;
        my $gterm = 2 / ( $gamma + 1 );
        if ( $pabs > 0 && $rho_amb > 0 ) {
            $ushock = sqrt( $pabs / ( $gterm * $rho_amb ) );
            $upshock = $gterm * $ushock;
        }
        $ushockx = $ushock;

  	# OLD:WRONG
        # my $dpdu=($gamma+1)*$rho_b*$upshock;
        # But dovp_du=dpdu/p_amb=infinity so du/dovp=0
        # $dup_dovp = 0;

        $dup_dovp =
          ( 1 + ( $sspd_amb / $ushock )**2 ) / ( 2 * $rho_amb * $ushock );
        $rhoshock = ( $gamma + 1 ) / ( $gamma - 1 ) * $rho_amb;
    }
    return ( $ushock, $upshock, $ushockx, $dup_dovp, $rhoshock );
}

sub strong_shock_from_pabs {

    # Might be useful in case you want to evaluate a strong shock
    # for an atmosphere with P_amb>0
    my ( $self, $p_shock_e ) = @_;
    my $gamma   = $self->{_gamma};
    my $rho_amb = $self->{_rho_amb};
    my $gterm   = 2 / ( $gamma + 1 );
    my $u_shock_e  = 0;
    my $up_shock_e = 0;
    if ( $rho_amb > 0 && $p_shock_e > 0 ) {
        $u_shock_e = sqrt( $p_shock_e / ( $gterm * $rho_amb ) );
        $up_shock_e = $gterm * $u_shock_e;
    }
    my $rho_shock_e      = ( $gamma + 1 ) / ( $gamma - 1 ) * $rho_amb;
    my $usmc_e           = $u_shock_e;
    my $dup_dovp_shock_e = 0;
    return ( $u_shock_e, $up_shock_e, $usmc_e, $dup_dovp_shock_e,
        $rho_shock_e );
}

sub ushock_from_ovprat {
    my ( $self, $ovprat ) = @_;

    # shock speed from overpressure
    my $gamma    = $self->{_gamma};
    my $sspd_amb = $self->{_sspd_amb};
    my $rho_amb  = $self->{_rho_amb};
    my $p_amb    = $self->{_p_amb};
    if ( $sspd_amb <=0 || $p_amb <=0 || $rho_amb <=0) {
	print STDERR "ushock from ovprat: sspd_amb=$sspd_amb, p_amb=$p_amb, rho_amb=$rho_amb\n";
        my ( $a, $b, $c ) = caller();
        die "$a, $b, $c\n";
    }
    my $term     = $ovprat * ( $gamma + 1 ) / ( 2 * $gamma );
    my $mshock   = sqrt( 1 + $term );
    my $ushock   = $sspd_amb * $mshock;
    my $ushockx  = $ushock - $sspd_amb;

    # Taylor series expansion of the square root for small overpressures
    # Switch at a point which keeps error < about 1.e-17
    if ( $term < 1.e-4 ) {
        $ushockx = $term * ( 0.5 + $term * ( -0.125 + $term / 16 ) );
    }

    # shock particle velocity
    my $upshock = $ovprat * $p_amb / ( $rho_amb * $ushock );

    # dup/dovp
    my $dup_dovp =
      $sspd_amb *
      $sspd_amb / $gamma / $ushock *
      ( 1 - ( $gamma + 1 ) / 4 * $upshock / $ushock );

    my $pabs = ( $ovprat + 1 ) * $p_amb;
    my $top  = ( $gamma + 1 ) * $pabs + ( $gamma - 1 ) * $p_amb;
    my $bot  = ( $gamma - 1 ) * $pabs + ( $gamma + 1 ) * $p_amb;
    my $rhoshock;
    if   ( $bot > 0 ) { $rhoshock = $rho_amb * $top / $bot }
    else              { $rhoshock = $rho_amb }
    return ( $ushock, $upshock, $ushockx, $dup_dovp, $rhoshock );
}

sub pdc_from_sigma {

    my ( $self, $sigma, $lnArat ) = @_;

    # Given sigma=2/(gamma-1)*(sspd-sspd_amb)
    #       lnArat = log(A)-log(Aref)
    # return primitive state variables
    #       overdensity, overpressure, oversspd and their absolute values
    $lnArat = 0 unless ( defined($lnArat) );

    my $gamma    = $self->{_gamma};
    my $sspd_amb = $self->{_sspd_amb};
    my $rho_amb  = $self->{_rho_amb};
    my $p_amb    = $self->{_p_amb};
    my $A_amb    = $self->{_A_amb};
    my $A_ref    = $self->{_A_ref};

    # We have
    #  sspd=$sspd_amb+$sspdx
    #  sspd**2=gamma*p/rho
    #  A=p/rho**gamma
    # so dividing these last to to eliminate p we get
    #  rho=(sspd*82/(gamma*A))**(1/(gamma-1))
    # Take log to get the result

    my $sspdx = $sigma * ( $gamma - 1 ) / 2;
    my $lnA   = $lnArat + log($A_ref);
    my $sspd  = $sspd_amb + $sspdx;
    my ( $ovd, $rho );
    my $eps = 1.e-4;
    if ( $sspd_amb > 0 && abs($sspdx) < $eps * $sspd_amb ) {
        my $x      = $sspdx / $sspd_amb;
        my $p=2/($gamma-1);
        my $series = $x * $p *
          ( 1 + $x * ( ( $p - 1 ) / 2 * ( 1 + $x * ( $p - 2 ) / 3 ) ) );
        my $A = exp($lnA);
        my $Aterm = ( $A_amb / $A )**( 1 / ( $gamma - 1 ) );
        $ovd = $rho_amb * ( ( $Aterm - 1 ) + $Aterm * $series );
        $rho = $ovd + $rho_amb;
    }
    else {
        my $ln_rho =
          ( log( $sspd**2 / $gamma ) - $lnA ) / ( $gamma - 1 );
        $rho = exp($ln_rho);
        $ovd = $rho - $rho_amb;
    }

    # pabs = $rho * $sspd**2 / $gamma 
    #  = ($rho+$ovd)*($sspd_amb+$sspdx)**2 /  $gamma;
    # $ovp = $pabs - $p_amb;
    # Expand to get:
    my $ovp =
      ( $rho_amb * $sspdx * ( 2 * $sspd_amb + $sspdx ) + $ovd * $sspd**2 ) /
      $gamma;
    my $pabs = $p_amb+$ovp;
    return wantarray ? ( $ovp, $pabs, $ovd, $rho, $sspdx, $sspd ) : $ovp;
}

sub ovd_from_sigma {

    # For convenience
    my ( $self, $sigma, $lnArat ) = @_;
    return ($self->pdc_from_sigma($sigma, $lnArat))[2];
}

sub ovp_from_sigma {

    # For convenience
    my ( $self, $sigma, $lnArat ) = @_;
    return $self->pdc_from_sigma($sigma, $lnArat);
}

sub lnArat_from_shock_ovp {

    # Temporary routine for adding entropy 
    my ( $self, $ovp ) = @_;

    my $gamma    = $self->{_gamma};
    my $sspd_amb = $self->{_sspd_amb};
    my $rho_amb  = $self->{_rho_amb};
    my $p_amb    = $self->{_p_amb};
    my $A_ref    = $self->{_A_ref};

    my $pabs    = $ovp + $p_amb;
    my $top     = ( $gamma + 1 ) * $pabs + ( $gamma - 1 ) * $p_amb;
    my $bot     = ( $gamma - 1 ) * $pabs + ( $gamma + 1 ) * $p_amb;
    my $rhoshock;
    if   ( $bot > 0 ) { $rhoshock = $rho_amb * $top / $bot }
    else              { $rhoshock = $rho_amb }

    # $A = $pabs / $rhoshock**$gamma;
    # $Arat = $A/$A_ref
    my $lnArat = 0;
    if ( $pabs > 0 && $rhoshock > 0 && $A_ref > 0 ) {
        $lnArat = log($pabs) - $gamma * log($rhoshock) - log($A_ref);
    }
    return ($lnArat);
}

sub sigma_from_ovp {

    # Returns riemann integral from overpressure ratio and entropy variable
    # FIXME: use Taylor Series for small ovp
    # Could also return sspd and rho
    my ( $self, $ovp, $lnArat ) = @_;
    $lnArat=0 unless defined($lnArat);
    my $gamma    = $self->{_gamma};
    my $sspd_amb = $self->{_sspd_amb};
    my $p_amb = $self->{_p_amb};
    my $A_ref    = $self->{_A_ref};
    my $pabs     = $ovp + $p_amb;
    my $A        = $A_ref * exp($lnArat);
    my $rho = ( $pabs / $A )**( 1 / $gamma );
    my $sspd = sqrt( $gamma * $pabs / $rho );
    my $sspdx = ( $sspd - $sspd_amb ); 
    my $sigma = ( 2 / ( $gamma - 1 ) ) * $sspdx;
    return ($sigma);
}

sub sigma_from_ovprat {

    # Returns riemann integral from overpressure ratio and entropy variable
    my ( $self, $ovp, $log_Arat ) = @_;
    my $gamma    = $self->{_gamma};
    my $sspd_amb = $self->{_sspd_amb};
    if ( defined($log_Arat) && $log_Arat != 0 ) {
        my $Arat = exp($log_Arat);
        my $Aterm = $Arat**( 1 / ( $gamma - 1 ) );
        $ovp = $Aterm * ( 1 + $ovp ) - 1;
    }
    my $qq = ( $gamma - 1 ) / ( 2. * $gamma );
    my $dsspd;
    if ( $ovp < 1.e-4 ) {
        $dsspd = $sspd_amb * (
            $ovp * $qq * (
                1 + $ovp * ( ( $qq - 1 ) / 2 * ( 1 + $ovp * ( $qq - 2 ) / 3 ) )
            )
        );
    }
    else {
        $dsspd = $sspd_amb * ( ( 1 + $ovp )**$qq - 1 );
    }
    my $sigma = ( 2 / ( $gamma - 1 ) ) * $dsspd;
    return ($sigma);
}

sub dsigma_dovprat {

    # Returns d(sigma)/d(ovp) a given overpressure and entropy variable
    my ( $self, $ovp, $log_Arat ) = @_;
    my $gamma    = $self->{_gamma};
    my $sspd_amb = $self->{_sspd_amb};
    my $qq       = ( $gamma - 1 ) / ( 2. * $gamma );
    my $Arat_pow = 1;
    if ( defined($log_Arat) && $log_Arat != 0 ) {
        my $Arat = exp($log_Arat);
        $Arat_pow = $Arat**( 1 / ( 2. * $gamma ) );
    }
    my $dsigma_dovp =
      $sspd_amb / $gamma * $Arat_pow * ( 1 + $ovp )**( $qq - 1 );
    return ($dsigma_dovp);
}

sub ovprat_from_sigma {
    my ( $self, $sigma, $log_Arat ) = @_;

    # Returns overpressure from riemann integral
    my $gamma    = $self->{_gamma};
    my $sspd_amb = $self->{_sspd_amb};
    my $gpow     = ( $gamma - 1 ) / ( 2. * $gamma );
    my $p        = 2 * $gamma / ( $gamma - 1 );

    # avoid roundoff, use taylor series for small x
    my $ovprat;
    my $x = ( $gamma - 1 ) / ( 2 * $sspd_amb ) * $sigma;
    if ( $x < 1.e-4 ) {
        $ovprat = $x * $p *
          ( 1 + $x * ( ( $p - 1 ) / 2 * ( 1 + $x * ( $p - 2 ) / 3 ) ) );
    }
    else {
        $ovprat =
          ( ( $gamma - 1 ) / ( 2 * $sspd_amb ) * $sigma + 1 )**( 1 / $gpow ) -
          1;
    }
    if ( defined($log_Arat) ) {
        my $Arat = exp($log_Arat);
        my $Aterm = $Arat**( 1 / ( $gamma - 1 ) );
        $ovprat = ( 1 + $ovprat ) / $Aterm - 1;
    }
    return ($ovprat);
}

sub sspd_from_sigma {
    my ( $self, $sigma ) = @_;
    my $gamma    = $self->{_gamma};
    my $sspd_amb = $self->{_sspd_amb};
    my $sspdx    = 0.5 * ( $gamma - 1 ) * $sigma;
    my $sspd     = $sspd_amb + $sspdx;
    return ( $sspd, $sspdx );
}

sub odrat_from_sigma {

    # Returns overdensity ratio [=(rho-rho_amb)/rho_amb] from sigma
    # Density can then be computed from:
    #   my $density = $rho_amb * ( $odrat + 1 );
    my ( $self, $sigma ) = @_;
    my $gamma    = $self->{_gamma};
    my $sspd_amb = $self->{_sspd_amb};
    my $sspdx    = 0.5 * ( $gamma - 1 ) * $sigma;
    my $p      = 2 / ( $gamma - 1 );

    # Here is what we are approximating
    # my $od      = ( $sspd / $sspd_amb )**$pow - 1;

    my $x = $sspdx / $sspd_amb;
    my $odrat;

    # avoid roundoff, use taylor series for small x
    if ( $x < 1.e-4 ) {
        $odrat = $x * $p *
          ( 1 + $x * ( ( $p - 1 ) / 2 * ( 1 + $x * ( $p - 2 ) / 3 ) ) );
    }
    else {
        $odrat = ( 1 + $x )**$p - 1;
    }
    return ($odrat);
}

## TO BE DELETED
sub Alim_from_Hmax {
    my ( $self, $HH_max ) = @_;
    my $gamma    = $self->{_gamma};
    my $rho_amb  = $self->{_rho_amb};
    my $sspd_amb = $self->{_sspd_amb};

    # FIXME: I also found this equation: which is right????
    #my $Alim = sqrt( 4 * $HH_max / ( $sspd_amb * ( $gamma + 1 ) ) );

    # Miles writes I+=($k+1)/2 * A * A * rho * sspd * Rd**2 / R
    # where Rd=(E0/P0)^1/3 and k=gamma-1)/2 -> k+1=(gamma+1)/2

    # Check:
    # We also have I+ R = Hmax * rho ??

    # If so, then
    # I+ R =($gamma+1)/4 * A * A * rho * sspd * Rd**2 = Hmax rho
    # A**2 = ( 4* Hmax / (sspd * Rd**2 * (gamma+1))

    #my $Alim = sqrt( 4 * $HH_max * $gamma / ( $rho_amb * ( $gamma + 1 ) ) );
    # FIXME: I also found this equation: which is right????
    # This is for Rd=1
    my $Alim=1.1111111111;
    if ( $sspd_amb > 0 ) {
        $Alim = 4 * $HH_max / ( $sspd_amb * ( $gamma + 1 ) );
        if ( $Alim < 0 ) { $Alim = 0 }
        else { $Alim = sqrt( 4 * $HH_max / ( $sspd_amb * ( $gamma + 1 ) ) ); }
    }
    return ($Alim);
}

sub miles_AB {
    my ( $self, $dlnp_dlnr, $lambda, $ovp ) = @_;
    my $gamma = $self->{_gamma};
    my $A     = 0;
    my $B     = 0;
    my $term  = ( -$dlnp_dlnr - 1 );
    if ( $term > 0 && $lambda > 0 ) {
        $A = $lambda * $ovp / ( $gamma * sqrt( 2 * $term ) );
        $B = 0.5 / $term - log($lambda);
    }
    return ( $A, $B );
}

sub sadek_dlnp_dlnr_from_dpdr {
    my ( $medium, $lambda, $ovp, $slope ) = @_;
    return $medium->sadek_dlnp_dlnr_from_slope( $lambda, $ovp, $slope, 0,
        0 );
}
sub sadek_dlnp_dlnr_from_dudr {
    my ( $medium, $lambda, $ovp, $slope ) = @_;
    return $medium->sadek_dlnp_dlnr_from_slope( $lambda, $ovp, $slope, 1,
        0 );
}
sub sadek_dlnp_dlnr_from_dpdt {
    my ( $medium, $lambda, $ovp, $slope ) = @_;
    return $medium->sadek_dlnp_dlnr_from_slope( $lambda, $ovp, $slope, 0,
        1 );
}
sub sadek_dlnp_dlnr_from_dudt {
    my ( $medium, $lambda, $ovp, $slope ) = @_;
    return $medium->sadek_dlnp_dlnr_from_slope( $lambda, $ovp, $slope, 1,
        1 );
}

sub sadek_dlnp_dlnr_from_slope {

    my ( $medium, $lambda, $ovp, $slope, $ipu, $irt ) = @_;

    # Given
    #   lambda = the range
    #   ovp = the overpressure at the shock front
    #   slope = df/dz = the slope of a wave profile at the shock front
    #   ipu = 0 f = ovp  [DEFAULT]
    #       = 1 f = up
    #   irt = 0 z = r  ( slope is df/dr ) [DEFAULT]
    #       = 1 z = t  ( slope is df/dt )
    #   
    # Compute
    #   dlnp_dlnr_s = the slope of the shock decay curve

    # The _ovp in the name refers to the call arg which is overpressure
    #
    # Reference: Initial decay of flow properties ... sadek 1983
    # This works even for zero initial pressure and sound speed
    # Note that sadek's equations assume positive zero sound speed and pressure,
    # So I had to redo them.

    my $p_amb    = $medium->{_p_amb};
    my $symmetry     = $medium->{_symmetry};
    my $gamma    = $medium->{_gamma};
    my $sspd_amb = $medium->{_sspd_amb};
    my $rho_amb  = $medium->{_rho_amb};
    my $pabs     = $ovp + $p_amb;
    if ($lambda<=0) {croak "Negative lambda in shock; check args\n"; }
    if ( $pabs < 0 ) { print STDERR "Negative abs p in sadek\n"; return }

    # Shock front properties ...
    my ( $D, $up, $D_minus_c0, $dudp_sf, $rho ) =
      $medium->ushock_from_ovp($ovp);

    #    Alternate coding, for reference...
    #    my $kk  = ( $gamma + 1 ) / 2;
    #    my $Dsq = ( $sspd_amb**2 + $kk * $ovp / $rho_amb );
    #    if ( $Dsq < 0 ) { $Dsq = 0 }
    #    my $D         = sqrt($Dsq);
    #    my $up        = $ovp / ( $rho_amb * $D );
    #    my $dlnu_dlnp = 1 - 0.5 * $kk * $up / $D;
    #    # dudp = u/ovp*dlnu_dlnp =
    #    my $dudp_sf = $dlnu_dlnp / ( $rho_amb * $D );
    #
    #    # This gives identical result:
    #    # my $dudp_old = ( 1 + ( $sspd_amb / $D )**2 ) / ( 2 * $rho_amb * $D );
    #
    #    # Density
    #    my $topr = $rho_amb * $sspd_amb**2 + $kk * $ovp;
    #    my $botr = $topr - $ovp;
    #    my $rho  = $rho_amb * $topr / $botr;

    # Sound speed
    my $rhocsq  = $gamma * $pabs;
    my $sspd_sq = $rhocsq / $rho;
    if ( $sspd_sq < 0 ) { $sspd_sq = 0 }
    my $aa        = sqrt($sspd_sq);

    #######################################################
    # FIXME: need to expand aa-c0 for small ovp
    #######################################################
    # Attempted Roundoff control for det...
    # sspd_sq = $gamma*$pabs/$rho = $gamma*($p_amb+$ovp)/$rho 
    #    = $sspd_amb_sq+$gamma*$ovp/rho;
    # sspd_sq - sspd_amb_sq = gamma*$ovp/$rho
    # (sspd-sspd_amb)*(sspd+sspd_amb)=gamma*ovp/rho
    # ... so ..
    # sspdx == sspd-sspd_amb=gamma*ovp/rho/(sspd+sspd_amb)
    #my $aa_minus_c0 = $gamma * $ovp / $rho / ( $aa + $sspd_amb );
    my $aa_minus_c0 = $aa - $sspd_amb;
    #######################################################

    my $D_minus_u = $D - $up;
    #my $det       = $sspd_sq - $D_minus_u**2;
    my $det = ( $aa + $D_minus_u ) * ( $aa_minus_c0 - $D_minus_c0 + $up );

    # From notes ...; sadek eq 3.29, p24 of pdf

    # Solve
    #      det * dpdr|t + (D-u)D dpdr|s + rho*a**2 * D * dudr|s +
    #              N*rho*a**2 * up * (D-u)/r = 0
    # where
    #       dudr|s = dpdr|s * dudp_s
    # So
    #
    #      det * dpdr|t + [ (D-u)D + rho*a**2 * D *dudp_s ]*dpdr|s +
    #              N*rho*a**2 * up * (D-u)/r = 0
    # Or
    #      det * dpdr|t + C1 * dpdr|s + C2 = 0
    # or
    #      dpdr|s = [- det * dpdr|t + C2 ] / C1

    # And in general, with different values for sign, C1 and C2:
    #     
    #      dpdr|s = [ $sign * det * slope + C2 ] / C1
    #

    my ($sign, $C1, $C2);

    # Handle slope = dpdr
    if ( !$ipu && !$irt ) {
        $sign = -1;
        $C1   = $D_minus_u * $D + $rhocsq * $D * $dudp_sf;
        $C2   = -$symmetry * $rhocsq * $up * $D_minus_u / $lambda;
    }

    # Handle slope = du/dr
    elsif ( $ipu && !$irt ) {
        $sign = -1;
        $C1   = $D_minus_u * $D * $dudp_sf + $D / $rho;
        $C2   = -$symmetry * $sspd_sq * $up / $lambda;
    }

    # Handle slope = dp/dt
    elsif ( !$ipu && $irt ) {
        $sign = 1;
        $C1 =
          ( $sspd_sq + $up * $D_minus_u ) * $D + $rhocsq * $D * $D * $dudp_sf;
        $C2 = -$symmetry * $rhocsq * $up * $D * $D_minus_u / $lambda;
    }

    # Handle slope = du/dt
    elsif ( $ipu && $irt ) {
        $sign = 1;
        $C1 = ( $sspd_sq + $up * $D_minus_u ) * $D * $dudp_sf + $D * $D / $rho;
        $C2 = -$symmetry * $sspd_sq * $up * $D / $lambda;
    }

    my $dpdr_sf      = ($sign* $det * $slope + $C2 ) / $C1;
    my $dlnp_dlnr_sf = $dpdr_sf * $lambda / $ovp;

    return ($dlnp_dlnr_sf);
}

sub sadek_slopes_from_dlnp_dlnr {

    my ( $medium, $lambda, $ovp, $dlnp_dlnr_sf ) = @_;

    # Given
    #   lambda = the range
    #   ovp = the overpressure at the shock front
    #   dlnp_dlnr_sf
    #
    # Compute wave front slopes:
    #   dp/dr, du/dr, dp/dt, du/dt
    #
    my $p_amb    = $medium->{_p_amb};
    my $symmetry     = $medium->{_symmetry};
    my $gamma    = $medium->{_gamma};
    my $sspd_amb = $medium->{_sspd_amb};
    my $rho_amb  = $medium->{_rho_amb};

    my $pabs     = $ovp + $p_amb;
    if ( $pabs < 0 ) { print STDERR "Negative abs p in sadek\n"; return }

    # Shock front properties ...
    my ( $D, $up, $D_minus_c0, $dudp_sf, $rho ) =
      $medium->ushock_from_ovp($ovp);

    # Sound speed
    my $rhocsq  = $gamma * $pabs;
    my $sspd_sq = $rhocsq / $rho;
    if ( $sspd_sq < 0 ) { $sspd_sq = 0 }
    my $aa = sqrt($sspd_sq);

    # From notes ...; sadek eq 3.29, p24 of pdf

    # Solve
    #      det * dpdr|t + (D-u)D dpdr|s + rho*a**2 * D * dudr|s +
    #              N*rho*a**2 * up * (D-u)/r = 0
    # where
    #       dudr|s = dpdr|s * dudp_s
    # So
    #
    #      det * dpdr|t + [ (D-u)D + rho*a**2 * D *dudp_s ]*dpdr|s +
    #              N*rho*a**2 * up * (D-u)/r = 0
    # Or
    #      det * dpdr|t + C1 * dpdr|s + C2 = 0
    # or
    #      dpdr|s = [- det * dpdr|t + C2 ] / C1

    # And in general, with different values for sign, C1 and C2:
    #
    #      dpdr|s = [ $sign * det * slope + C2 ] / C1
    #
    my ( $dpdr_t, $dudr_t, $dpdt_r, $dudt_r ) = ( 0, 0, 0, 0 );
    ##my $det = $sspd_sq - $D_minus_u**2;

    ###########################################################
    # FIXME: expand for small ovp
    ###########################################################
    # Roundoff control for det...
    # sspd_sq = $gamma*$pabs/$rho = $gamma*($p_amb+$ovp)/$rho 
    #    = $sspd_amb_sq+$gamma*$ovp/rho;
    # sspd_sq - sspd_amb_sq = gamma*$ovp/$rho
    # (sspd-sspd_amb)*(sspd+sspd_amb)=gamma*ovp/rho
    # ... so ..
    # sspdx == sspd-sspd_amb=gamma*ovp/rho/(sspd+sspd_amb)
    ##my $aa_minus_c0 = $gamma * $ovp / $rho / ( $aa + $sspd_amb );
    my $aa_minus_c0 = $aa - $sspd_amb;
    ###########################################################

    my $D_minus_u = $D - $up;
    #my $det       = $sspd_sq - $D_minus_u**2;
    my $det = ( $aa + $D_minus_u ) * ( $aa_minus_c0 - $D_minus_c0 + $up );
    if ( $det != 0 ) {

        my $D_minus_u = $D - $up;
        my $dpdt_sf   = $D * $dlnp_dlnr_sf * $ovp / $lambda;
        my $dudt_sf   = $dudp_sf * $dpdt_sf;
        $dpdr_t =
          -( $D_minus_u * $dpdt_sf +
              $rhocsq * $dudt_sf +
              $symmetry * $rhocsq * $up * $D_minus_u / $lambda ) /
          $det;

        $dudr_t =
          -( $D_minus_u * $dudt_sf +
              $dpdt_sf / $rho +
              $symmetry * $sspd_sq * $up / $lambda ) /
          $det;

        $dpdt_r =
          ( ( $sspd_sq + $up * $D_minus_u ) * $dpdt_sf +
              $rhocsq * $D * $dudt_sf +
              $symmetry * $rhocsq * $up * $D * $D_minus_u / $lambda ) /
          $det;

        $dudt_r =
          ( ( $sspd_sq + $up * $D_minus_u ) * $dudt_sf +
              $D / $rho * $dpdt_sf +
              $symmetry * $sspd_sq * $up * $D / $lambda ) /
          $det;

    }

    return wantarray ? ( $dpdr_t, $dudr_t, $dpdt_r, $dudt_r ) : $dpdr_t;
}

sub shock_front_values {

    my ( $self, $X, $Y, $dYdX ) = @_;

    # Given:
    #   symmetry = (0,1,2) 
    #   gamma
    #   X = ln(range)
    #   Y = ln(ovp) = ln(overpressure ratio) 
    #   dYdX = d(ln P)/d(ln R) along shock front
    #
    # Compute some shock front values and return in a hash ref
    #
    # FIXME: some of these expressions would benefit from Taylor series
    # expansions for small ovp

    my $symmetry = $self->{_symmetry};
    my $gamma    = $self->{_gamma};
    my $p_amb    = $self->{_p_amb};
    my $sspd_amb = $self->{_sspd_amb};
    my $rho_amb  = $self->{_rho_amb};

    my $lambda = exp($X);
    my $ovprat = exp($Y);
    my $ovp    = $p_amb * $ovprat;
    my $pabs   = $ovp + $p_amb;

    if ( $pabs < 0 ) { print STDERR "Negative abs p in sadek\n"; return }

    # Compute some shock front properties ...

    # D = Shock speed D
    my $term     = $ovprat * ( $gamma + 1 ) / ( 2 * $gamma );
    my $mshock   = sqrt( 1 + $term );
    my $D   = $sspd_amb * $mshock;
    my $D_minus_c0  = $D - $sspd_amb;

    # Taylor series expansion of the square root for small overpressures
    # Switch at a point which keeps error < about 1.e-17
    if ( $term < 1.e-4 ) {
        $D_minus_c0 = $term * ( 0.5 + $term * ( -0.125 + $term / 16 ) );
    }

    # up = shock particle velocity
    my $up = $ovprat * $p_amb / ( $rho_amb * $D );

    # dup/dovp along shock front
    my $dudp_sf =
      $sspd_amb *
      $sspd_amb / $gamma / $D *
      ( 1 - ( $gamma + 1 ) / 4 * $up / $D );

    # density ratio and density
    my $top  = ( $gamma + 1 ) * $pabs + ( $gamma - 1 ) * $p_amb;
    my $bot  = ( $gamma - 1 ) * $pabs + ( $gamma + 1 ) * $p_amb;
    my $rhorat = ( $bot > 0 ) ? $top / $bot : 1;
    my $rho = $rho_amb * $rhorat;

    # Sound speed
    my $rhocsq  = $gamma * $pabs;
    my $sspd_sq = $rhocsq / $rho;
    if ( $sspd_sq < 0 ) { $sspd_sq = 0 }
    my $aa = sqrt($sspd_sq);

    # sigma and Sigma, the C+ characteristic variable
    my $sspdx = ( $aa - $sspd_amb ); 
    my $sigma = ( 2 / ( $gamma - 1 ) ) * $sspdx;
    my $Sigma = ( $symmetry == 0 ) ? $sigma : $sigma * $lambda**( $symmetry / 2 );

    # FIXME: duplicate coding; call previous routine here

    # Compute the slopes of the wave profile at the shock front ...
    #   dp/dr, du/dr, dp/dt, du/dt

    # Reference for slopes:
    #  INITIAL DECAY OF FLOW PROPERTJES OF PLANAR, CYLINDRJCAL AND SPHERICAL BLAST WAVES
    #  H. S. I. Sadek and J. J. Gott1ieb
    #  UTIAS Technica1 Note No. 244, CN ISSN 0082-5263
    #  October, 1983
    #  See eq 3.29, p24 of pdf

    # Solve
    #      det * dpdr|t + (D-u)D dpdr|s + rho*a**2 * D * dudr|s +
    #              N*rho*a**2 * up * (D-u)/r = 0
    # where
    #       dudr|s = dpdr|s * dudp_s
    # So
    #
    #      det * dpdr|t + [ (D-u)D + rho*a**2 * D *dudp_s ]*dpdr|s +
    #              N*rho*a**2 * up * (D-u)/r = 0
    # Or
    #      det * dpdr|t + C1 * dpdr|s + C2 = 0
    # or
    #      dpdr|s = [- det * dpdr|t + C2 ] / C1

    # And in general, with different values for sign, C1 and C2:
    #
    #      dpdr|s = [ $sign * det * slope + C2 ] / C1
    #
    my ( $dpdr_t, $dudr_t, $dpdt_r, $dudt_r ) = ( 0, 0, 0, 0 );

    ###########################################################
    # FIXME: expand for small ovp
    ###########################################################
    # Roundoff control for det...
    # sspd_sq = $gamma*$pabs/$rho = $gamma*($p_amb+$ovp)/$rho 
    #    = $sspd_amb_sq+$gamma*$ovp/rho;
    # sspd_sq - sspd_amb_sq = gamma*$ovp/$rho
    # (sspd-sspd_amb)*(sspd+sspd_amb)=gamma*ovp/rho
    # ... so ..
    # sspdx == sspd-sspd_amb=gamma*ovp/rho/(sspd+sspd_amb)
    ##my $aa_minus_c0 = $gamma * $ovp / $rho / ( $aa + $sspd_amb );
    my $aa_minus_c0 = $aa - $sspd_amb;
    ###########################################################

    my $D_minus_u = $D - $up;
    #  $det = $sspd_sq - $D_minus_u**2;
    my $det = ( $aa + $D_minus_u ) * ( $aa_minus_c0 - $D_minus_c0 + $up );
    if ( $det != 0 ) {

        my $D_minus_u = $D - $up;
        my $dpdt_sf   = $D * $dYdX * $ovp / $lambda;
        my $dudt_sf   = $dudp_sf * $dpdt_sf;
        $dpdr_t =
          -( $D_minus_u * $dpdt_sf +
              $rhocsq * $dudt_sf +
              $symmetry * $rhocsq * $up * $D_minus_u / $lambda ) /
          $det;

        $dudr_t =
          -( $D_minus_u * $dudt_sf +
              $dpdt_sf / $rho +
              $symmetry * $sspd_sq * $up / $lambda ) /
          $det;

        $dpdt_r =
          ( ( $sspd_sq + $up * $D_minus_u ) * $dpdt_sf +
              $rhocsq * $D * $dudt_sf +
              $symmetry * $rhocsq * $up * $D * $D_minus_u / $lambda ) /
          $det;

        $dudt_r =
          ( ( $sspd_sq + $up * $D_minus_u ) * $dudt_sf +
              $D / $rho * $dpdt_sf +
              $symmetry * $sspd_sq * $up * $D / $lambda ) /
          $det;

    }

    my $rhash = {
        dpdr_t        => $dpdr_t,
        dudr_t        => $dudr_t,
        dpdt_r        => $dpdt_r,
        dudt_r        => $dudt_r,
        up            => $up,
        density_ratio => $rhorat,
        ushock        => $D,
        sspd          => $aa,
	sigma         => $sigma,
	Sigma         => $Sigma,
    };

    return $rhash;
}

1;
