use strict;
use warnings;
use Test;
use Blast::IPS;

BEGIN {

    # This compares the tables to Brodes suggestion that the shock particle
    # velocity varies like the similarity solution down to relatively low overpressure.
    # See the notes in the subroutine for more details.
    # Here we verify that the maximum error is below 4% to a scaled range of 2 
    # where the overpressure ratio is 0.167. The error rises rapidly beyond there.

    # Brode merely noted this correlation but did not suggest this as a shock
    # prediction method. Others, such as Korobeinikov, did this.

=pod
References:

RAND PROJECT AIR FORCE SANTA MONICA CA, and Brode, H L. 1954. Numerical Solutions of Spherical Blast Waves. http://oai.dtic.mil/oai/oai?&verb=getRecord&metadataPrefix=html&identifier=ADA595878. 

Brode, Harold L. 1955. "Numerical Solutions of Spherical Blast Waves". Journal of Applied Physics. 26 (6): 766-775. 
=cut

    plan tests => 1;
}

# max error <4% down to scaled range of 2. Actual max is about 3.6%.
my $TOL     = 0.04; 
my $VERBOSE = 0;

# Create a table for this case
my $symmetry    = 2;
my $gamma       = 1.4;
my %args        = ( 'symmetry' => $symmetry, 'gamma' => $gamma );
my $blast_table = Blast::IPS->new( \%args );
if ( !defined($blast_table) ) {
    die "missing table for sym=$symmetry, gamma=$gamma\n";
}

# Loop to calculate the maximum absolute relative error in overpressure
my $err_max;
my $rbounds = $blast_table->get_table_bounds();
my $Xmin= $rbounds->[0]->[0];
my $Xmax= log(2);  # max scaled range=2
my $rtable = $blast_table->table_gen( 50, $Xmin, $Xmax );
if ($VERBOSE) {print "X\tY\tdYdX\tY_b\tdYdX_b\terr\n"}
foreach my $point ( @{$rtable} ) {
    my ( $X, $Y,$dYdX,$Z,$dZdX ) = @{$point};
    my ($X_b, $Y_b, $dYdX_b)= brode_velocity_fit($X); 
    my $err    = abs( $Y_b-$Y);
    if ($VERBOSE) {print "$X\t$Y\t$dYdX\t$Y_b\t$dYdX_b\t$err\n"}

    my $lambda = exp($X);
    my $ovprat = exp($Y);

    if ( !defined($err_max) || $err > $err_max ) { $err_max = $err }
}

if ($VERBOSE) {
    my $err_max_pr = sprintf "%0.3g", $err_max;
    print "max error : $err_max_pr\n";
}
ok( $err_max <= $TOL );


sub brode_velocity_fit {

    # Brode uses us=0.30*lambda**(-3/2)
    # Says it is adequate down to 1/10 atmosphere

    # This uses Brode's velocity fit idea, and equals the first part of
    # the Korobeinikov model.  Error below about 2% down to lambda=1.7,
    # Then rises rapidly to 4% at 2 and 10% at lambda=3.

    my ( $X)=@_;
    my $lambda=exp($X);

    # define a sub to calculate overpressure ratio, given range in meters
    my $ovp_calc = sub {
        my ($lambda) = @_;

        # From Barenblatt, pp82-83, we have
        # D=2/5*xi*(E t**-3 / rho0)**(1/5)
        # R=xi*(E*t**2/$rho0)**(1/5
        # P=2/(gamma+1)*rho0*D**2
        # PR**3 = (2/(gamma+1)*xi**5 *(2/5)**2 * E
        # or P/P0 = const / lambda**3
        # where
        # const = (2/(gamma+1))* xi0**5 * (2/5)**2
        #   = 0.1568 for xi0=1.033 (gamma=1.4)
        # Note that this gives the pressure ratio, for overpressure we have
        #    ovp = P/P0 -1 = const / lambda**3 -1
        # The solution is only valid when P/P0>>1 anyway
        # In the following I must drop the 1 from the cubic term
        my $ovp_atm; 

        # Brode suggested 0.3, but 0.3055 works better
        # I get 0.305372 (for alpha=0.851079) see p 182 korob.
        #my $voverc = 0.3055 / $lambda**1.5;
        my $voverc = 0.305372 / $lambda**1.5;

        #my $voverc = 0.3/ $lambda**1.5;

        my $term1 = 0.25 * ( $gamma + 1 ) * $voverc;
        my $soverc = $term1 + sqrt( $term1**2 + 1 );

        $ovp_atm = $gamma * $soverc * $voverc;
        my $yy = 0.5 * $gamma * ( $gamma + 1 ) * $ovp_atm;
        my $dlnp_dlnr = -3 / ( 2 - ( $yy / ( $yy + $gamma * $gamma ) ) );

        # d2(lnp)/d2(lnr)
        my $kk = 2 * $gamma / ( $gamma + 1 );
        my $d2lnp_d2lnr =
          ( -$dlnp_dlnr )**3 * $kk /
          ( 3 * $ovp_atm * ( 1 + $kk / $ovp_atm )**2 );
        return ( $ovp_atm, $dlnp_dlnr, $d2lnp_d2lnr );
    };

    my ( $ovp_atm, $dYdX, $d2lnp_d2lnr ) = $ovp_calc->($lambda);
    my $Y=log($ovp_atm);
    return ( $X, $Y, $dYdX);
}


sub brode_ideal_point_source_fit {

    # This combines two analytical models for a spherical point source
    # with gamma=1.4 given by H.Brode. 

    # The high pressure part is a modified similarity model
    #
    #       $ovp_atm = 1 + 0.1567 / $lambda**3;
    #
    # The low pressure part is a cubic polynomial in 1/lambda which is suggested
    # to be used in the range 0.26 <= lambda <=2.8
    #
    # The high pressure fit is better than the second down to about lambda=0.4
    # instead of lambda =0.26 as suggested by Brode. 
    # The high pressure fit has an error less than about 4% down to lambda=.35
    # or .4.  The lower pressure fit has about an 8% error in its suggested
    # range. Note that there is a discontinuity at the junction of these two
    # models (they were not actually intended to be joined as is done here).

    # Despite the poor accuracy, these are of historical interest because they
    # have been widely used in the literature. 

    # Reference: JAP, June 1955, Numerical solutions of spherical blast waves

    my ( $X ) = @_;
    my $lambda=exp($X);

    # define a sub to calculate overpressure ratio, given scaled range
    my $ovp_calc = sub {
        my ($lambda) = @_;

        my ($ovp_atm);


        ## if ( $rsc < 0.26 ) {
        if ( $lambda < 0.4 ) {
            $ovp_atm = 1 + 0.1567 / $lambda**3;
        }
        else {
            # For .26 < lambda <2.8, or
            # 0.1 < ovp_atm < 10
            $ovp_atm =
              0.137 / $lambda**3 + 0.119 / $lambda**2 + 0.269 / $lambda - 0.019;
        }
        return $ovp_atm;
    };

    #################################
    # Check if out of bounds
    if ( $lambda > 2.8 ) { return ( $lambda, 0, 0 ) }
    #################################

    my $ovp_atm = $ovp_calc->($lambda);

    # differentiate to get d(ln P)/d(ln R)
    my $delta_r   = 0.001 * $lambda;
    my ($ovpm)    = $ovp_calc->( $lambda - $delta_r );
    my ($ovpp)    = $ovp_calc->( $lambda + $delta_r );
    my $dYdX = 0.5 * $lambda * ( $ovpp - $ovpm ) / ( $delta_r * $ovp_atm );
    my $Y=log($ovp_atm);
    return ( $X, $Y, $dYdX );
}


