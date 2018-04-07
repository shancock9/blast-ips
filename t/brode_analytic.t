use strict;
use warnings;
use Test;
use Blast::IPS;

BEGIN {

    # This compares the tables to two analytical models for a spherical point
    # source with gamma=1.4 given by H.Brode.  See the notes in the subroutine
    # for more details.

=pod
References:

RAND PROJECT AIR FORCE SANTA MONICA CA, and Brode, H L. 1954. Numerical Solutions of Spherical Blast Waves. http://oai.dtic.mil/oai/oai?&verb=getRecord&metadataPrefix=html&identifier=ADA595878. 

Brode, Harold L. 1955. "Numerical Solutions of Spherical Blast Waves". Journal of Applied Physics. 26 (6): 766-775. 
=cut

    plan tests => 1;
}

# The actual maximum error is about 8.5% near lambda=0.45
my $TOL     = 0.09;
my $VERBOSE = 0;

# Create a table for this case
my $symmetry    = 2;
my $gamma       = 1.4;
my $blast_table =
          Blast::IPS->new( 'symmetry' => $symmetry, 'gamma' => $gamma );
if ( !defined($blast_table) ) {
    die "missing table for sym=$symmetry, gamma=$gamma\n";
}

# Loop to calculate the maximum absolute relative error in overpressure
my $err_max;
my $rbounds = $blast_table->get_table_bounds();
my $Xmin= $rbounds->[0]->[0];
my $Xmax= log(2.7999);  # Brode fit stops at lambda=2.8
my $rtable = $blast_table->table_gen( 50, $Xmin, $Xmax );
if ($VERBOSE) {print "X\tY\tdYdX\tY_b\tdYdX_b\terr\n"}
foreach my $point ( @{$rtable} ) {
    my ( $X, $Y,$dYdX,$Z,$dZdX ) = @{$point};
    my ($X_b, $Y_b, $dYdX_b)= brode_ideal_point_source_fit($X); 
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


