use strict;
use warnings;
use Test;
use Blast::IPS;

BEGIN {

=pod

    This compares an approximation of Sakurai for a spherical blast with gamma=1.4.
    Although this is an improvement over the simple similarity theory, the error
    is still large.
    

    Reference: Sachdev, p71
    Sachdev, P. L. 2004. Shock waves and explosions. Boca Raton: Chapman & Hall/CRC. 
=cut

    plan tests => 1;
}

# Allow up to 9% error down to lambda=0.4 Actual error is 8.4% there.
my $TOL     = 0.09; 
my $VERBOSE = 1;

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
my $Xmax= log(0.4);  # max scaled range=0.4
my $rtable = $blast_table->table_gen( 50, $Xmin, $Xmax );
if ($VERBOSE) {print "X\tY\tdYdX\tY_b\tdYdX_b\terr\n"}
foreach my $point ( @{$rtable} ) {
    my ( $X, $Y,$dYdX,$Z,$dZdX ) = @{$point};
    my ($X_b, $Y_b, $dYdX_b)= sakurai_approximation($X); 
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

sub sakurai_approximation {

    # This is an approximation by Sakurai for the overpressure in atmospheres
    # for a point source in an ideal gas with constant gamma=1.4
    # Error exceeds 1% by lambda=0.3 and rises rapidly with range beyond that.
    # The velocity approximation is much more accurate.

    my ( $X)=@_;
    my $lambda=exp($X);

    my $ovp_calc = sub {
        my ($lambda) = @_;

        #my $rsc = $range_m / $r0d;
        my ($ovp_atm);

        # NOTE:
        # 0.156= 1.96/(2*pi), beause sakurai used different normalization
        # This exceeds 1% error by lambda=0.3
        # There is a second approximation that could be considered
        $ovp_atm = 1.07 + 0.157 / $lambda**3;    ## 0.157?
        return $ovp_atm;
    };

    my $ovp_atm = $ovp_calc->($lambda);

    # differentiate to get d(ln P)/d(ln R)
    my $delta_r   = 0.001 * $lambda;
    my ($ovpm)    = $ovp_calc->( $lambda - $delta_r );
    my ($ovpp)    = $ovp_calc->( $lambda + $delta_r );
    my $dYdX = 0.5 * $lambda * ( $ovpp - $ovpm ) / ( $delta_r * $ovp_atm );
    my $Y=log($ovp_atm);
    return ( $X, $Y, $dYdX);
}

