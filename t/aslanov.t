use strict;
use warnings;
use Test;
use Blast::IPS;

=pod

S.K Aslanov developed an approximate analytical model for a spherical point explosion 
in an ideal gas. His paper gives coefficients for the case of gamma=1.4, so this
is what is used here.

He shows a graph of the error out to a scaled range of about 20.  The error in
the model has a local peak of about 14% near a scaled range of about 0.8. 
A comparison here gives close to the same error profile over that range.

However if we look over a wider range, the error increases beyond 20 and
rapidly grows to about 16% at long distances. This is unexpected because the
model was presented as asymptotically accurate at long range.  

Reference:

Aslanov, S. K. 2006. "On point blast theory". Fluid Dynamics. 41 (1): 147-151. 

=cut


BEGIN {

    plan tests => 1;
}

my $VERBOSE = 1;

my $symmetry = 2;
my $gamma    = 1.4;

# Create a table for this case
my %args = ( 'symmetry' => $symmetry, 'gamma' => $gamma );
my $blast_table = Blast::IPS->new( \%args );
if ( !defined($blast_table) ) {
    die "missing table for sym=$symmetry, gamma=$gamma\n";
}
my $table_name = $blast_table->get_table_name();

my $rbounds = $blast_table->get_table_bounds();
my $Xmin    = $rbounds->[0]->[0];
my $Xmax    = $rbounds->[1]->[0];
if ( $Xmax > 20 ) { $Xmax = 20 }

my $TOL = 0.2;

# Generate a wide range of table of values for testing
my $rtable = $blast_table->table_gen( 51, $Xmin, $Xmax );

my $err_max;
if ($VERBOSE) {
        print "lambda\tX\tY\tY_k\terr\n";
}
foreach my $item ( @{$rtable} ) {
    my ( $X,   $Y,   $dYdX )   = @{$item};
    my ( $X_k, $Y_k, $dYdX_k ) = aslanov_model($X);
    my $err = abs( $Y - $Y_k );
    if ( !defined($err_max) || $err > $err_max ) { $err_max = $err }

    if ($VERBOSE) {
	my $lambda=exp($X);
        print "$lambda\t$X\t$Y\t$Y_k\t$err\n";
    }
}
my $err_max_pr = sprintf "%0.3g", $err_max;

if ( !ok( $err_max <= $TOL ) ) {
    my $text = "# Max error is $err_max_pr\n";
    print "$text";
}

sub aslanov_model {

    my ($X)    = @_;
    my $lambda = exp($X);
    my $a1     = 0.264;
    my $a2     = 0.0167;
    my $a3     = -0.00737;
    my $a4     = -0.0923;
    my $c0     = 1.36;
    my $c1     = $c0;
    my $c2     = $c0;
    my $Rstar0 = 3.06;

    #####################################################################
    # Here are equations for $c0 and $Rstar0 for any gamma from Aslanov 2006
    # my $alpha = ( 2 * $gamma + 1 ) / ( 8 * $gamma * ( $gamma - 1 ) );
    # my $c0 = 2 * $gamma / ( $alpha * ( $gamma + 1 ) );
    # my $Rstar0 = ( $gamma + 1 ) * ( 2 * $gamma - 1 ) / $gamma;
    #####################################################################

    # define a sub to calculate overpressure ratio, given range in meters
    my $ovp_calc = sub {
        my ($R) = @_;

        my $Rrat    = $R / $Rstar0;
        my $term1   = $a1 / ( $R * sqrt( log( $Rrat + $c0 ) ) );
        my $term2   = $a2 / ( $R * ( log( $Rrat + 1.0 ) )**2 );
        my $term3   = $a3 / ( $R * ( log( $Rrat + $c1 ) )**3.5 );
        my $term4   = $a4 * sqrt( log( $Rrat + $c2 ) ) / $R**2;
        my $ovp_atm = $term1 + $term2 + $term3 + $term4;
        return $ovp_atm;
    };

    my $ovp_atm = $ovp_calc->($lambda);

    # differentiate to get d(ln P)/d(ln R)
    my $delta_r = 0.001 * $lambda;
    my ($ovpm)  = $ovp_calc->( $lambda - $delta_r );
    my ($ovpp)  = $ovp_calc->( $lambda + $delta_r );
    my $dYdX    = 0.5 * $lambda * ( $ovpp - $ovpm ) / ( $delta_r * $ovp_atm );
    my $Y       = log($ovp_atm);
    return ( $X, $Y, $dYdX );
}

