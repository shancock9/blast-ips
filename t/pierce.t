use strict;
use warnings;
use Test;
use Blast::IPS;

=pod

Allan Pierce used the theory of N waves to derive equations for
overpressure ratio and positive phase duration at a long distance
from a point explosion. He evaluated the constants in the model by
fitting it to the values obtained by Brode at a scaled range of 3.

Reference:
Pierce, Allan D. Acoustics: An Introduction to Its Physical Principles and Applications. New York: Acoustical Society of America, 1995. 

=cut


BEGIN {

    plan tests => 1;
}

my $VERBOSE = 0;

my $symmetry = 2;
my $gamma    = 1.4;

# Create a table for this case
my $blast_table = Blast::IPS->new( 'symmetry' => $symmetry, 'gamma' => $gamma );

if ( !defined($blast_table) ) {
    die "missing table for sym=$symmetry, gamma=$gamma\n";
}
my $table_name = $blast_table->get_table_name();

my $rbounds = $blast_table->get_table_bounds();
my $Xmin    = $rbounds->[0]->[0];
my $Xmax    = $rbounds->[1]->[0];
if ( $Xmin < 1 ) {$Xmin=1}
if ( $Xmax > 20 ) { $Xmax = 20 }

# Brode's model has an error on the order of 3%.
# Pierce's model has an error on the order of 8% at long range.
# We will allow a maximum error of 10% in testing here.
my $TOL = 0.1;

# Generate a wide range of table of values for testing
my $rtable = $blast_table->table_gen( 51, $Xmin, $Xmax );

my $err_max;
my $err2_max;
if ($VERBOSE) {
        #print "lambda\tX\tY\tY_k\terr\tTpos\n";
        print "lambda\tX\tY\tY_k\terr\tTpos\tTpos_k\terr2\n";
}
foreach my $item ( @{$rtable} ) {
    my ( $X,   $Y,   $dYdX, $Z )   = @{$item};
    my ( $X_k, $Y_k, $dYdX_k, $Tpos_k ) = pierce($X);
    my $err = abs( $Y - $Y_k );
    if ( !defined($err_max) || $err > $err_max ) { $err_max = $err }


    ## OLD method, has been deleted:
    ##my $rs=exp($X);
    ##my $zs=exp($Z);
    #my ($Tpos, $Lpos)= $blast_table->get_positive_phase( $rs, $zs ); 

    my $ret = $blast_table->wavefront( 'X' => $X );
#    my $z_pose_rs = $ret->{z_pose_rs};
#    my $Tpos      = $zs - $z_pose_rs;
    my $Tpos = $ret->{Tpos};

    my $err2 = abs($Tpos-$Tpos_k)/$Tpos;
    if ( !defined($err2_max) || $err2 > $err2_max ) { $err2_max = $err2 }
    #print STDERR "$rs, $zs, $Tpos, $Tpos_k, $err2\n";

    if ($VERBOSE) {
	my $lambda=exp($X);
        print "$lambda\t$X\t$Y\t$Y_k\t$err\t$Tpos\t$Tpos_k\t$err2\n";
    }
}
my $err_max_pr = sprintf "%0.3g", $err_max;
my $err2_max_pr = sprintf "%0.3g", $err2_max;

if ( !ok( $err_max <= $TOL && $err2_max <= $TOL) ) {
    my $text = "# Max ovp error is $err_max_pr; max Tpos err is $err2_max_pr\n";
    print "$text";
}

sub pierce {

    my ($X)    = @_;
    my $lambda = exp($X);
    my $gamma=1.4;

    # define a sub to calculate overpressure ratio, given range in meters
    # Note that Pierce normalized to a distance (E/(rho*C**2)^(1/3) which
    # differs from the normalization used here by gamma**-(1/3).
    my $ovp_calc = sub {
        my ($r)     = @_;
        my $rstar   = 0.7 / $gamma**( 1 / 3 );
        my $rrat    = $r / $rstar;
        my $ovp_atm = 0;
        my $Tpos = 0;
        if ( $rrat > 1 ) {
            $ovp_atm = $gamma * 0.4 / $rrat / sqrt( log($rrat) );
            $Tpos    = 0.32 * sqrt( log($rrat) ) / $gamma**( 1 / 3 );
        }
        return ($ovp_atm, $Tpos);
    };

    my ($ovp_atm, $Tpos) = $ovp_calc->($lambda);

    # differentiate to get d(ln P)/d(ln R)
    my $delta_r = 0.001 * $lambda;
    my ($ovpm, $Tposm)  = $ovp_calc->( $lambda - $delta_r );
    my ($ovpp, $Tposp)  = $ovp_calc->( $lambda + $delta_r );
    my $dYdX    = 0.5 * $lambda * ( $ovpp - $ovpm ) / ( $delta_r * $ovp_atm );
    my $Y       = log($ovp_atm);
    return ( $X, $Y, $dYdX, $Tpos );
}

