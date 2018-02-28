use strict;
use Test;

my $rgoldstine;

BEGIN {

    # This test compares the table of values for this problem published in 1955 by
    # Herman Goldstine and John Von Neumann.  

=pod
Reference:

Goldstine, Herman H., and John Von Neumann. 1955. "Blast wave calculation". Communications on Pure and Applied Mathematics. 8 (2): 327-353. 
DOI: 10.1002/cpa.3160080207  
=cut

    # 
    #   [lambda, ovp_atm]
    $rgoldstine = [
        [ 0.116259524252,   100.006454 ],
        [ 0.1221334324475,  86.361648 ],
        [ 0.12789606950381, 75.346492 ],
        [ 0.13396092733693, 65.70844 ],
        [ 0.1403498420368,  57.274879 ],
        [ 0.14708674038283, 49.895875 ],
        [ 0.15419880133804, 43.440928 ],
        [ 0.16171668834789, 37.794133 ],
        [ 0.1696747816391,  32.855619 ],
        [ 0.17792789444341, 28.622687 ],
        [ 0.18668277271883, 24.911026 ],
        [ 0.19598773462182, 21.657636 ],
        [ 0.20589806727374, 18.806757 ],
        [ 0.21647788515146, 16.309079 ],
        [ 0.22780082698383, 14.122594 ],
        [ 0.2399530756365,  12.209147 ],
        [ 0.25303405500836, 10.536768 ],
        [ 0.26715991451403, 9.076439 ],
        [ 0.28246631666976, 7.801273 ],
        [ 0.29911168927709, 6.690302 ],
        [ 0.31728094220405, 5.723753 ],
        [ 0.33719173945363, 4.881288 ],
        [ 0.36059050362386, 4.141609 ],
        [ 0.38601212599888, 3.498135 ],
        [ 0.41555078052407, 2.929444 ],
        [ 0.45011359412749, 2.42999 ],
        [ 0.48368054265007, 2.063264 ],
        [ 0.52114685931983, 1.749343 ],
        [ 0.56404362572013, 1.476218 ],
        [ 0.61223796692065, 1.244782 ],
        [ 0.66646162423341, 1.049211 ],
        [ 0.72894048503996, 0.88737472 ],
        [ 0.80113245679242, 0.737472 ],
        [ 0.88475492473518, 0.615528 ],
        [ 0.95536424456421, 0.537405 ],
        [ 1.03311210741498, 0.469579 ],
        [ 1.11870586322245, 0.410702 ],
        [ 1.21291558260546, 0.35958 ],
        [ 1.31952425199777, 0.314056 ],
        [ 1.43710509199034, 0.274633 ],
        [ 1.5667501393793,  0.240462 ],
        [ 1.71360899461067, 0.210088 ],
        [ 1.87586043486341, 0.183801 ],
        [ 2.05990452518119, 0.160478 ],
        [ 2.26365359598588, 0.140311 ],
        [ 2.48916604720312, 0.122848 ],
        [ 2.7452838691693,  0.107352 ],
        [ 3.02936675339156, 0.093939 ],
        [ 3.34441576844453, 0.082312 ],
        [ 3.70273601561048, 0.071987 ],
        [ 4.10100817691879, 0.063039 ],
        [ 4.543637102769,   0.055274 ],
        [ 5.04787074893142, 0.048374 ],
        [ 5.60956885337298, 0.042387 ],
        [ 6.23523369262219, 0.037184 ],
        [ 6.94922249581862, 0.032559 ],
        [ 7.74642190113362, 0.028544 ],
        [ 8.63648671250697, 0.025051 ],
        [ 9.63019420182122, 0.022008 ],
        [ 10.7662351328749, 0.019314 ],
        [ 12.0376602861922, 0.016985 ],
    ];

    plan tests => 1;
}

# The paper includes an error estimate which indicates an error rising to about
# 7.9e-3 at the lowest overpressure.  The actual error is much greater, being
# below about 2 percent to a scaled range of about 7  and increases to about 4
# percent at the end of the table.  Here we will compare with their values
# below a scaled range of 6 and require that the difference be less than 2
# percent.

my $TOL        = 0.02;
my $lambda_stop = 6;
print "Comparison with Goldstine-Von Neuman Tables with tolerance $TOL\n";

# Create a table for this case
use Blast::IPS;
my $symmetry    = 2;
my $gamma       = 1.4;
my %args        = ( 'symmetry' => $symmetry, 'gamma' => $gamma );
my $blast_table = Blast::IPS->new( \%args );
if ( !defined($blast_table) ) {
    die "missing table for sym=$symmetry, gamma=$gamma\n";
}

# Loop to calculate the maximum absolute relative error in overpressure
my $iQ = 'X';
my $err_max;
my $lambda_max = 0;
foreach my $point ( @{$rgoldstine} ) {
    my ( $lambda_t, $ovprat_t ) = @{$point};
    last if ($lambda_t > $lambda_stop);

    # Lookup the overpressure at the given scaled range
    my $Q = log($lambda_t);
    my $ret = $blast_table->lookup( $Q, $iQ );
    my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};
    my $ovprat = exp($Y);
    my $lambda = exp($X);
    my $err    = abs( $ovprat - $ovprat_t ) / $ovprat_t;
    #print STDERR "lambda=$lambda, err=$err\n";
    if ( !defined($err_max) || $err > $err_max ) { $err_max = $err }
    if ( $lambda_t > $lambda_max ) { $lambda_max = $lambda_t }
}

     my $err_max_pr    = sprintf "%0.3g", $err_max;
     my $lambda_max_pr = sprintf "%0.5g", $lambda_max;
#    my $lambda_pr = sprintf "%0.7g", $lambda;
print "Goldstine-Von Neumann Table max error to lambda=$lambda_max_pr: $err_max_pr\n";
#print STDERR ""( $err_max <= $TOL )";
ok( $err_max <= $TOL );

