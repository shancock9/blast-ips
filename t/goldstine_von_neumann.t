use strict;
use warnings;
use Test;
use Blast::IPS;

my $rgoldstine;

BEGIN {

  # This test compares the table of values for this problem published in 1955 by
  # Herman Goldstine and John Von Neumann.

=pod
Reference:

Goldstine, Herman H., and John Von Neumann. 1955. "Blast wave calculation". Communications on Pure and Applied Mathematics. 8 (2): 327-353. 
DOI: 10.1002/cpa.3160080207  

The radii in this calculation must be scaled to the scaling method used here.
There are different ways to do this.  A very accurate way is to scale the radii
by the factor which makes the first overpressure nearly exact according to the
high resolution calculations.  By this method the scaled range to the first
overpressure should be lambda=0.1165525.  So we will scale the radii in their
calculation by 0.1165525/0.500474. 

=cut

    #   [range, ovp_atm]
    $rgoldstine = [
        [ 0.500474,  100.006454 ],
        [ 0.52576,   86.361648 ],
        [ 0.550567,  75.346492 ],
        [ 0.576675,  65.70844 ],
        [ 0.604178,  57.274879 ],
        [ 0.633179,  49.895875 ],
        [ 0.663795,  43.440928 ],
        [ 0.696158,  37.794133 ],
        [ 0.730416,  32.855619 ],
        [ 0.765944,  28.622687 ],
        [ 0.803632,  24.911026 ],
        [ 0.843688,  21.657636 ],
        [ 0.88635,   18.806757 ],
        [ 0.931894,  16.309079 ],
        [ 0.980637,  14.122594 ],
        [ 1.03295,   12.209147 ],
        [ 1.089261,  10.536768 ],
        [ 1.15007,   9.076439 ],
        [ 1.215961,  7.801273 ],
        [ 1.287616,  6.690302 ],
        [ 1.365831,  5.723753 ],
        [ 1.451543,  4.881288 ],
        [ 1.55227,   4.141609 ],
        [ 1.661705,  3.498135 ],
        [ 1.788863,  2.929444 ],
        [ 1.937649,  2.42999 ],
        [ 2.082148,  2.063264 ],
        [ 2.243433,  1.749343 ],
        [ 2.428095,  1.476218 ],
        [ 2.635562,  1.244782 ],
        [ 2.868984,  1.049211 ],
        [ 3.137943,  0.88737472 ],
        [ 3.448715,  0.737472 ],
        [ 3.808693,  0.615528 ],
        [ 4.112652,  0.537405 ],
        [ 4.447341,  0.469579 ],
        [ 4.815805,  0.410702 ],
        [ 5.221359,  0.35958 ],
        [ 5.680288,  0.314056 ],
        [ 6.18645,   0.274633 ],
        [ 6.744546,  0.240462 ],
        [ 7.376744,  0.210088 ],
        [ 8.075204,  0.183801 ],
        [ 8.867477,  0.160478 ],
        [ 9.744576,  0.140311 ],
        [ 10.715362, 0.122848 ],
        [ 11.817898, 0.107352 ],
        [ 13.040818, 0.093939 ],
        [ 14.397041, 0.082312 ],
        [ 15.939538, 0.071987 ],
        [ 17.65402,  0.063039 ],
        [ 19.559449, 0.055274 ],
        [ 21.730074, 0.048374 ],
        [ 24.148072, 0.042387 ],
        [ 26.841434, 0.037184 ],
        [ 29.915013, 0.032559 ],
        [ 33.346797, 0.028544 ],
        [ 37.178348, 0.025051 ],
        [ 41.45606,  0.022008 ],
        [ 46.346489, 0.019314 ],
        [ 51.81972,  0.016985 ],
    ];

    # Convert to scaled range (see note above)
    foreach my $item ( @{$rgoldstine} ) {
        $item->[0] *= 0.1165525 / 0.500474;
    }

    plan tests => 1;
}

# The paper includes an error estimate which indicates an error rising to about
# 7.9e-3 at the lowest overpressure.  The actual error is somewhat higher,
# being below 1 percent to a scaled range of about 3 (overpressure ratio about
# 0.1) and increasing to about 4 percent at the end of the table.

my $TOL     = 0.04;
my $VERBOSE = 0;
my $lambda_stop;    # Use to limit range of comparison;

# Get a built-in table for this case
my $symmetry    = 2;
my $gamma       = 1.4;
my $blast_table =
          Blast::IPS->new( 'symmetry' => $symmetry, 'gamma' => $gamma );
if ( !defined($blast_table) ) {
    die "missing table for sym=$symmetry, gamma=$gamma\n";
}

# Loop over the Goldstine-Von_Neumann points to calculate the maximum absolute
# relative error in overpressure
my $err_max;
my $lambda_max = 0;
foreach my $point ( @{$rgoldstine} ) {
    my ( $lambda_t, $ovprat_t ) = @{$point};
    last if ( defined($lambda_stop) && $lambda_t > $lambda_stop );

    # Lookup the overpressure at the given scaled range
    my $Q = log($lambda_t);
    my $ret = $blast_table->wavefront( 'X' => $Q )->{'TableVars'};
    my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};
    my $ovprat = exp($Y);
    my $lambda = exp($X);
    my $err    = abs( $ovprat - $ovprat_t ) / $ovprat_t;
    if ($VERBOSE) {
        print "lambda=$lambda, ovp=$ovprat, X=$X, err=$err\n";
    }

    if ( !defined($err_max) || $err > $err_max ) { $err_max = $err }
    if ( $lambda_t > $lambda_max ) { $lambda_max = $lambda_t }
}

if ( !ok( $err_max <= $TOL ) || $VERBOSE ) {
    my $err_max_pr    = sprintf "%0.3g", $err_max;
    my $lambda_max_pr = sprintf "%0.5g", $lambda_max;
    print
"Goldstine-Von Neumann Table max error to lambda=$lambda_max_pr: $err_max_pr\n";
}

