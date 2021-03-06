use strict;
use warnings;
use Test;
use Blast::IPS;

my @test_cases;

BEGIN {

    #  Some tests cases to check the tables
    #  SYM=0,1,2 for plane, cylindrical, and spherical symmetry
    #  gamma is the ideal gas gamma
    #  ovprat is the overpressure ratio, =(P-P0)/P0
    #  lambda is the scaled range
    #
    #   [SYM, gamma, ovprat, lambda, d(lnp)/d(lnr)]
    @test_cases = (
        [ 0, 1.1,   10,    0.02133992, -0.8944939 ],
        [ 0, 1.1,   1,     0.431521,   -0.6435972 ],
        [ 0, 1.1,   0.001, 233791,     -0.5003435 ],
        [ 0, 1.2,   10,    0.04092579, -0.8951952 ],
        [ 0, 1.2,   1,     0.8315136,  -0.6406585 ],
        [ 0, 1.2,   0.001, 466516.2,   -0.5003057 ],
        [ 0, 1.3,   10,    0.05964527, -0.8953936 ],
        [ 0, 1.3,   1,     1.218027,   -0.6381148 ],
        [ 0, 1.3,   0.001, 704954,     -0.5002374 ],
        [ 0, 1.4,   10,    0.07788091, -0.8953013 ],
        [ 0, 1.4,   1,     1.598639,   -0.6358781 ],
        [ 0, 1.4,   0.001, 951743.8,   -0.5000896 ],
        [ 0, 1.667, 10,    0.1254105,  -0.8943829 ],
        [ 0, 1.667, 1,     2.608393,   -0.6310154 ],
        [ 0, 1.667, 0.001, 1654987,    -0.4992258 ],
        [ 0, 2,     10,    0.1837311,  -0.8928181 ],
        [ 0, 2,     1,     3.876623,   -0.6264955 ],
        #[ 0, 2,     0.001, 2605587,    -0.4982453 ],
        [ 0, 2,     0.001, 2605558,    -0.4982463 ],
        [ 0, 3,     10,    0.3576851,  -0.8885918 ],
        [ 0, 3,     1,     7.780978,   -0.6182701 ],
        [ 0, 3,     0.001, 5748512,    -0.4968255 ],
        [ 1, 1.1,   10,    0.08132839, -1.8205614 ],
        [ 1, 1.1,   1,     0.351888,   -1.3087360 ],
        [ 1, 1.1,   0.001, 1009.464,   -0.7554066 ],
        [ 1, 1.2,   10,    0.1122667,  -1.8188836 ],
        [ 1, 1.2,   1,     0.4886248,  -1.2992444 ],
        [ 1, 1.2,   0.001, 1443.134,   -0.7550975 ],
        [ 1, 1.3,   10,    0.1349748,  -1.8171189 ],
        [ 1, 1.3,   1,     0.5906154,  -1.2909596 ],
        [ 1, 1.3,   0.001, 1788.812,   -0.7548495 ],
        [ 1, 1.4,   10,    0.1535647,  -1.8153120 ],
        [ 1, 1.4,   1,     0.6752847,  -1.2836106 ],
        [ 1, 1.4,   0.001, 2090.885,   -0.7546326 ],
        [ 1, 1.667, 10,    0.1926937,  -1.8105006 ],
        [ 1, 1.667, 1,     0.8572118,  -1.2673588 ],
        [ 1, 1.667, 0.001, 2784.219,   -0.7542023 ],
        [ 1, 2,     10,    0.2303877,  -1.8048921 ],
        [ 1, 2,     1,     1.037228,   -1.2517761 ],
        [ 1, 2,     0.001, 3522.062,   -0.7538409 ],
        [ 1, 3,     10,    0.312814,   -1.7917042 ],
        [ 1, 3,     1,     1.444709,   -1.2215822 ],
        [ 1, 3,     0.001, 5326.871,   -0.7532539 ],
        [ 2, 1.1,   10,    0.1699955,  -2.7492500 ],
        [ 2, 1.1,   1,     0.4451088,  -2.0056088 ],
        [ 2, 1.1,   0.001, 86.56106,   -1.1024200 ],
        [ 2, 1.2,   10,    0.2105483,  -2.7458680 ],
        [ 2, 1.2,   1,     0.5535127,  -1.9914313 ],
        [ 2, 1.2,   0.001, 109.9193,   -1.1011801 ],
        [ 2, 1.3,   10,    0.2377607,  -2.7427010 ],
        [ 2, 1.3,   1,     0.6272934,  -1.9790988 ],
        [ 2, 1.3,   0.001, 126.8763,   -1.1001032 ],
        [ 2, 1.4,   10,    0.2587651,  -2.7396801 ],
        [ 2, 1.4,   1,     0.6849294,  -1.9681878 ],
        [ 2, 1.4,   0.001, 140.8054,   -1.0991398 ],
        [ 2, 1.667, 10,    0.2999177,  -2.7321424 ],
        [ 2, 1.667, 1,     0.7998176,  -1.9441100 ],
        [ 2, 1.667, 0.001, 170.427,    -1.0970483 ],
        [ 2, 2,     10,    0.3363553,  -2.7237125 ],
        [ 2, 2,     1,     0.9038715,  -1.9210264 ],
        [ 2, 2,     0.001, 199.2837,   -1.0950729 ],
        [ 2, 3,     10,    0.4076786,  -2.7040766 ],
        [ 2, 3,     1,     1.113608,   -1.8757688 ],
        [ 2, 3,     0.001, 262.1752,   -1.0913330 ],
    );

    my $ncases = @test_cases;
    plan tests => 2 * $ncases;
}

my $VERBOSE = 0;
foreach my $rcase (@test_cases) {
    my ( $SYM, $gamma, $ovprat_t, $lambda_t, $dlnp_dlnr_t ) = @{$rcase};

    # Create a table for this case
    my %args = ( 'symmetry' => $SYM, 'gamma' => $gamma );
    my $blast_table = Blast::IPS->new( \%args );
    if ( !defined($blast_table) ) {
        die "missing table for sym=$SYM, gamma=$gamma\n";
    }

    # Lookup the overpressure at the given scaled range
    my $iQ  = 'X';
    my $Q   = log($lambda_t);
    #my $ret = $blast_table->lookup( $Q, $iQ );
    my $ret = $blast_table->wavefront( $iQ => $Q )->{'TableVars'};
    my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};
    my $ovprat    = exp($Y);
    my $lambda    = exp($X);
    my $err       = abs( $ovprat - $ovprat_t ) / $ovprat_t;
    my $ovprat_pr = sprintf "%0.7g", $ovprat;
    my $lambda_pr = sprintf "%0.7g", $lambda;
    my $err_pr    = sprintf "%0.3g", $err;
    if ($VERBOSE) {
        print
"err=$err_pr for SYM=$SYM, gamma=$gamma, ovprat=$ovprat_pr, lambda=$lambda_pr\n";
    }
    ok( $err <= 1.e-5 );

    # Lookup the scaled range at the given overpressure ratio
    $iQ  = 'Y';
    $Q   = log($ovprat_t);
    #$ret = $blast_table->lookup( $Q, $iQ );
    $ret = $blast_table->wavefront( $iQ => $Q )->{'TableVars'};
    ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};
    $ovprat = exp($Y);
    $lambda = exp($X);
    $err = abs( $lambda - $lambda_t ) / $lambda_t;
    $ovprat_pr = sprintf "%0.7g", $ovprat;
    $lambda_pr = sprintf "%0.7g", $lambda;
    $err_pr    = sprintf "%0.3g", $err;

    if ($VERBOSE) {
        print
"err=$err_pr for SYM=$SYM, gamma=$gamma, ovprat=$ovprat_pr, lambda=$lambda_pr\n";
    }
    ok( $err <= 1.e-5 );
}

