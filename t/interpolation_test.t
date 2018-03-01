use strict;
use Test;

my $rspherical_n16000;

BEGIN {

    # In this test we check that the cubic and linear interpolation errors
    # are about as expected.

    # The following table of points are for a spherical blast, gamma=1.4.
    # I took them from a calculation with 16000 points between the center of
    # explosion and the blast front. For this calculation the Energy error
    # to P=0.1 P0 was 5.05e-08. The maximum relative errors in interpolated
    # shock overpressure depend on these factors:

    # The maimumum finite difference Error (based on N/2 run) is 2.05e-7
    # The maximum interpolation error using cubic splines is 5.0741e-07
    # The maximum interpolation error using linear interpolation is 0.001

    # So if we use cubic interpolation the maximum error should be about 7.e-7
    # If we use linear interpolation the maximum error should be about 0.001

    #  [lambda,tau,(P-P0)/P0],
    $rspherical_n16000 = [
        [ 0.0053911547001131,  2.32944067480237E-06, 999832.708579685 ],
        [ 0.00810725153990657, 6.45997179576777E-06, 294003.937692866 ],
        [ 0.0121917346441529,  1.79146511287645E-05, 86453.2802532769 ],
        [ 0.0183340054157294,  4.96800295697442E-05, 25422.5283101106 ],
        [ 0.0275707899152082,  0.00013776554613493,  7476.30627845435 ],
        [ 0.0414611231598561,  0.0003819876864135,   2199.18001044867 ],
        [ 0.0623494915807821,  0.00105873758569414,  647.428105274647 ],
        [ 0.093761548269244,   0.00293058418814978,  191.127732703807 ],
        [ 0.140995078681226,   0.00807583452199088,  56.9424617330524 ],
        [ 0.212027283865301,   0.0219516345847545,   17.4388063881383 ],
        [ 0.318825404980837,   0.05746908087173,     5.72429692863669 ],
        [ 0.479454151258903,   0.139422513635584,    2.12162115731283 ],
        [ 0.72092210127718,    0.304187583530025,    0.904973965084675 ],
        [ 1.08412185712612,    0.598280378264377,    0.433953603756783 ],
        [ 1.63024765756362,    1.08396196742379,     0.226285981191879 ],
        [ 2.45142852435918,    1.8522654790345,      0.124944870699227 ],
        [ 3.68638867735128,    3.04042287882942,     0.0717495025237446 ],
        [ 5.54319772386355,    4.85522557550231,     0.0423612369279363 ],
        [ 8.33515598840773,    7.60897121403397,     0.0255176255375715 ],
        [ 12.5340489700162,    11.7726899897313,     0.0156019980842465 ],
        [ 18.8358437626056,    18.0418430964733,     0.0096565016653828 ],
    ];

    my $num = @{$rspherical_n16000};
    plan tests => 2;    #* $num;
}

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
my $iQ      = 'X';
my $VERBOSE = 0;
foreach my $interp ( 0 .. 1 ) {
    my $err_max;
    my $Z_err_max;
    my $lambda_max = 0;
    my $TOL = $interp == 0 ? 7.e-7 : 1.e-3;
    foreach my $point ( @{$rspherical_n16000} ) {
        my ( $lambda_t, $tau_t, $ovprat_t ) = @{$point};
        my $z_t = $lambda_t - $tau_t;
        my $Z_t = log($z_t);

        # Lookup the overpressure at the given scaled range
        my $Q = log($lambda_t);
        my $ret = $blast_table->lookup( $Q, $iQ, $interp );
        my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};
        my $ovprat = exp($Y);
        my $lambda = exp($X);
        my $err    = abs( $ovprat - $ovprat_t ) / $ovprat_t;
        my $Z_err  = abs( $Z - $Z_t );

        #print STDERR "lambda=$lambda, err=$err\n";
        if ( !defined($err_max) || $err > $err_max ) { $err_max = $err }
        if ( !defined($Z_err_max) || $Z_err > $Z_err_max ) {
            $Z_err_max = $Z_err;
        }
        if ( $lambda_t > $lambda_max ) { $lambda_max = $lambda_t }
    }

    my $err_max_pr    = sprintf "%0.3g", $err_max;
    my $Z_err_max_pr  = sprintf "%0.3g", $Z_err_max;
    my $lambda_max_pr = sprintf "%0.5g", $lambda_max;

    #    my $lambda_pr = sprintf "%0.7g", $lambda;
    my $interp_str = $interp == 0 ? 'cubic' : 'linear';
    if ($VERBOSE) {
        print
"Spherical blast data $interp_str interpolation test using tolerance $TOL..\n";
        print
"maximum Z relative error = $Z_err_max_pr; max ovp ratio relative error=$err_max_pr\n";
    }

    #print STDERR ""( $err_max <= $TOL )";
    ok( $err_max <= $TOL );
}

