use strict;
use Test;

my $rmodels;

BEGIN {


=pod

This test compares some tables of values for this problem published by
Korobeinikov, Chushkin and Sharovatova.  The table for spherical symmetry with
gamma=1.4 also appears in an appendix of Korobeinikov's book "Problems of
Point-Blast Theory" (but with an obvious typographical error).

These tables have errors which are typically on the order of 1 percent at high
overpressure but can become much larger at longer ranges.  Here we will compare
with the spherical and cylindrical blast tables for gamma=1.4.

REFERENCES

AIR FORCE SYSTEMS COMMAND WRIGHT-PATTERSON AFB OH FOREIGN TECHNOLOGY DIVISION, Korobeinikov, V P, Chushkin, P I, and Sharovatova, K V. 1971. Gas-Dynamic Functions of a Point Explosion. http://oai.dtic.mil/oai/oai?&verb=getRecord&metadataPrefix=html&identifier=AD0733923. 

Korobejnikov, Viktor Pavlovič. 1991. Problems of point-blast theory. New York: American Institute of Physics. 

=cut


    # The value of 'q' must be converted to overpressure
    # The value of 'Ts' must be multiplied by sqrt(gamma) to change to my scaling

    $rmodels = [
        {
            symmetry => 1,
            gamma    => 1.4,
            table    => [

                # [q,Ts,Rs],
                [ 0.05, 0.0097, 0.0999 ],
                [ 0.1,  0.0207, 0.148 ],
                [ 0.15, 0.0339, 0.193 ],
                [ 0.2,  0.0498, 0.238 ],
                [ 0.25, 0.0687, 0.285 ],
                [ 0.3,  0.0901, 0.333 ],
                [ 0.35, 0.116,  0.387 ],
                [ 0.4,  0.148,  0.448 ],
                [ 0.45, 0.186,  0.518 ],
                [ 0.5,  0.234,  0.599 ],
                [ 0.55, 0.295,  0.7 ],
                [ 0.6,  0.376,  0.826 ],
                [ 0.65, 0.485,  0.988 ],
                [ 0.7,  0.637,  1.21 ],
                [ 0.75, 0.861,  1.52 ],
                [ 0.8,  1.2,    1.98 ],
                [ 0.85, 1.86,   2.83 ],
                [ 0.9,  3.16,   4.48 ],
            ],
        },
        {
            symmetry => 2,
            gamma    => 1.4,
            table    => [

                # [q,Ts,Rs],
                [ 0.05, 0.0151, 0.194 ],
                [ 0.1,  0.0284, 0.252 ],
                [ 0.15, 0.0425, 0.3 ],
                [ 0.2,  0.058,  0.344 ],
                [ 0.25, 0.0755, 0.388 ],
                [ 0.3,  0.095,  0.432 ],
                [ 0.35, 0.116,  0.477 ],
                [ 0.4,  0.141,  0.524 ],
                [ 0.45, 0.17,   0.576 ],
                [ 0.5,  0.204,  0.635 ],
                [ 0.55, 0.245,  0.702 ],
                [ 0.6,  0.296,  0.781 ],
                [ 0.65, 0.359,  0.876 ],
                [ 0.7,  0.441,  0.994 ],
                [ 0.75, 0.557,  1.16 ],
                [ 0.8,  0.735,  1.39 ],
                [ 0.85, 1.02,   1.76 ],
                [ 0.9,  1.47,   2.33 ],
            ],
        },
    ];

    plan tests => 2;
}

# We will compare the shock radius and TOA at a given overpressure.
# The errors in TOA are comparable in magnitude to the errors in ovp.
# The error in these particular tables is on the order of 1 percent over the
# the initial part of the table then rises toward the end to about 3 percent.
my $TOL     = 0.04;
my $VERBOSE = 0;

# Create a table for this case
use Blast::IPS;

foreach my $model ( @{$rmodels} ) {
    my $symmetry    = $model->{symmetry};
    my $gamma    = $model->{gamma};
    my $rtable      = $model->{table};
    my %args        = ( 'symmetry' => $symmetry, 'gamma' => $gamma );
    my $blast_table = Blast::IPS->new( \%args );
    if ( !defined($blast_table) ) {
        die "missing table for sym=$symmetry, gamma=$gamma\n";
    }

    # Loop to calculate the maximum absolute relative error in range
    my $err_max;
    my $t_err_max;
    my $lambda_max = 0;
    foreach my $point ( @{$rtable} ) {
        my ( $q, $Ts, $Rs ) = @{$point};

        # Convert q=(c0/D)**2 to Y1=ovp ratio for comparison with tables in
        # Korobeinikov
        my $ovprat =
          ( 2 * $gamma - ( $gamma - 1 ) * $q ) / ( ( $gamma + 1 ) * $q ) - 1;
        my $YY = log($ovprat);
        my $ret = $blast_table->lookup( $YY, 'Y' );
        my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};

        my $ovp    = exp($Y);
        my $z      = exp($Z);
        my $lambda = exp($X);
        my $Toa    = $lambda - $z;

        my $err   = abs( $lambda - $Rs ) / $lambda;

	# multiply MIR times by sqrt of gamma to get same scaling as present work
        my $t_err = abs( $Ts * sqrt($gamma) - $Toa ) / $Toa;

	if ($VERBOSE) {
          print "sym=$symmetry, lambda=$lambda, err=$err, t_err=$t_err, Ts=$Ts, toa=$Toa\n";
	}
        if ( !defined($err_max) || $err > $err_max ) { $err_max = $err }
        if ( !defined($t_err_max) || $t_err > $t_err_max ) {
            $t_err_max = $t_err;
        }
        if ( $lambda > $lambda_max ) { $lambda_max = $lambda }
    }

    my $err_max_pr   = sprintf "%0.3g", $err_max;
    my $t_err_max_pr = sprintf "%0.3g", $t_err_max;
    if (!ok( $err_max <= $TOL && $t_err_max <= $TOL)) {
        print <<EOM;
Method of Integral Relations Table for symmetry $symmetry, gamma=$gamma, has these max error 
Rel error in shock range = $err_max_pr; Relative error in TOA=$t_err_max_pr\n";
EOM
    }
}

