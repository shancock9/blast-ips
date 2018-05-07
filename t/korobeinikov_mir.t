use strict;
use warnings;
use Test;
use Blast::IPS;

my $rmodels;

BEGIN {

=pod

This test compares some tables of values for this problem published by
Korobeinikov, Chushkin and Sharovatova.  The table for spherical symmetry with
gamma=1.4 also appears in an appendix of Korobeinikov's book "Problems of
Point-Blast Theory" (with an obvious typographical error).

These tables have errors which are typically on the order of 1 percent at high
overpressure but can become much larger at longer ranges.

REFERENCES

AIR FORCE SYSTEMS COMMAND WRIGHT-PATTERSON AFB OH FOREIGN TECHNOLOGY DIVISION, Korobeinikov, V P, Chushkin, P I, and Sharovatova, K V. 1971. Gas-Dynamic Functions of a Point Explosion. http://oai.dtic.mil/oai/oai?&verb=getRecord&metadataPrefix=html&identifier=AD0733923. 

Korobejnikov, Viktor Pavlovič. 1991. Problems of point-blast theory. New York: American Institute of Physics. 

=cut

   # The value of 'q' must be converted to overpressure
   # The value of 'Ts' must be multiplied by sqrt(gamma) to change to my scaling

    $rmodels = [
        {
            # The error in the MIR plane symmetry tables become
            # very large (to 27%) for q>0.75 so I have commented
            # out those points. The remaining error is up to about
            # 7%.
            symmetry => 0,
            gamma    => 1.4,
            table    => [

                # [q,Ts,Rs],
                [ 0.05, 0.0021, 0.0163 ],
                [ 0.1,  0.0066, 0.0358 ],
                [ 0.15, 0.0151, 0.0641 ],
                [ 0.2,  0.0277, 0.0999 ],
                [ 0.25, 0.0446, 0.142 ],
                [ 0.3,  0.0687, 0.196 ],
                [ 0.35, 0.102,  0.266 ],
                [ 0.4,  0.15,   0.358 ],
                [ 0.45, 0.219,  0.483 ],
                [ 0.5,  0.318,  0.653 ],
                [ 0.55, 0.465,  0.892 ],
                [ 0.6,  0.687,  1.24 ],
                [ 0.65, 1.04,   1.76 ],
                [ 0.7,  1.61,   2.59 ],
                [ 0.75, 2.67,   4.06 ],

                # Up to 28% error:
                #[ 0.8,  5.17,   7.4 ],
                #[ 0.85, 11.8,   16 ],
                #[ 0.9,  25.4,   33.2 ],
            ],
        },

        {
            symmetry => 0,
            gamma    => 1.667,
            table    => [

                # [q,Ts,Rs],
                [ 0.05, 0.0029, 0.0244 ],
                [ 0.1,  0.009,  0.0534 ],
                [ 0.15, 0.0196, 0.092 ],
                [ 0.2,  0.0356, 0.141 ],
                [ 0.25, 0.058,  0.202 ],
                [ 0.3,  0.0889, 0.278 ],
                [ 0.35, 0.132,  0.377 ],
                [ 0.4,  0.193,  0.506 ],
                [ 0.45, 0.283,  0.682 ],
                [ 0.5,  0.408,  0.916 ],
                [ 0.55, 0.599,  1.26 ],
                [ 0.6,  0.885,  1.74 ],
                [ 0.65, 1.34,   2.48 ],
                [ 0.7,  2.08,   3.64 ],
                [ 0.75, 3.49,   5.78 ],

                # up to 22% error:
                #[ 0.8,  6.64,   10.4 ],
                #[ 0.85, 15,     22.2 ],
                #[ 0.9,  33.2,   47.3 ],

            ],
        },

        {
            symmetry => 1,
            gamma    => 1.3,
            table    => [

                # [q,Ts,Rs],
                [ 0.05, 0.009,  0.0898 ],
                [ 0.1,  0.0192, 0.132 ],
                [ 0.15, 0.0318, 0.174 ],
                [ 0.2,  0.0468, 0.215 ],
                [ 0.25, 0.0639, 0.256 ],
                [ 0.3,  0.0846, 0.301 ],
                [ 0.35, 0.106,  0.349 ],
                [ 0.4,  0.137,  0.402 ],
                [ 0.45, 0.174,  0.466 ],
                [ 0.5,  0.22,   0.543 ],
                [ 0.55, 0.278,  0.633 ],
                [ 0.6,  0.352,  0.745 ],
                [ 0.65, 0.452,  0.89 ],
                [ 0.7,  0.596,  1.09 ],
                [ 0.75, 0.807,  1.37 ],
                [ 0.8,  1.13,   1.79 ],
                [ 0.85, 1.74,   2.55 ],
                [ 0.9,  2.96,   4.04 ],
            ],
        },

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
            symmetry => 1,
            gamma    => 1.667,
            table    => [

                # [q,Ts,Rs],
                [ 0.05, 0.0107, 0.121 ],
                [ 0.1,  0.0229, 0.178 ],
                [ 0.15, 0.0375, 0.232 ],
                [ 0.2,  0.0546, 0.286 ],
                [ 0.25, 0.0745, 0.34 ],
                [ 0.3,  0.0986, 0.399 ],
                [ 0.35, 0.128,  0.464 ],
                [ 0.4,  0.162,  0.538 ],
                [ 0.45, 0.204,  0.621 ],
                [ 0.5,  0.257,  0.72 ],
                [ 0.55, 0.324,  0.84 ],
                [ 0.6,  0.412,  0.989 ],
                [ 0.65, 0.531,  1.18 ],
                [ 0.7,  0.698,  1.44 ],
                [ 0.75, 0.946,  1.82 ],
                [ 0.8,  1.31,   2.36 ],
                [ 0.85, 2.05,   3.41 ],

                # 8.2% error:
                #[ 0.9,  3.32,   5.15 ],
            ],
        },

        {
            symmetry => 2,
            gamma    => 1.3,
            table    => [

                # [q,Ts,Rs],
                [ 0.05, 0.0146, 0.181 ],
                [ 0.1,  0.0274, 0.234 ],
                [ 0.15, 0.0415, 0.281 ],
                [ 0.2,  0.0569, 0.323 ],
                [ 0.25, 0.0737, 0.363 ],
                [ 0.3,  0.0918, 0.402 ],
                [ 0.35, 0.113,  0.444 ],
                [ 0.4,  0.137,  0.489 ],
                [ 0.45, 0.164,  0.536 ],
                [ 0.5,  0.196,  0.59 ],
                [ 0.55, 0.236,  0.653 ],
                [ 0.6,  0.287,  0.729 ],
                [ 0.65, 0.349,  0.819 ],
                [ 0.7,  0.429,  0.93 ],
                [ 0.75, 0.541,  1.08 ],
                [ 0.8,  0.713,  1.3 ],
                [ 0.85, 0.98,   1.64 ],
                [ 0.9,  1.41,   2.16 ],
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

        {
            symmetry => 2,
            gamma    => 1.667,
            table    => [

                # [q,Ts,Rs],
                [ 0.05, 0.0156, 0.22 ],
                [ 0.1,  0.0293, 0.285 ],
                [ 0.15, 0.0441, 0.34 ],
                [ 0.2,  0.06,   0.389 ],
                [ 0.25, 0.0772, 0.436 ],
                [ 0.3,  0.0978, 0.487 ],
                [ 0.35, 0.121,  0.54 ],
                [ 0.4,  0.148,  0.595 ],
                [ 0.45, 0.178,  0.655 ],
                [ 0.5,  0.213,  0.722 ],
                [ 0.55, 0.257,  0.799 ],
                [ 0.6,  0.31,   0.889 ],
                [ 0.65, 0.376,  0.996 ],
                [ 0.7,  0.461,  1.13 ],
                [ 0.75, 0.58,   1.31 ],
                [ 0.8,  0.763,  1.58 ],
                [ 0.85, 1.05,   1.99 ],
                [ 0.9,  1.52,   2.63 ],
            ],
        },
    ];

    my $ntests = @{$rmodels};
    plan tests => $ntests;
}

# We will compare the shock radius and TOA at a given overpressure.  The errors
# in TOA are comparable in magnitude to the errors in ovp.  The error in these
# particular tables is mostly on the order of 1 percent over the the initial
# part of the table then rises toward the end to about 3 percent.

# However, some of the MIR table points have particularly high errors and these
# have been commented out.
my $VERBOSE = 0;

# Create a table for this case
foreach my $model ( @{$rmodels} ) {
    my $symmetry = $model->{symmetry};
    my $gamma    = $model->{gamma};
    my $rtable   = $model->{table};
    my $blast_table =
      Blast::IPS->new( 'symmetry' => $symmetry, 'gamma' => $gamma );
    if ( !defined($blast_table) ) {
        die "missing table for sym=$symmetry, gamma=$gamma\n";
    }

    # Loop to calculate the maximum absolute relative error in range
    my $err_max;
    my $t_err_max;
    my $lambda_max = 0;

    # The plane symmetry MIR Tables have high errors
    my $TOL = ( $symmetry == 0 ) ? 0.08 : 0.04;
    foreach my $point ( @{$rtable} ) {
        my ( $q, $Ts, $Rs ) = @{$point};

        # the MIR tables use two times the energy for scaling in
        # plane symmetry compared to the present model
        if ( $symmetry == 0 ) {
            $Rs *= 2;
            $Ts *= 2;
        }

        # Convert q=(c0/D)**2 to Y1=ovp ratio for comparison with tables in
        # Korobeinikov
        my $ovprat =
          ( 2 * $gamma - ( $gamma - 1 ) * $q ) / ( ( $gamma + 1 ) * $q ) - 1;
        my $YY = log($ovprat);
        my $ret = $blast_table->wavefront( 'Y' => $YY )->{'TableVars'};
        my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};

        my $ovp    = exp($Y);
        my $z      = exp($Z);
        my $lambda = exp($X);
        my $Toa    = $lambda - $z;

        my $err = abs( $lambda - $Rs ) / $lambda;

        # multiply MIR times by sqrt of gamma to get same scaling as present
        # work
        my $t_err = abs( $Ts * sqrt($gamma) - $Toa ) / $Toa;

        if ($VERBOSE) {
            print
"sym=$symmetry, lambda=$lambda, err=$err, t_err=$t_err, Ts=$Ts, toa=$Toa\n";
        }
        if ( !defined($err_max) || $err > $err_max ) { $err_max = $err }
        if ( !defined($t_err_max) || $t_err > $t_err_max ) {
            $t_err_max = $t_err;
        }
        if ( $lambda > $lambda_max ) { $lambda_max = $lambda }
    }

    my $err_max_pr   = sprintf "%0.3g", $err_max;
    my $t_err_max_pr = sprintf "%0.3g", $t_err_max;
    if ( !ok( $err_max <= $TOL && $t_err_max <= $TOL ) ) {
        print <<EOM;
Method of Integral Relations Table for symmetry $symmetry, gamma=$gamma, has these max error 
Rel error in shock range = $err_max_pr; Relative error in TOA=$t_err_max_pr\n";
EOM
    }
}

