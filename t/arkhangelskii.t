use strict;
use warnings;
use Test;
use Blast::IPS;

my $rdata;

BEGIN {

=pod

This test compares some data points for a cylindrical blast wave with gamma=1.4 and
gamma=1.667 computed by Arkhangel'skii in the following reference.

Note that the title of the paper is misleading because it is for a cylindrical
explosion, not spherical.  There are some obvious typographical errors in the tables.

Also note that the distance scaling is based on E = E0/alpha rather than on E0, so this
has to be taken into account.  

Reference

Arkhangel'skii, N.A. 1974. "Numerical solution of the spherical explosion problem allowing for back pressure". USSR Computational Mathematics and Mathematical Physics. 14 (5): 183-192. 

=cut

    $rdata = [
        {
            symmetry => 1,
            gamma    => 1.4,
            table    => [

                # From Table 1
                # [R, T, $Pabs],
                [ 0.4040, 0.1259, 2.987 ],
                [ 2.051,  1.266,  1.277 ],
            ],
        },
        {
            symmetry => 1,
            gamma    => 1.667,
            table    => [

                # From Table 2
                # [R, T, $Pabs],
                [ 0.4040, 0.1244, 2.877 ],
                [ 2.051,  1.187,  1.271 ],
            ],
        },
    ];

    plan tests => 2;
}

# The error in these data points is about 1% so we allow 2% difference here
my $TOL     = 0.02;
my $VERBOSE = 0;

# Create a table for this case
foreach my $model ( @{$rdata} ) {
    my $symmetry    = $model->{symmetry};
    my $gamma       = $model->{gamma};
    my $rtable      = $model->{table};
    my $blast_table =
          Blast::IPS->new( 'symmetry' => $symmetry, 'gamma' => $gamma );
    if ( !defined($blast_table) ) {
        die "missing table for sym=$symmetry, gamma=$gamma\n";
    }
    my $alpha = $blast_table->get_alpha();

    # Loop to calculate the maximum absolute relative error in range
    my $err_max;
    my $t_err_max;
    my $lambda_max = 0;
    foreach my $point ( @{$rtable} ) {
        my ( $r_t, $tt, $pabs ) = @{$point};

	# Convert to my scaling
        my $lambda_t = $r_t / sqrt($alpha);
        my $tau      = sqrt($gamma) * $tt / sqrt($alpha);
        my $ovprat_t = $pabs - 1;

        my $Y_t      = log($ovprat_t);
        my $X_t      = log($lambda_t);
        my $ret      = $blast_table->wavefront( 'X' => $X_t )->{'TableVars'};
        my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};
        my $ovp    = exp($Y);
        my $z      = exp($Z);
        my $lambda = exp($X);
        my $Toa    = $lambda - $z;

        my $err   = abs( $Y - $Y_t );
        my $t_err = abs( $Toa - $tau );

        if ($VERBOSE) {
            print
"gamma=$gamma, sym=$symmetry, lambda=$lambda_t, ovp_t=$ovprat_t, ovp=$ovp, err=$err, t_err=$t_err, Ts=$tau, toa=$Toa\n";
        }
        if ( !defined($err_max) || $err > $err_max ) { $err_max = $err }
        if ( !defined($t_err_max) || $t_err > $t_err_max ) {
            $t_err_max = $t_err;
        }
        if ( $lambda > $lambda_max ) { $lambda_max = $lambda }
    }

    my $err_max_pr   = sprintf "%0.3g", $err_max;
    my $t_err_max_pr = sprintf "%0.3g", $t_err_max;
    if ( !ok( $err_max <= $TOL && $t_err_max <= $TOL ) || $VERBOSE ) {
        print <<EOM;
Arkhangelskii result for symmetry $symmetry gamma=$gamma:
Rel ovp error = $err_max_pr; Relative TOA error = $t_err_max_pr\n";
EOM
    }
}

