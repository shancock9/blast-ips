use strict;
use warnings;
use Test;
use Blast::IPS;

BEGIN {

=pod

John Miles wrote an important paper discussing the problem of the decay of
spherical blast waves at long range.  

He used the available numerical calculations to estimate values for the two
parameters A and B that define the asymptotic decay, and he concluded that the
analytic asymptote could be used to extrapolate blast strength beyond lambda=10
with about 1% accuracy.

In this test we verify this by comparing the asymptotic model with the builtin table
for a spherical blast with gamma=1.4.

Reference:

Miles, John W. 1967. "Decay of Spherical Blast Waves". Physics of Fluids. 10 (12): 2706

=cut

    plan tests => 1;
}

miles_test();

sub miles_test {

    my $VERBOSE  = 0;
    my $symmetry = 2;
    my $gamma    = 1.4;

    # Create a table for this case
    my %args = ( 'symmetry' => $symmetry, 'gamma' => $gamma );
    my $blast_table = Blast::IPS->new( \%args );
    if ( !defined($blast_table) ) {
        die "missing table for sym=$symmetry, gamma=$gamma\n";
    }
    my $table_name = $blast_table->get_table_name();
    my $rbounds    = $blast_table->get_table_bounds();
    my $Xmin_tab   = $rbounds->[0]->[0];
    my $Xmax_tab   = $rbounds->[1]->[0];

    # We are looking for about 1% max error at long range (lambda > 10)
    # Here we start closer in and accept a little more error
    my $TOL  = 0.02;
    my $Xmin = log(1.75);
    my $Xmax = $Xmax_tab + 10;

    # Generate a table of values for testing
    my $rtable = $blast_table->table_gen( 100, $Xmin, $Xmax );

    my $err_max;
    foreach my $item ( @{$rtable} ) {
        my ( $X,   $Y,   $dYdX )   = @{$item};
        my ( $X_k, $Y_k, $dYdX_k ) = miles_long_range($X);
        my $err = abs( $Y - $Y_k );
        if ( !defined($err_max) || $err > $err_max ) { $err_max = $err }

        if ($VERBOSE) {
            print
              "gamma=$gamma, sym=$symmetry, X=$X, Y=$Y, Yk=$Y_k, err=$err\n";
        }
    }
    my $err_max_pr = sprintf "%0.3g", $err_max;

    if ( !ok( $err_max <= $TOL ) ) {
        my $text =
"# Max error for Miles long range model using '$table_name' Xmin=$Xmin, Xmax=$Xmax, is $err_max_pr\n";
        print "$text";
    }
}

sub miles_long_range {

    my ( $X, $A, $B, $gamma ) = @_;

    # This is a model for the long range overpressure of a spherical blast. The
    # same type of model is used in the long range part of Korobeinikov's
    # analytical model but with different coefficients.

    # Reference:
    # Miles, John W. 1967. "Decay of Spherical Blast Waves". Physics of Fluids. 10 (12): 2706

    # Miles estimated that A falls in the range 0.23-0.24 and B in the range 0.25-0.63
    # Miles concluded that this equation could be used to extrapolate blast strength
    # beyond lambda=10 with about 1% accuracy.

    # {NOTE: The parameter A does have a limiting value of about 0.23, but B
    # varies over the entire history of a spherical shock. But using a single
    # fixed value of B can give a good approximation at long range}

    # The relative error in overpressure at very long range for these
    # parameters is about 0.006:
    $A     = 0.23 unless defined($A);
    $B     = 0.25 unless defined($B);
    $gamma = 1.4   unless defined($gamma);

    # Cutoff the solution where the decay rate equals -3 (about lambda=1)
    # From the exact relationship for the blast slope
    #    -$dlnp_dlnr = 1 + 0.5 / ( log($lambda) + $B );
    # we can derive the Value of lambda where the slope reaches -3. Since this
    # is the maximum slope (at the origin), this is a good cutoff point.
    # The prediction quickly goes bad closer than this.
    # This gives lambda_min = 1.03 for the default parameters.
    # But I suggest using this model only above lambda=1.75
    my $lambda_min = exp(0.25) - $B;
    my $X_min      = log($lambda_min);

    # print STDERR "lambda_min=$lambda_min, Xmin=$X_min\n";

    my $ovp_calc = sub {
        my ($XX) = @_;
        my $lambda = exp($XX);
        my $term = $XX + $B;
        my $ovp_atm = 1.e10;
        if ( $term > 0 ) {
            $ovp_atm = $gamma * $A / ( $lambda * sqrt($term) );
        }
        return ($ovp_atm);
    };

    # Clip the solution to the range above X_min to avoid numerical trouble
    my $X_clip = $X >= $X_min ? $X : $X_min;
    my $ovp_atm = $ovp_calc->($X_clip);
    my $Y = log($ovp_atm);

    # logarithmic derivative
    my $dYdX = -( 1 + 0.5 / ( $X + $B ) );

    # Here is the second derivative if needed:
    # d2(lnp)/d2(lnr)
    # my $d2lnp_d2lnr = 2 * ( -$dlnp_dlnr - 1 )**2;

    return ( $X, $Y, $dYdX );
}
