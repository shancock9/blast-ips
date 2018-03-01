use strict;
use Test;
use Blast::IPS;

my $ntests;
my @symmetries; 
my @gammas; 

BEGIN {

# This test compares the tables with a simple analytical model developed by
# Korobeinikov.

=pod

Reference:


=cut

    @symmetries = ( 0, 1, 2 );
    @gammas = ( 1.1, 1.2, 1.3, 1.4, 1.667, 2, 3 );

    my $nsym = @symmetries; 
    my $ngam = @gammas; 

    # Cases 1 and 2 long range are not yet programmed in Blast::IPS
    # So for now, we run case3 just for spherical symmetry
    #            case 1        case 2        case3 (sym=2 only)
    $ntests = $nsym*$ngam + $nsym*$ngam + $ngam;
    plan tests => $ntests;
}

my $VERBOSE = 0;

foreach my $symmetry (@symmetries) {
    foreach my $gamma (@gammas) {

        # Create a table for this case
        my %args = ( 'symmetry' => $symmetry, 'gamma' => $gamma );
        my $blast_table = Blast::IPS->new( \%args );
        if ( !defined($blast_table) ) {
            die "missing table for sym=$symmetry, gamma=$gamma\n";
        }
        my $table_name = $blast_table->get_table_name();

        # Use the blast table value of alpha. Another test verifies that this
        # value is accurate.
        my $alpha   = $blast_table->get_alpha();
        my $rbounds = $blast_table->table_bounds();

        my $Xmin_tab = $rbounds->[0]->[0];
        my $Xmax_tab = $rbounds->[1]->[0];

        my $TOL = 1;
        foreach my $case ( 1 .. 3 ) {

            my ( $Xmin, $Xmax );
            if ( $case == 1 ) {

                # high overpressure range: the error should be zero
                $Xmin = $Xmin_tab - 10;
                $Xmax = $Xmin_tab - 0.1;
                $TOL  = 0.0001;
            }
            elsif ( $case == 2 ) {

                # Within the bounds of the table
		# max errors depend strongly on the join point and symmetry
                $Xmin = $Xmin_tab + 0.1;
                $Xmax = $Xmax_tab - 0.1;
                $TOL  = ( $symmetry == 2 ) ? 0.2 : $symmetry == 1 ? 0.3 : 0.15;
            }
            elsif ( $case == 3 ) {

		# TODO: Blast::IPS not programmed for long range cyl and plane
		next unless ($symmetry == 2);

                $Xmin = $Xmax_tab + 0.1;
                $Xmax = $Xmin + 10;
                $TOL  = ( $symmetry == 2 ) ? 0.2 : $symmetry == 1 ? 0.3 : 0.15;
            }
            else {

		# undefined case
                ok(0);
            }

            # Generate a wide range of table of values for testing
            my $rtable = $blast_table->table_gen( 10, $Xmin, $Xmax );

            my $rjoin = 2;

	    # Korobeinikov suggest R*=2, but I find that these values are
	    # better.
	    # FIXME: can we improve these??
            $rjoin = ( $symmetry == 2 ) ? 1.75 : ( $symmetry == 1 ) ? 4 : 2;

            my $err_max;
            foreach my $item ( @{$rtable} ) {
                my ( $X, $Y, $dYdX ) = @{$item};
                my ( $X_k, $Y_k, $dYdX_k ) =
                  korobeinikov_analytic_model( $X, $symmetry, $gamma, $alpha,
                    $rjoin );
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
"# Max error for '$table_name' case=$case, Xmin=$Xmin, Xmax=$Xmax, sym=$symmetry, gamma=$gamma, alpha=$alpha, rjoin=$rjoin is $err_max_pr\n";
                print "$text";
            }
        }
    }
}

sub korobeinikov_analytic_model {

    my ( $X, $symmetry, $gamma, $alpha, $rjoin ) = @_;

    # This is the Korobeinikov simple analytical model for a point source
    # explosion

=pod

Lindberg, H. E., and R. D. Firth. 1967. STRUCTURAL RESPONSE OF SPINE VEHICLES. VOLUME 2. SIMULATION OF TRANSIENT SURFACE LOADS BY EXPLOSIVE BLAST WAVES. Ft. Belvoir: Defense Technical Information Center. http://handle.dtic.mil/100.2/AD816184. 

Korobejnikov, Viktor Pavlovič. 1991. Problems of point-blast theory. New York: American Institute of Physics. 

=cut

    my $lambda = exp($X);

    # The suggested rjoin=2 for spherical symmetry gives about 8.1% error at
    # lambda=7, and 5% error at long range.  This can be reduced by adjusting
    # rjoin.

    # a sub to calculate overpressure ratio, given range
    my $ovp_calc = sub {
        my ($rrat) = @_;

        my $am2;
        if ( $rrat <= $rjoin ) {

            if ( $symmetry == 0 ) {
                $am2 = 9 * $alpha * $gamma * $rrat;
            }
            elsif ( $symmetry == 1 ) {
                $am2 = 16 * $alpha * $gamma * $rrat**2;
            }
            elsif ( $symmetry == 2 ) {
                $am2 = 25 * $alpha * $gamma * $rrat**3;
            }
        }
        else {
            if ( $symmetry == 0 ) {
                $am2 = 9 * $alpha * $gamma * $rrat;
            }
            elsif ( $symmetry == 1 ) {

                # Need to review the cylindrical model ... discontinuous dYdX
                # and high error these are from p15 of sri_blast.pdf and are
                # for rjoin=2 Not sure if the cyl equation is correct for rjoin
                # other than 1
                $am2 = 16 * sqrt($rjoin) * $alpha * $gamma * $rrat**( 3 / 2 );
            }
            elsif ( $symmetry == 2 ) {
                $am2 =
                  ( 25 * $rjoin * $alpha * $gamma * $rrat**2 ) *
                  ( log( $rrat / $rjoin ) + 1 );
            }
        }

        my $up_over_c = 4 / ( $gamma + 1 ) / sqrt($am2);

        my $ovp_atm =
          4 * $gamma / $am2 / ( $gamma + 1 ) * ( 1 + sqrt( 1 + $am2 ) );

        # for reference:
        # This is the equation on page 182 of point blast theory
        # It is missing a factor of gamma.
        # It is a cleaner way to code and typeset if you include the
        # missing factor of gamma.
        # my $ovp_atm_check = 4 / ( $gamma + 1 ) / ( -1 + sqrt( 1 + $am2 ) );
        #  my $diff = ( $ovp_atm - $ovp_atm_check );

        return $ovp_atm;
    };

    # Eq 20:
    my $ovp_atm = $ovp_calc->($lambda);
    my $Y       = log($ovp_atm);

    # differentiate to get d(ln P)/d(ln R)
    my $delta_r = 0.001 * $lambda;
    my ($ovpm)  = $ovp_calc->( $lambda - $delta_r );
    my ($ovpp)  = $ovp_calc->( $lambda + $delta_r );
    my $dYdX    = 0.5 * $lambda * ( $ovpp - $ovpm ) / ( $delta_r * $ovp_atm );
    return ( $X, $Y, $dYdX );
}

