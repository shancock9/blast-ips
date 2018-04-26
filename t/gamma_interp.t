#!/usr/bin/perl 
use strict;
use warnings;
use Test;
use Blast::IPS;

# Test interpolation to arbitrary gamma by comparing each table to the value
# that would be interpolated to it if it were not in the set.

# 0 = quiet, 1=summary table, 2=summary table + individual tables
my $verbose = 2;

# Current allowable tolerances; these should be reduced as soon as possible
my $X_err_tol    = 1.5e-5;
my $dYdX_err_tol = 5e-5;


my $rgamma_table;

INIT {

    $rgamma_table = Blast::IPS::get_rgamma_table();

    my $n0 = @{ $rgamma_table->[0] } - 2;
    my $n1 = @{ $rgamma_table->[1] } - 2;
    my $n2 = @{ $rgamma_table->[2] } - 2;

    plan tests => $n0 + $n1 + $n2;

}

my @full_comparison;
my @summary_table;
push @summary_table, "symmetry\tgamma\trerr->{X_err}\trerr->{dYdX_err}\n";

# loop over symmetry
foreach my $symmetry ( 0, 1, 2 ) {
    my $rtab = $rgamma_table->[$symmetry];

    # loop over tables for each gamma
    my $imin = 1;
    my $imax = @{$rtab} - 2;
    for ( my $i = $imin ; $i <= $imax ; $i++ ) {
        my $gamma = $rtab->[$i]->[0];

        # find this table
        my $result = Blast::IPS::_gamma_lookup( $symmetry, $gamma );
        if ( !$result ) {
            print "No result for gamma=$gamma\n";
            next;
        }
        my $table_name = $result->{table_name};
        if ($table_name) {
            my $igam = $result->{igam};
            my $ngam = @{ $rgamma_table->[$symmetry] };

            push @full_comparison, "sym=$symmetry, gamma=$gamma\n";

            # get this table
            my $blast_table_mid =
              Blast::IPS->new( symmetry => $symmetry, gamma => $gamma );
            my $gamma_chk = $blast_table_mid->get_gamma();

            my $rtable_mid;

            my $DEBUG = 0;

            my $rilist_4 = set_interpolation_points_with_gap( $igam, $ngam, 4 );

            # At the low gamma range, do not make the points too lopsided
            my $six = $igam == 0 ? 4 : $igam == 1 ? 5 : 6;
            my $rilist_6 =
              set_interpolation_points_with_gap( $igam, $ngam, $six );

            if ($DEBUG) {

                # testing local routine
                $rtable_mid =
                  _make_intermediate_gamma_table( $symmetry, $gamma,
                    $rilist_6, $rilist_4 );

            }
            else {
                $rtable_mid =
                  Blast::IPS::_make_intermediate_gamma_table( $symmetry, $gamma,
                    $rilist_6, $rilist_4 );
            }
            my $rerr =
              compare_tables( $rtable_mid, $blast_table_mid, $verbose );
            push @summary_table,
              "$symmetry\t$gamma\t$rerr->{X_err}\t$rerr->{dYdX_err}\n";
            my $X_err    = $rerr->{X_err};
            my $dYdX_err = $rerr->{dYdX_err};
            if ( !ok( $X_err <= $X_err_tol && $dYdX_err <= $dYdX_err_tol ) ) {
                print STDERR
"X Error=$X_err dYdX Error=$dYdX_err at symmetry=$symmetry gamma=$gamma\n";
            }
            push @full_comparison, @{ $rerr->{full_comparison} };
            next;
        }
    }
}

if ($verbose) {
    my $fsummary = "gamma_interp.txt";
    open( my $fh, ">", $fsummary ) || die "cannot open $fsummary: $!\n";
    foreach my $line (@summary_table) {
        $fh->print($line);
    }
    $fh->close();
}

if ( $verbose > 1 ) {
    my $fall = "gamma_interp.full_comparison.txt";
    open( my $fh, ">", $fall ) || die "cannot open $fall: $!\n";
    foreach my $line (@full_comparison) {
        $fh->print($line);
    }
    $fh->close();
}

sub set_interpolation_points_with_gap {
    my ( $jmid, $ntab, $NLAG ) = @_;

    # Find the index range for NLAG lagrange interpolation points
    # Given:
    #   $jmid the table index to be skipped (must not be end point)
    #   $ntab is the number of points in the table
    #   $NLAG is the number of interpolation points desired
    # Return:
    #   a reference to a list of consecutive indexes

    return if ( $ntab <= 0 || $NLAG <= 0 || $jmid < 0 || $jmid > $ntab - 1 );
    my $j_lo = $jmid - int( $NLAG / 2 );
    my $j_hi = $j_lo + $NLAG;
    if ( $j_lo < 0 ) {
        my $jshift = 0 - $j_lo;
        $j_lo += $jshift;
        $j_hi += $jshift;    # FIXME: optional
    }
    if ( $j_hi > $ntab - 1 ) {
        my $jshift = $ntab - 1 - $j_hi;
        $j_lo += $jshift;    # FIXME: optional
        $j_hi += $jshift;
    }
    if ( $j_lo < 0 )         { $j_lo = 0 }
    if ( $j_hi > $ntab - 1 ) { $j_hi = $ntab - 1 }
    return [ $j_lo .. $jmid - 1, $jmid + 1 .. $j_hi ];
}

sub compare_tables {

    # Make a detailed comparison of the original and interpolated tables
    my ( $rtable_mid, $blast_table_mid, $verbose ) = @_;
    my $gamma_mid = $blast_table_mid->get_gamma();
    my ( $X_err_max, $dYdX_err_max, $Z_err_max, $dZdX_err_max ) =
      ( 0, 0, 0, 0 );
    my @comparison;
    push @comparison,
"Y_int\tX_mid\tX_int\tX_err\tdYdX_mid\tdYdX_int\tdYdX_err\tZ_mid\tZ_int\tZ_err\tdZdX_mid\tdZdX_err\n";
    foreach my $item ( @{$rtable_mid} ) {
        my ( $X_int, $Y_int, $dYdX_int, $Z_int, $dZdX_int ) = @{$item};
        my $item_mid = $blast_table_mid->lookup( $Y_int, 'Y' );
        my ( $X_mid, $Y_mid, $dYdX_mid, $Z_mid, $dZdX_mid ) = @{$item_mid};
        my $X_err    = $X_int - $X_mid;
        my $dYdX_err = $dYdX_int - $dYdX_mid;
        my $Z_err    = $Z_int - $Z_mid;
        my $dZdX_err = $dZdX_int - $dZdX_mid;
        if ( abs($X_err) > $X_err_max )       { $X_err_max    = abs($X_err) }
        if ( abs($dYdX_err) > $dYdX_err_max ) { $dYdX_err_max = abs($dYdX_err) }
        if ( abs($Z_err) > $Z_err_max )       { $Z_err_max    = abs($Z_err) }
        if ( abs($dZdX_err) > $dZdX_err_max ) { $dZdX_err_max = abs($dZdX_err) }
        push @comparison,
"$Y_int\t$X_mid\t$X_int\t$X_err\t$dYdX_mid\t$dYdX_int\t$dYdX_err\t$Z_mid\t$Z_int\t$Z_err\t$dZdX_mid\t$dZdX_err\n";
    }
    return {
        X_err           => $X_err_max,
        dYdX_err        => $dYdX_err_max,
        Z_err           => $Z_err_max,
        dZdX_err        => $dZdX_err_max,
        full_comparison => \@comparison,
    };
}

# Interface routines
sub polint {
    return Blast::IPS::polint(@_);
}

sub alpha_interpolate {
    return Blast::IPS::alpha_interpolate(@_);
}

sub get_builtin_table {
    return Blast::IPS::get_builtin_table(@_);
}

sub is_monotonic_list {
    return Blast::IPS::is_monotonic_list(@_);
}

sub _locate_2d {
    return Blast::IPS::_locate_2d(@_);
}

sub _interpolate_rows {
    return Blast::IPS::_interpolate_rows(@_);
}

sub _make_intermediate_gamma_table {

    # Given: a symmetry, gamma, and a list of the indexes of the builtin tables
    #     to be interpolated
    # Return: a table of X,Y,dYdX,Z,dZdX for an intermediate gamma
    #     by using cubic interpolation between the four builtin tables

    my ( $symmetry, $gamma_x, $rilist_6, $rilist_4 ) = @_;
    $rilist_4 = $rilist_6 unless defined($rilist_4);
    my $rtable_new;
    my $alpha_x = alpha_interpolate( $symmetry, $gamma_x );
    my $A_x = -log( ( $gamma_x + 1 ) * $alpha_x );

    my $rtable_closest;
    my $dA_min;
    my $rlag_points_6;
    my $rlag_points_4;

    my $rA_6;
    foreach my $igam ( @{$rilist_6} ) {
        my ( $gamma, $table_name ) = @{ $rgamma_table->[$symmetry]->[$igam] };
        my $alpha  = alpha_interpolate( $symmetry, $gamma );
        my $rtable = get_builtin_table($table_name);
        my $A      = -log( ( $gamma + 1 ) * $alpha );
        my $dA     = ( $A_x - $A );
        if ( !defined($dA_min) || abs($dA) < $dA_min ) {
            $rtable_closest = $rtable;
        }
        push @{$rlag_points_6}, [ $rtable, $A, $gamma, $alpha, $table_name ];
        push @{$rA_6}, $A;
    }

    my $rA_4;
    foreach my $igam ( @{$rilist_4} ) {
        my ( $gamma, $table_name ) = @{ $rgamma_table->[$symmetry]->[$igam] };
        my $alpha  = alpha_interpolate( $symmetry, $gamma );
        my $rtable = get_builtin_table($table_name);
        my $A      = -log( ( $gamma + 1 ) * $alpha );
        push @{$rlag_points_4}, [ $rtable, $A, $gamma, $alpha, $table_name ];
        push @{$rA_4}, $A;
    }

    # check that A values vary monotonically;
    # otherwise interpolation will fail
    # NOTE: assuming that rA_4 is monotonic if rA_6 is
    if ( !is_monotonic_list($rA_6) ) {

        print STDERR <<EOM;
Program error: non-monotonic A
symmetry=$symmetry gamma = $gamma_x alpha=$alpha_x 
indexes of tables = @{$rilist_6}
A values are (@{$rA_6})
EOM
        return;
    }

    my $lookup = sub {

        # This is a specialized version of lookup for gamma interpolation
        # for this routine which avoids extrapolation

        # Given:
        #     $Y = a value of 'Y' = ln(overpressure) to lookup
        #     $rtab = a builtin blast table

        # Returns undef if out of bounds of table; otheerwise
        # Returns
        # 	[X, Y, dY/dX, Z]

        my ( $Y, $rtab ) = @_;
        my $icol   = 1;
        my $interp = 0;    # use cubic interpolation

        my $ntab  = @{$rtab};
        my $rhash = {
            _rtable => $rtab,
            _jl     => undef,
            _ju     => undef,
        };

        # Locate this point in the table
        my ( $il, $iu ) = _locate_2d( $rhash, $Y, $icol );

        # Handle case before start of the table
        if ( $il < 0 ) {
            return;
        }

        # Handle case beyond end of the table
        elsif ( $iu >= $ntab ) {
            return;
        }

        # otherwise interpolate within table
        else {
            return _interpolate_rows( $Y, $icol, $rtab->[$il],
                $rtab->[$iu], $interp );
        }
    };

    # Loop over all Y values in the closest table to make the interpolated
    # table.  Use the Y values of the closest table so that we converge to it.
    # But clip the Y to the lower bounds of the other table to avoid
    # extrapolation at long range where the theory is not great.

    # FIXME: an improvement would be to start at the minimum Y value of the
    # the collection of 6 tables, or better to extrapolate start the tables
    # to the desired initial Y value if necessary.

    # We are interpolating X and Z with 6 interpolation points.
    # We are interpolating dYdX and dZdX with 4 interpolation points.
    my $rtab_x = [];

    foreach my $item_ref ( @{$rtable_closest} ) {
        my $Y = $item_ref->[1];
        my ( $rX, $rdYdX, $rZ, $rdZdX );

        my $missing_item;
        foreach my $rpoint ( @{$rlag_points_6} ) {
            my ( $rtable, $A ) = @{$rpoint};
            my $item = $lookup->( $Y, $rtable );
            if ( !defined($item) ) { $missing_item = 1; last }
            my ( $X, $YY, $dYdX, $Z, $dZdX ) = @{$item};
            push @{$rX}, $X;
            push @{$rZ}, $Z;
        }

        # All values must exist in order to interpolate
        # (tables all start and end at slightly different values)
        next if ($missing_item);

        foreach my $rpoint ( @{$rlag_points_4} ) {
            my ( $rtable, $A ) = @{$rpoint};
            my $item = $lookup->( $Y, $rtable );
            if ( !defined($item) ) { $missing_item = 1; last }
            my ( $X, $YY, $dYdX, $Z, $dZdX ) = @{$item};
            push @{$rdYdX}, $dYdX;
            push @{$rdZdX}, $dZdX;
        }
        next if ($missing_item);

        my $X_x    = polint( $A_x, $rA_6, $rX );
        my $Z_x    = polint( $A_x, $rA_6, $rZ );
        my $dYdX_x = polint( $A_x, $rA_4, $rdYdX );
        my $dZdX_x = polint( $A_x, $rA_4, $rdZdX );

        push @{$rtab_x}, [ $X_x, $Y, $dYdX_x, $Z_x, $dZdX_x ];
    }
    return $rtab_x;
}
