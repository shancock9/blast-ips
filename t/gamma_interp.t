#!/usr/bin/perl 
use strict;
use warnings;
use Test;
use Blast::IPS;

# Test interpolation to arbitrary gamma by comparing each table to the value
# that would be interpolated to it if it were not in the set.

# 0 = quiet, 1=summary table, 2=summary table + individual tables
my $verbose = $ARGV[0];

# Current allowable tolerances; these should be reduced as soon as possible
my $Y_err_tol    = 7.e-6;
my $dYdX_err_tol = 4e-5;

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
push @summary_table, "symmetry\tgamma\trerr->{Y_err}\trerr->{dYdX_err}\n";

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
            my $rilist_4 = set_interpolation_points_with_gap( $igam, $ngam, 4 );

            # At the low gamma range, do not make the points too lopsided
            my $six = $igam == 0 ? 4 : $igam == 1 ? 5 : 6;

            my $rilist_6 =
              set_interpolation_points_with_gap( $igam, $ngam, $six );
            $rtable_mid =
              Blast::IPS::_make_intermediate_gamma_table( $symmetry, $gamma,
                $rilist_6, $rilist_4 );
            my $rerr =
              compare_tables( $rtable_mid, $blast_table_mid, $verbose );
            push @summary_table,
              "$symmetry\t$gamma\t$rerr->{Y_err}\t$rerr->{dYdX_err}\n";
            my $Y_err    = $rerr->{Y_err};
            my $dYdX_err = $rerr->{dYdX_err};

            if ( !ok( $Y_err <= $Y_err_tol && $dYdX_err <= $dYdX_err_tol ) ) {
                print STDERR
"Y Error=$Y_err dYdX Error=$dYdX_err at symmetry=$symmetry gamma=$gamma\n";
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

if ( $verbose && $verbose > 1 ) {
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
    my ( $X_err_max, $Y_err_max, $dYdX_err_max, $Z_err_max, $dZdX_err_max ) =
      ( 0, 0, 0, 0, 0 );
    my @comparison;
    push @comparison,
"Y_int\tX_mid\tX_int\tX_err\tY_err\tdYdX_mid\tdYdX_int\tdYdX_err\tZ_mid\tZ_int\tZ_err\tdZdX_mid\tdZdX_err\n";
    foreach my $item ( @{$rtable_mid} ) {
        my ( $X_int, $Y_int, $dYdX_int, $Z_int, $dZdX_int ) = @{$item};
        my $item_mid = $blast_table_mid->lookup( $Y_int, 'Y' );
        my ( $X_mid, $Y_mid, $dYdX_mid, $Z_mid, $dZdX_mid ) = @{$item_mid};
        my $X_err    = $X_int - $X_mid;
        my $Y_err = $X_err * $dYdX_mid;
        my $dYdX_err = $dYdX_int - $dYdX_mid;
        my $Z_err    = $Z_int - $Z_mid;
        my $dZdX_err = $dZdX_int - $dZdX_mid;
        if ( abs($X_err) > $X_err_max )       { $X_err_max    = abs($X_err) }
        if ( abs($Y_err) > $Y_err_max )       { $Y_err_max    = abs($Y_err) }
        if ( abs($dYdX_err) > $dYdX_err_max ) { $dYdX_err_max = abs($dYdX_err) }
        if ( abs($Z_err) > $Z_err_max )       { $Z_err_max    = abs($Z_err) }
        if ( abs($dZdX_err) > $dZdX_err_max ) { $dZdX_err_max = abs($dZdX_err) }
        push @comparison,
"$Y_int\t$X_mid\t$X_int\t$X_err\t$Y_err\t$dYdX_mid\t$dYdX_int\t$dYdX_err\t$Z_mid\t$Z_int\t$Z_err\t$dZdX_mid\t$dZdX_err\n";
    }
    return {
        X_err           => $X_err_max,
	Y_err		=> $Y_err_max,
        dYdX_err        => $dYdX_err_max,
        Z_err           => $Z_err_max,
        dZdX_err        => $dZdX_err_max,
        full_comparison => \@comparison,
    };
}

sub old_compare_tables {

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
