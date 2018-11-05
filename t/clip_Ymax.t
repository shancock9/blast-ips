use strict;
use warnings;
use Test;
use Blast::IPS;

# In this test we will check the accuracy of extrapolation in the high overpressure
# range by comparing results for each table with the results that would be obtained
# by extrapolation if that table is clipped to a certain Ymax.

# This is a new routine and is being used for model development. Only gamma=1.4
# is being done now, and the error tolerances are still set fairly high.

my $rindex;

BEGIN {

    $rindex = Blast::IPS::get_index();

    my $ntests =
      @{ $rindex->[0] } +
      @{ $rindex->[1] } +
      @{ $rindex->[2] };    

    # TESTING: only doing three tests for now
    $ntests = 3;

    plan tests => $ntests;
}

# 0 = quiet, 1=summary table, 2=summary table + individual tables
my $VERBOSE = $ARGV[0];
my @summary_table;
my @full_comparison;

# Maximum allowable relative errors for success at any gamma...
# My goal is to reduce these if possible by improving the model.
# [$Y_err_tol, $dYdX_err_tol]
my $rtol = [ [ 0.002, 0.002 ], [ 2.e-4, 3.e-4 ], [ 8.e-4, 3.e-3 ], ];

# We are removing all high pressure points down to an overpressure ratio of 100
# for this test (Y=4.605). Normally the tables go up to an overpressure ratio
# of 162000 psi (Y=12).
my $Ymax_clip = log(100);
my $Ymin_clip;
foreach my $symmetry ( 0 .. 2 ) {

    foreach my $item(@{$rindex->[$symmetry]}) {
        my ( $gamma, $table_name ) = @{$item};

	# TESTING: Only doing gamma=1.4 for all symmetries for now
	next unless ($gamma == 1.4 );

        my ( $Y_err_tol, $dYdX_err_tol ) = @{ $rtol->[$symmetry] };

        # Create a table for this case
        my $blast_table = Blast::IPS->new( 'table_name' => $table_name );
        if ( !defined($blast_table) ) {
            die "missing table for sym=$symmetry, gamma=$gamma\n";
        }

        # Create a clipped table
        my $blast_table_clipped = Blast::IPS->new(
            'symmetry'  => $symmetry,
            'gamma'     => $gamma,
            'Ymax_clip' => $Ymax_clip,
            'Ymin_clip' => $Ymin_clip,
        );

        my $rerr =
          compare_tables( $blast_table, $blast_table_clipped, $VERBOSE );
        push @summary_table,
              "$symmetry\t$gamma\t$rerr->{Y_err}\t$rerr->{dYdX_err}\n";
        my $Y_err    = $rerr->{Y_err};
        my $dYdX_err = $rerr->{dYdX_err};

        if ( !ok( $Y_err <= $Y_err_tol && $dYdX_err <= $dYdX_err_tol ) ) {
            print STDERR
"Y Error=$Y_err dYdX Error=$dYdX_err at symmetry=$symmetry gamma=$gamma\n";
        }

        #print "$symmetry\t$gamma\t$Y_err\n" if ($VERBOSE);
        push @full_comparison, @{ $rerr->{full_comparison} };
    }

}

if ($VERBOSE) {
    my $fsummary = "clip_Ymax.txt";
    open( my $fh, ">", $fsummary ) || die "cannot open $fsummary: $!\n";
    foreach my $line (@summary_table) {
        $fh->print($line);
    }
    $fh->close();
}

if ( $VERBOSE && $VERBOSE > 1 ) {
    my $fall = "clip_Ymax.full_comparison.txt";
    open( my $fh, ">", $fall ) || die "cannot open $fall: $!\n";
    foreach my $line (@full_comparison) {
        $fh->print($line);
    }
    $fh->close();
}

sub compare_tables {

    # Make a detailed comparison of the two objects
    my ( $blast_table, $blast_table_clipped, $VERBOSE ) = @_;
    my ( $X_err_max, $Y_err_max, $dYdX_err_max, $Z_err_max, $dZdX_err_max ) =
      ( 0, 0, 0, 0, 0 );
    my @comparison;
    push @comparison,
"Y_int\tX\tX_int\tX_err\tY_err\tdYdX\tdYdX_int\tdYdX_err\tZ\tZ_int\tZ_err\tdZdX\tdZdX_err\n";

    my $rbounds = $blast_table->get_table_bounds();
    my $Xmin    = $rbounds->[0]->[0];
    my $Xmax    = $rbounds->[1]->[0];

    # Get the interpolated shock table
    my $rtable = $blast_table->get_table();

    foreach my $item ( @{$rtable} ) {
        my ( $X_int, $Y_int, $dYdX_int, $Z_int, $dZdX_int ) = @{$item};
        my $item = $blast_table_clipped->wavefront( 'Y' => $Y_int )->{'TableVars'};
        my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$item};
        my $X_err    = $X_int - $X;
        my $Y_err    = $X_err * $dYdX;
        my $dYdX_err = $dYdX_int - $dYdX;
        my $Z_err    = $Z_int - $Z;
        my $dZdX_err = $dZdX_int - $dZdX;

        # Tentative: do not include points out of bounds in comparison
        if ( $X >= $Xmin && $X <= $Xmax ) {
            if ( abs($X_err) > $X_err_max ) { $X_err_max = abs($X_err) }
            if ( abs($Y_err) > $Y_err_max ) { $Y_err_max = abs($Y_err) }
            if ( abs($dYdX_err) > $dYdX_err_max ) {
                $dYdX_err_max = abs($dYdX_err);
            }
            if ( abs($Z_err) > $Z_err_max ) { $Z_err_max = abs($Z_err) }
            if ( abs($dZdX_err) > $dZdX_err_max ) {
                $dZdX_err_max = abs($dZdX_err);
            }
        }
        if ($Y_err != 0) {
            push @comparison,
"$Y_int\t$X\t$X_int\t$X_err\t$Y_err\t$dYdX\t$dYdX_int\t$dYdX_err\t$Z\t$Z_int\t$Z_err\t$dZdX\t$dZdX_err\n";
	}
    }
    return {
        X_err           => $X_err_max,
        Y_err           => $Y_err_max,
        dYdX_err        => $dYdX_err_max,
        Z_err           => $Z_err_max,
        dZdX_err        => $dZdX_err_max,
        full_comparison => \@comparison,
    };
}
