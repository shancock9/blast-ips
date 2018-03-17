use strict;
use warnings;
use Test;
use Blast::IPS;

my $rindex;

BEGIN {

    # What is the error of the peak shock pressure using the Taylor-Sedov
    # analytic model, which ignores the ambient pressure P0?

    # Many people use the Taylor-Sedov analytic solution, sometimes well below
    # an overpressure at which it is accurate. I have seen physicists using
    # it well below shock overpressures of 1 atmosphere. 

    # This script find the range at which the error reaches 10% for each table

    # Remarkably, the error reaches 10% when the overpressure ratio is between
    # 21 and 26 atmospheres for ALL SYMMETRIES AND GAMMA values. The range at
    # which this occurs varies, of course.

=pod

Reference:


=cut

    $rindex = Blast::IPS->get_table_index();

    my $ntests = @{$rindex};
    plan tests => $ntests;
}

my $VERBOSE = 0;

# for each table, find the point of 10% overpressure error
my $err_stop=0.1;

print "symmetry\tgamma\tX\tY\n" if ($VERBOSE);
foreach my $item ( @{$rindex} ) {
    my ($table_name, $symmetry, $gamma)=@{$item};

    # Create a table for this case
    #my %args = ( 'symmetry' => $symmetry, 'gamma' => $gamma );
    my %args = ( 'table_name' => $table_name );
    my $blast_table = Blast::IPS->new( \%args );
    if ( !defined($blast_table) ) {
        die "missing table for sym=$symmetry, gamma=$gamma\n";
    }
    my $table_name_check = $blast_table->get_table_name();
    if ( $table_name ne $table_name_check ) {
        die <<EOM;
Asked for table '$table_name' but got '$table_name_check'
gamma=$gamma, symmetry=$symmetry
EOM
    }

    # Use the blast table value of alpha. Another test verifies that this
    # value is accurate.
    my $alpha   = $blast_table->get_alpha();
    my $rbounds = $blast_table->get_table_bounds();

    my $Xmin= $rbounds->[0]->[0];
    my $Xmax= $rbounds->[1]->[0];
    my $rtable = $blast_table->get_table( );

    my $TOL  = 1;
    my $dX   = 0.1;
    my $num  = int( ( $Xmax - $Xmin ) / $dX );

    my $err_max;
    my $X_stop;
    my $Y_stop;
    my $err_last;
    my ( $X0, $Y0, $dYdX0 ) = @{ $rtable->[0] };
    my $P0=exp($Y0)-1;
    my $R0=exp($X0);
    my $err;

    foreach my $item ( @{$rtable} ) {
        my ( $X, $Y, $dYdX ) = @{$item};
	my $R=exp($X);
	my $P=$P0*($R0/$R)**($symmetry+1);
        my $Y_k = log($P-1);
        $err = abs( $Y - $Y_k );
        if ( $err >= $err_stop ) {

	    # interpolate
            if ( defined($err_last) && $err_last < $err ) {
                my $ff = ( $err_stop - $err_last ) / ( $err - $err_last );
		if ($ff>1) {$ff=1}
		if ($ff<0) {$ff=0}
                my $X_int = $X_stop + $ff * ( $X - $X_stop );
                my $Y_int = $Y_stop + $ff * ( $Y - $Y_stop );
		$X_stop=$X_int;
		$Y_stop=$Y_int;
            }
            last;
        }
        $X_stop = $X;
        $Y_stop = $Y;
	$err_last=$err;
    }
    print "$symmetry\t$gamma\t$X_stop\t$Y_stop\n" if ($VERBOSE);
    my $ovp_stop=exp($Y_stop);

    # The overpressure will be between 20 and 30 atmospheres for all tables
    # for the 10% error case
    ok($ovp_stop > 20 && $ovp_stop < 30);
}

