use strict;
use warnings;
use Test;
use Blast::IPS;

my $rindex;

BEGIN {

 # This test checks the builtin tables by looking at the value of 'alpha'
 # which characterizes a point source at high pressure.  These values are
 # obtained by fitting the high pressure part of the shock tables.  These values
 # typically have errors below 1.e-6 but for low gamma (about 1.1) the errors
 # may be a little larger.

    # We compare these values with alpha values calculated analytically from
    # the similarity solution.  Tables of the analytical values are contained
    # in Blast::IPS::AlphaTable.pm and they can be extracted with the
    # appropriate call.  These tabulated values have estimated maximum errors of
    # about 1.e-7 or less.

    $rindex = Blast::IPS->get_index();

    ##my $ntests = 0 + keys %{$rindex};
    my $ntests =
      @{ $rindex->[0] } +
      @{ $rindex->[1] } +
      @{ $rindex->[2] };    #0 + keys %{$rindex};
    plan tests => $ntests;
}

my $VERBOSE = $ARGV[0];

# In most cases the differences are below 1.e-6. The current maximum error is
# 5.75-06 for plane symmetry, gamma = 1.1.  So we set the acceptable tolerance
# just above that.
my $tol = 6.e-6;

print "symmetry\tgamma\talpha1\talpha2\tdiff\n" if ($VERBOSE);
foreach my $symmetry ( 0 .. 2 ) {
    foreach my $item ( @{ $rindex->[$symmetry] } ) {
        my ( $gamma, $table_name ) = @{$item};

        # Create a table for this case
        my $blast_table = Blast::IPS->new( 'table_name' => $table_name );
        if ( !defined($blast_table) ) {
            die "missing table for sym=$symmetry, gamma=$gamma\n";
        }

        # Now compare these two different alpha values:
        # $alpha = the value estimated from the finite difference calculation
        # $alpha2 = the value computed by integrating the analytical solution

        # This gets the value extracted from the FD calculation:
        my $alpha = $blast_table->get_alpha();

       # The easiest way to get the analytical value is with the following call.
       # These values are in a table called AlphaTable.pm; no interpolation will
       # be actually needed.
        my $alpha2 = Blast::IPS::alpha_interpolate( $symmetry, $gamma );

        my $dalp = abs( $alpha2 - $alpha );
        ok( $dalp < $tol ) || print STDERR <<EOM;
Large d(alpha) for sym=$symmetry, gamma=$gamma: | $alpha - $alpha2 | = $dalp
EOM
        print "$symmetry\t$gamma\t$alpha\t$alpha2\t$dalp\n" if ($VERBOSE);
    }
}

