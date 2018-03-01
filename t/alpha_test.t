use strict;
use Test;

my $ralpha_table;

BEGIN {

  # This test checks the builtin tables by looking at the value of 'alpha' which
  # characterizes a point source at high pressure.  The value of alpha can be
  # obtained with available software and compared with the value extracted from
  # the builtin tables.

=pod

Reference:

'sedov3.f' computer program, obtained from http://cococubed.asu.edu/research_pages/sedov.shtml

=cut

    # The following alpha values are from running the sedov3 code to provide a
    # check on the value of alpha extracted from the blast tables in Blast::IPS.
    #   [symmetry, gamma, alpha ]
    $ralpha_table = [

 # NOTE: I have excluded the gamma 1.1 cases for now because the available codes
 # seem to have difficulty producing accurate values for gamma this low.
 # I got 3.41741, 3.99818, 4.50433/2=2.252165 by Aamer Haque program sedov.c

        #[2,1.1   , 3.41410750], # sedov3
        #[1,1.1   , 3.99292335], # sedov3
        #[0,1.1   , 2.24621135], # sedov3
        #[2,1.1   , 3.41741],  # sedov.c
        #[1,1.1   , 3.99818],  # sedov.c
        #[0,1.1   , 2.252165], # sedov.c

        [ 2, 1.2,   1.71980352 ],
        [ 1, 1.2,   2.00518582 ],
        [ 0, 1.2,   1.11957392 ],
        [ 2, 1.3,   1.14360009 ],
        [ 1, 1.3,   1.32797111 ],
        [ 0, 1.3,   0.734208034 ],
        [ 2, 1.4,   0.851072, ],
        [ 1, 1.4,   0.984074056 ],
        [ 0, 1.4,   0.53874280 ],
        [ 2, 1.667, 0.493319042 ],
        [ 1, 1.667, 0.563965674 ],
        [ 0, 1.667, 0.301289264 ],
        [ 2, 2,     3.11987797E-01 ],
        [ 1, 2,     0.351935918 ],
        [ 0, 2,     0.183165883 ],
        [ 2, 3,     0.132524085 ],
        [ 1, 3,     1.44403233e-01 ],
        [ 0, 3,     7.09617192e-02 ],

        # I am not doing spherical blasts for gamma greater than 3
        #[2,7 , 2.79252680e-02],
        #[1,7 , 2.74799190e-02],
        #[0,7 , 1.21555899e-02],
    ];

    my $ncases = @{$ralpha_table};
    plan tests => $ncases;
}

# Except for low values of gamma, as noted above, the program 'sedov3' appears
# to have an accuracy of better than 1.e-5.
my $TOL     = 1.e-5;
my $VERBOSE = 0;

foreach my $rcase ( @{$ralpha_table} ) {
    my ( $symmetry, $gamma, $alpha_t ) = @{$rcase};

    # Create a table for this case
    use Blast::IPS;
    my %args = ( 'symmetry' => $symmetry, 'gamma' => $gamma );
    my $blast_table = Blast::IPS->new( \%args );
    if ( !defined($blast_table) ) {
        die "missing table for sym=$symmetry, gamma=$gamma\n";
    }

    # get the value of alpha and compare
    my $alpha  = $blast_table->get_alpha();
    my $err    = abs( $alpha - $alpha_t ) / $alpha_t;
    my $err_pr = sprintf "%0.3g", $err;
    if ($VERBOSE) {
        print "alpha err=$err_pr for symmetry=$symmetry, gamma=$gamma\n";
    }
    ok( $err <= $TOL );
}
