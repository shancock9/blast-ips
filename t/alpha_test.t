use strict;
use warnings;
use Test;
use Blast::IPS;

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

       # NOTE: for low gamma values such as gamma 1.1 the available codes have
       # difficulty producing accurate values.  I haven't looked into this in
       # detail but I was able to get somewhat more accurate results by reducing
       # the tolerances in Aamer Haque's program sedov.c

        #[2,1.1   , 3.41410750], # sedov3.c
        #[1,1.1   , 3.99292335], # sedov3.c
        #[0,1.1   , 2.24621135], # sedov3.c

        [ 0, 1.1, 2.252165 ],    # sedov.c
        [ 1, 1.1, 3.99818 ],     # sedov.c
        [ 2, 1.1, 3.41741 ],     # sedov.c

        [ 0, 1.15, 1.50008064 ],
        [ 1, 1.15, 2.67423145 ],
        [ 2, 1.15, 2.289427137313403 ],

        [ 0, 1.2, 1.11957392 ],
        [ 1, 1.2, 2.00518582 ],
        [ 2, 1.2, 1.71980352 ],

        [ 0, 1.25, 0.889118185 ],
        [ 1, 1.25, 1.6002117 ],
        [ 2, 1.25, 1.37519129 ],

        [ 0, 1.3, 0.734208034 ],
        [ 1, 1.3, 1.32797111 ],
        [ 2, 1.3, 1.14360009 ],

        [ 0, 1.35, 6.22784684E-01 ],
        [ 1, 1.35, 1.13202879E+00 ],
        [ 2, 1.35, 9.76931970E-01 ],

        [ 0, 1.4, 0.53874280 ],
        [ 1, 1.4, 0.984074056 ],
        [ 2, 1.4, 0.851072, ],

        [ 0, 1.45, 4.73087918E-01 ],
        [ 1, 1.45, 8.68317887E-01 ],
        [ 2, 1.45, 7.52578653E-01 ],

        [ 0, 1.5, 0.420393169 ],
        [ 1, 1.5, 0.77524613 ],
        [ 2, 1.5, 0.673357454 ],

        [ 0, 1.55, 0.377185195 ],
        [ 1, 1.55, 0.698776277 ],
        [ 2, 1.55, 0.608236038 ],

        [ 0, 1.6, 0.341134760 ],
        [ 1, 1.6, 0.634833347 ],
        [ 2, 1.6, 0.55375103 ],

        # Note: The second shock vanishes in plane symmetry at gamma=1.6

        [ 0, 1.667, 0.301289264 ],
        [ 1, 1.667, 0.563965674 ],
        [ 2, 1.667, 0.493319042 ],

        [ 2, 1.8, 0.402941901 ],
        [ 1, 1.8, 0.45813829 ],
        [ 0, 1.8, 0.242077257 ],

        [ 0, 2, 0.183165883 ],
        [ 1, 2, 0.351935918 ],
        [ 2, 2, 3.11987797E-01 ],

        [ 2, 2.5, 0.191356397 ],
        [ 1, 2.5, 0.211963369 ],
        [ 0, 2.5, 0.106854459 ],

        [ 0, 3, 7.09617192e-02 ],
        [ 1, 3, 1.44403233e-01 ],
        [ 2, 3, 0.132524085 ],

        [ 2, 4, 0.0766253961 ],
        [ 1, 4, 0.0811006953 ],
        [ 0, 4, 0.0383442271 ],

        [ 2, 5, 0.050894114 ],
        [ 1, 5, 0.0524980647 ],
        [ 0, 5, 0.0241337021 ],

        [ 2, 6, 0.0366750787 ],
        [ 1, 6, 0.0369407716 ],
        [ 0, 6, 0.0166211108 ],

        [ 2, 6.5, 3.18160802E-02 ],

        [ 0, 7, 1.21555899e-02 ],
        [ 1, 7, 2.74799190e-02 ],

        # [2,7 , 2.79252680e-02],
    ];

    my $ncases = @{$ralpha_table};
    plan tests => 2 * $ncases;
}

my $VERBOSE = 0;

foreach my $rcase ( @{$ralpha_table} ) {
    my ( $symmetry, $gamma, $alpha_t ) = @{$rcase};

    # Except for low values of gamma, as noted above, the program 'sedov3'
    # appears have an accuracy of better than 1.e-5 over most of the
    # range but as high as 5.e-5 for cylinder gamma=1.5.
    my $TOL = $gamma < 1.2 ? 1.e-3 : 5.e-5;

    # The table_interpolate alphas get large (near .002) near gamma=1.1
    my $TOL_i = 2 * $TOL;

    # Create a table for this case
    my $blast_table =
      Blast::IPS->new( 'symmetry' => $symmetry, 'gamma' => $gamma );
    if ( !defined($blast_table) ) {
        die "missing table for sym=$symmetry, gamma=$gamma\n";
    }

    # get the value of alpha and compare
    my $alpha  = $blast_table->get_alpha();
    my $err    = abs( $alpha - $alpha_t ) / $alpha_t;
    my $err_pr = sprintf "%0.3g", $err;

    if ( !ok( $err <= $TOL ) || $VERBOSE ) {
        print STDERR
          "blast wave alpha err=$err_pr for symmetry=$symmetry, gamma=$gamma\n
alpha=$alpha but alpha_t=$alpha_t";
    }

    # Also get the alpha from the interpolation function
    my $alpha_i  = Blast::IPS::alpha_interpolate( $symmetry, $gamma );
    my $err_i    = abs( $alpha - $alpha_i ) / $alpha_i;
    my $err_i_pr = sprintf "%0.3g", $err_i;

    if ( !ok( $err_i <= $TOL_i ) || $VERBOSE ) {
        print STDERR
"table interp alpha err_i=$err_i_pr for symmetry=$symmetry, gamma=$gamma\nalpha_i=$alpha_i, alpha_t=$alpha_t, alpha=$alpha";
    }

}
