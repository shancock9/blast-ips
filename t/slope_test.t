use strict;
use warnings;
use Test;
use Blast::IPS;

BEGIN {

    # Make a simple check the shock front slopes. These have simple values in
    # spherical symmetry, gamma=1.4, in the similarity region. We will use a
    # small radius in the test case which should give slopes very close to these
    # theoretical limits.

# Reference:
#  INITIAL DECAY OF FLOW PROPERTJES OF PLANAR, CYLINDRJCAL AND SPHERICAL BLAST WAVES
#  H. S. I. Sadek and J. J. Gott1ieb
#  UTIAS Technica1 Note No. 244, CN ISSN 0082-5263
#  October, 1983
#  See p22, eq 4.12, 4.13

    plan tests => 1;
}

my $VERBOSE = 0;

my $symmetry = 2;
my $gamma    = 1.4;
my $rs       = 1.e-4;
my $tol      = 1.e-6;

# Create a table for this case
my $blast_table = Blast::IPS->new( 'symmetry' => $symmetry, 'gamma' => $gamma );

my $err = $blast_table->get_error();
if ($err) { die "$err\n"; }

my $ret    = $blast_table->wavefront( 'X' => log($rs) );
my $Y      = $ret->{Y};
my $dYdX   = $ret->{dYdX};
my $dpdr_t = $ret->{dpdr_t};
my $dudr_t = $ret->{dudr_t};
my $dpdt_r = $ret->{dpdt_r};
my $dudt_r = $ret->{dudt_r};

# Calculate shock front particle velocity
my $y    = exp($Y);
my $term = $y * ( $gamma + 1 ) / ( 2 * $gamma );
my $m    = sqrt( 1 + $term );
my $q    = 1 / $m**2;
my $up   = $y / ( $gamma * $m );
my $ps   = exp($Y);

# The theoretical values are given with these scale factors
my $dpdr_t_scaled = $rs / $ps * $dpdr_t;
my $dudr_t_scaled = $rs / $up * $dudr_t;
my $dpdt_r_scaled = $rs / ( $up * $ps ) * $dpdt_r;
my $dudt_r_scaled = $rs / $up**2 * $dudt_r;

my $dpdr_t_scaled_exact = 67 / 6;
my $dudr_t_scaled_exact = 13 / 6;
my $dpdt_r_scaled_exact = -17;
my $dudt_r_scaled_exact = -22 / 5;

my $err_1 = abs( $dpdr_t_scaled - $dpdr_t_scaled_exact );
my $err_2 = abs( $dudr_t_scaled - $dudr_t_scaled_exact );
my $err_3 = abs( $dpdt_r_scaled - $dpdt_r_scaled_exact );
my $err_4 = abs( $dudt_r_scaled - $dudt_r_scaled_exact );

if ( !ok( $err_1 < $tol && $err_2 < $tol && $err_3 < $tol && $err_4 < $tol ) ) {
    print STDERR <<EOM;
slope errors too large for symmetry=$symmetry, gamma=$gamma, rs=$rs:
err_1 = $err_1 = | $dpdr_t_scaled - $dpdr_t_scaled_exact |
err_2 = $err_2 = | $dudr_t_scaled - $dudr_t_scaled_exact |
err_3 = $err_3 = | $dpdt_r_scaled - $dpdt_r_scaled_exact |
err_4 = $err_4 = | $dudt_r_scaled - $dudt_r_scaled_exact |
EOM
}

