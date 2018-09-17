package Blast::IPS::ShockUtils;

# MIT License
# Copyright (c) 2018 Steven Hancock
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Some functions for shock waves in ideal gases

use strict;
use warnings;
use 5.006;

my $VERSION = 1.00;
use Carp;

our @EXPORT_OK = qw(
  shock_front_values
);
use Exporter;
our @ISA = qw(Exporter);

sub shock_front_values {

    my ( $symmetry, $gamma, $X, $Y, $dYdX ) = @_;

    # Given:
    #   symmetry = (0,1,2) 
    #   gamma
    #   X = ln(range)
    #   Y = ln(ovp) = ln(overpressure ratio) 
    #   dYdX = d(ln P)/d(ln R) along shock front
    #
    # Compute some shock front values and return in a hash ref
    #
    # FIXME: some of these expressions would benefit from Taylor series
    # expansions for small ovp

    my $p_amb    = 1; 
    my $sspd_amb = 1; 
    my $rho_amb  = $gamma; 

    my $lambda = exp($X);
    my $ovprat = exp($Y);
    my $ovp    = $p_amb * $ovprat;
    my $pabs   = $ovp + $p_amb;

    if ( $pabs < 0 ) { print STDERR "Negative abs p in sadek\n"; return }

    # Compute some shock front properties ...

    # D = Shock speed D
    my $term     = $ovprat * ( $gamma + 1 ) / ( 2 * $gamma );
    my $mshock   = sqrt( 1 + $term );
    my $D   = $sspd_amb * $mshock;
    my $D_minus_c0  = $D - $sspd_amb;

    # Taylor series expansion of the square root for small overpressures
    # Switch at a point which keeps error < about 1.e-17
    if ( $term < 1.e-4 ) {
        $D_minus_c0 = $term * ( 0.5 + $term * ( -0.125 + $term / 16 ) );
    }

    # up = shock particle velocity
    my $up = $ovprat * $p_amb / ( $rho_amb * $D );

    # dup/dovp along shock front
    my $dudp_sf =
      $sspd_amb *
      $sspd_amb / $gamma / $D *
      ( 1 - ( $gamma + 1 ) / 4 * $up / $D );

    # density ratio and density
    my $top  = ( $gamma + 1 ) * $pabs + ( $gamma - 1 ) * $p_amb;
    my $bot  = ( $gamma - 1 ) * $pabs + ( $gamma + 1 ) * $p_amb;
    my $rhorat = ( $bot > 0 ) ? $top / $bot : 1;
    my $rho = $rho_amb * $rhorat;

    # Sound speed
    my $rhocsq  = $gamma * $pabs;
    my $sspd_sq = $rhocsq / $rho;
    if ( $sspd_sq < 0 ) { $sspd_sq = 0 }
    my $aa = sqrt($sspd_sq);

    # sigma and Sigma, the C+ characteristic variable
    my $sspdx = ( $aa - $sspd_amb ); 
    my $sigma = ( 2 / ( $gamma - 1 ) ) * $sspdx;
    my $Sigma = ( $symmetry == 0 ) ? $sigma : $sigma * $lambda**( $symmetry / 2 );

    # Compute the slopes of the wave profile at the shock front ...
    #   dp/dr, du/dr, dp/dt, du/dt

    # Reference for slopes:
    #  INITIAL DECAY OF FLOW PROPERTJES OF PLANAR, CYLINDRJCAL AND SPHERICAL BLAST WAVES
    #  H. S. I. Sadek and J. J. Gott1ieb
    #  UTIAS Technica1 Note No. 244, CN ISSN 0082-5263
    #  October, 1983
    #  See eq 3.29, p24 of pdf

    # Solve
    #      det * dpdr|t + (D-u)D dpdr|s + rho*a**2 * D * dudr|s +
    #              N*rho*a**2 * up * (D-u)/r = 0
    # where
    #       dudr|s = dpdr|s * dudp_s
    # So
    #
    #      det * dpdr|t + [ (D-u)D + rho*a**2 * D *dudp_s ]*dpdr|s +
    #              N*rho*a**2 * up * (D-u)/r = 0
    # Or
    #      det * dpdr|t + C1 * dpdr|s + C2 = 0
    # or
    #      dpdr|s = [- det * dpdr|t + C2 ] / C1

    # And in general, with different values for sign, C1 and C2:
    #
    #      dpdr|s = [ $sign * det * slope + C2 ] / C1
    #
    my ( $dpdr_t, $dudr_t, $dpdt_r, $dudt_r ) = ( 0, 0, 0, 0 );

    ###########################################################
    # FIXME: expand for small ovp
    ###########################################################
    # Roundoff control for det...
    # sspd_sq = $gamma*$pabs/$rho = $gamma*($p_amb+$ovp)/$rho 
    #    = $sspd_amb_sq+$gamma*$ovp/rho;
    # sspd_sq - sspd_amb_sq = gamma*$ovp/$rho
    # (sspd-sspd_amb)*(sspd+sspd_amb)=gamma*ovp/rho
    # ... so ..
    # sspdx == sspd-sspd_amb=gamma*ovp/rho/(sspd+sspd_amb)
    ##my $aa_minus_c0 = $gamma * $ovp / $rho / ( $aa + $sspd_amb );
    my $aa_minus_c0 = $aa - $sspd_amb;
    ###########################################################

    my $D_minus_u = $D - $up;
    #  $det = $sspd_sq - $D_minus_u**2;
    my $det = ( $aa + $D_minus_u ) * ( $aa_minus_c0 - $D_minus_c0 + $up );
    if ( $det != 0 ) {

        my $D_minus_u = $D - $up;
        my $dpdt_sf   = $D * $dYdX * $ovp / $lambda;
        my $dudt_sf   = $dudp_sf * $dpdt_sf;
        $dpdr_t =
          -( $D_minus_u * $dpdt_sf +
              $rhocsq * $dudt_sf +
              $symmetry * $rhocsq * $up * $D_minus_u / $lambda ) /
          $det;

        $dudr_t =
          -( $D_minus_u * $dudt_sf +
              $dpdt_sf / $rho +
              $symmetry * $sspd_sq * $up / $lambda ) /
          $det;

        $dpdt_r =
          ( ( $sspd_sq + $up * $D_minus_u ) * $dpdt_sf +
              $rhocsq * $D * $dudt_sf +
              $symmetry * $rhocsq * $up * $D * $D_minus_u / $lambda ) /
          $det;

        $dudt_r =
          ( ( $sspd_sq + $up * $D_minus_u ) * $dudt_sf +
              $D / $rho * $dpdt_sf +
              $symmetry * $sspd_sq * $up * $D / $lambda ) /
          $det;

    }

    my $rhash = {
        dpdr_t        => $dpdr_t,
        dudr_t        => $dudr_t,
        dpdt_r        => $dpdt_r,
        dudt_r        => $dudt_r,
        up            => $up,
        density_ratio => $rhorat,
        ushock        => $D,
        sspd          => $aa,
	sigma         => $sigma,
	Sigma         => $Sigma,
    };

    return $rhash;
}

1;
