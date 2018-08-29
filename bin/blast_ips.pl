#!/usr/bin/perl 
use warnings;
use strict;

# This is a driver to illustrate usage of Blast::IPS.
use Blast::IPS;

# TODO:
# - implement the ground plane

my $audit_string = "";

my %symmetry_name = (
    0 => 'Plane',
    1 => 'Cylindrical',
    2 => 'Spherical',
);

# Setup a default blast table
my $gamma       = 1.4;
my $symmetry    = 2;
my $blast_table = Blast::IPS->new( symmetry => $symmetry, gamma => $gamma );

# Main Loop
my $selection = queryu(<<EOM);
Please select a starting mode (you can switch anytime):
D  - Dimensionless
SI - SI units
EOM
$selection = 'SI' unless $selection eq 'D';
while (1) {
    if ( $selection eq 'D' ) {
        ( $blast_table, $selection ) = eval_dimensionless( $blast_table );
    }
    elsif ( $selection eq 'SI' ) {
        ( $blast_table, $selection) = eval_si( $blast_table );
    }
    else { last; }
}

sub ask_ground_plane {
    my $ground_plane;
    if ( ifyes("Is this explosion on a hard half-plane?[Y/N]") ) {
        $ground_plane = "YES";
    }
    else {
        $ground_plane = "NO";
    }
    return $ground_plane;
}

{

    my $medium_stp;
    my $medium;

    BEGIN {
        my $gamma          = 1.4;
        my $symmetry       = 2;
        my $p_sea_level    = 1.01325e5;
        my $sspd_sea_level = 340.43;
        my $E0             = 1;
        my $medium_stp = {
            _gamma        => $gamma,
            _p_amb        => $p_sea_level,
            _sspd_amb     => $sspd_sea_level,
            _symmetry     => $symmetry,
            _rho_amb      => $gamma * $p_sea_level / $sspd_sea_level**2,
            _E0           => undef,
            _ground_plane => 'YES',
        };
        $medium = $medium_stp;
    }

    sub eval_si {
        my ($blast_table) = @_;
        my $symmetry      = $medium->{_symmetry};
        my $p_amb         = $medium->{_p_amb};
        my $sspd_amb      = $medium->{_sspd_amb};
        my $gamma         = $medium->{_gamma};
        my $ground_plane  = $medium->{_ground_plane}; 
        my $units         = 'SI';

        my $symmetry_2 = $blast_table->get_symmetry();
        my $gamma_2    = $blast_table->get_gamma();

        if ( $symmetry != $symmetry_2 || $gamma != $gamma_2 ) {
            $blast_table =
              Blast::IPS->new( symmetry => $symmetry, gamma => $gamma );
            query(
"Switching from symmetry=$symmetry_2, gamma=$gamma_2 to symmetry=$symmetry, gamma=$gamma"
            );
        }

        my $return_selection = 'D';
        while (1) {
            my $table_name = $blast_table->get_table_name();
            $symmetry = $blast_table->get_symmetry();
            $gamma    = $blast_table->get_gamma();
            print <<EOM;
====================
Main Menu - SI units
====================
First set these geometry and atmosphere parameters...
  S     - Symmetry: $symmetry_name{$symmetry}
  G     - Ground Plane?: $ground_plane
  A     - Atmosphere: gamma=$gamma, P0=$p_amb Pa, sspd=$sspd_amb m/s

Then evaluate the model with these settings...
  P     - do Point evaluations ...
  T     - do Table operations ...
  I     - show Global Information about this blast

Use one of these to go back...
  D     - switch to Dimensionless mode 
  Q     - Quit the program
EOM
            my $ans = queryu(":");
            if ( $ans eq 'S' ) {
                ( $blast_table, $medium ) =
                  select_geometry( $blast_table, $medium );
            }
            elsif ( $ans eq 'G' ) {
		$ground_plane = ask_ground_plane();
		$medium->{_ground_plane}=$ground_plane;
            }
            elsif ( $ans eq 'A' ) {
                ( $blast_table, $medium ) =
                  select_atmosphere_SI( $blast_table, $medium );
            }
            elsif ( $ans eq 'I' ) {
                show_summary_information( $blast_table, $medium, $units );
            }
            elsif ( $ans eq 'D' ) {
                $return_selection = 'D';
                last;
            }
            elsif ( $ans eq 'P' ) {
                point_evaluations_SI( $blast_table, $medium );
            }
            elsif ( $ans eq 'T' ) {
                table_operations( $blast_table, $medium );
            }
            elsif ( $ans eq 'Q' ) {
                if ( ifyes("Quit blast_ips? [Y/N]") ) {
                    $return_selection = 'Q';
                    last;
                }
            }
        }
        return ( $blast_table, $return_selection);
    }
}

sub eval_dimensionless {

    my ( $blast_table, $ground_plane ) = @_;

    # FIXME: eliminate $units, $medium

    my $symmetry = $blast_table->get_symmetry();
    my $gamma    = $blast_table->get_gamma();
    my $units    = 'D';
    my $medium   = {
        _gamma        => $gamma,
        _sspd_amb     => 1,
        _p_amb        => 1,
        _symmetry     => $symmetry,
        _rho_amb      => $gamma,
        _E0           => 1,
        _ground_plane => $ground_plane
    };

    # main loop
    my $return_selection;
    while (1) {
        my $table_name = $blast_table->get_table_name();
        $symmetry = $blast_table->get_symmetry();
        $gamma    = $blast_table->get_gamma();
        print <<EOM;
=========================
Main Menu - Dimensionless
=========================
First set the geometry and atmosphere parameters...
  A     - symmetry: $symmetry_name{$symmetry}; gamma=$gamma

Then evaluate the model with these settings...
  P     - do Point evaluations ...
  T     - do Table operations ...
  I     - show Global Information about this blast

Use one of these to go back...
  SI    - switch to SI Units
  Q     - Quit the program
EOM
        my $ans = queryu(":");
        if ( $ans eq 'A' ) {
            ( $blast_table, $medium ) =
              select_blast_table( $blast_table, $medium );
        }
##?        elsif ( $ans eq 'G' ) {
##?            $ground_plane = ask_ground_plane();
##?            $medium->{_ground_plane} = $ground_plane;
##?        }
        elsif ( $ans eq 'I' ) {
            show_summary_information( $blast_table, $medium, $units );
        }
        elsif ( $ans eq 'SI' ) {
            $return_selection = 'SI';
            last;
        }
        elsif ( $ans eq 'P' ) {
            my $vname = 'X';
            $vname = select_dimensionless_variable($vname);
            point_evaluations_dimensionless( $blast_table, $medium, $vname );
        }
        elsif ( $ans eq 'T' ) {
            table_operations( $blast_table, $medium );
        }
        elsif ( $ans eq 'Q' ) {
            if ( ifyes("Quit blast_ips? [Y/N]") ) {
                $return_selection = 'Q';
                last;
            }
        }
    }
    return ( $blast_table, $return_selection );
}

#sub get_medium {
#    my $patm_psi = 14.7;
#    my $patm_pa  = 1.e5 * $patm_psi / 14.5;
#    my $sspd     = 345;
#}
sub select_geometry {
    my ( $blast_table, $medium ) = @_;
    my $symmetry_now = $blast_table->get_symmetry();
    my $symmetry_name_now =
      ( $symmetry_now == 2 ? 'S' : $symmetry_now == 1 ? 'C' : 'P' );
    my $symmetry_name =
      queryu("Enter symmetry: 'S', 'C' or 'P', <cr>='$symmetry_name_now'");
    if ( !$symmetry_name ) { $symmetry_name = $symmetry_name_now }
    my $symmetry =
      ( $symmetry_name eq 'S' ? 2 : $symmetry_name eq 'C' ? 1 : 0 );
    if ( $symmetry != $symmetry_now ) {
        my $blast_table_old = $blast_table;
        $blast_table =
          Blast::IPS->new( 'symmetry' => $symmetry, 'gamma' => $gamma );
        my $err = $blast_table->get_error();
        if ($err) {
            query("Error: $err; no changes made");
            return ( $blast_table_old, $medium );
        }
        $medium->{_symmetry} = $symmetry;
    }
    return ( $blast_table, $medium );
}

{
    my $alt_m;

    BEGIN {
        $alt_m = 0;
    }

    sub select_atmosphere_SI {
        my ( $blast_table, $medium ) = @_;

        my $gamma           = $medium->{_gamma};
        my $p_amb           = $medium->{_p_amb};
        my $sspd_amb        = $medium->{_sspd_amb};
        my $symmetry        = $medium->{_symmetry};
        my $blast_table_old = $blast_table;
        my $touched;
        while (1) {
            print <<EOM;
==========================
Set atmospheric conditions
==========================
AL  - Set earth altitude; current value............: $alt_m m 
G   - change gamma; current value..................: $gamma 
P   - change ambient pressure; currrent value......: $p_amb Pa
C   - change ambient sound speed; current value....: $sspd_amb m/s
Y   - Yes, return with these values
Q   - Quit, keep original values
EOM
            my $ans = queryu('-->');
            if ( $ans ne 'Y' && $ans ne 'Q' ) { $touched = 1 }
            if ( $ans eq 'AL' ) {
                $gamma = 1.4;
                ( $p_amb, $sspd_amb, $alt_m ) =
                  altitude( $p_amb, $sspd_amb, $alt_m );
            }
            elsif ( $ans eq 'G' ) {
                $gamma = get_num( "Enter gamma", 1.4 );
                $alt_m = "?";
            }
            elsif ( $ans eq 'P' ) {
                $p_amb = get_num( "Enter ambient pressure, Pa", $p_amb );
                $alt_m = "?";
            }
            elsif ( $ans eq 'C' ) {
                $sspd_amb =
                  get_num( "Enter ambient sound speed, m/s", $sspd_amb );
                $alt_m = "?";
            }
            elsif ( $ans eq 'Y' ) {
                my $gamma_old    = $blast_table->get_gamma();
                my $symmetry_old = $blast_table->get_symmetry();
                if ( $gamma != $gamma_old || $symmetry != $symmetry_old ) {
                    my $blast_table_old = $blast_table;
                    $blast_table = Blast::IPS->new(
                        'symmetry' => $symmetry,
                        'gamma'    => $gamma
                    );
                    my $err = $blast_table->get_error();
                    if ($err) {
                        query("Error: $err; no changes made");
                        $blast_table = $blast_table_old;
                        last;
                    }
                }
                $medium->{_gamma}    = $gamma;
                $medium->{_p_amb}    = $p_amb;
                $medium->{_sspd_amb} = $sspd_amb;
                my $rho_amb = $gamma * $p_amb / $sspd_amb**2;
                $medium->{_rho_amb} = $rho_amb;
                last;
            }
            elsif ( $ans eq 'Q' ) {
                last if ( !$touched );
                last if ( ifyes("Return without keeping changes?[Y/N}") );
            }
        }
        return ( $blast_table, $medium );
    }
}

sub select_blast_table {
    my ( $blast_table, $medium ) = @_;
    my $symmetry_name = queryu("Enter symmetry: 'S', 'C' or 'P', <cr>='S'");
    if ( !$symmetry_name ) { $symmetry_name = 'S' }
    my $gamma = get_num( "Enter gamma", 1.4 );
    my $blast_table_old = $blast_table;

    my $symmetry =
      ( $symmetry_name eq 'S' ? 2 : $symmetry_name eq 'C' ? 1 : 0 );
    $blast_table =
      Blast::IPS->new( 'symmetry' => $symmetry, 'gamma' => $gamma );
    my $err = $blast_table->get_error();
    if ($err) {
        query("Error: $err; no changes made");
        return ( $blast_table_old, $medium );
    }
    my $p_amb    = $medium->{_p_amb};
    my $sspd_amb = $medium->{_sspd_amb};
    my $rho_amb  = $gamma * $p_amb / $sspd_amb**2;
    $medium->{_gamma}    = $gamma;
    $medium->{_symmetry} = $symmetry;
    $medium->{_rho_amb}  = $rho_amb;
    return ( $blast_table, $medium );
}

sub show_summary_information {
    my ( $blast_table, $medium, $units ) = @_;
    my $rinfo                   = $blast_table->get_info();
    my $table_name              = $rinfo->{table_name};
    my $alpha                   = $rinfo->{alpha};
    my $symmetry                = $rinfo->{symmetry};
    my $gamma                   = $rinfo->{gamma};
    my $Sint_pos                = $rinfo->{Sintegral_pos};
    my $Sint_neg                = $rinfo->{Sintegral_neg};
    my $Ixr_pos                 = $gamma * $Sint_pos;
    my $Ixr_neg                 = $gamma * $Sint_neg;
    my $r_tail_shock            = $rinfo->{r_tail_shock};
    my $z_tail_shock            = $rinfo->{z_tail_shock};
    my $r_tail_formed           = $rinfo->{r_tail_formed};
    my $z_tail_formed           = $rinfo->{z_tail_formed};
    my $KE_initial_ratio        = $rinfo->{KE_initial_ratio};
    my $pcenter_initial_ratio   = $rinfo->{pcenter_initial_ratio};
    my $t_u_plus_c_zero         = $rinfo->{t_u_plus_c_zero};
    my $t_pcenter_zero_1        = $rinfo->{t_pcenter_zero_1};
    my $t_pcenter_min           = $rinfo->{t_pcenter_min};
    my $pcenter_min             = $rinfo->{pcenter_min};
    my $t_dpdt_center_max       = $rinfo->{t_dpdt_center_max};
    my $dpdt_center_max         = $rinfo->{dpdt_center_max};
    my $pcenter_dpdt_center_max = $rinfo->{pcenter_dpdt_center_max};
    my %symmetry_name           = (
        0 => 'Plane',
        1 => 'Cylindrical',
        2 => 'Spherical',
    );

    foreach ( $alpha, $Sint_pos, $Sint_neg, $Ixr_pos, $Ixr_neg ) {
        $_ = sprintf( "%0.6g", $_ );
    }
    my $pstr =
      ( $symmetry == 2 ) ? "times r" : ( $symmetry == 1 ) ? "times r^1/2" : "";
    print <<EOM;

=====================
Blast Wave Parameters
=====================
symmetry  = $symmetry_name{$symmetry}
gamma     = $gamma
alpha     = $alpha = similarity solution parameter alpha
KE/E0     = $KE_initial_ratio = initial ratio of KE to total energy
t u+c=0   = $t_u_plus_c_zero = time of first zero velocity C+ characteristic
I+        = $Ixr_pos = positive impulse $pstr at long range
I-        = $Ixr_neg = negative impulse $pstr at long range
S+        = $Sint_pos = positive integral of Sigma x dr at long range
S-        = $Sint_neg = negative integral of Sigma x dr at long range
r TS      = $r_tail_shock = scaled radius at which Tail Shock first forms
z TS      = $z_tail_shock = scaled z at which Tail Shock first forms 
r TF      = $r_tail_formed = scaled radius at which Tail Shock becomes peak
z TF      = $z_tail_formed = scaled z at which Tail Shock becomes peak

Overpressure at r=0:
pc/ps0      = $pcenter_initial_ratio = initial ratio of overpressure to shock overpressure 
t_pc0       = $t_pcenter_zero_1 = time of first zero overpressure at r=0
t_pc_min    = $t_pcenter_min   = time of minimum overpressure at r=0
pc_min      = $pcenter_min   = minimum overpressure at r=0
t_dpdt_max  = $t_dpdt_center_max = time of maximum dp/dt at r=0
dpdt_max    = $dpdt_center_max   = maximum dp/dt at r=0
p_dpdt_max  = $pcenter_dpdt_center_max = ovp at time of maximum dp/dt at r=0

Note: zeros indicate undefined values.
EOM

    # DEBUG dump to check key names
    if (0) {
        print "\nKey dump\n";
        foreach my $key ( keys %{$rinfo} ) {
            print "$key => $rinfo->{$key}\n";
        }
    }

    hitcr();
    return;
}

{
    my $vname;

    BEGIN {
        $vname = 'x';
    }

    sub select_dimensionless_variable {

        # key => [ order, text ]
        my $i    = 0;
        my %menu = (
            'x'     => [ ++$i, 'scaled range, = r/d' ],
            'y'     => [ ++$i, 'overpressure ratio, =(P-P0)/P0' ],
            't'     => [ ++$i, 'scaled time of arrival, = c0*TOA/d' ],
            'z'     => [ ++$i, 'x - t' ],
            'X'     => [ ++$i, 'ln(x)' ],
            'Y'     => [ ++$i, 'ln(y)' ],
            'T'     => [ ++$i, 'ln(t)' ],
            'Z'     => [ ++$i, 'ln(z)' ],
            'Z-X'   => [ ++$i, 'ln(z/x)' ],
            'E1'    => [ ++$i, 'residual energy from primary shock' ],
            'E'     => [ ++$i, 'residual energy (total)' ],
            'dYdX'  => [ ++$i, 'dY/dX' ],
            'dZdX'  => [ ++$i, 'dZ/dX' ],
            'dWdX'  => [ ++$i, 'dW/dX' ],
            'dE1dX' => [ ++$i, 'dE1/dX' ],
            'M' => [ ++$i, '= (S/c0), = shock Mach number (S=shock speed)' ],
            'q' => [ ++$i, '= 1/M^2' ],
        );

        my $menu_text = <<EOM;
========================================
VARIABLE SELECTION MENU (Dimensionless): 
========================================
Let
 d = scaled distance (E/P0)^(1/N) 
 E is energy and N = 1,2, or 3 is symmetry
 P0, c0 = ambient atmospheric pressure, sound speed 
You may vary One of these variables: 
EOM
        foreach
          my $key ( sort { $menu{$a}->[0] <=> $menu{$b}->[0] } keys(%menu) )
        {
            $menu_text .= "    $key : $menu{$key}->[1]\n";
        }

        # Default; could also be previous value
        my $default = 'x';
        while (1) {
            print $menu_text;
            my $ans = query("Select one of these variables; <cr>='$default':");
            $ans = $default unless ($ans);
            if ( defined( $menu{$ans} ) ) {
                $vname = $ans;
                return $vname;
            }
            else {
                hitcr("error, try again");
            }
        }
        return $vname;
    }
}

sub ground_factor {
    my ( $symmetry, $ground_plane ) = @_;

    # by what factor should we scale the user's energy to
    # get the energy for use in the analytical model?
    my $escale = 1;
    if ( $ground_plane eq 'YES' ) {

        # in spherical and cylindrical symmetry with energy on the ground,
        # double the energy going from user to model.  Plane symmetry
        # already assumes a ground plane.
        if ( $symmetry > 0 ) { $escale = 2 }
    }
    elsif ( $ground_plane eq 'NO' ) {

        # in plane symmetry, halve the energy if no
        # going from user to model
        if ( $symmetry == 0 ) { $escale = 0.5 }
    }
    else {
	die "Ground plane is '$ground_plane'; expected to be YES or NO\n";
    }
    return $escale;
}

{

    my $range;
    BEGIN {
        $range        = undef;
    }

    sub point_evaluations_SI {
        my ( $blast_table, $medium ) = @_;

        # perform point evaluations with units
        # internal units are SI but other display units may be used
        my $gamma        = $medium->{_gamma};
        my $symmetry     = $medium->{_symmetry};
        my $p_amb        = $medium->{_p_amb};
        my $sspd_amb     = $medium->{_sspd_amb};
        my $ground_plane = $medium->{_ground_plane};
        my $E0           = $medium->{_E0};
	my $dscale;
        if ( $p_amb <= 0 || $sspd_amb <= 0 ) {
            query("Enter positive atmospheric pressure and sspd first");
            return;
        }
        my $ground_factor = ground_factor( $symmetry, $ground_plane );

        my $set_dscale = sub {

            # Set the scaling distance whenever E0 changes
            return unless defined($E0);
            $dscale =
              ( $ground_factor * $E0 / $p_amb )**( 1 / ( $symmetry + 1 ) );
        };
        my $ask_for_E0 = sub {

            # Get energy in Joules, but must be >0
            # Reset distance scale
            my $test = get_num( "Enter energy E0, Joules:", $E0 );
            if ( defined($test) && $test > 0 ) {
                $E0 = $test;
                $set_dscale->();
            }
        };
        my $ask_for_range = sub {

            # Get range in m, but must be >0
            my $test = get_num( "Enter range, m:", $range );
            if ( defined($test) && $test > 0 ) { $range = $test }
        };
        my $Nsym      = $symmetry + 1;    # = (1,2 or 3)
        my $to_string = sub {
            my ($val) = @_;
            return ( defined($val) ? $val : '?' );
        };
        my $dscale_to_string = sub {
            my ($dscale) = @_;
            my $str = '?';
            if ( defined($dscale) ) {
		$dscale = sprintf("%0.5g", $dscale);
                $str = "$dscale";
            }
            return $str;
        };
        my $range_to_string = sub {
            my ($range) = @_;
            my $str = '?';
            if ( defined($range) ) {
                my $range_ft = $range * 3.2808;
		$range = sprintf("%0.5g", $range);
		$range_ft = sprintf("%0.5g", $range_ft);
                $str = "$range  m = $range_ft ft";
            }
            return $str;
        };
        my $E0_to_string = sub {
            my ($E0) = @_;
            my $str = '?';
            if ( defined($E0) ) {
                my $E0_kt = $E0 / 4.184e12;
                my $E0_mj = $E0 / 1.e6;
                #if ($E0_kt > 0.1) {$E0_kt = sprintf( "%0.4f", $E0_kt );}
                #if ( $E0_mj > 0.1 ) { $E0_mj = sprintf( "%0.4f", $E0_mj ); }
		$E0_kt = sprintf("%0.5g", $E0_kt);
		$E0_mj = sprintf("%0.5g", $E0_mj);
                $str = "$E0_mj MJ = $E0_kt kt";
            }
            return $str;
        };

        $set_dscale->();
	my $pstr="1/".($symmetry+1);

        while (1) {
            my $range_str  = $range_to_string->($range);
            my $dscale_str = $dscale_to_string->($dscale);
            my $E0_str     = $E0_to_string->($E0);
            #my $gtext      = $ground_plane ? "on hard surface" : "in free air";

            print <<EOM;
 ----- Dimensional Solution Menu -------
     gamma....................... $gamma
     symmetry.................... $symmetry
     explosion on a ground plane? $ground_plane (g = $ground_factor)
     atmospheric Pressure, Pa.... $p_amb 
     ambient sound speed, m/s.... $sspd_amb
  E  Energy...................... $E0_str   
  R  Range....................... $range_str
     scale dist. (g*E0/P0)^$pstr... $dscale_str m 
   
  RE or C Calculate blast parameters, given: R, E
  RI calculate Energy, given: Range, Impulse               [TBD]
  RP calculate Energy, given: Range, measured Overpressure
  RT calculate Energy, given: Range, measured TOA
  PI calculate Energy, given: OVP, Impulse                 [TBD]
  EP calculate Range, given Energy and OVP
  EI calculate Range, given Energy and IMP                 [TBD]
  ET calculate Range, given Energy and TOA   
  
  Z Zoom  - view latest results, 1 screen per channel
  L List  - view latest results, 1 line per channel
  F files (read/write)...
  Q=Quit   CL=Clear  LG=List Gages LD=List Data

EOM
            my $ans = queryu(":");
            if ( $ans eq 'E' ) {
                $ask_for_E0->();
            }
            elsif ( $ans eq 'R' ) {
                $ask_for_range->();
            }
            elsif ( $ans eq 'G' ) {
                print <<EOM;
Select a ground plane option:
  0 = explosion is in free air
  1 = explosion is on a ground plane (rigid surface)
EOM
                $ground_plane = query(":");
            }
            elsif ( $ans eq 'AA' ) {
		query("Not programmed yet");
		next
            }
            elsif ( $ans eq 'RE' || $ans eq 'ER' || $ans eq 'C' ) {
                if ( !defined($range) ) { $ask_for_range->(); }
                if ( !defined($E0) )    { $ask_for_E0->(); }
                my $x    = $range / $dscale;
                my $X    = log($x);
                my $ret  = $blast_table->wavefront( 'X' => $X );
                my $Y    = $ret->{Y};
                my $y    = exp($Y);
                my $dYdX = $ret->{dYdX};
                query("ovp ratio=$y; dYdX=$dYdX");
            }
            elsif ( $ans eq 'RI' || $ans eq 'IR' ) {
		query("Not programmed yet");
		next
            }
            elsif ( $ans eq 'RT' || $ans eq 'TR' ) {
                if ( !defined($range) ) { $ask_for_range->(); }
                my $t = get_num("Enter toa, s:");
                next if ( $t <= 0 );
                my $z = ( $range - $t * $sspd_amb );
                if ( $z <= 0 ) {
                    my $toa_min = $range / $sspd_amb;
                    query("toa must exceed $toa_min at range $range");
                    next;
                }
                my $ff   = log( $z / $range );
                my $ret  = $blast_table->wavefront( 'Z-X' => $ff );
                my $X    = $ret->{X};
                my $Y    = $ret->{Y};
                my $T    = $ret->{T};
                my $dYdX = $ret->{dYdX};
                my $x    = exp($X);
                $dscale = $range / $x;
                my $ttest = $dscale * exp($T);
                $E0 = $p_amb * ($dscale)**$Nsym / $ground_factor;
                query("Energy is E0=$E0, dYdX=$dYdX, ttest=$ttest =? $t");
            }
            elsif ( $ans eq 'PI' || $ans eq 'IP' ) {
		query("Not programmed yet");
		next
            }
            elsif ( $ans eq 'RP' || $ans eq 'PR' ) {
                if ( !defined($range) ) { $ask_for_range->(); }
                my $y = get_num("Enter incident overpressure ratio:");
                next if ( $y <= 0 );
                my $Y    = log($y);
                my $ret  = $blast_table->wavefront( 'Y' => $Y );
                my $X    = $ret->{X};
                my $dYdX = $ret->{dYdX};
                my $x    = exp($X);
                $dscale = $range / $x;
                $E0 = $p_amb * ($dscale)**$Nsym;
                query("E0=$E0 gives ovp ratio=$y at range=$range, dYdX=$dYdX");
            }
            elsif ( $ans eq 'EP' || $ans eq 'PE' ) {
                if ( !defined($E0) ) { $ask_for_E0->(); }
                my $y = get_num("Enter incident overpressure ratio:");
                next if ( $y <= 0 );
                my $Y    = log($y);
                my $ret  = $blast_table->wavefront( 'Y' => $Y );
                my $X    = $ret->{X};
                my $dYdX = $ret->{dYdX};
                my $x    = exp($X);
                $range = $x * $dscale;
                query("range=$range gives ovp ratio=$y; dYdX=$dYdX");
            }
            elsif ( $ans eq 'ET' || $ans eq 'TE' ) {
                if ( !defined($E0) ) { $ask_for_E0->(); }
                my $t = get_num("Enter toa, s:");
                next if ( $t <= 0 );
                my $T    = log( $t * $sspd_amb / $dscale );
                my $ret  = $blast_table->wavefront( 'T' => $T );
                my $X    = $ret->{X};
                my $Y    = $ret->{Y};
                my $dYdX = $ret->{dYdX};
                my $y    = exp($Y);
                my $x    = exp($X);
                $range = $x * $dscale;
                query("t=$t gives range=$range, ovp ratio=$y; dYdX=$dYdX");
            }
            elsif ( $ans eq 'EI' || $ans eq 'IE' ) {
		query("Not programmed yet");
		next
                #if ( !defined($E0) ) { $ask_for_E0->(); }
            }
            elsif ( $ans eq 'Q' ) {
                last;
            }
            else {
            }
        }
	$medium->{$E0} = $E0;
        return;
    }
}

sub point_evaluations_dimensionless {
    my ( $blast_table, $medium, $vname ) = @_;
    my $gamma    = $medium->{_gamma};
    my $symmetry = $medium->{_symmetry};
    my $p_amb    = $medium->{_p_amb};
    my $pi       = 4 * atan2( 1, 1 );
    my $NN       = $symmetry + 1;
    my $another  = 'a';

    # Note that a ground_plane and ground factor are not currently used
    # in dimensionless mode. They have to be applied before converting
    # to an actual geometry.
    # my $ground_factor = ground_factor( $symmetry, $ground_plane );

    while (1) {
        my $val =
          get_num("Enter $another value for '$vname', or <cr> to quit:");
        $another = 'another';
        last if $val eq "" || $val !~ /\d/;    #^\s*[\-\+\d\.]/;

        my ( $iQ, $Q );
        if ( $vname =~ /^([XYZT]|E1|Er|dYdX|dZdX|dTdX|dE1dX|dErdX)/ ) {
            $Q  = $val;
            $iQ = $vname;
        }
        elsif ( $vname =~ /^[xyzt]$/ ) { $Q = log($val); $iQ = uc($vname) }
        elsif ( $vname eq 'q' ) {

            # Convert q=(c0/D)**2 to Y=ln(ovp ratio)
            $iQ = 'Y';
            my $q = $val;
            if ( $q > 0 && $q < 1 ) {
                my $ovprat =
                  ( 2 * $gamma - ( $gamma - 1 ) * $q ) /
                  ( ( $gamma + 1 ) * $q ) - 1;

                if ( $ovprat <= 0 ) {
                    hitcr("cant take log of '$ovprat'");
                    next;
                }
                $Q = log($ovprat);
            }
            else {
                hitcr("must have 0< q < 1");
                next;
            }
        }
        elsif ( $vname eq 'm' ) {

            # Convert m=(S/c0) to Y=ln(ovp ratio)
            $iQ = 'Y';
            my $m = $val;
            if ( $m > 1 ) {
                my $q = 1 / $m**2;
                my $ovprat =
                  ( 2 * $gamma - ( $gamma - 1 ) * $q ) /
                  ( ( $gamma + 1 ) * $q ) - 1;

                if ( $ovprat <= 0 ) {
                    hitcr("cant take log of '$ovprat'");
                    next;
                }

                $Q = log($ovprat);
            }
            else {
                hitcr("m must be > 1");
                next;
            }
        }
        else {
            die "coding incomplete for variable '$vname'";
        }
        my $ret  = $blast_table->wavefront( $iQ => $Q );
        my $X    = $ret->{X};
        my $Y    = $ret->{Y};
        my $Z    = $ret->{Z};
        my $dYdX = $ret->{dYdX};
        my $dZdX = $ret->{dZdX};

        #my $Lpos      = $ret->{Lpos};
        #my $Tpos      = $ret->{Tpos};
        #my $Lneg      = $ret->{Lneg};
        my $Ixr_pos     = $ret->{Ixr_pos};
        my $Ixr_neg     = $ret->{Ixr_neg};
        my $E1          = $ret->{E1};
        my $Er          = $ret->{Er};
        my $W_atm       = $ret->{W_atm};
        my $W_blast     = $ret->{W_blast};
        my $dErdX       = $ret->{dErdX};
        my $dE1dX       = $ret->{dE1dX};
        my $dEdE1       = $dE1dX && $dErdX ? $dErdX / $dE1dX : 1;
        my $ke_pos      = $ret->{KE_pos};
        my $work_pos    = $ret->{W_pos};
        my $qint_pos    = $ret->{qint_pos};
        my $z_pose_rs   = $ret->{z_pose_rs};
        my $z_nege_rs   = $ret->{z_nege_rs};
        my $z_pmin_rs   = $ret->{z_pmin_rs};
        my $rovp_min_rs = $ret->{rovp_min_rs};
        my $disp_pos    = $ret->{disp_pos};
        my $TableLoc    = $ret->{TableLoc};
        my $z_pose_ts   = $ret->{z_pose_ts};
        my $z_nege_ts   = $ret->{z_nege_ts};
        my $z_pmin_ts   = $ret->{z_pmin_ts};

        my ( $il, $iu, $ntab ) = @{$TableLoc};
        my $table_location =
            ( $il < 0 )     ? "Before Table Start"
          : ( $iu > $ntab ) ? "Extrapolation Beyond Table End"
          :                   "Table Interpolation";

        # To be deleted:
        my $Ixr_pos_lim = $ret->{Ixr_pos_lim};
        my $Ixr_neg_lim = $ret->{Ixr_neg_lim};

        my $E0 = 1;           ## For future use
        my $E2 = $Er - $E1;

        my $Tpmin = $z_pmin_rs - $z_nege_rs;

        my $x    = exp($X);
        my $y    = exp($Y);
        my $z    = exp($Z);
        my $t    = $x - $z;
        my $T    = $t > 0 ? log($t) : -999;
        my $term = $y * ( $gamma + 1 ) / ( 2 * $gamma );
        my $m    = sqrt( 1 + $term );
        my $q    = 1 / $m**2;
        my $up   = $y / ( $gamma * $m );

        my $dV_atm   = $W_atm / $p_amb;
        my $rbar     = $x;
        my $disp_end = displacement( $dV_atm, $x, $symmetry );

        my $dt_pose_rs = $z - $z_pose_rs;
        my $dt_pmin_rs = $z - $z_pmin_rs;
        my $dt_nege_rs = $z - $z_nege_rs;

        my $dr_pose_ts = $z - $z_pose_ts;
        my $dr_pmin_ts = $z - $z_pmin_ts;
        my $dr_nege_ts = $z - $z_nege_ts;

        # derived variables
        my ( $density_ratio_shock, $density_ratio_equilibrium ) =
          density_ratios( $Y, $gamma );

        foreach ( $x, $y, $z, $t, $X, $Y, $Z, $T, $dYdX, $dZdX ) {
            $_ = sprintf( "%0.8g", $_ );
        }

        #$Lpos,     $Lneg,
        #$Lpos,                $Tpos,
        foreach (
            $m,                   $q,
            $up,                  $Ixr_pos,
            $Ixr_neg,             $E1,
            $W_atm,               $Er,
            $E2,                  $W_blast,
            $dEdE1,               $ke_pos,
            $work_pos,            $dt_pose_rs,
            $dt_pmin_rs,          $dt_nege_rs,
            $dr_pose_ts,          $dr_pmin_ts,
            $dr_nege_ts,          $rovp_min_rs,
            $disp_pos,            $disp_end,
            $density_ratio_shock, $density_ratio_equilibrium,
          )
        {
            $_ = sprintf( "%0.6g", $_ );
        }

        # preparing for future version with units
        my $t_unit = "(scaled)";
        my $d_unit = "(scaled)";
        my $e_unit = "(scaled)";
        my $u_unit = "(scaled)";

        my $pstr =
          ( $symmetry == 2 ) ? "x r" : ( $symmetry == 1 ) ? "x r^1/2" : "";

#r^n/2 I+ = $Ixr_pos_lim = limiting positive impulse $pstr
#r^n/2 I- = $Ixr_neg_lim = limiting negative impulse $pstr
#L-    = $Lneg = (OLD) length of negative phase $d_unit
#t+   - toa = $Tpos = (OLD) time duration to end of positive overpressure phase at shock radius, $t_unit
#L+    = $Lpos = (OLD) length of positive phase $d_unit
        print <<EOM;
----------------------------------------------------------
Symmetry: $symmetry   Gamma: $gamma Method: $table_location  
x    = $x = scaled range r/d;       ln(x) = X = $X 
y    = $y = (P-P0)/P0;              ln(y) = Y = $Y 
z    = $z = (r-c0 t)/d;             ln(z) = Z = $Z 
toa  = $t = scaled time of arrival  ln(t) = T = $T 
dYdX = $dYdX
t+   - toa = $dt_pose_rs = time duration to end of positive overpressure phase at shock radius, $t_unit
tmin - toa = $dt_pmin_rs = time duration to minimum overpressure at shock radius, $t_unit
t-   - toa = $dt_nege_rs = time duration to end of negative phase at shock radius, $t_unit (spherical only)
L+    = $dr_pose_ts = distance from shock front to end of positive phase $d_unit
Lmin  = $dr_pmin_ts = distance from shock front to minimum overpressure $d_unit
L-    = $dr_nege_ts = distance from shock front to end of negative phase $d_unit
m  = $m = S/c0=shock Mach number;  q = $q = 1/m^2
up = $up = shock particle velocity $u_unit
rho_s/rho0 = $density_ratio_shock = density ratio at shock
rho_end/rho0 = $density_ratio_equilibrium = final density ratio
r^n/2 pmin = $rovp_min_rs = minimum overpressure $pstr
r^n/2 I+ = $Ixr_pos = positive phase overpressure impulse $pstr
r^n/2 I- = $Ixr_neg = negative phase overpressure impulse $pstr
qint+   = $qint_pos = positive phase dynamic impulse
E0      = $E0 = initial total energy $e_unit
E1      = $E1 = residual energy of main shock to this range $e_unit
E2      = $E2 = residual energy of tail shock to this range $e_unit
Er      = $Er  = E1+E2 = total residual energy (main shock+tail shock) to this range $e_unit
W_atm   = $W_atm = (gamma-1)*Er = work of thermal expansion against the atmosphere $e_unit
W_blast = $W_blast = (E0-Er-W_atm) = work of the blast $e_unit
W_pos   = $work_pos = work of positive phase of the blast $e_unit
ke_pos  = $ke_pos = kinetic energy in the positive phase $e_unit
dEr/dE1  = $dEdE1 = energy dissipation ratio (>1 if tail shock)
disp_pos = $disp_pos = maximum particle displacement
disp_end = $disp_end = equilibrium particle displacement

Note: zeros indicate undefined values.
EOM

    }
}

sub table_operations {
    my ( $blast_table, $medium ) = @_;

    # Default run
    while (1) {
        my $table_name = $blast_table->get_table_name();
        print <<EOM;
Current Table Name: $table_name
Make a Table of values.

  PLOT		- print uniformly spaced points in ln(range) for plotting
  CSV 		- Write out the current table as a .csv file
  INFO		- Display Info of the current table
  q or Q        - done
EOM

        my $ans = queryu(":");
        if ( $ans eq 'PLOT' ) {

            my $rsolution = get_table_points($blast_table);

            my $hdr = "X\tY\tdY/dX\tZ\tdZ/dX\tT\tdT/dX";
            ( $rsolution, $hdr ) =
              add_extra_variables( $rsolution, $hdr, $medium );
            my $ofile = get_output_filename( "Output file:", "profile.txt" );
            print_table( $rsolution, $ofile, $audit_string . $hdr );
        }
        elsif ( $ans eq 'INFO' ) {
            my $rtable = $blast_table->get_table();
            my $nrows  = @{$rtable};
            my $ncols  = @{ $rtable->[0] };
            my $Xb     = $rtable->[0]->[0];
            my $Xe     = $rtable->[-1]->[0];
            my $A_far  = $blast_table->{_A_far};
            my $B_far  = $blast_table->{_B_far};
            print <<EOM;
Table has $nrows rows and $ncols cols
$Xb <= ln(r) <= $Xe
A=$A_far, B=$B_far
EOM
            hitcr("");
        }
        elsif ( $ans eq 'CSV' ) {
            my $Id_max = 4;
            if ( ifyes( "Do you want T and dT/dX too? [Y/N], <cr>=N", "N" ) ) {
                $Id_max = 6;
            }
            my $rtable_raw = $blast_table->get_table();
            my $rtable_trim;
            foreach my $item ( @{$rtable_raw} ) {
                push @{$rtable_trim}, [ @{$item}[ 0 .. $Id_max ] ];
            }
            my $types = queryu(
                "Enter any combination of 'T' (.txt), 'P' (.pl), 'L' (.tex)");
            my $basename = query("Enter a ROOT filename (no extension");
            if ( $types =~ /T/ ) {
                print_table( $rtable_trim, $basename . ".txt" );
            }
            if ( $types =~ /P/ ) {
                print_table( $rtable_trim, $basename . ".pl" );
            }
            if ( $types =~ /L/ ) {
                print_table( $rtable_trim, $basename . ".tex" );
            }
        }
        elsif ( $ans eq 'Q' ) {
            last;
        }
    }
}

sub add_extra_variables {
    my ( $rtable, $hdr, $medium ) = @_;
    foreach my $item ( @{$rtable} ) {
        my ( $X, $Y, $dYdX, $Z, $dZX, $T, $dTdx ) = @{$item};
        my $ovprat = exp($Y);
        my $range  = exp($X);
        my ($sigma) = sigma_from_ovprat( $medium, $ovprat );
        my $SS = $sigma * $range;
        my ( $ushock, $upshock, $ushockx, $dup_dovp ) =
          ushock_from_ovprat( $medium, $ovprat );
        push @{$item}, ( $SS, $ushock, $upshock );
    }
    $hdr .= "\tSS\tUs\tUp";
    return ( $rtable, $hdr );
}

sub sigma_from_ovprat {

    # Returns riemann integral from overpressure ratio
    my ( $self, $ovp ) = @_;
    my $gamma    = $self->{_gamma};
    my $sspd_amb = $self->{_sspd_amb};
    my $pow      = ( $gamma - 1 ) / ( 2. * $gamma );
    my $dsspd;
    if ( $ovp < 1.e-4 ) {
        $dsspd =
          $sspd_amb * $ovp * $pow *
          ( 1 + $ovp * ( ( $pow - 1 ) / 2 * ( 1 + $ovp * ( $pow - 2 ) / 3 ) ) );
    }
    else {
        $dsspd = $sspd_amb * ( ( 1 + $ovp )**$pow - 1 );
    }
    my $sigma = ( 2 / ( $gamma - 1 ) ) * $dsspd;
    return ($sigma);
}

sub ushock_from_ovprat {
    my ( $self, $ovprat ) = @_;

    # shock speed from overpressure
    my $gamma    = $self->{_gamma};
    my $sspd_amb = $self->{_sspd_amb};
    my $rho_amb  = $self->{_rho_amb};
    my $p_amb    = $self->{_p_amb};
    my $term     = $ovprat * ( $gamma + 1 ) / ( 2 * $gamma );
    my $mshock   = sqrt( 1 + $term );
    my $ushock   = $sspd_amb * $mshock;
    my $ushockx  = $ushock - $sspd_amb;

    # Taylor series expansion of the square root for small overpressures
    # Switch at a point which keeps error < about 1.e-17
    if ( $term < 1.e-4 ) {
        $ushockx = $term * ( 0.5 + $term * ( -0.125 + $term / 16 ) );
    }

    # shock particle velocity
    my $upshock = $ovprat * $p_amb / ( $rho_amb * $ushock );

    # dup/dovp
    my $dup_dovp =
      $sspd_amb *
      $sspd_amb / $gamma / $ushock *
      ( 1 - ( $gamma + 1 ) / 4 * $upshock / $ushock );
    return ( $ushock, $upshock, $ushockx, $dup_dovp );
}

sub get_table_points {
    my ($blast_obj) = @_;

    # Evaluate the table at a selected list of points

    my $x1 = get_num( "Minimum value of X=ln(r):", -5 );
    my $x2 = get_num( "Maximum value of X=ln(r):", 15 );

# Future option to add, in case second arg to this routine is list of other table points
# 2 use X values from other Table (if appropriate)

    print <<EOM;
0 use Equally Spaced Points in X
1 use X values in Blast::IPS 
EOM
    my $opt = get_num( "Your selection:", 0 );
    my $rsolution;
    if ( $opt eq 0 ) {
        my $npts = get_num( "Total number of points", 400 );
        $npts = 100 unless ($npts);
        $rsolution = $blast_obj->table_gen( $npts, $x1, $x2 );
    }
    else {
        $rsolution = $blast_obj->get_table();

        # include one extra point on each edge
        my ( $j1, $j2 );
        my ( $jl, $ju );
        my $jmax = @{$rsolution} - 1;
        ( $jl, $ju ) = locate_2d( $rsolution, $x1 );
        $j1 = ( $jl > 0 ) ? $jl : 0;
        ( $jl, $ju ) = locate_2d( $rsolution, $x2 );
        $j2 = ( $ju < $jmax ) ? $ju : $jmax;
        if ( $j2 >= $j1 && ( $j1 > 0 || $j2 < $jmax ) ) {
            my @vals = @{$rsolution}[ $j1 .. $j2 ];
            $rsolution = \@vals;
        }

    }
    return $rsolution;
}

sub print_table {

    my $SETTING_string = "";
    my $ERROR_string   = "";
    my ( $rtable, $ofile, $header ) = @_;
    my $format = "";
    if    ( $ofile =~ /\.pl$/ )  { $format = "PERL" }
    elsif ( $ofile =~ /\.tex$/ ) { $format = "TEX" }
    if ( $ofile && open( OUT, ">$ofile" ) ) {
        if ($header) { chomp $header; print OUT "$header\n"; }
        if ( $format eq "PERL" ) {
            print OUT audit_string('#');
            print OUT string_to_comment($SETTING_string) if $SETTING_string;
            print OUT string_to_comment($ERROR_string)   if $ERROR_string;
            print OUT "\n";
            print OUT "\$rpoint_source_table = [\n";
        }
        elsif ( $format eq "TEX" ) {

   #$X$ & $Y$ & $dY/dX$ & $Z$ & $dZ/dX$ \\ \hline
   #$X=\ln(R_s)$ & $Y=\ln(ovp_s)$ & $dY/dX$ & $Z=R_s-c_0T_s$ & $dZ/dX$ \\ \hline
            print OUT <<'EOM';
\begin{center}
\begin{tabular}{ | c | c | c | c | c | } \hline
$X=\ln(R_s)$ & $Y=\ln(P_s/P_0-1)$ & $dY/dX$ & $Z=ln(R_s-c_0T_s)$ & $dZ/dX$ \\ \hline
EOM
        }
        foreach my $item ( @{$rtable} ) {

            #print STDERR "BOOGA: Newline";
            my $str;
            if ( $format eq "PERL" ) {
                $str = '[' . join( ', ', @{$item} ) . "],";
            }
            elsif ( $format eq "TEX" ) {
                $str = "";

                #$str = join( ' & ', @{$item} ) . "\\\\";
                my $count  = 0;
                my $format = "%0.7f";
                foreach my $val ( @{$item} ) {
                    $count++;
                    if ( $count == 1 ) {
                        $str = $val;
                    }
                    else {
                        $val = sprintf( $format, $val );
                        $str .= " & $val";
                    }
                }
                $str .= " \\\\";
            }
            else {
                $str = join( "\t", @{$item} );
            }

            #print STDERR ":$str\n";
            print OUT $str, "\n";
        }
        if ( $format eq "PERL" ) {
            print OUT "];\n";
        }
        elsif ( $format eq "TEX" ) {
            print OUT <<'EOM';
\hline
\end{tabular}
\end{center}
EOM
        }
        close OUT;
    }
    else {

        # Direct to terminal
        foreach my $item ( @{$rtable} ) {
            my $str = join( "\t", @{$item} );
            print $str, "\n";
        }
    }
    print "Wrote '$ofile'\n";
    return;
}

sub query {
    my ($msg) = @_;
    print $msg;
    my $ans = <>;
    chomp $ans;
    return $ans;
}

sub hitcr {
    my ($msg) = @_;
    if ($msg) { $msg .= ". hit <cr> to continue" }
    else      { $msg = "hit <cr> to continue" }
    query($msg);
    ##query( $msg . ". hit <cr>" );
}

sub queryu {
    return uc query(@_);
}

sub get_num {
    my ( $msg, $default ) = @_;
    if ( defined($default) ) {
        $msg =~ s/:$//;
        $msg .= " (<cr>=$default):";
    }
    my $ans = query($msg);
    $ans = $default if ( defined($default) && $ans eq "" );
    my $val = eval($ans);
    if ($@) { warn $@; $val = $ans; }
    return $val;
}

sub ifyes {

    # Updated to have default, which should be "Y" or "N"
    my ( $msg, $default ) = @_;
    my $count = 0;
  ASK:
    my $ans = query($msg);
    if ( defined($default) ) {
        $ans = $default unless ($ans);
    }
    if    ( $ans =~ /^Y/i ) { return 1 }
    elsif ( $ans =~ /^N/i ) { return 0 }
    else {
        $count++;
        if ( $count > 6 ) { die "error count exceeded in ifyes\n" }
        print STDERR "Please answer 'Y' or 'N'\n";
        goto ASK;
    }
}

sub displacement {
    my ( $dV, $rs, $symmetry ) = @_;

    # Given a $dV=a volume change, $rs=starting radius, and symmetry:
    # Compute the displacement of the radius to give the volume change

    # plane symmetry, volume = radius
    if ( $symmetry == 0 ) { return $dV }

    my $pi = 4 * atan2( 1, 1 );
    return   if ( $rs < 0 );
    return 0 if ( $dV == 0 );

    # cylindrical symmetry, volume = pi * r**2
    # so solve dV/pi = (rs+d)**2-rs**2 = (2*rs+d)*d
    if ( $symmetry == 1 ) {
        my $root = $rs**2 - $dV / $pi;
        return -$rs if ( $root <= 0 );
        my $dis = 2 * sqrt($root) - $rs;
        return $dis;
    }

    # spherical symmetry
    my $V0 = 4 / 3 * $pi * $rs**3;
    if ( $dV <= -$V0 ) { return -$rs }

    my $fofx = sub {
        my ($xx) = @_;
        my $rr = $rs + $xx;

        # Note (a-b)(a**2+ab+b**2)=a**3-b**3
        my $rsq      = ( $rs**2 + $rs * $rr + $rr**2 ) / 3;
        my $drsq_dxx = ( $rs + 2 * $rr ) / 3;
        my $ff       = 4 * $pi * $rsq * $xx - $dV;
        my $dfdx     = 4 * $pi * ( $rsq + $xx * $drsq_dxx );
        return ( $ff, $dfdx );
    };

    # newton iterations
    ##my $xx0 = $dV / ( 4 * $pi * $rs**2 );

    # Start X at a value which works for very small rs
    my $xx = ( 3 * $dV / ( 4 * $pi ) )**( 1 / 3 );
    my $tol = 1.e-8;
    for ( my $it = 0 ; $it < 10 ; $it++ ) {
        my ( $ff, $dfdx ) = $fofx->($xx);
        my $dx      = -$ff / $dfdx;
        my $V       = 4 / 3 * $pi * ( $rs + $xx )**3;
        my $dV_test = $V - $V0;

        #print "it=$it, x=$xx, dx=$dx, f=$ff, dfdx=$dfdx dV_test=$dV_test\n";
        $xx += $dx;
        last if ( abs($dx) < $tol * $rs );
    }
    return $xx;
}

sub density_ratios {
    my ( $Y, $gamma ) = @_;
    my $y                   = exp($Y);
    my $rho_amb             = $gamma;
    my $top                 = ( $gamma + 1 ) * ( $y + 1 ) + ( $gamma - 1 );
    my $bot                 = ( $gamma - 1 ) * ( $y + 1 ) + ( $gamma + 1 );
    my $density_ratio_shock = $top / $bot;
    my $density_ratio_ambient =
      $density_ratio_shock * ( 1 / ( 1 + $y ) )**( 1 / $gamma );
    return ( $density_ratio_shock, $density_ratio_ambient );
}

#($p_amb, $sspd_amb)=altitude($p_amb, $sspd_amb);
#print "got ==$p_amb, c=$sspd_amb\n";

sub altitude {

    # return pressure and sound speed at a given altitude
    # leave them unchanged if not selected
    my ( $p_amb, $sspd_amb, $alt_m ) = @_;
    while (1) {
        my $zz_m  = query("altitude, m:");
        my $zz_ft = 3.2808 * $zz_m;
        my $zz_km = $zz_m * 1000;
        my ( $pp, $tt, $rr, $ww ) = atmos($zz_m);
        my $ttf            = 32. + ( $tt - 273.15 ) * 9. / 5.;
        my $psi            = $pp * 14.5e-5;
        my $sspdms         = 331.48 * sqrt( $tt / 273.15 );
        my $sspdz          = 0.001 * 3.2808 * $sspdms;
        my $pressure_ratio = $pp / 1.01325e5;

        foreach ( $zz_km, $zz_ft, $zz_m, $pp, $psi, $pressure_ratio, $tt, $ttf,
            $rr, $ww, $sspdms, $sspdz )
        {
            $_ = sprintf "%0.6g", $_;
        }

        my $menu = <<EOM;
--------Atmosphere----------
altitude, m  .......... $zz_m
altitude, km .......... $zz_km
altitude, ft .......... $zz_ft 
pressure, Pa .......... $pp 
pressure, psi ......... $psi 
pressure ratio, P/P(z=0) $pressure_ratio
temperature, K......... $tt 
temperature, F......... $ttf 
density, kg/m3......... $rr 
molecular wt........... $ww 
sound speed, m/s....... $sspdms 
sound speed, ft/ms..... $sspdz 

Z  enter another altitude
Y  YES..change to these values
Q  Quit; return without changing values
EOM
        print $menu;
        my $ans = queryu('-->');
        next if ( $ans eq 'Z' );

        if ( $ans eq 'Y' ) {
            $p_amb    = $pp;
            $sspd_amb = $sspdms;
            $alt_m    = $zz_m;
            last;
        }
        last if ( $ans eq 'Q' );
    }
    return ( $p_amb, $sspd_amb, $alt_m );
}

{    # sub atmos

    #
    #       atmospheric model, based upon table 2.7 in Regan, Reentry Vehicle...
    #       compute pressure, temperature, density as function of altitude
    # 	    converted from the fortran version with f2pl
    #
    #       input parameter -
    #               zz = altitude in meters
    #       output parameters -
    #               pp = pressure in Pa
    #               tt = temperature in K
    #               rr = density
    #

    # saved arrays
    my ( $rptab, $rrtab, $rttab, $rwtab, $rxlapse, $rztab );

    # saved constants
    my ( $rgas, $grav, $w0, $pref, $b );

    # cached values
    my ( $ppsave, $ttsave, $zzsave, $rrsave, $wwsave );

    INIT {

        #       constants
        $rgas   = 287.;
        $grav   = 9.806;
        $w0     = 28.964;
        $pref   = 1.01325e5;
        $b      = 3.139e-7;
        $ppsave = 0.;
        $ttsave = 0.;
        $zzsave = -1.;
        $rrsave = 0.;
        $wwsave = 0.;

        #       altitude, km
        $rztab = [
            0.0,   11.019, 20.063, 32.162, 47.435, 52.43, 61.59, 80.0,
            90.0,  100.0,  110.0,  120.0,  150.0,  160.0, 170.0, 190.0,
            230.0, 300.0,  400.0,  500.0,  600.0,  700.0
        ];

        #       temp
        $rttab = [
            288.1,   216.65,  216.65,  228.65,  270.65,  270.65,
            252.65,  180.65,  180.65,  210.65,  260.65,  360.65,
            960.65,  1110.6,  1210.65, 1350.65, 1550.65, 1830.65,
            2160.65, 2420.65, 2590.65, 2700.65
        ];

        #       p/p0
        $rptab = [
            1.000,     2.284e-1, 5.462e-2, 8.567e-3, 1.095e-3,  5.823e-4,
            1.797e-4,  1.024e-5, 1.622e-6, 2.98e-7,  7.22e-8,   2.488e-8,
            5.0e-9,    3.64e-9,  2.756e-9, 1.66e-9,  6.869e-10, 1.86e-10,
            3.977e-11, 1.08e-11, 3.4e-12,  1.176e-12
        ];

        #       molecular weight
        $rwtab = [
            28.964, 28.964, 28.964, 28.964, 28.964, 28.964, 28.964, 28.964,
            28.964, 28.88,  28.56,  28.07,  26.92,  26.66,  26.5,   25.85,
            24.69,  22.66,  19.94,  17.94,  16.84,  16.17
        ];

        for ( my $n = 0 ; $n < @{$rztab} ; $n++ ) {
            $rztab->[$n] *= 1000.;
            $rptab->[$n] *= $pref;
            $rrtab->[$n] = $rptab->[$n] / ( $rgas * $rttab->[$n] );
        }
        for ( my $n = 1 ; $n < @{$rztab} ; $n++ ) {
            $rxlapse->[ $n - 1 ] =
              ( $rttab->[$n] - $rttab->[ $n - 1 ] ) /
              ( $rztab->[$n] - $rztab->[ $n - 1 ] );
        }

    }

    sub atmos {

        my ($zz) = @_;

        my ( $pp, $tt, $rr, $ww );

        #       fixup altitude
        if ( $zz < 0.0 )          { $zz = 0.0 }
        if ( $zz > $rztab->[-1] ) { $zz = $rztab->[1]; }

        #       check cache
        if ( $zz == $zzsave ) {
            $pp = $ppsave;
            $tt = $ttsave;
            $rr = $rrsave;
            $ww = $wwsave;
            goto L900;
        }

        #       lookup the altitude
        for ( my $n = 1 ; $n < @{$rztab} ; $n += 1 ) {
            if ( $zz <= $rztab->[$n] ) {
                my $i  = $n - 1;
                my $dz = $zz - $rztab->[$i];

                #         in isothermal layer
                if ( abs( $rxlapse->[$i] ) < 1.e-5 ) {
                    my $tt = $rttab->[$i];
                    my $q8 =
                      $grav * $dz *
                      ( 1. - ( $b / 2. ) * ( $zz + $rztab->[$i] ) ) /
                      ( $rgas * $rttab->[ $i - 1 ] );
                    my $expq8 = exp( -$q8 );
                    $pp = $expq8 * $rptab->[$i];
                    $rr = $expq8 * $rrtab->[$i];
                }

                #         in non-isothermal layer
                else {

                    my $q1 = 1. +
                      $b * ( $rttab->[$i] / $rxlapse->[$i] - $rztab->[$i] );
                    my $q2 = ( $q1 * $grav ) / ( $rgas * $rxlapse->[$i] );
                    $tt = $rttab->[$i] + $rxlapse->[$i] * $dz;
                    my $q3 = $tt / $rttab->[$i];
                    my $q4 = $q3**( -$q2 );
                    my $q5 =
                      exp( $b * $grav * $dz / ( $rgas * $rxlapse->[$i] ) );
                    my $q6 = $q4 * $q5;
                    $pp = $rptab->[$i] * $q6;
                    my $q7 = $q2 + 1.0;
                    $rr = $rrtab->[$i] * ( $q3**( -$q7 ) ) * $q5;
                }
                $ww =
                  $rwtab->[$i] +
                  ( $rwtab->[ $i + 1 ] - $rwtab->[$i] ) /
                  ( $rztab->[ $i + 1 ] - $rztab->[$i] ) *
                  $dz;

                # ? unused var
                my $t9 = ( $ww * $tt ) / $w0;
                goto L900;
            }
        }
        print STDERR "**error** in atmos\n";

        #       save results
      L900:
        $ppsave = $pp;
        $zzsave = $zz;
        $ttsave = $tt;
        $rrsave = $rr;
        $wwsave = $ww;
        return ( $pp, $tt, $rr, $ww );
    }
}
__END__
Example output, english units:

 Slant Range, ft.....................    328.0 
 Weight, lb..........................    1.0 
  
 Incident Overpressure, psi..........    8.333975241E-2 
 Ground (measured) overpressure, psi.    8.333975241E-2 
 Normal reflected overpressure, psi..    0.1666795048 
  
 Incident Impulse, psi-ms............    0.2713459657 
 Ground (measured) impulse, psi-ms...    0.2713459657 
 Normal reflected impulse, psi-ms....    0.5426919315 
  
 Time of arrival, ms.................    286.3715576 
 Positive Phase duration, ms.........    6.946690553 
 Shape factor, alpha.................    0.1971577722 
 Shock velocity, kfps................    1.118106658 
 hit <cr> to continue

