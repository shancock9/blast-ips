#!/usr/bin/perl 
use warnings;
use strict;

# This is a driver to illustrate usage of Blast::IPS.
use Blast::IPS;
use Blast::IPS::Medium;

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
print <<EOM;
Please select a starting mode (you can switch anytime):
D  - Dimensionless
SI - with units (default SI units)
EOM
my $selection = queryu("<cr>='SI':");
$selection = 'SI' unless $selection eq 'D';
while (1) {
    if ( $selection eq 'D' ) {
        ( $blast_table, $selection ) = eval_dimensionless( $blast_table );
    }
    elsif ( $selection eq 'SI' ) {
        #( $blast_table, $selection) = eval_si( $blast_table );
        ( $blast_table, $selection) = point_evaluations_SI( $blast_table );
    }
    else { last; }
}

sub ask_ground_plane {
    my $ground_plane;
    if ( ifyes("Is this explosion on a hard surface (ground plane)?[Y/N]") ) {
        $ground_plane = "YES";
    }
    else {
        $ground_plane = "NO";
    }
    return $ground_plane;
}

sub eval_dimensionless {

    my ( $blast_table, $ground_plane ) = @_;

    my $symmetry = $blast_table->get_symmetry();
    my $gamma    = $blast_table->get_gamma();
=pod
    my $medium   = {
        _gamma        => $gamma,
        _sspd_amb     => 1,
        _p_amb        => 1,
        _symmetry     => $symmetry,
        _rho_amb      => $gamma,
        _E0           => 1,
        _ground_plane => $ground_plane
    };
=cut
    my $medium   = Blast::IPS::Medium->new(
        gamma        => $gamma,
        sspd_amb     => 1,
        p_amb        => 1,
        symmetry     => $symmetry,
        rho_amb      => $gamma,
        E0           => 1,
        ground_plane => $ground_plane
    );

    # main loop
    my $return_selection;
    while (1) {
        my $table_name = $blast_table->get_table_name();
        $symmetry = $blast_table->get_symmetry();
        $gamma    = $blast_table->get_gamma();
        print <<EOM;
==============================================================
Main Menu - Dimensionless.  Menu Commands are Case Insensitive
==============================================================
First check/set the geometry and atmosphere parameters...
  S     - Symmetry: $symmetry_name{$symmetry}
        - gamma: $gamma

Then evaluate the model with these settings...
  C     - Calculate with these settings
  T     - do Table operations ...
  I     - show Global Information about this blast

Use one of these to go back...
  SI    - switch to Units mode (default SI units)
  Q     - Quit the program
EOM
        my $ans = queryu(":");
        if ( $ans eq 'S' ) {
	    print <<EOM;

=========================================================================
NOTE: In dimensionless mode you are working with one of three symmetries:

P : an infinte planar source on a rigid wall
C : an infinitely long cylinder source in air (no ground plane)
S : a sphere in free air (no ground plane)

- The scaling distance is (E/P0)^(1/n), where E is energy and n=1,2, or 3.  
- You must handle any ground plane yourself:  
  - For spherical and cylindrical geometries, if you have a ground plane, then 
   you need to use TWO TIMES your energy in this equation.
  - For planar geometry, if you do not have a ground plane (rigid wall),
   you need to use ONE HALF your energy in this equation.
- To avoid dealing with this issue, use the SI mode which does this for you

EOM
            ( $blast_table, $medium ) =
              select_blast_table( $blast_table, $medium );
        }
        elsif ( $ans eq 'I' ) {
            show_summary_information( $blast_table, $medium );
        }
        elsif ( $ans eq 'SI' ) {
            $return_selection = 'SI';
            last;
        }
        elsif ( $ans eq 'C' ) {
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

sub select_geometry {
    my ( $blast_table, $medium ) = @_;
    my $symmetry_now = $blast_table->get_symmetry();
    my $gamma = $blast_table->get_gamma();
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

You may either use an altitude model or adjust the parameters individually.

Z   - Set earth altitude; current value............: $alt_m m 

G   - change gamma; current value..................: $gamma 
P   - change ambient pressure; currrent value......: $p_amb Pa
C   - change ambient sound speed; current value....: $sspd_amb m/s
Y   - Yes, return with these values    Q - quit; use previous values
EOM
            my $ans = queryu('-->');
            if ( $ans ne 'Y' && $ans ne 'Q' ) { $touched = 1 }
            if ( $ans eq 'Z' ) {
                $gamma = 1.4;
                ( $p_amb, $sspd_amb, $alt_m ) =
                  altitude( $p_amb, $sspd_amb, $alt_m );
                $medium->{_gamma}    = $gamma;
                $medium->{_p_amb}    = $p_amb;
                $medium->{_sspd_amb} = $sspd_amb;
                my $rho_amb = $gamma * $p_amb / $sspd_amb**2;
                $medium->{_rho_amb} = $rho_amb;
                last;
            }
            elsif ( $ans eq 'G' ) {
                $gamma = get_num( "Enter gamma", 1.4 );
                $alt_m = "?";
            }
            elsif ( $ans eq 'P' ) {
                #$p_amb = get_num( "Enter ambient pressure, Pa", $p_amb );
                $p_amb = request_dimensional_value( "Enter ambient pressure, Pa:",
                    'P', $p_amb, 0 );
                $alt_m = "?";
            }
            elsif ( $ans eq 'C' ) {
                #$sspd_amb = get_num( "Enter ambient sound speed, m/s", $sspd_amb );
                $sspd_amb =
                  request_dimensional_value( "Enter ambient sound speed, m/s",
                   'L/T', $sspd_amb, 0 );
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
                last if ( ifyes("Return without keeping your changes?[Y/N]") );
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

sub get_timeline_dimensionless {
    my ( $blast_table, $medium, ) = @_;


    # FORMAT:
    # First line is a header, then come events:
    # [ti, ri, rsi, "description"]

    my $rinfo                   = $blast_table->get_info();
    my $rtimeline;
    my $symmetry                = $rinfo->{symmetry};
    my $gamma                   = $rinfo->{gamma};

    # First line is the header
    push @{$rtimeline}, ["T", "R", "Rs", "Event"];

# TBD
    my $t_pcenter_zero_1        = $rinfo->{t_pcenter_zero_1};
    my $t_pcenter_min           = $rinfo->{t_pcenter_min};

    my $pcenter_min             = $rinfo->{pcenter_min};
    my $t_dpdt_center_max       = $rinfo->{t_dpdt_center_max};
    my $dpdt_center_max         = $rinfo->{dpdt_center_max};
    my $pcenter_dpdt_center_max = $rinfo->{pcenter_dpdt_center_max};
    my $r_tail_shock            = $rinfo->{r_tail_shock};
    my $z_tail_shock            = $rinfo->{z_tail_shock};
    my $r_tail_formed           = $rinfo->{r_tail_formed};
    my $z_tail_formed           = $rinfo->{z_tail_formed};

#    foreach ( $alpha, $Sint_pos, $Sint_neg, $Ixr_pos, $Ixr_neg ) {
#        $_ = sprintf( "%0.6g", $_ );
#    }
=pod
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
=cut

    return $rtimeline;
}

sub show_summary_information {
    my ( $blast_table, $medium ) = @_;
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
            'I-X'   => [ ++$i, 'ln(i/x)'],
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

sub request_value {
    my ( $msg, $default ) = @_;
    my $test = get_num( $msg, $default );
    if ( !defined($test) ) { return $default }
    return $test;
}

sub request_positive_value {
    my ( $msg, $default) = @_;
    my $test = get_num( $msg, $default );
    if ( !defined($test) ) { return $default }
    if ( $test <= 0 ) { return $default }
    return $test;
}

{

    my $range;
    my $medium;

    BEGIN {
        $range = undef;
        my $gamma          = 1.4;
        my $symmetry       = 2;
        my $p_sea_level    = 1.01325e5;
        my $sspd_sea_level = 340.43;
        my $E0             = 1;
        $medium     = {
            _gamma        => $gamma,
            _p_amb        => $p_sea_level,
            _sspd_amb     => $sspd_sea_level,
            _symmetry     => $symmetry,
            _rho_amb      => $gamma * $p_sea_level / $sspd_sea_level**2,
            _E0           => undef,
            _ground_plane => 'YES',
        };
    }

    sub point_evaluations_SI {
        my ( $blast_table ) = @_;

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

	# Flags to avoid overwriting user settings
        my $user_set_E0;
        my $user_set_range;

        my $ground_factor = ground_factor( $symmetry, $ground_plane );
        my $set_dscale = sub {

            # Set the scaling distance whenever E0 changes
            return unless defined($E0);
            $dscale =
              ( $ground_factor * $E0 / $p_amb )**( 1 / ( $symmetry + 1 ) );
        };
        my $set_E_from_dscale = sub {
	    my $E0_save=$E0;
            $E0 = $p_amb * $dscale**( $symmetry + 1 ) / $ground_factor;
            if ($user_set_E0) {
                my $E0_str     = format_E($E0, $symmetry); 
                my $E0_save_str     = format_E($E0_save, $symmetry); 
                if ( !ifyes(<<EOM) )
This gives E0=$E0_str.  You previously entered E=$E0_save_str.
Change E0 to be $E0_str? [Y/N]
EOM
                {
                    $E0 = $E0_save;
                    $set_dscale->();
                }
		else {
                    $user_set_E0 = 0;
		}
            }
        };
        my $set_range_from_dscale = sub {
            my ($x) = @_;
            my $range_save = $range;
            $range = $x * $dscale;
            if ($user_set_range) {
                my $range_str  = format_value($range, 'L'); 
                my $range_save_str  = format_value($range_save, 'L'); 

		# Range may be well known, so ask before changing.
                if ( !ifyes(<<EOM) )
This gives range=$range_str. You previously entered range=$range_save_str
Change range to be $range_str? [Y/N]
EOM
                {
                    $range = $range_save;
                }
                else {
                    $user_set_range = 0;
                }
            }
        };
        my $ask_for_E0 = sub {
	    if ($symmetry == 0) {
            	$E0 = request_dimensional_value( 
                "Enter Energy E0, (J/m^2):", 'E/L^2',$E0,0);
	    }
	    elsif ($symmetry == 1) {
            	$E0 = request_dimensional_value( 
                "Enter Energy E0, (J/m):", 'E/L',$E0,0);
	    }
	    else {
            	$E0 = request_dimensional_value( 
                "Enter Energy E0, (J):", 'E',$E0,0);
            }
            $user_set_E0=1;;
            $set_dscale->();
        };
        my $ask_for_range = sub {
            $range =
              request_dimensional_value( "Enter range, m:", 'L', $range, 0 );
            $user_set_range=1;
        };
        my $ask_for_overpressure_ratio = sub {

	    my ($ans)=@_;
            my $y;
            if ($ans =~ /Y/) { 
                $y =
                  request_positive_value("Enter incident overpressure_ratio");
            }
            else {
                my $ovp = request_dimensional_value( "Enter Overpressure, Pa:",
                    'P', undef, 0 );
                $y = $ovp / $p_amb;
            }
            return $y;
        };
        my $ask_for_impulse = sub {

            my $imp = request_dimensional_value( 
                "Enter incident positive phase impulse, Pa-s:", 'P*T',undef,0);
            return $imp;
        };

        my $Nsym = $symmetry + 1;    # = (1,2 or 3)
        $set_dscale->();
        my $pstr = $symmetry == 0 ? "" : "^1/" . ( $symmetry + 1 );

        my $return_selection = 'D';
        while (1) {
            my $range_str  = format_value($range, 'L'); 
            my $dscale_str = format_value($dscale, 'L'); 
            my $p_amb_str   = format_value($p_amb, 'P'); 
            my $sspd_amb_str   = format_value($sspd_amb, 'L/T'); 
            my $E0_str     = format_E($E0, $symmetry); 
            $ground_factor = ground_factor( $symmetry, $ground_plane );

            print <<EOM;
=========================================================
Main Menu - SI units.  Menu Commands are Case Insensitive
=========================================================
1. Basic Settings:
   S  symmetry.................... $symmetry_name{$symmetry}
   G  explosion on a ground plane? $ground_plane (g = $ground_factor)
   A  atmospheric Pressure, Pa.... $p_amb_str
      ambient sound speed, m/s.... $sspd_amb_str
      gamma....................... $gamma
   E  Energy...................... $E0_str   
   R  Range....................... $range_str
      scale dist. (g*E0/P0)$pstr... $dscale_str m 

2. Commands to estimate E and/or R from two quantities:
   RI         (Range, Impulse)->Energy
   RP             (Range, OVP)->Energy 
   RY       (Range, OVP ratio)->Energy 
   RT             (Range, TOA)->Energy
   PI           (OVP, Impulse)->(Energy, Range)
   YI     (OVP ratio, Impulse)->(Energy, Range)
   EP            (Energy, OVP)->Range         
   EY      (Energy, OVP ratio)->Range         
   ET            (Energy, TOA)->Range         

3. Other Operations: 
   C Calculate blast parameters for the current value of Energy
   I show global Information 
   F Files (read/write) ...
   Q exit program    D switch to dimensionless mode
EOM
#  Z Zoom  - view latest results, 1 screen per channel
#  L List  - view latest results, 1 line per channel
#  F files (read/write)...
#  Q=Quit   CL=Clear  LG=List Gages LD=List Data
            my $ans = queryu(":");
	    $ans =~ s/\*$//; # in case user took the asterisk literally
            if ( $ans eq 'S' ) {
                ( $blast_table, $medium ) =
                  select_geometry( $blast_table, $medium );
		$symmetry=$medium->{_symmetry};
                $ground_plane = ask_ground_plane();
                $medium->{_ground_plane}=$ground_plane;
            }
            elsif ( $ans eq 'F' ) {
                query("Sorry, incomplete");
            }
            elsif ( $ans eq 'E' ) {
                $ask_for_E0->();
            }
            elsif ( $ans eq 'G' ) {
                $ground_plane = ask_ground_plane();
                $medium->{_ground_plane}=$ground_plane;
            }
            elsif ( $ans eq 'A' ) {
                ( $blast_table, $medium ) =
                  select_atmosphere_SI( $blast_table, $medium );
                $gamma        = $medium->{_gamma};
                $symmetry     = $medium->{_symmetry};
                $p_amb        = $medium->{_p_amb};
                $sspd_amb     = $medium->{_sspd_amb};
                $ground_plane = $medium->{_ground_plane};
                $E0           = $medium->{_E0};
            }
            elsif ( $ans eq 'R' ) {
                $ask_for_range->();
            }
            elsif ( $ans eq 'RE' || $ans eq 'ER' || $ans eq 'C' ) {
                if ( !defined($range) ) { $ask_for_range->(); }
                if ( !defined($E0) )    { $ask_for_E0->(); }
		$medium->{_E0}=$E0;
		display_result_SI($range, $medium, $blast_table); 
            }
            elsif ( $ans eq 'RI' || $ans eq 'IR' ) {
                if ( !defined($range) ) { $ask_for_range->(); }
                my $impulse = $ask_for_impulse->();
                next if ( !defined($impulse) || $impulse <= 0 );
                next if ( !defined($range)   || $range <= 0 );
                my $term =
                  log( ( $impulse * $sspd_amb ) / ( $p_amb * $range ) );
                my $ret = $blast_table->wavefront( 'I-X' => $term );
                my $X   = $ret->{X};
                my $x   = exp($X);
                $dscale = $range / $x;
                $set_E_from_dscale->();
                next;
            }
            elsif ( $ans eq 'RT' || $ans eq 'TR' ) {
                if ( !defined($range) ) { $ask_for_range->(); }
                my $t = request_dimensional_value( 
                  "Enter toa, s:", 'T',undef,0);
                next if ( $t <= 0 );
                my $z = ( $range - $t * $sspd_amb );
                if ( $z <= 0 ) {
                    my $toa_max = $range / $sspd_amb;
                    query("toa must not exceed $toa_max at range $range");
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
                $set_E_from_dscale->();
            }
            elsif ( $ans eq 'PI' || $ans eq 'YI' ) {
                my $y = $ask_for_overpressure_ratio->($ans);
                next if ( !defined($y) || $y <= 0 );
                my $Y       = log($y);
                my $ret     = $blast_table->wavefront( 'Y' => $Y );
                my $X       = $ret->{X};
                my $x       = exp($X);
                my $Ixr_pos = $ret->{Ixr_pos};
                my $I_pos   = $Ixr_pos / $x**( $symmetry / 2 );
                my $impulse = $ask_for_impulse->();
                next if ( !defined($impulse) || $impulse <= 0 );
                $dscale = $impulse * $sspd_amb / ( $I_pos * $p_amb );
                $set_E_from_dscale->();
		$set_range_from_dscale->($x);
            }
            elsif ( $ans eq 'RP' || $ans eq 'RY' ) {
                if ( !defined($range) ) { $ask_for_range->(); }
		my $y = $ask_for_overpressure_ratio->($ans);
                next if ( !defined($y) || $y <= 0 );
                my $Y    = log($y);
                my $ret  = $blast_table->wavefront( 'Y' => $Y );
                my $X    = $ret->{X};
                my $dYdX = $ret->{dYdX};
                my $x    = exp($X);
                $dscale = $range / $x;
                $set_E_from_dscale->();
            }
            elsif ( $ans eq 'EP' || $ans eq 'EY' ) {
                if ( !defined($E0) ) { $ask_for_E0->(); }
		my $y = $ask_for_overpressure_ratio->($ans);
                next if ( !defined($y) || $y <= 0 );
                my $Y    = log($y);
                my $ret  = $blast_table->wavefront( 'Y' => $Y );
                my $X    = $ret->{X};
                my $dYdX = $ret->{dYdX};
                my $x    = exp($X);
		$set_range_from_dscale->($x);
            }
            elsif ( $ans eq 'ET' || $ans eq 'TE' ) {
                if ( !defined($E0) ) { $ask_for_E0->(); }
                my $t = request_dimensional_value( 
                  "Enter toa, s:", 'T',undef,0);
                next if ( $t <= 0 );
                my $T    = log( $t * $sspd_amb / $dscale );
                my $ret  = $blast_table->wavefront( 'T' => $T );
                my $X    = $ret->{X};
                my $Y    = $ret->{Y};
                my $dYdX = $ret->{dYdX};
                my $y    = exp($Y);
                my $x    = exp($X);
		$set_range_from_dscale->($x);
            }
            elsif ( $ans eq 'D' ) {
                $return_selection = 'D';
                last;
            }
            elsif ( $ans eq 'I' ) {
                query("Coding incomplete");
            }
            elsif ( $ans eq 'Q' ) {
                if ( ifyes("Quit blast_ips? [Y/N]") ) {
                    $return_selection = 'Q';
                    last;
                }
            }
            else {
		# unknown response, keep going
            }
        }
	$medium->{_E0} = $E0;
        return ( $blast_table, $return_selection); 
    }
}

sub display_result_SI {
    my ( $range, $medium, $blast_table ) = @_;

    my $gamma         = $medium->{_gamma};
    my $symmetry      = $medium->{_symmetry};
    my $p_amb         = $medium->{_p_amb};
    my $sspd_amb      = $medium->{_sspd_amb};
    my $ground_plane  = $medium->{_ground_plane};
    my $E0            = $medium->{_E0};

    my $ground_factor = ground_factor( $symmetry, $ground_plane );
    my $dscale =
      ( $ground_factor * $E0 / $p_amb )**( 1 / ( $symmetry + 1 ) );
    my $x           = $range / $dscale;
    my $X           = log($x);
    my $ret         = $blast_table->wavefront( 'X' => $X );
    my $Y           = $ret->{Y};
    my $y           = exp($Y);
    my $dYdX        = $ret->{dYdX};
    my $T           = $ret->{T};
    my $toa         = exp($T) * $dscale / $sspd_amb;
    my $Z           = $ret->{Z};
    my $z           = exp($Z) * $dscale;
    my $z_pose_rs   = $ret->{z_pose_rs} * $dscale;
    my $tpos        = ( $z - $z_pose_rs ) / $sspd_amb;

    ###################################################################
    #my $imp         = $ret->{pint_pos} * $p_amb * $dscale / $sspd_amb;
    my $ImX         = $ret->{'I-X'};
    my $I           = $ImX + $X;
    my $imp         = exp($I) * $p_amb * $dscale / $sspd_amb;
    ###################################################################

    # FIXME: get these from the hash
    my $term        = $y * ( $gamma + 1 ) / ( 2 * $gamma );
    my $m           = sqrt( 1 + $term );
    my $up          = $y / ( $gamma * $m );
    my $shock_speed = $m * $sspd_amb;

    my $range_str = format_value( $range, 'L' );
    my $E0_str     = format_E($E0, $symmetry); 
    my $toa_str   = format_value( $toa, 'T' );
    my $tpos_str  = format_value( $tpos, 'T' );
    my $ovp       = $y * $p_amb;
    my $ovp_str   = format_value( $ovp, 'P' );
    my $y_str     = format_value($y);
    my $imp_str   = format_value( $imp, 'P*T' );
    my $shock_speed_str   = format_value( $shock_speed, 'L/T' );
    print <<EOM;
Energy...................... $E0_str   
Range....................... $range_str
Overpressure ratio.......... $y_str 
Overpressure ............... $ovp_str 
Incident Impulse............ $imp_str
Time of arrival............. $toa_str
Positive Phase duration..... $tpos_str 
Shock Speed................. $shock_speed_str 
EOM
    hitcr();
    return;
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
          query("Enter $another value for '$vname', or <cr> to quit:");
        $another = 'another';
        last if !$val  || $val !~ /\d/;    #^\s*[\-\+\d\.]/;

        my ( $iQ, $Q );
        if ( $vname =~ /^([XYZT]|E1|Er|dYdX|dZdX|dTdX|dE1dX|dErdX|I-X)/ ) {
            $Q  = $val;
            $iQ = $vname;
        }
        elsif ( $vname =~ /^[xyzt]$/ ) {
            next if ( $val <= 0 );
            $Q  = log($val);
            $iQ = uc($vname);
        }
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
        my $dpdr_t = $ret->{dpdr_t};
        my $dudr_t = $ret->{dudr_t};
        my $dpdt_r = $ret->{dpdt_r};
        my $dudt_r = $ret->{dudt_r};

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
	    $dpdr_t, $dudr_t, $dpdt_r, $dudt_r,
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
dpdt|r = $dpdt_r  dudt|r = $dudt_r  shock front slopes
dpdr|t = $dpdr_t  dudr|t = $dudr_t  shock front slopes

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
            ##print OUT audit_string_brief('#');
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
        my $root = $rs**2 + $dV / $pi;
        return -$rs if ( $root <= 0 );
        my $dis = sqrt($root) - $rs;
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

sub altitude {

    # return pressure and sound speed at a given altitude
    # leave them unchanged if not selected
    my ( $p_amb, $sspd_amb, $alt_m ) = @_;
    my $zz_m = $alt_m;
    while (1) {
        $zz_m =
              request_dimensional_value( "Altitude, m:", 'L', $zz_m );
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
Q  return 
EOM
        print $menu;
        my $ans = queryu('-->');
        next if ( $ans eq 'Z' );

        if ( $ans eq 'Q' ) {
            $p_amb    = $pp;
            $sspd_amb = $sspdms;
            $alt_m    = $zz_m;
            last;
        }
    }
    return ( $p_amb, $sspd_amb, $alt_m );
}

{    # sub atmos

    #
    #       atmospheric model, based upon table 2.7 in Regan, Reentry Vehicle...
    #       compute pressure, temperature, density as function of altitude
    # 	    converted from fortran with f2pl
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

{   # i/o with units

    my ( $kt_to_joules, $pa_to_psi, $m_to_ft, $rto_SI, $rSI_name );

    BEGIN {
        $kt_to_joules = 4.184e12;
        $pa_to_psi    = 0.00014503773800722;
        $m_to_ft      = 3.2808;
        $rto_SI       = {
            'P' => {
                'psi' => 1 / $pa_to_psi,
                'Pa'  => 1,
                'kPa' => 1.e3,
                'k'   => 1.e3,
                'MPa' => 1.e6,
                'M'   => 1.e6,
            },
            'L' => {
                'ft' => 1 / $m_to_ft,
                'kft' => 1000 / $m_to_ft,
                'm'  => 1,
                'km' => 1000,
            },
            'T' => {
                's' => 1, 
                'ms' => 0.001,
            },
            'E' => {
                'J'   => 1,
                'kJ'  => 1000,
                'MJ'  => 1.e6,
                't'   => 0.001 * $kt_to_joules,
                'kt'  => $kt_to_joules,
                'Mt'  => 1000 * $kt_to_joules,
                'lb'  => 0.001 * $kt_to_joules / 2000,
                'klb' => $kt_to_joules / 2000,
            },
            'E/L' => {
                'J/m'   => 1,
                'kJ/m'  => 1000,
                'MJ/m'  => 1.e6,
                't/m'   => 0.001 * $kt_to_joules,
                'kt/m'  => $kt_to_joules,
                'Mt/m'  => 1000 * $kt_to_joules,
                'lb/m'  => 0.001 * $kt_to_joules / 2000,
                'klb/m' => $kt_to_joules / 2000,
            },
            'E/L^2' => {
                'J/m^2'   => 1,
                'kJ/m^2'  => 1000,
                'MJ/m^2'  => 1.e6,
                't/m^2'   => 0.001 * $kt_to_joules,
                'kt/m^2'  => $kt_to_joules,
                'Mt/m^2'  => 1000 * $kt_to_joules,
                'lb/m^2'  => 0.001 * $kt_to_joules / 2000,
                'klb/m^2' => $kt_to_joules / 2000,
            },
            'L/T' => {
                'm/s'  => 1,
                'km/s' => 1000,
                'mps'  => 1,
                'fps'  => 1 / $m_to_ft,
                'kfps' => 1000 / $m_to_ft,
              },
            'P*T' => {
                'Pa-s'   => 1,
                'psi-s'  => 1 / $pa_to_psi,
                'psi-ms' => 0.001 / $pa_to_psi,
                'fps'    => 1 / $m_to_ft,
            },
        };
        $rSI_name       = {
            'P' => 'Pa',
            'T' => 's',
            'L' => 'm',
            'E' => 'J',
            'L/T' => 'm/s',
            'P*T' => 'Pa-s',
        };
    }

    sub request_dimensional_value {
        my ( $msg, $unit_type, $default, $min ) = @_;

	# Get a number with possible unit
        my $rconversion_factors = $rto_SI->{$unit_type};
        if ( !defined($rconversion_factors) ) {
            die "unknown unit type '$unit_type' in get_num_and_unit\n";
        }

        my @keys = sort ( keys %{$rconversion_factors} );
        if ( !defined($msg) ) { $msg=""; }
	if (defined($default)) {$msg .= "( (<cr>=$default):"}
	my $SI_name = $rSI_name->{$unit_type};
	$msg = <<EOM;
---
Default unit is $SI_name; for other units, append any of: (@keys)
$msg
EOM
        my $unit = "";
        my $val;
        while (1) {
            my $ans = query("$msg");

            $ans =~ s/\s+$//;
            $ans =~ s/^\s+//;

            # Examples of valid responses:
            # 4 psi, 4psi, 4.5 Pa-s, 10psi-ms,  3/3.2808,
            my $factor = 1;
            if ( $ans =~ /(.*[^A-Za-z\-])\s*([A-Za-z\-]+)$/ ) {
                $val  = $1;
                $unit = $2;
                $unit =~ s/\s*$//;
                $factor = $rconversion_factors->{$unit};
                if ( !defined($factor) ) {
                    local $" = ')(';
                    my @keys = sort ( keys %{$rconversion_factors} );
                    print <<EOM;
Unexpected unit: '$unit'
Expecting on of: (@keys)
EOM
                    next;
                }
            }
            else { $val = $ans }
            if ( !defined($val) ) { $val = $default }

            # optional, allow math in the numerical part
            else { $val = eval($val); }
            $val = $factor * $val;
	    if (defined($min) && $val<=$min) {
		print "***Value must be >$min; try again***\n";
		next;
	    }
            last;
        }
        return $val;
    }

    sub format_value {
        my ( $val, $unit_type ) = @_;

        # format $val, which is in SI units

        ##my $pa_to_psi = 0.00014503773800722;
        ##my $m_to_ft = 3.2808;

        # Format a numerical value in SI units as a string with
        # additional useful units
        my $str = '?';
        if ( defined($val) ) {
            $str = sprintf( "%0.6g", $val );
            if ( defined($unit_type) ) {
                if ( $unit_type eq 'L' ) {
                    my $val_ft = $val * $m_to_ft;
                    $val    = sprintf( "%0.6g", $val );
                    $val_ft = sprintf( "%0.6g", $val_ft );
                    $str    = "$val m = $val_ft ft";
                }
                elsif ( $unit_type eq 'T' ) {
                    $val = sprintf( "%0.6g", $val );
                    $str = "$val s";
                }
                elsif ( $unit_type eq 'L/T' ) {
                    $val = sprintf( "%0.6g", $val );
                    my $val_fps = sprintf( "%0.6g", $val * $m_to_ft );
                    $str = "$val m/s = $val_fps fps";
                }
                elsif ( $unit_type eq 'E' ) {
                    my $val_kt = $val / 4.184e12;
                    my $val_lb = $val_kt * 2.e6;
                    my $val_mj = $val / 1.e6;
                    $val_kt = sprintf( "%0.6g", $val_kt );
                    $val_mj = sprintf( "%0.6g", $val_mj );
                    $val_lb = sprintf( "%0.3g", $val_lb );
                    $str    = "$val_mj MJ = $val_kt kt =~ $val_lb lb TNT";
                }
                elsif ( $unit_type eq 'E/L' ) {
                    $val = sprintf( "%0.6g", $val );
                    $str = "$val J/m";
                }
                elsif ( $unit_type eq 'E/L^2' ) {
                    $val = sprintf( "%0.6g", $val );
                    $str = "$val J/m^2";
                }
                elsif ( $unit_type eq 'P' ) {
                    my $val_psi = $val * $pa_to_psi;
                    $val     = sprintf( "%0.6g", $val );
                    $val_psi = sprintf( "%0.6g", $val_psi );
                    $str     = "$val Pa = $val_psi psi";
                }
                elsif ( $unit_type eq 'P*T' ) {
                    $val = sprintf( "%0.6g", $val );
                    my $val_psi_ms =
                      sprintf( "%0.6g", $val * $pa_to_psi * 1000 );
                    $str = "$val Pa-s  =  $val_psi_ms psi-ms";
                }
            }
        }
        return $str;
    }

    sub format_E {
        my ( $val, $symmetry ) = @_;
        my $str;
        if ( $symmetry == 2 ) {
            $str = format_value( $val, 'E' );
        }
        elsif ( $symmetry == 1 ) {
            $str = format_value( $val, 'E/L' );
        }
        else {
            $str = format_value( $val, 'E/L^2' );
        }
        return $str;
    }

}
__END__

