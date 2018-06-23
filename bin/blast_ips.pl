#!/usr/bin/perl 
use warnings;
use strict;

# This is a driver to illustrate usage of Blast::IPS.
use Blast::IPS;

my $audit_string = "";

my %symmetry_name = (
    0 => 'Plane',
    1 => 'Cylindrical',
    2 => 'Spherical',
);

# Setup a default blast table
my $gamma       = 1.4;
my $symmetry    = 'S';
my $blast_table = Blast::IPS->new( symmetry => $symmetry, gamma => $gamma );

my $sspd_amb = 1;
my $p_amb    = 1;                                  # ambient pressure
my $rho_amb  = $gamma * $p_amb / $sspd_amb**2;
my $medium   = {
    _gamma    => $gamma,
    _sspd_amb => $sspd_amb,
    _rho_amb  => $rho_amb,
    _p_amb    => $p_amb,
    _symmetry => $symmetry,
};

# main loop
while (1) {
    my $table_name = $blast_table->get_table_name();
    $symmetry = $blast_table->get_symmetry();
    $gamma    = $blast_table->get_gamma();
    print <<EOM;
==============================================
MAIN MENU: Point Source Explosion in Ideal Gas 
==============================================
Enter one of the following:
  N     - New Table - Change Symmetry and/or Gamma
          Symmetry=$symmetry_name{$symmetry},  Gamma=$gamma
  P1    - Point evaluations with 1 dimensionless variable...
  P2    - Point evaluations with 2 dimensional variables...
  T     - Table operations ..
  Q     - Quit
EOM
    my $ans = queryu(":");
    if ( $ans eq 'N' ) {
        ($blast_table, $medium) = select_blast_table($blast_table, $medium);
    }
    elsif ( $ans eq 'P1' ) {
        my $vname = 'X';
        $vname = select_variable($vname);
        point_evaluations_dimensionless($blast_table, $medium, $vname);
    }
    elsif ( $ans eq 'P2' ) {
        point_evaluations_with_units($blast_table, $medium);
    }
    elsif ( $ans eq 'T' ) {
        table_operations($blast_table, $medium);
    }
    elsif ( $ans eq 'Q' ) {
        last;
    }
}

sub select_blast_table {
    my ($blast_table, $medium) = @_;
    my $symmetry = queryu("Enter symmetry: 'S', 'C' or 'P', <cr>='S'");
    if ( !$symmetry ) { $symmetry = 'S' }
    my $gamma = get_num( "Enter gamma", 1.4 );
    my $blast_table_old = $blast_table;
    $blast_table =
      Blast::IPS->new( 'symmetry' => $symmetry, 'gamma' => $gamma );
    my $err = $blast_table->get_error();
    if ($err) {
        query("Error: $err; no changes made");
        return ($blast_table_old, $medium);
    }

        $sspd_amb    = 1;
        $p_amb       = 1;                                  # ambient pressure
        $rho_amb     = $gamma * $p_amb / $sspd_amb**2;
        $medium      = {
            _gamma    => $gamma,
            _symmetry => $symmetry,
            _sspd_amb => $sspd_amb,
            _rho_amb  => $rho_amb,
            _p_amb    => $p_amb,
        };
    return ($blast_table, $medium);
}

sub select_variable {

    my ($vname) = @_;

    # key => [ order, text ]
    my %menu = (
        'x'    => [ 1,  'scaled range, = r/d'], 
        'y'    => [ 2,  'overpressure ratio, =(p-p0)/p0'], 
        'w'    => [ 3,  'scaled time of arrival, = c0 t / d'],
        'z'    => [ 4,  'x - w'],
        'X'    => [ 5,  'ln(x)'],
        'Y'    => [ 6,  'ln(y)'],
        'W'    => [ 7,  'ln(w)'],
        'Z'    => [ 8,  'ln(z)'],
        'dYdX' => [ 9,  'dY/dX'],
        'dZdX' => [ 10, 'dZ/dX'],
        'dWdX' => [ 11, 'dW/dX'],
        'm'    => [ 12, '= (S/c0), where S is the shock speed'],
        'q'    => [ 13, '= (c0/S)^2, where S is the shock speed'],
    );

    #my $menu_text = "Select a variable to evaluate:\n";
my $menu_text = <<EOM;
Point evaluation with one dimensionless variable

Let
 r = range; t = time of arrival; p = shock pressure; 
 d = scaled distance (E/p0)^(1/N) 
 E is energy and N = 1,2, or 3 is symmetry
 p0 = initial atmospheric pressure and c0 = sound speed 
 D = shock speed

Variables which may be selected to vary:
EOM
    foreach my $key ( sort { $menu{$a}->[0] <=> $menu{$b}->[0] } keys(%menu) ) {
	$menu_text .= "    $key : $menu{$key}->[1]\n";
    }


    while (1) {
	print $menu_text;
        my $ans = query("Select one of these; <cr>='x':");
	$ans = 'x' unless($ans);
	if (defined($menu{$ans}) ) {
	   return $ans;
	}
	else {
	   hitcr("error, try again");
        }
    }
    return $vname;
}

{

    my ($E0, $p_amb, $sspd_amb, $range, $ground_plane, $symmetry, $blast_table);

    BEGIN {
        $E0           = 1;
        $p_amb        = 1;
        $sspd_amb     = 1;
        $range        = 1;
        $ground_plane = 1;
    }

    sub point_evaluations_with_units {
        my ( $blast_table, $medium ) = @_;

	# perform point evaluations with units
	# internal units are SI but other display units may be used
#  AL ALtimeter reading, m............    0.0 
#  AT Atmospheric Temperature, K......    59.0 
	my ($p_amb, $sspd_amb);
        my $gamma    = $medium->{_gamma};
        my $symmetry = $medium->{_symmetry};
        while (1) {
	    $p_amb=$medium->{_p_amb};
	    $sspd_amb=$medium->{_sspd_amb};
	    my $gtext = $ground_plane ? "on hard surface" : "in free air";
            print <<EOM;
 ----- Dimensional Solution Menu -------
     gamma...........................    $gamma
     symmetry........................    $symmetry
  G  Ground plane option.............    $gtext
  AA Ambient Atmospheric conditions
     Atmospheric Pressure, Pa........    $p_amb 
     ambient sound speed, m/s........    $sspd_amb
  R  Range, m........................    $range
  E  Energy, Joules..................    $E0
   
  RE or C Calculate blast parameters, given: R, E
  RI calculate Energy, given: Range, Impulse
  RP  calculate Energy, given: Range, measured Overpressure
  RT calculate Energy, given: Range, measured TOA
  EP calculate Range, given Energy and OVP
  EI calculate Range, given Energy and IMP
  ET calculate Range, given Energy and TOA
  
  Z Zoom  - view latest results, 1 screen per channel
  L List  - view latest results, 1 line per channel
  F files (read/write)...
  X=eXit  CL=Clear  LG=List Gages LD=List Data

EOM
            my $ans = queryu(":");
	    if ($ans eq 'E') {
		$E0=get_num("Enter energy E0:");
	    }
	    elsif ($ans eq 'R') {
		$range=get_num("Enter range, m:");
	    }
            elsif ( $ans eq 'G' ) {
		print <<EOM;
Select a ground plane option:
  0 = explosion is in free air
  1 = explosion is on a ground plane (rigid surface)
EOM
                $ground_plane = query(":");
            }
	    elsif ($ans eq 'AA' ) {
	    }
	    elsif ($ans eq 'RE' || $ans eq 'C') {
	    }
	    elsif ($ans eq 'RI' ) {
	    }
	    elsif ($ans eq 'RT' ) {
	    }
	    elsif ($ans eq 'RP' ) {
	    }
	    elsif ($ans eq 'EP' ) {
	    }
	    elsif ($ans eq 'ET' ) {
	    }
	    elsif ($ans eq 'EI' ) {
	    }
	    elsif ($ans eq 'X') {
		return;
	    }
	    else {
	    }
        }
    }
}

sub point_evaluations_dimensionless {
    my ($blast_table, $medium, $vname) = @_;
    my $gamma    = $medium->{_gamma};
    my $symmetry = $medium->{_symmetry};
    while (1) {
        my $val = get_num("Enter a value for '$vname', or <cr> to quit:");
	last if $val eq "" || $val !~ /\d/; #^\s*[\-\+\d\.]/;
 
	my ($iQ, $Q);
	if ($vname =~ /^([XYZW]|dYdX|dZdX|dWdX)/) {$Q = $val; $iQ=$vname}
	elsif ($vname =~ /^[xyzw]$/) {$Q = log($val); $iQ=uc($vname)}
        elsif ( $vname eq 'q' ) {

            # Convert q=(c0/D)**2 to Y=ln(ovp ratio)
            $iQ = 'Y';
            my $q = $val;
            if ( $q > 0 && $q < 1) {
                my $ovprat =
                  ( 2 * $gamma - ( $gamma - 1 ) * $q ) /
                  ( ( $gamma + 1 ) * $q ) - 1;

		if ($ovprat<=0) {
		   hitcr("cant take log of '$ovprat'"); next;
		}
                $Q = log($ovprat);
            }
            else {
		hitcr("must have 0< q < 1"); next;
            }
        }
        elsif ( $vname eq 'm' ) {

            # Convert m=(S/c0) to Y=ln(ovp ratio)
            $iQ = 'Y';
            my $m = $val;
            if ( $m > 1 ) {
	        my $q = 1/$m**2;
                my $ovprat =
                  ( 2 * $gamma - ( $gamma - 1 ) * $q ) /
                  ( ( $gamma + 1 ) * $q ) - 1;

		if ($ovprat<=0) {
		   hitcr("cant take log of '$ovprat'"); next;
		}

                $Q = log($ovprat);
            }
            else {
		hitcr("m must be > 1"); next;
            }
        }
	else {
		die "coding incomplete for variable '$vname'";
        }
        my $ret  = $blast_table->wavefront( $iQ => $Q );
        my $X    = $ret->{X};
        my $Y    = $ret->{Y};
        my $Z    = $ret->{Z};
        my $dYdX    = $ret->{dYdX};
        my $dZdX    = $ret->{dZdX};
        my $Tpos = $ret->{Tpos};
        my $Lpos = $ret->{Lpos};
        my $Tneg = $ret->{Tneg};
        my $Lneg = $ret->{Lneg};
        my $Ixr_pos = $ret->{Ixr_pos};
        my $Ixr_neg = $ret->{Ixr_neg};
        my $x    = exp($X);
        my $y    = exp($Y);
        my $z    = exp($Z);
        my $w    = $x - $z;
        my $W    = $w > 0 ? log($w) : -999;
        my $term = $y * ( $gamma + 1 ) / ( 2 * $gamma );
        my $m    = sqrt( 1 + $term );
        my $q    = 1 / $m**2;
        my $up   = $y / ( $gamma * $m );
	foreach ($x, $y, $z, $w, $X, $Y, $Z, $W, $dYdX, $dZdX) {
           $_=sprintf("%0.8g", $_);
	}
	foreach ($Tpos, $Lpos, $Tneg, $Lneg, $m, $q, $up, $Ixr_pos, $Ixr_neg) {
           $_=sprintf("%0.6g", $_);
	}

	print <<EOM;

Results; all variables are dimensionless:

x    = $x = scaled range r/d;     ln(x) = X = $X 
y    = $y = (P-P0)/P0;            ln(y) = Y = $Y 
z    = $z = (r-c0 t)/d;           ln(z) = Z = $Z 
toa  = $w = scaled toa;           ln(t) = W = $W 
dYdX = $dYdX
T+ = $Tpos = time duration of positive phase 
L+ = $Lpos = length of positive phase 
T- = $Tneg = time duration of negative phase 
L- = $Lneg = length of negative phase 
m  = $m = S/c0=shock Mach number
q  = $q = 1/m^2
up = $up = (shock particle velocity/c0)
I+ = $Ixr_pos = limiting positive impulse x r
I- = $Ixr_neg = limiting negative impulse x r
EOM

    }
}

sub table_operations {
    my ($blast_table, $medium) = @_;

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
	    my $Id_max=4;
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

        #print "BOOGALOO: closed\n";
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
    query( $msg . ". hit <cr>" );
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
  my ($msg, $default)=@_;
    my $count = 0;
  ASK:
    my $ans   = query($msg);
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
