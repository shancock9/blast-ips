#!/usr/bin/perl 
use warnings;
use strict;

# This is a driver to illustrate usage of Blast::IPS.
use Blast::IPS;

my $audit_string = "";

# Setup a default blast table
my $gamma       = 1.4;
my $symmetry    = 'S';
my $blast_table = Blast::IPS->new( symmetry => $symmetry, gamma => $gamma );

my %symmetry_name = (
    0 => 'Plane',
    1 => 'Cylindrical',
    2 => 'Spherical',
);

my $sspd_amb = 1;
my $p_amb    = 1;                                  # ambient pressure
my $rho_amb  = $gamma * $p_amb / $sspd_amb**2;
my $medium   = {
    _gamma    => $gamma,
    _sspd_amb => $sspd_amb,
    _rho_amb  => $rho_amb,
    _p_amb    => $p_amb,
};

# main loop
while (1) {
    my $table_name = $blast_table->get_table_name();
    $symmetry = $blast_table->get_symmetry();
    $gamma    = $blast_table->get_gamma();
    print <<EOM;
Point Source Explosion in Ideal Gas 

  Symmetry=$symmetry_name{$symmetry},  Gamma=$gamma

Enter one of the following:
  N             - New Symmetry and/or Gamma
  P		- Point evaluations using current Table ..
  T		- Table operations ..
  q		- quit
EOM
    my $ans = queryu(":");
    if ( $ans eq 'N' ) {
        $blast_table = select_blast_table($blast_table);

        $sspd_amb    = 1;
        $p_amb       = 1;                                  # ambient pressure
        $rho_amb     = $gamma * $p_amb / $sspd_amb**2;
        $medium      = {
            _gamma    => $gamma,
            _sspd_amb => $sspd_amb,
            _rho_amb  => $rho_amb,
            _p_amb    => $p_amb,
        };
    }
    elsif ( $ans eq 'P' ) {
        my $vname = 'X';
        $vname = select_variable($vname);
        point_evaluations($blast_table, $medium, $vname, $gamma, $symmetry);
    }
    elsif ( $ans eq 'T' ) {
        table_operations($blast_table, $medium);
    }
    elsif ( $ans eq 'q' ) {
        exit;
    }
}

sub select_blast_table {
    my ($blast_table) = @_;
    my $symmetry = queryu("Enter symmetry: 'S', 'C' or 'P', <cr>='S'");
    if ( !$symmetry ) { $symmetry = 'S' }
    my $gamma = get_num( "Enter gamma", 1.4 );
    $blast_table =
      Blast::IPS->new( 'symmetry' => $symmetry, 'gamma' => $gamma );
    my $err = $blast_table->get_error();
    if ($err) {
        print "Error: $err\n";
        return;
    }
    return $blast_table;
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

    my $menu_text = "Select a variable to evaluate:\n";
    foreach my $key ( sort { $menu{$a}->[0] <=> $menu{$b}->[0] } keys(%menu) ) {
	$menu_text .= "    $key : $menu{$key}->[1]\n";
    }
    $menu_text .= <<EOM;

where
 r = range; t = time of arrival; p = shock pressure; 
 d = scaled distance (E/p0)^(1/N) 
 E is energy and N = 1,2, or 3 is symmetry
 p0 = initial atmospheric pressure and c0 = sound speed 
 D = shock speed
EOM

    while (1) {
	print $menu_text;
        my $ans = query(":");
	if (defined($menu{$ans}) ) {
	   return $ans;
	}
	else {
	   query("error, try again");
        }
    }
    return $vname;
}

sub point_evaluations {
    my ($blast_table, $vname, $medium, $gamma, $symmetry) = @_;
    while (1) {
        my $val = get_num("Enter a value for $vname, or 'q' to quit:");
        last if ( $val eq 'q' );
 
	my ($iQ, $Q);
	if ($vname =~ /^([XYZW]|dYdX|dZdX|dWdX)/) {$Q = $val; $iQ=$vname}
	elsif ($vname =~ /^[xyzw]$/) {$Q = log($val); $iQ=uc($vname)}
        elsif ( $vname eq 'q' ) {

            # Convert q=(c0/D)**2 to Y=ln(ovp ratio)
            $iQ = 'Y';
            my $q = $val;
            if ( $q > 0 ) {
                my $ovprat =
                  ( 2 * $gamma - ( $gamma - 1 ) * $q ) /
                  ( ( $gamma + 1 ) * $q ) - 1;

                $Q = log($ovprat);
            }
            else {
                $Q = 20;
            }
        }
        elsif ( $vname eq 's' ) {

            # Convert s=(D/c0) to Y=ln(ovp ratio)
            $iQ = 'Y';
            my $s = $val;
            if ( $s > 0 ) {
	        my $q = 1/$s**2;
                my $ovprat =
                  ( 2 * $gamma - ( $gamma - 1 ) * $q ) /
                  ( ( $gamma + 1 ) * $q ) - 1;

                $Q = log($ovprat);
            }
            else {
                $Q = 20;
            }
        }
	else {
		die "coding incomplete for $vname";
        }
        my $ret = $blast_table->wavefront( $iQ => $Q );
        my $X      = $ret->{X};
        my $Y      = $ret->{Y};
        my $Z      = $ret->{Z};
        my $x      = exp($X);
        my $y      = exp($Y);
        my $z      = exp($Z);
        my $w    = $x - $z;
        my $W    = log($w);

	print <<EOM;
X=$X   x=$x (=scaled range, r/d)
Y=$Y   y=$y (=overpressure ratio)
Z=$Z   z=$z 
W=$W   w=$w (=scaled toa, c0 t / d)
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
  q		- eXit
EOM

        my $ans = queryu(":");
        if ( $ans eq 'PLOT' ) {

        #            my $npts = get_num( "Total number of points",   400 );
        #            my $x1   = get_num( "Starting value of ln(r):", -5 );
        #            my $x2   = get_num( "Ending value of ln(r):",   15 );
        #            $npts = 100 unless ($npts);
        #            my $rsolution = $blast_table->table_gen( $npts, $x1, $x2 );

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
            query("Hit <cr>");
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
        elsif ( $ans eq 'q' ) {
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
