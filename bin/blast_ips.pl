#!/usr/bin/perl 
use warnings;
use strict;

# This is a driver to illustrate usage of Blast::IPS.
use Blast::IPS;

# Setup a default blast table
my $gamma       = 1.4;
my $symmetry    = 'S';
my $blast_table = Blast::IPS->new( symmetry => $symmetry, gamma => $gamma );

my %symmetry_name = (
    0 => 'Plane',
    1 => 'Cylindrical',
    2 => 'Spherical',
);

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
    }
    elsif ( $ans eq 'P' ) {
        my $vname = 'X';
        $vname = select_variable($vname);
        point_evaluations($blast_table, $vname, $gamma, $symmetry);
    }
    elsif ( $ans eq 'T' ) {
        table_operations($blast_table);
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

#my $ans =
#  queryu("Enter symmetry: S=spherical, C=cylindrical, P=plane; <cr>='S':");
#my $symmetry = $ans;
#if ( $ans !~ /^[012]$/ ) {
#    $symmetry = ( $ans =~ /^P/i ? 0 : $ans =~ /^C/i ? 1 : 2 );
#}
#
#my $gamma = get_num("Enter gamma; <cr>=1.4:");
#if ( !$gamma ) { $gamma = 1.4 }
#
## Create a blast object for this case
#my $blast_table = Blast::IPS->new( symmetry => $symmetry, gamma => $gamma );
#my $err=$blast_table->get_error();
#if ($err) {
#    print "Exiting due to error: $err\n";
#    exit;
#}
#
##my $alpha = $blast_table->get_alpha();
##my $alpha_i = Blast::IPS::alpha_interpolate($symmetry,$gamma);
##query("alpha=$alpha, alpha_i=$alpha_i");
#
#my %symmetry_name = (
#    0 => 'Plane',
#    1 => 'Cylindrical',
#    2 => 'Spherical',
#);
#
##my $iQ = query("Enter 'X' to enter range, 'Y' to enter overpressure");
#
#print <<EOM;
#Blast form a point source in an ideal homogeneous atmosphere
#Symmetry=$symmetry_name{$symmetry}
#Gamma=$gamma
#EOM

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
        's'    => [ 12, '= (D/c0), where D is the shock speed'],
        'q'    => [ 13, '= (c0/D)^2, where D is the shock speed'],
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
    my ($blast_table, $vname, $gamma, $symmetry) = @_;
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
    my ($blast_table) = @_;
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
