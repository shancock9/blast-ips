#!/usr/bin/perl 
use warnings;
use strict;

# This is a simple driver to illustrate usage of Blast::IPS.
# It will ask for a symmetry and gamma, 
# then go into a loop to evaluate the shock for any range
use Blast::IPS;

my $ans =
  queryu("Enter symmetry: S=spherical, C=cylindrical, P=plane; <cr>='S':");
my $symmetry = $ans;
if ( $ans !~ /^[012]$/ ) {
    $symmetry = ( $ans =~ /^P/i ? 0 : $ans =~ /^C/i ? 1 : 2 );
}

my $gamma = get_num("Enter gamma; <cr>=1.4:");
if ( !$gamma ) { $gamma = 1.4 }

# Create a blast object for this case
my $blast_table = Blast::IPS->new( symmetry => $symmetry, gamma => $gamma );
my $err=$blast_table->get_error();
if ($err) {
    print "Exiting due to error: $err\n";
    exit;
}

#my $alpha = $blast_table->get_alpha();
#my $alpha_i = Blast::IPS::alpha_interpolate($symmetry,$gamma);
#query("alpha=$alpha, alpha_i=$alpha_i");

my %symmetry_name = (
    0 => 'Plane',
    1 => 'Cylindrical',
    2 => 'Spherical',
);

#my $iQ = query("Enter 'X' to enter range, 'Y' to enter overpressure");

print <<EOM;
Blast form a point source in an ideal homogeneous atmosphere
Symmetry=$symmetry_name{$symmetry}
Gamma=$gamma
EOM

while (1) {
    my $lambda = query("Enter a scaled range, or <cr> to quit:");
    last unless ( $lambda && $lambda > 0 );
    my $iQ  = 'X';
    my $Q   = log($lambda);
    my $ret = $blast_table->wavefront( $iQ => $Q );
    #my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};
    my $Y=$ret->{Y};
    my $ovprat = exp($Y);
    print "The overpressure ratio is: $ovprat\n";
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
