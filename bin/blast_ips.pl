#!/usr/bin/perl 
use warnings;
use strict;

# This is a simple driver for Blast::IPS
use Blast::IPS;

my $ans =
  queryu("Enter symmetry: S=spherical, C=cylindrical, P=plane; <cr>='S':");
my $symmetry = $ans;
if ( $ans !~ /^[012]$/ ) {
    $symmetry = ( $ans =~ /^P/i ? 0 : $ans =~ /^C/i ? 1 : 2 );
}

my $gamma = get_num("Enter gamma; <cr>=1.4:");
if ( !$gamma ) { $gamma = 1.4 }

# Create a blast object with a table of values 
my $blast_table = Blast::IPS->new( symmetry => $symmetry, gamma => $gamma );

my %symmetry_name = (
    0 => 'Plane',
    1 => 'Cylindrical',
    2 => 'Spherical',
);

print <<EOM;
Blast form a point source in an ideal atmosphere
Symmetry=$symmetry_name{$symmetry}
Gamma=$gamma
EOM

while (1) {
    my $lambda = query("Enter a scaled range, or <cr> to quit:");
    last unless ( $lambda && $lambda > 0 );
    my $iQ  = 'X';
    my $Q   = log($lambda);
    my $ret = $blast_table->lookup( $Q, $iQ );
    my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};
    my $ovprat = exp($Y);
    print "Overpressure ratio=$ovprat\n";
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
