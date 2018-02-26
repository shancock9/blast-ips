#!/usr/bin/perl 
use warnings;
use strict;

# This is a simple driver for Blast::PointSource
use Blast::IPS;

my $ASYM=2; 
my $gamma=1.4;
my %args = ( 'ASYM' => $ASYM, 'gamma' => $gamma );
my $blast_table = Blast::PointSource->new( \%args );
my %symmetry_name = (
    0 => 'Plane',
    1 => 'Cylindrical',
    2 => 'Spherical',
);

print <<EOM;
Blast form a point source in an ideal atmosphere
Symmetry=$symmetry_name{$ASYM}
Gamma=$gamma
EOM

while (1) {
   my $lambda=query("Enter a scaled range, or <cr> to quit:");
   last unless ($lambda && $lambda>0);
   my $iQ          = 'X';
   my $Q           = log($lambda);
   my $ret         = $blast_table->lookup( $Q, $iQ );
   my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};
   my $ovprat = exp($Y);
   print "Overpressure ratio=$ovprat\n";
}

sub query {
    my ($msg) = @_;
    print $msg;
    my $ans = <STDIN>;
    chomp $ans;
    return $ans;
}
