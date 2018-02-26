use strict;
use Test;

BEGIN { plan tests => 2 }

# Create a point blast table for a spherical source with gamma=1.4
use Blast::IPS;
my %args = ( 'ASYM' => 2, 'gamma' => 1.4 );
my $blast_table = Blast::IPS->new( \%args );

# A test point: the overpressure ratio at a scaled range of 1
my $lambda_t = 1;
my $ovprat_t   = 0.4980896;

# Lookup the overpressure at a given scaled range 
my $iQ          = 'X';
my $Q           = log($lambda_t);
my $ret         = $blast_table->lookup( $Q, $iQ );
my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};
my $ovprat = exp($Y);
my $err = abs( $ovprat - $ovprat_t );

print STDERR "overpressure ratio=$ovprat at scaled range lambda=$lambda_t; error=$err\n";
ok( $err <= 1.e-5 );

# Lookup the scaled range to a given overpressure ratio
$iQ          = 'Y';
$Q           = log($ovprat_t);
$ret         = $blast_table->lookup( $Q, $iQ );
( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};
my $lambda = exp($X);
$err = abs($lambda-$lambda_t)/$lambda_t;

print STDERR "lambda=$lambda at overprressure ratio $ovprat_t; error=$err\n";
ok( $err <= 1.e-5 );

