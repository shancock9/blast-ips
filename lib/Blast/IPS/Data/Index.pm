package Blast::IPS::Data::Index;

########################################################################
# This file is created by 'Blast::IPS::Data::make_index.pl'.  It should be
# regenerated whenever there is a change in the list of data modules.
########################################################################
use strict;
use warnings;
use 5.006;
our $VERSION = 0.1.1;
use Carp;
my $rindex_table;

BEGIN {

   # This is a sorted table of available modules for a given symmetry and gamma:
   #    $rindex_table->[symmetry] = [
   #       [gamma, key, module],
   #       [gamma, key, module],
   #       ....
   #    ];
   # where
   #    symmetry = 0,1,2 for plane, cylindrical, spherical
   #    gamma    = the ideal gas gamma
   #    key      = a standard hash key for this case
   #    module   = the base name of the module file in Blast/IPS/Data/

    $rindex_table = [
        [
            [ '1.1',   'P1x1' ],
            [ '1.12',  'P1x12' ],
            [ '1.15',  'P1x15' ],
            [ '1.17',  'P1x17' ],
            [ '1.2',   'P1x2' ],
            [ '1.23',  'P1x23' ],
            [ '1.25',  'P1x25' ],
            [ '1.3',   'P1x3' ],
            [ '1.35',  'P1x35' ],
            [ '1.4',   'P1x4' ],
            [ '1.45',  'P1x45' ],
            [ '1.5',   'P1x5' ],
            [ '1.55',  'P1x55' ],
            [ '1.6',   'P1x6' ],
            [ '1.667', 'P1x667' ],
            [ '1.7',   'P1x7' ],
            [ '1.8',   'P1x8' ],
            [ '1.9',   'P1x9' ],
            [ 2,       'P2' ],
            [ '2.2',   'P2x2' ],
            [ '2.5',   'P2x5' ],
            [ 3,       'P3' ],
            [ '3.5',   'P3x5' ],
            [ 4,       'P4' ],
            [ '4.5',   'P4x5' ],
            [ 5,       'P5' ],
            [ '5.5',   'P5x5' ],
            [ 6,       'P6' ],
            [ '6.5',   'P6x5' ],
            [ 7,       'P7' ]
        ],
        [
            [ '1.1',   'C1x1' ],
            [ '1.12',  'C1x12' ],
            [ '1.15',  'C1x15' ],
            [ '1.17',  'C1x17' ],
            [ '1.2',   'C1x2' ],
            [ '1.23',  'C1x23' ],
            [ '1.25',  'C1x25' ],
            [ '1.3',   'C1x3' ],
            [ '1.35',  'C1x35' ],
            [ '1.4',   'C1x4' ],
            [ '1.45',  'C1x45' ],
            [ '1.49',  'C1x49' ],
            [ '1.55',  'C1x55' ],
            [ '1.6',   'C1x6' ],
            [ '1.667', 'C1x667' ],
            [ '1.7',   'C1x7' ],
            [ '1.8',   'C1x8' ],
            [ '1.9',   'C1x9' ],
            [ 2,       'C2' ],
            [ '2.2',   'C2x2' ],
            [ '2.5',   'C2x5' ],
            [ 3,       'C3' ],
            [ '3.5',   'C3x5' ],
            [ 4,       'C4' ],
            [ '4.5',   'C4x5' ],
            [ 5,       'C5' ],
            [ '5.5',   'C5x5' ],
            [ 6,       'C6' ],
            [ '6.5',   'C6x5' ],
            [ 7,       'C7' ]
        ],
        [
            [ '1.1',   'S1x1' ],
            [ '1.12',  'S1x12' ],
            [ '1.15',  'S1x15' ],
            [ '1.17',  'S1x17' ],
            [ '1.2',   'S1x2' ],
            [ '1.23',  'S1x23' ],
            [ '1.25',  'S1x25' ],
            [ '1.3',   'S1x3' ],
            [ '1.35',  'S1x35' ],
            [ '1.4',   'S1x4' ],
            [ '1.45',  'S1x45' ],
            [ '1.5',   'S1x5' ],
            [ '1.55',  'S1x55' ],
            [ '1.6',   'S1x6' ],
            [ '1.667', 'S1x667' ],
            [ '1.7',   'S1x7' ],
            [ '1.8',   'S1x8' ],
            [ '1.9',   'S1x9' ],
            [ 2,       'S2' ],
            [ '2.2',   'S2x2' ],
            [ '2.5',   'S2x5' ],
            [ 3,       'S3' ],
            [ '3.5',   'S3x5' ],
            [ 4,       'S4' ],
            [ '4.5',   'S4x5' ],
            [ 5,       'S5' ],
            [ '5.5',   'S5x5' ],
            [ 6,       'S6' ],
            [ '6.5',   'S6x5' ]
        ]
    ];
}

sub get_index {
    return $rindex_table;
}
1;
