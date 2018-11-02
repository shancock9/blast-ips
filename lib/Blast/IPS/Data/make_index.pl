#!/usr/bin/perl
use strict;

# TODO:
# Handle multiple data files; go by date

my @files = glob('*.pm');
my $rtmp  = [];
foreach my $file (@files) {
    next if ( $file eq 'Index.pm' );
    if ( $file =~ /^([SCP])(.*)(\.pm)$/ ) {
        my $module_name = $1 . $2;
        my $symmetry    = ( $1 eq 'S' ) ? 2 : ( $1 eq 'C' ) ? 1 : 0;
        my $N           = 0;
        my $gamma       = $2;

        # allow 'x' instead of '.' in file name
        $gamma =~ s/x/\./;

        # allow gamma to be preceded by _G
        # in this case, the number N of the FD run will be before the _
        if ( $gamma =~ /^(\d+)_G(.*)/ ) { $N = $1; $gamma = $2 }

        # allow trailing text
        if ( $gamma =~ /([^_]+)_/ ) { $gamma = $1 }

        # we should have a valid gamma
        if ( !is_float($gamma) ) {
            print STDERR "skipping $file: cannot parse a gamma value\n";
            next;
        }
        if ( $gamma <=1 ) {
            print STDERR "skipping $file: got invalid gamma='$gamma\n";
            next;
        }

        my $key = _make_table_name( $symmetry, $gamma );
        push @{ $rtmp->[$symmetry] }, [ $gamma, $key, $module_name ];
    }
    else {
        print STDERR "skipping $file: cannot parse symmetry from name\n";
        next;
    }
}

# Construct a sorted table of modules for a given symmetry and gamma:
#    $rindex_table->[symmetry] = [ [gamma, module], [gamma, module], .. ]
# This can be used by calling routines to see what is available and for
# loading tables for doing interpolation in gamma.
my $rindex_table;
foreach my $sym ( 0, 1, 2 ) {
    my @sorted = sort { $a->[0] <=> $b->[0] } @{ $rtmp->[$sym] };
    $rindex_table->[$sym] = \@sorted;

##   Reference code for selecting the best table, if necessary:
##    (would have to be modified)

##        # We require a table of unique gamma values.
##        # If there are multiple tables for a given gamma, use the
##        # table with least error.
##        my @unique;
##        my ( $gamma_last, $key_last );
##        foreach my $item (@sorted) {
##            my ( $gamma, $key ) = @{$item};
##            if ( !$gamma_last || $gamma != $gamma_last ) {
##                push @unique, $item;
##                $gamma_last = $gamma;
##                $key_last   = $key;
##                next;
##            }
##            my $err_last = $rshock_tables_info->{$key_last}->[2];
##            my $err      = $rshock_tables_info->{$key}->[2];
##            next if ( $err_last < $err );
##            $unique[-1] = $item;
##        }
##
##        $rgamma_table->[$sym] = \@unique;
##    $rindex_table->[$sym] = \@unique;

}

use Data::Dumper;
my $string=Data::Dumper->Dump( [$rindex_table],['rindex_table'] );
my $module=make_index_module($string);
$module = tidy($module);
my $findex = 'Index.pm';
open(my $fd, '>', $findex) || die "cannot open $findex: $!\n";
$fd->print($module);
$fd->close();
print "wrote $findex\n";

sub make_index_module {
    my ($code) = @_;
    my $header = <<'EOM';
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
    
EOM
    my $trailer = <<'EOM';
}

sub get_index {
   return $rindex_table;
}
1;
EOM
    return $header . $code . $trailer;
}

sub tidy {
    my ($source) = @_;
    my $errorfile_string;
    my $stderr_string;
    my $output;
    use Perl::Tidy;
    my $err = Perl::Tidy::perltidy(
        source      => \$source,
        destination => \$output,
        perltidyrc  => '',
        argv        => '',             # for safety; hide any ARGV from perltidy
        stderr      => \$stderr_string,
        errorfile => \$errorfile_string,    # not used when -se flag is set
    );
    if ( $err || $stderr_string || $errorfile_string ) {
        if ($err) {
            die <<EOM;
This error received calling Perl::Tidy
$err;
EOM
        }
        if ($stderr_string) {
            print STDERR "---------------------\n";
            print STDERR "<<STDERR>>\n$stderr_string\n";
            print STDERR "---------------------\n";
            print STDERR "This error received calling Perl::Tidy\n";
            die;
        }
        if ($errorfile_string) {
            print STDERR "---------------------\n";
            print STDERR "<<.ERR file>>\n$errorfile_string\n";
            print STDERR "---------------------\n";
            print STDERR "This error received calling Perl::Tidy:\n"; 
            print STDERR $errorfile_string;
            die;
        }
    }
    return $output;
}

sub is_float {
    defined( $_[0] )
      && $_[0] =~ /^\s*([+-]?)(?=\d|\.\d)\d*(\.\d*)?([EeDd]([+-]?\d+))?\s*$/;
}

sub _make_table_name {

    # Make a standard table name by combining the symmetry and gamma;
    # i.e. make the name 'S1.37' for spherical symmetry with gamma=1.37
    # This name will be used as the hash key for a module in the global table
    my ( $sym, $gamma ) = @_;
    my $table_name;
    my $gamma_pr = sprintf "%.4g", $gamma;
    if ( $sym =~ /\d/ ) {
        $sym = $sym == 0 ? 'P' : $sym == 1 ? 'C' : 'S';
    }
    $table_name = $sym . $gamma_pr;
    return $table_name;
}
