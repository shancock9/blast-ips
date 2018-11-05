#!/usr/bin/perl
use strict;
use warnings;
use Blast::IPS;

# Print a list of the builtin tables and their bounds
my $rindex = Blast::IPS::get_index();

my @output;
foreach my $symmetry ( 0 .. 2 ) {
    foreach my $item ( @{ $rindex->[$symmetry] } ) {
        my ( $gamma, $table_name ) = @{$item};

        # Create a table for this case
        my $blast_table = Blast::IPS->new( 'table_name' => $table_name );

        # do some checks
        if ( !defined($blast_table) ) {
            die "missing table for sym=$symmetry, gamma=$gamma\n";
        }
        my $table_name_check = $blast_table->get_table_name();
        if ( $table_name ne $table_name_check ) {
            die <<EOM;
Asked for table '$table_name' but got '$table_name_check'
gamma=$gamma, symmetry=$symmetry
EOM
        }

        my $rinfo = $blast_table->get_info();
##use Data::Dumper;
##print Data::Dumper->Dump([$rinfo]);
##exit 1;
        my $error_estimate = $rinfo->{Max_Error};

	# Patch until get_info has shock data
	$error_estimate=1.e-6 unless ($error_estimate);

        # my $alpha   = $blast_table->get_alpha();
        my $rbounds = $blast_table->get_table_bounds();
        my ( $Xu, $Yu, $dYdXu, $Zu, $dZdXu ) = @{ $rbounds->[0] };
        my ( $Xl, $Yl, $dYdXl, $Zl, $dZdXl ) = @{ $rbounds->[1] };
        push @output,
          [
            $symmetry, $gamma, $table_name, $error_estimate, $Xu, $Yu, $dYdXu,
            $Xl, $Yl, $dYdXl
          ];
    }
}

# sort by symmetry and gamma
my @sorted = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @output;

print
  "symmetry\tgamma\ttable_name\terror_estimate\tXu\tYu\tdYdXu\tXl\tYl\tdYdXl\n";
foreach my $item (@sorted) {
    my $str = join "\t", @{$item};
    print "$str\n";
}
