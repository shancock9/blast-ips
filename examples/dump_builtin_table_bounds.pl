#!/usr/bin/perl
use strict;
use warnings;
use Blast::IPS;

# Print a list of the builtin tables and their bounds
my $rindex = Blast::IPS->get_table_index();

my @output;
foreach my $key ( keys %{$rindex} ) {
    my $item           = $rindex->{$key};
    my $table_name     = $key;
    my $gamma          = $item->{gamma};
    my $symmetry       = $item->{symmetry};
    my $error_estimate = $item->{error_estimate};

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

    # my $alpha   = $blast_table->get_alpha();
    my $rbounds = $blast_table->get_table_bounds();
    my ( $Xu, $Yu, $dYdXu, $Zu, $dZdXu ) = @{ $rbounds->[0] };
    my ( $Xl, $Yl, $dYdXl, $Zl, $dZdXl ) = @{ $rbounds->[1] };
    push @output,
      [
        $symmetry, $gamma, $table_name, $error_estimate, $Xu, $Yu, $dYdXu, $Xl,
        $Yl, $dYdXl
      ];
}

# sort by symmetry and gamma
my @sorted = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @output;

print
  "symmetry\tgamma\ttable_name\terror_estimate\tXu\tYu\tdYdXu\tXl\tYl\tdYdXl\n";
foreach my $item (@sorted) {
    my $str = join "\t", @{$item};
    print "$str\n";
}
