#!/usr/bin/perl
use strict;
use warnings;
use Test;
use Blast::IPS;
my $rindex = Blast::IPS->get_table_index();
print "table_name\tsymmetry\tgamma\terr\tN\n";
foreach my $table_name ( sort keys %{$rindex} ) {
    my $symmetry = $rindex->{$table_name}->{symmetry};
    my $gamma    = $rindex->{$table_name}->{gamma};
    my $err      = $rindex->{$table_name}->{error_estimate};
    my $N        = $rindex->{$table_name}->{N};
    print "$table_name\t$symmetry\t$gamma\t$err\t$N\n";
}
