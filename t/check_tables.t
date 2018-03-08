use strict;
use warnings;
use Test;
use Blast::IPS;

my $rtables;

BEGIN {

    # use the builtin table check routine to be sure the tables are okay
    plan tests => 1;
}

my $err = Blast::IPS::check_tables();
if ($err) {
    print STDERR $err;
}

ok( !$err );
