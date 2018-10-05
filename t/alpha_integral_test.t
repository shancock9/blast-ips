#!/usr/bin/perl
use strict;
use warnings;
use Test;
use Blast::IPS::SimilaritySolution;
my $rtest_cases; 

BEGIN {

    # We test the integration routine which evaluates alpha

    # The alpha values in the following table were computed with
    # the sedov3.f program, downloaded from:
    # http://www.cococubed.com/code_pages/sedov.shtml

    # [symmetry, gamma, $alpha]
    $rtest_cases = [
        [ 0, 1.4,   0.53874280 ],
        [ 1, 1.4,   0.984074056 ],
        [ 2, 1.4,   0.85107188 ],
        [ 2, 1.667, 0.493319042 ],
    ];

    my $ntests = @{$rtest_cases};
    plan tests => $ntests;
}

my $VERBOSE = $ARGV[0];

# Settings..

# itmax = max iterations in case the tolerance is not reached;
my $itmax = 15;

# tol = stopping tolerance on absolute value of alpha
my $tol = 1.e-13;

my $alphatol = 5.e-8;
for ( my $icase = 0 ; $icase < @{$rtest_cases} ; $icase++ ) {
    my $case = $rtest_cases->[$icase];
    my ( $symmetry, $gamma, $alpha_ref ) = @{$case};

    my $obj = Blast::IPS::SimilaritySolution->new(
        symmetry => $symmetry,
        gamma    => $gamma
    );
    my ( $alpha, $err, $it, $alpha_trap, $err0 ) =
      $obj->alpha_integral( $tol, $itmax );
    my $alpha_err = abs( $alpha_ref - $alpha ) / $alpha_ref;
    if ( !ok( $alpha_err < $alphatol )
        || $VERBOSE )
    {
        print STDERR
"sym=$symmetry\tgamma=$gamma\talpha_ref=$alpha_ref\talpha=$alpha\t$err\tdiff=$alpha_err\n";
    }

}
