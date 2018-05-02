#!/usr/bin/perl

# Comparison of Blast::IPS with Korobeinikov small q theory as outlined in his book.
# His A1 values are in good agreement with my calculations
# The time values have about 3 times the error of the space values
# The space values are accurate to 1% to around q=0.1 to 0.2
# The maximum overall space errors are a little over 4%
# In this test we just check the space errors and accept results below 5% error

package main;
use strict;
use warnings;
use Test;
my $audit_string = ""; #audit_string('#');
use Blast::IPS;
my @gamma_list = ( 1.1, 1.2, 1.3, 1.4, 1.667, 2, 3 );
my @symmetry_list = ( 0, 1, 2 );
my @output;
my $verbose = $ARGV[0];
my $TOL = 0.05;

# Table of A1 values from Korob book near p135
#         gamma=1.1, 1.2, 1.3, 1.4, 1.667, 2, 3
# symmetry=0
# symmetry=1
# symmetry=2
my $rA1 = [
    [ 2.3257, 2.2437, 2.1862, 2.1433, 2.0683, 2.0143, 1.9407 ],
    [ 2.0866, 2.0424, 2.0092, 1.9836, 1.9374, 1.9043, 1.8632 ],
    [ 2.0010, 1.9666, 1.9396, 1.9182, 1.8785, 1.8498, 1.8141 ],
];

BEGIN {
    plan tests => 21;
}

foreach my $symmetry (@symmetry_list) {
    foreach my $igamma ( 0, 1, 2, 3, 4, 5, 6 ) {
        my $gamma = $gamma_list[$igamma];
        my $A1    = $rA1->[$symmetry]->[$igamma];
        my %hash  = ( 'symmetry' => $symmetry, 'gamma' => $gamma );
        push @output, "\nsym=$symmetry, gamma=$gamma, A1=$A1\n";
        my $rho_amb = $gamma;
        my $E0      = 1;

        my $blast_table = Blast::IPS->new( \%hash );
        my $error       = $blast_table->get_error();
        if ($error) { print "# error making blast table:\n$error"; next }

        my $err_R_max=0;
        for ( my $q = 0.009 ; $q <= 0.2 ; $q += 0.01 ) {

            # Convert q=(c0/D)**2 to Y1=ovp ratio for comparison with tables in
            # Korobeinikov
            if ( $q < 0.01 ) { $q = 0.001 }    ## One small value
            my $ovprat =
              ( 2 * $gamma - ( $gamma - 1 ) * $q ) / ( ( $gamma + 1 ) * $q ) -
              1;
            my $YY = log($ovprat);
            my $ret = $blast_table->lookup( $YY, 'Y' );
            my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};
            my $ovp    = exp($Y);
            my $z      = exp($Z);
            my $lambda = exp($X);
            my $Toa    = $lambda - $z;

            my $alpha = Blast::IPS::alpha_interpolate( $symmetry, $gamma );
            my $rho0  = $rho_amb;
            my $Eoad  = $E0 / ( $rho0 * $alpha );

            # NOTE: we are matching with the same q (or D) and not ovp
            #my $ret=$sedov_obj->shock_from_P($ovp+1);
            my $D = 1 / sqrt($q);
            my $retD = shock_from_D( $symmetry, $gamma, $rho0, $Eoad, $D );
            my ( $R_0, $T_0, $P_0, $U_0, $rho_0, $D_0, ) = @{$retD};

            my $delta = 2 / ( $symmetry + 3 );
            my $pow   = ( $symmetry + 1 );
            my $Rs    = ( $delta**2 * $Eoad / $D**2 )**( 1 / $pow );
            my $Ts    = $Rs**( 1 / $delta ) / sqrt($Eoad);

            my $Rs_cubed = $delta**2 * $q / ( $gamma * $alpha );
            $Rs = $Rs_cubed**( 1 / ( $symmetry + 1 ) );

            # Note the missing factor gamma in sqrt because of difference
            # in scaling
            $Ts = sqrt($q) * $delta * $Rs;

            ############################################
            # Base these on values returned from Sedov
            # small q theory:
            ############################################
            my $R_0_pow = $R_0**$pow;
            my $R_q_pow = ( $R_0_pow * exp( $A1 * $q ) );
            my $A1_r    = log( $lambda**$pow / $R_0_pow ) / $q;
            my $R_q     = $R_q_pow**( 1 / $pow );

            #my $R_q=($R_0_pow*exp($A1*$q))**(1/$pow);
            my $err_R = ( $R_q - $lambda ) / $lambda;
	    my $abs_err_R = abs($err_R);
	    if ($abs_err_R>$err_R_max) {$err_R_max=$abs_err_R}

            my $nu    = $symmetry + 1;
            my $top   = ( $nu * $delta + 2 ) * $A1 * $q;
            my $bot   = 2 * $nu * $delta * ( $nu * $delta + 1 );
            my $T_q   = $T_0 * ( 1 + $top / $bot );
            my $err_T = ( $T_q - $Toa ) / $Toa;
            my $A1_t =
              ( $Toa / $T_0 - 1 ) * $bot / ( ( $nu * $delta + 2 ) * $q );

            $lambda = sprintf "%0.6f", $lambda;
            $Toa    = sprintf "%0.6f", $Toa;
            push @output,
"$q\t$T_0\t$R_0\t$Toa\t$lambda\t$Ts\t$Rs\t$T_q\t$R_q\t$err_T\t$err_R\t$A1\t$A1_r\t$A1_t\n";
        }
	ok($err_R_max<$TOL) || print STDERR "sym=$symmetry, gam=$gamma, err_R_max=$err_R_max\n";
    }
}
if ($verbose) {
    my $ofile = "korob_small_q.txt";
    open( OUT, ">$ofile" ) || die "cannot open $ofile: $!\n";
    print OUT "$audit_string";
    print OUT
"q\tT_0\tR_0\tToa\tlambda\tTs\tRs\tT_q\tR_q\terr_T\terr_R\tA1\tA1_r\tA1_t\n";
    print OUT @output;
    close OUT;
    print "wrote $ofile\n";
}

sub shock_from_D {

    my ( $symmetry, $gamma, $rho0, $Eoad, $D ) = @_;

    # D is shock speed
    # See page 94 of Zeldovich
    my $U   = 2 / ( $gamma + 1 ) * $D;
    my $P   = $rho0 * $U * $D;
    my $rho = ( $gamma + 1 ) / ( $gamma - 1 ) * $rho0;
    my $R   = 1.e99;
    my $T   = 1.e99;
    if ( $D > 0 ) {
        my $delta = 2 / ( $symmetry + 3 );
        my $pow = ( $symmetry + 1 );
        $R = ( $delta**2 * $Eoad / $D**2 )**( 1 / $pow );
        $T = $R**( 1 / $delta ) / sqrt($Eoad);
    }
    return [ $R, $T, $P, $U, $rho, $D ];
}
