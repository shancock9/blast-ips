use strict;
use warnings;
use Test;
use Blast::IPS;

my $rokhotsimskii;

BEGIN {

# This test compares the table of values published in two papers by D.E. Okhotsimskii et al.

=pod

References:

D. E. Okhotsimskii, I. L. Kondrasheva, Z. I. Vlasova, R. K. Kazakova, "Computation Of Point Explosion Taking Into Account Counter Pressure", 1957 (In Russian)

D.Ye. Okhotsimskii; Z.P. Vlasova, "The Behaviour Of Shock Waves At Large Distances From The Point Of Explosion", 1963

=cut

# I have converted from the scaling used in the original papers to the more modern
# scaling that is used here.  This table combines three different calculations:
# The first part comes from the first reference above.
# The second part is the finite difference calculation in the second paper.
# The third part is the extrapolation to long range with a weak shock theory in the
# second paper.

    # NOTE: The second paper is in English, but the first paper is in Russian
    # and I did my best to understand it using google translate.  If anyone has
    # a clean translation in English I would appreciate seeing it.

    # [lambda, ovp]

    $rokhotsimskii = [
        [ 0.04479047645261, 1742.3237 ],
        [ 0.04736123243623, 1473.9341 ],
        [ 0.04994157304802, 1257.1257 ],
        [ 0.05423699122587, 981.35847 ],
        [ 0.05853740590426, 780.81703 ],
        [ 0.06283330632091, 631.52942 ],
        [ 0.0671352057039,  518.20386 ],
        [ 0.0714299495916,  430.47993 ],
        [ 0.07573042547316, 361.5077 ],
        [ 0.08002639447958, 306.59432 ],
        [ 0.08432572016355, 262.27742 ],
        [ 0.09033894894576, 213.51226 ],
        [ 0.09894186343139, 162.88294 ],
        [ 0.10753715390087, 127.11566 ],
        [ 0.11613870081211, 101.17722 ],
        [ 0.12473444397404, 81.910144 ],
        [ 0.13333745765101, 67.31941 ],
        [ 0.14193301087205, 56.025197 ],
        [ 0.15053516981504, 47.150192 ],
        [ 0.159130680827,   40.085577 ],
        [ 0.16773197448372, 34.385817 ],
        [ 0.17976414926882, 28.127146 ],
        [ 0.19697794309458, 21.624607 ],
        [ 0.2141776179811,  17.038194 ],
        [ 0.23139457748829, 13.712077 ],
        [ 0.2485985471493,  11.237532 ],
        [ 0.265819685356,   9.357089 ],
        [ 0.28302845630053, 7.9011431 ],
        [ 0.3002474946052,  6.7548102 ],
        [ 0.31745819661541, 5.8376267 ],
        [ 0.33468103373781, 5.0948185 ],
        [ 0.35877680797247, 4.2719535 ],
        [ 0.39324737502568, 3.4080611 ],
        [ 0.42778663736617, 2.787708 ],
        [ 0.46226397897767, 2.3286093 ],
        [ 0.49676956901989, 1.9795533 ],
        [ 0.56576353834958, 1.4925263 ],
        [ 0.60023257532342, 1.3185739 ],
        [ 0.63476234061959, 1.176004 ],
        [ 0.71761985255297, 0.9225047 ],
        [ 0.78653939671195, 0.7744325 ],
        [ 0.85562401005358, 0.6629836 ],
        [ 0.92452058191757, 0.5766384 ],
        [ 0.9935980197146,  0.5081888 ],
        [ 1.06247022638377, 0.4527996 ],
        [ 1.20040371532233, 0.3692012 ],
        [ 1.338327179603,   0.3095564 ],
        [ 1.43473959075177, 0.2772641 ],
        [ 1.57251943861769, 0.2404266 ],
        [ 1.71029116123459, 0.2115796 ],
        [ 1.84806531087392, 0.1884397 ],
        [ 1.98584621396699, 0.1695159 ],
        [ 2.2614313406731,  0.140528 ],
        [ 2.39922627828723, 0.1292188 ],
        ## END OF CLOSIN TABLE ##

        ## BEGINNING OF LONG RANGE TABLE ##
        [ 2.72997319307014, 0.1079008 ],
        [ 3.00562722610896, 0.094599695 ],
        [ 3.27999124396027, 0.08407005 ],
        [ 3.55331776248041, 0.075547085 ],
        [ 3.82579355687444, 0.068512888 ],
        [ 4.09755949996652, 0.062615611 ],
        [ 4.36872501881182, 0.057604889 ],
        [ 4.90958243179129, 0.049558569 ],
        [ 5.44887557604768, 0.043393995 ],
        [ 5.98694940533542, 0.038529658 ],
        [ 7.06034883032864, 0.031362485 ],
        [ 8.13106660657255, 0.026352864 ],
        [ 10.2672048132585, 0.019846321 ],
        [ 12.3987417019689, 0.015828619 ],
        [ 14.5273026390786, 0.013114765 ],
        [ 16.6537889996169, 0.011165597 ],
        [ 18.7787505569273, 0.0097015943 ],
        [ 23.0254267376162, 0.0076558389 ],
        [ 27.2690902447997, 0.0062996328 ],
        [ 31.5106855067746, 0.0053379625 ],
        [ 35.750779180519,  0.0046221739 ],
        [ 39.9897353193986, 0.0040698376 ],
        [ 48.4651509297275, 0.0032750895 ],
        [ 56.9381870028393, 0.0027323216 ],
        [ 65.4095473752413, 0.0023391945 ],
        [ 73.8796689109725, 0.0020418956 ],
        [ 82.3488375765897, 0.0018095306 ],
        [ 92.934229463511,  0.0015774986 ],

        # NOTE: The point at overpressure 0.001255858 causes a significant
        # distrubance, and the table is smoother without it.  This point is
        # the point at time 94.038755 in the published table.
        #[ 118.335799807659, 0.001255858 ],
        [ 135.268322837312, 0.0010386793 ],
        [ 186.060446954192, 0.0007318879 ],
        [ 253.775924176815, 0.0005211939 ],
        [ 389.194965107568, 0.0003274757 ],
        [ 626.162154951373, 0.0001959973 ],
        [ 1032.37374469223, 0.0001146161 ],
        [ 2250.96960768988, 4.99231E-05 ],
        [ 9291.58864500105, 1.11501E-05 ],
    ];

    plan tests => 1;
}

# This is the most accurate of the published tables that I have found.  I did
# not see an estimate of the error in the published papers.  The error is below
# about 1 percent over the entire table.
my $TOL     = 0.011;
my $VERBOSE = 0;
if ($VERBOSE) {
    print
"Comparison with combined Okhotsimskii et al tables with tolerance $TOL\n";
}

# Create a table for this case
my $symmetry    = 2;
my $gamma       = 1.4;
my $blast_table =
          Blast::IPS->new( 'symmetry' => $symmetry, 'gamma' => $gamma );
if ( !defined($blast_table) ) {
    die "missing table for sym=$symmetry, gamma=$gamma\n";
}

# Loop to calculate the maximum absolute relative error in overpressure
my $err_max;
my $lambda_max = 0;
foreach my $point ( @{$rokhotsimskii} ) {
    my ( $lambda_t, $ovprat_t ) = @{$point};

    #last if ($lambda_t > $lambda_max);

    # Lookup the overpressure at the given scaled range
    my $Q = log($lambda_t);
    my $ret = $blast_table->wavefront( 'X' => $Q )->{'TableVars'};
    my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};
    my $ovprat = exp($Y);
    my $lambda = exp($X);
    my $err    = abs( $ovprat - $ovprat_t ) / $ovprat_t;

    #print STDERR "lambda=$lambda, err=$err\n";
    if ( !defined($err_max) || $err > $err_max ) { $err_max = $err }
    if ( $lambda_t > $lambda_max ) { $lambda_max = $lambda_t }
}

my $err_max_pr    = sprintf "%0.3g", $err_max;
my $lambda_max_pr = sprintf "%0.5g", $lambda_max;
if ($VERBOSE) {
    print
"Okhotsimskii combined tables error to lambda = $lambda_max_pr: $err_max_pr\n";
}
ok( $err_max <= $TOL );

