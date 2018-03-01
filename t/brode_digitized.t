use strict;
use Test;

my $rbrode_ideal_point_source_digitized;

BEGIN {

# This test compares the graphical results for this problem published in 1954 by
# Hal Brode. He used a finite difference method in which the shock location is
# obtained with an artificial viscosity.  This has limited accuracy but can be
# easily applied to a wide range of problems.

=pod
References:

RAND PROJECT AIR FORCE SANTA MONICA CA, and Brode, H L. 1954. Numerical Solutions of Spherical Blast Waves. http://oai.dtic.mil/oai/oai?&verb=getRecord&metadataPrefix=html&identifier=ADA595878. 

Brode, Harold L. 1955. "Numerical Solutions of Spherical Blast Waves". Journal of Applied Physics. 26 (6): 766-775. 
=cut

    # Unfortunately Brode did not give digital results for this problem, but
    # rather only graphs, some of which are quite small in final published
    # form.  This table is digitized from fig 1 of Brode's 1954 Rand report.
    # The digitization error is estimated to be about 0.5 percent. The same
    # figure in his JAP paper was smaller and would have had a higher
    # digitization error.

    # [lambda, ovp]
    $rbrode_ideal_point_source_digitized = [
        [ 0.04610627114753, 1576.4129722386 ],
        [ 0.04739145036237, 1457.2028098546 ],
        [ 0.04852214597092, 1357.65680257604 ],
        [ 0.04948573213737, 1269.87215398409 ],
        [ 0.05066522783988, 1187.79091722597 ],
        [ 0.051875225788,   1102.30048881121 ],
        [ 0.05311045230247, 1039.15420655858 ],
        [ 0.05480750094261, 941.88959652284 ],
        [ 0.05633522069257, 870.662821733091 ],
        [ 0.05790552453401, 804.822297588503 ],
        [ 0.05975579105232, 729.491103787862 ],
        [ 0.06190560894995, 658.612546926041 ],
        [ 0.06413424699197, 592.283998992266 ],
        [ 0.06618201717891, 541.078154409523 ],
        [ 0.06856302967204, 486.58643774215 ],
        [ 0.07075383753047, 441.042051916673 ],
        [ 0.07330101696404, 398.189679909789 ],
        [ 0.075641476925,   360.910901212014 ],
        [ 0.07867196541032, 323.295753368095 ],
        [ 0.08214667699392, 282.845973076354 ],
        [ 0.08543975611366, 251.379582148897 ],
        [ 0.08851563718496, 225.174911123567 ],
        [ 0.09206191404476, 201.701905007839 ],
        [ 0.09612580725792, 177.860701138405 ],
        [ 0.09997927888128, 157.452659278142 ],
        [ 0.10358002225654, 141.039253291959 ],
        [ 0.10857779056275, 123.879658653037 ],
        [ 0.11116805313224, 115.417058033173 ],
        [ 0.1142646779684,  107.530087698037 ],
        [ 0.1165330195933,  100.97410211209 ],
        [ 0.120249976452,   93.3447721607543 ],
        [ 0.12458105375922, 83.2952047347421 ],
        [ 0.129068124699,   74.3275811938797 ],
        [ 0.1342406356684,  66.0647810497343 ],
        [ 0.14238885574824, 55.5789067689304 ],
        [ 0.15162820718095, 46.0288368077246 ],
        [ 0.16273726017825, 37.5249394378084 ],
        [ 0.17672988709338, 29.6455833103693 ],
        [ 0.18672221905518, 25.2359698574064 ],
        [ 0.19419542756635, 22.8754686381026 ],
        [ 0.21172541851307, 17.7905756177931 ],
        [ 0.22635078239569, 14.5609597038129 ],
        [ 0.24197527940297, 12.1537416391932 ],
        [ 0.26381227554797, 9.52671120257758 ],
        [ 0.29678776958093, 6.9576407390974 ],
        [ 0.32102954260055, 5.64992912646128 ],
        [ 0.35273861065185, 4.44604649431352 ],
        [ 0.38910068915112, 3.47129082164912 ],
        [ 0.41429910967078, 3.01347782490887 ],
        [ 0.4359764191244,  2.66795807082515 ],
        [ 0.49040286643973, 2.05858118927773 ],
        [ 0.55162380546936, 1.58839312993538 ],
        [ 0.61081824379879, 1.27971423239637 ],
        [ 0.66322517041752, 1.08079169791419 ],
        [ 0.68166657622078, 1.02697152023933 ],
        [ 0.71732009941273, 0.93817768751816 ],
        [ 0.7608218910157,  0.83706233389729 ],
        [ 0.80063929267519, 0.76167201633284 ],
        [ 0.84255197315838, 0.69035424550997 ],
        [ 0.88317283016218, 0.63313744599841 ],
        [ 0.94409152708854, 0.56046952084252 ],
        [ 1.00133219753948, 0.50203269315203 ],
        [ 1.09160621591325, 0.42895822543316 ],
        [ 1.18066448599778, 0.37674543612732 ],
        [ 1.25223200659475, 0.33879266569664 ],
        [ 1.33855431991929, 0.30346457436144 ],
        [ 1.4764202731956,  0.25827180635078 ],
        [ 1.60311249827632, 0.2268339210192 ],
        [ 1.78200123218987, 0.19533991751684 ],
        [ 1.98087467646391, 0.16755870685239 ],
        [ 2.14236764102001, 0.14949511444848 ],
        [ 2.34445013500602, 0.13078140538291 ],
        [ 2.49617483173356, 0.1190005585953 ],
        [ 2.72089009710977, 0.10533919617287 ],
        [ 3.02450982754867, 0.09071309819625 ],
        [ 3.30971044095179, 0.07998526720834 ],
        [ 3.57944752494144, 0.07192502143661 ],
        [ 3.85604417007523, 0.06493218240719 ],
        [ 4.07343696364219, 0.06025595860744 ],

    ];

    plan tests => 1;
}

# Brode did not give an estimate of his computational error, but it is in the
# range of 1 to 2 percent over most of range.  The error rises toward
# the end but remains below 5%, so this is what is used here.
my $TOL     = 0.05;
my $VERBOSE = 0;
if ($VERBOSE) {
    print "Comparison with digitized Brode Table with tolerance $TOL\n";
}

# Create a table for this case
use Blast::IPS;
my $symmetry    = 2;
my $gamma       = 1.4;
my %args        = ( 'symmetry' => $symmetry, 'gamma' => $gamma );
my $blast_table = Blast::IPS->new( \%args );
if ( !defined($blast_table) ) {
    die "missing table for sym=$symmetry, gamma=$gamma\n";
}

# Loop to calculate the maximum absolute relative error in overpressure
my $iQ = 'X';
my $err_max;
my $lambda_max = 0;
foreach my $point ( @{$rbrode_ideal_point_source_digitized} ) {
    my ( $lambda_t, $ovprat_t ) = @{$point};

    #last if ($lambda_t > $lambda_max);

    # Lookup the overpressure at the given scaled range
    my $Q = log($lambda_t);
    my $ret = $blast_table->lookup( $Q, $iQ );
    my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$ret};
    my $ovprat = exp($Y);
    my $lambda = exp($X);
    my $err    = abs( $ovprat - $ovprat_t ) / $ovprat_t;

    # print STDERR "lambda=$lambda, err=$err\n";
    if ( !defined($err_max) || $err > $err_max ) { $err_max = $err }
    if ( $lambda_t > $lambda_max ) { $lambda_max = $lambda_t }
}

my $err_max_pr = sprintf "%0.3g", $err_max;
if ($VERBOSE) {
    print "Brode digitized table error to lambda = $lambda_max: $err_max_pr\n";
}
ok( $err_max <= $TOL );

