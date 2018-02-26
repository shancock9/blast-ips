use ExtUtils::MakeMaker;
WriteMakefile(
NAME  => "Blast::IPS",
VERSION_FROM => "lib/Blast/IPS.pm",
(
  $] >= 5.005
? (
  ABSTRACT => 'Estimate Blast Strength for a Point Source Explosion in Ideal Gas Atmosphere',
  LICENSE => 'MIT',
  AUTHOR => 'Steve Hancock <perltidy@perltidy.sourceforge.net>'
  )
: ()
),
#EXE_FILES => ['bin/blast_ideal_point_source.pl'],
dist  => { COMPRESS => 'gzip', SUFFIX => 'gz' },
);
