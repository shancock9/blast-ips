# NAME

Blast::IPS

## SYNOPSIS

    # Installation
    perl Makefile.PL
    make
    make test
    make install

    # Run the example driver
    perl bin/blast_ips.pl


## DESCRIPTION

Blast::IPS evaluates the shock overpressure for the classical problem of a
point source explosion in an ideal homogeneous atmosphere.  This problem does
not have an analytic solution except for the special case that the ambient
atmospheric pressure is zero, so results are obtained by interpolating builtin
tables of pre-computed values.

This problem was one of the first problems to be addressed when digital
computers became available, and it remains an important theoretical model for
understanding explosions.  My own interest in this problem came from working
with models of the damaging effects of blast waves at long distances from
accidental rocket explosions.  Although many complex factors are involved in a
real explosion, a precise solution of this ideal problem would have provided a
very useful point of reference.  The numerical solutions and analytical fits
available were of unknown accuracy and typically defined over a limited range.

The builtin tables cover the three one-dimensional symmetries (plane,
cylindrical, spherical) and values of the ideal gas gamma from 1.1 to 7 (6.5
for spherical symmetry).  They were prepared with calculations using a finite
difference method for the initial blast wave formation and the method of
characteristics for propagation to very long distances.

The emphasis of this project is on getting an accurate estimate of the maximum
error in the final interpolated values.  Accuracy information is largely
unavailable in existing published results, and when it is given it is usually
highly approximate.  The tables used here were created with numerical methods
which produced very smooth results and which converged as 1/N^2, where N is the
number of spatial points in the finite difference calculation. This made it
relatively easy to obtain an accurate error estimate, and the maximum error was
driven to very small values by increasing N.

The estimated maximum relative error in shock overpressure obtained after
interpolating the builtin tables with cubic interpolation for the gamma value
of one of the builtin tables is about 1.e-6 at all ranges. 

The software can also perform accurate interpolation to values of gamma between
the tabulated vaules, not just the gamma values of the builtin tables.  This
additional interpolation increases the relative error in overpressure ratio by
about 1.e-6 over most of the range (above about gamma=1.2), so the maximum
estimated relative error in overpressure ratio for arbitrary gamma and range is
about 2.e-6.

Besides shock strength values, the tables include impulses, residual energy
distribution, and information on the tail shock. 

Test problems and example scripts illustrate how to use the module.  A detailed 
description of the numerical methods used to construct the tables is being prepared.

# AUTHOR

Steven Hancock  `<perltidy@users.sourceforge.net>`

# COPYRIGHT

Copyright (c) Steven Hancock `<perltidy@users.sourceforge.new>`.

### LICENCE

The MIT License

## DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.
