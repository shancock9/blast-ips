# NAME

Blast::IPS

## SYNOPSIS

    # Installation
    perl Makefile.PL
    make
    make test
    make install

    # Run the example driver
    perl blast_ideal_point_source.pl


## DESCRIPTION

Blast::IPS evaluates the shock overpressure for a point source explosion in an
ideal homogeneous atmosphere.  This problem has a precise mathematical
definition but it does not have an analytic solution except for the special
case that the ambient atmospheric pressure is zero, so results are obtained by
interpolating a builtin table of values.

The builtin tables cover the three one-dimensional symmetries (plane,
cylindrical, spherical) and several values of the ideal gas gamma ( from 1.1 to
about 7).  They were prepared with calculations using the finite difference
method and the method of characteristics.  

The emphasis of this project is on getting an accurate estimate of the maximum
error in the final interpolated values.  This is needed in in statistical
models and is largely unavailable in existing published results. The individual
point values of shock overpressure in the builtin tables have an estimated
maximum error of about 1.e-7 in most cases.

The estimated maximum relative error in shock overpressure obtained after
interpolating the builtin tables with cubic interpolation for one of the gamma
values is about 1.e-6 at all ranges. This is orders of magnitude lower than any
published tables that I have found.

Interpolation to arbitrary gamma values will soon be added and should have 
a maximum interpolation error for any gamma and any distance of about 1.e-4.

I will be adding documentation in the future. 

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
