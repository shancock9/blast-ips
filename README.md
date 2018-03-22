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
interpolating a table of values.

A table of values may be supplied as a call argument, or alternatively one of
the builtin tables may be used.  

The builtin tables cover the three one-dimensional symmetries (plane,
cylindrical, spherical) and several values of the ideal gas gamma (1.1, 1.2,
1.3, 1.4, 1.667, 2 and 3).  They were prepared with calculations using the
finite difference method and the method of characteristics.  The estimated
maximum relative error in shock overpressure obtained by interpolating
the builtin tables is about 1.e-6 at all ranges. This is orders of magnitude
lower than any published tables that I have found.

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
