package Blast::IPS;

# MIT License
# Copyright (c) 2018 Steven Hancock
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

use strict;
use warnings;
use 5.006;
our $VERSION = 20180225;

use Carp;
my $rtables;
my $rtables_info;

# Evaluate the shock overpressure for point source explosion in an ideal
# homogeneous atmosphere.  This does not have an analytic solution, so
# we interpolate tables of values.

# A table of values may be supplied as a call argument, or alternatively one of
# the builtin tables may be used.

# The builtin tables cover the three one-dimensional symmetries (plane,
# cylindrical, spherical) and several values of the ideal gas gamma (1.1, 1.2,
# 1.3, 1.4, 1.667, 2 and 3).  They were prepared with calculations using the
# finite difference method and the method of characteristics.  The estimated
# relative accuracy of interpolated shock overpressures depends on the table
# but is below 1.e-5 in all cases.

# TODO:
# Needs documentation and test cases
# Work is needed to handle interpolations beyond the ranges of the tables for some variables.
#    Right now the program writes warnings to STDERR when this happens.
# Question: how accurate would it be to interpolate on gamma (or gamma-1) between tables?
# Are there any theoretical results to guide this (such as low gamma theories)?

sub new {
    my ( $class, $rtable, $options ) = @_;
    my $self = {};
    bless $self, $class;
    $self->_setup( $rtable, $options );
    return $self;
}

sub _setup {
    my ( $self, $arg ) = @_;
    my ( $rtable, $options, $table_name );

    my $gamma;
    my $symmetry;
    my $alpha;

    # Convert missing arg to default spherical table
    if ( !defined($arg) ) { $arg = { symmetry => 2 } }

    # Handle table name abbreviations
    elsif ( !ref($arg) ) {
        if    ( $arg eq 'S' ) { $arg = { symmetry => 2 } }
        elsif ( $arg eq 'C' ) { $arg = { symmetry => 1 } }
        elsif ( $arg eq 'P' ) { $arg = { symmetry => 0 } }
    }

    my $reftype = ref($arg);

    # A reference to an array gives a specific table
    if ( $reftype eq "ARRAY" ) {
        $rtable = $arg;
    }

    # A reference to a hash gives a hash of call args
    # This should eventually be expanded to include things like Blast::Medium
    elsif ( $reftype eq "HASH" ) {
        $options    = $arg;
        $rtable     = $options->{rtable};
        $table_name = $options->{table_name};
        $gamma      = $options->{gamma};
        $symmetry   = $options->{symmetry};

        # Allow either 'symmetry' or the older 'ASYM' for the symmetry
        if ( !defined($symmetry) ) { $symmetry = $options->{ASYM}; }

        # Option to lookup table from symmetry and/or gamma
        if (   !defined($rtable)
            && !defined($table_name)
            && ( defined($symmetry) || defined($gamma) ) )
        {

            if ( !defined($gamma) )    { $gamma    = 1.4 }
            if ( !defined($symmetry) ) { $symmetry = 1.4 }

            my @keys = reverse sort keys %{$rtables_info};

            # Look through our list of tables for the best match
            my $err_est;
            foreach my $key (@keys) {
                my $item = $rtables_info->{$key};
                my ( $symmetry_t, $gamma_t, $err_est_t ) = @{$item};

                next if $symmetry_t != $symmetry;
                next if $gamma_t != $gamma;

                # use table with min error estimate
                if ( !defined($table_name) || $err_est_t < $err_est ) {
                    $table_name = $key;
                    $err_est    = $err_est_t;
                }
            }
            if ( defined($table_name) ) {
                ( $rtable, $table_name ) = get_builtin_table($table_name);

                #print "Found Table $table_name\n";
            }
        }
    }

    # A scalar value is the name of a pre-defined table
    elsif ( !$reftype ) {
        ( $rtable, $table_name ) = get_builtin_table($arg);
    }
    else {
        croak "Unexpected argument type $reftype in Blast::Table";
    }

    if ( defined($table_name) ) {
        my $item = $rtables_info->{$table_name};
        if ( defined($item) ) {
            ( $symmetry, $gamma, my $err_est ) = @{$item};
        }
    }
    else {
        $table_name = "NONAME";
    }

    my $error = "";
    if ( !defined($rtable) ) {
        $error .=
"Undefined Table in Blast::Table\nfor symmetry=$symmetry, gamma=$gamma, name=$table_name\n";
    }

    $gamma    = 1.4 unless defined($gamma);
    $symmetry = 2   unless defined($symmetry);
    if ( $gamma <= 1.0 ) { $error .= "Bad gamma='$gamma'\n" }
    if ( $symmetry != 0 && $symmetry != 1 && $symmetry != 2 ) {
        $error .= "Bad symmetry='$symmetry'\n";
    }

    my $table_error = _check_table($rtable);
    $error .= $table_error;

    my $num_table;
    if ( !$error ) {
        $rtable    = _setup_toa_table($rtable);
        $num_table = @{$rtable};
    }

    $self->{_rtable}     = $rtable;
    $self->{_table_name} = $table_name;
    $self->{_gamma}      = $gamma;
    $self->{_symmetry}   = $symmetry;
    $self->{_jl}         = -1;
    $self->{_ju}         = $num_table;    #@{$rtable};
    $self->{_error}      = $error;

    if ($error) { carp "$error\n" }
    else {
        $self->_end_model_setup();
    }
    return;
}

sub get_builtin_table {
    my ($table_name) = @_;
    return ( $rtables->{$table_name}, $table_name );
}

sub _check_table {

    # Do some simple checks on the table
    my ($rtable) = @_;
    my $error = "";

    if ( !defined($rtable) ) {
        return "No Table Defined";
    }

    my $ic_X    = 0;
    my $ic_Y    = 1;
    my $ic_dYdX = 2;
    my $ic_Z    = 3;
    my $ic_dZdX = 4;

    my $nrows = @{$rtable};
    if ( $nrows < 2 ) {
        $error .= "Table must have at least 2 rows\n";
        return $error;
    }
    my $ncols = @{ $rtable->[0] };
    if ( $ncols < 5 ) {
        $error .= "Table must have at least 5 cols\n";
        return $error;
    }

    my $is_monotonic = sub {

        # returns 1 if table is monotonic increasing
        # returns -1 if table is monotonic decreasing
        # returns 0 if table is non-monotonic

        my ($icol) = @_;
        my $mono = 0;
        my $i_non_mono;    # could be returned if needed
        my $num = @{$rtable};
        for ( my $i = 1 ; $i < $num ; $i++ ) {
            my $dF = $rtable->[$i]->[$icol] - $rtable->[ $i - 1 ]->[$icol];
            if ( $i == 1 ) { $mono = $dF > 0 ? 1 : $dF < 0 ? -1 : 0 }
            elsif ( $dF * $mono < 0 ) {
                $i_non_mono = $i;
                $mono       = 0;
                last;
            }
        }
        return ($mono);
    };

    # X must be increasing
    if ( !$is_monotonic->($ic_X) > 0 ) {
        $error .= "X is not monotonic increasing\n";
    }

    # Y must be decreasing
    if ( !$is_monotonic->($ic_Y) < 0 ) {
        $error .= "Y is not monotonic increasing\n";
    }

    # dYdX must be negative
    foreach my $row ( @{$rtable} ) {
        if ( $row->[$ic_dYdX] >= 0 ) {
            $error .= "dYdX not everywhere negative\n";
            last;
        }
    }

    # Z must be increasing
    if ( !$is_monotonic->($ic_Z) > 0 ) {
        $error .= "Z is not monotonic increasing\n";
    }

    # dZdX must be positive
    foreach my $row ( @{$rtable} ) {
        if ( $row->[$ic_dZdX] <= 0 ) {
            $error .= "dZdX not everywhere positive\n";
            last;
        }
    }

    return $error;
}

sub _setup_toa_table {
    my ($rtable) = @_;

    # Make a table of X, ln(T), dln(T)/dX, where T=time of arrival
    # for solving inverse problems involving T
    foreach my $row ( @{$rtable} ) {
        my ( $X, $Y, $dYdX, $Z, $dZdX ) = @{$row};
        my $z    = exp($Z);
        my $dzdX = $z * $dZdX;
        my $R    = exp($X);
        my $W    = log( $R - $z );
        my $dWdX = ( $R - $dzdX ) / ( $R - $z );
        $row->[5] = $W;
        $row->[6] = $dWdX;
    }
    return $rtable;
}

sub get_ASYM       { my $self = shift; return $self->{_symmetry}; }
sub get_symmetry   { my $self = shift; return $self->{_symmetry} }
sub get_gamma      { my $self = shift; return $self->{_gamma} }
sub get_alpha      { my $self = shift; return $self->{_alpha} }
sub get_error      { my $self = shift; return $self->{_error} }
sub get_table      { my $self = shift; return $self->{_rtable} }
sub get_table_name { my $self = shift; return $self->{_table_name} }

sub set_table {

    # install a new table, or
    # reset to defaults if called without arguments
    my ( $self, $rtable, $options ) = @_;
    $self->_setup( $rtable, $options );
    return;
}

sub clone {

    my ($self) = @_;
    my $class = ref($self);
    my $newobj = bless { %{$self} }, $class;
    my $rtable = $newobj->{_rtable};
    my $rcopy;
    foreach my $item ( @{$rtable} ) {
        push @{$rcopy}, [ @{$item} ];
    }
    $newobj->{_rtable} = $rcopy;
    return $newobj;
}

sub lookup {

    # This is the main entry for external calls.

    # Given:
    #     $Q = a value to lookup

    #     $id_Q defines what Q contains:
    #       0 or 'X' => ln(R)    [this is the default if $id_Q not specified]
    #       1 or 'Y' => ln(overpressure)
    #       3 or 'Z' => R-cT, where T=time of arrival
    #       5 or 'W' => ln(T)

    #      $interp = optional interpolation flag
    #	     = 0 => cubic [This is the default and recommended method]
    #	     = 1 => linear [This is only intended for error checking]

    # Returns
    # 	[X, Y, dY/dX, Z]

    # Positive phase duration and length can be computed from Z.
    # Note that W is not returned but Toa can be computed from
    #   Toa = exp(X)-Z

    my ( $self, $Q, $id_Q, $interp ) = @_;

    # do not proceed unless the table is good
    my $error = $self->{_error};
    if ($error) {
        carp "$error";
        return;
    }

    my $rtab = $self->{_rtable};
    my $ntab = @{$rtab};

    # Set table column number to search
    my $icol;
    if    ( !defined($id_Q) )     { $icol = 0 }
    elsif ( $id_Q =~ /^[0135]$/ ) { $icol = $id_Q }
    else {
        $id_Q = uc($id_Q);
        if    ( $id_Q eq 'X' ) { $icol = 0 }
        elsif ( $id_Q eq 'Y' ) { $icol = 1 }
        elsif ( $id_Q eq 'Z' ) { $icol = 3 }
        elsif ( $id_Q eq 'W' ) { $icol = 5 }
        else                   { croak "unexpected id_Q=$id_Q\n"; }
    }

    # Locate this point in the table
    my ( $il, $iu ) = $self->_locate_2d( $Q, $icol );

    # Handle case before start of the table
    if ( $il < 0 ) {
        return $self->_short_range_calc( $Q, $icol );
    }

    # Handle case beyond end of the table
    elsif ( $iu >= $ntab ) {
        return $self->_long_range_calc( $Q, $icol );
    }

    # otherwise interpolate within table
    else {
        return _interpolate_rows( $Q, $icol, $rtab->[$il], $rtab->[$iu],
            $interp );
    }
}

sub get_table_bounds {

    my ( $self ) = @_;

    # return the first and last rows of the table
    return $self->table_gen(2);
}

sub table_gen {

    my ( $self, $num, $X_b, $X_e ) = @_;

    # Generate a table of values uniformly spaced in ln(r)
    # between any two values of ln(r).

    # This is convenient for generating plot data

    # num is the total number of points [ default = 200 ]
    # $X_b = first value of ln(r) [ default is first table value ]
    # $X_e = last value of ln(r)  [ default is last table value ]

    my $rtable = $self->{_rtable};
    my $Xtmin  = $rtable->[0]->[0];
    my $Xtmax  = $rtable->[-1]->[0];

    my $NUM_LIM = 100000;    # a limit to avoid problems with a bad call
    if ( !defined($num) )  { $num = 200; }
    if ( $num < 1 )        { $num = 1 }
    if ( $num > $NUM_LIM ) { $num = $NUM_LIM }

    if ( !defined($X_b) ) { $X_b = $Xtmin; }
    if ( !defined($X_e) ) { $X_e = $Xtmax; }

    my $rtable_gen;
    my $dX = ( $X_e - $X_b );
    if ( $num > 1 ) { $dX = ( $X_e - $X_b ) / ( $num - 1 ) }
    my $X = $X_b;
    for ( my $n = 1 ; $n <= $num ; $n++ ) {
        push @{$rtable_gen}, $self->lookup($X);
        $X += $dX;
    }
    return ($rtable_gen);
}

sub _short_range_calc {
    my ( $self, $Q, $icol ) = @_;
    my $rtab     = $self->{_rtable};
    my $symmetry = $self->{_symmetry};
    my $gamma    = $self->{_gamma};
    my $alpha    = $self->{_alpha};
    my $delta    = ( 3 + $symmetry ) / 2;
    my $p_amb    = 1;                       ## PATCH!! FIXME

    # $icol is the column number of the variable $Q which we are given
    # 0 = X
    # 1 = Y = ln(ovprat)
    # 3 = Z
    # 5 = W = ln(TOA)

    my $ovprat_from_lambda = sub {
        my ($lambda) = @_;
        my $am2;
        if ( $symmetry == 0 ) {
            $am2 = 9 * $alpha * $gamma * $lambda;
        }
        elsif ( $symmetry == 1 ) {
            $am2 = 16 * $alpha * $gamma * $lambda**2;
        }
        elsif ( $symmetry == 2 ) {
            $am2 = 25 * $alpha * $gamma * $lambda**3;
        }
        my $up_over_c = 4 / ( $gamma + 1 ) / sqrt($am2);
        my $ovp_atm =
          4 * $gamma / $am2 / ( $gamma + 1 ) * ( 1 + sqrt( 1 + $am2 ) );

        # For future reference:
	# This is the equation on page 182 of Korobeinikov's book "point blast
	# theory".  It is missing a factor of gamma.
        # It is a cleaner way to code (and typeset) if you include the
        # missing factor of gamma.
        #  my $ovp_atm_check = 4 / ( $gamma + 1 ) / ( -1 + sqrt( 1 + $am2 ) );
        #  my $diff = ( $ovp_atm - $ovp_atm_check );

        my $dlnp_dlnr = -( $symmetry + 1 ) *
          ( 1 - 2 * $gamma / ( ( $gamma + 1 ) * sqrt( 1 + $am2 ) * $ovp_atm ) );
        return ( $ovp_atm, $dlnp_dlnr );
    };

    my $lambda_from_ovprat = sub {
        my ($ovprat) = @_;
        my $prat = $ovprat + 1;
        my $q = 2 * $gamma / ( ( $gamma + 1 ) * $prat + ( $gamma - 1 ) );
        my $uovcsq = ( 2 / ( $gamma + 1 ) * ( 1 - $q ) )**2 / $q;
        my $C =
          ( 4 / ( ( $gamma + 1 ) * ( $symmetry + 3 ) ) )**2 /
          ( $gamma * $alpha );
        my $lambda = ( $C / ($uovcsq) )**( 1 / ( $symmetry + 1 ) );
        return $lambda;
    };

    my ( $X_0, $Y_0, $dY_dX_0, $Z_0, $dZ_dX_0 ) = @{ $rtab->[0] };

    my $z_0     = exp($Z_0);
    my $dz_dX_0 = $z_0 * $dZ_dX_0;
    my $dY_dX_i = -( 1 + $symmetry );
    my $X_i;
    my $Y_i;
    my $r_0   = exp($X_0);
    my $lnt_0 = log( $r_0 - $z_0 );

    if ( $icol == 0 ) {
        $X_i = $Q;
        my $lambda = exp($X_i);
        ( my $ovprat, $dY_dX_i ) = $ovprat_from_lambda->($lambda);
        $Y_i = log($ovprat);
    }
    elsif ( $icol == 1 ) {

        #$X_i = $X_0 + ( $Q - $Y_0 ) / $dY_dX_i;
        $Y_i = $Q;
        my $ovprat = exp($Y_i);
        my $lambda = $lambda_from_ovprat->($ovprat);
        ( my $ovprat2, $dY_dX_i ) = $ovprat_from_lambda->($lambda);
        $X_i = log($lambda);
    }
    elsif ( $icol == 3 ) {

        # FIXME: old coding using simple blast theory
        # Need to use the better theory
        print STDERR "FIXME: Table.pm using old blast theory\n";
        my $dX_newton = sub {
            my $lnT_i = $lnt_0 + $delta * ( $X_i - $X_0 );
            my $T     = exp($lnT_i);
            my $R     = exp($X_i);
            my $fofx  = $R - $T - $Q;
            my $dfdx  = $R - $T * $delta;
            return -$fofx / $dfdx;
        };

        # First guess: Close-in, Toa greatly exceeds R because of the very
        # high shock speed, so an accurate first guess is z=R
        $X_i = log($Q);
        for ( my $it = 0 ; $it < 10 ; $it++ ) {
            my $dX_i = $dX_newton->();
            $X_i = $X_i + $dX_i;
            last if ( abs($dX_i) < 1.e-6 );
        }
        $Y_i = $Y_0 + $dY_dX_i * ( $X_i - $X_0 );
    }
    elsif ( $icol == 5 ) {

        # FIXME: old coding using simple blast theory
        print STDERR "FIXME: Table.pm using old blast theory\n";
        $X_i = $X_0 + ( 1 / $delta ) * ( $Q - $lnt_0 );
        $Y_i = $Y_0 + $dY_dX_i * ( $X_i - $X_0 );
    }
    else {
        carp "Unexpected column id=$icol";    # shouldn't happen
        return;
    }

    #my $Y_i       = $Y_0 + $dY_dX_i * ( $X_i - $X_0 );
    my $d2Y_dX2_i = 0;

    # FIXME: This still uses simple theory
    my $lnT_i     = $lnt_0 + $delta * ( $X_i - $X_0 );
    my $T         = exp($lnT_i);
    my $R         = exp($X_i);
    my $z_i       = $R - $T;
    my $dz_dX_i   = $R - $T * $delta;
    my $d2z_dX2_i = $R - $T * ($delta)**2;
    my $Z_i       = $z_i;
    my $dZ_dX_i   = $dz_dX_i;
    $Z_i     = log($z_i);
    $dZ_dX_i = $dz_dX_i / $z_i;
    return [ $X_i, $Y_i, $dY_dX_i, $Z_i, $dZ_dX_i, $d2Y_dX2_i, $d2z_dX2_i ];
}

sub _end_model_setup {
    my ($self) = @_;

    # Set parameters for evaluation beyond the ends of the table
    # Set asymptotic wave parameters, given a table and gamma
    my $rtable   = $self->{_rtable};
    my $gamma    = $self->{_gamma};
    my $symmetry = $self->{_symmetry};

    # Near region: we find the value of alpha at the first table point
    # based on R^(symmetry+1)*u**2=constant
    my ( $X_near, $Y_near, $dYdX_near, $Z_near, $dZdX_near ) =
      @{ $rtable->[0] };
    my $lambda = exp($X_near);
    my $Prat   = exp($Y_near) + 1;
    my $q      = 2 * $gamma / ( ( $gamma + 1 ) * $Prat + ( $gamma - 1 ) );
    my $uovcsq = ( 2 / ( $gamma + 1 ) * ( 1 - $q ) )**2 / $q;
    my $C      = $lambda**( $symmetry + 1 ) * $uovcsq;
    my $alpha =
      ( 4 / ( ( $gamma + 1 ) * ( $symmetry + 3 ) ) )**2 / ( $gamma * $C );
    $self->{_alpha} = $alpha;

    # Near region
    my ( $X_far, $Y_far, $dYdX_far, $Z_far, $dZdX_far ) = @{ $rtable->[-1] };
    my ( $A_far, $B_far, $Z_zero ) = ( 0, 0, 0 );
    if ( $symmetry == 2 ) {
        my $mu = -( $dYdX_far + 1 );
        if ( $mu > 0 ) {
            $A_far = exp($X_far) * exp($Y_far) / ( $gamma * sqrt( 2 * $mu ) );
            $B_far = 0.5 / $mu - $X_far;
        }
        else { $self->{_error} .= "ending dYdX=$dYdX_far is bad\n" }
        $Z_zero =
          exp($Z_far) - 0.5 * ( $gamma + 1 ) * $A_far * sqrt( $X_far + $B_far );
	$Z_zero = log($Z_zero);
    }
    $self->{_A_far}  = $A_far;
    $self->{_B_far}  = $B_far;
    $self->{_Z_zero} = $Z_zero;
    return;
}

sub _long_range_calc {
    my ( $self, $Q, $icol ) = @_;
    my $symmetry = $self->{_symmetry};
    return $self->_long_range_sphere($Q,$icol) if ( $symmetry == 2 );
    return $self->_long_range_non_sphere($Q,$icol);
}

sub _long_range_sphere {
    my ( $self, $Q, $icol ) = @_;
    my $symmetry = $self->{_symmetry};
    my $A_far    = $self->{_A_far};
    my $B_far    = $self->{_B_far};
    my $Z_zero   = $self->{_Z_zero};
    my $gamma    = $self->{_gamma};
    my $log_ga   = 0;
    if ( $symmetry == 2 ) {
        $log_ga = log( $gamma * $A_far );
    }
    my $kA = 0.5 * ( $gamma + 1 ) * $A_far;
    my $rtab = $self->{_rtable};
    my ( $X_e, $Y_e, $dYdX_e, $Z_e ) = @{ $rtab->[-1] };

    my $X_i;
    if ( $icol == 0 ) {
        $X_i = $Q;
    }
    else {
        my $dX_newton;
        if ( $icol == 1 ) {

            # first guess
            $X_i = $X_e + $Y_e - $Q;

            $dX_newton = sub {
                my $term = $X_i + $B_far;
                return 0 if ( $term <= 0 );    # shouldn't happen
                my $fofx = $log_ga - $X_i - 0.5 * log($term) - $Q;
                my $mu   = 1 + 0.5 / $term;
                return $fofx / $mu;
            };
        }
        elsif ( $icol == 3 ) {

	    # FIXME: definition of Z has changed
            $dX_newton = sub {
                my $term = $X_i + $B_far;
                return 0 if ( $term <= 0 );    # shouldn't happen
                my $Z_i  = $Z_zero + $kA * sqrt($term);
                my $fofx = $Z_i - $Q;
                my $dfdx = 0.5 * $kA / sqrt($term);
                return -$fofx / $dfdx;
            };

            # Iterating for Z at long range is a little tricky..
            # Evaluate Z at the end of the table and at a large X (X=100),
            # and only iterate if there is a solution in that range.
            $X_i = $X_e;
            my $dX_e = $dX_newton->();
            $X_i = 100;
            my $dX_i = $dX_newton->();
            if ( $dX_i > 0 ) {
                $dX_newton = sub { 0 };
            }
            else {
                $X_i = abs($dX_i) < abs($dX_e) ? $X_i + $dX_i : $X_e + $dX_e;
            }
        }
        elsif ( $icol == 5 ) {

            # Since Z hardly changes with distance, an accurate first guess
            # is made using the last Z in the table
            ##$X_i       = log( exp($Q) - $Z_e );
            $X_i       = log( exp($Q) + exp($Z_e) );
            $dX_newton = sub {
                my $term = $X_i + $B_far;
                return 0 if ( $term <= 0 );    # shouldn't happen
                my $R    = exp($X_i);
                my $Toa  = $R - $Z_zero - $kA * sqrt($term);
                my $dTdX = $R - 0.5 * $kA / sqrt($term);
                my $fofx = log($Toa) - $Q;
                my $dfdx = $dTdX / $Toa;
                return -$fofx / $dfdx;
            };
        }
        else {
            carp "Unexpected column id=$icol";    # shouldn't happen
            return;
        }

        # Newton iterations
        for ( my $it = 0 ; $it < 10 ; $it++ ) {
            my $dX_i = $dX_newton->();
            $X_i = $X_i + $dX_i;
            last if ( abs($dX_i) < 1.e-6 );
        }
    }

    my $term      = $X_i + $B_far;
    my $Y_i       = $log_ga - $X_i - 0.5 * log($term);
    my $mum_i     = 1 / ( 2 * $term );
    my $dY_dX_i   = -( 1 + $mum_i );
    my $d2Y_dX2_i = 2 * $mum_i**2;

    #    my $d3Y_d3X = -8 * $mum**3;
    my $Z_i = log(exp($Z_zero) + $kA * sqrt($term));

    # FIXME: These are incorrect because of the change in definition of Z
    # FIXME: these give slight discontinuities; can it be fixed?
    # For example, should the Z table be splined?
    my $dZ_dX_i   = $kA * 0.5 / sqrt($term);
    my $d2Z_dX2_i = -0.25 * $kA * $term**-1.5;

    return [ $X_i, $Y_i, $dY_dX_i, $Z_i, $dZ_dX_i, $d2Y_dX2_i, $d2Z_dX2_i ];
}

sub _long_range_non_sphere {
    my ( $self, $Q, $icol ) = @_;
    my $symmetry = $self->{_symmetry};
    my $gamma    = $self->{_gamma};
    my $rtab     = $self->{_rtable};
    my ( $X_e, $Y_e, $dY_dX_e, $Z_e, $dZ_dX_e ) = @{ $rtab->[-1] };
    my ( $X_i, $Y_i, $dY_dX_i, $Z_i, $dZ_dX_i, $d2Y_dX2_i, $d2Z_dX2_i );
    $dY_dX_i   = $dY_dX_e;
    $dZ_dX_i   = $dZ_dX_e;
    $d2Y_dX2_i = 0;
    $d2Z_dX2_i = 0;

    if ( $icol == 0 ) {
        $X_i = $Q;
        $Y_i = $Y_e + $dY_dX_i * ( $X_i - $X_e );
        $Z_i = $Z_e + $dZ_dX_i * ( $X_i - $X_e );
    }
    elsif ( $icol == 1 ) {
        $Y_i = $Q;
        $X_i = $X_e + ( $Y_i - $Y_e ) / $dY_dX_i;
        $Z_i = $Z_e + $dZ_dX_i * ( $X_i - $X_e );
    }
    elsif ( $icol == 3 ) {
        $Z_i = $Q;
        $X_i = $X_e + ( $Z_i - $Z_e ) / $dZ_dX_i;
        $Y_i = $Y_e + $dY_dX_i * ( $X_i - $X_e );
    }
    elsif ( $icol == 5 ) {

        # Since Z hardly changes with distance, an accurate first guess
        # is made using the last Z in the table
        $X_i = log( exp($Q) + exp($Z_e) );
        $Y_i = $Y_e + $dY_dX_i * ( $X_i - $X_e );
        $Z_i = $Z_e + $dZ_dX_i * ( $X_i - $X_e );
    }

    return [ $X_i, $Y_i, $dY_dX_i, $Z_i, $dZ_dX_i, $d2Y_dX2_i, $d2Z_dX2_i ];
}

sub _locate_2d {
    my ( $self, $xx, $icol ) = @_;

    # Binary search for two consecutive table row indexes, jl and ju, of a
    # 2D matrix such that the value of x lies between these two table values.
    # $icol is the column of the variable x
    # If x is out of bounds, returns either jl<0 or ju>=N

    my $rx             = $self->{_rtable};
    my $num            = @{$rx};
    my $dx_is_positive = $rx->[-1]->[$icol] > $rx->[0]->[$icol];

    # Set search bounds to previous location ...
    my $jl = $self->{_jl};
    my $ju = $self->{_ju};

    # ... but reset to beyond end of table if no longer valid
    $jl = -1
      if ( !defined($jl)
        || $jl < 0
        || $jl >= $num
        || ( $xx > $rx->[$jl]->[$icol] ne $dx_is_positive ) );
    $ju = $num
      if ( !defined($ju)
        || $ju < 0
        || $ju >= $num
        || ( $xx > $rx->[$ju]->[$icol] eq $dx_is_positive ) );

    # Loop until the requested point lies in a single interval
    while ( $ju - $jl > 1 ) {
        my $jm = int( ( $jl + $ju ) / 2 );
        if ( $xx > $rx->[$jm]->[$icol] eq $dx_is_positive ) {
            $jl = $jm;
        }
        else {
            $ju = $jm;
        }
    }
    $self->{_jl} = $jl;
    $self->{_ju} = $ju;
    return ( $jl, $ju );
}

sub _interpolate_rows {
    my ( $Q, $icol, $row_b, $row_e, $interp ) = @_;

    # Interpolate all of the variables between two table rows. Each row has
    # these values
    # 	[X,Y,Z,dYdX,dZdX,W,dWdX]

    my ( $X_b, $Y_b, $dY_dX_b, $Z_b, $dZ_dX_b, $W_b, $dW_dX_b ) = @{$row_b};
    my ( $X_e, $Y_e, $dY_dX_e, $Z_e, $dZ_dX_e, $W_e, $dW_dX_e ) = @{$row_e};

    # If this is an inverse problem, first find X=X_i
    my $X_i;
    if ( $icol == 0 ) {
        $X_i = $Q;
    }
    elsif ( $icol == 1 ) {
        ( $X_i, my $dX_dY_i, my $d2X_dY2_i ) =
          _interpolate_scalar( $Q, $Y_b, $X_b, 1 / $dY_dX_b,
            $Y_e, $X_e, 1 / $dY_dX_e, $interp );
    }
    elsif ( $icol == 3 ) {
        ( $X_i, my $dX_dZ_i, my $d2X_dZ2_i ) =
          _interpolate_scalar( $Q, $Z_b, $X_b, 1 / $dZ_dX_b,
            $Z_e, $X_e, 1 / $dZ_dX_e, $interp );
    }
    elsif ( $icol == 5 ) {
        ( $X_i, my $dX_dW_i, my $d2X_dW2_i ) =
          _interpolate_scalar( $Q, $W_b, $X_b, 1 / $dW_dX_b,
            $W_e, $X_e, 1 / $dW_dX_e, $interp );
    }
    else {
        carp "unknown call type icol='$icol'";    # Shouldn't happen
    }

    # Step 2: now given X, do a complete normal interpolation
    my ( $Y_i, $dY_dX_i, $d2Y_dX2_i ) =
      _interpolate_scalar( $X_i, $X_b, $Y_b, $dY_dX_b, $X_e, $Y_e, $dY_dX_e,
        $interp );
    my ( $Z_i, $dZ_dX_i, $d2Z_dX2_i ) =
      _interpolate_scalar( $X_i, $X_b, $Z_b, $dZ_dX_b, $X_e, $Z_e, $dZ_dX_e,
        $interp );

    return [ $X_i, $Y_i, $dY_dX_i, $Z_i, $dZ_dX_i, $d2Y_dX2_i, $d2Z_dX2_i ];
}

sub _interpolate_scalar {
    my ( $xx, $x1, $y1, $dydx1, $x2, $y2, $dydx2, $interp ) = @_;
    my ( $yy, $dydx, $d2ydx2, $d3ydx3 );
    if ( defined($interp) && $interp == 1 ) {
        ( $yy, $dydx ) = _linear_interpolation( $xx, $x1, $y1, $x2, $y2 );
        ( $dydx, $d2ydx2 ) =
          _linear_interpolation( $xx, $x1, $dydx1, $x2, $dydx2 );
        $d3ydx3 = 0;
    }
    else {
        ( $yy, $dydx, $d2ydx2, $d3ydx3 ) =
          _cubic_interpolation( $xx, $x1, $y1, $dydx1, $x2, $y2, $dydx2 );
    }
    return ( $yy, $dydx, $d2ydx2, $d3ydx3 );
}

sub _linear_interpolation {
    my ( $xx, $x1, $y1, $x2, $y2 ) = @_;
    my $dx21 = $x2 - $x1;
    my $dy21 = $y2 - $y1;
    my $dydx = ( $dx21 == 0 ) ? 0 : $dy21 / $dx21;
    my $yy   = $y1 + $dydx * ( $xx - $x1 );
    return ( $yy, $dydx );
}

sub _cubic_interpolation {
    my ( $xx, $x1, $y1, $dydx1, $x2, $y2, $dydx2 ) = @_;

    # Given the value of a function and its slope at two points,
    # Use a cubic polynomial to find the value of the function and its
    # derivatives at an intermediate point

    # let z=(x-x1)/(x2-x1)
    #  y=a+b*z+c*z**2+d*z**3
    my $dx21   = $x2 - $x1;
    my $dy21   = $y2 - $y1;
    my $deriv1 = $dydx1 * $dx21;
    my $deriv2 = $dydx2 * $dx21;
    my $aa     = $y1;
    my $bb     = $deriv1;
    my $cc     = 3 * $dy21 - 2 * $deriv1 - $deriv2;
    my $dd     = -2 * $dy21 + $deriv1 + $deriv2;
    my $zz     = ( $xx - $x1 ) / $dx21;
    my $yy     = $aa + $zz * ( $bb + $zz * ( $cc + $zz * $dd ) );
    my $yp1    = $bb + $zz * ( 2 * $cc + $zz * 3 * $dd );
    my $yp2    = 2 * $cc + $zz * 6 * $dd;
    my $yp3    = 6 * $dd;
    my $dydx   = $yp1 / $dx21;
    my $d2ydx2 = $yp2 / $dx21**2;
    my $d3ydx3 = $yp3 / $dx21**3;
    return ( $yy, $dydx, $d2ydx2, $d3ydx3 );
}

##################################################################
# TABLE DATA
##################################################################

BEGIN {

    # name => [symmetry, gamma, error]
    # symmetry = 0,1,2 for Plane, Cylindrical, Spherical
    # Error = the sum of all run errors plus log interpolation error;

    # name -> [symmetry, gamma, error ]
    $rtables_info = {
        P4000_G1X1   => [ 0, 1.1,   3.9e-6 ],
        P4000_G1X2   => [ 0, 1.2,   1.5e-6 ],
        P4000_G1X3   => [ 0, 1.3,   2.7e-6 ],
        P4000_G1X4_A => [ 0, 1.4,   1.7e-6 ],
        P4000_G1X667 => [ 0, 1.667, 1.5e-6 ],
        P4000_G2     => [ 0, 2,     6.7e-6 ],
        P4000_G3     => [ 0, 3,     3.4e-6 ],
        P4000_G7     => [ 0, 7,     2.5e-6 ],
        C4000_G1X1   => [ 1, 1.1,   3.5e-6 ],
        C4000_G1X2   => [ 1, 1.2,   3.4e-6 ],
        C4000_G1X3   => [ 1, 1.3,   1.9e-6 ],
        C8000_G1X4   => [ 1, 1.4,   9.e-7 ],
        C4000_G1X667 => [ 1, 1.667, 2.2e-6 ],
        C4000_G2     => [ 1, 2,     2.1e-6 ],
        C4000_G3     => [ 1, 3,     1.8e-6 ],
        C4000_G7     => [ 1, 7,     1.6e-6 ],
        S4000_G1X1   => [ 2, 1.1,   6.6e-6 ],
        S4000_G1X2   => [ 2, 1.2,   1.8e-6 ],
        S4000_G1X3   => [ 2, 1.3,   1.6e-6 ],
        S8000_G1X4   => [ 2, 1.4,   1.4e-6 ],
        S16000_G1X4  => [ 2, 1.4,   9.e-7 ],
        S4000_G1X667 => [ 2, 1.667, 1.7e-6 ],
        S4000_G2     => [ 2, 2,     1.9e-6 ],
        S4000_G3     => [ 2, 3,     3.7e-6 ],
    };

    # All tables

    # Let.. 
    #  P0 = ambient atmospheric pressure
    #  E  = explosion energy
    #  n  = symmetry = 0,1, or 2  for plane, cylindrical or spherical
    #  d = (E/P0)^(1/(n+1)) = scaling distance
    #  lambda = scaled range = r/d, where r is distance
    #  tau = scaled time = c0 t / d, where t is time and c0 is ambient sound speed

    # The table format is:
    #  [ X, Y, dY/dX, Z, dZ/dX ]
    # where
    #  X = ln(lambda) where lambda = scaled range
    #  Y = ln(overpressure ratio) = ln (P/P0-1)
    #  dYdX = derivative of Y with respect to X 
    #  Z = ln ( lambda-tau) 
    #  dZdX = derivative of X with respect to X 

    # Comments above each table give their estimated errors.  The maximum error
    # for a table is the estimated maximum error after interpolation to any
    # range within the bounds of the table.  It is estimated as the sum of a
    # finite difference error, a method of characteristics error, and an
    # interpolation error.

    # The finite difference (FD) calculations converged like 1/N^2, where N is
    # the number of points between the origin and the shock front.  The digits
    # after the first letter in a table name give the number of points used for
    # that table.  Thus S4000 is for a table based on a calculation with 4000
    # points between the origin and the shock front.  The error was accurately
    # estimated by calculating with N and N/2 and comparing the results.  The
    # error was found to vary roughly as 10/N^2.

    # The method of characteristics (MOC) calculations which carried the wave
    # to long range were found to have errors of about the same order of
    # magnitude as the FD calculations.  

    # Cubic interpolation among the table points has a maximum error typically
    # between 5.e-7 and 1.e-6, depending on the number of table points used.  

    # These errors can be made almost arbitrarily small by increasing the number 
    # of points in the FD calculation and in the tables.  My goal for these tables
    # was to achieve a maximum error on the order of 1.e-6.

    $rtables = {

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 8.55e-07
        #    Est Max FD Error (based on N/2 run) to R=30.01: 4.5e-7
        #    Est Max MOC error for R>30.01: -2.26e-07
        #    Max interpolation error with cubic splines: 1.0041e-06
        #    Est Max overall error after cubic interpolation: 1.7e-6
        S4000_G1X667 => [
            [ -4.47137282, 12.00033596,  -2.99997697, -4.47248183, 0.99833559 ],
            [ -4.31618890, 11.53478855,  -2.99997501, -4.31758877, 0.99789875 ],
            [ -3.84862177, 10.13211965,  -2.99988249, -3.85144653, 0.99575704 ],
            [ -3.64686028, 9.52687163,   -2.99975988, -3.65068525, 0.99425200 ],
            [ -3.28588411, 8.44409941,   -2.99929151, -3.29246597, 0.99009693 ],
            [ -3.00126829, 7.59057033,   -2.99833772, -3.01137141, 0.98477672 ],
            [ -2.76703686, 6.88844032,   -2.99664924, -2.78141994, 0.97829310 ],
            [ -2.56729075, 6.29011738,   -2.99391620, -2.58673888, 0.97060045 ],
            [ -2.39128128, 5.76349639,   -2.98973031, -2.41666316, 0.96156859 ],
            [ -2.23473417, 5.29590047,   -2.98367910, -2.26691076, 0.95121086 ],
            [ -2.09253823, 4.87219400,   -2.97521960, -2.13246210, 0.93939688 ],
            [ -1.96131488, 4.48248428,   -2.96371443, -2.01004202, 0.92598947 ],
            [ -1.83837697, 4.11902067,   -2.94839153, -1.89710783, 0.91080622 ],
            [ -1.72103733, 3.77418031,   -2.92821827, -1.79121880, 0.89353744 ],
            [ -1.60697473, 3.44161262,   -2.90185362, -1.69039718, 0.87377311 ],
            [ -1.49254619, 3.11146170,   -2.86703808, -1.59170247, 0.85065355 ],
            [ -1.37173893, 2.76786011,   -2.81938086, -1.49060439, 0.82237457 ],
            [ -1.24566599, 2.41624996,   -2.75610393, -1.38901511, 0.78847470 ],
            [ -1.13141934, 2.10525248,   -2.68624277, -1.30087144, 0.75399553 ],
            [ -1.02527623, 1.82402814,   -2.61118909, -1.22266984, 0.71909589 ],
            [ -0.92450745, 1.56482554,   -2.53220152, -1.15196661, 0.68389168 ],
            [ -0.82828605, 1.32502313,   -2.45148814, -1.08783475, 0.64894823 ],
            [ -0.73311554, 1.09564476,   -2.36852715, -1.02774807, 0.61370600 ],
            [ -0.63806030, 0.87449401,   -2.28454646, -0.97109024, 0.57843465 ],
            [ -0.54266586, 0.66055872,   -2.20099370, -0.91758125, 0.54354116 ],
            [ -0.44588580, 0.45155918,   -2.11857811, -0.86665054, 0.50917516 ],
            [ -0.34682205, 0.24571357,   -2.03796069, -0.81789114, 0.47551580 ],
            [ -0.24470969, 0.04165019,   -1.95976573, -0.77102501, 0.44276799 ],
            [ -0.13872176, -0.16201584,  -1.88444759, -0.72579682, 0.41109753 ],
            [ -0.02810004, -0.36643045,  -1.81241357, -0.68202820, 0.38067324 ],
            [ 0.08822861,  -0.57320594,  -1.74380561, -0.63946627, 0.35156599 ],
            [ 0.21116652,  -0.78351626,  -1.67884786, -0.59797826, 0.32388827 ],
            [ 0.34171807,  -0.99861614,  -1.61767489, -0.55743698, 0.29771675 ],
            [ 0.48100292,  -1.21985213,  -1.56035046, -0.51772147, 0.27309643 ],
            [ 0.63027732,  -1.44868634,  -1.50688029, -0.47871572, 0.25004400 ],
            [ 0.78764585,  -1.68189434,  -1.45817287, -0.44106540, 0.22896358 ],
            [ 0.95111520,  -1.91661175,  -1.41462206, -0.40522155, 0.21003832 ],
            [ 1.12234888,  -2.15540280,  -1.37544624, -0.37075624, 0.19293871 ],
            [ 1.30256329,  -2.40001282,  -1.34012687, -0.33741701, 0.17744520 ],
            [ 1.49288528,  -2.65195504,  -1.30825048, -0.30501751, 0.16338159 ],
            [ 1.69430756,  -2.91249133,  -1.27948422, -0.27342836, 0.15060561 ],
            [ 1.90865120,  -3.18387385,  -1.25344147, -0.24242836, 0.13894969 ],
            [ 2.13699165,  -3.46732289,  -1.22989403, -0.21194602, 0.12831604 ],
            [ 2.38051289,  -3.76416618,  -1.20862333, -0.18191071, 0.11861142 ],
            [ 2.64116814,  -4.07662144,  -1.18937919, -0.15218203, 0.10972783 ],
            [ 2.92056030,  -4.40642476,  -1.17198018, -0.12269172, 0.10158853 ],
            [ 3.22078249,  -4.75584855,  -1.15624208, -0.09334381, 0.09411554 ],
            [ 3.54898208,  -5.13288761,  -1.14180590, -0.06362938, 0.08714572 ],
            [ 3.96597996,  -5.60575784,  -1.12669738, -0.02889018, 0.07970235 ],
            [ 4.42272079,  -6.11719875,  -1.11329970, 0.00592118,  0.07294231 ],
            [ 4.92427564,  -6.67249200,  -1.10141119, 0.04091387,  0.06678648 ],
            [ 5.47590072,  -7.27704381,  -1.09086079, 0.07615793,  0.06117104 ],
            [ 6.08463942,  -7.93813391,  -1.08147663, 0.11178182,  0.05603067 ],
            [ 6.75692320,  -8.66228579,  -1.07312861, 0.14781803,  0.05132126 ],
            [ 7.50077412,  -9.45766686,  -1.06568913, 0.18433559,  0.04699844 ],
            [ 8.32720671,  -10.33553710, -1.05903202, 0.22147939,  0.04301549 ],
            [ 9.24642976,  -11.30618386, -1.05306795, 0.25927992,  0.03934418 ],
            [ 10.26984828, -12.38107711, -1.04771815, 0.29775903,  0.03595971 ],
            [ 11.40966486, -13.57244846, -1.04291569, 0.33691310,  0.03284147 ],
            [ 12.68068072, -14.89516836, -1.03859766, 0.37677083,  0.02996802 ],
            [ 14.09789637, -16.36423499, -1.03471434, 0.41730789,  0.02732339 ],
            [ 15.67871199, -17.99706674, -1.03121987, 0.45851676,  0.02489129 ],
            [ 17.19473728, -19.55827846, -1.02847678, 0.49473603,  0.02294528 ],
            [ 19.07787838, -21.49231998, -1.02567953, 0.53598303,  0.02092543 ],
            [ 20.80513759, -23.26205125, -1.02356036, 0.57074773,  0.01937005 ],
            [ 22.76618451, -25.26726529, -1.02154460, 0.60722112,  0.01786936 ],
            [ 25.00413757, -27.55122541, -1.01963044, 0.64554805,  0.01642416 ],
            [ 27.57244981, -30.16754913, -1.01781611, 0.68589294,  0.01503525 ],
            [ 30.53801771, -33.18332079, -1.01609989, 0.72844356,  0.01370339 ],
            [ 33.98543568, -36.68335067, -1.01448011, 0.77341568,  0.01242933 ],
            [ 36.95119866, -39.69030173, -1.01332758, 0.80888257,  0.01151216 ],
            [ 40.30516087, -43.08706649, -1.01222769, 0.84597127,  0.01062822 ],
            [ 44.11696742, -46.94342558, -1.01117978, 0.88481723,  0.00977780 ],
            [ 48.47251372, -51.34542786, -1.01018320, 0.92557381,  0.00896118 ],
            [ 53.47870276, -56.40015113, -1.00923730, 0.96841568,  0.00817864 ],
            [ 59.26990725, -62.24216829, -1.00834145, 1.01354295,  0.00743044 ],
            [ 66.01687133, -69.04245435, -1.00749501, 1.06118635,  0.00671686 ],
            [ 71.15175539, -74.21441107, -1.00695786, 1.09447702,  0.00626048 ],
            [ 76.88820631, -79.98925844, -1.00644222, 1.12909440,  0.00581967 ],
            [ 83.32289182, -86.46376486, -1.00594790, 1.16513911,  0.00539449 ],
            [ 90.57262031, -93.75485141, -1.00547473, 1.20272380,  0.00498501 ],
            [ 98.77958315, -102.00483615, -1.00502252, 1.24197512, 0.00459129 ],
            [
                108.11825879, -111.38834066, -1.00459110, 1.28303613,
                0.00421341
            ],
            [
                118.80460842, -122.12148869, -1.00418028, 1.32606931,
                0.00385141
            ],
            [
                131.10847507, -134.47430886, -1.00378990, 1.37126023,
                0.00350537
            ],
            [
                145.37052683, -148.78768259, -1.00341978, 1.41882213,
                0.00317533
            ],
            [
                153.36729955, -156.81108007, -1.00324226, 1.44356780,
                0.00301633
            ],
            [
                162.02575067, -165.49684361, -1.00306973, 1.46900172,
                0.00286135
            ],
            [
                171.42040609, -174.91953638, -1.00290219, 1.49516117,
                0.00271041
            ],
            [
                181.63655735, -185.16449034, -1.00273960, 1.52208657,
                0.00256349
            ],
            [
                192.77218041, -196.32972536, -1.00258194, 1.54982181,
                0.00242062
            ],
            [
                204.94026522, -208.52827922, -1.00242919, 1.57841470,
                0.00228179
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 7.22e-07
        #    Est Max FD Error (based on N/2 run) to R=30.00: 6.7e-7
        #    Est Max MOC error for R>30.00: -1.71e-06
        #    Max interpolation error with cubic splines: 1.0066e-06
        #    Est Max overall error after cubic interpolation: 3.4e-6
        C4000_G1X2 => [
            [ -7.08883517, 12.00033665,  -1.99998660, -7.09013028, 0.99870406 ],
            [ -6.27216020, 10.36701524,  -1.99992757, -6.27509333, 0.99706266 ],
            [ -5.60663839, 9.03607121,   -1.99972964, -5.61235248, 0.99427031 ],
            [ -5.15788142, 8.13875394,   -1.99933494, -5.16684544, 0.99099861 ],
            [ -4.75032769, 7.32406408,   -1.99850278, -4.76382962, 0.98641649 ],
            [ -4.42070182, 6.66551026,   -1.99711218, -4.43951893, 0.98103164 ],
            [ -4.13586207, 6.09693725,   -1.99491422, -4.16094492, 0.97466338 ],
            [ -3.88662729, 5.60010654,   -1.99167279, -3.91889900, 0.96733630 ],
            [ -3.66090522, 5.15102874,   -1.98702390, -3.70146869, 0.95886757 ],
            [ -3.45421726, 4.74095598,   -1.98059468, -3.50424718, 0.94919022 ],
            [ -3.26208478, 4.36120333,   -1.97192298, -3.32290111, 0.93817180 ],
            [ -3.08025312, 4.00363589,   -1.96040614, -3.15342095, 0.92559287 ],
            [ -2.90536330, 3.66203670,   -1.94529617, -2.99276889, 0.91118137 ],
            [ -2.73271279, 3.32781649,   -1.92543106, -2.83686239, 0.89440684 ],
            [ -2.55558367, 2.98902393,   -1.89879305, -2.68017529, 0.87425921 ],
            [ -2.35695721, 2.61555247,   -1.86010071, -2.50907994, 0.84783233 ],
            [ -2.17584042, 2.28252561,   -1.81594644, -2.35798345, 0.82008810 ],
            [ -2.00945997, 1.98428567,   -1.76794113, -2.22386599, 0.79163609 ],
            [ -1.85293247, 1.71148075,   -1.71690125, -2.10220296, 0.76254469 ],
            [ -1.70312288, 1.45821004,   -1.66373018, -1.99016326, 0.73297210 ],
            [ -1.55726214, 1.21949402,   -1.60911445, -1.88542574, 0.70300634 ],
            [ -1.41327475, 0.99177984,   -1.55371994, -1.78637295, 0.67278336 ],
            [ -1.26962657, 0.77258133,   -1.49821623, -1.69190534, 0.64249971 ],
            [ -1.12470143, 0.55946379,   -1.44308091, -1.60098618, 0.61230375 ],
            [ -0.97735508, 0.35085421,   -1.38884481, -1.51297946, 0.58242425 ],
            [ -0.82649462, 0.14535954,   -1.33594870, -1.42734580, 0.55308494 ],
            [ -0.67075418, -0.05866151,  -1.28465740, -1.34346161, 0.52444634 ],
            [ -0.50903424, -0.26236855,  -1.23527138, -1.26091934, 0.49671259 ],
            [ -0.34010414, -0.46698863,  -1.18800410, -1.17929502, 0.47005401 ],
            [ -0.16252716, -0.67388461,  -1.14298957, -1.09812276, 0.44460231 ],
            [ 0.02508017,  -0.88424325,  -1.10036826, -1.01701720, 0.42049273 ],
            [ 0.22415449,  -1.09922251,  -1.06025741, -0.93561052, 0.39784291 ],
            [ 0.43680615,  -1.32059897,  -1.02264960, -0.85331018, 0.37669272 ],
            [ 0.66498656,  -1.54985122,  -0.98759892, -0.76964796, 0.35710876 ],
            [ 0.90632907,  -1.78425183,  -0.95569280, -0.68565256, 0.33944317 ],
            [ 1.15638633,  -2.01960051,  -0.92741038, -0.60276272, 0.32396427 ],
            [ 1.41764458,  -2.25851689,  -0.90223986, -0.51995135, 0.31038178 ],
            [ 1.69072466,  -2.50176994,  -0.87993342, -0.43685822, 0.29854703 ],
            [ 1.97905182,  -2.75253262,  -0.86007264, -0.35231626, 0.28822042 ],
            [ 2.28366620,  -3.01176494,  -0.84249101, -0.26592606, 0.27929698 ],
            [ 2.60825033,  -3.28261617,  -0.82691145, -0.17656307, 0.27161626 ],
            [ 2.94944735,  -3.56237494,  -0.81339521, -0.08502866, 0.26518328 ],
            [ 3.31304190,  -3.85590460,  -0.80160736, 0.01037221,  0.25980876 ],
            [ 3.69635757,  -4.16117304,  -0.79153290, 0.10908836,  0.25545464 ],
            [ 4.18944995,  -4.54885341,  -0.78137650, 0.23398498,  0.25138035 ],
            [ 4.72253867,  -4.96308664,  -0.77311828, 0.36714679,  0.24841718 ],
            [ 5.29921555,  -5.40693249,  -0.76655582, 0.50977632,  0.24641708 ],
            [ 5.92597639,  -5.88568965,  -0.76146230, 0.66380347,  0.24522391 ],
            [ 6.61013570,  -6.40525340,  -0.75762570, 0.83135655,  0.24468873 ],
            [ 7.36296489,  -6.97448828,  -0.75483280, 1.01552940,  0.24466772 ],
            [ 8.20022435,  -7.60559570,  -0.75288352, 1.22051119,  0.24502791 ],
            [ 9.14460284,  -8.31593542,  -0.75159329, 1.45219362,  0.24564868 ],
            [ 10.23413446, -9.13433565,  -0.75079385, 1.72026015,  0.24642587 ],
            [ 11.53105459, -10.10772212, -0.75034167, 2.04041735,  0.24727073 ],
            [ 13.14739246, -11.32031560, -0.75011772, 2.44080080,  0.24810808 ],
            [ 15.32296547, -12.95213077, -0.75002768, 2.98147882,  0.24887505 ],
            [ 18.71298738, -15.49468442, -0.75000291, 3.82639429,  0.24951156 ],
            [ 26.02173536, -20.97624980, -0.75000002, 5.65193815,  0.24992082 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 9e-07
        #    Est Max FD Error (based on N/2 run) to R=50.01: 5.e-7
        #    Est Max MOC error for R>50.01: -5.76e-07
        #    Max interpolation error with cubic splines: 1.0024e-06
        #    Est Max overall error after cubic interpolation: 2.2e-6
        C4000_G1X667 => [
            [ -6.55083431, 12.00033643,  -1.99998464, -6.55222075, 0.99861261 ],
            [ -6.08707399, 11.07282703,  -1.99996331, -6.08927938, 0.99779221 ],
            [ -5.31249797, 9.52374222,   -1.99982841, -5.31728897, 0.99519793 ],
            [ -4.77775165, 8.45441348,   -1.99950101, -4.78594317, 0.99177696 ],
            [ -4.35429884, 7.60784016,   -1.99883749, -4.36683414, 0.98739338 ],
            [ -4.00765424, 6.91513187,   -1.99767905, -4.02542409, 0.98209284 ],
            [ -3.70982651, 6.32041991,   -1.99580154, -3.73382392, 0.97576512 ],
            [ -3.44872341, 5.79965009,   -1.99295339, -3.47997009, 0.96837593 ],
            [ -3.21538770, 5.33506724,   -1.98883415, -3.25496704, 0.95986099 ],
            [ -3.00310147, 4.91343459,   -1.98308024, -3.05219803, 0.95012134 ],
            [ -2.80738611, 4.52603231,   -1.97527724, -2.86729108, 0.93906127 ],
            [ -2.62376300, 4.16422657,   -1.96488737, -2.69597613, 0.92649452 ],
            [ -2.44864925, 3.82128182,   -1.95123838, -2.53495160, 0.91218249 ],
            [ -2.27840065, 3.49053159,   -1.93341913, -2.38101439, 0.89576080 ],
            [ -2.10776552, 3.16253781,   -1.90993504, -2.22976623, 0.87651902 ],
            [ -1.92772530, 2.82144075,   -1.87783835, -2.07402899, 0.85292800 ],
            [ -1.74007144, 2.47289684,   -1.83530631, -1.91657270, 0.82458823 ],
            [ -1.57002619, 2.16469034,   -1.78836388, -1.77877258, 0.79566139 ],
            [ -1.41193117, 1.88586736,   -1.73788436, -1.65527764, 0.76625304 ],
            [ -1.26178900, 1.62887293,   -1.68470908, -1.54244713, 0.73646607 ],
            [ -1.11843138, 1.39122016,   -1.63033648, -1.43898475, 0.70679056 ],
            [ -0.97668329, 1.16406839,   -1.57442700, -1.34092194, 0.67675895 ],
            [ -0.83545185, 0.94569765,   -1.51793062, -1.24746839, 0.64666432 ],
            [ -0.69370376, 0.73453262,   -1.46168123, -1.15792954, 0.61677969 ],
            [ -0.54983165, 0.52825726,   -1.40613946, -1.07132983, 0.58722849 ],
            [ -0.40266720, 0.32535581,   -1.35181869, -0.98706306, 0.55820859 ],
            [ -0.25097976, 0.12434473,   -1.29910765, -0.90455953, 0.52989195 ],
            [ -0.09371172, -0.07592005,  -1.24837189, -0.82340756, 0.50247126 ],
            [ 0.07076201,  -0.27718288,  -1.19973957, -0.74297008, 0.47603535 ],
            [ 0.24341252,  -0.48025751,  -1.15350299, -0.66299913, 0.45077833 ],
            [ 0.42591597,  -0.68670424,  -1.10972784, -0.58296159, 0.42678011 ],
            [ 0.61975257,  -0.89773340,  -1.06853088, -0.50247485, 0.40415164 ],
            [ 0.82683369,  -1.11491924,  -1.02994004, -0.42102784, 0.38295412 ],
            [ 1.04923436,  -1.33988163,  -0.99397934, -0.33810547, 0.36324276 ],
            [ 1.28330524,  -1.56863564,  -0.96142747, -0.25521535, 0.34547815 ],
            [ 1.52647486,  -1.79881613,  -0.93250410, -0.17316467, 0.32979951 ],
            [ 1.78052084,  -2.03235780,  -0.90677275, -0.09118683, 0.31597776 ],
            [ 2.04735065,  -2.27117063,  -0.88387019, -0.00854573, 0.30381877 ],
            [ 2.32833705,  -2.51658743,  -0.86354168, 0.07528168,  0.29318324 ],
            [ 2.62526905,  -2.77024806,  -0.84554198, 0.16091829,  0.28393460 ],
            [ 2.93981833,  -3.03363897,  -0.82967277, 0.24893114,  0.27595932 ],
            [ 3.27363011,  -3.30819634,  -0.81576325, 0.33987152,  0.26915655 ],
            [ 3.62871082,  -3.59563407,  -0.80364959, 0.43438569,  0.26342779 ],
            [ 4.00634682,  -3.89707843,  -0.79320607, 0.53293213,  0.25869139 ],
            [ 4.48848962,  -4.27688144,  -0.78274806, 0.65651783,  0.25421773 ],
            [ 5.00941084,  -4.68230446,  -0.77422977, 0.78801675,  0.25087223 ],
            [ 5.57996809,  -5.12198362,  -0.76736521, 0.93041956,  0.24848538 ],
            [ 6.19423559,  -5.59163290,  -0.76207795, 1.08254142,  0.24695761 ],
            [ 6.86479443,  -6.10122349,  -0.75807712, 1.24781993,  0.24611470 ],
            [ 7.60254436,  -6.65933705,  -0.75514872, 1.42924940,  0.24581487 ],
            [ 8.42127761,  -7.27669311,  -0.75309447, 1.63052881,  0.24592596 ],
            [ 9.34405361,  -7.97093883,  -0.75172353, 1.85763463,  0.24632764 ],
            [ 10.40530898, -8.76820501,  -0.75086723, 2.11935833,  0.24691443 ],
            [ 11.66205604, -9.71150685,  -0.75037813, 2.43010336,  0.24759616 ],
            [ 13.21789561, -10.87875161, -0.75013252, 2.81589289,  0.24829793 ],
            [ 15.28662052, -12.43044188, -0.75003195, 3.33028488,  0.24895718 ],
            [ 18.44421755, -14.79868035, -0.75000347, 4.11738923,  0.24951956 ],
            [ 25.46974999, -20.06783443, -0.75000001, 5.87217686,  0.24991685 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 4.15e-07
        #    Est Max FD Error (based on N/2 run) to R=100.03: 6.3e-7
        #    Est Max MOC error for R>100.03: -5.79e-07
        #    Max interpolation error with cubic splines: 1.0022e-06
        #    Est Max overall error after cubic interpolation: 1.5e-6
        'P4000_G1X667' => [
            [
                -11.89915783, 12.00011376, -0.99997729, -11.90100705,
                0.99907454
            ],
            [
                -11.61073518, 11.71169387, -0.99998167, -11.61287158,
                0.99893067
            ],
            [
                -11.19332851, 11.29429258, -0.99998283, -11.19596136,
                0.99868186
            ],
            [ -9.62894918, 9.72997606,  -0.99992051, -9.63471407, 0.99710950 ],
            [ -8.54109423, 8.64227741,  -0.99976436, -8.55104563, 0.99500083 ],
            [ -7.68814249, 7.78964279,  -0.99944767, -7.70342432, 0.99230535 ],
            [ -6.98343282, 7.08549781,  -0.99888444, -7.00523349, 0.98899427 ],
            [ -6.38217910, 6.48516249,  -0.99797040, -6.41172334, 0.98504291 ],
            [ -5.85594214, 5.96032818,  -0.99657918, -5.89452013, 0.98041246 ],
            [ -5.38635963, 5.49278960,  -0.99456227, -5.43534128, 0.97505867 ],
            [ -4.96021966, 5.06952652,  -0.99174403, -5.02108874, 0.96892297 ],
            [ -4.56722142, 4.68047955,  -0.98791213, -4.64163260, 0.96192258 ],
            [ -4.19930341, 4.31789497,  -0.98280983, -4.28914444, 0.95395083 ],
            [ -3.84885198, 3.97458273,  -0.97610291, -3.95638017, 0.94484379 ],
            [ -3.50847744, 3.64376273,  -0.96734020, -3.63651476, 0.93435955 ],
            [ -3.16855076, 3.31681007,  -0.95581631, -3.32093599, 0.92206696 ],
            [ -2.81115748, 2.97789684,  -0.94010631, -2.99402201, 0.90698749 ],
            [ -2.43338653, 2.62659643,  -0.91893753, -2.65479446, 0.88851843 ],
            [ -2.09164245, 2.31643900,  -0.89553840, -2.35432449, 0.86959058 ],
            [ -1.77410883, 2.03599076,  -0.87036022, -2.08122651, 0.85027136 ],
            [ -1.48022994, 1.78394964,  -0.84456206, -1.83413420, 0.83114619 ],
            [ -1.20083778, 1.55160997,  -0.81840897, -1.60455881, 0.81213677 ],
            [ -0.92540311, 1.32985176,  -0.79174484, -1.38350813, 0.79291414 ],
            [ -0.64567389, 1.11219486,  -0.76447936, -1.16446022, 0.77322814 ],
            [ -0.35667494, 0.89527729,  -0.73682946, -0.94392086, 0.75305712 ],
            [ -0.06146564, 0.68178799,  -0.70977657, -0.72459469, 0.73295530 ],
            [ 0.24386045,  0.46912422,  -0.68358864, -0.50387662, 0.71299165 ],
            [ 0.56287426,  0.25511104,  -0.65854405, -0.27960122, 0.69327293 ],
            [ 0.90005907,  0.03713666,  -0.63484408, -0.04915388, 0.67387487 ],
            [ 1.26046439,  -0.18757234, -0.61267024, 0.19023703,  0.65488231 ],
            [ 1.65041722,  -0.42237388, -0.59216551, 0.44193273,  0.63636970 ],
            [ 2.05966690,  -0.66092077, -0.57417763, 0.69876516,  0.61913040 ],
            [ 2.46896554,  -0.89279664, -0.55935016, 0.94900435,  0.60397168 ],
            [ 2.88629898,  -1.12356709, -0.54700608, 1.19817739,  0.59045688 ],
            [ 3.31707272,  -1.35690628, -0.53672080, 1.44985724,  0.57833478 ],
            [ 3.76575643,  -1.59573492, -0.52818751, 1.70684057,  0.56743802 ],
            [ 4.23627455,  -1.84253626, -0.52116828, 1.97146677,  0.55764997 ],
            [ 4.73346992,  -2.10017401, -0.51545382, 2.24648276,  0.54886506 ],
            [ 5.36414395,  -2.42348112, -0.51013244, 2.58962776,  0.53964643 ],
            [ 6.04679383,  -2.77027715, -0.50615551, 2.95518199,  0.53164552 ],
            [ 6.79386816,  -3.14726127, -0.50328414, 3.34967992,  0.52475640 ],
            [ 7.62287796,  -3.56359906, -0.50130747, 3.78216139,  0.51888406 ],
            [ 8.55140575,  -4.02843224, -0.50004908, 4.26156544,  0.51398094 ],
            [ 9.61748529,  -4.56109380, -0.49933825, 4.80723120,  0.50994549 ],
            [ 10.82591642, -5.16429340, -0.49904267, 5.42144613,  0.50681377 ],
            [ 12.34261799, -5.92115225, -0.49903138, 6.18806980,  0.50430334 ],
            [ 14.26643862, -6.88136373, -0.49921942, 7.15631543,  0.50246864 ],
            [ 17.04199333, -8.26741315, -0.49952509, 8.54888611,  0.50116210 ],
            [ 20.00009956, -9.74542025, -0.49974825, 10.03033999, 0.50054312 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 2.91e-07
        #    Est Max FD Error (based on N/2 run) to R=100.04: 6.9e-7
        #    Est Max MOC error for R>100.04: 1.49e-06
        #    Max interpolation error with cubic splines: 1.0037e-06
        #    Est Max overall error after cubic interpolation: 2.7e-6
        'P4000_G1X3' => [
            [
                -12.64183430, 12.00011381, -0.99999501, -12.64359272,
                0.99912003
            ],
            [
                -11.47627542, 10.83457292, -0.99997394, -11.47942692,
                0.99842181
            ],
            [ -10.12239088, 9.48076425, -0.99989767, -10.12860179, 0.99688523 ],
            [ -9.10081939,  8.45937461, -0.99971588, -9.11119098,  0.99478886 ],
            [ -8.28593585,  7.64484846, -0.99935922, -8.30156185,  0.99213122 ],
            [ -7.60429104,  6.96382915, -0.99873558, -7.62632548,  0.98887610 ],
            [ -7.01983212,  6.38037260, -0.99773881, -7.04944032,  0.98501221 ],
            [ -6.50568741,  5.86774426, -0.99623676, -6.54411077,  0.98049624 ],
            [ -6.04506213,  5.40931316, -0.99407586, -6.09362169,  0.97528490 ],
            [ -5.62541614,  4.99274135, -0.99107295, -5.68555401,  0.96931691 ],
            [ -5.23689651,  4.60843253, -0.98700574, -5.31023423,  0.96250614 ],
            [ -4.87144174,  4.24866010, -0.98159969, -4.95985781,  0.95473480 ],
            [ -4.52166211,  3.90649562, -0.97449934, -4.62742309,  0.94583114 ],
            [ -4.17927219,  3.57435546, -0.96519433, -4.30529854,  0.93551074 ],
            [ -3.83307545,  3.24225188, -0.95284924, -3.98349366,  0.92325504 ],
            [ -3.45821559,  2.88817588, -0.93553019, -3.64023762,  0.90772313 ],
            [ -3.08805827,  2.54573757, -0.91393487, -3.30744926,  0.88996477 ],
            [ -2.75037089,  2.24100697, -0.89025626, -3.00995075,  0.87169135 ],
            [ -2.43462872,  1.96383466, -0.86496428, -2.73763863,  0.85296431 ],
            [ -2.14269863,  1.71501690, -0.83937361, -2.49130934,  0.83445759 ],
            [ -1.86144177,  1.48257553, -0.81333130, -2.25922007,  0.81580339 ],
            [ -1.58144658,  1.25856328, -0.78672837, -2.03346258,  0.79671704 ],
            [ -1.29470674,  1.03689028, -0.75949104, -1.80784420,  0.77695272 ],
            [ -1.00109166,  0.81791000, -0.73229249, -1.58268205,  0.75681395 ],
            [ -0.70019768,  0.60160139, -0.70574659, -1.35801696,  0.73659926 ],
            [ -0.38833557,  0.38555585, -0.68012390, -1.13147359,  0.71639125 ],
            [ -0.06146564,  0.16731227, -0.65565412, -0.90062733,  0.69627245 ],
            [ 0.28453340,  -0.05546480, -0.63256248, -0.66320652, 0.67635391 ],
            [ 0.65535761,  -0.28593507, -0.61098937, -0.41609677, 0.65670391 ],
            [ 1.05729593,  -0.52739996, -0.59108973, -0.15608611, 0.63742914 ],
            [ 1.47161193,  -0.76864296, -0.57398656, 0.10426684,  0.61971011 ],
            [ 1.88782571,  -1.00450154, -0.55983362, 0.35885576,  0.60397867 ],
            [ 2.31310983,  -1.23999330, -0.54803225, 0.61265120,  0.58987000 ],
            [ 2.75328419,  -1.47897390, -0.53817516, 0.86943126,  0.57715118 ],
            [ 3.21204925,  -1.72391905, -0.52998960, 1.13151582,  0.56570229 ],
            [ 3.69432257,  -1.97782159, -0.52323188, 1.40178908,  0.55540417 ],
            [ 4.20402854,  -2.24304741, -0.51771484, 1.68246606,  0.54618911 ],
            [ 4.74629569,  -2.52252301, -0.51326532, 1.97635417,  0.53798741 ],
            [ 5.44193612,  -2.87805198, -0.50915277, 2.34751961,  0.52947152 ],
            [ 6.20316322,  -3.26439548, -0.50611036, 2.74768999,  0.52223136 ],
            [ 7.04817345,  -3.69107007, -0.50392582, 3.18629144,  0.51616813 ],
            [ 7.99420578,  -4.16702943, -0.50242579, 3.67213680,  0.51123374 ],
            [ 9.08034898,  -4.71214539, -0.50143584, 4.22513965,  0.50731383 ],
            [ 10.28919468, -5.31790517, -0.50083787, 4.83655099,  0.50446466 ],
            [ 11.83920472, -6.09387162, -0.50044857, 5.61663406,  0.50230574 ],
            [ 13.93580441, -7.14281551, -0.50019804, 6.66807428,  0.50089479 ],
            [ 16.00015763, -8.17526659, -0.50008249, 7.70143149,  0.50033535 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 4.8e-07
        #    Est Max FD Error (based on N/2 run) to R=70.01: 7.e-7
        #    Est Max MOC error for R>70.01: 3.11e-07
        #    Max interpolation error with cubic splines: 1.0023e-06
        #    Est Max overall error after cubic interpolation: 1.5e-6
        'P4000_G1X2' => [
            [
                -13.01929359, 12.00011385, -0.99999415, -13.02102096,
                0.99913557
            ],
            [
                -11.82708696, 10.80792625, -0.99997322, -11.83022438,
                0.99842887
            ],
            [ -10.53366369, 9.51457558, -0.99990003, -10.53966215, 0.99699209 ],
            [ -9.50940503,  8.49049553, -0.99972143, -9.51943470,  0.99496145 ],
            [ -8.68802751,  7.66947259, -0.99936741, -8.70318665,  0.99236792 ],
            [ -8.00436732,  6.98643240, -0.99874935, -8.02576338,  0.98920136 ],
            [ -7.41766910,  6.40073062, -0.99775856, -7.44644943,  0.98543625 ],
            [ -6.90240474,  5.88697350, -0.99626587, -6.93977128,  0.98104050 ],
            [ -6.44065979,  5.42741220, -0.99411604, -6.48790569,  0.97596505 ],
            [ -6.01989416,  5.00970951, -0.99112564, -6.07843211,  0.97014941 ],
            [ -5.63025488,  4.62427006, -0.98707218, -5.70167424,  0.96350868 ],
            [ -5.26368045,  4.26336834, -0.98168138, -5.34982310,  0.95592684 ],
            [ -4.91233330,  3.91964044, -0.97458858, -5.01544552,  0.94722284 ],
            [ -4.56815194,  3.58573168, -0.96528356, -4.69112024,  0.93711780 ],
            [ -4.21949196,  3.25123833, -0.95291219, -4.36642765,  0.92508141 ],
            [ -3.84015349,  2.89293136, -0.93546304, -4.01834414,  0.90972533 ],
            [ -3.46910045,  2.54968018, -0.91392625, -3.68394590,  0.89230365 ],
            [ -3.13029342,  2.24392888, -0.89033489, -3.38461686,  0.87433608 ],
            [ -2.81320766,  1.96553616, -0.86515178, -3.11026519,  0.85587241 ],
            [ -2.51963461,  1.71524833, -0.83966676, -2.86166625,  0.83756215 ],
            [ -2.23630163,  1.48099670, -0.81370865, -2.62696597,  0.81902930 ],
            [ -1.95392512,  1.25496278, -0.78717837, -2.39836885,  0.79999586 ],
            [ -1.66460675,  1.03115575, -0.76001491, -2.16977312,  0.78022085 ],
            [ -1.36905170,  0.81055223, -0.73296205, -1.94215933,  0.76005844 ],
            [ -1.06625547,  0.59265026, -0.70656699, -1.71510230,  0.73976186 ],
            [ -0.75259745,  0.37507841, -0.68109934, -1.48628223,  0.71941691 ],
            [ -0.42394054,  0.15529653, -0.65677467, -1.25321059,  0.69909965 ],
            [ -0.07602113, -0.06912890, -0.63380462, -1.01353328, 0.67891525 ],
            [ 0.29662565,  -0.30121656, -0.61234238, -0.76431223, 0.65894953 ],
            [ 0.70031296,  -0.54429847, -0.59253586, -0.50233589, 0.63931308 ],
            [ 1.11764639,  -0.78790533, -0.57544624, -0.23939164, 0.62115938 ],
            [ 1.53610395,  -1.02565185, -0.56131630, 0.01709407,  0.60503855 ],
            [ 1.96472032,  -1.26361945, -0.54949104, 0.27324412,  0.59052488 ],
            [ 2.40693676,  -1.50435350, -0.53962887, 0.53142972,  0.57747180 ],
            [ 2.86905632,  -1.75175150, -0.53139867, 0.79549952,  0.56568897 ],
            [ 3.35424520,  -2.00786151, -0.52459289, 1.06733021,  0.55511154 ],
            [ 3.86742094,  -2.27557463, -0.51900961, 1.34970291,  0.54565427 ],
            [ 4.41327256,  -2.55758341, -0.51448200, 1.64518656,  0.53726081 ],
            [ 5.11344718,  -2.91624377, -0.51025984, 2.01820098,  0.52858358 ],
            [ 5.88013257,  -3.30615687, -0.50708774, 2.42051918,  0.52125102 ],
            [ 6.73073710,  -3.73642564, -0.50475714, 2.86117878,  0.51517255 ],
            [ 7.68877224,  -4.21913755, -0.50308659, 3.35224257,  0.51026836 ],
            [ 8.78566491,  -4.77027432, -0.50191873, 3.90971996,  0.50646709 ],
            [ 9.83848027,  -5.29831255, -0.50123120, 4.44159275,  0.50408322 ],
            [ 11.35838846, -6.05964917, -0.50064793, 5.20603958,  0.50204640 ],
            [ 13.26379683, -7.01319121, -0.50027885, 6.16131860,  0.50082927 ],
            [ 15.54447178, -8.15392038, -0.50009524, 7.30279854,  0.50027196 ],
            [ 16.00014827, -8.38179758, -0.50007638, 7.53074775,  0.50021712 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 7.97e-07
        #    Est Max FD Error (based on N/2 run) to R=20.00: 6.7e-7
        #    Est Max MOC error for R>20.00: -2.16e-07
        #    Max interpolation error with cubic splines: 1.0114e-06
        #    Est Max overall error after cubic interpolation: 1.6e-6
        S4000_G1X3 => [
            [ -4.70228458, 12.00033614,  -2.99997917, -4.70333915, 0.99841733 ],
            [ -4.25375155, 10.65475547,  -2.99992756, -4.25581922, 0.99689536 ],
            [ -3.80521852, 9.30922741,   -2.99971455, -3.80927445, 0.99390429 ],
            [ -3.49149171, 8.36819737,   -2.99926081, -3.49799250, 0.99021932 ],
            [ -3.21090665, 7.52676738,   -2.99828628, -3.22082468, 0.98505695 ],
            [ -2.97645129, 6.82398369,   -2.99654454, -2.99057510, 0.97868704 ],
            [ -2.78140771, 6.23976761,   -2.99381401, -2.80036921, 0.97134206 ],
            [ -2.60472645, 5.71116151,   -2.98954180, -2.62949686, 0.96250390 ],
            [ -2.44887074, 5.24566535,   -2.98342132, -2.48023721, 0.95245391 ],
            [ -2.30667480, 4.82200342,   -2.97484498, -2.34559071, 0.94094792 ],
            [ -2.17516862, 4.43151367,   -2.96316875, -2.22268239, 0.92786090 ],
            [ -2.05223071, 4.06812736,   -2.94767611, -2.10949471, 0.91307077 ],
            [ -1.93444320, 3.72207001,   -2.92724104, -2.00291306, 0.89617963 ],
            [ -1.81889500, 3.38531490,   -2.90030956, -1.90045982, 0.87664197 ],
            [ -1.70267501, 3.05022110,   -2.86466283, -1.79987742, 0.85367413 ],
            [ -1.57808767, 2.69627758,   -2.81509292, -1.69525966, 0.82504283 ],
            [ -1.45223866, 2.34584821,   -2.75163142, -1.59347788, 0.79174577 ],
            [ -1.33792540, 2.03517645,   -2.68192826, -1.50487966, 0.75778112 ],
            [ -1.23133442, 1.75321028,   -2.60723272, -1.42592564, 0.72323148 ],
            [ -1.12985182, 1.49254939,   -2.52877310, -1.35429187, 0.68821339 ],
            [ -1.03221807, 1.24955344,   -2.44827061, -1.28880377, 0.65311015 ],
            [ -0.93570397, 1.01722649,   -2.36577156, -1.22747744, 0.61763951 ],
            [ -0.83970651, 0.79410332,   -2.28274151, -1.16988765, 0.58220340 ],
            [ -0.74319241, 0.57778986,   -2.20004862, -1.11540109, 0.54700383 ],
            [ -0.64522608, 0.36628118,   -2.11843836, -1.06352562, 0.51224626 ],
            [ -0.54512909, 0.15825762,   -2.03869830, -1.01396920, 0.47820154 ],
            [ -0.44207897, -0.04780100,  -1.96135129, -0.96641571, 0.44506511 ],
            [ -0.33517765, -0.25343257,  -1.88676458, -0.92057361, 0.41299229 ],
            [ -0.22343627, -0.46020232,  -1.81518904, -0.87617631, 0.38210486 ],
            [ -0.10625182, -0.66885467,  -1.74706659, -0.83315608, 0.35261559 ],
            [ 0.01735789,  -0.88074626,  -1.68254356, -0.79133198, 0.32461236 ],
            [ 0.14858123,  -1.09745872,  -1.62165742, -0.75050760, 0.29813323 ],
            [ 0.28834017,  -1.32002256,  -1.56456533, -0.71061646, 0.27326390 ],
            [ 0.43814975,  -1.55031644,  -1.51117719, -0.67146389, 0.24997617 ],
            [ 0.59713560,  -1.78657780,  -1.46213510, -0.63346543, 0.22855673 ],
            [ 0.76211776,  -2.02409600,  -1.41827986, -0.59737904, 0.20937516 ],
            [ 0.93432041,  -2.26485262,  -1.37891796, -0.56284615, 0.19212845 ],
            [ 1.11630276,  -2.51246132,  -1.34324298, -0.52934403, 0.17646043 ],
            [ 1.30726584,  -2.76583285,  -1.31121092, -0.49702835, 0.16234902 ],
            [ 1.51057581,  -3.02937679,  -1.28209730, -0.46536409, 0.14947146 ],
            [ 1.72632965,  -3.30307981,  -1.25579117, -0.43440957, 0.13777559 ],
            [ 1.95531519,  -3.58784720,  -1.23206429, -0.40410822, 0.12715855 ],
            [ 2.19969224,  -3.88623873,  -1.21058754, -0.37424811, 0.11747242 ],
            [ 2.46135185,  -4.20038345,  -1.19112846, -0.34469992, 0.10861276 ],
            [ 2.74208761,  -4.53223152,  -1.17350338, -0.31537760, 0.10049759 ],
            [ 3.04088924,  -4.88044636,  -1.15768941, -0.28647963, 0.09312155 ],
            [ 3.42815127,  -5.32537750,  -1.14072670, -0.25202334, 0.08507902 ],
            [ 3.85068842,  -5.80409202,  -1.12570309, -0.21765824, 0.07780941 ],
            [ 4.31338545,  -6.32176012,  -1.11237880, -0.18322996, 0.07121261 ],
            [ 4.82090973,  -6.88321486,  -1.10056541, -0.14865749, 0.06521428 ],
            [ 5.37871461,  -7.49408838,  -1.09008654, -0.11385302, 0.05974682 ],
            [ 5.99304054,  -8.16079608,  -1.08078196, -0.07873126, 0.05475080 ],
            [ 6.67211693,  -8.89181463,  -1.07249558, -0.04315575, 0.05016776 ],
            [ 7.42436384,  -9.69571552,  -1.06510372, -0.00705099, 0.04595507 ],
            [ 8.25839432,  -10.58119306, -1.05850395, 0.02961077,  0.04208029 ],
            [ 9.18621638,  -11.56045419, -1.05259047, 0.06694384,  0.03850614 ],
            [ 10.21783446, -12.64349474, -1.04729276, 0.10491575,  0.03521343 ],
            [ 11.36785087, -13.84506308, -1.04253234, 0.14360937,  0.03217480 ],
            [ 12.65026668, -15.17918018, -1.03825220, 0.18301720,  0.02937297 ],
            [ 14.07948230, -16.66021871, -1.03440483, 0.22309560,  0.02679389 ],
            [ 15.67469792, -18.30745547, -1.03094056, 0.26388254,  0.02441910 ],
            [ 17.41563155, -20.09950035, -1.02788356, 0.30447302,  0.02227943 ],
            [ 19.13056235, -21.86006669, -1.02541484, 0.34112873,  0.02051944 ],
            [ 21.09604290, -23.87312773, -1.02307709, 0.37974452,  0.01882482 ],
            [ 23.36199118, -26.18877967, -1.02086763, 0.42049867,  0.01719672 ],
            [ 25.99128542, -28.87010500, -1.01878388, 0.46359550,  0.01563623 ],
            [ 28.40942322, -31.33170053, -1.01720569, 0.49991682,  0.01443725 ],
            [ 31.16222141, -34.12973393, -1.01570507, 0.53802537,  0.01328277 ],
            [ 34.31305658, -37.32773797, -1.01428081, 0.57807731,  0.01217332 ],
            [ 37.94093279, -41.00489472, -1.01293169, 0.62025061,  0.01110942 ],
            [ 42.14532619, -45.26088434, -1.01165654, 0.66474922,  0.01009155 ],
            [ 45.75187495, -48.90779021, -1.01074801, 0.69978938,  0.00935865 ],
            [ 49.81956251, -53.01738214, -1.00987996, 0.73638215,  0.00865212 ],
            [ 54.42920697, -57.67060533, -1.00905190, 0.77465457,  0.00797214 ],
            [ 59.68003197, -62.96682699, -1.00826336, 0.81475015,  0.00731890 ],
            [ 65.69492021, -69.02909285, -1.00751388, 0.85683185,  0.00669259 ],
            [ 72.62750045, -76.01121783, -1.00680301, 0.90108579,  0.00609338 ],
            [ 80.67183632, -84.10747949, -1.00613027, 0.94772582,  0.00552143 ],
            [ 90.07586757, -93.56606509, -1.00549523, 0.99699936,  0.00497691 ],
            [ 97.25526514, -100.78343304, -1.00509259, 1.03145190, 0.00462921 ],
            [
                105.29779605, -108.86532726, -1.00470636, 1.06729950,
                0.00429381
            ],
            [
                114.34609987, -117.95449631, -1.00433642, 1.10465114,
                0.00397076
            ],
            [
                124.57342501, -128.22431101, -1.00398264, 1.14362916,
                0.00366010
            ],
            [
                136.19184568, -139.88698363, -1.00364489, 1.18437163,
                0.00336187
            ],
            [
                149.46317038, -153.20447925, -1.00332303, 1.22703506,
                0.00307610
            ],
            [
                164.71359568, -168.50317292, -1.00301694, 1.27179791,
                0.00280283
            ],
            [
                182.35365074, -186.19379856, -1.00272648, 1.31886489,
                0.00254211
            ],
            [
                202.90573193, -206.79898894, -1.00245152, 1.36847237,
                0.00229395
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 8e-07
        #    Est Max FD Error (based on N/2 run) to R=15.00: 9.3e-7
        #    Est Max MOC error for R>15.00: -2.25e-07
        #    Max interpolation error with cubic splines: 1.0028e-06
        #    Est Max overall error after cubic interpolation: 1.8e-6
        S4000_G1X2 => [
            [ -4.82347692, 12.00033639,  -2.99997990, -4.82451288, 0.99844528 ],
            [ -4.30866043, 10.45591053,  -2.99991193, -4.31090415, 0.99663073 ],
            [ -3.87445895, 9.15338727,   -2.99966874, -3.87876678, 0.99352499 ],
            [ -3.56050822, 8.21170962,   -2.99914308, -3.56741551, 0.98960595 ],
            [ -3.29649402, 7.42001234,   -2.99810853, -3.30677322, 0.98451072 ],
            [ -3.07144374, 6.74546919,   -2.99629043, -3.08587634, 0.97821909 ],
            [ -2.87528052, 6.15796862,   -2.99334080, -2.89468992, 0.97066203 ],
            [ -2.70352573, 5.64419676,   -2.98890895, -2.72869354, 0.96190049 ],
            [ -2.55013326, 5.18617189,   -2.98255523, -2.58188443, 0.95187104 ],
            [ -2.40905696, 4.76599097,   -2.97363791, -2.44838304, 0.94032949 ],
            [ -2.27862456, 4.37887240,   -2.96155893, -2.32655993, 0.92723140 ],
            [ -2.15591058, 4.01637929,   -2.94549228, -2.21366098, 0.91235399 ],
            [ -2.03767521, 3.66931107,   -2.92423233, -2.10676971, 0.89527327 ],
            [ -1.92235094, 3.33360292,   -2.89646074, -2.00462514, 0.87563968 ],
            [ -1.80501130, 2.99581285,   -2.85939564, -1.90321263, 0.85229776 ],
            [ -1.67714691, 2.63340275,   -2.80705312, -1.79608006, 0.82267665 ],
            [ -1.55278351, 2.28815902,   -2.74288379, -1.69578587, 0.78953391 ],
            [ -1.43925069, 1.98063465,   -2.67267223, -1.60803957, 0.75565614 ],
            [ -1.33299892, 1.70056908,   -2.59767421, -1.52956139, 0.72113947 ],
            [ -1.23155835, 1.44099432,   -2.51910785, -1.45817044, 0.68611054 ],
            [ -1.13343470, 1.19773992,   -2.43839985, -1.39256828, 0.65084500 ],
            [ -1.03647275, 0.96528925,   -2.35599923, -1.33118292, 0.61525761 ],
            [ -0.93995865, 0.74189741,   -2.27323128, -1.27351867, 0.57971105 ],
            [ -0.84284153, 0.52513957,   -2.19090159, -1.21893929, 0.54440264 ],
            [ -0.74431206, 0.31329109,   -2.10981862, -1.16702512, 0.50958532 ],
            [ -0.64365856, 0.10494774,   -2.03069774, -1.11746243, 0.47551726 ],
            [ -0.53975475, -0.10201012,  -1.95380801, -1.06979851, 0.44229842 ],
            [ -0.43200209, -0.30849421,  -1.87973599, -1.02389178, 0.41018834 ],
            [ -0.31950251, -0.51591076,  -1.80876849, -0.97950722, 0.37933111 ],
            [ -0.20149107, -0.72531241,  -1.74122362, -0.93650810, 0.34989016 ],
            [ -0.07676171, -0.93842215,  -1.67713952, -0.89464419, 0.32190640 ],
            [ 0.05553475,  -1.15622435,  -1.61673660, -0.85383950, 0.29549772 ],
            [ 0.19661104,  -1.38022030,  -1.56004234, -0.81393996, 0.27069132 ],
            [ 0.34768179,  -1.61180330,  -1.50709092, -0.77483842, 0.24751107 ],
            [ 0.50868790,  -1.85042633,  -1.45826461, -0.73675150, 0.22612727 ],
            [ 0.67461559,  -2.08870636,  -1.41491049, -0.70084596, 0.20712894 ],
            [ 0.84810712,  -2.33071348,  -1.37592006, -0.66643163, 0.19002656 ],
            [ 1.03042421,  -2.57827958,  -1.34075824, -0.63323036, 0.17457971 ],
            [ 1.22373972,  -2.83430390,  -1.30884920, -0.60087444, 0.16052895 ],
            [ 1.42704969,  -3.09741548,  -1.28018751, -0.56955762, 0.14786671 ],
            [ 1.64392318,  -3.37215266,  -1.25413003, -0.53877576, 0.13630356 ],
            [ 1.87461375,  -3.65867659,  -1.23056593, -0.50857663, 0.12578636 ],
            [ 2.12154045,  -3.95982195,  -1.20917441, -0.47873527, 0.11616901 ],
            [ 2.38342400,  -4.27389976,  -1.18996659, -0.44948324, 0.10745635 ],
            [ 2.66669624,  -4.60842803,  -1.17241759, -0.42021366, 0.09941018 ],
            [ 2.96763204,  -4.95881822,  -1.15670148, -0.39142430, 0.09211354 ],
            [ 3.35765633,  -5.40655847,  -1.13984155, -0.35709766, 0.08415945 ],
            [ 3.78356588,  -5.88873397,  -1.12489627, -0.32283337, 0.07696567 ],
            [ 4.24944242,  -6.40960378,  -1.11165613, -0.28854305, 0.07044609 ],
            [ 4.76035254,  -6.97445508,  -1.09991873, -0.25411271, 0.06451927 ],
            [ 5.32194811,  -7.58913588,  -1.08950556, -0.21944338, 0.05911618 ],
            [ 5.94106829,  -8.26070288,  -1.08025109, -0.18442076, 0.05417422 ],
            [ 6.62474116,  -8.99632566,  -1.07201871, -0.14897949, 0.04964544 ],
            [ 7.38118607,  -9.80437413,  -1.06468284, -0.11304725, 0.04548661 ],
            [ 8.22101556,  -10.69567257, -1.05812368, -0.07650515, 0.04165533 ],
            [ 9.15483714,  -11.68092978, -1.05224952, -0.03930754, 0.03812254 ],
            [ 10.19365502, -12.77119237, -1.04698447, -0.00144980, 0.03486583 ],
            [ 11.35167136, -13.98078005, -1.04225379, 0.03713018,  0.03186019 ],
            [ 12.64188713, -15.32267247, -1.03800392, 0.07639338,  0.02909068 ],
            [ 14.08210276, -16.81476805, -1.03417782, 0.11639252,  0.02653697 ],
            [ 15.68691836, -18.47157433, -1.03073851, 0.15703507,  0.02418919 ],
            [ 17.28440808, -20.11586598, -1.02794542, 0.19407807,  0.02224374 ],
            [ 19.13602403, -22.01668249, -1.02528730, 0.23346570,  0.02035797 ],
            [ 21.28120075, -24.21332286, -1.02278136, 0.27513312,  0.01854762 ],
            [ 23.33854390, -26.31544220, -1.02080689, 0.31175865,  0.01709763 ],
            [ 25.68985635, -28.71340645, -1.01893387, 0.35027155,  0.01570180 ],
            [ 28.39251850, -31.46476495, -1.01716043, 0.39084082,  0.01436090 ],
            [ 31.51850125, -34.64168175, -1.01548478, 0.43365962,  0.01307571 ],
            [ 35.15903550, -38.33561000, -1.01390517, 0.47895022,  0.01184697 ],
            [ 38.51947891, -41.74071144, -1.01270946, 0.51713404,  0.01090510 ],
            [ 42.35891827, -45.62669849, -1.01157325, 0.55721503,  0.01000016 ],
            [ 46.77119878, -50.08758849, -1.01049568, 0.59936533,  0.00913252 ],
            [ 51.87416707, -55.24142515, -1.00947593, 0.64378197,  0.00830250 ],
            [ 56.24316617, -59.65019319, -1.00874855, 0.67871696,  0.00770487 ],
            [ 61.16212477, -64.61042719, -1.00805288, 0.71516253,  0.00712872 ],
            [ 66.72600420, -70.21721067, -1.00738857, 0.75324083,  0.00657418 ],
            [ 73.05106803, -76.58694503, -1.00675530, 0.79308968,  0.00604137 ],
            [ 80.28086829, -83.86333855, -1.00615272, 0.83486526,  0.00553041 ],
            [ 88.59428702, -92.22545138, -1.00558050, 0.87874565,  0.00504141 ],
            [ 98.21648031, -101.89864386, -1.00503831, 0.92493502, 0.00457450 ],
            [
                109.43398571, -113.16968931, -1.00452583, 0.97366910,
                0.00412978
            ],
            [
                117.97664826, -121.74959036, -1.00420051, 1.00770573,
                0.00384567
            ],
            [
                127.52570267, -131.33722586, -1.00388815, 1.04308661,
                0.00357150
            ],
            [
                138.24411523, -142.09566525, -1.00358865, 1.07991487,
                0.00330730
            ],
            [
                150.32909858, -154.22223673, -1.00330192, 1.11830606,
                0.00305308
            ],
            [
                164.02110847, -167.95752616, -1.00302786, 1.15839024,
                0.00280889
            ],
            [
                179.61572001, -183.59725576, -1.00276638, 1.20031452,
                0.00257473
            ],
            [
                197.47948679, -201.50814643, -1.00251737, 1.24424616,
                0.00235065
            ],
            [
                207.40138730, -211.45441978, -1.00239752, 1.26702346,
                0.00224239
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 7.47e-07
        #    Est Max FD Error (based on N/2 run) to R=50.01: 6.3e-7
        #    Est Max MOC error for R>50.01: -2.99e-07
        #    Max interpolation error with cubic splines: 1.0031e-06
        #    Est Max overall error after cubic interpolation: 1.9e-6
        C4000_G1X3 => [
            [ -6.90501884, 12.00033649,  -1.99998611, -6.90633721, 0.99868076 ],
            [ -5.96898884, 10.12831491,  -1.99991123, -5.97235386, 0.99662946 ],
            [ -5.50142170, 9.19325111,   -1.99976624, -5.50679784, 0.99461000 ],
            [ -4.98145480, 8.15352950,   -1.99934409, -4.99051300, 0.99090362 ],
            [ -4.58800870, 7.36702956,   -1.99856066, -4.60146052, 0.98646701 ],
            [ -4.25591959, 6.70352910,   -1.99720992, -4.27471344, 0.98105482 ],
            [ -3.96951232, 6.13179139,   -1.99506984, -3.99460403, 0.97465335 ],
            [ -3.71714251, 5.62866601,   -1.99187547, -3.74952866, 0.96721736 ],
            [ -3.49052472, 5.17775102,   -1.98731219, -3.53127007, 0.95867745 ],
            [ -3.28316497, 4.76627495,   -1.98099057, -3.33345518, 0.94891639 ],
            [ -3.09080856, 4.38599014,   -1.97246442, -3.15195774, 0.93781868 ],
            [ -2.90920083, 4.02874705,   -1.96115053, -2.98275596, 0.92517693 ],
            [ -2.73453494, 3.68743463,   -1.94628405, -2.82238832, 0.91069382 ],
            [ -2.56278016, 3.35474830,   -1.92678556, -2.66737562, 0.89390629 ],
            [ -2.38766641, 3.01952025,   -1.90077838, -2.51254832, 0.87388815 ],
            [ -2.19442667, 2.65565162,   -1.86363247, -2.34610470, 0.84811516 ],
            [ -2.01237107, 2.32023063,   -1.81971085, -2.19418964, 0.82018297 ],
            [ -1.84554276, 2.02054934,   -1.77179145, -2.05970772, 0.79157612 ],
            [ -1.68901526, 1.74713998,   -1.72073832, -1.93806140, 0.76238848 ],
            [ -1.53979556, 1.49429168,   -1.66757023, -1.82648500, 0.73283219 ],
            [ -1.39513638, 1.25697859,   -1.61304484, -1.72261968, 0.70301866 ],
            [ -1.25171212, 1.02960934,   -1.55737994, -1.62395045, 0.67282562 ],
            [ -1.10880236, 0.81103461,   -1.50158359, -1.52995712, 0.64261729 ],
            [ -0.96481497, 0.59883060,   -1.44617229, -1.43960049, 0.61254519 ],
            [ -0.81843312, 0.39115596,   -1.39163460, -1.35212487, 0.58279751 ],
            [ -0.66851488, 0.18655091,   -1.33841532, -1.26696065, 0.55358375 ],
            [ -0.51377882, -0.01651249,  -1.28682418, -1.18353005, 0.52507577 ],
            [ -0.35309212, -0.21924329,  -1.23715979, -1.10140389, 0.49746580 ],
            [ -0.18507541, -0.42304901,  -1.18960049, -1.02008668, 0.47089679 ],
            [ -0.00861808, -0.62890270,  -1.14437280, -0.93926819, 0.44554739 ],
            [ 0.17791612,  -0.83829631,  -1.10154715, -0.85844313, 0.42151130 ],
            [ 0.37612095,  -1.05254737,  -1.06121459, -0.77718700, 0.39889087 ],
            [ 0.58773049,  -1.27302338,  -1.02344300, -0.69506487, 0.37776985 ],
            [ 0.81496537,  -1.50148709,  -0.98823375, -0.61150418, 0.35818743 ],
            [ 1.05460870,  -1.73438467,  -0.95628833, -0.52783655, 0.34056022 ],
            [ 1.30390686,  -1.96914929,  -0.92786732, -0.44492500, 0.32503986 ],
            [ 1.56387128,  -2.20699273,  -0.90263377, -0.36224655, 0.31143722 ],
            [ 1.83676946,  -2.45017135,  -0.88018873, -0.27893135, 0.29952649 ],
            [ 2.12377974,  -2.69985684,  -0.86029373, -0.19449900, 0.28916696 ],
            [ 2.42776224,  -2.95860780,  -0.84264485, -0.10800948, 0.28018379 ],
            [ 2.74947730,  -3.22712225,  -0.82711078, -0.01915251, 0.27249137 ],
            [ 3.09000700,  -3.50639182,  -0.81353931, 0.07248921,  0.26599107 ],
            [ 3.45202957,  -3.79870141,  -0.80173460, 0.16776039,  0.26056319 ],
            [ 3.83808596,  -4.10617779,  -0.79154319, 0.26745352,  0.25610951 ],
            [ 4.24814863,  -4.42891632,  -0.78288132, 0.37171028,  0.25256091 ],
            [ 4.77306246,  -4.83750752,  -0.77432889, 0.50338772,  0.24936824 ],
            [ 5.33955023,  -5.27412501,  -0.76751954, 0.64397713,  0.24716654 ],
            [ 5.95608830,  -5.74559516,  -0.76219957, 0.79589604,  0.24579175 ],
            [ 6.62902671,  -6.25706563,  -0.75816910, 0.96102840,  0.24510034 ],
            [ 7.36768426,  -6.81592620,  -0.75522059, 1.14198783,  0.24494876 ],
            [ 8.18856680,  -7.43495278,  -0.75314520, 1.34314461,  0.24520232 ],
            [ 9.11117343,  -8.12911121,  -0.75176030, 1.56960503,  0.24573835 ],
            [ 10.17229547, -8.92630931,  -0.75089160, 1.83073988,  0.24645118 ],
            [ 11.42694735, -9.86806281,  -0.75039308, 2.14046149,  0.24724864 ],
            [ 12.97672167, -11.03077455, -0.75014052, 2.52429369,  0.24805449 ],
            [ 15.03140530, -12.57194391, -0.75003526, 3.03479501,  0.24880510 ],
            [ 18.13977974, -14.90327044, -0.75000428, 3.80928610,  0.24944258 ],
            [ 24.75488575, -19.86460618, -0.75000004, 5.46125340,  0.24989278 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 2.04e-07
        #    Est Max FD Error (based on N/2 run) to R=20.01: 1.6e-7
        #    Est Max MOC error for R>20.01: 2.87e-07
        #    Max interpolation error with cubic splines: 1.0055e-06
        #    Est Max overall error after cubic interpolation: 1.4e-6
        S8000_G1X4 => [
            [ -4.61798880, 12.00032860,  -2.99997850, -4.61906014, 0.99839214 ],
            [ -3.99663876, 10.13631445,  -2.99987271, -3.99936179, 0.99591005 ],
            [ -3.58468222, 8.90054843,   -2.99956211, -3.58973934, 0.99239616 ],
            [ -3.26779258, 7.95011114,   -2.99886804, -3.27593888, 0.98773504 ],
            [ -3.01304914, 7.18631353,   -2.99757216, -3.02500685, 0.98196964 ],
            [ -2.79749700, 6.54039198,   -2.99537598, -2.81405134, 0.97499853 ],
            [ -2.61094642, 5.98189417,   -2.99193794, -2.63289302, 0.96680137 ],
            [ -2.44578309, 5.48812279,   -2.98684127, -2.47396394, 0.95730684 ],
            [ -2.29698634, 5.04419070,   -2.97960010, -2.33229790, 0.94643689 ],
            [ -2.16071259, 4.63878487,   -2.96963322, -2.20413672, 0.93407530 ],
            [ -2.03362081, 4.26217268,   -2.95619879, -2.08628816, 0.92002285 ],
            [ -1.91347149, 3.90799763,   -2.93843433, -1.97667798, 0.90407679 ],
            [ -1.79746525, 3.56840774,   -2.91510193, -1.87282876, 0.88583216 ],
            [ -1.68280270, 3.23583193,   -2.88445102, -1.77243822, 0.86468442 ],
            [ -1.56475983, 2.89766894,   -2.84327587, -1.67182677, 0.83934943 ],
            [ -1.43162139, 2.52294997,   -2.78318908, -1.56223135, 0.80616047 ],
            [ -1.31250275, 2.19528733,   -2.71614403, -1.46818134, 0.77229407 ],
            [ -1.20274108, 1.90105882,   -2.64339737, -1.38527774, 0.73782750 ],
            [ -1.09927605, 1.63148582,   -2.56624395, -1.31072783, 0.70288338 ],
            [ -1.00084990, 1.38277298,   -2.48670935, -1.24325402, 0.66794325 ],
            [ -0.90442137, 1.14690646,   -2.40485455, -1.18053955, 0.63268216 ],
            [ -0.80857065, 0.92038256,   -2.32158887, -1.12159496, 0.59722586 ],
            [ -0.71271993, 0.70185490,   -2.23830392, -1.06604435, 0.56196036 ],
            [ -0.61588434, 0.48912064,   -2.15581627, -1.01332343, 0.52708256 ],
            [ -0.51723424, 0.28047176,   -2.07488090, -0.96302913, 0.49281582 ],
            [ -0.41587071, 0.07419068,   -1.99604178, -0.91478851, 0.45933769 ],
            [ -0.31098909, -0.13110982,  -1.91980264, -0.86833654, 0.42684177 ],
            [ -0.20186828, -0.33654999,  -1.84662603, -0.82349131, 0.39552765 ],
            [ -0.08754166, -0.54360791,  -1.77673566, -0.78001513, 0.36550689 ],
            [ 0.03292069,  -0.75356637,  -1.71034727, -0.73773911, 0.33689399 ],
            [ 0.16046037,  -0.96762519,  -1.64763833, -0.69653383, 0.30978807 ],
            [ 0.29618075,  -1.18715697,  -1.58868580, -0.65625933, 0.28424225 ],
            [ 0.44118856,  -1.41344000,  -1.53355696, -0.61681707, 0.26030086 ],
            [ 0.59691753,  -1.64816251,  -1.48219715, -0.57806291, 0.23794986 ],
            [ 0.75865264,  -1.88407299,  -1.43619232, -0.54124032, 0.21788627 ],
            [ 0.92739917,  -2.12285350,  -1.39490058, -0.50603365, 0.19983489 ],
            [ 1.10443499,  -2.36642549,  -1.35772263, -0.47213480, 0.18353531 ],
            [ 1.29103812,  -2.61656879,  -1.32417134, -0.43929881, 0.16877379 ],
            [ 1.48841569,  -2.87485976,  -1.29385792, -0.40734253, 0.15537839 ],
            [ 1.69767432,  -3.14266804,  -1.26646776, -0.37613437, 0.14320933 ],
            [ 1.92041890,  -3.42193060,  -1.24167562, -0.34550223, 0.13212172 ],
            [ 2.15816006,  -3.71438590,  -1.21922718, -0.31532583, 0.12200216 ],
            [ 2.41187283,  -4.02107492,  -1.19894719, -0.28557374, 0.11277343 ],
            [ 2.68363649,  -4.34434108,  -1.18059607, -0.25610357, 0.10432952 ],
            [ 2.97511090,  -4.68596732,  -1.16400385, -0.22685085, 0.09659661 ],
            [ 3.28522453,  -5.04456971,  -1.14913890, -0.19801328, 0.08956756 ],
            [ 3.68784605,  -5.50391394,  -1.13318232, -0.16354629, 0.08188559 ],
            [ 4.12862956,  -6.00017088,  -1.11901932, -0.12903668, 0.07491552 ],
            [ 4.61158163,  -6.53746468,  -1.10646578, -0.09443243, 0.06858517 ],
            [ 5.14104788,  -7.12025303,  -1.09534894, -0.05969004, 0.06282937 ],
            [ 5.72409682,  -7.75591246,  -1.08547194, -0.02463826, 0.05757018 ],
            [ 6.36800767,  -8.45192717,  -1.07668018, 0.01083180,  0.05275027 ],
            [ 7.07989568,  -9.21551547,  -1.06885035, 0.04676115,  0.04832807 ],
            [ 7.86941934,  -10.05653762, -1.06185640, 0.08326086,  0.04425864 ],
            [ 8.74627984,  -10.98479576, -1.05560003, 0.12037483,  0.04051028 ],
            [ 9.72163805,  -12.01155126, -1.04999325, 0.15814878,  0.03705456 ],
            [ 10.80840999, -13.14982035, -1.04495848, 0.19663149,  0.03386609 ],
            [ 12.01914422, -14.41215170, -1.04043641, 0.23579916,  0.03092764 ],
            [ 13.36780218, -15.81250678, -1.03637345, 0.27562787,  0.02822265 ],
            [ 14.87299221, -17.36959393, -1.03271494, 0.31617277,  0.02573065 ],
            [ 16.51777462, -19.06543376, -1.02948245, 0.35658587,  0.02348149 ],
            [ 18.32734551, -20.92565442, -1.02659740, 0.39716649,  0.02143408 ],
            [ 20.24527573, -22.89212464, -1.02410230, 0.43649423,  0.01963108 ],
            [ 22.46314552, -25.16074824, -1.02174710, 0.47805288,  0.01790002 ],
            [ 24.79341497, -27.53926107, -1.01972496, 0.51795194,  0.01638980 ],
            [ 27.48545611, -30.28174241, -1.01781391, 0.56006205,  0.01494100 ],
            [ 30.28095583, -33.12469425, -1.01618722, 0.60002924,  0.01369033 ],
            [ 33.50262344, -36.39594824, -1.01464737, 0.64214217,  0.01249078 ],
            [ 37.23985507, -40.18510750, -1.01319290, 0.68660779,  0.01134296 ],
            [ 41.08608460, -44.07965283, -1.01197058, 0.72830023,  0.01036660 ],
            [ 45.53081239, -48.57492804, -1.01081363, 0.77223321,  0.00943200 ],
            [ 50.01126122, -53.10161067, -1.00985419, 0.81268663,  0.00864881 ],
            [ 55.15652041, -58.29515325, -1.00894343, 0.85519633,  0.00789815 ],
            [ 61.10284752, -64.29200911, -1.00808074, 0.89995802,  0.00718028 ],
            [ 68.02307732, -71.26523853, -1.00726551, 0.94719711,  0.00649543 ],
            [ 74.89632940, -78.18608261, -1.00660404, 0.98985664,  0.00593488 ],
            [ 82.82467294, -86.16421507, -1.00597660, 1.03471778,  0.00539891 ],
            [ 92.03290078, -95.42464751, -1.00538282, 1.08199919,  0.00488765 ],
            [
                100.88931861, -104.32659604, -1.00491344, 1.12343069,
                0.00448057
            ],
            [
                111.04329492, -114.52812446, -1.00446695, 1.16688668,
                0.00409080
            ],
            [
                122.75853987, -126.29313585, -1.00404315, 1.21256004,
                0.00371843
            ],
            [
                136.36904580, -139.95584466, -1.00364181, 1.26067282,
                0.00336351
            ],
            [
                148.90568093, -152.53616799, -1.00333677, 1.30109104,
                0.00309219
            ],
            [
                163.20262019, -166.87867008, -1.00304587, 1.34338574,
                0.00283211
            ],
            [
                179.60268249, -183.32634215, -1.00276898, 1.38772918,
                0.00258330
            ],
            [
                198.53610751, -202.30962118, -1.00250601, 1.43431863,
                0.00234580
            ],
            [
                203.72621983, -207.51257358, -1.00244243, 1.44634284,
                0.00228819
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 8.3e-07
        #    Est Max FD Error (based on N/2 run) to R=20.01: 1.4e-6
        #    Est Max MOC error for R>20.01: 1.13e-06
        #    Max interpolation error with cubic splines: 1.007e-06
        #    Est Max overall error after cubic interpolation: 3.5e-6
        C4000_G1X1 => [
            [ -7.41093286, 12.00033716,  -1.99998713, -7.41220199, 0.99873007 ],
            [ -7.06092930, 11.30033656,  -1.99997409, -7.06273076, 0.99819693 ],
            [ -6.52103257, 10.22057052,  -1.99991837, -6.52412553, 0.99690237 ],
            [ -5.95359237, 9.08577753,   -1.99974263, -5.95905381, 0.99452430 ],
            [ -5.48624917, 8.15129037,   -1.99934374, -5.49497764, 0.99123608 ],
            [ -5.08451763, 7.34823152,   -1.99853802, -5.09758730, 0.98685379 ],
            [ -4.75825072, 6.69637013,   -1.99719776, -4.77640257, 0.98170709 ],
            [ -4.47453062, 6.12999822,   -1.99507516, -4.49869731, 0.97559693 ],
            [ -4.21723433, 5.61705797,   -1.99181015, -4.24857912, 0.96828463 ],
            [ -3.99554302, 5.17595164,   -1.98733838, -4.03477771, 0.96023054 ],
            [ -3.78773540, 4.76358043,   -1.98102784, -3.83617738, 0.95082275 ],
            [ -3.59403541, 4.38063446,   -1.97247596, -3.65301068, 0.94006698 ],
            [ -3.41108410, 4.02074986,   -1.96112135, -3.48211280, 0.92779403 ],
            [ -3.23529856, 3.67725730,   -1.94622694, -3.32022112, 0.91372944 ],
            [ -3.06108054, 3.33982452,   -1.92655001, -3.16242828, 0.89726656 ],
            [ -2.88283176, 2.99866480,   -1.90021801, -3.00420667, 0.87750785 ],
            [ -2.68263779, 2.62192466,   -1.86189992, -2.83106973, 0.85148896 ],
            [ -2.50013433, 2.28598192,   -1.81817213, -2.67811597, 0.82410161 ],
            [ -2.33241029, 1.98492463,   -1.77059921, -2.54221875, 0.79592061 ],
            [ -2.17456536, 1.70937119,   -1.71997766, -2.41884031, 0.76700817 ],
            [ -2.02341219, 1.45333596,   -1.66717341, -2.30511438, 0.73750967 ],
            [ -1.87628978, 1.21202305,   -1.61289552, -2.19880278, 0.70753607 ],
            [ -1.73118273, 0.98196538,   -1.55782147, -2.09832660, 0.67724097 ],
            [ -1.58666557, 0.76082299,   -1.50264990, -2.00265002, 0.64685037 ],
            [ -1.44111067, 0.54610891,   -1.44784035, -1.91071162, 0.61651649 ],
            [ -1.29320781, 0.33598802,   -1.39384663, -1.82176316, 0.58643982 ],
            [ -1.14178868, 0.12896389,   -1.34107657, -1.73522351, 0.55683843 ],
            [ -0.98570903, -0.07631002,  -1.28986901, -1.65059166, 0.52792498 ],
            [ -0.82387597, -0.28100756,  -1.24052011, -1.56745098, 0.49991035 ],
            [ -0.65492363, -0.48653854,  -1.19320220, -1.48530102, 0.47294847 ],
            [ -0.47757058, -0.69408886,  -1.14809659, -1.40374381, 0.44720152 ],
            [ -0.29046007, -0.90483550,  -1.10534910, -1.32239306, 0.42280980 ],
            [ -0.09198879, -1.12013375,  -1.06504337, -1.24080173, 0.39987288 ],
            [ 0.11951037,  -1.34130508,  -1.02725938, -1.15854405, 0.37847812 ],
            [ 0.34626966,  -1.57015111,  -0.99199126, -1.07502508, 0.35865725 ],
            [ 0.59057307,  -1.80839731,  -0.95926503, -0.98968962, 0.34044947 ],
            [ 0.83975315,  -2.04377485,  -0.93069651, -0.90686722, 0.32475522 ],
            [ 1.10212659,  -2.28450955,  -0.90505270, -0.82353484, 0.31087876 ],
            [ 1.37191325,  -2.52557783,  -0.88266478, -0.74131828, 0.29897897 ],
            [ 1.65908192,  -2.77607726,  -0.86252834, -0.65701468, 0.28849849 ],
            [ 1.96845856,  -3.04002734,  -0.84435921, -0.56923627, 0.27927796 ],
            [ 2.28471734,  -3.30454210,  -0.82888808, -0.48216008, 0.27165995 ],
            [ 2.62542446,  -3.58452337,  -0.81508941, -0.39076395, 0.26510446 ],
            [ 2.98968634,  -3.87915387,  -0.80300903, -0.29523896, 0.25961377 ],
            [ 3.37276715,  -4.18472730,  -0.79270550, -0.19667386, 0.25518099 ],
            [ 3.86371459,  -4.57124405,  -0.78233955, -0.07246948, 0.25105079 ],
            [ 4.39386250,  -4.98365344,  -0.77390065, 0.05977266,  0.24805137 ],
            [ 4.96778961,  -5.42578230,  -0.76717102, 0.20150465,  0.24602807 ],
            [ 5.59024529,  -5.90159004,  -0.76194127, 0.35422906,  0.24482864 ],
            [ 6.27015608,  -6.41820999,  -0.75798316, 0.52047410,  0.24429899 ],
            [ 7.01621957,  -6.98255872,  -0.75509476, 0.70270699,  0.24429534 ],
            [ 7.84551413,  -7.60784498,  -0.75306594, 0.90544066,  0.24468243 ],
            [ 8.77929700,  -8.31035451,  -0.75171350, 1.13421709,  0.24533814 ],
            [ 9.85361077,  -9.11742636,  -0.75086762, 1.39822767,  0.24615620 ],
            [ 11.12708577, -10.07328769, -0.75038261, 1.71228239,  0.24704576 ],
            [ 12.70454394, -11.25675902, -0.75013715, 2.10271748,  0.24792996 ],
            [ 14.80320296, -12.83090987, -0.75003477, 2.62395664,  0.24874359 ],
            [ 17.99028420, -15.22126750, -0.75000441, 3.41794523,  0.24942528 ],
            [ 24.75045336, -20.29140125, -0.75000006, 5.10610598,  0.24989334 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 9.28e-07
        #    Est Max FD Error (based on N/2 run) to R=30.01: 2.7e-7
        #    Est Max MOC error for R>30.01: -3.98e-07
        #    Max interpolation error with cubic splines: 1.0031e-06
        #    Est Max overall error after cubic interpolation: 3.7e-6
        S4000_G3 => [
            [ -4.16835563, 12.00033618,  -2.99997236, -4.16957049, 0.99817660 ],
            [ -3.43610249, 9.80363731,   -2.99980735, -3.43975058, 0.99451826 ],
            [ -3.14991916, 8.94517905,   -2.99952488, -3.15552837, 0.99156390 ],
            [ -2.87134947, 8.10967777,   -2.99889985, -2.87987971, 0.98715468 ],
            [ -2.61338139, 7.33620164,   -2.99761851, -2.62596476, 0.98102098 ],
            [ -2.38788325, 6.66047453,   -2.99532526, -2.40556936, 0.97327716 ],
            [ -2.20450354, 6.11147685,   -2.99192043, -2.22784161, 0.96467800 ],
            [ -2.03969060, 5.61875616,   -2.98682017, -2.06964653, 0.95458995 ],
            [ -1.89144853, 5.17648154,   -2.97958720, -1.92895822, 0.94306405 ],
            [ -1.75574657, 4.77277742,   -2.96963682, -1.80184048, 0.92997176 ],
            [ -1.62944970, 4.39851750,   -2.95624982, -1.68529547, 0.91513462 ],
            [ -1.51046063, 4.04774975,   -2.93860638, -1.57737225, 0.89838327 ],
            [ -1.39603209, 3.71274501,   -2.91552395, -1.47563187, 0.87933966 ],
            [ -1.28351026, 3.38630497,   -2.88535311, -1.37788746, 0.85745225 ],
            [ -1.16926361, 3.05884760,   -2.84540329, -1.28136148, 0.83172059 ],
            [ -1.04341460, 2.70418046,   -2.78861211, -1.17869726, 0.79906435 ],
            [ -0.92450744, 2.37645800,   -2.72137832, -1.08572421, 0.76407636 ],
            [ -0.81552201, 2.08376453,   -2.64804612, -1.00434923, 0.72876148 ],
            [ -0.71340964, 1.81727768,   -2.57006024, -0.93172624, 0.69331966 ],
            [ -0.61678689, 1.57279677,   -2.48952731, -0.86642084, 0.65824177 ],
            [ -0.52318389, 1.34360024,   -2.40712581, -0.80643297, 0.62342260 ],
            [ -0.42953886, 1.12214197,   -2.32240364, -0.74969443, 0.58837004 ],
            [ -0.33584945, 0.90854609,   -2.23736928, -0.69620099, 0.55366883 ],
            [ -0.24112680, 0.70063314,   -2.15296037, -0.64538372, 0.51949289 ],
            [ -0.14476789, 0.49719188,   -2.07025080, -0.59694680, 0.48611671 ],
            [ -0.04579055, 0.29630693,   -1.98979710, -0.55045536, 0.45364916 ],
            [ 0.05687832,  0.09606317,   -1.91197106, -0.50551586, 0.42216264 ],
            [ 0.16391716,  -0.10453838,  -1.83735334, -0.46197375, 0.39184450 ],
            [ 0.27633033,  -0.30701338,  -1.76619893, -0.41958552, 0.36277399 ],
            [ 0.39478963,  -0.51217481,  -1.69890567, -0.37827956, 0.33510586 ],
            [ 0.52063864,  -0.72190398,  -1.63543772, -0.33779260, 0.30882808 ],
            [ 0.65490969,  -0.93741200,  -1.57593834, -0.29802737, 0.28400816 ],
            [ 0.79884238,  -1.16015122,  -1.52043860, -0.25886677, 0.26067025 ],
            [ 0.95066734,  -1.38705878,  -1.46989107, -0.22095616, 0.23923234 ],
            [ 1.10824125,  -1.61503137,  -1.42478509, -0.18481646, 0.21992904 ],
            [ 1.27300445,  -1.84636508,  -1.38433280, -0.15005484, 0.20245151 ],
            [ 1.44608378,  -2.08273572,  -1.34797512, -0.11642154, 0.18658210 ],
            [ 1.62891982,  -2.32611776,  -1.31519547, -0.08366271, 0.17211632 ],
            [ 1.82243791,  -2.57769650,  -1.28565567, -0.05166336, 0.15892418 ],
            [ 2.02809263,  -2.83928052,  -1.25898902, -0.02025152, 0.14686039 ],
            [ 2.24716552,  -3.11238060,  -1.23491384, 0.01068114,  0.13581534 ],
            [ 2.48128167,  -3.39887548,  -1.21315738, 0.04126121,  0.12568148 ],
            [ 2.73178679,  -3.70024717,  -1.19351373, 0.07155025,  0.11638062 ],
            [ 3.00043476,  -4.01843054,  -1.17577244, 0.10163773,  0.10783130 ],
            [ 3.28908133,  -4.35543341,  -1.15974886, 0.13159793,  0.09996295 ],
            [ 3.60287878,  -4.71699989,  -1.14514201, 0.16178783,  0.09264502 ],
            [ 4.00465186,  -5.17388554,  -1.12974696, 0.19737513,  0.08475017 ],
            [ 4.44532198,  -5.66861562,  -1.11609146, 0.23308840,  0.07755881 ],
            [ 4.92956804,  -6.20603998,  -1.10397488, 0.26900820,  0.07099822 ],
            [ 5.46325486,  -6.79223981,  -1.09320522, 0.30524799,  0.06499751 ],
            [ 6.05283306,  -7.43384191,  -1.08361854, 0.34189737,  0.05949792 ],
            [ 6.70514070,  -8.13780791,  -1.07507607, 0.37901118,  0.05445194 ],
            [ 7.42860513,  -8.91272786,  -1.06744672, 0.41667466,  0.04981311 ],
            [ 8.23224476,  -9.76772907,  -1.06062171, 0.45493738,  0.04554463 ],
            [ 9.12647120,  -10.71333142, -1.05450435, 0.49385134,  0.04161332 ],
            [ 10.12289121, -11.76122597, -1.04901212, 0.53345372,  0.03799105 ],
            [ 11.23290838, -12.92281211, -1.04408049, 0.57371539,  0.03465779 ],
            [ 12.47112446, -14.21276472, -1.03964477, 0.61466877,  0.03158978 ],
            [ 13.85054017, -15.64401960, -1.03565885, 0.65623847,  0.02877264 ],
            [ 15.38915580, -17.23464213, -1.03207155, 0.69845350,  0.02618543 ],
            [ 17.10477145, -19.00241235, -1.02884356, 0.74127598,  0.02381287 ],
            [ 18.78517892, -20.72904933, -1.02625977, 0.77962161,  0.02188171 ],
            [ 20.81234205, -22.80676432, -1.02370270, 0.82194961,  0.01994086 ],
            [ 22.75382200, -24.79224390, -1.02168350, 0.85911346,  0.01838633 ],
            [ 24.96987039, -27.05415113, -1.01976412, 0.89815220,  0.01688980 ],
            [ 27.51372705, -29.64589352, -1.01794302, 0.93923327,  0.01545207 ],
            [ 30.45210200, -32.63437274, -1.01621872, 0.98254810,  0.01407393 ],
            [ 33.86946826, -36.10428234, -1.01458976, 1.02831690,  0.01275613 ],
            [ 37.01995573, -39.29873086, -1.01335428, 1.06686887,  0.01174584 ],
            [ 40.61537213, -42.93998199, -1.01217819, 1.10730369,  0.01077501 ],
            [ 44.74234286, -47.11483184, -1.01106079, 1.14979187,  0.00984400 ],
            [ 49.50953705, -51.93214485, -1.01000140, 1.19452853,  0.00895316 ],
            [ 55.05474030, -57.52993160, -1.00899934, 1.24173840,  0.00810283 ],
            [ 59.82871889, -62.34512067, -1.00828502, 1.27892411,  0.00749186 ],
            [ 65.23197685, -67.79124725, -1.00760229, 1.31777252,  0.00690399 ],
            [ 71.37893453, -73.98287390, -1.00695088, 1.35842570,  0.00633936 ],
            [ 78.41114519, -81.06171526, -1.00633051, 1.40104494,  0.00579809 ],
            [ 86.50538262, -89.20472939, -1.00574093, 1.44581438,  0.00528029 ],
            [ 95.88468434, -98.63516548, -1.00518185, 1.49294553,  0.00478610 ],
            [
                106.83365131, -109.63786927, -1.00465301, 1.54268305,
                0.00431562
            ],
            [
                115.18231060, -118.02393782, -1.00431712, 1.57742821,
                0.00401520
            ],
            [
                124.52515759, -127.40557279, -1.00399448, 1.61355372,
                0.00372540
            ],
            [
                135.02535605, -137.94604561, -1.00368501, 1.65116719,
                0.00344625
            ],
            [
                146.88087279, -149.84344396, -1.00338862, 1.69038944,
                0.00317778
            ],
            [
                160.33376814, -163.33996442, -1.00310524, 1.73135667,
                0.00292002
            ],
            [
                175.68251161, -178.73423094, -1.00283480, 1.77422327,
                0.00267300
            ],
            [
                193.29850260, -196.39781889, -1.00257722, 1.81916517,
                0.00243674
            ],
            [
                203.09808524, -206.22203963, -1.00245322, 1.84247624,
                0.00232266
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 1.01e-06
        #    Est Max FD Error (based on N/2 run) to R=27.93: TBD
        #    Est Max MOC error for R>27.93: 2.15e-06
        #    Max interpolation error with cubic splines: 1.0042e-06
        #    Est Max overall error after cubic interpolation: TBD
        'P4000_G1X1' => [
            [
                -13.67229811, 12.00011429, -0.99999927, -13.67399083,
                0.99915293
            ],
            [
                -13.23787270, 11.56569344, -0.99998840, -13.23997652,
                0.99894700
            ],
            [
                -11.89070608, 10.21856416, -0.99994935, -11.89483630,
                0.99793073
            ],
            [ -10.69805159, 9.02602597, -0.99983304, -10.70556197, 0.99623134 ],
            [ -9.77523350,  8.10346094, -0.99958030, -9.78717204,  0.99399761 ],
            [ -9.02640950,  7.35510436, -0.99911393, -9.04381287,  0.99123035 ],
            [ -8.39335764,  6.72283418, -0.99833562, -8.41730936,  0.98790083 ],
            [ -7.84293617,  6.17363182, -0.99712588, -7.87457656,  0.98397564 ],
            [ -7.35409562,  5.68659766, -0.99534246, -7.39463704,  0.97941470 ],
            [ -6.91250444,  5.24758107, -0.99281911, -6.96324975,  0.97417078 ],
            [ -6.50674215,  4.84539097, -0.98935351, -6.56914277,  0.96817054 ],
            [ -6.12762761,  4.47114467, -0.98469555, -6.20335355,  0.96131055 ],
            [ -5.76754717,  4.11762612, -0.97853242, -5.85857679,  0.95345083 ],
            [ -5.41866326,  3.77757545, -0.97043714, -5.52747203,  0.94437090 ],
            [ -5.07112294,  3.44208263, -0.95975014, -5.20107191,  0.93368273 ],
            [ -4.70813140,  3.09622809, -0.94520670, -4.86447123,  0.92056045 ],
            [ -4.31401351,  2.72754095, -0.92492559, -4.50486722,  0.90385715 ],
            [ -3.95908347,  2.40313310, -0.90241432, -4.18706826,  0.88655732 ],
            [ -3.63035332,  2.11039048, -0.87811461, -3.89852068,  0.86868648 ],
            [ -3.32401623,  1.84519704, -0.85290307, -3.63514898,  0.85058807 ],
            [ -3.03425001,  1.60172534, -0.82734266, -3.39129067,  0.83239680 ],
            [ -2.74918633,  1.36958893, -0.80121421, -3.15664978,  0.81373398 ],
            [ -2.46098762,  1.14253176, -0.77448904, -2.92491376,  0.79437720 ],
            [ -2.16540836,  0.91762309, -0.74744238, -2.69307443,  0.77432771 ],
            [ -1.86444560,  0.69669624, -0.72090635, -2.46309194,  0.75402969 ],
            [ -1.55421978,  0.47709700, -0.69513409, -2.23236542,  0.73354261 ],
            [ -1.23120318,  0.25661472, -0.67038941, -1.99876600,  0.71296909 ],
            [ -0.89130324,  0.03281930, -0.64688482, -1.75995699,  0.69240616 ],
            [ -0.52954663, -0.19710817,  -0.62478478, -1.51322216, 0.67194365 ],
            [ -0.14015803, -0.43628669,  -0.60423981, -1.25557901, 0.65169025 ],
            [ 0.28385129,  -0.68836978,  -0.58538299, -0.98355717, 0.63177353 ],
            [ 0.70559456,  -0.93186886,  -0.56983951, -0.72090629, 0.61412661 ],
            [ 1.13143735,  -1.17169152,  -0.55693322, -0.46280911, 0.59837903 ],
            [ 1.57049653,  -1.41375479,  -0.54609185, -0.20327953, 0.58415344 ],
            [ 2.02710423,  -1.66095396,  -0.53700801, 0.06044731,  0.57131868 ],
            [ 2.50268633,  -1.91448376,  -0.52947219, 0.32935505,  0.55984010 ],
            [ 3.00594035,  -2.17930004,  -0.52320272, 0.60842925,  0.54953227 ],
            [ 3.53964180,  -2.45709948,  -0.51805217, 0.89919841,  0.54038246 ],
            [ 4.22240985,  -2.80904488,  -0.51315716, 1.26479390,  0.53092148 ],
            [ 4.96844189,  -3.19038173,  -0.50937631, 1.65776080,  0.52292654 ],
            [ 5.79197662,  -3.60861148,  -0.50650320, 2.08555025,  0.51632455 ],
            [ 6.71500630,  -4.07506594,  -0.50434964, 2.55953982,  0.51102271 ],
            [ 7.77083423,  -4.60666807,  -0.50275756, 3.09677426,  0.50692479 ],
            [ 9.00204704,  -5.22490430,  -0.50161231, 3.71891429,  0.50394814 ],
            [ 10.43133177, -5.94124997,  -0.50084768, 4.43766799,  0.50201104 ],
            [ 12.10263413, -6.77788607,  -0.50038737, 5.27562595,  0.50089446 ],
            [ 14.08601655, -7.77007192,  -0.50014811, 6.26845003,  0.50033631 ],
            [ 16.60563114, -9.03009287,  -0.50004258, 7.52874056,  0.50009599 ],
            [ 20.00009424, -10.72739419, -0.50000784, 9.22612909,  0.50001762 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 8.1e-07
        #    Est Max FD Error (based on N/2 run) to R=10.00: 3.e-6
        #    Est Max MOC error for R>10.00: -2.09e-07
        #    Max interpolation error with cubic splines: 2.0058e-06
        #    Est Max overall error after cubic interpolation: 6.6e-6
        S4000_G1X1 => [
            [ -5.03706899, 12.00033733,  -2.99998070, -5.03808417, 0.99847647 ],
            [ -4.39797100, 10.08307984,  -2.99987272, -4.40062087, 0.99602008 ],
            [ -3.95145334, 8.74365006,   -2.99950190, -3.95663693, 0.99220561 ],
            [ -3.74252607, 8.11701386,   -2.99904148, -3.74962383, 0.98931853 ],
            [ -3.38670031, 7.05013257,   -2.99728215, -3.39883125, 0.98170789 ],
            [ -3.14843813, 6.33629149,   -2.99445585, -3.16581848, 0.97374587 ],
            [ -2.93771942, 5.70576010,   -2.98962864, -2.96162057, 0.96383087 ],
            [ -2.75185700, 5.15074201,   -2.98205204, -2.78352838, 0.95199703 ],
            [ -2.58838765, 4.66409639,   -2.97107378, -2.62896944, 0.93842619 ],
            [ -2.43745841, 4.21677796,   -2.95537013, -2.48848635, 0.92255302 ],
            [ -2.29371011, 3.79344459,   -2.93316467, -2.35717429, 0.90376370 ],
            [ -2.15621671, 3.39212345,   -2.90281069, -2.23437454, 0.88179709 ],
            [ -2.01626007, 2.98869496,   -2.86000942, -2.11279192, 0.85481154 ],
            [ -1.87675129, 2.59350669,   -2.80281715, -1.99571084, 0.82278217 ],
            [ -1.75269373, 2.24965821,   -2.73836846, -1.89563977, 0.78981567 ],
            [ -1.63916090, 1.94265549,   -2.66804890, -1.80785519, 0.75605195 ],
            [ -1.53271195, 1.66256208,   -2.59310001, -1.72918635, 0.72158962 ],
            [ -1.43104745, 1.40286785,   -2.51479768, -1.65759013, 0.68659678 ],
            [ -1.33218539, 1.15820932,   -2.43410641, -1.59145420, 0.65116611 ],
            [ -1.23477557, 0.92508596,   -2.35207849, -1.52975778, 0.61549905 ],
            [ -1.13777158, 0.70092089,   -2.26972489, -1.47178265, 0.57984223 ],
            [ -1.04013783, 0.48333145,   -2.18782460, -1.41690644, 0.54440273 ],
            [ -0.94116049, 0.27080014,   -2.10723060, -1.36476183, 0.50947237 ],
            [ -0.83978870, 0.06121883,   -2.02837164, -1.31486734, 0.47519973 ],
            [ -0.73543703, -0.14641194,  -1.95192254, -1.26703677, 0.44187332 ],
            [ -0.62716989, -0.35369749,  -1.87820242, -1.22096341, 0.40964415 ],
            [ -0.51408492, -0.56203570,  -1.80750442, -1.17641658, 0.37866366 ],
            [ -0.39535966, -0.77256420,  -1.74012752, -1.13324476, 0.34908980 ],
            [ -0.27018245, -0.98632022,  -1.67632437, -1.09133384, 0.32105825 ],
            [ -0.13746036, -1.20473236,  -1.61617331, -1.05051241, 0.29462175 ],
            [ 0.00394851,  -1.42919487,  -1.55972465, -1.01064263, 0.26981714 ],
            [ 0.15532561,  -1.66121618,  -1.50698535, -0.97159304, 0.24665172 ],
            [ 0.31804075,  -1.90233534,  -1.45793266, -0.93325468, 0.22511724 ],
            [ 0.49382629,  -2.15451700,  -1.41246135, -0.89548224, 0.20516208 ],
            [ 0.68438635,  -2.41956747,  -1.37049886, -0.85818916, 0.18674823 ],
            [ 0.89169140,  -2.69956889,  -1.33195729, -0.82128048, 0.16981641 ],
            [ 1.11823594,  -2.99719688,  -1.29665313, -0.78462160, 0.15427914 ],
            [ 1.36674919,  -3.31530431,  -1.26444298, -0.74810274, 0.14005512 ],
            [ 1.63018122,  -3.64455666,  -1.23615851, -0.71291156, 0.12750122 ],
            [ 1.92737919,  -4.00790712,  -1.20989104, -0.67681971, 0.11575468 ],
            [ 2.24193293,  -4.38478252,  -1.18711771, -0.64207659, 0.10547133 ],
            [ 2.58499322,  -4.78843241,  -1.16679911, -0.60753706, 0.09618248 ],
            [ 3.03239420,  -5.30552624,  -1.14561026, -0.56679096, 0.08633205 ],
            [ 3.52850957,  -5.86913831,  -1.12725049, -0.52620574, 0.07760601 ],
            [ 4.08008283,  -6.48633191,  -1.11135433, -0.48561865, 0.06985120 ],
            [ 4.69562640,  -7.16599964,  -1.09757845, -0.44483198, 0.06292818 ],
            [ 5.38482362,  -7.91815376,  -1.08562749, -0.40367952, 0.05672283 ],
            [ 6.15853325,  -8.75392654,  -1.07524685, -0.36203101, 0.05114234 ],
            [ 7.03099438,  -9.68792040,  -1.06619683, -0.31969215, 0.04609985 ],
            [ 8.01743082,  -10.73558760, -1.05828764, -0.27655369, 0.04153222 ],
            [ 9.13625598,  -11.91558053, -1.05135225, -0.23249291, 0.03738433 ],
            [ 10.40707657, -13.24762358, -1.04525922, -0.18746780, 0.03361618 ],
            [ 11.85229560, -14.75421180, -1.03989664, -0.14145130, 0.03019339 ],
            [ 13.49771422, -16.46122231, -1.03516939, -0.09442406, 0.02708565 ],
            [ 15.36993277, -18.39521963, -1.03100298, -0.04644974, 0.02427024 ],
            [ 17.36973555, -20.45340967, -1.02752902, -0.00039270, 0.02186324 ],
            [ 19.64264781, -22.78522416, -1.02442553, 0.04673444,  0.01966376 ],
            [ 22.35527078, -25.56001707, -1.02153429, 0.09715164,  0.01756995 ],
            [ 25.18113563, -28.44327031, -1.01917414, 0.14427498,  0.01582638 ],
            [ 28.53763309, -31.86027262, -1.01696857, 0.19453233,  0.01416701 ],
            [ 31.94055749, -35.31781686, -1.01519845, 0.24038263,  0.01281282 ],
            [ 35.94824707, -39.38295917, -1.01353716, 0.28907884,  0.01152245 ],
            [ 40.70977314, -44.20509075, -1.01198261, 0.34094085,  0.01029676 ],
            [ 45.39424780, -48.94274902, -1.01076719, 0.38684815,  0.00932537 ],
            [ 50.88788646, -54.49227868, -1.00962329, 0.43547184,  0.00839994 ],
            [ 57.38497789, -61.04826828, -1.00854976, 0.48711164,  0.00752090 ],
            [ 65.14142491, -68.86697363, -1.00754546, 0.54212050,  0.00668870 ],
            [ 72.47967140, -76.25772547, -1.00679112, 0.58883260,  0.00605692 ],
            [ 81.07089917, -84.90412894, -1.00607981, 0.63822062,  0.00545557 ],
            [ 91.21347887, -95.10484972, -1.00541096, 0.69057953,  0.00488485 ],
            [
                103.29969931, -107.25252382, -1.00478401, 0.74625647,
                0.00434492
            ],
            [
                113.94866948, -117.95000333, -1.00434097, 0.79043679,
                0.00396029
            ],
            [
                126.26889082, -130.32103121, -1.00392095, 0.83691728,
                0.00359315
            ],
            [
                140.62549977, -144.73097882, -1.00352373, 0.88593283,
                0.00324358
            ],
            [
                157.48865673, -161.65027968, -1.00314908, 0.93775651,
                0.00291162
            ],
            [
                177.47191612, -181.69280831, -1.00279675, 0.99270838,
                0.00259734
            ],
            [
                201.38798495, -205.67165022, -1.00246653, 1.05116687,
                0.00230080
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 1.2e-06
        #    Est Max FD Error (based on N/2 run) to R=100.04: 3.5e-7
        #    Est Max MOC error for R>100.04: 4.32e-07
        #    Max interpolation error with cubic splines: 1.003e-06
        #    Est Max overall error after cubic interpolation: 1.8e-6
        C4000_G3 => [
            [ -6.07231246, 12.00033674,  -1.99998157, -6.07383128, 0.99848004 ],
            [ -5.03193080, 9.91962969,   -1.99987756, -5.03623536, 0.99568646 ],
            [ -4.49136228, 8.83861646,   -1.99963547, -4.49876408, 0.99257224 ],
            [ -4.20517894, 8.26639483,   -1.99933524, -4.21504459, 0.99008909 ],
            [ -3.81509180, 7.48661339,   -1.99855029, -3.82969584, 0.98530031 ],
            [ -3.48322662, 6.82356512,   -1.99719048, -3.50362970, 0.97941863 ],
            [ -3.19793901, 6.25406887,   -1.99504466, -3.22515442, 0.97248477 ],
            [ -2.94691278, 5.75362837,   -1.99185119, -2.98200058, 0.96444683 ],
            [ -2.72074285, 5.30361118,   -1.98727590, -2.76487813, 0.95518747 ],
            [ -2.51472668, 4.89480767,   -1.98096584, -2.56914064, 0.94465657 ],
            [ -2.32348993, 4.51674023,   -1.97244847, -2.38959514, 0.93268668 ],
            [ -2.14381834, 4.16330333,   -1.96119634, -2.22319879, 0.91913728 ],
            [ -1.97206355, 3.82765761,   -1.94650155, -2.06661759, 0.90374812 ],
            [ -1.80389165, 3.50185003,   -1.92730921, -1.91608209, 0.88605075 ],
            [ -1.63460010, 3.17761962,   -1.90203791, -1.76779270, 0.86532371 ],
            [ -1.45380619, 2.83676038,   -1.86721043, -1.61360590, 0.83973272 ],
            [ -1.27253766, 2.50214265,   -1.82317010, -1.46399071, 0.81042205 ],
            [ -1.10745876, 2.20505890,   -1.77483373, -1.33262245, 0.78070447 ],
            [ -0.95339449, 1.93553079,   -1.72306979, -1.21462943, 0.75070566 ],
            [ -0.80701264, 1.68721958,   -1.66886455, -1.10692583, 0.72062342 ],
            [ -0.66694735, 1.45730779,   -1.61361894, -1.00806595, 0.69089097 ],
            [ -0.52788644, 1.23684972,   -1.55684584, -0.91406882, 0.66096951 ],
            [ -0.38878349, 1.02427731,   -1.49950470, -0.82420423, 0.63114865 ],
            [ -0.24874044, 0.81828911,   -1.44249648, -0.73788670, 0.60171483 ],
            [ -0.10647575, 0.61708883,   -1.38642579, -0.65435233, 0.57283936 ],
            [ 0.03941173,  0.41885504,   -1.33173074, -0.57285781, 0.54464818 ],
            [ 0.19034097,  0.22190453,   -1.27875705, -0.49274550, 0.51725602 ],
            [ 0.34745785,  0.02504552,   -1.22787616, -0.41358089, 0.49082212 ],
            [ 0.51204685,  -0.17298910,  -1.17934717, -0.33491591, 0.46547314 ],
            [ 0.68559308,  -0.37359274,  -1.13333026, -0.25626938, 0.44130333 ],
            [ 0.86982839,  -0.57831011,  -1.08991096, -0.17711911, 0.41838084 ],
            [ 1.06625310,  -0.78831021,  -1.04922644, -0.09710354, 0.39680935 ],
            [ 1.27703520,  -1.00537169,  -1.01127636, -0.01564134, 0.37662400 ],
            [ 1.50072285,  -1.22760497,  -0.97661102, 0.06648629,  0.35815434 ],
            [ 1.73237577,  -1.45018948,  -0.94590566, 0.14750863,  0.34179447 ],
            [ 1.97390142,  -1.67526772,  -0.91863720, 0.22826116,  0.32729221 ],
            [ 2.22701177,  -1.90463234,  -0.89440721, 0.30943091,  0.31445454 ],
            [ 2.49359786,  -2.14011661,  -0.87287858, 0.39170342,  0.30311576 ],
            [ 2.77522933,  -2.38318026,  -0.85380128, 0.47562325,  0.29315224 ],
            [ 3.07343614,  -2.63520156,  -0.83696423, 0.56170397,  0.28445747 ],
            [ 3.39007404,  -2.89779838,  -0.82216883, 0.65054036,  0.27692869 ],
            [ 3.72679245,  -3.17239015,  -0.80925291, 0.74266120,  0.27047942 ],
            [ 4.08556918,  -3.46065248,  -0.79805953, 0.83868524,  0.26502352 ],
            [ 4.46839506,  -3.76426216,  -0.78844798, 0.93923619,  0.26048067 ],
            [ 4.87667774,  -4.08444257,  -0.78029735, 1.04479465,  0.25677757 ],
            [ 5.40068689,  -4.49112209,  -0.77228877, 1.17839180,  0.25333965 ],
            [ 5.97021596,  -4.92905270,  -0.76592679, 1.32191118,  0.25083274 ],
            [ 6.59162190,  -5.40338080,  -0.76099619, 1.47720244,  0.24911984 ],
            [ 7.27223637,  -5.91998174,  -0.75729117, 1.64635938,  0.24806646 ],
            [ 8.02356876,  -6.48787411,  -0.75460164, 1.83250927,  0.24753892 ],
            [ 8.86270579,  -7.12023687,  -0.75272948, 2.04015037,  0.24741259 ],
            [ 9.81331307,  -7.83514367,  -0.75149553, 2.27540266,  0.24757460 ],
            [ 10.91503974, -8.66262007,  -0.75073604, 2.54834755,  0.24792559 ],
            [ 12.23252109, -9.65138187,  -0.75031110, 2.87528736,  0.24838071 ],
            [ 13.88478279, -10.89089091, -0.75010411, 3.28609788,  0.24886925 ],
            [ 16.13344579, -12.57750937, -0.75002317, 3.84628437,  0.24933399 ],
            [ 19.71933687, -15.26695922, -0.75000213, 4.74116557,  0.24972366 ],
            [ 26.02169453, -19.99373069, -0.75000002, 6.31587683,  0.24994228 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 9.97e-07
        #    Est Max FD Error (based on N/2 run) to R=50.01: 4.5e-7
        #    Est Max MOC error for R>50.01: -6.31e-07
        #    Max interpolation error with cubic splines: 1.0044e-06
        #    Est Max overall error after cubic interpolation: 2.1e-6
        C4000_G2 => [
            [ -6.37389118, 12.00033646,  -1.99998362, -6.37532307, 0.99856709 ],
            [ -5.85616358, 10.96489509,  -1.99995396, -5.85856775, 0.99759300 ],
            [ -5.45577563, 10.16414460,  -1.99991307, -5.45936570, 0.99640366 ],
            [ -4.98350595, 9.21967819,   -1.99976104, -4.98926912, 0.99422092 ],
            [ -4.48772355, 8.22831587,   -1.99935646, -4.49720198, 0.99047974 ],
            [ -4.09002277, 7.43330511,   -1.99857615, -4.10416037, 0.98577267 ],
            [ -3.75748580, 6.76890299,   -1.99723674, -3.77724909, 0.98006919 ],
            [ -3.47085460, 6.19670833,   -1.99511349, -3.49725021, 0.97332188 ],
            [ -3.21848479, 5.69356901,   -1.99194376, -3.25255907, 0.96548601 ],
            [ -2.99186700, 5.24263524,   -1.98741080, -3.03474320, 0.95648311 ],
            [ -2.78495511, 4.83202266,   -1.98113817, -2.83786034, 0.94621373 ],
            [ -2.59327050, 4.45303128,   -1.97268139, -2.65756674, 0.93455702 ],
            [ -2.41278242, 4.09794187,   -1.96147855, -2.49004878, 0.92132390 ],
            [ -2.23988458, 3.76001171,   -1.94680398, -2.33202107, 0.90624513 ],
            [ -2.07081696, 3.43241438,   -1.92765190, -2.18023188, 0.88890514 ],
            [ -1.90018183, 3.10555368,   -1.90236145, -2.03025018, 0.86851245 ],
            [ -1.71655872, 2.75932892,   -1.86720890, -1.87304984, 0.84307712 ],
            [ -1.53464514, 2.42350940,   -1.82331589, -1.72225249, 0.81422220 ],
            [ -1.36882783, 2.12505539,   -1.77521202, -1.58963468, 0.78488366 ],
            [ -1.21386784, 1.85388050,   -1.72374517, -1.47028598, 0.75516059 ],
            [ -1.06670339, 1.60411104,   -1.66998970, -1.36133429, 0.72529349 ],
            [ -0.92540317, 1.37198848,   -1.61511052, -1.26093983, 0.69559030 ],
            [ -0.78484352, 1.14892142,   -1.55868203, -1.16527621, 0.66555299 ],
            [ -0.64433036, 0.93390534,   -1.50178521, -1.07386852, 0.63554737 ],
            [ -0.50298811, 0.72564675,   -1.44530329, -0.98614478, 0.60586960 ],
            [ -0.35922465, 0.52189013,   -1.38969015, -0.90115807, 0.57663537 ],
            [ -0.21194708, 0.32125044,   -1.33547735, -0.81835669, 0.54804514 ],
            [ -0.05967426, 0.12194044,   -1.28295717, -0.73704708, 0.52021160 ],
            [ 0.09875334,  -0.07725895,  -1.23246947, -0.65679107, 0.49330379 ],
            [ 0.26446200,  -0.27743338,  -1.18429713, -0.57721854, 0.46748494 ],
            [ 0.43912789,  -0.48021942,  -1.13852641, -0.49775577, 0.44282956 ],
            [ 0.62400726,  -0.68663756,  -1.09534393, -0.41808667, 0.41947594 ],
            [ 0.82101151,  -0.89834103,  -1.05477069, -0.33766177, 0.39747756 ],
            [ 1.03199171,  -1.11678345,  -1.01686366, -0.25602311, 0.37690769 ],
            [ 1.25889681,  -1.34342104,  -0.98167227, -0.17272064, 0.35783322 ],
            [ 1.49355707,  -1.57002604,  -0.95048548, -0.09078117, 0.34098475 ],
            [ 1.73797403,  -1.79886225,  -0.92276334, -0.00930974, 0.32608974 ],
            [ 1.99383124,  -2.03171689,  -0.89810804, 0.07239226,  0.31294526 ],
            [ 2.26306459,  -2.27048111,  -0.87617583, 0.15504235,  0.30137267 ],
            [ 2.54711727,  -2.51651723,  -0.85672383, 0.23916359,  0.29124360 ],
            [ 2.84778285,  -2.77143979,  -0.83952294, 0.32536167,  0.28243429 ],
            [ 3.16659128,  -3.03659837,  -0.82439291, 0.41415114,  0.27484423 ],
            [ 3.50532506,  -3.31353304,  -0.81116314, 0.50611306,  0.26837594 ],
            [ 3.86562943,  -3.60365921,  -0.79968673, 0.60179122,  0.26294189 ],
            [ 4.24929976,  -3.90851368,  -0.78982280, 0.70177559,  0.25845565 ],
            [ 4.73994514,  -4.29351045,  -0.77998281, 0.82748780,  0.25422679 ],
            [ 5.27097299,  -4.70548007,  -0.77200997, 0.96159670,  0.25107454 ],
            [ 5.83914435,  -5.14224333,  -0.76576144, 1.10357820,  0.24887926 ],
            [ 6.46625131,  -5.62082664,  -0.76085299, 1.25915605,  0.24743799 ],
            [ 7.15199484,  -6.14123123,  -0.75717249, 1.42852624,  0.24664589 ],
            [ 7.90887965,  -6.71323924,  -0.75450525, 1.61507229,  0.24636425 ],
            [ 8.75258578,  -7.34897441,  -0.75265696, 1.82295233,  0.24646582 ],
            [ 9.70837459,  -8.06771900,  -0.75144365, 2.05868495,  0.24683506 ],
            [ 10.81549095, -8.89919653,  -0.75070209, 2.33225365,  0.24737120 ],
            [ 12.14056722, -9.89362324,  -0.75029113, 2.66045747,  0.24798894 ],
            [ 13.80782700, -11.14435942, -0.75009407, 3.07446926,  0.24861766 ],
            [ 16.08908999, -12.85541489, -0.75001944, 3.64234926,  0.24919865 ],
            [ 19.79838395, -15.63741100, -0.75000141, 4.56772241,  0.24967909 ],
            [ 25.30332165, -19.76611620, -0.75000004, 5.94299590,  0.24991868 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 8.85e-07
        #    Est Max FD Error (based on N/2 run) to R=20.00: 3.5e-7
        #    Est Max MOC error for R>20.00: -5.9e-07
        #    Max interpolation error with cubic splines: 1.003e-06
        #    Est Max overall error after cubic interpolation: 1.9e-6
        S4000_G2 => [
            [ -4.35786149, 12.00033595,  -2.99997543, -4.35900684, 0.99828100 ],
            [ -3.89566869, 10.61377889,  -2.99991317, -3.89796101, 0.99655768 ],
            [ -3.67218593, 9.94335716,   -2.99984785, -3.67539259, 0.99518255 ],
            [ -3.45004675, 9.27699085,   -2.99968424, -3.45452418, 0.99326951 ],
            [ -3.13296099, 8.32590197,   -2.99917746, -3.14017454, 0.98914351 ],
            [ -2.85819812, 7.50196417,   -2.99812652, -2.86910918, 0.98355397 ],
            [ -2.62978889, 6.81735052,   -2.99628967, -2.64518825, 0.97675070 ],
            [ -2.43407352, 6.23119093,   -2.99334590, -2.45477183, 0.96869798 ],
            [ -2.26205371, 5.71662615,   -2.98890242, -2.28890747, 0.95932355 ],
            [ -2.10821337, 5.25726751,   -2.98250958, -2.14211966, 0.94856842 ],
            [ -1.96803280, 4.83975967,   -2.97361058, -2.00997762, 0.93630804 ],
            [ -1.83837697, 4.45494820,   -2.96155134, -1.88945188, 0.92240500 ],
            [ -1.71678264, 4.09575936,   -2.94555404, -1.77821972, 0.90668764 ],
            [ -1.60048074, 3.75434597,   -2.92453966, -1.67378043, 0.88882045 ],
            [ -1.48694792, 3.42380200,   -2.89705241, -1.57400301, 0.86833291 ],
            [ -1.37285858, 3.09525712,   -2.86079442, -1.47626828, 0.84438319 ],
            [ -1.25126425, 2.75031379,   -2.81078355, -1.37534742, 0.81487730 ],
            [ -1.12738859, 2.40597625,   -2.74625410, -1.27648873, 0.78049760 ],
            [ -1.01497542, 2.10113636,   -2.67537828, -1.19067624, 0.74569319 ],
            [ -0.91002073, 1.82425933,   -2.59924736, -1.11423922, 0.71048430 ],
            [ -0.81081947, 1.57029753,   -2.51982287, -1.04549024, 0.67530767 ],
            [ -0.71564895, 1.33431910,   -2.43858431, -0.98287543, 0.64039760 ],
            [ -0.62126551, 1.10808112,   -2.35512250, -0.92408934, 0.60525234 ],
            [ -0.52676679, 0.88951576,   -2.27067628, -0.86855563, 0.57014167 ],
            [ -0.43177817, 0.67782708,   -2.18676502, -0.81605153, 0.53549084 ],
            [ -0.33517766, 0.47060701,   -2.10403498, -0.76598005, 0.50141144 ],
            [ -0.23620033, 0.26638981,   -2.02326858, -0.71801425, 0.46811906 ],
            [ -0.13401922, 0.06369390,   -1.94504317, -0.67185186, 0.43578531 ],
            [ -0.02787611, -0.13871313,  -1.86986445, -0.62727414, 0.40458467 ],
            [ 0.08330213,  -0.34253843,  -1.79793353, -0.58398602, 0.37458888 ],
            [ 0.20019392,  -0.54863981,  -1.72965446, -0.54190090, 0.34596985 ],
            [ 0.32402755,  -0.75875054,  -1.66505548, -0.50077524, 0.31875144 ],
            [ 0.45569876,  -0.97390761,  -1.60433786, -0.46053269, 0.29303388 ],
            [ 0.59646380,  -1.19565077,  -1.54750894, -0.42102442, 0.26883601 ],
            [ 0.74756226,  -1.42538150,  -1.49459208, -0.38215548, 0.24618154 ],
            [ 0.90453763,  -1.65618234,  -1.44718961, -0.34515123, 0.22577362 ],
            [ 1.06826847,  -1.88956325,  -1.40467591, -0.30972907, 0.20736234 ],
            [ 1.23986271,  -2.12723409,  -1.36645795, -0.27561118, 0.19070617 ],
            [ 1.42055569,  -2.37095066,  -1.33201935, -0.24255144, 0.17559224 ],
            [ 1.61167716,  -2.62247762,  -1.30092260, -0.21033940, 0.16183875 ],
            [ 1.81437641,  -2.88325015,  -1.27283365, -0.17883714, 0.14930717 ],
            [ 2.03010801,  -3.15502572,  -1.24743179, -0.14789288, 0.13786325 ],
            [ 2.26015350,  -3.43927725,  -1.22446796, -0.11741215, 0.12740374 ],
            [ 2.50591404,  -3.73758145,  -1.20371447, -0.08730766, 0.11783445 ],
            [ 2.76918320,  -4.05194170,  -1.18494662, -0.05747029, 0.10906180 ],
            [ 3.04999549,  -4.38225245,  -1.16807221, -0.02799773, 0.10105523 ],
            [ 3.41339325,  -4.80333063,  -1.15000259, 0.00708964,  0.09232345 ],
            [ 3.81018693,  -5.25635582,  -1.13398633, 0.04210479,  0.08441327 ],
            [ 4.24427253,  -5.74541388,  -1.11979349, 0.07714080,  0.07723458 ],
            [ 4.72053091,  -6.27562027,  -1.10720340, 0.11232012,  0.07070094 ],
            [ 5.24422717,  -6.85242746,  -1.09602745, 0.14773677,  0.06474100 ],
            [ 5.82161248,  -7.48228634,  -1.08609243, 0.18349461,  0.05929002 ],
            [ 6.45892518,  -8.17154988,  -1.07725670, 0.21964089,  0.05429877 ],
            [ 7.16479308,  -8.92906865,  -1.06937474, 0.25629988,  0.04971348 ],
            [ 7.94803488,  -9.76379001,  -1.06233184, 0.29353390,  0.04549510 ],
            [ 8.81826254,  -10.68541959, -1.05602965, 0.33138212,  0.04161139 ],
            [ 9.78708318,  -11.70568636, -1.05037741, 0.36990722,  0.03803136 ],
            [ 10.86570062, -12.83580849, -1.04530600, 0.40909396,  0.03473413 ],
            [ 12.06851679, -14.09027905, -1.04074661, 0.44898707,  0.03169555 ],
            [ 13.40893252, -15.48246752, -1.03664886, 0.48953886,  0.02890030 ],
            [ 14.90394814, -17.02941908, -1.03296178, 0.53076169,  0.02632959 ],
            [ 16.47505678, -18.64976029, -1.02981686, 0.57031720,  0.02409261 ],
            [ 18.21549999, -20.43953885, -1.02697197, 0.61040043,  0.02203180 ],
            [ 20.23433552, -22.51002262, -1.02428889, 0.65281789,  0.02005385 ],
            [ 22.23219253, -24.55416578, -1.02211533, 0.69120892,  0.01842565 ],
            [ 24.52764430, -26.89794129, -1.02005610, 0.73165320,  0.01686054 ],
            [ 27.18146119, -29.60231157, -1.01810918, 0.77434273,  0.01535948 ],
            [ 29.72166361, -32.18650340, -1.01657113, 0.81179079,  0.01415822 ],
            [ 32.61955386, -35.13023031, -1.01510862, 0.85109680,  0.01300264 ],
            [ 35.94426038, -38.50278668, -1.01372056, 0.89242663,  0.01189327 ],
            [ 39.78233961, -42.39091859, -1.01240589, 0.93596978,  0.01083062 ],
            [ 44.24330879, -46.90436130, -1.01116356, 0.98194418,  0.00981517 ],
            [ 48.35347250, -51.05841516, -1.01022110, 1.02064239,  0.00903712 ],
            [ 53.04145444, -55.79214555, -1.00932378, 1.06120538,  0.00828983 ],
            [ 58.41927354, -61.21773992, -1.00847107, 1.10380155,  0.00757352 ],
            [ 64.62732315, -67.47578436, -1.00766247, 1.14862355,  0.00688842 ],
            [ 71.84343005, -74.74432746, -1.00689749, 1.19589308,  0.00623473 ],
            [ 78.05178432, -80.99376417, -1.00635208, 1.23310448,  0.00576521 ],
            [ 85.07411331, -88.05881659, -1.00583071, 1.27196020,  0.00531356 ],
            [ 93.05761375, -96.08682136, -1.00533318, 1.31260011,  0.00487987 ],
            [
                102.18418409, -105.25983566, -1.00485928, 1.35518285,
                0.00446421
            ],
            [
                112.68070419, -115.80492043, -1.00440882, 1.39988945,
                0.00406665
            ],
            [
                124.83304697, -128.00815641, -1.00398157, 1.44692766,
                0.00368728
            ],
            [
                139.00545319, -142.23402448, -1.00357735, 1.49653752,
                0.00332616
            ],
            [
                149.80257687, -153.06835257, -1.00332056, 1.53117556,
                0.00309558
            ],
            [
                161.87620580, -165.18054509, -1.00307386, 1.56717406,
                0.00287316
            ],
            [
                175.43425313, -178.77862053, -1.00283717, 1.60463834,
                0.00265892
            ],
            [
                190.72862580, -194.11460397, -1.00261044, 1.64368649,
                0.00245288
            ],
            [
                208.06687002, -211.49617470, -1.00239362, 1.68445159,
                0.00225504
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 7.64e-07
        #    Est Max FD Error (based on N/2 run) to R=100.03: 5.8e-7
        #    Est Max MOC error for R>100.03: -5.12e-06
        #    Max interpolation error with cubic splines: 1.0016e-06
        #    Est Max overall error after cubic interpolation: 6.7e-6
        P4000_G2 => [
            [
                -11.51913736, 12.00011367, -0.99999181, -11.52104722,
                0.99904417
            ],
            [
                -11.07687438, 11.55785545, -0.99998206, -11.07925746,
                0.99880706
            ],
            [
                -10.54548702, 11.02647687, -0.99997798, -10.54859646,
                0.99844290
            ],
            [ -9.11569000, 9.59674978,  -0.99990807, -9.12205573, 0.99680733 ],
            [ -8.06455966, 8.54579045,  -0.99973722, -8.07534946, 0.99457759 ],
            [ -7.23198561, 7.71355759,  -0.99939651, -7.24838895, 0.99173664 ],
            [ -6.54138358, 7.02355543,  -0.99879827, -6.56462265, 0.98826125 ],
            [ -5.94998280, 6.43312248,  -0.99783526, -5.98132582, 0.98412143 ],
            [ -5.43135949, 5.91596739,  -0.99637951, -5.47213637, 0.97928017 ],
            [ -4.96759917, 5.45433665,  -0.99427910, -5.01922940, 0.97368818 ],
            [ -4.54593781, 5.03566276,  -0.99135441, -4.60996367, 0.96728269 ],
            [ -4.15674639, 4.65056030,  -0.98739226, -4.23487823, 0.95998332 ],
            [ -3.79173948, 4.29106300,  -0.98212692, -3.88594750, 0.95167129 ],
            [ -3.44375129, 3.95043481,  -0.97522164, -3.55638025, 0.94218311 ],
            [ -3.10516819, 3.62169536,  -0.96620864, -3.23917213, 0.93125910 ],
            [ -2.76613723, 3.29604264,  -0.95435080, -2.92556383, 0.91843879 ],
            [ -2.40762430, 2.95668374,  -0.93811555, -2.59905164, 0.90264714 ],
            [ -2.03498810, 2.61095214,  -0.91667626, -2.26614353, 0.88369286 ],
            [ -1.69707679, 2.30508047,  -0.89303030, -1.97075128, 0.86431489 ],
            [ -1.38248764, 2.02805884,  -0.86762515, -1.70191084, 0.84458936 ],
            [ -1.09155973, 1.77936586,  -0.84168806, -1.45899772, 0.82516744 ],
            [ -0.81417847, 1.54951828,  -0.81537945, -1.23277115, 0.80589667 ],
            [ -0.54020266, 1.32978731,  -0.78855776, -1.01463182, 0.78646326 ],
            [ -0.26150453, 1.11384370,  -0.76113742, -0.79821248, 0.76662670 ],
            [ 0.02642372,  0.89870843,  -0.73338274, -0.58040043, 0.74640808 ],
            [ 0.32081006,  0.68684154,  -0.70625367, -0.36363914, 0.72635153 ],
            [ 0.62557474,  0.47565113,  -0.68001807, -0.14532086, 0.70652852 ],
            [ 0.94439724,  0.26291018,  -0.65494814, 0.07679365,  0.68704179 ],
            [ 1.28173770,  0.04604926,  -0.63125103, 0.30529787,  0.66797096 ],
            [ 1.64280348,  -0.17777722, -0.60910255, 0.54306829,  0.64939381 ],
            [ 2.03391479,  -0.41189006, -0.58865247, 0.79346179,  0.63138636 ],
            [ 2.44030233,  -0.64739491, -0.57092452, 1.04662271,  0.61488328 ],
            [ 2.84755888,  -0.87683045, -0.55629657, 1.29402483,  0.60041361 ],
            [ 3.26355320,  -1.10562400, -0.54411275, 1.54105425,  0.58754705 ],
            [ 3.69342675,  -1.33726220, -0.53396698, 1.79109088,  0.57603876 ],
            [ 4.14155180,  -1.57458824, -0.52556095, 2.04685685,  0.56571969 ],
            [ 4.61221692,  -1.82025764, -0.51865658, 2.31088411,  0.55646159 ],
            [ 5.20579233,  -2.12608033, -0.51214336, 2.63819757,  0.54672008 ],
            [ 5.84404861,  -2.45128253, -0.50719098, 2.98433394,  0.53821601 ],
            [ 6.53645609,  -2.80111569, -0.50354208, 3.35434445,  0.53083320 ],
            [ 7.29537815,  -3.18221065, -0.50096650, 3.75469036,  0.52447059 ],
            [ 8.13706128,  -3.60308010, -0.49925979, 4.19373959,  0.51904310 ],
            [ 9.09195150,  -4.07927014, -0.49823529, 4.68705951,  0.51444332 ],
            [ 10.05059549, -4.55664602, -0.49777969, 5.17851921,  0.51105275 ],
            [ 11.32796187, -5.19237414, -0.49765758, 5.82915617,  0.50787570 ],
            [ 12.87428094, -5.96203024, -0.49784834, 6.61241604,  0.50538156 ],
            [ 14.82518095, -6.93367057, -0.49825703, 7.59634339,  0.50348055 ],
            [ 17.73088557, -8.38233612, -0.49883850, 9.05683579,  0.50195524 ],
            [ 20.00011492, -9.51472254, -0.49918041, 10.19508048, 0.50129585 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 1.47e-06
        #    Est Max FD Error (based on N/2 run) to R=70.01: 2.9e-7
        #    Est Max MOC error for R>70.01: -3.19e-07
        #    Max interpolation error with cubic splines: 3.5313e-06
        #    Est Max overall error after cubic interpolation: 1.6e-6
        C4000_G7 => [
            [ -5.58930901, 12.00033753,  -1.99997850, -5.59094962, 0.99835806 ],
            [ -4.56281105, 9.94740342,   -1.99986470, -4.56739702, 0.99540385 ],
            [ -4.09905073, 9.01999018,   -1.99964644, -4.10635209, 0.99267334 ],
            [ -3.69418417, 8.21047846,   -1.99920122, -3.70514836, 0.98898030 ],
            [ -3.30700812, 7.43659221,   -1.99827266, -3.32319434, 0.98369753 ],
            [ -2.98857878, 6.80050177,   -1.99674209, -3.01089296, 0.97747533 ],
            [ -2.70956122, 6.24368205,   -1.99432711, -2.73914398, 0.97006885 ],
            [ -2.47510586, 5.77646470,   -1.99097429, -2.51261780, 0.96196521 ],
            [ -2.25511186, 5.33896229,   -1.98608278, -2.30200849, 0.95235412 ],
            [ -2.05245465, 4.93710875,   -1.97932871, -2.11008303, 0.94135509 ],
            [ -1.86457686, 4.56603641,   -1.97029159, -1.93435470, 0.92891761 ],
            [ -1.68722381, 4.21759965,   -1.95837496, -1.77082159, 0.91482342 ],
            [ -1.51717856, 3.88584286,   -1.94284691, -1.61658407, 0.89881833 ],
            [ -1.35024159, 3.56312273,   -1.92260359, -1.46803466, 0.88042107 ],
            [ -1.18157980, 3.24100453,   -1.89595598, -1.32131562, 0.85886502 ],
            [ -0.99885242, 2.89781714,   -1.85878657, -1.16677979, 0.83194419 ],
            [ -0.82179207, 2.57256171,   -1.81361774, -1.02204864, 0.80230905 ],
            [ -0.66000551, 2.28303186,   -1.76428208, -0.89463187, 0.77239070 ],
            [ -0.50858638, 2.01979856,   -1.71163081, -0.77993201, 0.74231228 ],
            [ -0.36459899, 1.77724882,   -1.65672699, -0.67519464, 0.71231856 ],
            [ -0.22612345, 1.55168263,   -1.60070951, -0.57860096, 0.68269979 ],
            [ -0.08833738, 1.33507631,   -1.54321035, -0.48658117, 0.65299498 ],
            [ 0.04948860,  1.12637468,   -1.48531773, -0.39861638, 0.62355140 ],
            [ 0.18832559,  0.92416065,   -1.42790002, -0.31406339, 0.59462388 ],
            [ 0.32984974,  0.72610590,   -1.37140253, -0.23193056, 0.56628786 ],
            [ 0.47518072,  0.53083714,   -1.31639019, -0.15165524, 0.53871481 ],
            [ 0.62557477,  0.33690515,   -1.26327612, -0.07266670, 0.51203325 ],
            [ 0.78249544,  0.14272890,   -1.21233196, 0.00563629,  0.48633004 ],
            [ 0.94730837,  -0.05301157,  -1.16381785, 0.08372780,  0.46171094 ],
            [ 1.12145315,  -0.25161120,  -1.11793326, 0.16205394,  0.43827106 ],
            [ 1.30659403,  -0.45450741,  -1.07479704, 0.24109985,  0.41607966 ],
            [ 1.50475360,  -0.66339758,  -1.03445176, 0.32143286,  0.39517748 ],
            [ 1.71804420,  -0.87993510,  -0.99694702, 0.40358359,  0.37561584 ],
            [ 1.93932969,  -1.09676156,  -0.96362720, 0.48471844,  0.35813206 ],
            [ 2.16901373,  -1.31461493,  -0.93415119, 0.56514366,  0.34258764 ],
            [ 2.40895208,  -1.53553213,  -0.90802385, 0.64563901,  0.32875570 ],
            [ 2.66087402,  -1.76128210,  -0.88486043, 0.72686726,  0.31646052 ],
            [ 2.92660642,  -1.99361122,  -0.86434103, 0.80946939,  0.30555512 ],
            [ 3.20757058,  -2.23383786,  -0.84622974, 0.89392598,  0.29593263 ],
            [ 3.50554899,  -2.48354741,  -0.83030116, 0.98080869,  0.28748797 ],
            [ 3.82241082,  -2.74435636,  -0.81636061, 1.07069706,  0.28012904 ],
            [ 4.15991414,  -3.01776516,  -0.80424347, 1.16413093,  0.27377664 ],
            [ 4.52020218,  -3.30557281,  -0.79378999, 1.26175454,  0.26835130 ],
            [ 4.98118758,  -3.66897690,  -0.78333004, 1.38416800,  0.26300621 ],
            [ 5.48138024,  -4.05855306,  -0.77479410, 1.51459916,  0.25874693 ],
            [ 6.02577017,  -4.47838348,  -0.76795783, 1.65450694,  0.25544866 ],
            [ 6.62012665,  -4.93314207,  -0.76260732, 1.80555390,  0.25298784 ],
            [ 7.27400122,  -5.43037019,  -0.75852066, 1.97035848,  0.25123512 ],
            [ 7.99471444,  -5.97588533,  -0.75551275, 2.15097024,  0.25007569 ],
            [ 8.79777188,  -6.58168005,  -0.75337816, 2.35148587,  0.24938670 ],
            [ 9.70525533,  -7.26464239,  -0.75193219, 2.57762406,  0.24905688 ],
            [ 10.75002654, -8.04970666,  -0.75101001, 2.83777541,  0.24898683 ],
            [ 11.98513106, -8.97690465,  -0.75046690, 3.14535188,  0.24908961 ],
            [ 13.50480189, -10.11711559, -0.75018042, 3.52403724,  0.24929179 ],
            [ 15.49366360, -11.60896758, -0.75005242, 4.02009415,  0.24953278 ],
            [ 18.39533922, -13.78529548, -0.75000901, 4.74452424,  0.24976273 ],
            [ 23.76067375, -17.80931110, -0.75000042, 6.08514812,  0.24993649 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 6.44e-07
        #    Est Max FD Error (based on N/2 run) to R=100.03: 4.8e-7
        #    Est Max MOC error for R>100.03: -1.89e-06
        #    Max interpolation error with cubic splines: 1.0031e-06
        #    Est Max overall error after cubic interpolation: 3.4e-6
        P4000_G3 => [
            [
                -10.85856780, 12.00011389, -0.99999078, -10.86059363,
                0.99898607
            ],
            [ -9.22656143, 10.36814314, -0.99995575, -9.23114849, 0.99770133 ],
            [ -8.00994637, 9.15163344,  -0.99985034, -8.01838995, 0.99576112 ],
            [ -7.06428735, 8.20620996,  -0.99961511, -7.07786812, 0.99316656 ],
            [ -6.29978821, 7.44215238,  -0.99917434, -6.31974985, 0.98992921 ],
            [ -5.65665947, 6.79976793,  -0.99843261, -5.68428478, 0.98602199 ],
            [ -5.09952008, 6.24379532,  -0.99727308, -5.13615715, 0.98140485 ],
            [ -4.60664878, 5.75265676,  -0.99555824, -4.65371713, 0.97603604 ],
            [ -4.16281829, 5.31129891,  -0.99312584, -4.22183991, 0.96986072 ],
            [ -3.75660814, 4.90851661,  -0.98978209, -3.82925526, 0.96280416 ],
            [ -3.37906112, 4.53562538,  -0.98529244, -3.46721943, 0.95476638 ],
            [ -3.02278750, 4.18559068,  -0.97936813, -3.12864255, 0.94561457 ],
            [ -2.68062151, 3.85174324,  -0.97162875, -2.80682567, 0.93514885 ],
            [ -2.34450164, 3.52678008,  -0.96152570, -2.49448438, 0.92304357 ],
            [ -2.00251817, 3.20014848,  -0.94811908, -2.18121494, 0.90866991 ],
            [ -1.62586686, 2.84647777,  -0.92905658, -1.84234193, 0.89027293 ],
            [ -1.27074627, 2.52041528,  -0.90653992, -1.52963937, 0.87043816 ],
            [ -0.94510912, 2.22911093,  -0.88199229, -1.24942865, 0.85026525 ],
            [ -0.63918001, 1.96321179,  -0.85586235, -0.99239928, 0.82984540 ],
            [ -0.35667505, 1.72509158,  -0.82963290, -0.76074380, 0.81004228 ],
            [ -0.08318702, 1.50182289,  -0.80296480, -0.54189275, 0.79034370 ],
            [ 0.19011699,  1.28608243,  -0.77575983, -0.32860187, 0.77048970 ],
            [ 0.47129143,  1.07188396,  -0.74792285, -0.11481399, 0.75024850 ],
            [ 0.75965447,  0.86023139,  -0.72023403, 0.09859824,  0.73003345 ],
            [ 1.05617617,  0.65070751,  -0.69327440, 0.31209289,  0.71013184 ],
            [ 1.36436809,  0.44110109,  -0.66734173, 0.52791383,  0.69065120 ],
            [ 1.68826138,  0.22902209,  -0.64267478, 0.74849575,  0.67168076 ],
            [ 2.03257119,  0.01182705,  -0.61946712, 0.97654414,  0.65329796 ],
            [ 2.40290590,  -0.21348137, -0.59788170, 1.21513681,  0.63557513 ],
            [ 2.80635563,  -0.45057349, -0.57804962, 1.46805175,  0.61857316 ],
            [ 3.20936198,  -0.68011897, -0.56164380, 1.71429257,  0.60379455 ],
            [ 3.61751425,  -0.90647975, -0.54801497, 1.95802299,  0.59083208 ],
            [ 4.03675225,  -1.13376208, -0.53665578, 2.20325622,  0.57935679 ],
            [ 4.47175396,  -1.36507637, -0.52721404, 2.45299923,  0.56914889 ],
            [ 4.92708022,  -1.60328301, -0.51941684, 2.71001794,  0.56004110 ],
            [ 5.49875764,  -1.89799427, -0.51201436, 3.02736361,  0.55050515 ],
            [ 6.11127126,  -2.20976795, -0.50633198, 3.36192010,  0.54219258 ],
            [ 6.77262538,  -2.54314278, -0.50210502, 3.71802365,  0.53496774 ],
            [ 7.49486686,  -2.90461137, -0.49908866, 4.10204569,  0.52869710 ],
            [ 8.28669592,  -3.29893822, -0.49708493, 4.51846255,  0.52331419 ],
            [ 9.17520489,  -3.74000539, -0.49588565, 4.98127588,  0.51867369 ],
            [ 10.18064394, -4.23825962, -0.49534042, 5.50069131,  0.51473453 ],
            [ 11.42393886, -4.85404845, -0.49532078, 6.13834462,  0.51121982 ],
            [ 12.82282075, -5.54718757, -0.49571438, 6.85143507,  0.50846390 ],
            [ 14.53199938, -6.39502457, -0.49639942, 7.71840275,  0.50617404 ],
            [ 16.83527551, -7.53944918, -0.49731707, 8.88178038,  0.50417485 ],
            [ 17.98128672, -8.10961270, -0.49772102, 9.45915028,  0.50346612 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 5.12e-07
        #    Est Max FD Error (based on N/2 run) to R=500.13: 3.8e-7
        #    Est Max MOC error for R>500.13: 1.15e-06
        #    Max interpolation error with cubic splines: 1.0008e-06
        #    Est Max overall error after cubic interpolation: 2.5e-6
        P4000_G7 => [
            [ -9.78736350, 12.00011442, -0.99998925, -9.78955181, 0.99890466 ],
            [ -8.89410425, 11.10686911, -0.99998041, -8.89752675, 0.99828587 ],
            [ -7.81095183, 10.02376131, -0.99993242, -7.81684121, 0.99704688 ],
            [ -6.69286624, 8.90581540,  -0.99979281, -6.70318865, 0.99481345 ],
            [ -5.83230086, 8.04553296,  -0.99950985, -5.84821500, 0.99198434 ],
            [ -5.09265802, 7.30642634,  -0.99897553, -5.11576814, 0.98832581 ],
            [ -4.49341967, 6.70802829,  -0.99813814, -4.52471354, 0.98414377 ],
            [ -3.95755370, 6.17347976,  -0.99683016, -3.99862736, 0.97912121 ],
            [ -3.48080540, 5.69866237,  -0.99492198, -3.53316286, 0.97329989 ],
            [ -3.04929110, 5.26987746,  -0.99224249, -3.11455773, 0.96661621 ],
            [ -2.65248604, 4.87683067,  -0.98858665, -2.73246135, 0.95898510 ],
            [ -2.28265532, 4.51207068,  -0.98371502, -2.37935242, 0.95031022 ],
            [ -1.93175603, 4.16795130,  -0.97730848, -2.04757324, 0.94042080 ],
            [ -1.59264317, 3.83788295,  -0.96894767, -1.73053412, 0.92908224 ],
            [ -1.25708644, 3.51449777,  -0.95801902, -1.42092648, 0.91591502 ],
            [ -0.91114038, 3.18549856,  -0.94338996, -1.10673436, 0.90011921 ],
            [ -0.52766250, 2.82756182,  -0.92254739, -0.76535246, 0.87984526 ],
            [ -0.18306002, 2.51352477,  -0.89934580, -0.46565065, 0.85919075 ],
            [ 0.13525405,  2.23115770,  -0.87423387, -0.19544469, 0.83827222 ],
            [ 0.43210363,  1.97547178,  -0.84802869, 0.05033654,  0.81748402 ],
            [ 0.71038979,  1.74312121,  -0.82159313, 0.27502515,  0.79722589 ],
            [ 0.98224156,  1.52341083,  -0.79468491, 0.48901822,  0.77708334 ],
            [ 1.25620966,  1.30945786,  -0.76719422, 0.69913175,  0.75680466 ],
            [ 1.54013465,  1.09564471,  -0.73904657, 0.91106704,  0.73618951 ],
            [ 1.82937976,  0.88590808,  -0.71141356, 1.12105778,  0.71595735 ],
            [ 2.12758658,  0.67780604,  -0.68459721, 1.33158537,  0.69620746 ],
            [ 2.43828700,  0.46915869,  -0.65888327, 1.54488180,  0.67705165 ],
            [ 2.76582425,  0.25742365,  -0.63448052, 1.76356539,  0.65856737 ],
            [ 3.11486331,  0.04005296,  -0.61158684, 1.99027906,  0.64084348 ],
            [ 3.49144137,  -0.18614858, -0.59034252, 2.22835370,  0.62394035 ],
            [ 3.89653187,  -0.42128270, -0.57115095, 2.47782800,  0.60815021 ],
            [ 4.29955137,  -0.64818508, -0.55537211, 2.72013493,  0.59465233 ],
            [ 4.70952181,  -0.87308444, -0.54222551, 2.96145294,  0.58290315 ],
            [ 5.13207889,  -1.09980320, -0.53125138, 3.20552497,  0.57259189 ],
            [ 5.57173661,  -1.33128883, -0.52212756, 3.45521566,  0.56350482 ],
            [ 6.03236189,  -1.56999086, -0.51461074, 3.71287946,  0.55548453 ],
            [ 6.49243206,  -1.80534933, -0.50877829, 3.96684633,  0.54874219 ],
            [ 7.11920825,  -2.12229679, -0.50293053, 4.30833278,  0.54120052 ],
            [ 7.78551206,  -2.45589430, -0.49869412, 4.66672757,  0.53481348 ],
            [ 8.50877071,  -2.81542405, -0.49573677, 5.05147205,  0.52932270 ],
            [ 9.30188036,  -3.20776995, -0.49384071, 5.46933011,  0.52459119 ],
            [ 10.18502143, -3.64338142, -0.49281507, 5.93072926,  0.52048764 ],
            [ 11.17564929, -4.13136727, -0.49250247, 6.44450215,  0.51693530 ],
            [ 12.30484525, -4.68759785, -0.49275426, 7.02639521,  0.51383765 ],
            [ 13.62605149, -5.33905525, -0.49344496, 7.70338202,  0.51109311 ],
            [ 15.23865091, -6.13560729, -0.49447647, 8.52545507,  0.50859589 ],
            [ 17.40237648, -7.20699013, -0.49580627, 9.62315531,  0.50618675 ],
            [ 19.89180096, -8.44287765, -0.49704928, 10.88072484, 0.50426325 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 8.97e-07
        #    Est Max FD Error (based on N/2 run) to R=16.41: 6.9e-7
        #    Est Max MOC error for R>16.41: 4.69e-07
        #    Est Max MOC error for R>16.41 based on N/2: 6.9e-07:
        #    Max interpolation error with cubic splines: 5.0491e-07
        #    Est Max overall error after cubic interpolation: 1.7e-6
        'P4000_G1X4_A' => [
            [
                -12.37483994, 12.00011384, -0.99999283, -12.37662634,
                0.99910601
            ],
            [
                -11.25518684, 10.88047738, -0.99997280, -11.25831576,
                0.99843313
            ],
            [ -9.69371861, 9.31910369,   -0.99988037, -9.70056161, 0.99656723 ],
            [ -8.70640850, 8.33199501,   -0.99967906, -8.71764267, 0.99435329 ],
            [ -7.90854370, 7.53452171,   -0.99928828, -7.92532794, 0.99154383 ],
            [ -7.24324582, 6.86989576,   -0.99861855, -7.26672342, 0.98814079 ],
            [ -6.66975950, 6.29747752,   -0.99755680, -6.70113793, 0.98410593 ],
            [ -6.16457202, 5.79389200,   -0.99597080, -6.20511513, 0.97940598 ],
            [ -5.71066465, 5.34229102,   -0.99370098, -5.76173689, 0.97398715 ],
            [ -5.29639300, 4.93123655,   -0.99056133, -5.35947942, 0.96778899 ],
            [ -4.91212805, 4.55136435,   -0.98632358, -4.98890535, 0.96072014 ],
            [ -4.55003224, 4.19518255,   -0.98070594, -4.64244757, 0.95265774 ],
            [ -4.20271585, 3.85578273,   -0.97333934, -4.31313237, 0.94341889 ],
            [ -3.86189344, 3.52561495,   -0.96369178, -3.99337078, 0.93270314 ],
            [ -3.51569670, 3.19411208,   -0.95086494, -3.67262383, 0.91994113 ],
            [ -3.13613430, 2.83650820,   -0.93266165, -3.32647265, 0.90358297 ],
            [ -2.77202311, 2.50077696,   -0.91069991, -3.00069287, 0.88547510 ],
            [ -2.43903827, 2.20142015,   -0.88671848, -2.70888444, 0.86689923 ],
            [ -2.12657579, 1.92827852,   -0.86114564, -2.44094406, 0.84789080 ],
            [ -1.83815315, 1.68357234,   -0.83542749, -2.19906217, 0.82922749 ],
            [ -1.55905369, 1.45403757,   -0.80923977, -1.97023612, 0.81041987 ],
            [ -1.28015141, 1.23206521,   -0.78247625, -1.74688328, 0.79118993 ],
            [ -0.99370211, 1.01185922,   -0.75508398, -1.52309546, 0.77131045 ],
            [ -0.70042178, 0.79443240,   -0.72782357, -1.29984989, 0.75115304 ],
            [ -0.39953227, 0.57947435,   -0.70127424, -1.07688535, 0.73099498 ],
            [ -0.08721783, 0.36450757,   -0.67569504, -0.85174644, 0.71091313 ],
            [ 0.24072525,  0.14698891,   -0.65130500, -0.62191025, 0.69098328 ],
            [ 0.58862612,  -0.07551622,  -0.62832090, -0.38498446, 0.67130902 ],
            [ 0.96186376,  -0.30593164,  -0.60690968, -0.13809153, 0.65198096 ],
            [ 1.36750306,  -0.54799972,  -0.58719108, 0.12246940,  0.63307132 ],
            [ 1.77917716,  -0.78619840,  -0.57055821, 0.37950838,  0.61602812 ],
            [ 2.19431779,  -1.02010224,  -0.55676748, 0.63203484,  0.60088126 ],
            [ 2.61967069,  -1.25439247,  -0.54526468, 0.88466602,  0.58729387 ],
            [ 3.06010253,  -1.49235526,  -0.53568096, 1.14057146,  0.57506569 ],
            [ 3.52033712,  -1.73699203,  -0.52773050, 1.40263693,  0.56404877 ],
            [ 4.00230354,  -1.98970486,  -0.52121848, 1.67204914,  0.55418724 ],
            [ 4.51548905,  -2.25575738,  -0.51589155, 1.95410460,  0.54530291 ],
            [ 5.06010592,  -2.53550536,  -0.51163796, 2.24887262,  0.53741837 ],
            [ 5.64469201,  -2.83356870,  -0.50828344, 2.56093046,  0.53043361 ],
            [ 6.27524823,  -3.15320773,  -0.50570059, 2.89339944,  0.52431322 ],
            [ 6.95758035,  -3.49756269,  -0.50377179, 3.24928177,  0.51902964 ],
            [ 7.71141091,  -3.87675173,  -0.50236421, 3.63875166,  0.51447832 ],
            [ 8.54522850,  -4.29519268,  -0.50139241, 4.06606459,  0.51066335 ],
            [ 9.48481951,  -4.76596918,  -0.50075832, 4.54431966,  0.50751823 ],
            [ 10.56570495, -5.30700217,  -0.50037962, 5.09143571,  0.50499281 ],
            [ 11.83255611, -5.94076541,  -0.50018340, 5.72986345,  0.50305308 ],
            [ 13.38925146, -6.71931999,  -0.50009863, 6.51174974,  0.50163028 ],
            [ 15.63786212, -7.84379498,  -0.50005745, 7.63842183,  0.50062166 ],
            [ 22.41672584, -11.23338380, -0.50000408, 11.02913458, 0.50002414 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 5.05e-08
        #    Est Max FD Error (based on N/2 run) to R=10.00: 2.05e-7
        #    Est Max MOC error for R>10.00: 1.95e-07
        #    Max interpolation error with cubic splines: 5.0741e-07
        #    Est Max overall error after cubic interpolation: 9.e-7
        S16000_G1X4 => [
            [ -4.61791146, 12.00009666,  -2.99997850, -4.61898292, 0.99839195 ],
            [ -4.54937966, 11.79450298,  -2.99997802, -4.55056720, 0.99821765 ],
            [ -4.34378429, 11.17772357,  -2.99994874, -4.34540115, 0.99757278 ],
            [ -4.18482860, 10.70086569,  -2.99991727, -4.18688126, 0.99691793 ],
            [ -3.66552111, 9.14303362,   -2.99966246, -3.66999950, 0.99326808 ],
            [ -3.38405925, 8.29879992,   -2.99920191, -3.39089772, 0.98970975 ],
            [ -3.11799464, 7.50093133,   -2.99822707, -3.12820282, 0.98461796 ],
            [ -2.88697980, 6.80847721,   -2.99646079, -2.90144221, 0.97817326 ],
            [ -2.68916703, 6.21599335,   -2.99361329, -2.70866537, 0.97052533 ],
            [ -2.51638181, 5.69908055,   -2.98932415, -2.54170503, 0.96166061 ],
            [ -2.36078553, 5.23440142,   -2.98308939, -2.39284074, 0.95140127 ],
            [ -2.21955191, 4.81366049,   -2.97441234, -2.25926581, 0.93972688 ],
            [ -2.08875921, 4.42535389,   -2.96261006, -2.13719631, 0.92644813 ],
            [ -1.96586110, 4.06216744,   -2.94688375, -2.02423581, 0.91137496 ],
            [ -1.84833804, 3.71699497,   -2.92618516, -1.91810895, 0.89420219 ],
            [ -1.73344651, 3.38228778,   -2.89902611, -1.81647862, 0.87442835 ],
            [ -1.61771512, 3.04876910,   -2.86306272, -1.71659139, 0.85117985 ],
            [ -1.49353859, 2.69621866,   -2.81303708, -1.61264908, 0.82221330 ],
            [ -1.36837758, 2.34798711,   -2.74917470, -1.51179426, 0.78866036 ],
            [ -1.25454986, 2.03894144,   -2.67899260, -1.42393730, 0.75446170 ],
            [ -1.14831875, 1.75826764,   -2.60375550, -1.34561359, 0.71971136 ],
            [ -1.04742472, 1.49948598,   -2.52497292, -1.27475339, 0.68464662 ],
            [ -0.95044999, 1.25850817,   -2.44428727, -1.21005094, 0.64960235 ],
            [ -0.85414183, 1.02708506,   -2.36129012, -1.14919558, 0.61409523 ],
            [ -0.75834290, 0.80487364,   -2.27783721, -1.09206366, 0.57868597 ],
            [ -0.66207275, 0.58959362,   -2.19486016, -1.03804883, 0.54359218 ],
            [ -0.56437015, 0.37917007,   -2.11308917, -0.98663853, 0.50900260 ],
            [ -0.46435393, 0.17186024,   -2.03315274, -0.93743922, 0.47511105 ],
            [ -0.36122030, -0.03378208,  -1.95560809, -0.89015880, 0.44211937 ],
            [ -0.25415724, -0.23910216,  -1.88089892, -0.84455428, 0.41020915 ],
            [ -0.14234997, -0.44534047,  -1.80938442, -0.80042968, 0.37954429 ],
            [ -0.02488290, -0.65381544,  -1.74130075, -0.75759539, 0.35024590 ],
            [ 0.09919100,  -0.86578969,  -1.67683173, -0.71589758, 0.32241597 ],
            [ 0.23087954,  -1.08252789,  -1.61610430, -0.67520493, 0.29613084 ],
            [ 0.37137682,  -1.30549620,  -1.55915095, -0.63537311, 0.27142210 ],
            [ 0.52187801,  -1.53605627,  -1.50600721, -0.59630307, 0.24831746 ],
            [ 0.68084946,  -1.77151489,  -1.45746712, -0.55854950, 0.22717020 ],
            [ 0.84616789,  -2.00877769,  -1.41399703, -0.52260179, 0.20818905 ],
            [ 1.01911914,  -2.24986538,  -1.37492528, -0.48811142, 0.19108373 ],
            [ 1.20107401,  -2.49675139,  -1.33969268, -0.45478618, 0.17560992 ],
            [ 1.39318111,  -2.75098028,  -1.30788422, -0.42243185, 0.16158470 ],
            [ 1.59676084,  -3.01423378,  -1.27912709, -0.39086689, 0.14884275 ],
            [ 1.81303565,  -3.28798874,  -1.25312325, -0.35996191, 0.13725137 ],
            [ 2.04344475,  -3.57393496,  -1.22959564, -0.32958674, 0.12668713 ],
            [ 2.28940411,  -3.87367549,  -1.20831290, -0.29964341, 0.11704723 ],
            [ 2.55226652,  -4.18869650,  -1.18907861, -0.27006333, 0.10824517 ],
            [ 2.83402426,  -4.52120775,  -1.17167983, -0.24072947, 0.10018730 ],
            [ 3.13661932,  -4.87330289,  -1.15594369, -0.21156062, 0.09279869 ],
            [ 3.46191882,  -5.24695077,  -1.14172460, -0.18250523, 0.08601778 ],
            [ 3.81218248,  -5.64453954,  -1.12887837, -0.15349676, 0.07978429 ],
            [ 4.18986168,  -6.06863720,  -1.11727426, -0.12447622, 0.07404490 ],
            [ 4.59800016,  -6.52243518,  -1.10678396, -0.09536504, 0.06874768 ],
            [ 5.03943358,  -7.00885135,  -1.09730372, -0.06612615, 0.06385309 ],
            [ 5.51799110,  -7.53185953,  -1.08872438, -0.03668213, 0.05931883 ],
            [ 6.03749541,  -8.09537750,  -1.08095551, -0.00698706, 0.05511194 ],
            [ 6.60236377,  -8.70392553,  -1.07391232, 0.02301076,  0.05120198 ],
            [ 7.21760882,  -9.36261795,  -1.06751762, 0.05336343,  0.04756179 ],
            [ 7.88883945,  -10.07715759, -1.06170241, 0.08411940,  0.04416766 ],
            [ 8.62206160,  -10.85362119, -1.05640720, 0.11531232,  0.04100004 ],
            [ 9.42387905,  -11.69867276, -1.05157929, 0.14696969,  0.03804192 ],
            [ 10.30149406, -12.61956231, -1.04717241, 0.17911119,  0.03527858 ],
            [ 11.26290792, -13.62433411, -1.04314539, 0.21175442,  0.03269676 ],
            [ 12.31612128, -14.72099198, -1.03946456, 0.24488787,  0.03028639 ],
            [ 13.47033446, -15.91875434, -1.03609808, 0.27851184,  0.02803702 ],
            [ 14.73534758, -17.22742539, -1.03301825, 0.31261648,  0.02593940 ],
            [ 16.11816068, -18.65389454, -1.03020696, 0.34710125,  0.02398965 ],
            [ 17.62757389, -20.20690959, -1.02764282, 0.38190916,  0.02218048 ],
            [ 19.24509083, -21.86722366, -1.02534153, 0.41641641,  0.02053043 ],
            [ 21.01788460, -23.68301225, -1.02322562, 0.45141067,  0.01899012 ],
            [ 22.95942935, -25.66770860, -1.02128220, 0.48684987,  0.01755482 ],
            [ 25.08379632, -27.83533729, -1.01949920, 0.52268620,  0.01621976 ],
            [ 27.40549258, -30.20035241, -1.01786530, 0.55886550,  0.01498024 ],
            [ 29.93923198, -32.77740781, -1.01636987, 0.59532664,  0.01383159 ],
            [ 32.69962331, -35.58104472, -1.01500298, 0.63200090,  0.01276919 ],
            [ 35.70075976, -38.62528044, -1.01375526, 0.66881130,  0.01178851 ],
            [ 38.95569285, -41.92308169, -1.01261800, 0.70567195,  0.01088512 ],
            [ 42.47577527, -45.48570676, -1.01158302, 0.74248740,  0.01005468 ],
            [ 46.47219808, -49.52639091, -1.01059681, 0.78102862,  0.00925567 ],
            [ 50.80188930, -53.89996955, -1.00970233, 0.81947489,  0.00852414 ],
            [ 55.47514690, -58.61662138, -1.00889255, 0.85770596,  0.00785600 ],
            [ 60.49723667, -63.68148044, -1.00816090, 0.89558914,  0.00724727 ],
            [ 66.20384951, -69.43261488, -1.00746342, 0.93522978,  0.00666228 ],
            [ 72.33738091, -75.60994189, -1.00683571, 0.97442031,  0.00613173 ],
            [ 78.89175415, -82.20721910, -1.00627218, 1.01299424,  0.00565194 ],
            [ 85.84880552, -89.20610214, -1.00576763, 1.05077014,  0.00521944 ],
            [ 93.73338442, -97.13420180, -1.00528586, 1.09023863,  0.00480377 ],
            [
                102.07903034, -105.52211959, -1.00485656, 1.12873151,
                0.00443105
            ],
            [
                111.55130442, -115.03839828, -1.00444674, 1.16896015,
                0.00407310
            ],
            [
                121.53864613, -125.06828559, -1.00408349, 1.20800215,
                0.00375398
            ],
            [
                132.88531077, -136.45925880, -1.00373675, 1.24880967,
                0.00344767
            ],
            [
                144.78595141, -148.40249937, -1.00343121, 1.28818010,
                0.00317634
            ],
            [
                158.31050503, -161.97142765, -1.00313952, 1.32932739,
                0.00291598
            ],
            [
                172.39530335, -176.09859656, -1.00288420, 1.36874142,
                0.00268699
            ],
            [
                188.39308393, -192.14050991, -1.00264034, 1.40991938,
                0.00246727
            ],
            [
                204.89778348, -208.68699076, -1.00242852, 1.44901601,
                0.00227558
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 1.96e-07
        #    Est Max FD Error (based on N/2 run) to R=30.01: 1.5e-7
        #    Est Max MOC error for R>30.01: 2.71e-07
        #    Max interpolation error with cubic splines: 5.0651e-07
        #    Est Max overall error after cubic interpolation: 9.e-7
        C8000_G1X4 => [
            [ -6.77636571, 12.00017690,  -1.99998567, -6.77770516, 0.99865965 ],
            [ -6.50370625, 11.45486317,  -1.99996835, -6.50546592, 0.99823879 ],
            [ -6.29028163, 11.02802023,  -1.99995807, -6.29246041, 0.99781889 ],
            [ -5.43602330, 9.31958777,   -1.99979533, -5.44114995, 0.99486072 ],
            [ -5.01074173, 8.46916254,   -1.99951558, -5.01859575, 0.99211697 ],
            [ -4.55791427, 7.56386204,   -1.99881071, -4.57029241, 0.98755240 ],
            [ -4.21101529, 6.87065893,   -1.99762458, -4.22856630, 0.98231530 ],
            [ -3.91584883, 6.28128025,   -1.99572587, -3.93948619, 0.97613271 ],
            [ -3.65539465, 5.76182930,   -1.99283751, -3.68615078, 0.96887919 ],
            [ -3.42248637, 5.29813042,   -1.98866403, -3.46142523, 0.96052123 ],
            [ -3.21074138, 4.87761509,   -1.98284815, -3.25901382, 0.95097462 ],
            [ -3.01489686, 4.49001171,   -1.97494739, -3.07379990, 0.94010306 ],
            [ -2.83114566, 4.12802473,   -1.96443845, -2.90215495, 0.92775055 ],
            [ -2.65545667, 3.78404721,   -1.95061274, -2.74036323, 0.91364318 ],
            [ -2.48413471, 3.45134112,   -1.93252350, -2.58519044, 0.89739642 ],
            [ -2.31169299, 3.12006589,   -1.90859633, -2.43205008, 0.87825462 ],
            [ -2.12790388, 2.77219498,   -1.87555754, -2.27276517, 0.85449359 ],
            [ -1.94056947, 2.42468982,   -1.83286839, -2.11524802, 0.82654389 ],
            [ -1.77036726, 2.11662036,   -1.78588681, -1.97696096, 0.79793614 ],
            [ -1.61181047, 1.83737150,   -1.73549797, -1.85272509, 0.76877259 ],
            [ -1.46095888, 1.57950990,   -1.68252672, -1.73896893, 0.73914547 ],
            [ -1.31630990, 1.34002365,   -1.62829485, -1.63418728, 0.70945775 ],
            [ -1.17295569, 1.11058180,   -1.57252927, -1.53463990, 0.67929246 ],
            [ -1.03052332, 0.89059695,   -1.51643335, -1.44003749, 0.64909559 ],
            [ -0.88740115, 0.67757137,   -1.46057945, -1.34929566, 0.61902189 ],
            [ -0.74228137, 0.46963261,   -1.40551126, -1.26163407, 0.58926369 ],
            [ -0.59382523, 0.26500940,   -1.35164297, -1.17634392, 0.55999075 ],
            [ -0.44095314, 0.06241942,   -1.29938902, -1.09294329, 0.53141266 ],
            [ -0.28232345, -0.13965254,  -1.24899990, -1.01087208, 0.50368029 ],
            [ -0.11676717, -0.34237589,  -1.20073411, -0.92972765, 0.47697014 ],
            [ 0.05713022,  -0.54711289,  -1.15474329, -0.84904291, 0.45141250 ],
            [ 0.24065747,  -0.75496778,  -1.11119454, -0.76846431, 0.42715178 ],
            [ 0.43547783,  -0.96737063,  -1.07015512, -0.68752106, 0.40427675 ],
            [ 0.64338721,  -1.18577628,  -1.03167808, -0.60574502, 0.38286436 ],
            [ 0.86637272,  -1.41172852,  -0.99579907, -0.52264478, 0.36297585 ],
            [ 1.10320327,  -1.64357842,  -0.96298010, -0.43887909, 0.34489916 ],
            [ 1.34881729,  -1.87642863,  -0.93385164, -0.35617460, 0.32899638 ],
            [ 1.60498180,  -2.11224281,  -0.90796322, -0.27373949, 0.31502139 ],
            [ 1.87383724,  -2.35316962,  -0.88492002, -0.19074343, 0.30275518 ],
            [ 2.15681636,  -2.60060535,  -0.86446113, -0.10663327, 0.29204882 ],
            [ 2.45576708,  -2.85624668,  -0.84633789, -0.02076048, 0.28275859 ],
            [ 2.77198800,  -3.12127340,  -0.83037043, 0.06734652,  0.27477554 ],
            [ 3.10770741,  -3.39761642,  -0.81635776, 0.15840985,  0.26797956 ],
            [ 3.46408023,  -3.68629841,  -0.80416838, 0.25285420,  0.26228419 ],
            [ 3.84335450,  -3.98923131,  -0.79364232, 0.35140183,  0.25758832 ],
            [ 4.24768376,  -4.30823576,  -0.78464198, 0.45474929,  0.25380098 ],
            [ 4.67899229,  -4.64495401,  -0.77704344, 0.56354152,  0.25083605 ],
            [ 5.13979349,  -5.00149624,  -0.77071662, 0.67858061,  0.24860378 ],
            [ 5.63298697,  -5.38027137,  -0.76553306, 0.80076941,  0.24701464 ],
            [ 6.16225924,  -5.78428818,  -0.76136327, 0.93120759,  0.24597900 ],
            [ 6.74069101,  -6.22367017,  -0.75803918, 1.07330025,  0.24540377 ],
            [ 7.33591360,  -6.67410445,  -0.75560084, 1.21929903,  0.24522003 ],
            [ 8.00761631,  -7.18096375,  -0.75370228, 1.38403517,  0.24532547 ],
            [ 8.74622006,  -7.73710721,  -0.75233118, 1.56534675,  0.24565868 ],
            [ 9.56915578,  -8.35580331,  -0.75138053, 1.76770638,  0.24615279 ],
            [ 10.50104685, -9.05568849,  -0.75075534, 1.99737190,  0.24674831 ],
            [ 11.58131105, -9.86647133,  -0.75037124, 2.26428106,  0.24739413 ],
            [ 12.87476188, -10.83688138, -0.75015654, 2.58471043,  0.24804616 ],
            [ 14.50041048, -12.05627232, -0.75005208, 2.98847902,  0.24866655 ],
            [ 16.71546856, -13.71762531, -0.75001142, 3.53995703,  0.24922172 ],
            [ 20.25635767, -16.37330726, -0.75000099, 4.42334270,  0.24967625 ],
        ],
    };
}

1;
