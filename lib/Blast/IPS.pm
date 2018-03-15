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
our $VERSION = 0.1.1;

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
# Allow lookup on slope dYdX or dZdX, with quadratic interpolation
# Work is needed to handle interpolations beyond the ranges of the tables for some
# spherical symmetry variables.
# Question: how accurate would it be to interpolate on gamma (or gamma-1) between tables?
# Are there any theoretical results to guide this (such as low gamma theories)?
# An interesting way to test is to have the driver get a builtin table, truncate at
# both ends, reinstall it and look at the errors at both missing ends.

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

sub get_info {
    my ($self)     = @_;

    # Return some information about this case
    my $rinfo      = {};
    my $table_name = $self->{_table_name};
    my $item       = $rtables_info->{$table_name};
    my ( $symmetry, $gamma, $err_est ) = @{$item};
    $rinfo->{table_name}      = $table_name;
    $rinfo->{alpha}           = $self->{_alpha};
    $rinfo->{symmetry}        = $self->{_symmetry};
    $rinfo->{gamma}           = $self->{_gamma};
    $rinfo->{estimated_error} = $err_est;
    return ($rinfo);
}

sub get_table_index {

    # returns a list of references of the form
    #   [NAME, symmetry, gamma, error]
    # one per built-in table

    my @table_list;
    foreach my $key(sort keys %{$rtables_info}) {
        my $item = $rtables_info->{$key};
        my( $symmetry, $gamma, $err_est ) = @{$item};
	push @table_list, [$key, $symmetry, $gamma, $err_est];
    }
    return ( \@table_list );
}


sub check_tables {

    # Should be called at program installation to check the tables

    my @info_keys  = keys %{$rtables_info};
    my @table_keys = keys %{$rtables};

    # Check that the keys of the two table hashes are the same
    my @missing_info_keys =
      grep { !exists $rtables_info->{$_} } @table_keys;
    my @missing_table_keys =
      grep { !exists $rtables->{$_} } @info_keys;
    my $error = @missing_info_keys || @missing_table_keys;
    if ($error) {
        local $" = ')(';
        $error =<<EOM;
------------------------------------------------------------------------
Program error detected checking hash keys
Have info but not table for: (@missing_table_keys)
Have table but no info for: (@missing_info_keys)
------------------------------------------------------------------------
EOM
	croak $error;
    }

    return $error if ($error); 

    foreach my $key (@info_keys) {
	my $rtable=$rtables->{$key};
	#print STDERR "Checking table $key\n";
	my $error=_check_table($rtable);
        if ($error) {
            return "Table $key:\n" . $error;
        }
    } 
    return; 
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

    my ($self) = @_;

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

    # Distant region
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
    return $self->_long_range_sphere( $Q, $icol ) if ( $symmetry == 2 );
    return $self->_long_range_non_sphere( $Q, $icol );
}

sub _long_range_sphere {

    # FIXME: This routine works but needs review because the definition of Z
    # has changed.
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
    my $Z_i = log( exp($Z_zero) + $kA * sqrt($term) );

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

    # For planar and cylindrical waves at long range it is a good approximation
    # that dY_dX and dZ_dX become constant.  We can use the slopes at the end
    # point of the table. But in fact, at long range:
    # dY_dX is about -0.5 for sym=0 and -0.75 for sym=1
    # dZ_dX is about +0.5 for sym=0 and +0.25 for sym=1

    if ( $icol == 0 ) {

        # Given X
        $X_i = $Q;
        $Y_i = $Y_e + $dY_dX_i * ( $X_i - $X_e );
        $Z_i = $Z_e + $dZ_dX_i * ( $X_i - $X_e );
    }
    elsif ( $icol == 1 ) {

        # Given Y
        $Y_i = $Q;
        $X_i = $X_e + ( $Y_i - $Y_e ) / $dY_dX_i;
        $Z_i = $Z_e + $dZ_dX_i * ( $X_i - $X_e );
    }
    elsif ( $icol == 3 ) {

        # Given Z
        $Z_i = $Q;
        $X_i = $X_e + ( $Z_i - $Z_e ) / $dZ_dX_i;
        $Y_i = $Y_e + $dY_dX_i * ( $X_i - $X_e );
    }
    elsif ( $icol == 5 ) {

        # Given W=ln(TOA). We have to iterate in this case because assuming
        # constant dW/dX is not sufficiently accurate.  Since Z hardly changes
        # with distance, an accurate first guess is made using the last Z in
        # the table. Simple iteration converges in just a couple of steps.
        $Z_i = $Z_e;
        $X_i = $X_e;
        foreach my $it ( 0 .. 5 ) {
            my $X_last = $X_i;
            $X_i = log( exp($Q) + exp($Z_i) );
            $Z_i = $Z_e + $dZ_dX_i * ( $X_i - $X_e );
            my $dX = $X_i - $X_last;
            ##print STDERR "it=$it, X_i=$X_i, Z_i=$Z_i, dX=$dX\n";
            last if ( abs($dX) < 1.e-8 );
        }
        $Y_i = $Y_e + $dY_dX_i * ( $X_i - $X_e );
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
        P8000_G1X1   => [ 0, 1.1,   1.0e-6 ],
        P8000_G1X2   => [ 0, 1.2,   1.e-6 ],
        P8000_G1X3   => [ 0, 1.3,   1.1e-6 ],
        P8000_G1X4   => [ 0, 1.4,   1.1e-6 ],
        P8000_G1X667 => [ 0, 1.667, 1.1e-6 ],
        P8000_G2     => [ 0, 2,     1.3e-6 ],
        P8000_G3     => [ 0, 3,     2.7e-6 ],
        P4000_G7     => [ 0, 7,     2.5e-6 ],
        C8000_G1X1   => [ 1, 1.1,   1.8e-6 ],
        C8000_G1X2   => [ 1, 1.2,   1.0e-6 ],
        C8000_G1X3   => [ 1, 1.3,   8.e-7 ],
        C16000_G1X4  => [ 1, 1.4,   8.e-7 ],
        C8000_G1X667 => [ 1, 1.667, 6.8e-7 ],
        C8000_G2     => [ 1, 2,     8.e-7 ],
        C8000_G3     => [ 1, 3,     1.6e-6 ],
        C4000_G7     => [ 1, 7,     1.6e-6 ],
        S8000_G1X1   => [ 2, 1.1,   2.6e-6 ],
        S8000_G1X15  => [ 2, 1.15,  1.0e-6 ],
        S8000_G1X2   => [ 2, 1.2,   1.1e-6 ],
        S8000_G1X3   => [ 2, 1.3,   8.9e-7 ],
        S16000_G1X4  => [ 2, 1.4,   9.e-7 ],
        S32000_G1X4  => [ 2, 1.4,   7.2e-7 ],
        S8000_G1X667 => [ 2, 1.667, 8.2e-7 ],
        S8000_G2     => [ 2, 2,     8.3e-7 ],
        S8000_G3     => [ 2, 3,     1.5e-6 ],
    };

    # Let..
    #  P0 = ambient atmospheric pressure
    #  E  = explosion energy
    #  n  = symmetry = 0,1, or 2  for plane, cylindrical or spherical
    #  d = (E/P0)^(1/(n+1)) = scaling distance
    #  lambda = scaled range = r/d, where r is distance
    #  tau = scaled time = c0 t / d, where t is time and c0 is ambient sound
    #  speed

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

    # These errors can be made almost arbitrarily small by increasing the
    # number of points in the FD calculation and in the tables.  My goal for
    # these tables was to achieve a maximum error on the order of 1.e-6.

    $rtables = {

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 2.14e-07
        #    Est Max FD Error (based on N/2 run) to R=10.00: 1.2e-7
        #    Est Max MOC error for R>10.00: 1.97e-07
        #    Max interpolation error with cubic splines: 5.0522e-07
        #    Est Max overall error after cubic interpolation: 8.2e-7

        S8000_G1X667 => [
            [ -4.47131967, 12.00017656,  -2.99997696, -4.47242876, 0.99833545 ],
            [ -4.34086863, 11.60882679,  -2.99995356, -4.34221760, 0.99797521 ],
            [ -4.18432739, 11.13920907,  -2.99992945, -4.18603368, 0.99743842 ],
            [ -4.00830247, 10.61114536,  -2.99990805, -4.01052493, 0.99666270 ],
            [ -3.90394164, 10.29807267,  -2.99990634, -3.90654119, 0.99609575 ],
            [ -3.80898673, 10.01321978,  -2.99989199, -3.81198478, 0.99549639 ],
            [ -3.75423088, 9.84896131,   -2.99982618, -3.75748598, 0.99510968 ],
            [ -3.36690888, 8.68712263,   -2.99944452, -3.37273544, 0.99123623 ],
            [ -3.07655735, 7.81632530,   -2.99867172, -3.08557716, 0.98641492 ],
            [ -2.82203786, 7.05327322,   -2.99715821, -2.83527547, 0.98003001 ],
            [ -2.61566768, 6.43497396,   -2.99473421, -2.63374459, 0.97268482 ],
            [ -2.43594758, 5.89706555,   -2.99100619, -2.45966996, 0.96409622 ],
            [ -2.27481527, 5.41552665,   -2.98550030, -2.30509490, 0.95410350 ],
            [ -2.12902364, 4.98079254,   -2.97772943, -2.16679676, 0.94267606 ],
            [ -1.99510137, 4.58267559,   -2.96709304, -2.04139114, 0.92969842 ],
            [ -1.87013712, 4.21273373,   -2.95285423, -1.92610224, 0.91499650 ],
            [ -1.75144347, 3.86330586,   -2.93405429, -1.81846069, 0.89829516 ],
            [ -1.63666896, 3.52789490,   -2.90945043, -1.71642463, 0.87922099 ],
            [ -1.52265719, 3.19794423,   -2.87710678, -1.61741379, 0.85707375 ],
            [ -1.40441145, 2.86020889,   -2.83346983, -1.51760582, 0.83042664 ],
            [ -1.27454318, 2.49606791,   -2.77187453, -1.41190242, 0.79663364 ],
            [ -1.15783906, 2.17645070,   -2.70344383, -1.32090150, 0.76227052 ],
            [ -1.04989503, 1.88853752,   -2.62940478, -1.24047575, 0.72741128 ],
            [ -0.94788567, 1.62424540,   -2.55110531, -1.16805218, 0.69220861 ],
            [ -0.85100924, 1.38094967,   -2.47091140, -1.10267575, 0.65728555 ],
            [ -0.75583037, 1.14967237,   -2.38850323, -1.04178420, 0.62214644 ],
            [ -0.66067440, 0.92638273,   -2.30452529, -0.98426559, 0.58680099 ],
            [ -0.56538356, 0.71078418,   -2.22073602, -0.93002273, 0.55177668 ],
            [ -0.46894693, 0.50063867,   -2.13792894, -0.87848562, 0.51724522 ],
            [ -0.37052079, 0.29423388,   -2.05684762, -0.82925362, 0.48340904 ],
            [ -0.26922237, 0.08991267,   -1.97804738, -0.78197219, 0.45043537 ],
            [ -0.16413266, -0.11390880,  -1.90194577, -0.73633655, 0.41846772 ],
            [ -0.05450900, -0.31834833,  -1.82900029, -0.69217375, 0.38769107 ],
            [ 0.06046656,  -0.52457702,  -1.75952995, -0.64931888, 0.35824852 ],
            [ 0.18184761,  -0.73407809,  -1.69366519, -0.60756663, 0.33021201 ],
            [ 0.31061902,  -0.94809318,  -1.63155988, -0.56678778, 0.30366628 ],
            [ 0.44790054,  -1.16798790,  -1.57328453, -0.52685325, 0.27865976 ],
            [ 0.59490197,  -1.39516909,  -1.51886836, -0.48765220, 0.25522014 ],
            [ 0.75092634,  -1.62815515,  -1.46889949, -0.44955775, 0.23361331 ],
            [ 0.91295437,  -1.86244672,  -1.42420978, -0.41331600, 0.21421169 ],
            [ 1.08238397,  -2.10026209,  -1.38406074, -0.37854189, 0.19670579 ],
            [ 1.26035692,  -2.34328759,  -1.34791333, -0.34497817, 0.18086815 ],
            [ 1.44821660,  -2.59336101,  -1.31528679, -0.31238468, 0.16649363 ],
            [ 1.64717352,  -2.85203589,  -1.28580978, -0.28059304, 0.15342303 ],
            [ 1.85838466,  -3.12072444,  -1.25917527, -0.24947724, 0.14152449 ],
            [ 2.08336066,  -3.40122415,  -1.23508204, -0.21889045, 0.13066789 ],
            [ 2.32339214,  -3.69499592,  -1.21329977, -0.18874674, 0.12075446 ],
            [ 2.58016892,  -4.00394230,  -1.19359658, -0.15893468, 0.11168451 ],
            [ 2.85527744,  -4.32979096,  -1.17578408, -0.12938179, 0.10337820 ],
            [ 3.15059686,  -4.67457670,  -1.15968168, -0.10000671, 0.09575934 ],
            [ 3.46799715,  -5.04028422,  -1.14513448, -0.07075197, 0.08876397 ],
            [ 3.80994126,  -5.42954181,  -1.13198398, -0.04152921, 0.08232638 ],
            [ 4.17868276,  -5.84469602,  -1.12010355, -0.01229407, 0.07639637 ],
            [ 4.57706811,  -6.28872605,  -1.10936515, 0.01702229,  0.07092308 ],
            [ 5.00813597,  -6.76478221,  -1.09965573, 0.04647574,  0.06586293 ],
            [ 5.47551790,  -7.27662574,  -1.09086735, 0.07613458,  0.06117462 ],
            [ 5.98283886,  -7.82796598,  -1.08290955, 0.10603769,  0.05682553 ],
            [ 6.53451816,  -8.42333411,  -1.07569399, 0.13624303,  0.05278351 ],
            [ 7.13537008,  -9.06763910,  -1.06914293, 0.16679822,  0.04902121 ],
            [ 7.79080481,  -9.76638120,  -1.06318628, 0.19774905,  0.04551444 ],
            [ 8.50682928,  -10.52564772, -1.05776188, 0.22913618,  0.04224216 ],
            [ 9.28964796,  -11.35168960, -1.05281748, 0.26097706,  0.03918773 ],
            [ 10.14666359, -12.25197783, -1.04830341, 0.29330650,  0.03633444 ],
            [ 11.08527773, -13.23393648, -1.04417965, 0.32612673,  0.03366986 ],
            [ 12.11369120, -14.30578802, -1.04040991, 0.35943996,  0.03118236 ],
            [ 13.24050442, -15.47613513, -1.03696295, 0.39323418,  0.02886196 ],
            [ 14.47551756, -16.75479211, -1.03380948, 0.42750705,  0.02669851 ],
            [ 15.82873069, -18.15174622, -1.03092482, 0.46223551,  0.02468359 ],
            [ 17.31134383, -19.67819384, -1.02828581, 0.49740384,  0.02280856 ],
            [ 18.88832591, -21.29787507, -1.02593576, 0.53199868,  0.02111200 ],
            [ 20.62829918, -23.08102816, -1.02376103, 0.56730907,  0.01951828 ],
            [ 22.47666448, -24.97146735, -1.02182006, 0.60201764,  0.01807571 ],
            [ 24.57274049, -27.11128586, -1.01997232, 0.63840684,  0.01668379 ],
            [ 26.79231697, -29.37331798, -1.01833054, 0.67401031,  0.01543100 ],
            [ 29.12104014, -31.74298790, -1.01687645, 0.70860383,  0.01430827 ],
            [ 31.75254986, -34.41703052, -1.01548928, 0.74479270,  0.01322510 ],
            [ 34.74059872, -37.44932923, -1.01416806, 0.78270772,  0.01218190 ],
            [ 37.87220238, -40.62343037, -1.01300624, 0.81936585,  0.01125481 ],
            [ 41.42775859, -44.22320424, -1.01189914, 0.85775203,  0.01036248 ],
            [ 45.12637165, -47.96398797, -1.01093179, 0.89458407,  0.00957532 ],
            [ 49.32169678, -52.20318341, -1.01000930, 0.93312234,  0.00881787 ],
            [ 53.64336503, -56.56633145, -1.00920917, 0.96976052,  0.00815526 ],
            [ 58.53542848, -61.50152425, -1.00844530, 1.00805366,  0.00751755 ],
            [ 63.51102490, -66.51746290, -1.00778849, 1.04404789,  0.00696504 ],
            [ 69.12428190, -72.17262660, -1.00716047, 1.08161094,  0.00643295 ],
            [ 75.48753090, -78.57947626, -1.00656093, 1.12087224,  0.00592139 ],
            [ 81.88429598, -85.01653747, -1.00605171, 1.15731476,  0.00548399 ],
            [ 89.10090155, -92.27500893, -1.00556459, 1.19533038,  0.00506296 ],
            [ 97.28202476, -100.49969830, -1.00509937, 1.23504959, 0.00465837 ],
            [
                105.37058788, -108.62789211, -1.00471013, 1.27131691,
                0.00431790
            ],
            [
                114.47879249, -117.77725333, -1.00433738, 1.30911264,
                0.00399010
            ],
            [
                124.78374849, -128.12501709, -1.00398100, 1.34856201,
                0.00367502
            ],
            [
                136.50260420, -139.88847317, -1.00364086, 1.38980643,
                0.00337270
            ],
            [
                147.87455236, -151.30019806, -1.00336213, 1.42670857,
                0.00312375
            ],
            [
                160.68879541, -164.15576660, -1.00309516, 1.46516124,
                0.00288423
            ],
            [
                175.19786073, -178.70783464, -1.00283987, 1.50529288,
                0.00265415
            ],
            [
                188.80505695, -192.35225188, -1.00263598, 1.54012265,
                0.00246964
            ],
            [
                204.02114087, -207.60691643, -1.00244010, 1.57631322,
                0.00229172
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 1.8e-07
        #    Est Max FD Error (based on N/2 run) to R=30.01: 1.8e-7
        #    Est Max MOC error for R>30.01: 3.52e-07
        #    Max interpolation error with cubic splines: 5.061e-07
        #    Est Max overall error after cubic interpolation: 1.e-6

        C8000_G1X2 => [
            [ -7.08875546, 12.00017695,  -1.99998660, -7.09005067, 0.99870396 ],
            [ -7.00902916, 11.84072538,  -1.99998428, -7.01043193, 0.99859625 ],
            [ -6.18981905, 10.20233910,  -1.99991639, -6.19300432, 0.99680978 ],
            [ -5.59041613, 9.00363091,   -1.99972076, -5.59622393, 0.99417610 ],
            [ -5.10063688, 8.02430482,   -1.99925741, -5.11013134, 0.99046380 ],
            [ -4.72619201, 7.27582944,   -1.99842900, -4.74002583, 0.98608079 ],
            [ -4.40112819, 6.62642031,   -1.99699737, -4.42032029, 0.98065109 ],
            [ -4.12119034, 6.06766916,   -1.99476351, -4.14664772, 0.97428213 ],
            [ -3.87137381, 5.56972816,   -1.99141854, -3.90414766, 0.96682387 ],
            [ -3.64675168, 5.12290750,   -1.98665967, -3.68790156, 0.95826822 ],
            [ -3.44105335, 4.71488651,   -1.98009409, -3.49175666, 0.94850188 ],
            [ -3.24979982, 4.33698200,   -1.97125753, -3.31138050, 0.93739190 ],
            [ -3.06851207, 3.98062340,   -1.95952740, -3.14255877, 0.92470015 ],
            [ -2.89371888, 3.63939124,   -1.94412337, -2.98216482, 0.91013310 ],
            [ -2.72094124, 3.30515993,   -1.92386983, -2.82634115, 0.89316328 ],
            [ -2.54323669, 2.96559223,   -1.89667152, -2.66939024, 0.87273700 ],
            [ -2.34235329, 2.58841093,   -1.85685932, -2.49671367, 0.84572416 ],
            [ -2.16261618, 2.25853420,   -1.81238619, -2.34715271, 0.81792772 ],
            [ -1.99711693, 1.96248713,   -1.76411123, -2.21410850, 0.78941927 ],
            [ -1.84113556, 1.69125014,   -1.71285193, -2.09322074, 0.76027157 ],
            [ -1.69176072, 1.43932982,   -1.65955981, -1.98184822, 0.73067306 ],
            [ -1.54594802, 1.20131240,   -1.60479819, -1.87748525, 0.70064796 ],
            [ -1.40194799, 0.97420576,   -1.54933671, -1.77876610, 0.67039370 ],
            [ -1.25808284, 0.75531166,   -1.49377937, -1.68450256, 0.64007514 ],
            [ -1.11304900, 0.54267359,   -1.43871426, -1.59386541, 0.60990517 ],
            [ -0.96546577, 0.33436698,   -1.38457138, -1.50606897, 0.58006143 ],
            [ -0.81416942, 0.12891930,   -1.33176155, -1.42054338, 0.55075392 ],
            [ -0.65798700, -0.07503734,  -1.28061253, -1.33678048, 0.52218034 ],
            [ -0.49570899, -0.27880312,  -1.23138248, -1.25431522, 0.49452298 ],
            [ -0.32599377, -0.48372552,  -1.18425439, -1.17267745, 0.46793559 ],
            [ -0.14761737, -0.69089989,  -1.13942119, -1.09150901, 0.44258350 ],
            [ 0.04087092,  -0.90159243,  -1.09700053, -1.01039247, 0.41858886 ],
            [ 0.24108246,  -1.11714369,  -1.05707351, -0.92889117, 0.39604829 ],
            [ 0.45496149,  -1.33913857,  -1.01966937, -0.84648651, 0.37502190 ],
            [ 0.68448208,  -1.56907812,  -0.98483477, -0.76270107, 0.35557153 ],
            [ 0.92695126,  -1.80393456,  -0.95318713, -0.67866686, 0.33806421 ],
            [ 1.17833858,  -2.03993433,  -0.92513359, -0.59566470, 0.32272748 ],
            [ 1.44105019,  -2.27961029,  -0.90017772, -0.51269973, 0.30927900 ],
            [ 1.71671009,  -2.52461028,  -0.87799839, -0.42911370, 0.29753125 ],
            [ 2.00682461,  -2.77639515,  -0.85833485, -0.34432412, 0.28732831 ],
            [ 2.31389457,  -3.03720823,  -0.84091247, -0.25749545, 0.27850798 ],
            [ 2.63851050,  -3.30761886,  -0.82560268, -0.16835364, 0.27098303 ],
            [ 2.98317311,  -3.58978708,  -0.81219675, -0.07609470, 0.26462565 ],
            [ 3.34901383,  -3.88472120,  -0.80056728, 0.01970960,  0.25934773 ],
            [ 3.73865008,  -4.19462873,  -0.79054982, 0.11988332,  0.25504369 ],
            [ 4.15368787,  -4.52089939,  -0.78202387, 0.22499049,  0.25162735 ],
            [ 4.59655718,  -4.86558133,  -0.77485540, 0.33581405,  0.24900595 ],
            [ 5.06956475,  -5.23062831,  -0.76891941, 0.45311155,  0.24708975 ],
            [ 5.57620571,  -5.61891103,  -0.76408207, 0.57793928,  0.24578601 ],
            [ 6.12036019,  -6.03358387,  -0.76021594, 0.71144808,  0.24500504 ],
            [ 6.70749600,  -6.47899642,  -0.75719218, 0.85517714,  0.24465872 ],
            [ 7.34386811,  -6.96007454,  -0.75488954, 1.01085612,  0.24466306 ],
            [ 8.03912150,  -7.48428164,  -0.75318707, 1.18104269,  0.24493864 ],
            [ 8.80709139,  -8.06220406,  -0.75197198, 1.36932313,  0.24541343 ],
            [ 9.66680474,  -8.70829657,  -0.75114272, 1.58056798,  0.24602258 ],
            [ 10.64748204, -9.44463602,  -0.75060769, 1.82217731,  0.24670938 ],
            [ 11.79533878, -10.30601987, -0.75028731, 2.10578614,  0.24742475 ],
            [ 13.18918694, -11.35166803, -0.75011451, 2.45116920,  0.24812642 ],
            [ 14.97903703, -12.69417547, -0.75003482, 2.89589861,  0.24877716 ],
            [ 17.51170195, -14.59371676, -0.75000646, 3.52675467,  0.24934237 ],
            [ 21.91761237, -17.89815885, -0.75000037, 4.62646982,  0.24977997 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 2.22e-07
        #    Est Max FD Error (based on N/2 run) to R=50.01: 1.3e-7
        #    Est Max MOC error for R>50.01: -5.38e-08
        #    Max interpolation error with cubic splines: 5.0857e-07
        #    Est Max overall error after cubic interpolation: 7.e-7
        #
        C8000_G1X667 => [
            [ -6.55075459, 12.00017704,  -1.99998464, -6.55214114, 0.99861249 ],
            [ -5.93701885, 10.77272290,  -1.99994759, -5.93958175, 0.99743389 ],
            [ -5.70814167, 10.31498268,  -1.99992219, -5.71136475, 0.99677185 ],
            [ -5.07940127, 9.05759958,   -1.99972676, -5.08545358, 0.99393019 ],
            [ -4.58984597, 8.07871574,   -1.99927433, -4.59973865, 0.99006195 ],
            [ -4.21159395, 7.32262203,   -1.99845413, -4.22606456, 0.98543583 ],
            [ -3.88216309, 6.66448012,   -1.99702002, -3.90232950, 0.97966028 ],
            [ -3.60032167, 6.10192346,   -1.99478174, -3.62712652, 0.97290608 ],
            [ -3.35117699, 5.60531502,   -1.99145518, -3.38566734, 0.96506345 ],
            [ -3.12711474, 5.15959719,   -1.98672176, -3.17040218, 0.95606720 ],
            [ -2.92220023, 4.75311425,   -1.98019764, -2.97550383, 0.94581508 ],
            [ -2.73184250, 4.37695475,   -1.97141909, -2.79653343, 0.93416967 ],
            [ -2.55212240, 4.02363873,   -1.95980489, -2.62979837, 0.92093394 ],
            [ -2.37945674, 3.68649557,   -1.94459981, -2.47205297, 0.90582446 ],
            [ -2.20964574, 3.35789353,   -1.92469793, -2.31967828, 0.88835821 ],
            [ -2.03686810, 3.02753489,   -1.89825637, -2.16793478, 0.86764400 ],
            [ -1.84673432, 2.67003595,   -1.86066457, -2.00542091, 0.84116832 ],
            [ -1.66712619, 2.33970955,   -1.81614337, -1.85685875, 0.81254537 ],
            [ -1.50261364, 2.04482503,   -1.76760190, -1.72554605, 0.78339658 ],
            [ -1.34822283, 1.77584562,   -1.71587782, -1.60685575, 0.75380822 ],
            [ -1.20162134, 1.52818062,   -1.66224202, -1.49850589, 0.72412950 ],
            [ -1.05952489, 1.29585986,   -1.60728372, -1.39771577, 0.69436410 ],
            [ -0.91776436, 1.07199782,   -1.55087913, -1.30141813, 0.66419599 ],
            [ -0.77632182, 0.85664001,   -1.49434507, -1.20960206, 0.63413212 ],
            [ -0.63380037, 0.64767419,   -1.43831147, -1.12135528, 0.60435631 ],
            [ -0.48865456, 0.44293720,   -1.38321087, -1.03578027, 0.57499694 ],
            [ -0.33976671, 0.24102961,   -1.32952955, -0.95232900, 0.54625531 ],
            [ -0.18585584, 0.04044560,   -1.27760080, -0.87043069, 0.51828855 ],
            [ -0.02561942, -0.16021613,  -1.22768743, -0.78957769, 0.49124507 ],
            [ 0.14209644,  -0.36205851,  -1.18004815, -0.70939824, 0.46529219 ],
            [ 0.31879321,  -0.56649935,  -1.13480704, -0.62940819, 0.44053728 ],
            [ 0.50601564,  -0.77488002,  -1.09207165, -0.54916835, 0.41708496 ],
            [ 0.70535246,  -0.98848281,  -1.05193695, -0.46827352, 0.39503433 ],
            [ 0.91866510,  -1.20878166,  -1.01444700, -0.38625505, 0.37445415 ],
            [ 1.14821726,  -1.43754929,  -0.97959904, -0.30254413, 0.35538290 ],
            [ 1.38610500,  -1.66680560,  -0.94863878, -0.22006169, 0.33853059 ],
            [ 1.63384848,  -1.89832409,  -0.92111472, -0.13808644, 0.32366411 ],
            [ 1.89320893,  -2.13396205,  -0.89662853, -0.05588845, 0.31057267 ],
            [ 2.16599831,  -2.37549841,  -0.87485355, 0.02721666,  0.29908025 ],
            [ 2.45375152,  -2.62438012,  -0.85554184, 0.11178815,  0.28905011 ],
            [ 2.75810308,  -2.88208938,  -0.83847608, 0.19839336,  0.28035959 ],
            [ 3.08072139,  -3.15009933,  -0.82346828, 0.28759557,  0.27289975 ],
            [ 3.42348042,  -3.43002708,  -0.81034679, 0.38000742,  0.26656878 ],
            [ 3.78818345,  -3.72341488,  -0.79896354, 0.47622033,  0.26127535 ],
            [ 4.17669322,  -4.03184991,  -0.78918132, 0.57684630,  0.25693191 ],
            [ 4.59098988,  -4.35701607,  -0.78086872, 0.68253538,  0.25345201 ],
            [ 5.03350854,  -4.70095755,  -0.77389301, 0.79406061,  0.25074751 ],
            [ 5.50695662,  -5.06592971,  -0.76812559, 0.91226833,  0.24873108 ],
            [ 6.01462862,  -5.45463861,  -0.76343792, 1.03815490,  0.24731505 ],
            [ 6.56040579,  -5.87023349,  -0.75970349, 1.17286311,  0.24641258 ],
            [ 7.15075820,  -6.31781728,  -0.75679003, 1.31817220,  0.24593733 ],
            [ 7.79134123,  -6.80185264,  -0.75458237, 1.47565738,  0.24580882 ],
            [ 8.49260154,  -7.33040223,  -0.75295886, 1.64806985,  0.24594943 ],
            [ 9.26877549,  -7.91434766,  -0.75180860, 1.83909278,  0.24628893 ],
            [ 10.14009105, -8.56904080,  -0.75103096, 2.05389184,  0.24676441 ],
            [ 11.13756952, -9.31790498,  -0.75053576, 2.30031370,  0.24732086 ],
            [ 12.31082700, -10.19828444, -0.75024484, 2.59084005,  0.24791098 ],
            [ 13.74647604, -11.27524610, -0.75009240, 2.94719089,  0.24849457 ],
            [ 15.61512788, -12.67683216, -0.75002544, 3.41208208,  0.24903696 ],
            [ 18.33439732, -14.71631517, -0.75000376, 4.08998728,  0.24950631 ],
            [ 23.48552631, -18.57966703, -0.75000009, 5.37633571,  0.24986309 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 2.08e-07
        #    Est Max FD Error (based on N/2 run) to R=24.70: 1.57e-7
        #    Est Max MOC error for R>24.70: 4.69e-07
        #    Max interpolation error with cubic splines: 5.0545e-07
        #    Est Max overall error after cubic interpolation: 1.1e-6
        'P8000_G1X667' => [
            [
                -11.89911045, 12.00006630, -0.99999232, -11.90095972,
                0.99907452
            ],
            [
                -11.08146799, 11.18243389, -0.99998260, -11.08425250,
                0.99860583
            ],
            [ -9.72600894, 9.82702826,  -0.99992814, -9.73150001, 0.99724714 ],
            [ -8.63033221, 8.73149503,  -0.99978450, -8.63984739, 0.99522090 ],
            [ -7.76073329, 7.86219463,  -0.99948631, -7.77546670, 0.99258319 ],
            [ -7.04420438, 7.14620326,  -0.99895005, -7.06534645, 0.98932947 ],
            [ -6.43427580, 6.53715579,  -0.99807288, -6.46305111, 0.98543607 ],
            [ -5.90216994, 6.00640104,  -0.99673226, -5.93985311, 0.98087198 ],
            [ -5.42795523, 5.53416320,  -0.99478041, -5.47591048, 0.97558767 ],
            [ -4.99819465, 5.10719333,  -0.99204444, -5.05789509, 0.96952674 ],
            [ -4.60247451, 4.71531326,  -0.98831637, -4.67555552, 0.96261029 ],
            [ -4.23239668, 4.35042778,  -0.98334189, -4.32072675, 0.95473070 ],
            [ -3.88068277, 4.00566357,  -0.97680056, -3.98646959, 0.94573820 ],
            [ -3.53971847, 3.67399729,  -0.96825163, -3.66572136, 0.93539489 ],
            [ -3.20020985, 3.34708905,  -0.95702186, -3.35014730, 0.92329380 ],
            [ -2.84569657, 3.01039616,  -0.94180091, -3.02537553, 0.90854528 ],
            [ -2.46408530, 2.65483553,  -0.92084162, -2.68209551, 0.89011696 ],
            [ -2.11984166, 2.34172143,  -0.89761973, -2.37886945, 0.87122886 ],
            [ -1.80060054, 2.05907705,  -0.87257657, -2.10377386, 0.85194082 ],
            [ -1.50474115, 1.80467776,  -0.84679044, -1.85452671, 0.83277983 ],
            [ -1.22460225, 1.57108557,  -0.82067841, -1.62387841, 0.81377732 ],
            [ -0.94911738, 1.34865442,  -0.79405751, -1.40233137, 0.79458000 ],
            [ -0.66996832, 1.13079568,  -0.76683715, -1.18326622, 0.77493611 ],
            [ -0.38127038, 0.91342800,  -0.73914567, -0.96246368, 0.75475902 ],
            [ -0.08664585, 0.69968820,  -0.71202355, -0.74307203, 0.73464275 ],
            [ 0.21767968,  0.48704903,  -0.68575283, -0.52256532, 0.71466436 ],
            [ 0.53544319,  0.27320339,  -0.66059871, -0.29864114, 0.69491806 ],
            [ 0.87107568,  0.05556408,  -0.63676795, -0.06870848, 0.67548110 ],
            [ 1.22945181,  -0.16854460, -0.61445302, 0.16990306,  0.65644500 ],
            [ 1.61671631,  -0.40239017, -0.59380221, 0.42046077,  0.63788749 ],
            [ 2.02653226,  -0.64187388, -0.57551000, 0.67822853,  0.62044620 ],
            [ 2.43538762,  -0.87399655, -0.56045880, 0.92870442,  0.60514144 ],
            [ 2.85182629,  -1.10469471, -0.54793085, 1.17780451,  0.59150415 ],
            [ 3.27929866,  -1.33661712, -0.53753432, 1.42799222,  0.57932934 ],
            [ 3.72801836,  -1.57579032, -0.52883033, 1.68541023,  0.56829246 ],
            [ 4.19544411,  -1.82124596, -0.52170966, 1.94868138,  0.55843910 ],
            [ 4.68884723,  -2.07716326, -0.51590547, 2.22197446,  0.54959465 ],
            [ 5.21476040,  -2.34719519, -0.51122877, 2.50886438,  0.54165625 ],
            [ 5.77458880,  -2.63231346, -0.50755731, 2.81006686,  0.53461758 ],
            [ 6.37712996,  -2.93723533, -0.50472796, 3.13025511,  0.52838554 ],
            [ 7.02816668,  -3.26510028, -0.50261877, 3.47241241,  0.52292820 ],
            [ 7.74115401,  -3.62287939, -0.50110075, 3.84349040,  0.51817155 ],
            [ 8.52783109,  -4.01664352, -0.50007217, 4.24944716,  0.51408751 ],
            [ 9.40875574,  -4.45685750, -0.49943457, 4.70072019,  0.51062645 ],
            [ 10.41164153, -4.95754115, -0.49910362, 5.21129662,  0.50774855 ],
            [ 11.55874472, -5.52998864, -0.49900614, 5.79233147,  0.50544447 ],
            [ 12.95640742, -6.22746671, -0.49907908, 6.49738000,  0.50359206 ],
            [ 14.69741181, -7.09652532, -0.49926981, 7.37280379,  0.50218856 ],
            [ 17.18460718, -8.33865386, -0.49953858, 8.62035517,  0.50111934 ],
            [ 19.23703232, -9.36409667, -0.49970334, 9.64834868,  0.50065893 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 2.32e-07
        #    Est Max FD Error (based on N/2 run) to R=34.09: 1.7e-7
        #    Est Max MOC error for R>34.09: 4.44e-07
        #    Max interpolation error with cubic splines: 5.0693e-07
        #    Est Max overall error after cubic interpolation: 1.1e-6

        'P8000_G1X3' => [
            [
                -12.64178695, 12.00006627, -0.99999305, -12.64354541,
                0.99912001
            ],
            [
                -11.55327663, 10.91157197, -0.99997938, -11.55630893,
                0.99848159
            ],
            [
                -10.65825935, 10.01659000, -0.99994115, -10.66300708,
                0.99762065
            ],
            [ -9.59080036, 8.94924522,   -0.99982596, -9.59890966, 0.99592965 ],
            [ -8.68066644, 8.03936958,   -0.99956787, -8.69347716, 0.99355652 ],
            [ -7.93379227, 7.29297481,   -0.99908955, -7.95245205, 0.99059200 ],
            [ -7.30382014, 6.66380066,   -0.99829487, -7.32946634, 0.98703558 ],
            [ -6.75693356, 6.11815278,   -0.99706518, -6.79075930, 0.98285381 ],
            [ -6.27107344, 5.63412334,   -0.99525652, -6.31435974, 0.97799810 ],
            [ -5.83201892, 5.19767299,   -0.99270136, -5.88614365, 0.97241778 ],
            [ -5.42879645, 4.79805498,   -0.98919855, -5.49528385, 0.96604081 ],
            [ -5.05267196, 4.42682535,   -0.98450352, -5.13325580, 0.95876920 ],
            [ -4.69591917, 4.07664834,   -0.97830548, -4.79264998, 0.95046025 ],
            [ -4.35092377, 3.74047025,   -0.97018510, -4.46634895, 0.94089630 ],
            [ -4.00839182, 3.40989748,   -0.95950522, -4.14592793, 0.92970202 ],
            [ -3.65298273, 3.07133529,   -0.94507651, -3.81785851, 0.91609640 ],
            [ -3.26174160, 2.70541667,   -0.92464052, -3.46277884, 0.89859882 ],
            [ -2.91002769, 2.38408319,   -0.90191320, -3.14983660, 0.88056989 ],
            [ -2.58473991, 2.09461111,   -0.87733764, -2.86636424, 0.86205145 ],
            [ -2.28254156, 1.83327410,   -0.85185747, -2.60863366, 0.84344921 ],
            [ -1.99745286, 1.59406447,   -0.82605265, -2.37079919, 0.82490742 ],
            [ -1.71773896, 1.36667354,   -0.79971884, -2.14268615, 0.80605299 ],
            [ -1.43520462, 1.14452966,   -0.77279280, -1.91768532, 0.78664439 ],
            [ -1.14373020, 0.92329270,   -0.74537858, -1.69132750, 0.76656467 ],
            [ -0.84664221, 0.70588011,   -0.71846658, -1.46660131, 0.74636561 ],
            [ -0.54041311, 0.48991048,   -0.69235534, -1.24115884, 0.72613491 ],
            [ -0.22123999, 0.27299050,   -0.66729712, -1.01264381, 0.70595787 ],
            [ 0.11499845,  0.05269475,   -0.64351390, -0.77867875, 0.68592773 ],
            [ 0.47321343,  -0.17373041,  -0.62118424, -0.53656383, 0.66613196 ],
            [ 0.85920629,  -0.40939656,  -0.60046444, -0.28326103, 0.64666655 ],
            [ 1.27681726,  -0.65610105,  -0.58162474, -0.01722770, 0.62777585 ],
            [ 1.69095579,  -0.89367272,  -0.56617919, 0.23924942,  0.61117135 ],
            [ 2.11129741,  -1.12886858,  -0.55332639, 0.49296300,  0.59632952 ],
            [ 2.54389150,  -1.36582983,  -0.54259365, 0.74797681,  0.58297627 ],
            [ 2.99347481,  -1.60768449,  -0.53365117, 1.00730178,  0.57094171 ],
            [ 3.46423937,  -1.85709520,  -0.52624595, 1.27346514,  0.56011071 ],
            [ 3.96072428,  -2.11679377,  -0.52016501, 1.54907159,  0.55039031 ],
            [ 4.48662276,  -2.38899207,  -0.51523694, 1.83617492,  0.54172591 ],
            [ 5.04836821,  -2.67726026,  -0.51129007, 2.13825821,  0.53403750 ],
            [ 5.65317286,  -2.98549919,  -0.50817929, 2.45912740,  0.52726605 ],
            [ 6.30543724,  -3.31614102,  -0.50578921, 2.80105756,  0.52139761 ],
            [ 7.01807119,  -3.67589994,  -0.50398736, 3.17075062,  0.51635438 ],
            [ 7.80467699,  -4.07178280,  -0.50266718, 3.57516273,  0.51209219 ],
            [ 8.68181627,  -4.51224958,  -0.50173312, 4.02271354,  0.50857387 ],
            [ 9.67724551,  -5.01134309,  -0.50109395, 4.52746779,  0.50574488 ],
            [ 10.83007401, -5.58875040,  -0.50066965, 5.10915199,  0.50355856 ],
            [ 12.24501854, -6.29694707,  -0.50038347, 5.82039808,  0.50192868 ],
            [ 14.09849459, -7.22419190,  -0.50018548, 6.74955952,  0.50082944 ],
            [ 16.50506588, -8.42775810,  -0.50006569, 7.95403569,  0.50026241 ],
            [ 19.33335373, -9.84200507,  -0.50001718, 9.36858084,  0.50006512 ],
            [ 23.28159166, -11.81615394, -0.50000244, 11.34281230, 0.50000910 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 1.99e-07
        #    Est Max FD Error (based on N/2 run) to R=10.00: 2.2e-7
        #    Est Max MOC error for R>10.00: 1.75e-07
        #    Max interpolation error with cubic splines: 5.083e-07
        #    Est Max overall error after cubic interpolation: 9.e-7
        S8000_G1X3 => [
            [ -4.70223143, 12.00017664,  -2.99997917, -4.70328608, 0.99841721 ],
            [ -4.45655798, 11.26316300,  -2.99995647, -4.45808291, 0.99771089 ],
            [ -3.96633083, 9.79252728,   -2.99982473, -3.96951469, 0.99521686 ],
            [ -3.59782065, 8.68711660,   -2.99946310, -3.60336057, 0.99166845 ],
            [ -3.29907098, 7.79112733,   -2.99868401, -3.30775561, 0.98692171 ],
            [ -3.05037420, 7.04552293,   -2.99723068, -3.06300758, 0.98094626 ],
            [ -2.84019686, 6.41580123,   -2.99481177, -2.85754665, 0.97379094 ],
            [ -2.66316417, 5.88591324,   -2.99120708, -2.68583804, 0.96569558 ],
            [ -2.50091219, 5.40098954,   -2.98578028, -2.52989989, 0.95607931 ],
            [ -2.35433674, 4.96386969,   -2.97811598, -2.39053789, 0.94508483 ],
            [ -2.21961152, 4.56330574,   -2.96760331, -2.26402517, 0.93257617 ],
            [ -2.09375147, 4.19063983,   -2.95349480, -2.14751722, 0.91836908 ],
            [ -1.97460992, 3.83980427,   -2.93492173, -2.03903302, 0.90226365 ],
            [ -1.85837973, 3.50003452,   -2.91038224, -1.93521383, 0.88367154 ],
            [ -1.74293336, 3.16582162,   -2.87808721, -1.83441451, 0.86202626 ],
            [ -1.62267206, 2.82222056,   -2.83428678, -1.73228513, 0.83578019 ],
            [ -1.49197600, 2.45562314,   -2.77317689, -1.62515954, 0.80273422 ],
            [ -1.37431303, 2.13319451,   -2.70533846, -1.53265795, 0.76896018 ],
            [ -1.26536122, 1.84234990,   -2.63201884, -1.45072793, 0.73453550 ],
            [ -1.16231806, 1.57506888,   -2.55454680, -1.37682103, 0.69960832 ],
            [ -1.06377994, 1.32724346,   -2.47471494, -1.30959833, 0.66457330 ],
            [ -0.96680947, 1.09123175,   -2.39255768, -1.24686805, 0.62911758 ],
            [ -0.87060488, 0.86504880,   -2.30946051, -1.18805265, 0.59358972 ],
            [ -0.77430626, 0.64665301,   -2.22652368, -1.13259559, 0.55827113 ],
            [ -0.67691079, 0.43381506,   -2.14449697, -1.07993160, 0.52335149 ],
            [ -0.57758884, 0.22484403,   -2.06409310, -1.02966720, 0.48905549 ],
            [ -0.47544149, 0.01804122,   -1.98582576, -0.98143875, 0.45556359 ],
            [ -0.36973696, -0.18782376,  -1.91022771, -0.93502020, 0.42309443 ],
            [ -0.25959251, -0.39416833,  -1.83763286, -0.89016627, 0.39180135 ],
            [ -0.14420109, -0.60215298,  -1.76836555, -0.84671180, 0.36184462 ],
            [ -0.02259609, -0.81312368,  -1.70260690, -0.80447462, 0.33332678 ],
            [ 0.10617532,  -1.02829361,  -1.64052132, -0.76332318, 0.30634216 ],
            [ 0.24323289,  -1.24905499,  -1.58217498, -0.72311458, 0.28093846 ],
            [ 0.38979197,  -1.47684534,  -1.52759502, -0.68372418, 0.25714082 ],
            [ 0.54720056,  -1.71320357,  -1.47676662, -0.64503688, 0.23495021 ],
            [ 0.71050332,  -1.95055887,  -1.43130134, -0.60833182, 0.21507382 ],
            [ 0.88081751,  -2.19076871,  -1.39051473, -0.57326062, 0.19721338 ],
            [ 1.05954105,  -2.43591980,  -1.35378355, -0.53949077, 0.18109427 ],
            [ 1.24792769,  -2.68775069,  -1.32063391, -0.50678425, 0.16650555 ],
            [ 1.44720883,  -2.94786448,  -1.29067892, -0.47495520, 0.15327361 ],
            [ 1.65881894,  -3.21804003,  -1.26357071, -0.44382702, 0.14124171 ],
            [ 1.88391500,  -3.49963018,  -1.23905171, -0.41329843, 0.13029332 ],
            [ 2.12400764,  -3.79438033,  -1.21686365, -0.38324561, 0.12031188 ],
            [ 2.38037768,  -4.10370479,  -1.19680555, -0.35359899, 0.11120726 ],
            [ 2.65508639,  -4.42991320,  -1.17864849, -0.32422374, 0.10287711 ],
            [ 2.94960391,  -4.77456009,  -1.16223795, -0.29507684, 0.09525378 ],
            [ 3.26600246,  -5.13987596,  -1.14740460, -0.26607370, 0.08826365 ],
            [ 3.60634462,  -5.52803912,  -1.13400423, -0.23715483, 0.08184555 ],
            [ 3.97328522,  -5.94186341,  -1.12189148, -0.20823474, 0.07593830 ],
            [ 4.36927013,  -6.38388462,  -1.11094941, -0.17927038, 0.07049503 ],
            [ 4.79733801,  -6.85726647,  -1.10106061, -0.15019768, 0.06546905 ],
            [ 5.26072023,  -7.36534479,  -1.09212084, -0.12096572, 0.06082052 ],
            [ 5.76344177,  -7.91228173,  -1.08402780, -0.09150167, 0.05651039 ],
            [ 6.30952168,  -8.50218516,  -1.07669644, -0.06176375, 0.05250831 ],
            [ 6.90397417,  -9.14019253,  -1.07004314, -0.03168595, 0.04878399 ],
            [ 7.55200936,  -9.83159984,  -1.06399759, -0.00122611, 0.04531353 ],
            [ 8.25923415,  -10.58208161, -1.05849797, 0.02964619,  0.04207675 ],
            [ 9.03245304,  -11.39853555, -1.05348521, 0.06098122,  0.03905358 ],
            [ 9.87846879,  -12.28780751, -1.04891180, 0.09279497,  0.03622933 ],
            [ 10.80468299, -13.25733209, -1.04473565, 0.12509689,  0.03359102 ],
            [ 11.81949648, -14.31554853, -1.04091827, 0.15790157,  0.03112626 ],
            [ 12.93150970, -15.47106412, -1.03742765, 0.19120075,  0.02882521 ],
            [ 14.15032284, -16.73349001, -1.03423441, 0.22498974,  0.02667826 ],
            [ 15.48573596, -18.11261135, -1.03131364, 0.25924410,  0.02467736 ],
            [ 16.94894910, -19.61962937, -1.02864154, 0.29395207,  0.02281398 ],
            [ 18.46398782, -21.17625260, -1.02632006, 0.32723756,  0.02116827 ],
            [ 20.16438713, -22.91948902, -1.02412869, 0.36185492,  0.01959058 ],
            [ 21.94999295, -24.74640072, -1.02219122, 0.39553974,  0.01817520 ],
            [ 23.96945697, -26.80876150, -1.02034592, 0.43082668,  0.01680833 ],
            [ 26.07699099, -28.95741732, -1.01872299, 0.46493412,  0.01559024 ],
            [ 28.45908300, -31.38221348, -1.01717606, 0.50063369,  0.01441460 ],
            [ 31.16445423, -34.13200129, -1.01570398, 0.53805551,  0.01328191 ],
            [ 33.97931731, -36.98920907, -1.01441927, 0.57399703,  0.01228178 ],
            [ 37.17273303, -40.22665963, -1.01319558, 0.61163703,  0.01131854 ],
            [ 40.46265958, -43.55820572, -1.01213533, 0.64745032,  0.01047519 ],
            [ 44.18697233, -47.32577652, -1.01112416, 0.68490902,  0.00966289 ],
            [ 47.97483625, -51.15408741, -1.01025553, 0.72014251,  0.00895858 ],
            [ 52.24786766, -55.46911805, -1.00942568, 0.75693384,  0.00827982 ],
            [ 57.09122872, -60.35615806, -1.00863415, 0.79541048,  0.00762678 ],
            [ 61.95925407, -65.26453530, -1.00796237, 0.83114276,  0.00706803 ],
            [ 67.44996457, -70.79715241, -1.00732021, 0.86843414,  0.00652985 ],
            [ 73.67273999, -77.06351736, -1.00670735, 0.90741278,  0.00601235 ],
            [ 79.82307523, -83.25348737, -1.00619490, 0.94301661,  0.00557659 ],
            [ 86.74613536, -90.21769071, -1.00570444, 0.98013215,  0.00515681 ],
            [ 94.57538072, -98.08971042, -1.00523578, 1.01888214,  0.00475309 ],
            [
                103.47423954, -107.03311322, -1.00478870, 1.05940501,
                0.00436548
            ],
            [
                112.10488614, -115.70346731, -1.00442252, 1.09566840,
                0.00404611
            ],
            [
                121.82443781, -125.46425318, -1.00407193, 1.13346087,
                0.00373867
            ],
            [
                132.82216570, -136.50486741, -1.00373680, 1.17290781,
                0.00344318
            ],
            [
                143.12921188, -146.84901586, -1.00346923, 1.20714570,
                0.00320610
            ],
            [
                154.64688100, -158.40512454, -1.00321225, 1.24272245,
                0.00297737
            ],
            [
                167.57073757, -171.36886052, -1.00296576, 1.27974066,
                0.00275700
            ],
            [
                182.13734906, -185.97690621, -1.00272970, 1.31831522,
                0.00254501
            ],
            [
                198.63503301, -202.51770847, -1.00250399, 1.35857544,
                0.00234141
            ],
            [
                217.41803654, -221.34566066, -1.00228857, 1.40066755,
                0.00214623
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 2e-07
        #    Est Max FD Error (based on N/2 run) to R=10.00: 2.9e-7
        #    Est Max MOC error for R>10.00: 2.89e-07
        #    Max interpolation error with cubic splines: 5.0604e-07
        #    Est Max overall error after cubic interpolation: 1.1e-6
        S8000_G1X2 => [
            [ -4.75175969, 11.78518518,  -2.99997507, -4.75291337, 0.99826850 ],
            [ -4.61850927, 11.38543838,  -2.99996282, -4.61991838, 0.99788487 ],
            [ -4.18874869, 10.09618801,  -2.99986478, -4.19143514, 0.99596507 ],
            [ -3.74319962, 8.75966276,   -2.99950480, -3.74844720, 0.99210914 ],
            [ -3.42395850, 7.80220490,   -2.99871650, -3.43244200, 0.98722572 ],
            [ -3.20963809, 7.15963188,   -2.99754664, -3.22135505, 0.98233457 ],
            [ -2.99509372, 6.51673685,   -2.99534135, -3.01128943, 0.97554388 ],
            [ -2.80999882, 5.96260295,   -2.99191366, -2.83142173, 0.96760018 ],
            [ -2.64584326, 5.47184660,   -2.98684452, -2.67330817, 0.95840193 ],
            [ -2.49736423, 5.02886012,   -2.97963202, -2.53176004, 0.94784088 ],
            [ -2.36041864, 4.62145284,   -2.96963821, -2.40275718, 0.93574283 ],
            [ -2.23316051, 4.24434646,   -2.95621682, -2.28452113, 0.92203123 ],
            [ -2.11245132, 3.88852069,   -2.93842074, -2.17413856, 0.90641024 ],
            [ -1.99577322, 3.54696897,   -2.91503028, -2.06939655, 0.88849822 ],
            [ -1.87987895, 3.21084108,   -2.88415785, -1.96760414, 0.86759948 ],
            [ -1.75995358, 2.86735024,   -2.84245906, -1.86503165, 0.84237040 ],
            [ -1.62647922, 2.49178740,   -2.78253041, -1.75472297, 0.80968912 ],
            [ -1.50675672, 2.16253043,   -2.71573128, -1.65975086, 0.77620115 ],
            [ -1.39634923, 1.86659565,   -2.64339891, -1.57591179, 0.74202182 ],
            [ -1.29212328, 1.59501592,   -2.56673616, -1.50036744, 0.70723783 ],
            [ -1.19243937, 1.34306752,   -2.48738009, -1.43160235, 0.67218095 ],
            [ -1.09468508, 1.10388543,   -2.40567738, -1.36762251, 0.63667676 ],
            [ -0.99805053, 0.87539935,   -2.32303871, -1.30781563, 0.60108791 ],
            [ -0.90139805, 0.65487485,   -2.24034036, -1.25143696, 0.56560570 ],
            [ -0.80386770, 0.44038507,   -2.15846387, -1.19799371, 0.53048275 ],
            [ -0.70454575, 0.23002478,   -2.07805837, -1.14703333, 0.49592292 ],
            [ -0.60255930, 0.02212778,   -1.99966726, -1.09819539, 0.46212498 ],
            [ -0.49716466, -0.18458133,  -1.92382298, -1.05123892, 0.42931266 ],
            [ -0.38754100, -0.39142483,  -1.85090996, -1.00593441, 0.39766494 ],
            [ -0.27280557, -0.59972708,  -1.78121042, -0.96207578, 0.36732699 ],
            [ -0.15209637, -0.81066667,  -1.71497569, -0.91951045, 0.33843362 ],
            [ -0.02444472, -1.02551176,  -1.65236193, -0.87808863, 0.31107676 ],
            [ 0.11130328,  -1.24573175,  -1.59342464, -0.83764692, 0.28529950 ],
            [ 0.25619911,  -1.47252523,  -1.53825913, -0.79809677, 0.26115655 ],
            [ 0.41162712,  -1.70752051,  -1.48682738, -0.75929781, 0.23863763 ],
            [ 0.57418657,  -1.94534899,  -1.44037207, -0.72219979, 0.21828829 ],
            [ 0.74342401,  -2.18550195,  -1.39873375, -0.68684030, 0.20003598 ],
            [ 0.92079263,  -2.43018422,  -1.36124813, -0.65285603, 0.18358473 ],
            [ 1.10779424,  -2.68149115,  -1.32738157, -0.61995324, 0.16869430 ],
            [ 1.30593081,  -2.94137547,  -1.29671241, -0.58790288, 0.15517303 ],
            [ 1.51517774,  -3.20974143,  -1.26908868, -0.55674551, 0.14294874 ],
            [ 1.73798536,  -3.48963962,  -1.24405394, -0.52616821, 0.13181481 ],
            [ 1.97591948,  -3.78286807,  -1.22135697, -0.49604534, 0.12165576 ],
            [ 2.22960019,  -4.09003034,  -1.20085768, -0.46638763, 0.11240703 ],
            [ 2.50102598,  -4.41338673,  -1.18231496, -0.43705412, 0.10395974 ],
            [ 2.79219931,  -4.75513262,  -1.16553449, -0.40793942, 0.09622667 ],
            [ 3.10464219,  -5.11685911,  -1.15037538, -0.37900943, 0.08914624 ],
            [ 3.44082020,  -5.50121910,  -1.13666983, -0.35016203, 0.08264534 ],
            [ 3.80258811,  -5.91012612,  -1.12429820, -0.32137210, 0.07667445 ],
            [ 4.19319444,  -6.34703442,  -1.11311126, -0.29252581, 0.07117054 ],
            [ 4.61507827,  -6.81444288,  -1.10300670, -0.26359911, 0.06609392 ],
            [ 5.07167217,  -7.31592079,  -1.09387152, -0.23452069, 0.06140012 ],
            [ 5.56640190,  -7.85498636,  -1.08561035, -0.20524749, 0.05705371 ],
            [ 6.10368744,  -8.43619792,  -1.07812702, -0.17570567, 0.05301860 ],
            [ 6.68834364,  -9.06448676,  -1.07133799, -0.14583404, 0.04926471 ],
            [ 7.32518115,  -9.74473268,  -1.06517396, -0.11560259, 0.04576902 ],
            [ 8.02000732,  -10.48283678, -1.05956797, -0.08496335, 0.04250881 ],
            [ 8.77942698,  -11.28549524, -1.05446006, -0.05386782, 0.03946395 ],
            [ 9.61004314,  -12.15935362, -1.04980175, -0.02230081, 0.03661959 ],
            [ 10.51965753, -13.11227491, -1.04554719, 0.00976773,  0.03396091 ],
            [ 11.51567111, -14.15166026, -1.04166079, 0.04232382,  0.03147782 ],
            [ 12.60748436, -15.28696022, -1.03810597, 0.07539149,  0.02915800 ],
            [ 13.80349750, -16.52654653, -1.03485590, 0.10893598,  0.02699392 ],
            [ 15.11431062, -17.88104172, -1.03188234, 0.14296155,  0.02497572 ],
            [ 16.55092375, -19.36144424, -1.02916140, 0.17745456,  0.02309513 ],
            [ 18.12453690, -20.97892918, -1.02667263, 0.21238289,  0.02134509 ],
            [ 19.75522059, -22.65129290, -1.02450864, 0.24589639,  0.01979893 ],
            [ 21.48222690, -24.41890719, -1.02257194, 0.27884614,  0.01839484 ],
            [ 23.43112719, -26.40994776, -1.02072610, 0.31333871,  0.01703785 ],
            [ 25.64061869, -28.66323529, -1.01896958, 0.34949767,  0.01572863 ],
            [ 27.93523452, -30.99957060, -1.01743662, 0.38422591,  0.01457102 ],
            [ 30.53327206, -33.64095654, -1.01597627, 0.42059216,  0.01345459 ],
            [ 33.20424713, -36.35286706, -1.01471075, 0.45518810,  0.01247577 ],
            [ 36.22126125, -39.41239038, -1.01350349, 0.49136535,  0.01153173 ],
            [ 39.64554535, -42.88088698, -1.01235366, 0.52925349,  0.01062279 ],
            [ 43.13749700, -46.41420686, -1.01136726, 0.56493627,  0.00983503 ],
            [ 47.08735515, -50.40705174, -1.01042615, 0.60224371,  0.00907621 ],
            [ 51.04886676, -54.40824189, -1.00962717, 0.63687914,  0.00842617 ],
            [ 55.50727733, -58.90782404, -1.00886312, 0.67301244,  0.00779933 ],
            [ 60.54796598, -63.99129755, -1.00813362, 0.71076385,  0.00719583 ],
            [ 66.27535772, -69.76322299, -1.00743826, 0.75026888,  0.00661581 ],
            [ 71.95065698, -75.47904400, -1.00685754, 0.78639448,  0.00612767 ],
            [ 78.35556160, -81.92604265, -1.00630241, 0.82409599,  0.00565771 ],
            [ 85.61903720, -89.23331521, -1.00577262, 0.86350405,  0.00520601 ],
            [ 92.64749568, -96.30078086, -1.00533849, 0.89875216,  0.00483344 ],
            [
                100.54445752, -104.23821405, -1.00492262, 0.93546707,
                0.00447440
            ],
            [
                109.45747964, -113.19328916, -1.00452486, 0.97376622,
                0.00412894
            ],
            [
                119.56669894, -123.34627539, -1.00414506, 1.01378188,
                0.00379710
            ],
            [
                129.06326354, -132.88072909, -1.00384215, 1.04854689,
                0.00353101
            ],
            [
                139.69965774, -143.55640366, -1.00355151, 1.08470491,
                0.00327444
            ],
            [
                151.66409816, -155.56162464, -1.00327304, 1.12236487,
                0.00302740
            ],
            [
                165.18516895, -169.12509847, -1.00300665, 1.16164904,
                0.00278993
            ],
            [
                180.54265466, -184.52674769, -1.00275225, 1.20269539,
                0.00256205
            ],
            [
                194.38028121, -198.40107703, -1.00255730, 1.23690566,
                0.00238666
            ],
            [
                209.82746038, -213.88627542, -1.00236993, 1.27243344,
                0.00221743
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 1.87e-07
        #    Est Max FD Error (based on N/2 run) to R=50.01: 1.6e-7
        #    Est Max MOC error for R>50.01: 1.06e-07
        #    Max interpolation error with cubic splines: 5.0633e-07
        #    Est Max overall error after cubic interpolation: 8.e-7
        #
        C8000_G1X3 => [
            [ -6.90493912, 12.00017707,  -1.99998611, -6.90625760, 0.99868066 ],
            [ -6.61413969, 11.41858341,  -1.99997516, -6.61590354, 0.99823461 ],
            [ -6.05404004, 10.29841004,  -1.99992013, -6.05713028, 0.99690510 ],
            [ -5.36348937, 8.91742312,   -1.99969414, -5.36966297, 0.99380824 ],
            [ -4.89498539, 7.98065237,   -1.99922022, -4.90486546, 0.99007481 ],
            [ -4.51236634, 7.21586210,   -1.99832654, -4.52688201, 0.98539065 ],
            [ -4.19245337, 6.57678513,   -1.99683439, -4.21248917, 0.97979407 ],
            [ -3.91441910, 6.02189167,   -1.99450068, -3.94094715, 0.97319074 ],
            [ -3.66784985, 5.53050097,   -1.99104591, -3.70189329, 0.96552547 ],
            [ -3.44557920, 5.08845574,   -1.98614554, -3.48822516, 0.95673415 ],
            [ -3.24178445, 4.68433258,   -1.97940460, -3.29423390, 0.94670840 ],
            [ -3.05187461, 4.30923482,   -1.97034012, -3.11549351, 0.93529796 ],
            [ -2.87204254, 3.95592515,   -1.95835114, -2.94843114, 0.92229820 ],
            [ -2.69836910, 3.61711090,   -1.94261772, -2.78951192, 0.90737920 ],
            [ -2.52648726, 3.28490640,   -1.92194126, -2.63500356, 0.89000732 ],
            [ -2.34945457, 2.94701274,   -1.89417229, -2.47924637, 0.86910512 ],
            [ -2.15085879, 2.57466635,   -1.85392690, -2.30929201, 0.84175287 ],
            [ -1.97293030, 2.24866874,   -1.80903080, -2.16196890, 0.81367151 ],
            [ -1.80888672, 1.95581121,   -1.76033140, -2.03081312, 0.78493123 ],
            [ -1.65424905, 1.68752472,   -1.70869941, -1.91167320, 0.75563684 ],
            [ -1.50653277, 1.43902791,   -1.65524541, -1.80222148, 0.72606030 ],
            [ -1.36199577, 1.20373237,   -1.60027692, -1.69943617, 0.69607955 ],
            [ -1.21877957, 0.97853249,   -1.54450756, -1.60190739, 0.66585708 ],
            [ -1.07564931, 0.76146571,   -1.48871570, -1.50876813, 0.63564388 ],
            [ -0.93108940, 0.55027277,   -1.43341615, -1.41905933, 0.60560279 ],
            [ -0.78393613, 0.34336464,   -1.37912968, -1.33213849, 0.57595056 ],
            [ -0.63290457, 0.13910686,   -1.32622725, -1.24736733, 0.54686565 ],
            [ -0.47678519, -0.06389819,  -1.27503027, -1.16422703, 0.51853328 ],
            [ -0.31434835, -0.26695482,  -1.22579664, -1.08225346, 0.49112810 ],
            [ -0.14442504, -0.47118568,  -1.17875941, -1.00106849, 0.46482637 ],
            [ 0.03439926,  -0.67790785,  -1.13405480, -0.92022710, 0.43975756 ],
            [ 0.22363724,  -0.88843628,  -1.09178859, -0.83929706, 0.41603488 ],
            [ 0.42484018,  -1.10402485,  -1.05205723, -0.75787904, 0.39376266 ],
            [ 0.64002795,  -1.32632171,  -1.01487582, -0.67543420, 0.37299347 ],
            [ 0.87129961,  -1.55693310,  -0.98027529, -0.59145097, 0.35378136 ],
            [ 1.11305707,  -1.79007027,  -0.94921425, -0.50804559, 0.33668073 ],
            [ 1.36460578,  -2.02527857,  -0.92160150, -0.42529930, 0.32164453 ],
            [ 1.62757787,  -2.26431756,  -0.89705052, -0.34250171, 0.30845570 ],
            [ 1.90384657,  -2.50904433,  -0.87522412, -0.25892811, 0.29692179 ],
            [ 2.19511185,  -2.76106427,  -0.85586077, -0.17395409, 0.28689003 ],
            [ 2.50279663,  -3.02168824,  -0.83875580, -0.08705983, 0.27823675 ],
            [ 2.82898334,  -3.29274597,  -0.82369962, 0.00244569,  0.27083572 ],
            [ 3.17511643,  -3.57550313,  -0.81054159, 0.09506726,  0.26459005 ],
            [ 3.54295790,  -3.87148242,  -0.79913097, 0.19139952,  0.25940196 ],
            [ 3.93437435,  -4.18228651,  -0.78932714, 0.29206843,  0.25517779 ],
            [ 4.35196925,  -4.51009612,  -0.78098547, 0.39789237,  0.25182201 ],
            [ 4.79760295,  -4.85650566,  -0.77398712, 0.50950582,  0.24924922 ],
            [ 5.27378185,  -5.22362406,  -0.76820405, 0.62771475,  0.24736922 ],
            [ 5.78420016,  -5.61447201,  -0.76350180, 0.75362253,  0.24608988 ],
            [ 6.33253655,  -6.03204748,  -0.75975575, 0.88832834,  0.24532292 ],
            [ 6.92485862,  -6.48115327,  -0.75683452, 1.03351688,  0.24498045 ],
            [ 7.56742031,  -6.96670979,  -0.75461904, 1.19091558,  0.24497973 ],
            [ 8.27046607,  -7.49662841,  -0.75298859, 1.36322808,  0.24524243 ],
            [ 9.04823051,  -8.08179104,  -0.75183203, 1.55413877,  0.24569757 ],
            [ 9.92094017,  -8.73754957,  -0.75104867, 1.76881429,  0.24628161 ],
            [ 10.91941520, -9.48717699,  -0.75054841, 2.01505252,  0.24693896 ],
            [ 12.09307084, -10.36786737, -0.75025312, 2.30528602,  0.24762175 ],
            [ 13.52751898, -11.44393735, -0.75009717, 2.66098714,  0.24828894 ],
            [ 15.38997032, -12.84087899, -0.75002767, 3.12402582,  0.24890456 ],
            [ 18.08263900, -14.86041475, -0.75000446, 3.79503289,  0.24943464 ],
            [ 23.05336348, -18.58846448, -0.75000016, 5.03610029,  0.24983594 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 2.08e-07
        #    Est Max FD Error (based on N/2 run) to R=20.00: 7.e-7
        #    Est Max MOC error for R>20.00: 5.55e-07
        #    Max interpolation error with cubic splines: 5.0717e-07
        #    Est Max overall error after cubic interpolation: 1.8e-6
        #
        C8000_G1X1 => [
            [ -7.41085315, 12.00017707,  -1.99998713, -7.41212238, 0.99872997 ],
            [ -6.42603181, 10.03057725,  -1.99990492, -6.42943352, 0.99659265 ],
            [ -6.01541897, 9.20941537,   -1.99977291, -6.02055217, 0.99485417 ],
            [ -5.54243599, 8.26362876,   -1.99941326, -5.55068573, 0.99171847 ],
            [ -5.13137526, 7.44188093,   -1.99866787, -5.14384329, 0.98746194 ],
            [ -4.79780132, 6.77536431,   -1.99740910, -4.81524400, 0.98242623 ],
            [ -4.50341868, 6.18763563,   -1.99534974, -4.52689068, 0.97630353 ],
            [ -4.24889919, 5.68013577,   -1.99230426, -4.27925576, 0.96929237 ],
            [ -4.01912621, 5.22282573,   -1.98790932, -4.05743422, 0.96117708 ],
            [ -3.80939678, 4.80650023,   -1.98180679, -3.85678522, 0.95189926 ],
            [ -3.61411213, 4.42024527,   -1.97350996, -3.67189640, 0.94128213 ],
            [ -3.43013698, 4.05812744,   -1.96248224, -3.49980321, 0.92917858 ],
            [ -3.25344021, 3.71258055,   -1.94797569, -3.33681194, 0.91529380 ],
            [ -3.07920690, 3.37476641,   -1.92885762, -3.17870920, 0.89910574 ],
            [ -2.90161433, 3.03438459,   -1.90332069, -3.02070954, 0.87973718 ],
            [ -2.70509796, 2.66379688,   -1.86669252, -2.85022958, 0.85461679 ],
            [ -2.52033898, 2.32277037,   -1.82343278, -2.69479913, 0.82730520 ],
            [ -2.35103257, 2.01794925,   -1.77621425, -2.55707110, 0.79918499 ],
            [ -2.19224607, 1.73983368,   -1.72589975, -2.43243119, 0.77035093 ],
            [ -2.04040778, 1.48172218,   -1.67328429, -2.31767769, 0.74090108 ],
            [ -1.89293652, 1.23892418,   -1.61913937, -2.21060965, 0.71097455 ],
            [ -1.74770477, 1.00775557,   -1.56413017, -2.10954474, 0.68071095 ],
            [ -1.60301181, 0.78543620,   -1.50887128, -2.01325172, 0.65028301 ],
            [ -1.45755611, 0.56996934,   -1.45396493, -1.92087852, 0.61991514 ],
            [ -1.30988380, 0.35928134,   -1.39982465, -1.83157055, 0.58978015 ],
            [ -1.15891528, 0.15198155,   -1.34689992, -1.74478839, 0.56011513 ],
            [ -1.00338182, -0.05346505,  -1.29549282, -1.65994976, 0.53110894 ],
            [ -0.84223165, -0.25818771,  -1.24591965, -1.57665542, 0.50298188 ],
            [ -0.67417987, -0.46351241,  -1.19837836, -1.49443671, 0.47590133 ],
            [ -0.49793100, -0.67066307,  -1.15304279, -1.41287782, 0.45002537 ],
            [ -0.31206532, -0.88090376,  -1.11004134, -1.33155692, 0.42548527 ],
            [ -0.11504403, -1.09552815,  -1.06947406, -1.25004995, 0.40238969 ],
            [ 0.09499843,  -1.31607471,  -1.03138399, -1.16784987, 0.38080697 ],
            [ 0.31996300,  -1.54400497,  -0.99582752, -1.08448843, 0.36080427 ],
            [ 0.56101703,  -1.77999095,  -0.96296021, -0.99978220, 0.34249453 ],
            [ 0.81098854,  -2.01695978,  -0.93377114, -0.91623285, 0.32643345 ],
            [ 1.07089406,  -2.25619842,  -0.90789124, -0.83326822, 0.31240319 ],
            [ 1.34247823,  -2.49956342,  -0.88493230, -0.75013638, 0.30017353 ],
            [ 1.63030880,  -2.75123298,  -0.86439383, -0.66532958, 0.28945874 ],
            [ 1.93132401,  -3.00863538,  -0.84636884, -0.57962588, 0.28028501 ],
            [ 2.25091899,  -3.27650151,  -0.83040958, -0.49135428, 0.27239784 ],
            [ 2.59074064,  -3.55623080,  -0.81637485, -0.39996926, 0.26570365 ],
            [ 2.95085582,  -3.84794995,  -0.80418183, -0.30533007, 0.26013455 ],
            [ 3.33325561,  -4.15338733,  -0.79366921, -0.20676447, 0.25558363 ],
            [ 3.74105812,  -4.47514368,  -0.78467103, -0.10331648, 0.25194383 ],
            [ 4.17586175,  -4.81460360,  -0.77707382, 0.00558313,  0.24913060 ],
            [ 4.63997770,  -5.17372566,  -0.77075035, 0.12069445,  0.24705258 ],
            [ 5.13650429,  -5.55507802,  -0.76556879, 0.24297742,  0.24561681 ],
            [ 5.66892514,  -5.96151807,  -0.76140129, 0.37348823,  0.24473226 ],
            [ 6.24271186,  -6.39740594,  -0.75811494, 0.51376932,  0.24430878 ],
            [ 6.86312179,  -6.86691804,  -0.75559006, 0.66530868,  0.24426110 ],
            [ 7.53940273,  -7.37722853,  -0.75370186, 0.83056804,  0.24450767 ],
            [ 8.28299181,  -7.93712815,  -0.75233732, 1.01254598,  0.24497350 ],
            [ 9.11171809,  -8.56018427,  -0.75138976, 1.21581468,  0.24559162 ],
            [ 10.05020341, -9.26503110,  -0.75076518, 1.44663451,  0.24630192 ],
            [ 11.13846439, -10.08182595, -0.75037980, 1.71509347,  0.24705305 ],
            [ 12.44111365, -11.05914684, -0.75016287, 2.03742222,  0.24779986 ],
            [ 14.07736168, -12.28649657, -0.75005589, 2.44349304,  0.24850348 ],
            [ 16.29901962, -13.95280553, -0.75001312, 2.99633559,  0.24912778 ],
            [ 19.80350775, -16.58118986, -0.75000140, 3.87041184,  0.24963356 ],
            [ 27.62990345, -22.45098887, -0.75000035, 5.82574929,  0.24994812 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 2.32e-07
        #    Est Max FD Error (based on N/2 run) to R=10.00: 6.4e-8
        #    Est Max MOC error for R>10.00: -8.88e-07
        #    Max interpolation error with cubic splines: 5.079e-07
        #    Est Max overall error after cubic interpolation: 1.5e-6
        #
        S8000_G3 => [
            [ -4.16830247, 12.00017702,  -2.99997236, -4.16951744, 0.99817646 ],
            [ -3.53922615, 10.11299026,  -2.99985977, -3.54235061, 0.99530621 ],
            [ -3.20990726, 9.12511712,   -2.99960569, -3.21503262, 0.99229327 ],
            [ -2.91037377, 8.22671013,   -2.99902372, -2.91841723, 0.98789017 ],
            [ -2.66570809, 7.49306634,   -2.99796156, -2.67733670, 0.98246735 ],
            [ -2.41298021, 6.73565222,   -2.99567051, -2.43000828, 0.97427697 ],
            [ -2.25152444, 6.25218504,   -2.99297595, -2.27325971, 0.96711850 ],
            [ -2.08132223, 5.74313484,   -2.98834991, -2.10944644, 0.95738364 ],
            [ -1.92903604, 5.28851771,   -2.98172485, -1.96446587, 0.94623790 ],
            [ -1.79052280, 4.87610190,   -2.97256371, -1.83424481, 0.93358738 ],
            [ -1.66175139, 4.49407299,   -2.96013271, -1.71492203, 0.91919919 ],
            [ -1.54113309, 4.13796331,   -2.94370731, -1.60499893, 0.90298259 ],
            [ -1.42557475, 3.79897641,   -2.92214824, -1.50168744, 0.88455289 ],
            [ -1.31283868, 3.47105646,   -2.89402558, -1.40312409, 0.86347707 ],
            [ -1.19938183, 3.14472174,   -2.85694198, -1.30651960, 0.83886050 ],
            [ -1.07833671, 2.80186830,   -2.80581161, -1.20676951, 0.80858385 ],
            [ -0.95594788, 2.46232193,   -2.74046994, -1.10989917, 0.77370527 ],
            [ -0.84462665, 2.16113718,   -2.66868789, -1.02570133, 0.73847054 ],
            [ -0.74082570, 1.88803782,   -2.59179654, -0.95086769, 0.70302645 ],
            [ -0.64264640, 1.63746110,   -2.51161681, -0.88356564, 0.66774203 ],
            [ -0.54847531, 1.40476576,   -2.42970142, -0.82231982, 0.63287738 ],
            [ -0.45506202, 1.18171320,   -2.34560409, -0.76483321, 0.59791109 ],
            [ -0.36133883, 0.96586884,   -2.26042244, -0.71043299, 0.56304191 ],
            [ -0.26698287, 0.75659444,   -2.17575236, -0.65893464, 0.52870409 ],
            [ -0.17107512, 0.55194569,   -2.09244293, -0.60985265, 0.49506811 ],
            [ -0.07276095, 0.35026054,   -2.01121665, -0.56280675, 0.46229837 ],
            [ 0.02888957,  0.14986436,   -1.93258981, -0.51744816, 0.43051573 ],
            [ 0.13470608,  -0.05058085,  -1.85704776, -0.47353658, 0.39986161 ],
            [ 0.24556147,  -0.25238230,  -1.78495153, -0.43086527, 0.37045262 ],
            [ 0.36246351,  -0.45697263,  -1.71651467, -0.38922897, 0.34236405 ],
            [ 0.48630801,  -0.66547610,  -1.65197486, -0.34851208, 0.31569356 ],
            [ 0.61830477,  -0.87944408,  -1.59137459, -0.30854098, 0.29046604 ],
            [ 0.75966040,  -1.10029854,  -1.53476566, -0.26919802, 0.26671365 ],
            [ 0.91015499,  -1.32725330,  -1.48264801, -0.23075730, 0.24466080 ],
            [ 1.06614758,  -1.55481885,  -1.43617185, -0.19417636, 0.22481915 ],
            [ 1.22900388,  -1.78523075,  -1.39453213, -0.15905951, 0.20687463 ],
            [ 1.39984334,  -2.02019433,  -1.35713782, -0.12514141, 0.19059751 ],
            [ 1.58001134,  -2.26159283,  -1.32345919, -0.09216964, 0.17577899 ],
            [ 1.77068210,  -2.51096515,  -1.29309492, -0.05997442, 0.16226213 ],
            [ 1.97312010,  -2.76988765,  -1.26569520, -0.02840806, 0.14990976 ],
            [ 2.18872927,  -3.04004170,  -1.24095032, 0.00266378,  0.13860012 ],
            [ 2.41885565,  -3.32297399,  -1.21860581, 0.03333627,  0.12823457 ],
            [ 2.66497462,  -3.62034278,  -1.19843059, 0.06369687,  0.11872377 ],
            [ 2.92877813,  -3.93402070,  -1.18021025, 0.09383428,  0.10998475 ],
            [ 3.21214664,  -4.26605549,  -1.16375123, 0.12383168,  0.10194300 ],
            [ 3.51695251,  -4.61844225,  -1.14888857, 0.15374652,  0.09453664 ],
            [ 3.84566235,  -4.99382003,  -1.13545448, 0.18366898,  0.08770101 ],
            [ 4.20053390,  -5.39454262,  -1.12331599, 0.21364222,  0.08138773 ],
            [ 4.58441845,  -5.82359631,  -1.11233944, 0.24373536,  0.07554682 ],
            [ 5.00035960,  -6.28413638,  -1.10240767, 0.27400363,  0.07013558 ],
            [ 5.45179388,  -6.77970946,  -1.09341330, 0.30450240,  0.06511513 ],
            [ 5.94235081,  -7.31403113,  -1.08526232, 0.33527149,  0.06045252 ],
            [ 6.47645386,  -7.89163735,  -1.07786385, 0.36637062,  0.05611493 ],
            [ 7.05892076,  -8.51743936,  -1.07113879, 0.39784806,  0.05207455 ],
            [ 7.69476446,  -9.19651187,  -1.06501991, 0.42973055,  0.04830883 ],
            [ 8.38979402,  -9.93473595,  -1.05944460, 0.46205362,  0.04479620 ],
            [ 9.15001541,  -10.73815725, -1.05436018, 0.49483000,  0.04151938 ],
            [ 9.98243239,  -11.61383374, -1.04971730, 0.52808466,  0.03846130 ],
            [ 10.89424716, -12.56898896, -1.04547494, 0.56181869,  0.03560837 ],
            [ 11.89346090, -13.61164548, -1.04159594, 0.59603399,  0.03294778 ],
            [ 12.98847423, -14.75020606, -1.03804829, 0.63071806,  0.03046849 ],
            [ 14.18848740, -15.99387191, -1.03480298, 0.66585780,  0.02815995 ],
            [ 15.50350055, -17.35264054, -1.03183385, 0.70143773,  0.02601219 ],
            [ 16.94431370, -18.83730460, -1.02911739, 0.73743875,  0.02401583 ],
            [ 18.45007586, -20.38506711, -1.02673698, 0.77222955,  0.02224066 ],
            [ 20.18212331, -22.16137713, -1.02444223, 0.80920670,  0.02050528 ],
            [ 22.05632948, -24.07939121, -1.02236778, 0.84610676,  0.01891536 ],
            [ 24.07542792, -26.14170748, -1.02049612, 0.88279224,  0.01746277 ],
            [ 26.23948437, -28.34824290, -1.01881098, 0.91911239,  0.01613954 ],
            [ 28.69893120, -30.85192037, -1.01720499, 0.95719449,  0.01486425 ],
            [ 31.33289533, -33.52924973, -1.01576488, 0.99478217,  0.01370837 ],
            [ 34.13469912, -36.37336791, -1.01447698, 1.03168828,  0.01266429 ],
            [ 37.31698478, -39.59971084, -1.01324855, 1.07034504,  0.01165892 ],
            [ 40.69282789, -43.01838019, -1.01215513, 1.10813778,  0.01075588 ],
            [ 44.24158952, -46.60850878, -1.01118529, 1.14483669,  0.00994815 ],
            [ 48.26046879, -50.67042484, -1.01025879, 1.18321108,  0.00917035 ],
            [ 52.46128663, -54.91257420, -1.00944170, 1.22024378,  0.00847916 ],
            [ 57.21818094, -59.71247107, -1.00866091, 1.25895217,  0.00781387 ],
            [ 62.15340664, -64.68870575, -1.00797678, 1.29602884,  0.00722694 ],
            [ 67.73661891, -70.31457579, -1.00732264, 1.33475809,  0.00666206 ],
            [ 74.08516567, -76.70756971, -1.00669825, 1.37528041,  0.00611937 ],
            [ 80.64212037, -83.30661749, -1.00615623, 1.41380860,  0.00564538 ],
            [ 88.08544877, -90.79378519, -1.00563843, 1.45408593,  0.00518993 ],
            [ 95.67848133, -98.42789288, -1.00519296, 1.49195774,  0.00479596 ],
            [
                104.27165039, -107.06380276, -1.00476684, 1.53149751,
                0.00441714
            ],
            [
                114.04710498, -116.88380722, -1.00435992, 1.57284855,
                0.00405355
            ],
            [
                123.91000396, -126.78795594, -1.00401423, 1.61125667,
                0.00374318
            ],
            [
                135.07945489, -138.00034383, -1.00368354, 1.65135369,
                0.00344492
            ],
            [
                146.11157231, -149.07152985, -1.00340640, 1.68793866,
                0.00319391
            ],
            [
                158.52261288, -161.52315325, -1.00314060, 1.72603907,
                0.00295224
            ],
            [
                172.55075943, -175.59352110, -1.00288609, 1.76577913,
                0.00271992
            ],
            [
                188.48779585, -191.57455755, -1.00264279, 1.80729907,
                0.00249698
            ],
            [
                203.93845363, -207.06446545, -1.00244314, 1.84442428,
                0.00231337
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 2.52e-07
        #    Est Max FD Error (based on N/2 run) to R=24.00: 1.9e-7
        #    Est Max MOC error for R>24.00: 3.05e-07
        #    Max interpolation error with cubic splines: 5.0528e-07
        #    Est Max overall error after cubic interpolation: 1.e-6

        'P8000_G1X1' => [
            [
                -13.67225064, 12.00006664, -0.99999356, -13.67394340,
                0.99915291
            ],
            [
                -12.84419449, 11.17202136, -0.99998527, -12.84675657,
                0.99871734
            ],
            [
                -12.30156297, 10.62940384, -0.99997465, -12.30492494,
                0.99831624
            ],
            [
                -11.89666087, 10.22451846, -0.99994939, -11.90077879,
                0.99793690
            ],
            [ -10.68900900, 9.01698466, -0.99983145, -10.69655354, 0.99621414 ],
            [ -9.76711769,  8.09534826, -0.99957689, -9.77910505,  0.99397294 ],
            [ -9.01878784,  7.34748916, -0.99910720, -9.03625818,  0.99119637 ],
            [ -8.38680016,  6.71628731, -0.99832476, -8.41083137,  0.98786034 ],
            [ -7.83700223,  6.16771461, -0.99710897, -7.86873786,  0.98392692 ],
            [ -7.34834273,  5.68087124, -0.99531606, -7.38900275,  0.97935383 ],
            [ -6.90716069,  5.24227537, -0.99278161, -6.95804422,  0.97409966 ],
            [ -6.50169871,  4.84040092, -0.98930163, -6.56426008,  0.96808776 ],
            [ -6.12288682,  4.46647615, -0.98462657, -6.19879641,  0.96121607 ],
            [ -5.76299872,  4.11317506, -0.97844151, -5.85424033,  0.95334218 ],
            [ -5.41408419,  3.77313156, -0.97031459, -5.52314801,  0.94424117 ],
            [ -5.06640139,  3.43755105, -0.95958437, -5.19666388,  0.93352530 ],
            [ -4.70281811,  3.09120610, -0.94496544, -4.85958062,  0.92035259 ],
            [ -4.30922550,  2.72311260, -0.92464879, -4.50054012,  0.90363799 ],
            [ -3.95460024,  2.39908760, -0.90210444, -4.18309420,  0.88632524 ],
            [ -3.62606518,  2.10662528, -0.87777749, -3.89479622,  0.86844226 ],
            [ -3.32003713,  1.84180350, -0.85256186, -3.63176496,  0.85034474 ],
            [ -3.03013350,  1.59831992, -0.82697061, -3.38786474,  0.83213203 ],
            [ -2.74493282,  1.36618139, -0.80082030, -3.15318925,  0.81345111 ],
            [ -2.45648486,  1.13904497, -0.77407255, -2.92133767,  0.79407243 ],
            [ -2.16091693,  0.91426653, -0.74703743, -2.68959739,  0.77402311 ],
            [ -1.85970380,  0.69327844, -0.72049922, -2.45951736,  0.75371243 ],
            [ -1.54928765,  0.47366912, -0.69473941, -2.22874842,  0.73322197 ],
            [ -1.22603832,  0.25315289, -0.67001238, -1.99508460,  0.71264772 ],
            [ -0.88570190,  0.02919657, -0.64651954, -1.75607967,  0.69207771 ],
            [ -0.52348523, -0.20089452,  -0.62443965, -1.50915041, 0.67161405 ],
            [ -0.13365580, -0.44021489,  -0.60392390, -1.25134280, 0.65136798 ],
            [ 0.29115498,  -0.69264447,  -0.58508785, -0.97894427, 0.63144980 ],
            [ 0.71220257,  -0.93563393,  -0.56961894, -0.71684920, 0.61386671 ],
            [ 1.13939086,  -1.17612052,  -0.55671596, -0.45805121, 0.59810372 ],
            [ 1.57835134,  -1.41804384,  -0.54591797, -0.19869226, 0.58391636 ],
            [ 2.03416617,  -1.66474612,  -0.53688268, 0.06448105,  0.57113470 ],
            [ 2.51166229,  -1.91923599,  -0.52934620, 0.33437902,  0.55964045 ],
            [ 3.01523556,  -2.18416311,  -0.52310091, 0.61353622,  0.54935802 ],
            [ 3.54986666,  -2.46239630,  -0.51796610, 0.90472270,  0.54022324 ],
            [ 4.12076624,  -2.75685424,  -0.51378639, 1.21076490,  0.53218821 ],
            [ 4.73473907,  -3.07121982,  -0.51041792, 1.53528926,  0.52520055 ],
            [ 5.40071496,  -3.41020365,  -0.50773018, 1.88298209,  0.51920704 ],
            [ 6.12895863,  -3.77913985,  -0.50561069, 2.25917042,  0.51416519 ],
            [ 6.93552187,  -4.18623872,  -0.50395419, 2.67211554,  0.51001656 ],
            [ 7.84099586,  -4.64193955,  -0.50267517, 3.13233321,  0.50671032 ],
            [ 8.87232966,  -5.15983044,  -0.50170710, 3.65352761,  0.50419281 ],
            [ 10.05162083, -5.75104318,  -0.50100811, 4.24697514,  0.50241052 ],
            [ 11.39076584, -6.42162402,  -0.50054279, 4.91893019,  0.50126578 ],
            [ 12.91546915, -7.18456453,  -0.50026206, 5.68264282,  0.50059995 ],
            [ 14.71238785, -8.08333763,  -0.50010885, 6.58181654,  0.50024643 ],
            [ 16.96149921, -9.20804102,  -0.50003570, 7.70670569,  0.50008037 ],
            [ 20.01992708, -10.73731089, -0.50000777, 9.23604570,  0.50001747 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 2.02e-07
        #    Est Max FD Error (based on N/2 run) to R=10.00: 1.35e-6
        #    Est Max MOC error for R>10.00: 2.83e-07
        #    Max interpolation error with cubic splines: 1.0054e-06
        #    Est Max overall error after cubic interpolation: 2.6e-6
        S8000_G1X1 => [
            [ -5.03701586, 12.00017661,  -2.99998069, -5.03803112, 0.99847635 ],
            [ -4.48654608, 10.34879417,  -2.99990195, -4.48886590, 0.99651634 ],
            [ -4.29742007, 9.78144172,   -2.99979105, -4.30050197, 0.99537027 ],
            [ -4.08164398, 9.13416734,   -2.99965002, -4.08590614, 0.99359377 ],
            [ -3.81211982, 8.32573565,   -2.99926050, -3.81851192, 0.99038336 ],
            [ -3.57170920, 7.60477191,   -2.99843618, -3.58088848, 0.98617426 ],
            [ -3.29602641, 6.77839378,   -2.99643308, -3.30993530, 0.97901365 ],
            [ -3.06838095, 6.09662065,   -2.99296514, -3.08799629, 0.97035004 ],
            [ -2.87242446, 5.51062184,   -2.98741544, -2.89880864, 0.96005067 ],
            [ -2.69516782, 4.98178100,   -2.97880552, -2.72968142, 0.94766700 ],
            [ -2.53717090, 4.51204252,   -2.96645950, -2.58103245, 0.93343639 ],
            [ -2.39025951, 4.07743480,   -2.94898198, -2.44507678, 0.91681276 ],
            [ -2.24880722, 3.66193054,   -2.92435996, -2.31674184, 0.89704986 ],
            [ -2.11130175, 3.26201471,   -2.89051464, -2.19494931, 0.87367236 ],
            [ -1.96842148, 2.85229600,   -2.84214399, -2.07214500, 0.84442362 ],
            [ -1.83461119, 2.47582031,   -2.78239456, -1.96126299, 0.81205766 ],
            [ -1.71446186, 2.14539103,   -2.71587174, -1.86565271, 0.77881650 ],
            [ -1.60347343, 1.84786852,   -2.64382893, -1.78107390, 0.74477953 ],
            [ -1.49877667, 1.57499476,   -2.56760816, -1.70489324, 0.71010859 ],
            [ -1.39811102, 1.32047391,   -2.48831952, -1.63516809, 0.67491714 ],
            [ -1.29981975, 1.07986509,   -2.40703622, -1.57056957, 0.63936831 ],
            [ -1.20271127, 0.85010450,   -2.32484469, -1.51021147, 0.60369859 ],
            [ -1.10562883, 0.62840214,   -2.24257773, -1.45333389, 0.56809789 ],
            [ -1.00765058, 0.41269072,   -2.16104042, -1.39940909, 0.53280310 ],
            [ -0.90797477, 0.20130707,   -2.08095864, -1.34804443, 0.49806582 ],
            [ -0.80574147, -0.00740614,  -2.00286996, -1.29887785, 0.46409418 ],
            [ -0.70014891, -0.21485469,  -1.92725272, -1.25163484, 0.43109756 ],
            [ -0.59032420, -0.42246175,  -1.85445814, -1.20606232, 0.39924637 ],
            [ -0.47552365, -0.63129655,  -1.78486235, -1.16200745, 0.36873009 ],
            [ -0.35481445, -0.84268054,  -1.71867090, -1.11928311, 0.33966539 ],
            [ -0.22716279, -1.05799551,  -1.65601562, -1.07771534, 0.31213328 ],
            [ -0.09154296, -1.27850011,  -1.59702230, -1.03717881, 0.28620601 ],
            [ 0.05310602,  -1.50542269,  -1.54177896, -0.99757442, 0.26193338 ],
            [ 0.20819159,  -1.74043817,  -1.49023763, -0.95874957, 0.23929894 ],
            [ 0.37515303,  -1.98515280,  -1.44236333, -0.92059456, 0.21828576 ],
            [ 0.55574049,  -2.24151842,  -1.39805649, -0.88297726, 0.19884279 ],
            [ 0.75173998,  -2.51142488,  -1.35723891, -0.84580904, 0.18092468 ],
            [ 0.96527657,  -2.79712866,  -1.31978024, -0.80898350, 0.16446262 ],
            [ 1.19885984,  -3.10128499,  -1.28552914, -0.77238283, 0.14937400 ],
            [ 1.45530085,  -3.42682032,  -1.25433251, -0.73590108, 0.13557622 ],
            [ 1.73308827,  -3.77126130,  -1.22646116, -0.70001593, 0.12317577 ],
            [ 2.03305661,  -4.13532897,  -1.20172787, -0.66478330, 0.11208144 ],
            [ 2.35795520,  -4.52208374,  -1.17977833, -0.63003555, 0.10213061 ],
            [ 2.71106216,  -4.93511439,  -1.16029139, -0.59560297, 0.09317764 ],
            [ 3.09554010,  -5.37778419,  -1.14300924, -0.56137832, 0.08510816 ],
            [ 3.51488987,  -5.85378241,  -1.12769676, -0.52726421, 0.07782077 ],
            [ 3.97359637,  -6.36784140,  -1.11412285, -0.49312916, 0.07121794 ],
            [ 4.47612766,  -6.92459387,  -1.10209646, -0.45889331, 0.06522296 ],
            [ 5.02813774,  -7.52991304,  -1.09143028, -0.42444373, 0.05976214 ],
            [ 5.63586735,  -8.19022677,  -1.08196089, -0.38968729, 0.05477410 ],
            [ 6.30674613,  -8.91316429,  -1.07353825, -0.35452107, 0.05020400 ],
            [ 7.04899452,  -9.70711052,  -1.06603310, -0.31886321, 0.04600692 ],
            [ 7.87222589,  -10.58184370, -1.05932977, -0.28262854, 0.04214351 ],
            [ 8.78684841,  -11.54788863, -1.05333142, -0.24576231, 0.03858267 ],
            [ 9.80466672,  -12.61715334, -1.04795354, -0.20821631, 0.03529761 ],
            [ 10.93848322, -13.80250488, -1.04312526, -0.16996767, 0.03226662 ],
            [ 12.20229905, -15.11798287, -1.03878615, -0.13101095, 0.02947124 ],
            [ 13.61151468, -16.57900864, -1.03488401, -0.09135255, 0.02689517 ],
            [ 15.18293028, -18.20238483, -1.03137374, -0.05101196, 0.02452383 ],
            [ 16.86595752, -19.93556367, -1.02832792, -0.01155223, 0.02242181 ],
            [ 18.75198101, -21.87230742, -1.02555366, 0.02885761,  0.02046885 ],
            [ 20.77856513, -23.94813426, -1.02312464, 0.06854161,  0.01872721 ],
            [ 23.12635291, -26.34743080, -1.02083411, 0.11050651,  0.01705603 ],
            [ 25.63863790, -28.90946777, -1.01884013, 0.15146212,  0.01557701 ],
            [ 28.55347800, -31.87638615, -1.01695937, 0.19475694,  0.01416003 ],
            [ 31.65437642, -35.02726693, -1.01533290, 0.23669919,  0.01291641 ],
            [ 35.25474759, -38.67998161, -1.01379803, 0.28101355,  0.01172637 ],
            [ 39.05403745, -42.52910800, -1.01248084, 0.32356034,  0.01069159 ],
            [ 43.46366308, -46.99093029, -1.01123681, 0.36847632,  0.00970212 ],
            [ 48.06542747, -51.64188895, -1.01017878, 0.41112568,  0.00885074 ],
            [ 53.39658276, -57.02454959, -1.00917821, 0.45609339,  0.00803677 ],
            [ 58.87612750, -62.55200624, -1.00833640, 0.49820015,  0.00734492 ],
            [ 65.20022876, -68.92622079, -1.00753875, 0.54251396,  0.00668310 ],
            [ 72.54913252, -76.32765789, -1.00678470, 0.58925346,  0.00605152 ],
            [ 79.99991984, -83.82659518, -1.00616023, 0.63233839,  0.00552385 ],
            [ 88.60760448, -92.48467178, -1.00556832, 0.67766954,  0.00501961 ],
            [ 98.62158275, -102.55150958, -1.00500861, 0.72547245, 0.00453892 ],
            [
                108.56449968, -112.54189596, -1.00455420, 0.76861094,
                0.00414574
            ],
            [
                120.03056471, -124.05762525, -1.00412294, 0.81394723,
                0.00377002
            ],
            [
                133.34410523, -137.42324377, -1.00371458, 0.86170023,
                0.00341183
            ],
            [
                146.14604982, -150.27061230, -1.00339163, 0.90352222,
                0.00312678
            ],
            [
                160.82068915, -164.99269602, -1.00308428, 0.94736666,
                0.00285398
            ],
            [
                177.74849410, -181.97015901, -1.00279243, 0.99342672,
                0.00259347
            ],
            [
                197.41134579, -201.68510480, -1.00251592, 1.04192479,
                0.00234528
            ],
            [
                215.52011054, -219.83747073, -1.00230568, 1.08265052,
                0.00215562
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 3e-07
        #    Est Max FD Error (based on N/2 run) to R=16.69: 8.8e-8
        #    Est Max MOC error for R>16.69: -9.8e-07
        #    Max interpolation error with cubic splines: 5.0475e-07
        #    Est Max overall error after cubic interpolation: 1.6e-6
        C8000_G3 => [
            [ -6.07223275, 12.00017755,  -1.99998157, -6.07375169, 0.99847992 ],
            [ -5.54796515, 10.95165727,  -1.99994742, -5.55053230, 0.99742962 ],
            [ -5.24406463, 10.34387526,  -1.99992188, -5.24754501, 0.99651372 ],
            [ -4.74790279, 9.35162225,   -1.99977973, -4.75362517, 0.99426191 ],
            [ -4.50032576, 8.85654002,   -1.99964553, -4.50766128, 0.99263897 ],
            [ -4.24502244, 8.34605629,   -1.99938593, -4.25450104, 0.99047950 ],
            [ -3.84538318, 7.54715344,   -1.99863540, -3.85954870, 0.98574418 ],
            [ -3.50766616, 6.87237708,   -1.99732423, -3.52757244, 0.97992333 ],
            [ -3.22056190, 6.29920493,   -1.99526235, -3.24716196, 0.97311203 ],
            [ -2.96626636, 5.79218068,   -1.99215725, -3.00067285, 0.96514332 ],
            [ -2.73817301, 5.33825332,   -1.98770427, -2.78153413, 0.95598034 ],
            [ -2.53057110, 4.92619939,   -1.98154400, -2.58411526, 0.94554774 ],
            [ -2.33886967, 4.54708201,   -1.97325102, -2.40394772, 0.93373733 ],
            [ -2.15802519, 4.19117295,   -1.96222495, -2.23626511, 0.92029876 ],
            [ -1.98625533, 3.85529165,   -1.94788780, -2.07945314, 0.90511821 ],
            [ -1.81795670, 3.52897037,   -1.92912399, -1.92855562, 0.88763781 ],
            [ -1.64853831, 3.20414678,   -1.90437182, -1.77986650, 0.86714626 ],
            [ -1.46935702, 2.86582297,   -1.87055246, -1.62668273, 0.84207744 ],
            [ -1.28641254, 2.52746453,   -1.82687662, -1.47525170, 0.81279234 ],
            [ -1.12010343, 2.22752609,   -1.77880720, -1.34250923, 0.78307531 ],
            [ -0.96512984, 1.95577582,   -1.72721883, -1.22345305, 0.75305737 ],
            [ -0.81786459, 1.70535244,   -1.67302392, -1.11475835, 0.72289498 ],
            [ -0.67735869, 1.47412939,   -1.61780947, -1.01527069, 0.69312146 ],
            [ -0.53839754, 1.25323649,   -1.56117310, -0.92102828, 0.66323357 ],
            [ -0.39918641, 1.03989865,   -1.50378297, -0.83078158, 0.63336418 ],
            [ -0.25914460, 0.83331871,   -1.44668174, -0.74415829, 0.60387169 ],
            [ -0.11699112, 0.63168882,   -1.39048654, -0.66038698, 0.57493042 ],
            [ 0.02866563,  0.43318686,   -1.33564812, -0.57872155, 0.54666970 ],
            [ 0.17916022,  0.23622304,   -1.28254625, -0.49853983, 0.51921933 ],
            [ 0.33581995,  0.03935632,   -1.23149220, -0.41930403, 0.49270539 ],
            [ 0.49986353,  -0.15859999,  -1.18277360, -0.34059789, 0.46726761 ],
            [ 0.67286512,  -0.35914750,  -1.13653152, -0.26189704, 0.44298898 ],
            [ 0.85629496,  -0.56353964,  -1.09292027, -0.18279205, 0.41997305 ],
            [ 1.05181477,  -0.77314098,  -1.05203369, -0.10284361, 0.39830034 ],
            [ 1.26147668,  -0.98961750,  -1.01389247, -0.02151193, 0.37801698 ],
            [ 1.48472051,  -1.21195855,  -0.97891739, 0.06074510,  0.35938355 ],
            [ 1.71570232,  -1.43440097,  -0.94795734, 0.14180058,  0.34288702 ],
            [ 1.95654779,  -1.65931037,  -0.92045233, 0.22257304,  0.32825617 ],
            [ 2.20888480,  -1.88840501,  -0.89601120, 0.30372308,  0.31530231 ],
            [ 2.47457890,  -2.12350199,  -0.87429458, 0.38593136,  0.30385891 ],
            [ 2.75519173,  -2.36605970,  -0.85504902, 0.46974263,  0.29380077 ],
            [ 3.05233630,  -2.61753031,  -0.83805604, 0.55569600,  0.28501761 ],
            [ 3.36741446,  -2.87915757,  -0.82313575, 0.64425970,  0.27741656 ],
            [ 3.70274495,  -3.15291976,  -0.81009177, 0.73615181,  0.27089377 ],
            [ 4.06007304,  -3.44029613,  -0.79877940, 0.83192368,  0.26536958 ],
            [ 4.44113522,  -3.74276109,  -0.78906388, 0.93213154,  0.26076670 ],
            [ 4.84826304,  -4.06226379,  -0.78080473, 1.03749508,  0.25700283 ],
            [ 5.28397615,  -4.40089683,  -0.77386997, 1.14878633,  0.25399782 ],
            [ 5.75078195,  -4.76074554,  -0.76813513, 1.26678019,  0.25167303 ],
            [ 6.25237899,  -5.14481311,  -0.76346818, 1.39255680,  0.24994539 ],
            [ 6.79285219,  -5.55638839,  -0.75974417, 1.52729269,  0.24873424 ],
            [ 7.37767496,  -5.99980595,  -0.75683812, 1.67250911,  0.24795920 ],
            [ 8.01430894,  -6.48088717,  -0.75462805, 1.83021675,  0.24754247 ],
            [ 8.71220411,  -7.00693048,  -0.75299942, 2.00291422,  0.24741071 ],
            [ 9.48540083,  -7.58866478,  -0.75184255, 2.19423292,  0.24749567 ],
            [ 10.35453055, -8.24174158,  -0.75105724, 2.40943622,  0.24773609 ],
            [ 11.35001721, -8.98913275,  -0.75055468, 2.65622216,  0.24807793 ],
            [ 12.52087901, -9.86773269,  -0.75025725, 2.94692295,  0.24847456 ],
            [ 13.95193001, -10.94125883, -0.75009956, 3.30280874,  0.24888640 ],
            [ 15.80878237, -12.33400408, -0.75002880, 3.76534279,  0.24927988 ],
            [ 18.48625101, -14.34214144, -0.75000482, 4.43329193,  0.24962513 ],
            [ 23.37397346, -18.00794031, -0.75000020, 5.65416041,  0.24988873 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 2.49e-07
        #    Est Max FD Error (based on N/2 run) to R=42.61: 1.1e-7
        #    Est Max MOC error for R>42.61: 1.9e-07
        #    Max interpolation error with cubic splines: 5.0862e-07
        #    Est Max overall error after cubic interpolation: 8.e-7
        #
        C8000_G2 => [
            [ -6.37381146, 12.00017718,  -1.99998362, -6.37524347, 0.99856698 ],
            [ -5.77060139, 10.79377453,  -1.99994526, -5.77322060, 0.99737743 ],
            [ -4.96034929, 9.17337044,   -1.99974978, -4.96624785, 0.99408479 ],
            [ -4.46967424, 8.19222894,   -1.99933290, -4.47932607, 0.99030484 ],
            [ -4.07507385, 7.40342872,   -1.99853308, -4.08942575, 0.98555578 ],
            [ -3.74463522, 6.74323763,   -1.99716516, -3.76465630, 0.97980736 ],
            [ -3.45965849, 6.17437125,   -1.99500372, -3.48635450, 0.97301576 ],
            [ -3.20838628, 5.67345397,   -1.99178145, -3.24281091, 0.96512801 ],
            [ -2.98264440, 5.22430700,   -1.98718094, -3.02592384, 0.95607033 ],
            [ -2.77649817, 4.81526941,   -1.98082561, -2.82986025, 0.94574593 ],
            [ -2.58524463, 4.43720029,   -1.97225748, -2.65006828, 0.93401789 ],
            [ -2.40518861, 4.08304867,   -1.96092290, -2.48305475, 0.92071351 ],
            [ -2.23260064, 3.74583369,   -1.94608479, -2.32542256, 0.90555318 ],
            [ -2.06351818, 3.41834811,   -1.92670091, -2.17374693, 0.88809418 ],
            [ -1.89253214, 3.09100591,   -1.90106934, -2.02361012, 0.86752631 ],
            [ -1.70777316, 2.74293259,   -1.86530196, -1.86564871, 0.84176699 ],
            [ -1.52668829, 2.40900978,   -1.82118257, -1.71577920, 0.81287845 ],
            [ -1.36132391, 2.11174293,   -1.77285823, -1.58375022, 0.78349238 ],
            [ -1.20668624, 1.84151009,   -1.72122265, -1.46486787, 0.75373621 ],
            [ -1.05997278, 1.59287936,   -1.66743768, -1.35645735, 0.72389740 ],
            [ -0.91866016, 1.36110658,   -1.61243437, -1.25625433, 0.69415635 ],
            [ -0.77788946, 1.13809186,   -1.55586680, -1.16065315, 0.66406334 ],
            [ -0.63738357, 0.92348229,   -1.49898216, -1.06945867, 0.63407316 ],
            [ -0.49593293, 0.71545951,   -1.44252317, -0.98187547, 0.60440931 ],
            [ -0.35208397, 0.51197625,   -1.38699158, -0.89704563, 0.57521492 ],
            [ -0.20455568, 0.31138899,   -1.33284141, -0.81431108, 0.54665171 ],
            [ -0.05204554, 0.11216263,   -1.28042723, -0.73308373, 0.51886692 ],
            [ 0.10671230,  -0.08705862,  -1.23004750, -0.65287011, 0.49200908 ],
            [ 0.27299538,  -0.28752952,  -1.18194340, -0.57323479, 0.46621988 ],
            [ 0.44812450,  -0.49045241,  -1.13630300, -0.49377728, 0.44162911 ],
            [ 0.63375735,  -0.69730696,  -1.09320848, -0.41400246, 0.41831921 ],
            [ 0.83143647,  -0.90932665,  -1.05276873, -0.33352382, 0.39639132 ],
            [ 1.04319268,  -1.12816295,  -1.01499761, -0.25180713, 0.37589537 ],
            [ 1.27054667,  -1.35484777,  -0.98000603, -0.16855728, 0.35693136 ],
            [ 1.50565985,  -1.58152069,  -0.94900515, -0.08665922, 0.34018706 ],
            [ 1.75063856,  -1.81054033,  -0.92144436, -0.00518453, 0.32538373 ],
            [ 2.00727252,  -2.04378071,  -0.89692172, 0.07659433,  0.31231601 ],
            [ 2.27730161,  -2.28294776,  -0.87511666, 0.15932895,  0.30081748 ],
            [ 2.56225543,  -2.52947936,  -0.85578009, 0.24356871,  0.29075626 ],
            [ 2.86380761,  -2.78488633,  -0.83869135, 0.32988415,  0.28201283 ],
            [ 3.18362655,  -3.05063593,  -0.82366252, 0.41883002,  0.27448257 ],
            [ 3.52347423,  -3.32824919,  -0.81052572, 0.51098100,  0.26806931 ],
            [ 3.88518448,  -3.61929164,  -0.79912959, 0.60693046,  0.26268343 ],
            [ 4.27068258,  -3.92539705,  -0.78933423, 0.70729972,  0.25823910 ],
            [ 4.68206754,  -4.24833743,  -0.78100591, 0.81276151,  0.25465173 ],
            [ 5.12165964,  -4.59006107,  -0.77401482, 0.92405206,  0.25183692 ],
            [ 5.59216720,  -4.95282000,  -0.76823211, 1.04201172,  0.24970942 ],
            [ 6.09708645,  -5.33947057,  -0.76352760, 1.16768132,  0.24818289 ],
            [ 6.64029991,  -5.75315717,  -0.75977622, 1.30219839,  0.24717279 ],
            [ 7.22727849,  -6.19822218,  -0.75685089, 1.44709254,  0.24659492 ],
            [ 7.86388025,  -6.67928454,  -0.75463261, 1.60398576,  0.24636875 ],
            [ 8.56215412,  -7.20561347,  -0.75299590, 1.77602199,  0.24641748 ],
            [ 9.33513707,  -7.78718232,  -0.75183458, 1.96658752,  0.24667208 ],
            [ 10.20265827, -8.43904298,  -0.75104825, 2.18074863,  0.24706984 ],
            [ 11.19554006, -9.18446982,  -0.75054642, 2.42630178,  0.24755599 ],
            [ 12.36259929, -10.06020774, -0.75025079, 2.71552915,  0.24808325 ],
            [ 13.78884905, -11.13012455, -0.75009529, 3.06975071,  0.24861153 ],
            [ 15.64150089, -12.51971326, -0.75002655, 3.53083053,  0.24910664 ],
            [ 18.32596960, -14.53309695, -0.75000403, 4.20018889,  0.24953758 ],
            [ 23.34909541, -18.30044680, -0.75000011, 5.45464354,  0.24986758 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 2.21e-07
        #    Est Max FD Error (based on N/2 run) to R=10.00: 9.7e-8
        #    Est Max MOC error for R>10.00: -2.3e-07
        #    Max interpolation error with cubic splines: 5.0701e-07
        #    Est Max overall error after cubic interpolation: 8.3e-7
        #
        S8000_G2 => [
            [ -4.35780834, 12.00017668,  -2.99997543, -4.35895378, 0.99828087 ],
            [ -4.01113331, 10.96016446,  -2.99993048, -4.01306074, 0.99710613 ],
            [ -3.62605081, 9.80495970,   -2.99980493, -3.62948762, 0.99483623 ],
            [ -3.19953751, 8.52558181,   -2.99932630, -3.20606340, 0.99018134 ],
            [ -2.91220930, 7.66390422,   -2.99840635, -2.92226743, 0.98484465 ],
            [ -2.67526991, 6.95363573,   -2.99676138, -2.68964754, 0.97830100 ],
            [ -2.47326676, 6.34852441,   -2.99407960, -2.49277410, 0.97050953 ],
            [ -2.29676240, 5.82038611,   -2.98998867, -2.32224102, 0.96141921 ],
            [ -2.13932535, 5.35008341,   -2.98404295, -2.17166888, 0.95095229 ],
            [ -1.99655705, 4.92461010,   -2.97572164, -2.03672405, 0.93901876 ],
            [ -1.86498626, 4.53379118,   -2.96439525, -1.91403765, 0.92548300 ],
            [ -1.74181361, 4.16953664,   -2.94929207, -1.80095868, 0.91015671 ],
            [ -1.62457564, 3.82487188,   -2.92943086, -1.69524431, 0.89276918 ],
            [ -1.51067585, 3.49261958,   -2.90344907, -1.59466100, 0.87288058 ],
            [ -1.39702109, 3.16448443,   -2.86928041, -1.49673587, 0.84975059 ],
            [ -1.27790244, 2.82534970,   -2.82282818, -1.39714562, 0.82170117 ],
            [ -1.15123253, 2.47162002,   -2.75978328, -1.29518207, 0.78744932 ],
            [ -1.03679393, 2.15966972,   -2.69003396, -1.20702281, 0.75269976 ],
            [ -0.93052953, 1.87772691,   -2.61481863, -1.12888296, 0.71754578 ],
            [ -0.82984594, 1.61839059,   -2.53553458, -1.05840433, 0.68216971 ],
            [ -0.73421917, 1.37975418,   -2.45473020, -0.99483161, 0.64727029 ],
            [ -0.63984703, 1.15199677,   -2.37168539, -0.93540027, 0.61218393 ],
            [ -0.54534001, 0.93184338,   -2.28725215, -0.87920875, 0.57700892 ],
            [ -0.45047104, 0.71885667,   -2.20312162, -0.82612437, 0.54223453 ],
            [ -0.35421150, 0.51080716,   -2.12006405, -0.77558659, 0.50801340 ],
            [ -0.25578536, 0.30616814,   -2.03888193, -0.72724534, 0.47456042 ],
            [ -0.15427885, 0.10325196,   -1.96010733, -0.68074375, 0.44202120 ],
            [ -0.04891025, -0.09923121,  -1.88425778, -0.63584702, 0.41056931 ],
            [ 0.06125039,  -0.30274013,  -1.81164981, -0.59230937, 0.38032043 ],
            [ 0.17703268,  -0.50842947,  -1.74260272, -0.54997678, 0.35140863 ],
            [ 0.29942151,  -0.71763021,  -1.67728838, -0.50868178, 0.32391676 ],
            [ 0.42942464,  -0.93160460,  -1.61583925, -0.46829578, 0.29791588 ],
            [ 0.56825189,  -1.15184134,  -1.55828047, -0.42867337, 0.27343257 ],
            [ 0.71722182,  -1.37988374,  -1.50460517, -0.38968963, 0.25047805 ],
            [ 0.87320319,  -1.61069662,  -1.45609485, -0.35228566, 0.22961661 ],
            [ 1.03557836,  -1.84351443,  -1.41265647, -0.31656415, 0.21082719 ],
            [ 1.20566051,  -2.08037645,  -1.37361306, -0.28218698, 0.19383311 ],
            [ 1.38464933,  -2.32300797,  -1.33844094, -0.24890684, 0.17841908 ],
            [ 1.57374073,  -2.57301604,  -1.30670887, -0.21652746, 0.16440667 ],
            [ 1.77426530,  -2.83209150,  -1.27803916, -0.18487255, 0.15163849 ],
            [ 1.98745291,  -3.10171665,  -1.25213092, -0.15381859, 0.13998945 ],
            [ 2.21459552,  -3.38339658,  -1.22871915, -0.12326051, 0.12934947 ],
            [ 2.45717424,  -3.67881936,  -1.20755879, -0.09309410, 0.11961666 ],
            [ 2.71687029,  -3.98986301,  -1.18842688, -0.06321823, 0.11069843 ],
            [ 2.99531298,  -4.31829579,  -1.17113758, -0.03356350, 0.10251887 ],
            [ 3.29437885,  -4.66613941,  -1.15551537, -0.00405626, 0.09500661 ],
            [ 3.61633693,  -5.03582692,  -1.14139241, 0.02538989,  0.08809340 ],
            [ 3.96344727,  -5.42973528,  -1.12862768, 0.05483355,  0.08172354 ],
            [ 4.33816235,  -5.85042362,  -1.11709248, 0.08432678,  0.07584716 ],
            [ 4.74332726,  -6.30085457,  -1.10666457, 0.11392859,  0.07041728 ],
            [ 5.18217971,  -6.78438480,  -1.09723106, 0.14369963,  0.06539104 ],
            [ 5.65835423,  -7.30476315,  -1.08868920, 0.17369844,  0.06073017 ],
            [ 6.17586712,  -7.86610980,  -1.08094684, 0.20397760,  0.05640118 ],
            [ 6.73914087,  -8.47294015,  -1.07392140, 0.23458322,  0.05237482 ],
            [ 7.35318913,  -9.13036026,  -1.06753756, 0.26556229,  0.04862470 ],
            [ 8.02342155,  -9.84385315,  -1.06172951, 0.29694971,  0.04512854 ],
            [ 8.75584465,  -10.61949197, -1.05643835, 0.32877682,  0.04186659 ],
            [ 9.55726257,  -11.46414766, -1.05161090, 0.36107669,  0.03882088 ],
            [ 10.43467781, -12.38485486, -1.04720298, 0.39385805,  0.03597731 ],
            [ 11.39589177, -13.38944649, -1.04317462, 0.42713005,  0.03332275 ],
            [ 12.44930517, -14.48634190, -1.03949092, 0.46089294,  0.03084566 ],
            [ 13.60331836, -15.68392626, -1.03612287, 0.49512110,  0.02853687 ],
            [ 14.86813150, -16.99242068, -1.03304132, 0.52981755,  0.02638562 ],
            [ 16.25414464, -18.42221695, -1.03022193, 0.56496283,  0.02438311 ],
            [ 17.77235779, -19.98429700, -1.02764302, 0.60052951,  0.02252121 ],
            [ 19.36628238, -21.62042176, -1.02537363, 0.63506581,  0.02085765 ],
            [ 21.15253047, -23.45003047, -1.02323879, 0.67086642,  0.01927019 ],
            [ 23.05956571, -25.39950620, -1.02132576, 0.70620375,  0.01782821 ],
            [ 25.22524864, -27.60934292, -1.01950462, 0.74326622,  0.01643750 ],
            [ 27.53291554, -29.96009015, -1.01787946, 0.77970926,  0.01518096 ],
            [ 29.97241574, -32.44139368, -1.01643342, 0.81532726,  0.01404998 ],
            [ 32.73708992, -35.24953923, -1.01505475, 0.85262242,  0.01295983 ],
            [ 35.64709726, -38.20152971, -1.01383411, 0.88887872,  0.01198451 ],
            [ 38.94675308, -41.54485648, -1.01267012, 0.92683068,  0.01104516 ],
            [ 42.40004628, -45.04008392, -1.01164532, 0.96349282,  0.01021029 ],
            [ 46.31521870, -48.99887979, -1.01066766, 1.00185083,  0.00940662 ],
            [ 50.38201442, -53.10727106, -1.00981241, 1.03862555,  0.00869757 ],
            [ 54.98818325, -57.75670251, -1.00899593, 1.07707285,  0.00801515 ],
            [ 59.72602317, -62.53543816, -1.00828696, 1.11359578,  0.00741804 ],
            [ 65.08111849, -67.93304503, -1.00760939, 1.15173903,  0.00684323 ],
            [ 71.16411895, -74.06030888, -1.00696292, 1.19163858,  0.00629086 ],
            [ 77.37332821, -80.31097995, -1.00640744, 1.22917669,  0.00581300 ],
            [ 84.40455099, -87.38533511, -1.00587669, 1.26838892,  0.00535351 ],
            [ 91.46519704, -94.48583913, -1.00542551, 1.30476524,  0.00496056 ],
            [ 99.42219530, -102.48424206, -1.00499356, 1.34269100, 0.00458224 ],
            [
                108.43236162, -111.53748753, -1.00458067, 1.38229393,
                0.00421861
            ],
            [
                117.33047972, -120.47478649, -1.00423493, 1.41843530,
                0.00391253
            ],
            [
                127.33897586, -130.52396469, -1.00390358, 1.45607987,
                0.00361777
            ],
            [
                138.64890262, -141.87619698, -1.00358651, 1.49534991,
                0.00333436
            ],
            [
                151.49409615, -154.76545752, -1.00328364, 1.53638338,
                0.00306235
            ],
            [
                163.94319780, -167.25385087, -1.00303525, 1.57307646,
                0.00283828
            ],
            [
                177.95355597, -181.30502062, -1.00279715, 1.61129081,
                0.00262261
            ],
            [
                193.79558501, -197.18950605, -1.00256927, 1.65115135,
                0.00241537
            ],
            [
                211.80070675, -215.23887004, -1.00235156, 1.69279912,
                0.00221658
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 4.79e-07
        #    Est Max FD Error (based on N/2 run) to R=57.86: 3.3e-7
        #    Est Max MOC error for R>57.86: 5.03e-07
        #    Max interpolation error with cubic splines: 5.0444e-07
        #    Est Max overall error after cubic interpolation: 1.3e-6
        P8000_G2 => [
            [
                -11.51920152, 12.00017799, -0.99999181, -11.52111132,
                0.99904420
            ],
            [ -9.73385547, 10.21487281, -0.99995049, -9.73852484, 0.99765999 ],
            [ -8.53341803, 9.01455035,  -0.99983551, -8.54194386, 0.99571970 ],
            [ -7.61074969, 8.09213123,  -0.99958656, -7.62430524, 0.99317944 ],
            [ -6.85867970, 7.34052453,  -0.99912417, -6.87847943, 0.99001192 ],
            [ -6.22377890, 6.70639955,  -0.99835125, -6.25106519, 0.98619621 ],
            [ -5.67292746, 6.15675735,  -0.99715001, -5.70899732, 0.98169836 ],
            [ -5.18487872, 5.67049693,  -0.99538127, -5.23110066, 0.97647719 ],
            [ -4.74432255, 5.23248804,  -0.99287778, -4.80217964, 0.97047225 ],
            [ -4.34063562, 4.83232504,  -0.98944486, -4.41174813, 0.96361596 ],
            [ -3.96506930, 4.46153650,  -0.98484602, -4.05126465, 0.95581194 ],
            [ -3.60981231, 4.11268347,  -0.97878018, -3.71323856, 0.94691660 ],
            [ -3.26767827, 3.77909949,  -0.97085228, -3.39096094, 0.93672445 ],
            [ -2.93023099, 3.45315879,  -0.96048019, -3.07680978, 0.92489104 ],
            [ -2.58465998, 3.12353349,  -0.94664355, -2.75957750, 0.91074126 ],
            [ -2.19541442, 2.75881108,  -0.92650049, -2.40859761, 0.89217123 ],
            [ -1.84359439, 2.43672052,  -0.90376897, -2.09802665, 0.87296280 ],
            [ -1.51976450, 2.14795641,  -0.87908552, -1.81845985, 0.85337949 ],
            [ -1.21757080, 1.88615482,  -0.85317138, -1.56351726, 0.83369887 ],
            [ -0.93588977, 1.64947236,  -0.82706989, -1.33137770, 0.81442219 ],
            [ -0.66158969, 1.42623311,  -0.80049899, -1.11062285, 0.79510117 ],
            [ -0.38615768, 1.20948349,  -0.77337365, -0.89432675, 0.77548597 ],
            [ -0.10128767, 0.99314606,  -0.74557866, -0.67629329, 0.75532041 ],
            [ 0.18965779,  0.78024452,  -0.71814862, -0.45947971, 0.73519126 ],
            [ 0.48920374,  0.56916359,  -0.69149906, -0.24226406, 0.71526501 ],
            [ 0.80125521,  0.35743738,  -0.66588822, -0.02216167, 0.69562152 ],
            [ 1.12985733,  0.14269881,  -0.64155647, 0.20321374,  0.67635406 ],
            [ 1.47971718,  -0.07766714, -0.61870297, 0.43650121,  0.65754912 ],
            [ 1.85643678,  -0.30664301, -0.59749594, 0.68070955,  0.63929061 ],
            [ 2.26597482,  -0.54724879, -0.57811843, 0.93884319,  0.62169802 ],
            [ 2.67158986,  -0.77842233, -0.56226210, 1.18784682,  0.60642522 ],
            [ 3.08321852,  -1.00706047, -0.54908239, 1.43462114,  0.59290211 ],
            [ 3.50649335,  -1.23706575, -0.53810304, 1.68296636,  0.58083551 ],
            [ 3.94639256,  -1.47169264, -0.52897786, 1.93603648,  0.57001910 ],
            [ 4.40950657,  -1.71484618, -0.52141598, 2.19770228,  0.56026670 ],
            [ 4.89044485,  -1.96408871, -0.51533586, 2.46502722,  0.55165003 ],
            [ 5.40426685,  -2.22755216, -0.51040580, 2.74642218,  0.54387695 ],
            [ 5.95103473,  -2.50550903, -0.50652515, 3.04184795,  0.53696082 ],
            [ 6.53645446,  -2.80111498, -0.50354207, 3.35434349,  0.53083321 ],
            [ 7.17265734,  -3.12071129, -0.50130440, 3.69027049,  0.52539459 ],
            [ 7.86169021,  -3.46553778, -0.49971537, 4.05059112,  0.52065720 ],
            [ 8.62160762,  -3.84483759, -0.49864736, 4.44461547,  0.51653154 ],
            [ 9.46952921,  -4.26734688, -0.49800477, 4.88102019,  0.51298175 ],
            [ 10.43457354, -4.74776977, -0.49770204, 5.37453695,  0.50995877 ],
            [ 11.55059966, -5.30317766, -0.49766986, 5.94217678,  0.50743933 ],
            [ 12.86906376, -5.95943817, -0.49784750, 6.60977546,  0.50538748 ],
            [ 14.47221058, -6.75781955, -0.49817969, 7.41857759,  0.50375309 ],
            [ 16.59848679, -7.81757360, -0.49862878, 8.48816121,  0.50242761 ],
            [ 19.01901911, -9.02504878, -0.49904669, 9.70313859,  0.50154336 ],
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
        #    Energy error to P=0.1 P0: 3.98e-07
        #    Est Max FD Error (based on N/2 run) to R=100.04: 2.7e-7
        #    Est Max MOC error for R>100.04: -1.88e-06
        #    Max interpolation error with cubic splines: 5.0321e-07
        #    Est Max overall error after cubic interpolation: 2.7e-6
        P8000_G3 => [
            [
                -10.85863198, 12.00017826, -0.99999079, -10.86065774,
                0.99898610
            ],
            [ -9.09234545, 10.23393348, -0.99994949, -9.09725166, 0.99754102 ],
            [ -7.89815703, 9.03986173,  -0.99983293, -7.90708808, 0.99551540 ],
            [ -6.98361242, 8.12556728,  -0.99958303, -6.99775589, 0.99288170 ],
            [ -6.23872881, 7.38114484,  -0.99912252, -6.25931501, 0.98961152 ],
            [ -5.60445291, 6.74764523,  -0.99834905, -5.63281778, 0.98564393 ],
            [ -5.05391392, 6.19831624,  -0.99714692, -5.09140906, 0.98096411 ],
            [ -4.56586518, 5.71205785,  -0.99537584, -4.61392122, 0.97552658 ],
            [ -4.12593391, 5.27467256,  -0.99287287, -4.18607789, 0.96927998 ],
            [ -3.72255944, 4.87482145,  -0.98943893, -3.79648427, 0.96214199 ],
            [ -3.34730557, 4.50434387,  -0.98483848, -3.43691218, 0.95401632 ],
            [ -2.99267348, 4.15610667,  -0.97877550, -3.10017910, 0.94476291 ],
            [ -2.65147679, 3.82343637,  -0.97085898, -2.77958520, 0.93417388 ],
            [ -2.31559177, 3.49899674,  -0.96052255, -2.46781565, 0.92191125 ],
            [ -1.97232430, 3.17154127,  -0.94676484, -2.15379939, 0.90729577 ],
            [ -1.58912869, 2.81238453,  -0.92693591, -1.80967053, 0.88833184 ],
            [ -1.23756770, 2.49037598,  -0.90420626, -1.50079217, 0.86846580 ],
            [ -0.91441613, 2.20207798,  -0.87949650, -1.22336200, 0.84827535 ],
            [ -0.61159753, 1.93963898,  -0.85337990, -0.96953629, 0.82794631 ],
            [ -0.33042915, 1.70334991,  -0.82711715, -0.73950812, 0.80816882 ],
            [ -0.05754453, 1.48126524,  -0.80042489, -0.52165027, 0.78848263 ],
            [ 0.21621611,  1.26586954,  -0.77315859, -0.30851751, 0.76859760 ],
            [ 0.49857725,  1.05151250,  -0.74525544, -0.09436939, 0.74830699 ],
            [ 0.78750733,  0.84020696,  -0.71762832, 0.11890506,  0.72812158 ],
            [ 1.08479987,  0.63089915,  -0.69077067, 0.33239272,  0.70826765 ],
            [ 1.39410554,  0.42129131,  -0.66496362, 0.54842505,  0.68884329 ],
            [ 1.71973308,  0.20883136,  -0.64042442, 0.76960691,  0.66992350 ],
            [ 2.06620375,  -0.00897192, -0.61736509, 0.99848757,  0.65160153 ],
            [ 2.43943861,  -0.23528806, -0.59593337, 1.23832606,  0.63393898 ],
            [ 2.84615352,  -0.47354356, -0.57628572, 1.49263858,  0.61701935 ],
            [ 3.24944253,  -0.70260066, -0.56018013, 1.73846562,  0.60243633 ],
            [ 3.65834460,  -0.92883059, -0.54679867, 1.98212232,  0.58963673 ],
            [ 4.07887183,  -1.15634454, -0.53564449, 2.22763602,  0.57829695 ],
            [ 4.51564090,  -1.38819584, -0.52637664, 2.47795657,  0.56820499 ],
            [ 4.97278920,  -1.62700948, -0.51873516, 2.73559761,  0.55920584 ],
            [ 5.45410897,  -1.87512257, -0.51251140, 3.00276910,  0.55118292 ],
            [ 5.96348947,  -2.13485449, -0.50752527, 3.28165841,  0.54404066 ],
            [ 6.50597906,  -2.40906315, -0.50361040, 3.57501899,  0.53768871 ],
            [ 7.08777430,  -2.70113866, -0.50061755, 3.88614998,  0.53204794 ],
            [ 7.71561374,  -3.01470761, -0.49841575, 4.21857065,  0.52705696 ],
            [ 8.39481396,  -3.35267148, -0.49688948, 4.57500760,  0.52267864 ],
            [ 9.14481241,  -3.72493392, -0.49591334, 4.96550989,  0.51881231 ],
            [ 9.97967552,  -4.13870630, -0.49540040, 5.39717655,  0.51543051 ],
            [ 10.81993680, -4.55489548, -0.49526738, 5.82911089,  0.51277259 ],
            [ 11.71329827, -4.99738219, -0.49537781, 6.28617444,  0.51056287 ],
            [ 12.92039665, -5.59555924, -0.49575009, 6.90104102,  0.50830636 ],
            [ 14.35193348, -6.30564688, -0.49632430, 7.62723995,  0.50637483 ],
            [ 16.16094808, -7.20418000, -0.49706195, 8.54163858,  0.50466803 ],
            [ 17.97413328, -8.10605247, -0.49771495, 9.45554860,  0.50347012 ],
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
        #    Energy error to P=0.1 P0: 1.27e-08
        #    Est Max FD Error (based on N/2 run) to R=10.00: 5.4e-8
        #    Est Max MOC error for R>10.00: -1.62e-07
        #    Max interpolation error with cubic splines: 5.0857e-07
        #    Est Max overall error after cubic interpolation: 7.2e-7
        S32000_G1X4 => [
            [ -4.61789818, 12.00005700,  -2.99997850, -4.61896966, 0.99839192 ],
            [ -4.28161470, 10.99121774,  -2.99994103, -4.28338973, 0.99733513 ],
            [ -4.20868609, 10.77243643,  -2.99992661, -4.21066652, 0.99702648 ],
            [ -3.71151174, 9.28099075,   -2.99971068, -3.71569101, 0.99371858 ],
            [ -3.43113554, 8.43999369,   -2.99930987, -3.43750638, 0.99041533 ],
            [ -3.16713677, 7.64827666,   -2.99846957, -3.17661645, 0.98571980 ],
            [ -2.92945389, 6.93575830,   -2.99688234, -2.94301846, 0.97953476 ],
            [ -2.72385401, 6.31984376,   -2.99424058, -2.74235645, 0.97203917 ],
            [ -2.54540988, 5.78586779,   -2.99020621, -2.56964442, 0.96331858 ],
            [ -2.38871835, 5.31774642,   -2.98442571, -2.41944456, 0.95342726 ],
            [ -2.24488609, 4.88903814,   -2.97623741, -2.28310211, 0.94200956 ],
            [ -2.11252276, 4.49578575,   -2.96508677, -2.15924321, 0.92905860 ],
            [ -1.98841814, 4.12867789,   -2.95017507, -2.04482747, 0.91434989 ],
            [ -1.87005263, 3.78058313,   -2.93049897, -1.93756323, 0.89759655 ],
            [ -1.75496261, 3.44472517,   -2.90470932, -1.83533555, 0.87837254 ],
            [ -1.63970461, 3.11181032,   -2.87063845, -1.73536019, 0.85586838 ],
            [ -1.51803217, 2.76525434,   -2.82392072, -1.63286250, 0.82826655 ],
            [ -1.39023213, 2.40820272,   -2.76133470, -1.52909776, 0.79483249 ],
            [ -1.27461978, 2.09284157,   -2.69218514, -1.43914239, 0.76073294 ],
            [ -1.16721453, 1.80760040,   -2.61776782, -1.35927330, 0.72607170 ],
            [ -1.06528245, 1.54470494,   -2.53935786, -1.28703619, 0.69097415 ],
            [ -0.96785766, 1.30118607,   -2.45904221, -1.22141441, 0.65596327 ],
            [ -0.87146971, 1.06813169,   -2.37634493, -1.15989212, 0.62050795 ],
            [ -0.77566868, 0.84446950,   -2.29291336, -1.10214513, 0.58507002 ],
            [ -0.67951329, 0.62800295,   -2.20975492, -1.04758424, 0.54988932 ],
            [ -0.58208851, 0.41673960,   -2.12767789, -0.99571186, 0.51517794 ],
            [ -0.48255002, 0.20898476,   -2.04736872, -0.94613920, 0.48114621 ],
            [ -0.38003003, 0.00313137,   -1.96935328, -0.89852992, 0.44797661 ],
            [ -0.27375188, -0.20211745,  -1.89411810, -0.85264746, 0.41586468 ],
            [ -0.16284095, -0.40813560,  -1.82199439, -0.80826224, 0.38495971 ],
            [ -0.04646313, -0.61610893,  -1.75328058, -0.76520930, 0.35540826 ],
            [ 0.07635369,  -0.82736659,  -1.68814907, -0.72331636, 0.32730720 ],
            [ 0.20661735,  -1.04318922,  -1.62673247, -0.68244537, 0.30073577 ],
            [ 0.34542731,  -1.26490838,  -1.56910542, -0.64247224, 0.27574464 ],
            [ 0.49405610,  -1.49402798,  -1.51527037, -0.60326756, 0.25234815 ],
            [ 0.65210130,  -1.72949726,  -1.46572170, -0.56513174, 0.23076970 ],
            [ 0.81621598,  -1.96631546,  -1.42139639, -0.52888568, 0.21142329 ],
            [ 0.98777318,  -2.20666334,  -1.38157212, -0.49414658, 0.19399725 ],
            [ 1.16804146,  -2.45239933,  -1.34569060, -0.46063040, 0.17824815 ],
            [ 1.35828904,  -2.70525152,  -1.31329670, -0.42811140, 0.16397571 ],
            [ 1.55975529,  -2.96680888,  -1.28401967, -0.39641492, 0.15101564 ],
            [ 1.77368828,  -3.23859502,  -1.25754724, -0.36540113, 0.13922899 ],
            [ 2.00148149,  -3.52225353,  -1.23359894, -0.33494061, 0.12849081 ],
            [ 2.24458777,  -3.81944261,  -1.21193271, -0.30492576, 0.11869346 ],
            [ 2.50434694,  -4.13163839,  -1.19234897, -0.27528626, 0.10974886 ],
            [ 2.78277749,  -4.46108776,  -1.17462869, -0.24589877, 0.10156049 ],
            [ 3.08143108,  -4.80943502,  -1.15861595, -0.21671670, 0.09406124 ],
            [ 3.40257840,  -5.17912911,  -1.14413791, -0.18764379, 0.08717672 ],
            [ 3.74827989,  -5.57233201,  -1.13105748, -0.15862907, 0.08084988 ],
            [ 4.12098842,  -5.99161951,  -1.11924059, -0.12960955, 0.07502569 ],
            [ 4.52354898,  -6.43996798,  -1.10856114, -0.10051694, 0.06965327 ],
            [ 4.95879863,  -6.92030583,  -1.09891138, -0.07130857, 0.06469114 ],
            [ 5.43056776,  -7.43661636,  -1.09017888, -0.04190176, 0.06009531 ],
            [ 5.94248006,  -7.99260782,  -1.08227390, -0.01225765, 0.05583328 ],
            [ 6.49895373,  -8.59281068,  -1.07510858, 0.01768142,  0.05187307 ],
            [ 7.10500217,  -9.24234754,  -1.06860344, 0.04797265,  0.04818642 ],
            [ 7.76583484,  -9.94650281,  -1.06269108, 0.07865088,  0.04475073 ],
            [ 8.48765816,  -10.71157622, -1.05730779, 0.10976538,  0.04154427 ],
            [ 9.27687624,  -11.54402759, -1.05240064, 0.14134024,  0.03855013 ],
            [ 10.14069157, -12.45111444, -1.04792176, 0.17340035,  0.03575291 ],
            [ 11.08670558, -13.44046940, -1.04383039, 0.20595432,  0.03313991 ],
            [ 12.12311899, -14.52031283, -1.04009066, 0.23900281,  0.03069993 ],
            [ 13.25893219, -15.69966066, -1.03667047, 0.27254419,  0.02842260 ],
            [ 14.50354532, -16.98790900, -1.03354218, 0.30656219,  0.02629902 ],
            [ 15.86735844, -18.39545725, -1.03068038, 0.34104335,  0.02432039 ],
            [ 17.35403478, -19.92575035, -1.02807440, 0.37580018,  0.02248710 ],
            [ 18.96327961, -21.57821864, -1.02571426, 0.41059295,  0.02079943 ],
            [ 20.73034343, -23.38874714, -1.02354429, 0.44591679,  0.01922357 ],
            [ 22.66960032, -25.37167253, -1.02155126, 0.48173310,  0.01775473 ],
            [ 24.79622727, -27.54212884, -1.01972275, 0.51799776,  0.01638814 ],
            [ 27.12607839, -29.91592105, -1.01804721, 0.55466046,  0.01511903 ],
            [ 29.67549429, -32.50933384, -1.01651370, 0.59166410,  0.01394267 ],
            [ 32.46103042, -35.33885924, -1.01511201, 0.62894410,  0.01285438 ],
            [ 35.49908643, -38.42082525, -1.01383254, 0.66642771,  0.01184956 ],
            [ 38.80541696, -41.77090598, -1.01266634, 0.70403326,  0.01092370 ],
            [ 42.39450405, -45.40349327, -1.01160499, 0.74166950,  0.01007240 ],
            [ 46.37981353, -49.43302640, -1.01061770, 0.78017273,  0.00927267 ],
            [ 50.69697674, -53.79403804, -1.00972221, 0.81857972,  0.00854047 ],
            [ 55.35616761, -58.49658292, -1.00891148, 0.85677030,  0.00787168 ],
            [ 60.51080326, -63.69515772, -1.00815909, 0.89568743,  0.00724575 ],
            [ 66.21946063, -69.44834248, -1.00746167, 0.93533374,  0.00666081 ],
            [ 72.35529514, -75.62797854, -1.00683404, 0.97453012,  0.00613030 ],
            [ 79.13083143, -82.44779365, -1.00625338, 1.01434355,  0.00563587 ],
            [ 86.36744802, -89.72772703, -1.00573325, 1.05346946,  0.00518987 ],
            [ 94.32295787, -97.72688195, -1.00525305, 1.09306239,  0.00477537 ],
            [
                103.06762307, -106.51549052, -1.00481029, 1.13309201,
                0.00439074
            ],
            [
                112.67681522, -116.16888900, -1.00440260, 1.17352260,
                0.00403441
            ],
            [
                123.23073263, -126.76723432, -1.00402775, 1.21431243,
                0.00370485
            ],
            [
                134.81385274, -138.39495572, -1.00368359, 1.25541299,
                0.00340056
            ],
            [
                146.97310079, -150.59709748, -1.00338042, 1.29507748,
                0.00313110
            ],
            [
                160.19059853, -163.85738912, -1.00310285, 1.33477872,
                0.00288316
            ],
            [
                175.22189833, -178.93327819, -1.00283789, 1.37627726,
                0.00264534
            ],
            [
                191.61405531, -195.36991416, -1.00259615, 1.41780170,
                0.00242734
            ],
            [
                209.44629116, -213.24642424, -1.00237600, 1.45925730,
                0.00222792
            ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 4.91e-08
        #    Est Max FD Error (based on N/2 run) to R=50.01: 5.2e-8
        #    Est Max MOC error for R>50.01: 1.79e-07
        #    Max interpolation error with cubic splines: 5.0703e-07
        #    Est Max overall error after cubic interpolation: 7.3.e-7
        C16000_G1X4 => [
            [ -6.77632586, 12.00009717,  -1.99998567, -6.77766537, 0.99865960 ],
            [ -6.42946765, 11.30638774,  -1.99996678, -6.43136306, 0.99810282 ],
            [ -6.33204500, 11.11154539,  -1.99996249, -6.33413456, 0.99790829 ],
            [ -5.65663731, 9.76077908,   -1.99986548, -5.66074700, 0.99588213 ],
            [ -5.03895534, 8.52557650,   -1.99954647, -5.04659009, 0.99233780 ],
            [ -4.62390128, 7.69576248,   -1.99895764, -4.63548483, 0.98835523 ],
            [ -4.26892451, 6.98634739,   -1.99788318, -4.28548090, 0.98332363 ],
            [ -3.96394684, 6.37728016,   -1.99611554, -3.98646387, 0.97727223 ],
            [ -3.69704238, 5.84483835,   -1.99340439, -3.72652955, 0.97017361 ],
            [ -3.46115640, 5.37504830,   -1.98949375, -3.49859841, 0.96205131 ],
            [ -3.24609871, 4.94774365,   -1.98398969, -3.29266873, 0.95271629 ],
            [ -3.04783803, 4.55509430,   -1.97648592, -3.10480084, 0.94208621 ],
            [ -2.86228708, 4.18923205,   -1.96647364, -2.93108174, 0.93000675 ],
            [ -2.68535860, 3.84241432,   -1.95327261, -2.76772161, 0.91621769 ],
            [ -2.51352521, 3.50819033,   -1.93599796, -2.61160930, 0.90037188 ],
            [ -2.34180380, 3.17760557,   -1.91323728, -2.45854881, 0.88181069 ],
            [ -2.16132224, 2.83498466,   -1.88219355, -2.30139795, 0.85908138 ],
            [ -1.97039625, 2.47946973,   -1.84029926, -2.13997154, 0.83124667 ],
            [ -1.79783499, 2.16578623,   -1.79398692, -1.99894477, 0.80275121 ],
            [ -1.63759154, 1.88222550,   -1.74409097, -1.87260806, 0.77366447 ],
            [ -1.48564395, 1.62115384,   -1.69148099, -1.75727603, 0.74410012 ],
            [ -1.34032654, 1.37924013,   -1.63747710, -1.65128611, 0.71445289 ],
            [ -1.19680607, 1.14819882,   -1.58189007, -1.55090151, 0.68434076 ],
            [ -1.05431145, 0.92678136,   -1.52579703, -1.45553818, 0.65413316 ],
            [ -0.91141959, 0.71276364,   -1.46987067, -1.36422372, 0.62403001 ],
            [ -0.76668537, 0.50404371,   -1.41462408, -1.27607463, 0.59419857 ],
            [ -0.61890439, 0.29901880,   -1.36053997, -1.19044873, 0.56483824 ],
            [ -0.46681745, 0.09613807,   -1.30798203, -1.10674876, 0.53612503 ],
            [ -0.30919539, -0.10597893,  -1.25725211, -1.02446801, 0.50823331 ],
            [ -0.14481349, -0.30858976,  -1.20859961, -0.94316594, 0.48133161 ],
            [ 0.02763579,  -0.51294463,  -1.16221175, -0.86241820, 0.45556840 ],
            [ 0.20954738,  -0.72028937,  -1.11823077, -0.78181388, 0.43107381 ],
            [ 0.40245132,  -0.93191855,  -1.07675728, -0.70093346, 0.40795557 ],
            [ 0.60806247,  -1.14922387,  -1.03785510, -0.61933003, 0.38629750 ],
            [ 0.82841910,  -1.37382580,  -1.00154402, -0.53648109, 0.36615309 ],
            [ 1.06402318,  -1.60574956,  -0.96806773, -0.45244676, 0.34769194 ],
            [ 1.30806270,  -1.83827819,  -0.93837298, -0.36963258, 0.33145392 ],
            [ 1.56255098,  -2.07363264,  -0.91196695, -0.28715151, 0.31717073 ],
            [ 1.82921681,  -2.31360500,  -0.88848293, -0.20429434, 0.30463894 ],
            [ 2.10978815,  -2.55987736,  -0.86762121, -0.12040619, 0.29368897 ],
            [ 2.40599868,  -2.81405657,  -0.84913439, -0.03486807, 0.28417790 ],
            [ 2.71941619,  -3.07755519,  -0.83282344, 0.05286939,  0.27598719 ],
            [ 3.05183429,  -3.35194431,  -0.81850505, 0.14340849,  0.26900565 ],
            [ 3.40483675,  -3.63860194,  -0.80602674, 0.23729046,  0.26313672 ],
            [ 3.78047383,  -3.93927677,  -0.79523521, 0.33518282,  0.25828274 ],
            [ 4.18063073,  -4.25557799,  -0.78599692, 0.43771274,  0.25435450 ],
            [ 4.60743149,  -4.58930774,  -0.77818072, 0.54557635,  0.25126280 ],
            [ 5.06330757,  -4.94251154,  -0.77165705, 0.65955401,  0.24891822 ],
            [ 5.55096048,  -5.31744619,  -0.76629854, 0.78049887,  0.24723159 ],
            [ 6.07427958,  -5.71727707,  -0.76197252, 0.90956068,  0.24611227 ],
            [ 6.63733948,  -6.14529929,  -0.75855435, 1.04793393,  0.24547297 ],
            [ 7.24600336,  -6.60615406,  -0.75591699, 1.19725095,  0.24522768 ],
            [ 7.90912290,  -7.10671761,  -0.75393557, 1.35987377,  0.24529516 ],
            [ 8.63753858,  -7.65533379,  -0.75249654, 1.53865140,  0.24560056 ],
            [ 9.44728231,  -8.26422315,  -0.75149270, 1.73771162,  0.24607575 ],
            [ 10.36197849, -8.95127721,  -0.75082697, 1.96306316,  0.24666034 ],
            [ 11.41824558, -9.74410823,  -0.75041355, 2.22394716,  0.24730191 ],
            [ 12.67609758, -10.68784982, -0.75017889, 2.53544152,  0.24795543 ],
            [ 14.24414585, -11.86405919, -0.75006200, 2.92476532,  0.24858236 ],
            [ 16.35000171, -13.44352035, -0.75001468, 3.44888797,  0.24914869 ],
            [ 19.61228403, -15.89025115, -0.75000155, 4.26255051,  0.24961999 ],
            [ 27.12887201, -21.52769436, -0.75000001, 6.14040822,  0.24994149 ],
        ],

        #    Overall Error Estimate:
        #    Energy error to P=0.1 P0: 2.24e-07
        #    Est Max FD Error (based on N/2 run) to R=20.70: 2.8e-7
        #    Est Max MOC error for R>20.70: 3.86e-07
        #    Max interpolation error with cubic splines: 5.0464e-07
        #    Est Max overall error after cubic interpolation: 1.1e-6
        'P8000_G1X4' => [
            [
                -12.37479245, 12.00006633, -0.99999283, -12.37657889,
                0.99910599
            ],
            [
                -11.01496637, 10.64026364, -0.99997073, -11.01849531,
                0.99823248
            ],
            [ -9.99397706, 9.61933095,   -0.99991177, -9.99986341, 0.99704844 ],
            [ -8.95507173, 8.58058739,   -0.99974986, -8.96498634, 0.99501944 ],
            [ -8.11861745, 7.74446034,   -0.99942292, -8.13371661, 0.99239807 ],
            [ -7.42101233, 7.04743694,   -0.99884277, -7.44247438, 0.98916712 ],
            [ -6.82463275, 6.45199975,   -0.99790509, -6.85364460, 0.98531635 ],
            [ -6.30159688, 5.93040018,   -0.99648120, -6.33941511, 0.98080528 ],
            [ -5.83454858, 5.46543984,   -0.99442241, -5.88249895, 0.97559551 ],
            [ -5.41005083, 5.04387839,   -0.99154922, -5.46958155, 0.96962442 ],
            [ -5.01802587, 4.65588466,   -0.98764612, -5.09075518, 0.96280994 ],
            [ -4.65041149, 4.29371356,   -0.98245056, -4.73819538, 0.95504302 ],
            [ -4.29959338, 3.95018904,   -0.97562005, -4.40466217, 0.94615722 ],
            [ -3.95784526, 3.61822783,   -0.96668259, -4.08301976, 0.93589722 ],
            [ -3.61486541, 3.28860940,   -0.95488413, -3.76404600, 0.92380084 ],
            [ -3.25027436, 2.94330625,   -0.93862434, -3.42990599, 0.90877072 ],
            [ -2.87459777, 2.59453509,   -0.91734383, -3.09179569, 0.89081506 ],
            [ -2.53363347, 2.28564074,   -0.89389606, -2.79114857, 0.87236320 ],
            [ -2.21602835, 2.00565096,   -0.86873942, -2.51704066, 0.85347415 ],
            [ -1.92231756, 1.75420915,   -0.84310431, -2.26908735, 0.83476807 ],
            [ -1.64148392, 1.52106603,   -0.81706720, -2.03727097, 0.81603370 ],
            [ -1.36333950, 1.29749045,   -0.79048002, -1.81294076, 0.79695283 ],
            [ -1.07990439, 1.07730147,   -0.76327608, -1.58984180, 0.77728727 ],
            [ -0.78830319, 0.85874702,   -0.73587217, -1.36612576, 0.75715353 ],
            [ -0.48999829, 0.64326718,   -0.70907651, -1.14328603, 0.73698315 ],
            [ -0.18137686, 0.42848055,   -0.68317460, -0.91896484, 0.71686302 ],
            [ 0.14142456,  0.21201453,   -0.65841299, -0.69081738, 0.69688197 ],
            [ 0.48284327,  -0.00870037,  -0.63499466, -0.45630388, 0.67712432 ],
            [ 0.84800876,  -0.23648290,  -0.61309556, -0.21264576, 0.65767957 ],
            [ 1.24322470,  -0.47467581,  -0.59286105, 0.04344843,  0.63863616 ],
            [ 1.65601953,  -0.71564675,  -0.57520427, 0.30334102,  0.62090796 ],
            [ 2.06930671,  -0.95026125,  -0.56063300, 0.55664734,  0.60523748 ],
            [ 2.49092724,  -1.18398788,  -0.54849317, 0.80880510,  0.59121228 ],
            [ 2.92640168,  -1.42055701,  -0.53836677, 1.06345040,  0.57859315 ],
            [ 3.38023141,  -1.66290040,  -0.52995326, 1.32338941,  0.56722625 ],
            [ 3.86026431,  -1.91554885,  -0.52296922, 1.59313937,  0.55692990 ],
            [ 4.35785127,  -2.17431832,  -0.51737510, 1.86794350,  0.54787012 ],
            [ 4.89307828,  -2.44995138,  -0.51281143, 2.15892133,  0.53968338 ],
            [ 5.46449932,  -2.74189835,  -0.50920459, 2.46517109,  0.53243990 ],
            [ 6.07976569,  -3.05428465,  -0.50640583, 2.79073529,  0.52607084 ],
            [ 6.74671997,  -3.39128393,  -0.50428950, 3.13968238,  0.52053134 ],
            [ 7.47740041,  -3.75915092,  -0.50273647, 3.51820958,  0.51576536 ],
            [ 8.28618788,  -4.16528028,  -0.50164349, 3.93364564,  0.51173103 ],
            [ 9.19343961,  -4.62003626,  -0.50091619, 4.39631489,  0.50838374 ],
            [ 10.22567553, -5.13684384,  -0.50047039, 4.91960790,  0.50568419 ],
            [ 11.42780311, -5.73830649,  -0.50022761, 5.52614728,  0.50357831 ],
            [ 12.86088433, -6.45507998,  -0.50011727, 6.24660401,  0.50202377 ],
            [ 14.78381232, -7.41671609,  -0.50006989, 7.21075190,  0.50090447 ],
            [ 19.28989861, -9.66994172,  -0.50001699, 9.46554066,  0.50011249 ],
            [ 23.80817302, -11.92911182, -0.50000208, 11.72488212, 0.50001208 ],
        ],

        'P8000_G1X2' => [
            [
                -13.01924624, 12.00006662, -0.99999330, -13.02097366,
                0.99913555
            ],
            [
                -11.74508114, 10.72592253, -0.99996892, -11.74835009,
                0.99836291
            ],
            [ -10.29668273, 9.27762108, -0.99987332, -10.30343822, 0.99661128 ],
            [ -9.32619421,  8.30734036, -0.99966554, -9.33719102,  0.99447325 ],
            [ -8.53822516,  7.51977205, -0.99926547, -8.55457212,  0.99176591 ],
            [ -7.87689999,  6.85913461, -0.99858016, -7.89971815,  0.98847752 ],
            [ -7.30638666,  6.28971146, -0.99749691, -7.33683466,  0.98458384 ],
            [ -6.80339434,  5.78835102, -0.99588243, -6.84268657,  0.98005247 ],
            [ -6.35112676,  5.33842929, -0.99357633, -6.40057493,  0.97483176 ],
            [ -5.93726665,  4.92784496, -0.99038480, -5.99832438,  0.96885054 ],
            [ -5.55296797,  4.54801969, -0.98607953, -5.62726458,  0.96202576 ],
            [ -5.19005655,  4.19114039, -0.98036889, -5.27950625,  0.95422694 ],
            [ -4.84103005,  3.85020918, -0.97287071, -4.94797488,  0.94526466 ],
            [ -4.49715440,  3.51727787, -0.96302317, -4.62466833,  0.93482003 ],
            [ -4.14521654,  3.18057287, -0.94984088, -4.29782089,  0.92226554 ],
            [ -3.75050417,  2.80927985, -0.93067491, -3.93696574,  0.90573567 ],
            [ -3.38804065,  2.47581068, -0.90861694, -3.61178237,  0.88818029 ],
            [ -3.05502659,  2.17713001, -0.88461111, -3.31896766,  0.87009022 ],
            [ -2.74351176,  1.90544249, -0.85926738, -3.05076222,  0.85162210 ],
            [ -2.45327220,  1.65972368, -0.83368230, -2.80622530,  0.83328727 ],
            [ -2.17090235,  1.42797989, -0.80759774, -2.57354487,  0.81466042 ],
            [ -1.88749328,  1.20287667, -0.78091465, -2.34537410,  0.79546825 ],
            [ -1.59532906,  0.97872609, -0.75358747, -2.11588562,  0.77547947 ],
            [ -1.29839390,  0.75898523, -0.72666438, -1.88862442,  0.75527858 ],
            [ -0.99334759,  0.54135865, -0.70046354, -1.66134323,  0.73496132 ],
            [ -0.67646290,  0.32344662, -0.67524167, -1.43169323,  0.71461112 ],
            [ -0.34368585,  0.10281075, -0.65122076, -1.19729717,  0.69431694 ],
            [ 0.00962984,  -0.12319069,  -0.62858881, -0.95558754, 0.67417146 ],
            [ 0.38900813,  -0.35756176,  -0.60750478, -0.70365393, 0.65427261 ],
            [ 0.80131514,  -0.60392059,  -0.58810048, -0.43799701, 0.63472242 ],
            [ 1.21820172,  -0.84558479,  -0.57179349, -0.17713607, 0.61710057 ],
            [ 1.63899933,  -1.08325003,  -0.55825044, 0.07916011,  0.60137663 ],
            [ 2.07020253,  -1.32144560,  -0.54693971, 0.33535961,  0.58724189 ],
            [ 2.51724141,  -1.56375825,  -0.53748790, 0.59496215,  0.57449542 ],
            [ 2.98429285,  -1.81288532,  -0.52962751, 0.86053249,  0.56301962 ],
            [ 3.46442832,  -2.06558940,  -0.52326700, 1.12837348,  0.55294000 ],
            [ 3.99507719,  -2.34175390,  -0.51783417, 1.41922338,  0.54354768 ],
            [ 4.54953187,  -2.62762149,  -0.51353667, 1.71826540,  0.53540149 ],
            [ 5.14440164,  -2.93203645,  -0.51010602, 2.03455751,  0.52824713 ],
            [ 5.78698754,  -3.25890949,  -0.50740980, 2.37193053,  0.52203816 ],
            [ 6.48716390,  -3.61341271,  -0.50532516, 2.73551016,  0.51672440 ],
            [ 7.25735988,  -4.00196343,  -0.50374212, 3.13168829,  0.51226039 ],
            [ 8.11497494,  -4.43343846,  -0.50256016, 3.56935432,  0.50859714 ],
            [ 9.08677925,  -4.92137379,  -0.50168868, 4.06210340,  0.50567936 ],
            [ 10.21650530, -5.48775538,  -0.50105075, 4.63202571,  0.50344858 ],
            [ 11.56222480, -6.16169386,  -0.50059354, 5.30835570,  0.50186117 ],
            [ 13.17437370, -6.96845425,  -0.50029045, 6.11653122,  0.50086580 ],
            [ 15.04799172, -7.90562699,  -0.50012090, 7.05440548,  0.50034735 ],
            [ 17.30573323, -9.03466385,  -0.50004032, 8.18374888,  0.50011360 ],
            [ 20.36617188, -10.56494670, -0.50000885, 9.71414663,  0.50002472 ],
        ],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 2.01e-07
#    Est Max FD Error (based on N/2 run) to R=10.00: 4.7e-7
#    Est Max MOC error for R>10.00: -6.23e-08 
#    Max interpolation error with cubic splines: 5.0436e-07
#    Est Max overall error after cubic interpolation: 1.e-6

        'S8000_G1X15' => [
[-4.91112446, 12.00017628, -2.99998029, -4.91215040, 0.99846031],
[-4.41115547, 10.50029178, -2.99991166, -4.41332851, 0.99673698],
[-3.90334821, 8.97697380, -2.99960949, -3.90800830, 0.99299442],
[-3.64524552, 8.20282083, -2.99913850, -3.65211590, 0.98966167],
[-3.37986444, 7.40702752, -2.99809014, -3.39010971, 0.98456213],
[-3.15882552, 6.74450949, -2.99629856, -3.17312355, 0.97842338],
[-2.96387680, 6.16064224, -2.99337924, -2.98306961, 0.97099172],
[-2.79087522, 5.64313503, -2.98893349, -2.81580862, 0.96225842],
[-2.63601360, 5.18072350, -2.98252231, -2.66753854, 0.95221782],
[-2.49828418, 4.77050102, -2.97384641, -2.53713158, 0.94106227],
[-2.36704932, 4.38097231, -2.96178107, -2.41445792, 0.92803832],
[-2.24413407, 4.01785308, -2.94578838, -2.30126630, 0.91329964],
[-2.12611228, 3.67136692, -2.92470868, -2.19444416, 0.89643456],
[-2.01077788, 3.33556315, -2.89714171, -2.09214487, 0.87701050],
[-1.89219621, 2.99412429, -2.85993107, -1.98949710, 0.85364785],
[-1.76331283, 2.62876577, -2.80743698, -1.88133646, 0.82402295],
[-1.63857252, 2.28242129, -2.74338098, -1.78056300, 0.79099548],
[-1.52467273, 1.97383904, -2.67330370, -1.69236193, 0.75719566],
[-1.41796043, 1.69248831, -2.59838014, -1.61337866, 0.72268508],
[-1.31608595, 1.43172709, -2.51989768, -1.54152599, 0.68762831],
[-1.21732388, 1.18681690, -2.43908653, -1.47535377, 0.65222129],
[-1.12010343, 0.95367021, -2.35688089, -1.41367256, 0.61659633],
[-1.02335691, 0.72964539, -2.27430894, -1.35574242, 0.58099385],
[-0.92605051, 0.51235016, -2.19217883, -1.30093481, 0.54562045],
[-0.82738247, 0.30006753, -2.11129331, -1.24883069, 0.51073439],
[-0.72638090, 0.09085511, -2.03215869, -1.19898839, 0.47650992],
[-0.62237889, -0.11645777, -1.95540265, -1.15118050, 0.44320645],
[-0.51443486, -0.32348219, -1.88135030, -1.10510078, 0.41097436],
[-0.40176184, -0.53140043, -1.81037657, -1.06056597, 0.37999808],
[-0.28355518, -0.74133590, -1.74278182, -1.01742420, 0.35043463],
[-0.15881488, -0.95466021, -1.67871458, -0.97549354, 0.32237379],
[-0.02646026, -1.17276439, -1.61827294, -0.93461482, 0.29587924],
[0.11455056, -1.39687231, -1.56156326, -0.89468435, 0.27101210],
[0.26549957, -1.62849655, -1.50858512, -0.85556984, 0.24778035],
[0.42573599, -1.86622951, -1.45988719, -0.81761921, 0.22642649],
[0.59255050, -2.10602439, -1.41619270, -0.78148611, 0.20726428],
[0.76716275, -2.34978907, -1.37687503, -0.74683946, 0.19001276],
[0.94923432, -2.59719424, -1.34169521, -0.71368592, 0.17456032],
[1.14239457, -2.85318933, -1.30973755, -0.68136101, 0.16049680],
[1.34680174, -3.11787923, -1.28085719, -0.64989072, 0.14775115],
[1.56366298, -3.39274040, -1.25474690, -0.61913638, 0.13618177],
[1.79590255, -3.68130722, -1.23098217, -0.58877132, 0.12559480],
[2.04220877, -3.98180088, -1.20960761, -0.55904835, 0.11600691],
[2.30620812, -4.29850517, -1.19021586, -0.52961128, 0.10723354],
[2.58893609, -4.63246308, -1.17267864, -0.50045668, 0.09921634],
[2.89230559, -4.98574365, -1.15682090, -0.47150019, 0.09187704],
[3.21818422, -5.36032475, -1.14249402, -0.44268439, 0.08515070],
[3.56883154, -5.75860170, -1.12955178, -0.41393794, 0.07897452],
[3.94689869, -6.18337196, -1.11785643, -0.38518306, 0.07329052],
[4.35502826, -6.63738196, -1.10729098, -0.35636827, 0.06805139],
[4.79625551, -7.12377965, -1.09774424, -0.32743747, 0.06321335],
[5.27420893, -7.64632465, -1.08910975, -0.29832262, 0.05873518],
[5.79271064, -8.20894264, -1.08129484, -0.26897299, 0.05458270],
[6.35637747, -8.81637563, -1.07421072, -0.23932261, 0.05072340],
[6.97002167, -9.47352645, -1.06778185, -0.20932742, 0.04713122],
[7.63885184, -10.18567870, -1.06194110, -0.17895327, 0.04378393],
[8.36947377, -10.95955427, -1.05662236, -0.14813540, 0.04065857],
[9.16849112, -11.80181782, -1.05177297, -0.11684575, 0.03773849],
[10.04290611, -12.71951050, -1.04734752, -0.08507107, 0.03500992],
[11.00031995, -13.72026228, -1.04330576, -0.05280477, 0.03246065],
[12.04933330, -14.81270753, -1.03961090, -0.02003556, 0.03007902],
[13.19894647, -16.00585703, -1.03623172, 0.01323128, 0.02785535],
[14.45875959, -17.30930973, -1.03314074, 0.04698198, 0.02578092],
[15.83917271, -18.73346049, -1.03031316, 0.08119977, 0.02384750],
[17.24307967, -20.17818086, -1.02789702, 0.11346982, 0.02216680],
[18.88624914, -21.86517869, -1.02552029, 0.14847547, 0.02048641],
[20.60571292, -23.62667892, -1.02343433, 0.18237941, 0.01898827],
[22.55449584, -25.61913934, -1.02145004, 0.21793561, 0.01754182],
[24.57765276, -27.68389236, -1.01971884, 0.25209712, 0.01626198],
[26.86743409, -30.01688264, -1.01807011, 0.28788109, 0.01502683],
[29.22047725, -32.41071861, -1.01664169, 0.32193231, 0.01394322],
[31.87710798, -35.10969870, -1.01527930, 0.35754812, 0.01289743],
[34.89073161, -38.16735547, -1.01398196, 0.39485576, 0.01188988],
[37.96212927, -41.27994091, -1.01286915, 0.42999750, 0.01101610],
[41.43413882, -44.79473628, -1.01180756, 0.46674383, 0.01017390],
[45.37825464, -48.78336889, -1.01079648, 0.50522768, 0.00936355],
[49.35121853, -52.79747762, -1.00993962, 0.54101540, 0.00867020],
[53.84232368, -57.33133318, -1.00912167, 0.57841382, 0.00800242],
[58.94426259, -62.47776244, -1.00834216, 0.61755734, 0.00736037],
[63.99871275, -67.57268809, -1.00769128, 0.65335853, 0.00681982],
[69.70162667, -73.31764086, -1.00706918, 0.69072697, 0.00629922],
[76.16728688, -79.82703268, -1.00647555, 0.72979205, 0.00579867],
[82.42199174, -86.12067948, -1.00598917, 0.76473715, 0.00538567],
[89.44753400, -93.18661818, -1.00552330, 0.80113957, 0.00498754],
[97.37442141, -101.15547270, -1.00507776, 0.83911545, 0.00460434],
[106.36186551, -110.18658539, -1.00465235, 0.87879554, 0.00423614],
[114.80149468, -118.66401092, -1.00431310, 0.91327024, 0.00394079],
[124.25053910, -128.15223185, -1.00398759, 0.94912698, 0.00365592],
[134.87497442, -138.81733136, -1.00367572, 0.98647326, 0.00338156],
[146.87620919, -150.86083866, -1.00337739, 1.02542973, 0.00311774],
[160.50054627, -164.52919267, -1.00309249, 1.06613241, 0.00286448],
[172.77043062, -176.83564942, -1.00287417, 1.10005472, 0.00266951],
[186.46084482, -190.56393935, -1.00266431, 1.13528151, 0.00248132],
[201.79802500, -205.94039736, -1.00246288, 1.17191152, 0.00229994],
],
    };
}

1;
