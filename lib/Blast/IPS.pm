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
my $ralpha_table;

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
    my ( $self, @args ) = @_;

    my @input_keys = qw(
      rtable
      table_name
      gamma
      symmetry
    );
    my %valid_input_keys;
    @valid_input_keys{@input_keys} = (1) x scalar(@input_keys);

    # Convert the various input methods to a hash ref
    my $rinput_hash;
    my $reftype;
    if (@args) {
        my $arg0    = $args[0];
        my $reftype = ref($arg0);
        if ( !$reftype ) {

            if ( defined($arg0) ) {

                # simple hash of named values
                my %input_hash = @args;
                $rinput_hash = \%input_hash;
            }

        }
        elsif ( $reftype eq 'HASH' ) {
            $rinput_hash = $arg0;
        }
        elsif ( $reftype eq "ARRAY" ) {

            # useful? maybe delete this option?
            $rinput_hash = { rtable => $arg0 };
        }
        else {
            carp "Unexpected ref type: $reftype\n";
            return;
        }
    }
    else {

        # default to spherical table with gamma 1.4 if no arg given
        $rinput_hash = { symmetry => 2, gamma => 1.4 };
    }

    # Validate input keys
    _check_keys( $rinput_hash, \%valid_input_keys,
        "Checking for valid input keys" );

    # The following four quantities are needed:
    my $rtable     = $rinput_hash->{rtable};
    my $table_name = $rinput_hash->{table_name};
    my $gamma      = $rinput_hash->{gamma};
    my $symmetry   = $rinput_hash->{symmetry};

    # Allow input of the older 'ASYM' keyword for the symmetry
    if ( !defined($symmetry) ) { $symmetry = $rinput_hash->{ASYM}; }

    # There are three input options:

    # Option 1: a user-supplied table is given
    if ( defined($rtable) ) {

        if ( !defined($gamma) )      { $gamma      = 1.4 }
        if ( !defined($symmetry) )   { $symmetry   = 2 }
        if ( !defined($table_name) ) { $table_name = 'NONAME' }
    }

    # Option 2: a builtin table name is given
    elsif ( defined($table_name) ) {
        ( $rtable, $table_name ) = get_builtin_table($table_name);
        my $item = $rtables_info->{$table_name};
        if ( defined($item) ) {
            ( $symmetry, $gamma, my $err_est ) = @{$item};
        }
        else {
            carp("Unknown table name: '$table_name'");
            return;
        }
    }

    # Option 3: lookup a table given symmetry and gamma
    else {
        if ( !defined($gamma) )    { $gamma    = 1.4 }
        if ( !defined($symmetry) ) { $symmetry = 2 }

	# allow simple abbreviations for symmetry
        $symmetry = uc $symmetry;
        if    ( $symmetry eq 'S' ) { $symmetry = 2 }
        elsif ( $symmetry eq 'C' ) { $symmetry = 1 }
        elsif ( $symmetry eq 'P' ) { $symmetry = 0 }
        my @keys = reverse sort keys %{$rtables_info};

        # Look through our list of tables for the best match
        my $err_est;

        # FIXME: interpolate if necessary

        # allow a small roundoff error when comparing gamma
        my $eps = 1.e-6;
        foreach my $key (@keys) {
            my $item = $rtables_info->{$key};
            my ( $symmetry_t, $gamma_t, $err_est_t ) = @{$item};
            next if $symmetry_t != $symmetry;
            next if abs( $gamma_t - $gamma ) > $eps;

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
        else {

            # Make a table by interpolation

        }
    }

##    # A scalar value is the name of a pre-defined table
##    elsif ( !$reftype ) {
##        ( $rtable, $table_name ) = get_builtin_table($rinput_hash);
##    }
##    else {
##        croak "Unexpected argument type $reftype in Blast::Table";
##    }
##
##?    if ( defined($table_name) ) {
##?        my $item = $rtables_info->{$table_name};
##?        if ( defined($item) ) {
##?            ( $symmetry, $gamma, my $err_est ) = @{$item};
##?        }
##?    }
##?    else {
##?        $table_name = "NONAME";
##?    }

    # Now check the results

    my $error = "";
##    if ( !defined($rtable) ) {
##        $error .=
##"Undefined Table in Blast::Table\nfor symmetry=$symmetry, gamma=$gamma, name=$table_name\n";
##    }
##
##    $gamma    = 1.4 unless defined($gamma);
##    $symmetry = 2   unless defined($symmetry);

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

sub _check_keys {
    my ( $rtest, $rvalid, $msg, $exact_match ) = @_;

    # Check the keys of a hash:
    # $rtest  = ref to hash to test
    # $rvalid = ref to has with valid keys

    # $msg = a message to write in case of error
    # $exact_match defines the type of check:
    #     = false: test hash must not have unknown key
    #     = true:  test hash must have exactly same keys as known hash
    my @unknown_keys =
      grep { $_ && !exists $rvalid->{$_} } keys %{$rtest};
    my @missing_keys =
      grep { !exists $rtest->{$_} } keys %{$rvalid};
    my $error = @unknown_keys;
    if ($exact_match) { $error ||= @missing_keys }
    if ($error) {
        local $" = ')(';
        my @expected_keys = sort keys %{$rvalid};
        @unknown_keys = sort @unknown_keys;
        croak(<<EOM);
------------------------------------------------------------------------
Program error detected checking hash keys
Message is: '$msg'
Expected keys: (@expected_keys)
Unknown key(s): (@unknown_keys)
Missing key(s): (@missing_keys)
------------------------------------------------------------------------
EOM
    }
    return;
}

sub get_builtin_table {
    my ($table_name) = @_;
    return ( $rtables->{$table_name}, $table_name );
}

sub get_info {
    my ($self) = @_;
    if (!ref($self)) {
	my $what=ref $self;
	print STDERR "$what\n";
	croak("get_info not called by a Blast::IPS object");
    }

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

    # returns hash ref with index to all built-in tables
    # NOTE:
    # -  The hash key is the name of the builtin table ('S1.4', 'C1.667', etc);
    # -  symmetry = 0, 1, or 2 for plane, cylindrical or spherical
    # -  gamma is the ideal gas gamma
    # -  error_estimate = the estimated maximum relative error in shock
    #      overpressure after interpolation using the standard cubic
    #      interpolation
    # -  N = the number of points across the blast wave in the finite difference
    #      calculation used to make the table
    # Additional specific parameters can be obtained by loading a table and
    # using the 'get_info()' method.
    my $rtable_index = {};
    foreach my $key ( sort keys %{$rtables_info} ) {
        my $item = $rtables_info->{$key};
        my ( $symmetry, $gamma, $err_est, $N ) = @{$item};
        $rtable_index->{$key} = {
            symmetry       => $symmetry,
            gamma          => $gamma,
            error_estimate => $err_est,
            N              => $N,
        };
    }
    return ($rtable_index);
}

sub get_table_index_as_array_ref {

    # OLD: to be deleted
    # returns a list of references of the form
    #   [NAME, symmetry, gamma, error]
    # one per built-in table

    my @table_list;
    foreach my $key ( sort keys %{$rtables_info} ) {
        my $item = $rtables_info->{$key};
        my ( $symmetry, $gamma, $err_est ) = @{$item};
        push @table_list, [ $key, $symmetry, $gamma, $err_est ];
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
        $error = <<EOM;
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
        my $rtable = $rtables->{$key};

        #print STDERR "Checking table $key\n";

 	# Check for a naming error
        my $item = $rtables_info->{$key};
        my ( $symmetry, $gamma, $err_est ) = @{$item};
	my $letter=$symmetry eq 0 ? 'P' : $symmetry eq 1 ? 'C' : 'S';
	my $str=$letter.$gamma;
        #print STDERR "Checking if table $key similar to $str\n";
	if ($key !~ /^$str/i) {
            return "Table $key: expected to be like $str\n";
	}

 	# Check for numerical problems
        my $error = _check_table($rtable);
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
    my ( $X_e, $Y_e, $dYdX_e, $Z_e, $dZdX_e ) = @{ $rtab->[-1] };

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

    ## PATCH: Return ending dZdX
    $dZ_dX_i = $dZdX_e;

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

sub alpha_interpolate {
    my ( $sym, $gamma ) = @_;

    # Given: a 1d symmetry (0, 1, or 2) and an ideal gas gamma
    # return: parameter alpha for the similarity solution
    # returns undef if out of bounds of table 
    # (currently gamma<1.1 || $gamma>7 )

    # alpha is obtained by interpolating a table of pre-computed values.
    # The estimated accuracy over most of the range of the table is about 1.e-6
    # The table is spaced closely enough that cubic interpolation of the 
    # interpolated values have comparable accuracy.

    return if ( $sym != 0 && $sym != 1 && $sym != 2 );

    my $rtab = $ralpha_table->[$sym];
    my $ntab = @{$rtab};
    my ( $jl, $ju );
    my $icol = 0;
    my $eps  = 1.e-5;

    my $rhash = {
        _jl     => $jl,
        _ju     => $ju,
        _rtable => $rtab,
    };
    ( $jl, $ju ) = _locate_2d( $rhash, $gamma, $icol );
    my ( $gamma_min, $alpha_min ) = @{ $rtab->[0] };
    my ( $gamma_max, $alpha_max ) = @{ $rtab->[1] };
    if ( $jl < 0 ) {
        return if ( $gamma + $eps < $gamma_min );
        return $alpha_min;
    }
    if ( $ju >= $ntab ) {
        return if ( $gamma - $eps > $gamma_max );
        return $alpha_max;
    }

    # (gamma, alpha, dalpha/dgamma) =
    my ( $x1, $y1, $dydx1 ) = @{ $rtab->[$jl] };
    my ( $x2, $y2, $dydx2 ) = @{ $rtab->[$ju] };
    my ( $alpha, $dydx, $d2ydx2, $d3ydx3 ) =
      _cubic_interpolation( $gamma, $x1, $y1, $dydx1, $x2, $y2, $dydx2 );
    return ($alpha);
}

sub _make_intermediate_table {

    # Given: a symmetry, gamma, and two tables with same symmetry and
    # bounding gamma values
    # Return: a table of X,Y,dYdX,Z,dZdX interpolated from the two 
    # bounding tables.

    # The error depends on the difference between the gamma value and the gamma
    # values of the bounding tables. For the current set of tables, the maximum
    # relative error in overpressure ratio for the generated table is about
    # 5.e-4.  I would like to reduce that.

    my ( $symmetry, $gamma_mid, $blast_table_lo, $blast_table_hi) = @_;
    my $rtable_new;
    my $alpha_mid = Blast::IPS::alpha_interpolate( $symmetry, $gamma_mid );

    # Evaluate the interpolation fraction
    my $gamma_lo = $blast_table_lo->get_gamma();
    my $gamma_hi = $blast_table_hi->get_gamma();
    my $alpha_lo = $blast_table_lo->get_alpha();
    my $alpha_hi = $blast_table_hi->get_alpha();

    my $A_lo  = 1 / ( ( $gamma_lo + 1 ) * $alpha_lo );
    my $A_mid = 1 / ( ( $gamma_mid + 1 ) * $alpha_mid );
    my $A_hi  = 1 / ( ( $gamma_hi + 1 ) * $alpha_hi );
    my $ff = log( $A_mid / $A_lo ) / log( $A_hi / $A_lo );

    # Small errors in alpha can cause ff to go out of bounds
    if ( $ff < 0 ) { $ff = 0 }
    if ( $ff > 1 ) { $ff = 1 }

    # Use the Y values of the closest table so that we converge to it. But clip
    # the Y to the lower bounds of the other table to avoid extrapolation at
    # long range where the theory is not great. 
    my $rtable_closest;
    my $rbounds;
    if ( $ff <= 0.5 ) {
        $rtable_closest = $blast_table_lo->get_table();
        $rbounds        = $blast_table_hi->get_table_bounds();
    }
    else {
        $rtable_closest = $blast_table_hi->get_table();
        $rbounds        = $blast_table_lo->get_table_bounds();
    }
    my $Ymax = $rbounds->[0]->[1];
    my $Ymin = $rbounds->[1]->[1];

    my $iQ = 'Y';

    my $rtab_int = [];
    my $stop;
    foreach my $item_ref(@{$rtable_closest}){
	last if ($stop);
	my $Y=$item_ref->[1];
        if ( $Y < $Ymin ) {
            $Y    = $Ymin;
            $stop = 1;
        }
        my $item_lo = $blast_table_lo->lookup( $Y, $iQ );
        my $item_hi = $blast_table_hi->lookup( $Y, $iQ );
        my ( $X_lo, $YY_lo, $dYdX_lo, $Z_lo, $dZdX_lo ) = @{$item_lo};
        my ( $X_hi, $YY_hi, $dYdX_hi, $Z_hi, $dZdX_hi ) = @{$item_hi};

        my $X_int = $X_lo + $ff * ( $X_hi - $X_lo );
        my $Z_int = $Z_lo + $ff * ( $Z_hi - $Z_lo );

        my $dXdY_lo  = 1 / $dYdX_lo;
        my $dXdY_hi  = 1 / $dYdX_hi;
        my $dXdY_int = $dXdY_lo + $ff * ( $dXdY_hi - $dXdY_lo );
        my $dYdX_int = 1 / $dXdY_int;

        my $dXdZ_lo  = 1 / $dZdX_lo;
        my $dXdZ_hi  = 1 / $dZdX_hi;
        my $dXdZ_int = $dXdZ_lo + $ff * ( $dXdZ_hi - $dXdZ_lo );
        my $dZdX_int = 1 / $dXdZ_int;

        #print STDERR "$Y\t$Z_lo\t$Z_int\t$Z_hi\t$dZdX_lo\t$dZdX_int\t$dZdX_hi\n";

        push @{$rtab_int}, [ $X_int, $Y, $dYdX_int, $Z_int, $dZdX_int ];
    }
    return $rtab_int;
}

##################################################################
# ALPHA TABLE
##################################################################
BEGIN {
    $ralpha_table = [
        [
            [ 1.1,  2.24937958,   -22.17380450 ],
            [ 1.11, 2.04650122,   -18.64478050 ],
            [ 1.12, 1.87648397,   -15.72303900 ],
            [ 1.13, 1.73204044,   -13.41389050 ],
            [ 1.14, 1.60820616,   -11.59797750 ],
            [ 1.15, 1.50008089,   -10.14617050 ],
            [ 1.16, 1.40528275,   -8.93044450 ],
            [ 1.17, 1.321472,     -7.92274150 ],
            [ 1.18, 1.24682792,   -7.07801750 ],
            [ 1.19, 1.17991165,   -6.36271350 ],
            [ 1.2,  1.11957365,   -5.75154950 ],
            [ 1.21, 1.06488066,   -5.22521600 ],
            [ 1.22, 1.01506933,   -4.76860105 ],
            [ 1.23, 0.969508639,  -4.36984775 ],
            [ 1.24, 0.927672375,  -4.01952270 ],
            [ 1.25, 0.889118185,  -3.71004460 ],
            [ 1.26, 0.853471483,  -3.43526330 ],
            [ 1.27, 0.820412919,  -3.19014725 ],
            [ 1.28, 0.789668538,  -2.97054770 ],
            [ 1.29, 0.761001965,  -2.77301900 ],
            [ 1.3,  0.734208158,  -2.59467980 ],
            [ 1.31, 0.709108369,  -2.43311930 ],
            [ 1.32, 0.685545772,  -2.28625760 ],
            [ 1.33, 0.663383217,  -2.15234935 ],
            [ 1.34, 0.642498785,  -2.02992895 ],
            [ 1.35, 0.622784638,  -1.91769810 ],
            [ 1.36, 0.604144823,  -1.81454855 ],
            [ 1.37, 0.586493667,  -1.71951925 ],
            [ 1.38, 0.569754438,  -1.63177345 ],
            [ 1.39, 0.553858198,  -1.55058015 ],
            [ 1.4,  0.538742835,  -1.47529815 ],
            [ 1.41, 0.524352235,  -1.40536305 ],
            [ 1.42, 0.510635574,  -1.34027675 ],
            [ 1.43, 0.4975467,    -1.27960535 ],
            [ 1.44, 0.485043467,  -1.22294185 ],
            [ 1.45, 0.473087863,  -1.16993765 ],
            [ 1.46, 0.461644714,  -1.12029595 ],
            [ 1.47, 0.450681944,  -1.07372985 ],
            [ 1.48, 0.440170117,  -1.02998880 ],
            [ 1.49, 0.430082168,  -0.98884740 ],
            [ 1.5,  0.420393169,  -0.95010245 ],
            [ 1.51, 0.411080119,  -0.91357025 ],
            [ 1.52, 0.402121764,  -0.87908380 ],
            [ 1.53, 0.393498443,  -0.84649355 ],
            [ 1.54, 0.385191893,  -0.81566240 ],
            [ 1.55, 0.377185195,  -0.78646405 ],
            [ 1.56, 0.369462612,  -0.75878955 ],
            [ 1.57, 0.362009404,  -0.73252595 ],
            [ 1.58, 0.354812093,  -0.70757805 ],
            [ 1.59, 0.347857843,  -0.68386745 ],
            [ 1.6,  0.341134744,  -0.66130930 ],
            [ 1.61, 0.334631657,  -0.63982980 ],
            [ 1.62, 0.328338148,  -0.61936115 ],
            [ 1.63, 0.322244434,  -0.59984065 ],
            [ 1.64, 0.316341335,  -0.58121040 ],
            [ 1.65, 0.310620226,  -0.56341695 ],
            [ 1.66, 0.305072996,  -0.54641090 ],
            [ 1.67, 0.299692008,  -0.53014635 ],
            [ 1.68, 0.294470069,  -0.51458365 ],
            [ 1.69, 0.289400335,  -0.49967725 ],
            [ 1.7,  0.284476524,  -0.48539050 ],
            [ 1.71, 0.279692525,  -0.47169520 ],
            [ 1.72, 0.27504262,   -0.45855640 ],
            [ 1.73, 0.270521397,  -0.44594435 ],
            [ 1.74, 0.266123733,  -0.43383125 ],
            [ 1.75, 0.261844772,  -0.42219115 ],
            [ 1.76, 0.25767991,   -0.41099975 ],
            [ 1.77, 0.253624777,  -0.40023435 ],
            [ 1.78, 0.249675223,  -0.38987355 ],
            [ 1.79, 0.245827306,  -0.37989740 ],
            [ 1.8,  0.242077275,  -0.37028715 ],
            [ 1.81, 0.238421563,  -0.36102725 ],
            [ 1.82, 0.23485673,   -0.35209685 ],
            [ 1.83, 0.231379626,  -0.34348045 ],
            [ 1.84, 0.227987121,  -0.33516760 ],
            [ 1.85, 0.224676274,  -0.32714225 ],
            [ 1.86, 0.221444276,  -0.31939140 ],
            [ 1.87, 0.218288446,  -0.31190265 ],
            [ 1.88, 0.215206223,  -0.30466445 ],
            [ 1.89, 0.212195157,  -0.29766585 ],
            [ 1.9,  0.209252906,  -0.29089635 ],
            [ 1.91, 0.20637723,   -0.28434610 ],
            [ 1.92, 0.203565984,  -0.27800580 ],
            [ 1.93, 0.200817114,  -0.27186830 ],
            [ 1.94, 0.198128618,  -0.26592200 ],
            [ 1.95, 0.195498674,  -0.26015895 ],
            [ 1.96, 0.192925439,  -0.25457485 ],
            [ 1.97, 0.190407177,  -0.24916095 ],
            [ 1.98, 0.18794222,   -0.24391045 ],
            [ 1.99, 0.185528968,  -0.23881685 ],
            [ 2,    0.183165883,  -0.23389773 ],
            [ 2.02, 0.178584359,  -0.22448485 ],
            [ 2.04, 0.174186489,  -0.21556038 ],
            [ 2.06, 0.169961944,  -0.20713350 ],
            [ 2.08, 0.165901149,  -0.19916832 ],
            [ 2.1,  0.161995211,  -0.19163210 ],
            [ 2.12, 0.158235865,  -0.18449592 ],
            [ 2.14, 0.154615374,  -0.17773037 ],
            [ 2.16, 0.15112665,   -0.17131042 ],
            [ 2.18, 0.147762957,  -0.16521520 ],
            [ 2.2,  0.144518042,  -0.15942257 ],
            [ 2.22, 0.141386054,  -0.15391310 ],
            [ 2.24, 0.138361518,  -0.14866888 ],
            [ 2.26, 0.135439299,  -0.14367348 ],
            [ 2.28, 0.132614579,  -0.13891158 ],
            [ 2.3,  0.129882836,  -0.13436913 ],
            [ 2.32, 0.127239814,  -0.13003310 ],
            [ 2.34, 0.124681512,  -0.12589135 ],
            [ 2.36, 0.12220416,   -0.12193337 ],
            [ 2.38, 0.119804177,  -0.11814730 ],
            [ 2.4,  0.117478268,  -0.11452350 ],
            [ 2.42, 0.115223237,  -0.11105433 ],
            [ 2.44, 0.113036095,  -0.10773062 ],
            [ 2.46, 0.110914012,  -0.10454460 ],
            [ 2.48, 0.108854311,  -0.10148882 ],
            [ 2.5,  0.106854459,  -0.09855652 ],
            [ 2.52, 0.10491205,   -0.09574125 ],
            [ 2.54, 0.103024809,  -0.09303700 ],
            [ 2.56, 0.10119057,   -0.09043817 ],
            [ 2.58, 0.0994072823, -0.08793943 ],
            [ 2.6,  0.097672993,  -0.08553588 ],
            [ 2.62, 0.095985847,  -0.08322333 ],
            [ 2.64, 0.0943440599, -0.08099642 ],
            [ 2.66, 0.0927459903, -0.07885105 ],
            [ 2.68, 0.0911900179, -0.07678432 ],
            [ 2.7,  0.0896746174, -0.07479210 ],
            [ 2.72, 0.0881983341, -0.07287093 ],
            [ 2.74, 0.0867597804, -0.07101756 ],
            [ 2.76, 0.0853576317, -0.06922893 ],
            [ 2.78, 0.0839906231, -0.06750215 ],
            [ 2.8,  0.0826575457, -0.06583448 ],
            [ 2.82, 0.081357244,  -0.06422332 ],
            [ 2.84, 0.080088613,  -0.06266622 ],
            [ 2.86, 0.0788505951, -0.06116123 ],
            [ 2.88, 0.0776421638, -0.05970541 ],
            [ 2.9,  0.0764623785, -0.05829671 ],
            [ 2.92, 0.0753102956, -0.05693384 ],
            [ 2.94, 0.0741850249, -0.05561458 ],
            [ 2.96, 0.0730857124, -0.05433714 ],
            [ 2.98, 0.0720115392, -0.05309983 ],
            [ 3,    0.0709617192, -0.05191005 ],
            [ 3.05, 0.0684388348, -0.04909201 ],
            [ 3.1,  0.0660525181, -0.04646047 ],
            [ 3.15, 0.0637927881, -0.04401906 ],
            [ 3.2,  0.0616506126, -0.04175046 ],
            [ 3.25, 0.0596177422, -0.03963947 ],
            [ 3.3,  0.0576866654, -0.03767222 ],
            [ 3.35, 0.0558505204, -0.03583646 ],
            [ 3.4,  0.0541030191, -0.03412099 ],
            [ 3.45, 0.052438421,  -0.03251589 ],
            [ 3.5,  0.0508514304, -0.03101237 ],
            [ 3.55, 0.0493371841, -0.02960228 ],
            [ 3.6,  0.0478912021, -0.02827840 ],
            [ 3.65, 0.0465093445, -0.02703397 ],
            [ 3.7,  0.0451878056, -0.02586299 ],
            [ 3.75, 0.0439230454, -0.02476015 ],
            [ 3.8,  0.0427117905, -0.02372041 ],
            [ 3.85, 0.0415510046, -0.02273928 ],
            [ 3.9,  0.0404378623, -0.02181253 ],
            [ 3.95, 0.039369752,  -0.02093635 ],
            [ 4,    0.0383442271, -0.02011447 ],
            [ 4.1,  0.036411985,  -0.01859768 ],
            [ 4.2,  0.0346246919, -0.01722008 ],
            [ 4.3,  0.0329679691, -0.01597734 ],
            [ 4.4,  0.0314292245, -0.01485310 ],
            [ 4.5,  0.0299973487, -0.01383346 ],
            [ 4.6,  0.0286625332, -0.01290633 ],
            [ 4.7,  0.0274160836, -0.01206134 ],
            [ 4.8,  0.0262502644, -0.01128945 ],
            [ 4.9,  0.0251581945, -0.01058281 ],
            [ 5,    0.0241337021, -0.00993466 ],
            [ 5.1,  0.0231712625, -0.00933896 ],
            [ 5.2,  0.0222659091, -0.00879050 ],
            [ 5.3,  0.0214131635, -0.00828458 ],
            [ 5.4,  0.0206089938, -0.00781712 ],
            [ 5.5,  0.0198497397, -0.00738454 ],
            [ 5.6,  0.0191320859, -0.00698360 ],
            [ 5.7,  0.0184530204, -0.00661144 ],
            [ 5.8,  0.017809798,  -0.00626548 ],
            [ 5.9,  0.0171999241, -0.00594344 ],
            [ 6,    0.0166211108, -0.00564328 ],
            [ 6.1,  0.0160712674, -0.00536317 ],
            [ 6.2,  0.0155484775, -0.00510144 ],
            [ 6.3,  0.0150509792, -0.00485659 ],
            [ 6.4,  0.0145771597, -0.00462726 ],
            [ 6.5,  0.0141255272, -0.00441226 ],
            [ 6.6,  0.0136947081, -0.00421046 ],
            [ 6.7,  0.0132834346, -0.00402088 ],
            [ 6.8,  0.0128905323, -0.00384257 ],
            [ 6.9,  0.01251492,   -0.00367471 ],
            [ 7,    0.0121555899, -0.00351525 ],
        ],
        [
            [ 1.1,  3.99653115,   -39.36364617 ],
            [ 1.11, 3.63788337,   -32.88747050 ],
            [ 1.12, 3.33878174,   -27.71519650 ],
            [ 1.13, 3.08357944,   -23.71600950 ],
            [ 1.14, 2.86446155,   -20.46737800 ],
            [ 1.15, 2.67423188,   -17.84854400 ],
            [ 1.16, 2.50749067,   -15.70620500 ],
            [ 1.17, 2.36010778,   -13.93103450 ],
            [ 1.18, 2.22886998,   -12.44346600 ],
            [ 1.19, 2.11123846,   -11.18423200 ],
            [ 1.2,  2.00518534,   -10.10868750 ],
            [ 1.21, 1.90906471,   -9.18273750 ],
            [ 1.22, 1.82153059,   -8.37970400 ],
            [ 1.23, 1.74147063,   -7.67865650 ],
            [ 1.24, 1.66795746,   -7.06294650 ],
            [ 1.25, 1.6002117,    -6.51919450 ],
            [ 1.26, 1.53757357,   -6.03654900 ],
            [ 1.27, 1.47948072,   -5.60613600 ],
            [ 1.28, 1.42545085,   -5.22063700 ],
            [ 1.29, 1.37506798,   -4.87397600 ],
            [ 1.3,  1.32797133,   -4.56107400 ],
            [ 1.31, 1.2838465,    -4.27768100 ],
            [ 1.32, 1.24241771,   -4.02013200 ],
            [ 1.33, 1.20344386,   -3.78535150 ],
            [ 1.34, 1.16671068,   -3.57075750 ],
            [ 1.35, 1.13202871,   -3.37406450 ],
            [ 1.36, 1.09922939,   -3.19332050 ],
            [ 1.37, 1.0681623,    -3.02683250 ],
            [ 1.38, 1.03869274,   -2.87313200 ],
            [ 1.39, 1.01069966,   -2.73093110 ],
            [ 1.4,  0.984074118,  -2.59909935 ],
            [ 1.41, 0.958717673,  -2.47664625 ],
            [ 1.42, 0.934541193,  -2.36269555 ],
            [ 1.43, 0.911463762,  -2.25648470 ],
            [ 1.44, 0.889411499,  -2.15729865 ],
            [ 1.45, 0.868317789,  -2.06452465 ],
            [ 1.46, 0.848121006,  -1.97764140 ],
            [ 1.47, 0.828764961,  -1.89614520 ],
            [ 1.48, 0.810198102,  -1.81959610 ],
            [ 1.49, 0.792373039,  -1.74759860 ],
            [ 1.5,  0.77524613,   -1.67979585 ],
            [ 1.51, 0.758777122,  -1.61586565 ],
            [ 1.52, 0.742928817,  -1.55551635 ],
            [ 1.53, 0.727666795,  -1.49848325 ],
            [ 1.54, 0.712959152,  -1.44452590 ],
            [ 1.55, 0.698776277,  -1.39342540 ],
            [ 1.56, 0.685090644,  -1.34498980 ],
            [ 1.57, 0.671876481,  -1.29902080 ],
            [ 1.58, 0.659110228,  -1.25535180 ],
            [ 1.59, 0.646769445,  -1.21384545 ],
            [ 1.6,  0.634833319,  -1.17435325 ],
            [ 1.61, 0.62328238,   -1.13674610 ],
            [ 1.62, 0.612098397,  -1.10090505 ],
            [ 1.63, 0.601264279,  -1.06672050 ],
            [ 1.64, 0.590763987,  -1.03409125 ],
            [ 1.65, 0.580582454,  -1.00292370 ],
            [ 1.66, 0.570705513,  -0.97313140 ],
            [ 1.67, 0.561119826,  -0.94463410 ],
            [ 1.68, 0.551812831,  -0.91736255 ],
            [ 1.69, 0.542772575,  -0.89123705 ],
            [ 1.7,  0.53398809,   -0.86619345 ],
            [ 1.71, 0.525448706,  -0.84218265 ],
            [ 1.72, 0.517144437,  -0.81914350 ],
            [ 1.73, 0.509065836,  -0.79702395 ],
            [ 1.74, 0.501203958,  -0.77577550 ],
            [ 1.75, 0.493550326,  -0.75535290 ],
            [ 1.76, 0.4860969,    -0.73571365 ],
            [ 1.77, 0.478836053,  -0.71681795 ],
            [ 1.78, 0.471760541,  -0.69862870 ],
            [ 1.79, 0.464863479,  -0.68111095 ],
            [ 1.8,  0.458138322,  -0.66423195 ],
            [ 1.81, 0.45157884,   -0.64796470 ],
            [ 1.82, 0.445179028,  -0.63227260 ],
            [ 1.83, 0.438933388,  -0.61712770 ],
            [ 1.84, 0.432836474,  -0.60251375 ],
            [ 1.85, 0.426883113,  -0.58840260 ],
            [ 1.86, 0.421068422,  -0.57476985 ],
            [ 1.87, 0.415387716,  -0.56159490 ],
            [ 1.88, 0.409836524,  -0.54885730 ],
            [ 1.89, 0.40441057,   -0.53653800 ],
            [ 1.9,  0.399105764,  -0.52461890 ],
            [ 1.91, 0.393918192,  -0.51308280 ],
            [ 1.92, 0.388844108,  -0.50191340 ],
            [ 1.93, 0.383879924,  -0.49109825 ],
            [ 1.94, 0.379022143,  -0.48061695 ],
            [ 1.95, 0.374267585,  -0.47045575 ],
            [ 1.96, 0.369613028,  -0.46060730 ],
            [ 1.97, 0.365055439,  -0.45105620 ],
            [ 1.98, 0.360591904,  -0.44179060 ],
            [ 1.99, 0.356219627,  -0.43279930 ],
            [ 2,    0.351935918,  -0.42411302 ],
            [ 2.02, 0.343623973,  -0.40748383 ],
            [ 2.04, 0.335636565,  -0.39170852 ],
            [ 2.06, 0.327955632,  -0.37680352 ],
            [ 2.08, 0.320564424,  -0.36270615 ],
            [ 2.1,  0.313447386,  -0.34935925 ],
            [ 2.12, 0.306590054,  -0.33671265 ],
            [ 2.14, 0.29997888,   -0.32471480 ],
            [ 2.16, 0.293601462,  -0.31332210 ],
            [ 2.18, 0.287445996,  -0.30249842 ],
            [ 2.2,  0.281501525,  -0.29220495 ],
            [ 2.22, 0.275757798,  -0.28240782 ],
            [ 2.24, 0.270205212,  -0.27307593 ],
            [ 2.26, 0.264834761,  -0.26418045 ],
            [ 2.28, 0.259637994,  -0.25569480 ],
            [ 2.3,  0.254606969,  -0.24759443 ],
            [ 2.32, 0.249734217,  -0.23985655 ],
            [ 2.34, 0.245012707,  -0.23246005 ],
            [ 2.36, 0.240435815,  -0.22538665 ],
            [ 2.38, 0.235997241,  -0.21861550 ],
            [ 2.4,  0.231691195,  -0.21212967 ],
            [ 2.42, 0.227512054,  -0.20591610 ],
            [ 2.44, 0.223454551,  -0.19995870 ],
            [ 2.46, 0.219513706,  -0.19424375 ],
            [ 2.48, 0.215684801,  -0.18875842 ],
            [ 2.5,  0.211963369,  -0.18349077 ],
            [ 2.52, 0.20834517,   -0.17842957 ],
            [ 2.54, 0.204826186,  -0.17356430 ],
            [ 2.56, 0.201402598,  -0.16888515 ],
            [ 2.58, 0.19807078,   -0.16438287 ],
            [ 2.6,  0.194827283,  -0.16004883 ],
            [ 2.62, 0.191668827,  -0.15587572 ],
            [ 2.64, 0.188592254,  -0.15185407 ],
            [ 2.66, 0.185594664,  -0.14797672 ],
            [ 2.68, 0.182673185,  -0.14423873 ],
            [ 2.7,  0.179825115,  -0.14063255 ],
            [ 2.72, 0.177047883,  -0.13715252 ],
            [ 2.74, 0.174339014,  -0.13379295 ],
            [ 2.76, 0.171696165,  -0.13054808 ],
            [ 2.78, 0.169117091,  -0.12741303 ],
            [ 2.8,  0.166599644,  -0.12438298 ],
            [ 2.82, 0.164141772,  -0.12145340 ],
            [ 2.84, 0.161741508,  -0.11862000 ],
            [ 2.86, 0.159396972,  -0.11587930 ],
            [ 2.88, 0.157106336,  -0.11322615 ],
            [ 2.9,  0.154867926,  -0.11065690 ],
            [ 2.92, 0.15268006,   -0.10816940 ],
            [ 2.94, 0.15054115,   -0.10575965 ],
            [ 2.96, 0.148449674,  -0.10342450 ],
            [ 2.98, 0.14640417,   -0.10116102 ],
            [ 3,    0.144403233,  -0.09898250 ],
            [ 3.05, 0.139587152,  -0.09381583 ],
            [ 3.1,  0.13502165,   -0.08898305 ],
            [ 3.15, 0.130688847,  -0.08449106 ],
            [ 3.2,  0.126572544,  -0.08030926 ],
            [ 3.25, 0.122657921,  -0.07641080 ],
            [ 3.3,  0.118931464,  -0.07277109 ],
            [ 3.35, 0.115380812,  -0.06936845 ],
            [ 3.4,  0.111994619,  -0.06618299 ],
            [ 3.45, 0.108762513,  -0.06319707 ],
            [ 3.5,  0.105674912,  -0.06039511 ],
            [ 3.55, 0.102723002,  -0.05776259 ],
            [ 3.6,  0.0998986533, -0.05528661 ],
            [ 3.65, 0.0971943407, -0.05295458 ],
            [ 3.7,  0.0946031951, -0.05075694 ],
            [ 3.75, 0.0921186468, -0.04868463 ],
            [ 3.8,  0.0897347319, -0.04672647 ],
            [ 3.85, 0.087446,     -0.04487506 ],
            [ 3.9,  0.0852472257, -0.04312378 ],
            [ 3.95, 0.0831336218, -0.04146530 ],
            [ 4,    0.0811006953, -0.03990662 ],
            [ 4.1,  0.0772604149, -0.03702306 ],
            [ 4.2,  0.073696084,  -0.03439642 ],
            [ 4.3,  0.0703811308, -0.03201931 ],
            [ 4.4,  0.0672922223, -0.02986209 ],
            [ 4.5,  0.0644087137, -0.02789951 ],
            [ 4.6,  0.0617123207, -0.02610958 ],
            [ 4.7,  0.0591867986, -0.02447335 ],
            [ 4.8,  0.0568176508, -0.02297426 ],
            [ 4.9,  0.0545919475, -0.02159793 ],
            [ 5,    0.0524980647, -0.02033192 ],
            [ 5.1,  0.0505255638, -0.01916513 ],
            [ 5.2,  0.0486650383, -0.01808790 ],
            [ 5.3,  0.0469079841, -0.01709156 ],
            [ 5.4,  0.0452467264, -0.01616852 ],
            [ 5.5,  0.0436742792, -0.01531213 ],
            [ 5.6,  0.0421842995, -0.01451635 ],
            [ 5.7,  0.04077101,   -0.01377584 ],
            [ 5.8,  0.0394291322, -0.01308575 ],
            [ 5.9,  0.0381538595, -0.01244180 ],
            [ 6,    0.0369407716, -0.01184020 ],
            [ 6.1,  0.0357858201, -0.01127743 ],
            [ 6.2,  0.0346852854, -0.01075040 ],
            [ 6.3,  0.0336357403, -0.01025622 ],
            [ 6.4,  0.032634041,  -0.00979234 ],
            [ 6.5,  0.031677272,  -0.00935649 ],
            [ 6.6,  0.030762744,  -0.00894652 ],
            [ 6.7,  0.0298879686, -0.00856054 ],
            [ 6.8,  0.0290506363, -0.00819676 ],
            [ 6.9,  0.0282486166, -0.00785359 ],
            [ 7,    0.027479919,  -0.00752700 ],
        ],
        [
            [ 1.1,  3.41668546,   -33.55360867 ],
            [ 1.11, 3.11077043,   -28.05790200 ],
            [ 1.12, 2.85552742,   -23.63345700 ],
            [ 1.13, 2.63810129,   -20.20395550 ],
            [ 1.14, 2.45144831,   -17.43368950 ],
            [ 1.15, 2.2894275,    -15.20080050 ],
            [ 1.16, 2.1474323,    -13.37440500 ],
            [ 1.17, 2.0219394,    -11.86121900 ],
            [ 1.18, 1.91020792,   -10.59335300 ],
            [ 1.19, 1.81007234,   -9.52024050 ],
            [ 1.2,  1.71980311,   -8.60379050 ],
            [ 1.21, 1.63799653,   -7.81491550 ],
            [ 1.22, 1.5635048,    -7.13085450 ],
            [ 1.23, 1.49537944,   -6.53375250 ],
            [ 1.24, 1.43282975,   -6.00940750 ],
            [ 1.25, 1.37519129,   -5.54640700 ],
            [ 1.26, 1.32190161,   -5.13549450 ],
            [ 1.27, 1.2724814,    -4.76910100 ],
            [ 1.28, 1.22651959,   -4.44098550 ],
            [ 1.29, 1.18366169,   -4.14596600 ],
            [ 1.3,  1.14360027,   -3.87971100 ],
            [ 1.31, 1.10606747,   -3.63859650 ],
            [ 1.32, 1.07082834,   -3.41949850 ],
            [ 1.33, 1.0376775,    -3.21979200 ],
            [ 1.34, 1.0064325,    -3.03727995 ],
            [ 1.35, 0.976931901,  -2.87001400 ],
            [ 1.36, 0.94903222,   -2.71632665 ],
            [ 1.37, 0.922605368,  -2.57477800 ],
            [ 1.38, 0.89753666,   -2.44411330 ],
            [ 1.39, 0.873723102,  -2.32323625 ],
            [ 1.4,  0.851071935,  -2.21118575 ],
            [ 1.41, 0.829499387,  -2.10711640 ],
            [ 1.42, 0.808929607,  -2.01028185 ],
            [ 1.43, 0.78929375,   -1.92003235 ],
            [ 1.44, 0.77052896,   -1.83575900 ],
            [ 1.45, 0.75257857,   -1.75693990 ],
            [ 1.46, 0.735390162,  -1.68313105 ],
            [ 1.47, 0.718915949,  -1.61390350 ],
            [ 1.48, 0.703112092,  -1.54888270 ],
            [ 1.49, 0.687938295,  -1.48773190 ],
            [ 1.5,  0.673357454,  -1.43014740 ],
            [ 1.51, 0.659335347,  -1.37585485 ],
            [ 1.52, 0.645840357,  -1.32460605 ],
            [ 1.53, 0.632843226,  -1.27617580 ],
            [ 1.54, 0.620316841,  -1.23035940 ],
            [ 1.55, 0.608236038,  -1.18697055 ],
            [ 1.56, 0.59657743,   -1.14584600 ],
            [ 1.57, 0.585319118,  -1.10681710 ],
            [ 1.58, 0.574441088,  -1.06974205 ],
            [ 1.59, 0.563924277,  -1.03450410 ],
            [ 1.6,  0.553751006,  -1.00097690 ],
            [ 1.61, 0.543904739,  -0.96905065 ],
            [ 1.62, 0.534369993,  -0.93862430 ],
            [ 1.63, 0.525132253,  -0.90960460 ],
            [ 1.64, 0.516177901,  -0.88190545 ],
            [ 1.65, 0.507494144,  -0.85544745 ],
            [ 1.66, 0.499068952,  -0.83015695 ],
            [ 1.67, 0.490891005,  -0.80596490 ],
            [ 1.68, 0.482949654,  -0.78281445 ],
            [ 1.69, 0.475234716,  -0.76063750 ],
            [ 1.7,  0.467736904,  -0.73937785 ],
            [ 1.71, 0.460447159,  -0.71899485 ],
            [ 1.72, 0.453357007,  -0.69943640 ],
            [ 1.73, 0.446458431,  -0.68065840 ],
            [ 1.74, 0.439743839,  -0.66261960 ],
            [ 1.75, 0.433206039,  -0.64528145 ],
            [ 1.76, 0.42683821,   -0.62860795 ],
            [ 1.77, 0.42063388,   -0.61256540 ],
            [ 1.78, 0.414586902,  -0.59712220 ],
            [ 1.79, 0.408691436,  -0.58224870 ],
            [ 1.8,  0.402941928,  -0.56791700 ],
            [ 1.81, 0.397333096,  -0.55410425 ],
            [ 1.82, 0.391859843,  -0.54077945 ],
            [ 1.83, 0.386517507,  -0.52791945 ],
            [ 1.84, 0.381301454,  -0.51550905 ],
            [ 1.85, 0.376207326,  -0.50352450 ],
            [ 1.86, 0.371230964,  -0.49194640 ],
            [ 1.87, 0.366368398,  -0.48075650 ],
            [ 1.88, 0.361615834,  -0.46993765 ],
            [ 1.89, 0.356969645,  -0.45947355 ],
            [ 1.9,  0.352426363,  -0.44934875 ],
            [ 1.91, 0.34798267,   -0.43954880 ],
            [ 1.92, 0.343635387,  -0.43005985 ],
            [ 1.93, 0.339381473,  -0.42087130 ],
            [ 1.94, 0.335217961,  -0.41196590 ],
            [ 1.95, 0.331142155,  -0.40333185 ],
            [ 1.96, 0.327151324,  -0.39496300 ],
            [ 1.97, 0.323242895,  -0.38684640 ],
            [ 1.98, 0.319414396,  -0.37897195 ],
            [ 1.99, 0.315663456,  -0.37132995 ],
            [ 2,    0.311987797,  -0.36394660 ],
            [ 2.02, 0.304853637,  -0.34981048 ],
            [ 2.04, 0.297995378,  -0.33639833 ],
            [ 2.06, 0.291397704,  -0.32372415 ],
            [ 2.08, 0.285046412,  -0.31173482 ],
            [ 2.1,  0.278928311,  -0.30038195 ],
            [ 2.12, 0.273031134,  -0.28962273 ],
            [ 2.14, 0.267343402,  -0.27941385 ],
            [ 2.16, 0.26185458,   -0.26971842 ],
            [ 2.18, 0.256554665,  -0.26050535 ],
            [ 2.2,  0.251434366,  -0.25174195 ],
            [ 2.22, 0.246484987,  -0.24339957 ],
            [ 2.24, 0.241698383,  -0.23545178 ],
            [ 2.26, 0.237066916,  -0.22787420 ],
            [ 2.28, 0.232583415,  -0.22064430 ],
            [ 2.3,  0.228241144,  -0.21374120 ],
            [ 2.32, 0.224033767,  -0.20714570 ],
            [ 2.34, 0.219955316,  -0.20083987 ],
            [ 2.36, 0.216000172,  -0.19480820 ],
            [ 2.38, 0.212162988,  -0.18903303 ],
            [ 2.4,  0.208438851,  -0.18350000 ],
            [ 2.42, 0.204822988,  -0.17819805 ],
            [ 2.44, 0.201310929,  -0.17311357 ],
            [ 2.46, 0.197898445,  -0.16823490 ],
            [ 2.48, 0.194581533,  -0.16355120 ],
            [ 2.5,  0.191356397,  -0.15905232 ],
            [ 2.52, 0.18821944,   -0.15472877 ],
            [ 2.54, 0.185167246,  -0.15057163 ],
            [ 2.56, 0.182196575,  -0.14657255 ],
            [ 2.58, 0.179304344,  -0.14272372 ],
            [ 2.6,  0.176487626,  -0.13901780 ],
            [ 2.62, 0.173743632,  -0.13544870 ],
            [ 2.64, 0.171069678,  -0.13200827 ],
            [ 2.66, 0.168463301,  -0.12869042 ],
            [ 2.68, 0.165922061,  -0.12549105 ],
            [ 2.7,  0.163443659,  -0.12240260 ],
            [ 2.72, 0.161025957,  -0.11942265 ],
            [ 2.74, 0.158666753,  -0.11654727 ],
            [ 2.76, 0.156364066,  -0.11376710 ],
            [ 2.78, 0.154116069,  -0.11107935 ],
            [ 2.8,  0.151920892,  -0.10848203 ],
            [ 2.82, 0.149776788,  -0.10597017 ],
            [ 2.84, 0.147682085,  -0.10354012 ],
            [ 2.86, 0.145635183,  -0.10118903 ],
            [ 2.88, 0.143634524,  -0.09891240 ],
            [ 2.9,  0.141678687,  -0.09670717 ],
            [ 2.92, 0.139766237,  -0.09457160 ],
            [ 2.94, 0.137895823,  -0.09250220 ],
            [ 2.96, 0.136066149,  -0.09049638 ],
            [ 2.98, 0.134275968,  -0.08855160 ],
            [ 3,    0.132524085,  -0.08667916 ],
            [ 3.05, 0.128304501,  -0.08223641 ],
            [ 3.1,  0.124300444,  -0.07807830 ],
            [ 3.15, 0.120496671,  -0.07421078 ],
            [ 3.2,  0.116879366,  -0.07060789 ],
            [ 3.25, 0.113435882,  -0.06724682 ],
            [ 3.3,  0.110154684,  -0.06410668 ],
            [ 3.35, 0.107025214,  -0.06116906 ],
            [ 3.4,  0.104037778,  -0.05841703 ],
            [ 3.45, 0.101183511,  -0.05583559 ],
            [ 3.5,  0.0984542193, -0.05341150 ],
            [ 3.55, 0.0958423613, -0.05113241 ],
            [ 3.6,  0.0933409782, -0.04898736 ],
            [ 3.65, 0.0909436256, -0.04696607 ],
            [ 3.7,  0.088644371,  -0.04505943 ],
            [ 3.75, 0.0864376822, -0.04325936 ],
            [ 3.8,  0.0843184349, -0.04155816 ],
            [ 3.85, 0.0822818665, -0.03994899 ],
            [ 3.9,  0.080323536,  -0.03842534 ],
            [ 3.95, 0.0784393328, -0.03698140 ],
            [ 4,    0.0766253961, -0.03562321 ],
            [ 4.1,  0.0731941797, -0.03310799 ],
            [ 4.2,  0.0700037985, -0.03081395 ],
            [ 4.3,  0.067031389,  -0.02873482 ],
            [ 4.4,  0.0642568343, -0.02684535 ],
            [ 4.5,  0.0616623193, -0.02512369 ],
            [ 4.6,  0.0592320964, -0.02355160 ],
            [ 4.7,  0.0569520001, -0.02211246 ],
            [ 4.8,  0.0548096037, -0.02079166 ],
            [ 4.9,  0.0527936682, -0.01957745 ],
            [ 5,    0.050894114,  -0.01845891 ],
            [ 5.1,  0.0491018861, -0.01742656 ],
            [ 5.2,  0.0474088028, -0.01647206 ],
            [ 5.3,  0.0458074744, -0.01558790 ],
            [ 5.4,  0.0442912227, -0.01476756 ],
            [ 5.5,  0.042853962,  -0.01400530 ],
            [ 5.6,  0.041490163,  -0.01329588 ],
            [ 5.7,  0.0401947856, -0.01263471 ],
            [ 5.8,  0.0389632209, -0.01201756 ],
            [ 5.9,  0.0377912733, -0.01144071 ],
            [ 6,    0.0366750787, -0.01090087 ],
            [ 6.1,  0.0356110989, -0.01039498 ],
            [ 6.2,  0.034596083,  -0.00992031 ],
            [ 6.3,  0.0336270366, -0.00947434 ],
            [ 6.4,  0.0327012159, -0.00905478 ],
            [ 6.5,  0.0318160802, -0.00865962 ],
            [ 6.6,  0.0309692927, -0.00828689 ],
            [ 6.7,  0.0301587021, -0.00793482 ],
            [ 6.8,  0.0293823293, -0.00760160 ],
            [ 6.9,  0.0286383823, -0.00728531 ],
            [ 7,    0.027925268,  -0.00698229 ],
        ],
    ];
}

##################################################################
# BLAST TABLES
##################################################################

BEGIN {

    # name => [symmetry, gamma, error, N]
    # symmetry = 0,1,2 for Plane, Cylindrical, Spherical
    # Error = the sum of all run errors plus log interpolation error;
    # N=number of points from center to shock in Finite Difference run
    $rtables_info = {
        'P1.1'   => [ 0, 1.1,   1.0e-6, 8000 ],
        'P1.15'  => [ 0, 1.15,  1.3e-6, 8000 ],
        'P1.2'   => [ 0, 1.2,   1.e-6,  8000 ],
        'P1.25'  => [ 0, 1.25,  1.3e-6, 8000 ],
        'P1.3'   => [ 0, 1.3,   1.1e-6, 8000 ],
        'P1.35'  => [ 0, 1.35,  1.4e-6, 8000 ],
        'P1.4'   => [ 0, 1.4,   1.1e-6, 8000 ],
        'P1.45'  => [ 0, 1.45,  1.5e-6, 8000 ],
        'P1.5'   => [ 0, 1.5,   1.7e-6, 8000 ],
        'P1.6'   => [ 0, 1.6,   1.7e-6, 8000 ],
        'P1.667' => [ 0, 1.667, 1.1e-6, 8000 ],
        'P1.8'   => [ 0, 1.8,   1.4e-6, 8000 ],
        'P2'     => [ 0, 2,     1.3e-6, 8000 ],
        'P2.5'   => [ 0, 2.5,   2.1e-6, 4000 ],
        'P3'     => [ 0, 3,     2.7e-6, 8000 ],
        'P4'     => [ 0, 4,     3.4e-6, 4000 ],
        'P5'     => [ 0, 5,     1.7e-6, 4000 ],
        'P6'     => [ 0, 6,     1.5e-6, 4000 ],
        'P7'     => [ 0, 7,     2.5e-6, 4000 ],
        'C1.1'   => [ 1, 1.1,   1.8e-6, 8000 ],
        'C1.15'  => [ 1, 1.15,  1.1e-6, 8000 ],
        'C1.2'   => [ 1, 1.2,   1.0e-6, 8000 ],
        'C1.25'  => [ 1, 1.25,  2.0e-6, 8000 ],
        'C1.3'   => [ 1, 1.3,   8.e-7,  8000 ],
        'C1.35'  => [ 1, 1.35,  8.7e-7, 8000 ],
        'C1.4'   => [ 1, 1.4,   8.e-7,  16000 ],
        'C1.45'  => [ 1, 1.45,  7.0e-7, 8000 ],
        'C1.5'   => [ 1, 1.5,   1.2e-6, 8000 ],
        'C1.6'   => [ 1, 1.6,   7.6e-7, 8000 ],
        'C1.667' => [ 1, 1.667, 6.8e-7, 8000 ],
        'C1.8'   => [ 1, 1.8,   6.8e-7, 8000 ],
        'C2'     => [ 1, 2,     8.e-7,  8000 ],
        'C2.5'   => [ 1, 2.5,   9.3e-7, 4000 ],
        'C3'     => [ 1, 3,     1.6e-6, 8000 ],
        'C4'     => [ 1, 4,     8.7e-7, 4000 ],
        'C5'     => [ 1, 5,     9.e-7,  4000 ],
        'C6'     => [ 1, 6,     9.e-7,  4000 ],
        'C7'     => [ 1, 7,     1.6e-6, 4000 ],
        'S1.1'   => [ 2, 1.1,   2.6e-6, 8000 ],
        'S1.15'  => [ 2, 1.15,  1.0e-6, 8000 ],
        'S1.2'   => [ 2, 1.2,   1.1e-6, 8000 ],
        'S1.25'  => [ 2, 1.25,  1.1e-6, 8000 ],
        'S1.3'   => [ 2, 1.3,   8.9e-7, 8000 ],
        'S1.35'  => [ 2, 1.35,  1.e-6,  8000 ],
        'S1.4a'  => [ 2, 1.4,   9.e-7,  16000 ],
        'S1.4'   => [ 2, 1.4,   7.2e-7, 32000 ],
        'S1.45'  => [ 2, 1.45,  9e-7,   8000 ],
        'S1.5'   => [ 2, 1.5,   1.e-6,  8000 ],
        'S1.6'   => [ 2, 1.6,   1.1e-6, 8000 ],
        'S1.667' => [ 2, 1.667, 8.2e-7, 8000 ],
        'S1.8'   => [ 2, 1.8,   7.5e-7, 8000 ],
        'S2'     => [ 2, 2,     8.3e-7, 8000 ],
        'S2.5'   => [ 2, 2.5,   9.8e-7, 4000 ],
        'S3'     => [ 2, 3,     1.5e-6, 8000 ],
        'S4'     => [ 2, 4,     9.1e-7, 4000 ],
        'S5'     => [ 2, 5,     1.1e-6, 4000 ],
        'S6'     => [ 2, 6,     1.0e-6, 4000 ],
        'S6.5'   => [ 2, 6.5,   4.1e-6, 4000 ],
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

        # Run: S8000_G1x667
        'S1.667' => [
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

        # Run: C8000_G1x2
        'C1.2' => [
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
        # Run: C8000_G1x667
        'C1.667' => [
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
        # Run: P8000_G1x667
        'P1.667' => [
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

        # Run: P8000_G1x3
        'P1.3' => [
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
        # Run: S8000_G1x3
        'S1.3' => [
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
        # Run: S8000_G1x2
        'S1.2' => [
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
        # Run: C8000_G1x3
        'C1.3' => [
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
        # Run: C8000_G1x1
        'C1.1' => [
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
        # Run: S8000_G3
        'S3' => [
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

        # Run: P8000_G1x1
        'P1.1' => [
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
        # Run: S8000_G1x1
        'S1.1' => [
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
        # Run: C8000_G3
        'C3' => [
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
        # Run: C8000_G2
        'C2' => [
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
        # Run: S8000_G2
        'S2' => [
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
        # Run: P8000_G2
        'P2' => [
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
        # Run: C4000_G7
        'C7' => [
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
        # Run: P8000_G3
        'P3' => [
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
        # Run: P4000_G7
        'P7' => [
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
        # Run: S16000_G1x4
        'S1.4a' => [
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
        # Run: S32000_G1x4
        'S1.4' => [
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
        # Run: C16000_G1x4
        'C1.4' => [
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
        # Run: P8000_G1x4
        'P1.4' => [
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

        # Run: P8000_G1x2
        'P1.2' => [
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

        # Run: S8000_G1x15
        'S1.15' => [
            [ -4.91112446, 12.00017628,  -2.99998029, -4.91215040, 0.99846031 ],
            [ -4.41115547, 10.50029178,  -2.99991166, -4.41332851, 0.99673698 ],
            [ -3.90334821, 8.97697380,   -2.99960949, -3.90800830, 0.99299442 ],
            [ -3.64524552, 8.20282083,   -2.99913850, -3.65211590, 0.98966167 ],
            [ -3.37986444, 7.40702752,   -2.99809014, -3.39010971, 0.98456213 ],
            [ -3.15882552, 6.74450949,   -2.99629856, -3.17312355, 0.97842338 ],
            [ -2.96387680, 6.16064224,   -2.99337924, -2.98306961, 0.97099172 ],
            [ -2.79087522, 5.64313503,   -2.98893349, -2.81580862, 0.96225842 ],
            [ -2.63601360, 5.18072350,   -2.98252231, -2.66753854, 0.95221782 ],
            [ -2.49828418, 4.77050102,   -2.97384641, -2.53713158, 0.94106227 ],
            [ -2.36704932, 4.38097231,   -2.96178107, -2.41445792, 0.92803832 ],
            [ -2.24413407, 4.01785308,   -2.94578838, -2.30126630, 0.91329964 ],
            [ -2.12611228, 3.67136692,   -2.92470868, -2.19444416, 0.89643456 ],
            [ -2.01077788, 3.33556315,   -2.89714171, -2.09214487, 0.87701050 ],
            [ -1.89219621, 2.99412429,   -2.85993107, -1.98949710, 0.85364785 ],
            [ -1.76331283, 2.62876577,   -2.80743698, -1.88133646, 0.82402295 ],
            [ -1.63857252, 2.28242129,   -2.74338098, -1.78056300, 0.79099548 ],
            [ -1.52467273, 1.97383904,   -2.67330370, -1.69236193, 0.75719566 ],
            [ -1.41796043, 1.69248831,   -2.59838014, -1.61337866, 0.72268508 ],
            [ -1.31608595, 1.43172709,   -2.51989768, -1.54152599, 0.68762831 ],
            [ -1.21732388, 1.18681690,   -2.43908653, -1.47535377, 0.65222129 ],
            [ -1.12010343, 0.95367021,   -2.35688089, -1.41367256, 0.61659633 ],
            [ -1.02335691, 0.72964539,   -2.27430894, -1.35574242, 0.58099385 ],
            [ -0.92605051, 0.51235016,   -2.19217883, -1.30093481, 0.54562045 ],
            [ -0.82738247, 0.30006753,   -2.11129331, -1.24883069, 0.51073439 ],
            [ -0.72638090, 0.09085511,   -2.03215869, -1.19898839, 0.47650992 ],
            [ -0.62237889, -0.11645777,  -1.95540265, -1.15118050, 0.44320645 ],
            [ -0.51443486, -0.32348219,  -1.88135030, -1.10510078, 0.41097436 ],
            [ -0.40176184, -0.53140043,  -1.81037657, -1.06056597, 0.37999808 ],
            [ -0.28355518, -0.74133590,  -1.74278182, -1.01742420, 0.35043463 ],
            [ -0.15881488, -0.95466021,  -1.67871458, -0.97549354, 0.32237379 ],
            [ -0.02646026, -1.17276439,  -1.61827294, -0.93461482, 0.29587924 ],
            [ 0.11455056,  -1.39687231,  -1.56156326, -0.89468435, 0.27101210 ],
            [ 0.26549957,  -1.62849655,  -1.50858512, -0.85556984, 0.24778035 ],
            [ 0.42573599,  -1.86622951,  -1.45988719, -0.81761921, 0.22642649 ],
            [ 0.59255050,  -2.10602439,  -1.41619270, -0.78148611, 0.20726428 ],
            [ 0.76716275,  -2.34978907,  -1.37687503, -0.74683946, 0.19001276 ],
            [ 0.94923432,  -2.59719424,  -1.34169521, -0.71368592, 0.17456032 ],
            [ 1.14239457,  -2.85318933,  -1.30973755, -0.68136101, 0.16049680 ],
            [ 1.34680174,  -3.11787923,  -1.28085719, -0.64989072, 0.14775115 ],
            [ 1.56366298,  -3.39274040,  -1.25474690, -0.61913638, 0.13618177 ],
            [ 1.79590255,  -3.68130722,  -1.23098217, -0.58877132, 0.12559480 ],
            [ 2.04220877,  -3.98180088,  -1.20960761, -0.55904835, 0.11600691 ],
            [ 2.30620812,  -4.29850517,  -1.19021586, -0.52961128, 0.10723354 ],
            [ 2.58893609,  -4.63246308,  -1.17267864, -0.50045668, 0.09921634 ],
            [ 2.89230559,  -4.98574365,  -1.15682090, -0.47150019, 0.09187704 ],
            [ 3.21818422,  -5.36032475,  -1.14249402, -0.44268439, 0.08515070 ],
            [ 3.56883154,  -5.75860170,  -1.12955178, -0.41393794, 0.07897452 ],
            [ 3.94689869,  -6.18337196,  -1.11785643, -0.38518306, 0.07329052 ],
            [ 4.35502826,  -6.63738196,  -1.10729098, -0.35636827, 0.06805139 ],
            [ 4.79625551,  -7.12377965,  -1.09774424, -0.32743747, 0.06321335 ],
            [ 5.27420893,  -7.64632465,  -1.08910975, -0.29832262, 0.05873518 ],
            [ 5.79271064,  -8.20894264,  -1.08129484, -0.26897299, 0.05458270 ],
            [ 6.35637747,  -8.81637563,  -1.07421072, -0.23932261, 0.05072340 ],
            [ 6.97002167,  -9.47352645,  -1.06778185, -0.20932742, 0.04713122 ],
            [ 7.63885184,  -10.18567870, -1.06194110, -0.17895327, 0.04378393 ],
            [ 8.36947377,  -10.95955427, -1.05662236, -0.14813540, 0.04065857 ],
            [ 9.16849112,  -11.80181782, -1.05177297, -0.11684575, 0.03773849 ],
            [ 10.04290611, -12.71951050, -1.04734752, -0.08507107, 0.03500992 ],
            [ 11.00031995, -13.72026228, -1.04330576, -0.05280477, 0.03246065 ],
            [ 12.04933330, -14.81270753, -1.03961090, -0.02003556, 0.03007902 ],
            [ 13.19894647, -16.00585703, -1.03623172, 0.01323128,  0.02785535 ],
            [ 14.45875959, -17.30930973, -1.03314074, 0.04698198,  0.02578092 ],
            [ 15.83917271, -18.73346049, -1.03031316, 0.08119977,  0.02384750 ],
            [ 17.24307967, -20.17818086, -1.02789702, 0.11346982,  0.02216680 ],
            [ 18.88624914, -21.86517869, -1.02552029, 0.14847547,  0.02048641 ],
            [ 20.60571292, -23.62667892, -1.02343433, 0.18237941,  0.01898827 ],
            [ 22.55449584, -25.61913934, -1.02145004, 0.21793561,  0.01754182 ],
            [ 24.57765276, -27.68389236, -1.01971884, 0.25209712,  0.01626198 ],
            [ 26.86743409, -30.01688264, -1.01807011, 0.28788109,  0.01502683 ],
            [ 29.22047725, -32.41071861, -1.01664169, 0.32193231,  0.01394322 ],
            [ 31.87710798, -35.10969870, -1.01527930, 0.35754812,  0.01289743 ],
            [ 34.89073161, -38.16735547, -1.01398196, 0.39485576,  0.01188988 ],
            [ 37.96212927, -41.27994091, -1.01286915, 0.42999750,  0.01101610 ],
            [ 41.43413882, -44.79473628, -1.01180756, 0.46674383,  0.01017390 ],
            [ 45.37825464, -48.78336889, -1.01079648, 0.50522768,  0.00936355 ],
            [ 49.35121853, -52.79747762, -1.00993962, 0.54101540,  0.00867020 ],
            [ 53.84232368, -57.33133318, -1.00912167, 0.57841382,  0.00800242 ],
            [ 58.94426259, -62.47776244, -1.00834216, 0.61755734,  0.00736037 ],
            [ 63.99871275, -67.57268809, -1.00769128, 0.65335853,  0.00681982 ],
            [ 69.70162667, -73.31764086, -1.00706918, 0.69072697,  0.00629922 ],
            [ 76.16728688, -79.82703268, -1.00647555, 0.72979205,  0.00579867 ],
            [ 82.42199174, -86.12067948, -1.00598917, 0.76473715,  0.00538567 ],
            [ 89.44753400, -93.18661818, -1.00552330, 0.80113957,  0.00498754 ],
            [ 97.37442141, -101.15547270, -1.00507776, 0.83911545, 0.00460434 ],
            [
                106.36186551, -110.18658539, -1.00465235, 0.87879554,
                0.00423614
            ],
            [
                114.80149468, -118.66401092, -1.00431310, 0.91327024,
                0.00394079
            ],
            [
                124.25053910, -128.15223185, -1.00398759, 0.94912698,
                0.00365592
            ],
            [
                134.87497442, -138.81733136, -1.00367572, 0.98647326,
                0.00338156
            ],
            [
                146.87620919, -150.86083866, -1.00337739, 1.02542973,
                0.00311774
            ],
            [
                160.50054627, -164.52919267, -1.00309249, 1.06613241,
                0.00286448
            ],
            [
                172.77043062, -176.83564942, -1.00287417, 1.10005472,
                0.00266951
            ],
            [
                186.46084482, -190.56393935, -1.00266431, 1.13528151,
                0.00248132
            ],
            [
                201.79802500, -205.94039736, -1.00246288, 1.17191152,
                0.00229994
            ],
        ],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 1.99e-07
#    Est Max FD Error (based on N/2 run) to R=10.00: 3.e-7
#    Est Max MOC error for R>10.00: 3.19e-07 
#    Max interpolation error with cubic splines: 5.0539e-07
#    Est Max overall error after cubic interpolation: 1.1e-6
#    
        # Run: S8000_G1x15
        'S1.25' => [
[-4.75637566, 12.00017662, -2.99997952, -4.75742125, 0.99843081],
[-4.60173799, 11.53626676, -2.99996744, -4.60305673, 0.99802062],
[-4.15182189, 10.18654836, -2.99987782, -4.15441321, 0.99610812],
[-3.92899136, 9.51809462, -2.99977702, -3.93261294, 0.99455818],
[-3.68880469, 8.79761671, -2.99952271, -3.69400095, 0.99218649],
[-3.39296638, 7.91033005, -2.99883675, -3.40107595, 0.98779060],
[-3.15289168, 7.19051489, -2.99760912, -3.16453497, 0.98244603],
[-2.92770968, 6.51573610, -2.99532001, -2.94406405, 0.97530276],
[-2.74261477, 5.96160747, -2.99187610, -2.76424797, 0.96727969],
[-2.57711552, 5.46684574, -2.98673080, -2.60490708, 0.95790307],
[-2.42886044, 5.02454782, -2.97946807, -2.46365430, 0.94723189],
[-2.29295524, 4.62025773, -2.96948268, -2.33571708, 0.93509424],
[-2.16552753, 4.24267429, -2.95596482, -2.21741592, 0.92122319],
[-2.04493031, 3.88721437, -2.93807448, -2.10724174, 0.90545693],
[-1.92836419, 3.54603773, -2.91455895, -2.00272053, 0.88738452],
[-1.81291782, 3.21126719, -2.88361911, -1.90145661, 0.86637298],
[-1.69344035, 2.86912818, -2.84184281, -1.79942076, 0.84102759],
[-1.56039281, 2.49485677, -2.78178174, -1.68965265, 0.80821098],
[-1.44102731, 2.16668036, -2.71481946, -1.59514701, 0.77461218],
[-1.33086668, 1.87152042, -2.64224908, -1.51167640, 0.74032831],
[-1.22695375, 1.60088576, -2.56539344, -1.43653752, 0.70550006],
[-1.12760577, 1.34992763, -2.48588117, -1.36817686, 0.67044798],
[-1.03029937, 1.11198926, -2.40414325, -1.30465500, 0.63502807],
[-0.93388878, 0.88419010, -2.32131274, -1.24514413, 0.59947656],
[-0.83746024, 0.66435075, -2.23846118, -1.18904843, 0.56406242],
[-0.74004187, 0.45030053, -2.15638400, -1.13581441, 0.52899452],
[-0.64085480, 0.24043899, -2.07585184, -1.08506667, 0.49452241],
[-0.53906940, 0.03317799, -1.99743890, -1.03646112, 0.46085312],
[-0.43378674, -0.17307441, -1.92155687, -0.98968256, 0.42815430],
[-0.32420216, -0.37959339, -1.84860867, -0.94451554, 0.39660988],
[-0.20948259, -0.58760247, -1.77891103, -0.90077796, 0.36637788],
[-0.08877339, -0.79826641, -1.71271465, -0.85832069, 0.33759194],
[0.03899024, -1.01301221, -1.65012537, -0.81696452, 0.31032170],
[0.17481608, -1.23305983, -1.59126783, -0.77659556, 0.28464068],
[0.31991948, -1.45986942, -1.53616474, -0.73707888, 0.26057244],
[0.47567690, -1.69504065, -1.48479075, -0.69828397, 0.23811331],
[0.63812438, -1.93239132, -1.43853955, -0.66128736, 0.21787524],
[0.80736181, -2.17224829, -1.39706604, -0.62599047, 0.19970651],
[0.98484241, -2.41680123, -1.35971635, -0.59203751, 0.18331769],
[1.17184402, -2.66783626, -1.32600257, -0.55917857, 0.16849019],
[1.36964466, -2.92702668, -1.29552466, -0.52721523, 0.15504345],
[1.57956344, -3.19600956, -1.26794606, -0.49598448, 0.14282458],
[1.80273303, -3.47611792, -1.24300422, -0.46538255, 0.13171375],
[2.04041716, -3.76880682, -1.22045442, -0.43530974, 0.12159980],
[2.29454577, -4.07629209, -1.20003294, -0.40561227, 0.11236337],
[2.56629292, -4.39981912, -1.18157615, -0.37625377, 0.10393110],
[2.85786131, -4.74182489, -1.16487358, -0.34710623, 0.09620906],
[3.17090066, -5.10404545, -1.14977961, -0.31812566, 0.08913392],
[3.50747520, -5.48867052, -1.13614483, -0.28924685, 0.08264166],
[3.87024115, -5.89852376, -1.12381974, -0.26037922, 0.07666877],
[4.26164548, -6.33614784, -1.11268488, -0.23147581, 0.07116680],
[4.68452780, -6.80449195, -1.10262510, -0.20248226, 0.06609001],
[5.14212050, -7.30690225, -1.09353250, -0.17334191, 0.06139663],
[5.63824938, -7.84733273, -1.08530494, -0.14398828, 0.05704776],
[6.17733433, -8.43033371, -1.07784890, -0.11435189, 0.05300853],
[6.76359008, -9.06018723, -1.07108900, -0.08440484, 0.04925332],
[7.40242729, -9.74241795, -1.06494880, -0.05408657, 0.04575513],
[8.09965328, -10.48292153, -1.05936280, -0.02335231, 0.04249187],
[8.86167284, -11.28817880, -1.05427330, 0.00783581, 0.03944461],
[9.69528893, -12.16504445, -1.04963091, 0.03949964, 0.03659783],
[10.60790329, -13.12096006, -1.04539212, 0.07165354, 0.03393802],
[11.60771686, -14.16416220, -1.04151811, 0.10430988, 0.03145292],
[12.70313011, -15.30305665, -1.03797639, 0.13745907, 0.02913272],
[13.90374325, -16.54726127, -1.03473650, 0.17110141, 0.02696738],
[15.21915637, -17.90635969, -1.03177327, 0.20521132, 0.02494898],
[16.66056950, -19.39155861, -1.02906215, 0.23978125, 0.02306871],
[18.23998266, -21.01485508, -1.02658147, 0.27479650, 0.02131857],
[19.92330987, -22.74101040, -1.02436829, 0.30931204, 0.01973187],
[21.74253522, -24.60268687, -1.02235923, 0.34384366, 0.01826985],
[23.80711031, -26.71139286, -1.02044861, 0.38006679, 0.01685923],
[25.96939958, -28.91603180, -1.01877074, 0.41512753, 0.01560338],
[28.42241470, -31.41307484, -1.01717386, 0.45187669, 0.01439255],
[30.97172356, -34.00433815, -1.01578005, 0.48716797, 0.01332264],
[33.86030276, -36.93652505, -1.01445222, 0.52412117, 0.01229145],
[36.83234013, -39.94975749, -1.01330152, 0.55927424, 0.01138795],
[40.19224792, -43.35245988, -1.01220379, 0.59603406, 0.01051709],
[43.60470911, -46.80491338, -1.01126057, 0.63060123, 0.00976146],
[47.44867307, -50.69038176, -1.01035917, 0.66668634, 0.00903271],
[51.79880402, -55.08365155, -1.00949910, 0.70441172, 0.00833104],
[56.16397992, -59.48865633, -1.00876888, 0.73943353, 0.00773022],
[61.07929204, -64.44530650, -1.00807058, 0.77596908, 0.00715108],
[66.63981443, -70.04879844, -1.00740384, 0.81414073, 0.00659375],
[72.12543261, -75.57346142, -1.00684607, 0.84899108, 0.00612407],
[78.28863463, -81.77716668, -1.00631198, 0.88530334, 0.00567126],
[85.24430977, -88.77491998, -1.00580134, 0.92319397, 0.00523540],
[93.13264692, -96.70704130, -1.00531393, 0.96279405, 0.00481656],
[100.76607492, -104.37946835, -1.00491445, 0.99820943, 0.00447115],
[109.34320288, -112.99706255, -1.00453175, 1.03509448, 0.00413833],
[119.02452354, -122.72043453, -1.00416565, 1.07356709, 0.00381816],
[130.00600908, -133.74568932, -1.00381603, 1.11376005, 0.00351067],
[140.32285106, -144.10042582, -1.00353717, 1.14867619, 0.00326414],
[151.87901901, -155.69588316, -1.00326958, 1.18498886, 0.00302646],
[164.87942401, -168.73708182, -1.00301316, 1.22280746, 0.00279766],
[179.57299115, -183.47306979, -1.00276785, 1.26225486, 0.00257774],
[196.26449051, -200.20875578, -1.00253356, 1.30346967, 0.00236674],
[211.30618294, -215.28717384, -1.00235402, 1.33781902, 0.00220436],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 1.92e-07
#    Est Max FD Error (based on N/2 run) to R=30.01: 3.2e-7
#    Est Max MOC error for R>30.01: 2.38e-07 
#    Max interpolation error with cubic splines: 5.0787e-07
#    Est Max overall error after cubic interpolation: 1.06e-6
# Run:C8000_G1x15
'C1.15' => [
[-7.22122335, 12.00017689, -1.99998686, -7.22250594, 0.99871660],
[-6.46359956, 10.48495394, -1.99994142, -6.46633754, 0.99725835],
[-5.99587941, 9.54956288, -1.99984212, -6.00025366, 0.99561651],
[-5.52670358, 8.61133716, -1.99958690, -5.53370534, 0.99297508],
[-5.10377348, 7.76575196, -1.99903642, -5.11447952, 0.98924148],
[-4.73335973, 7.02545252, -1.99798417, -4.74889909, 0.98435461],
[-4.42889932, 6.41737587, -1.99630350, -4.45001907, 0.97869402],
[-4.16049492, 5.88187740, -1.99370495, -4.18818985, 0.97200408],
[-3.92064417, 5.40410713, -1.98989650, -3.95594627, 0.96424455],
[-3.70374833, 4.97304640, -1.98455256, -3.74773188, 0.95537483],
[-3.50320086, 4.57573852, -1.97722237, -3.55711708, 0.94522400],
[-3.31575448, 4.20598138, -1.96742691, -3.38098635, 0.93367863],
[-3.13693018, 3.85525515, -1.95447962, -3.21516889, 0.92046342],
[-2.96303279, 3.51678287, -1.93748859, -3.05639164, 0.90521872],
[-2.78879948, 3.18107853, -1.91501260, -2.90019577, 0.88725894],
[-2.60493630, 2.83169158, -1.88418869, -2.73904729, 0.86508530],
[-2.41166721, 2.47137733, -1.84284690, -2.57440876, 0.83797854],
[-2.23685570, 2.15310608, -1.79715238, -2.43031326, 0.81007618],
[-2.07449175, 1.86522436, -1.74794951, -2.30107678, 0.78144774],
[-1.92052593, 1.60003444, -1.69608503, -2.18298859, 0.75219750],
[-1.77215887, 1.35233552, -1.64241706, -2.07357597, 0.72248201],
[-1.62670317, 1.11740900, -1.58752864, -1.97066859, 0.69236194],
[-1.48245812, 0.89240726, -1.53212701, -1.87298413, 0.66202679],
[-1.33792111, 0.67496433, -1.47681715, -1.77949411, 0.63166853],
[-1.19176753, 0.46313953, -1.42212558, -1.68938784, 0.60149153],
[-1.04272862, 0.25521598, -1.36848257, -1.60197723, 0.57169586],
[-0.88964066, 0.04975515, -1.31626517, -1.51671301, 0.54249158],
[-0.73130781, -0.15460629, -1.26576976, -1.43309401, 0.51407311],
[-0.56650331, -0.35915441, -1.21723491, -1.35066523, 0.48662134],
[-0.39403557, -0.56502579, -1.17087465, -1.26904313, 0.46031345],
[-0.21250592, -0.77350433, -1.12682510, -1.18779504, 0.43528461],
[-0.02046856, -0.98582174, -1.08520413, -1.10651832, 0.41165657],
[0.18386317, -1.20347852, -1.04606067, -1.02471707, 0.38950530],
[0.40233318, -1.42791972, -1.00945131, -0.94192704, 0.36890294],
[0.63711660, -1.66082462, -0.97539602, -0.85760559, 0.34989241],
[0.88227318, -1.89611022, -0.94486654, -0.77395164, 0.33303082],
[1.13735568, -2.13357632, -0.91773307, -0.69094244, 0.31824065],
[1.40409840, -2.37506918, -0.89360897, -0.60783270, 0.30529703],
[1.68412532, -2.62221918, -0.87218525, -0.52397088, 0.29401641],
[1.97983860, -2.87723998, -0.85315809, -0.43852336, 0.28421919],
[2.29141047, -3.14037044, -0.83640320, -0.35132345, 0.27581898],
[2.62249805, -3.41476987, -0.82163070, -0.26123540, 0.26864598],
[2.97354324, -3.70086234, -0.80874443, -0.16802755, 0.26262772],
[3.34722222, -4.00091148, -0.79756607, -0.07085967, 0.25765142],
[3.74384022, -4.31527374, -0.78800068, 0.03049607, 0.25364151],
[4.16708495, -4.64700352, -0.77987283, 0.13714503, 0.25048701],
[4.61899735, -4.99783493, -0.77306389, 0.24977098, 0.24810146],
[5.10187922, -5.36971768, -0.76745137, 0.36913271, 0.24639557],
[5.61942140, -5.76567185, -0.76290156, 0.49633646, 0.24527650],
[6.17569982, -6.18899949, -0.75928745, 0.63258240, 0.24465461],
[6.77577794, -6.64374168, -0.75648359, 0.77931142, 0.24444229],
[7.42830852, -7.13663404, -0.75436146, 0.93883955, 0.24455599],
[8.14333252, -7.67542524, -0.75280625, 1.11382209, 0.24491880],
[8.93548224, -8.27129125, -0.75170912, 1.30804320, 0.24546032],
[9.82658248, -8.94078020, -0.75097044, 1.52706574, 0.24611812],
[10.84965197, -9.70881005, -0.75050224, 1.77923577, 0.24683754],
[12.05750483, -10.61511993, -0.75022863, 2.07783651, 0.24757115],
[13.54315231, -11.72957248, -0.75008608, 2.44619113, 0.24827820],
[15.49180491, -13.19115631, -0.75002378, 2.93067258, 0.24892248],
[18.35847761, -15.34119151, -0.75000363, 3.64512261, 0.24946788],
[23.86801548, -19.47335041, -0.75000011, 5.02091096, 0.24986504],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 1.82e-07
#    Est Max FD Error (based on N/2 run) to R=30.01: 1.1e-6
#    Est Max MOC error for R>30.01: 4.15e-07 
#    Max interpolation error with cubic splines: 5.0693e-07
#    Est Max overall error after cubic interpolation: 2.e-6
#    

# Run:C8000_G1x25
'C1.25' => [
[-6.98719148, 12.00017704, -1.99998635, -6.98849863, 0.99869200],
[-6.80500793, 11.63581266, -1.99998035, -6.80657650, 0.99843021],
[-6.07705755, 10.17994506, -1.99991520, -6.08030848, 0.99674392],
[-5.50150534, 9.02893394, -1.99972709, -5.50729292, 0.99419642],
[-5.02516310, 8.07646671, -1.99929261, -5.03449764, 0.99062504],
[-4.63280222, 7.29216601, -1.99845250, -4.64665001, 0.98606655],
[-4.30684259, 6.64095723, -1.99703671, -4.32607155, 0.98061344],
[-4.02488919, 6.07817076, -1.99481105, -4.05044767, 0.97417864],
[-3.77529661, 5.58066060, -1.99149862, -3.80819365, 0.96669694],
[-3.55067448, 5.13381758, -1.98678133, -3.59197975, 0.95810705],
[-3.34475220, 4.72532278, -1.98026165, -3.39565957, 0.94828903],
[-3.15327472, 4.34693882, -1.97147720, -3.21511875, 0.93711582],
[-2.97221092, 3.99097303, -1.95983530, -3.04655933, 0.92438144],
[-2.79775365, 3.65033299, -1.94455009, -2.88653236, 0.90977815],
[-2.62553589, 3.31709690, -1.92446925, -2.73127606, 0.89279405],
[-2.44883912, 2.97932459, -1.89755817, -2.57527685, 0.87241233],
[-2.25006876, 2.60587645, -1.85835690, -2.40446048, 0.84561949],
[-2.07001273, 2.27513678, -1.81395313, -2.25466205, 0.81771825],
[-1.90440150, 1.97862823, -1.76567960, -2.12156939, 0.78912467],
[-1.74853211, 1.70734277, -1.71438818, -2.00081802, 0.75992942],
[-1.59958410, 1.45591264, -1.66109576, -1.88981329, 0.73035134],
[-1.45435235, 1.21861983, -1.60633401, -1.78590838, 0.70038955],
[-1.31071115, 0.99186743, -1.55074086, -1.68747048, 0.67016263],
[-1.16724496, 0.77338939, -1.49503353, -1.59349841, 0.63988873],
[-1.02257308, 0.56111368, -1.43978562, -1.50311128, 0.60976354],
[-0.87541981, 0.35326836, -1.38547961, -1.41558675, 0.57998464],
[-0.72458930, 0.14832757, -1.33252241, -1.33033045, 0.55075216],
[-0.56885478, -0.05515116, -1.28122953, -1.24680243, 0.52225043],
[-0.40691270, -0.25858148, -1.23183742, -1.16449237, 0.49464431],
[-0.23764538, -0.46303037, -1.18459944, -1.08304381, 0.46812343],
[-0.05965985, -0.66980227, -1.13965995, -1.00201507, 0.44282353],
[0.12845838, -0.88011700, -1.09714769, -0.92100848, 0.41887209],
[0.32831762, -1.09531026, -1.05714412, -0.83959044, 0.39636594],
[0.54182576, -1.31692879, -1.01968334, -0.75725762, 0.37536866],
[0.77108187, -1.54659947, -0.98479075, -0.67348765, 0.35593132],
[1.01273544, -1.78065843, -0.95316902, -0.58964341, 0.33846892],
[1.26338026, -2.01595831, -0.92512907, -0.50678157, 0.32315719],
[1.52536748, -2.25497338, -0.90018234, -0.42393100, 0.30972088],
[1.80038155, -2.49940073, -0.87800395, -0.34041920, 0.29797259],
[2.09019115, -2.75091936, -0.85831993, -0.25559347, 0.28775146],
[2.39619631, -3.01082751, -0.84092285, -0.16893619, 0.27893119],
[2.72003154, -3.28059199, -0.82561792, -0.07987400, 0.27139154],
[3.06392513, -3.56213625, -0.81221413, 0.01231629, 0.26501478],
[3.42926434, -3.85667101, -0.80057752, 0.10812664, 0.25971084],
[3.81821825, -4.16603955, -0.79055885, 0.20826129, 0.25538182],
[4.23267753, -4.49185839, -0.78202944, 0.31335646, 0.25193852],
[4.67476367, -4.83593329, -0.77486104, 0.42411577, 0.24929094],
[5.14738593, -5.20068427, -0.76891997, 0.54144566, 0.24734667],
[5.65343777, -5.58851528, -0.76408085, 0.66625140, 0.24601627],
[6.19700048, -6.00273609, -0.76021334, 0.79973298, 0.24520929],
[6.78094044, -6.44572608, -0.75719988, 0.94279190, 0.24483863],
[7.41651752, -6.92620689, -0.75489376, 1.09838246, 0.24481758],
[8.11097435, -7.44981534, -0.75318865, 1.26847288, 0.24506936],
[8.87794651, -8.02698749, -0.75197201, 1.45660017, 0.24552188],
[9.73646129, -8.67217887, -0.75114192, 1.66763444, 0.24611033],
[10.71573945, -9.40746680, -0.75060656, 1.90897530, 0.24677809],
[11.86159663, -10.26734889, -0.75028630, 2.19215835, 0.24747598],
[13.25284499, -11.31104546, -0.75011376, 2.53695697, 0.24816203],
[15.04009513, -12.65160171, -0.75003437, 2.98109130, 0.24879955],
[17.57076004, -14.54964220, -0.75000628, 3.61149099, 0.24935413],
[21.98847076, -17.86293416, -0.75000034, 4.71418529, 0.24978458],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 6.83e-07
#    Est Max FD Error (based on N/2 run) to R=100.05: 5.e-7
#    Est Max MOC error for R>100.05: 3.19e-07 
#    Max interpolation error with cubic splines: 5.0302e-07
#    Est Max overall error after cubic interpolation: 1.32e-6
# Run:P8000_G1x15
'P1.15' => [
[-13.28893903, 12.00017822, -0.99999343, -13.29064951, 0.99914403],
[-12.34002479, 11.05127715, -0.99998303, -12.34277519, 0.99862294],
[-11.68950143, 10.40077365, -0.99995930, -11.69331104, 0.99809165],
[-10.56217756, 9.27353668, -0.99987119, -10.56888078, 0.99663759],
[-9.56983261, 8.28140992, -0.99965356, -9.58086466, 0.99445549],
[-8.79432878, 7.50631194, -0.99924836, -8.81062554, 0.99179147],
[-8.13943110, 6.85210875, -0.99855632, -8.16210517, 0.98855117],
[-7.57045750, 6.28423963, -0.99745916, -7.60068885, 0.98469529],
[-7.06959826, 5.78502783, -0.99582985, -7.10856738, 0.98021971],
[-6.61810630, 5.33590634, -0.99350087, -6.66712643, 0.97505484],
[-6.20504584, 4.92615183, -0.99028254, -6.26554739, 0.96914235],
[-5.82104333, 4.54666519, -0.98594137, -5.89464789, 0.96239130],
[-5.45828751, 4.18999663, -0.98018756, -5.54689135, 0.95467652],
[-5.10865464, 3.84854763, -0.97262348, -5.21461353, 0.94579236],
[-4.76339608, 3.51437990, -0.96267290, -4.88981363, 0.93541223],
[-4.40907644, 3.17554907, -0.94932184, -4.56053801, 0.92289286],
[-4.00695177, 2.79758584, -0.92966428, -4.19268799, 0.90618212],
[-3.64544576, 2.46537564, -0.90754886, -3.86817327, 0.88879270],
[-3.31268525, 2.16728344, -0.88352137, -3.57535410, 0.87084752],
[-3.00179632, 1.89646941, -0.85826288, -3.30742569, 0.85255144],
[-2.71090427, 1.65048619, -0.83271549, -3.06205583, 0.83430397],
[-2.42719860, 1.41791854, -0.80665073, -2.82798165, 0.81570672],
[-2.14229773, 1.19189463, -0.78001322, -2.59831088, 0.79651759],
[-1.84890607, 0.96705133, -0.75279648, -2.36754933, 0.77653249],
[-1.55069718, 0.74658419, -0.72600781, -2.13900186, 0.75631003],
[-1.24412916, 0.52805298, -0.69994121, -1.91027949, 0.73593114],
[-0.92557889, 0.30913986, -0.67485852, -1.67912787, 0.71548678],
[-0.59097574, 0.08739799, -0.65097727, -1.44317308, 0.69506673],
[-0.23575643, -0.13976032, -0.62848810, -1.19992093, 0.67477164],
[0.14591461, -0.37553383, -0.60753014, -0.94627089, 0.65468742],
[0.56062347, -0.62336483, -0.58825049, -0.67893571, 0.63493846],
[0.97930806, -0.86616193, -0.57207141, -0.41689469, 0.61714810],
[1.40191683, -1.10498978, -0.55862898, -0.15950975, 0.60126321],
[1.83550254, -1.34468020, -0.54738242, 0.09801942, 0.58696284],
[2.28503433, -1.58855422, -0.53797683, 0.35891011, 0.57406926],
[2.75491779, -1.83942840, -0.53014252, 0.62585879, 0.56246359],
[3.24944254, -2.09992810, -0.52366066, 0.90136728, 0.55206018],
[3.77332665, -2.37281152, -0.51834312, 1.18808157, 0.54278878],
[4.33191962, -2.66109186, -0.51402523, 1.48891542, 0.53459056],
[4.93128026, -2.96808834, -0.51056144, 1.80710164, 0.52741639],
[5.58017144, -3.29844816, -0.50781469, 2.14724412, 0.52120770],
[6.28723580, -3.65670423, -0.50567049, 2.51382418, 0.51593001],
[7.06430941, -4.04896709, -0.50402066, 2.91294950, 0.51153849],
[7.93503951, -4.48724749, -0.50275875, 3.35671316, 0.50795615],
[8.92315799, -4.98353006, -0.50180877, 3.85715403, 0.50514849],
[10.07059776, -5.55888934, -0.50110426, 4.43548199, 0.50305342],
[11.43844570, -6.24394829, -0.50060258, 5.12249134, 0.50160472],
[12.92246978, -6.98660865, -0.50030336, 5.86620578, 0.50078381],
[14.70386796, -7.87767191, -0.50012910, 6.75783613, 0.50032655],
[16.94461877, -8.99822272, -0.50004272, 7.87865330, 0.50010741],
[19.96185964, -10.50691042, -0.50000963, 9.38744123, 0.50002377],
[24.57791348, -12.81495464, -0.50000095, 11.69551093, 0.50000236],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 6.31e-07
#    Est Max FD Error (based on N/2 run) to R=100.05: 4.54e-7
#    Est Max MOC error for R>100.05: 3.19e-07 
#    Max interpolation error with cubic splines: 5.0528e-07
#    Est Max overall error after cubic interpolation: 1.3e-6
# Run:P8000_G1x25
'P1.25' => [
[-12.81135743, 12.00017781, -0.99999317, -12.81310068, 0.99912762],
[-11.73808762, 10.92692396, -0.99998004, -11.74107082, 0.99850621],
[-10.88322119, 10.07209015, -0.99994461, -10.88779895, 0.99770601],
[-9.84932023, 9.03829229, -0.99983992, -9.85700821, 0.99614187],
[-8.91446629, 8.10368575, -0.99959289, -8.92676199, 0.99381696],
[-8.16083404, 7.35051107, -0.99913624, -8.17880298, 0.99094290],
[-7.52468343, 6.71512759, -0.99837213, -7.54945481, 0.98748201],
[-6.97101994, 6.16266576, -0.99717901, -7.00380087, 0.98338952],
[-6.48109649, 5.67452069, -0.99542151, -6.52312837, 0.97864418],
[-6.03866561, 5.23462746, -0.99293066, -6.09130741, 0.97318453],
[-5.63279153, 4.83227276, -0.98950880, -5.69753929, 0.96694293],
[-5.25472560, 4.45899094, -0.98491699, -5.33327278, 0.95982693],
[-4.89634409, 4.10704514, -0.97884592, -4.99070242, 0.95169174],
[-4.55014818, 3.76948334, -0.97088579, -4.66280563, 0.94232815],
[-4.20707679, 3.43811395, -0.96042101, -4.34134760, 0.93137856],
[-3.85213224, 3.09961511, -0.94630863, -4.01305557, 0.91809687],
[-3.45719395, 2.72971503, -0.92605411, -3.65379064, 0.90079682],
[-3.10287431, 2.40546865, -0.90349840, -3.33771703, 0.88295596],
[-2.77573792, 2.11380496, -0.87908825, -3.05182339, 0.86461435],
[-2.47141047, 1.85008032, -0.85369020, -2.79147994, 0.84611775],
[-2.18479109, 1.60904896, -0.82796318, -2.55158973, 0.82766637],
[-1.90420993, 1.38040245, -0.80172639, -2.32198171, 0.80890202],
[-1.62112916, 1.15724692, -0.77489059, -2.09572832, 0.78956280],
[-1.32911589, 0.93497962, -0.74752907, -1.86809507, 0.76950713],
[-1.03181230, 0.71676283, -0.72065730, -1.64233038, 0.74930344],
[-0.72569561, 0.50019681, -0.69457156, -1.41607542, 0.72903745],
[-0.40677948, 0.28274464, -0.66950669, -1.18683016, 0.70878031],
[-0.07097992, 0.06199941, -0.64568758, -0.95224216, 0.68862669],
[0.28611046, -0.16448601, -0.62332261, -0.70994728, 0.66869082],
[0.67060629, -0.40004957, -0.60254056, -0.45667612, 0.64904690],
[1.08979909, -0.64850904, -0.58346750, -0.18871406, 0.62979115],
[1.50502572, -0.88743018, -0.56783330, 0.06920484, 0.61285747],
[1.92580093, -1.12353230, -0.55482890, 0.32382771, 0.59772449],
[2.35851376, -1.36117928, -0.54396181, 0.57945642, 0.58410363],
[2.80803448, -1.60358698, -0.53489516, 0.83919550, 0.57182191],
[3.27850050, -1.85339727, -0.52737460, 1.10555037, 0.56076715],
[3.77426400, -2.11325122, -0.52118686, 1.38103191, 0.55085063],
[4.29973714, -2.38573773, -0.51615163, 1.66809608, 0.54200609],
[4.86073072, -2.67410429, -0.51210325, 1.96988803, 0.53416595],
[5.46298327, -2.98150460, -0.50890130, 2.28944642, 0.52728556],
[6.11591065, -3.31291948, -0.50640811, 2.63169662, 0.52130077],
[6.82570334, -3.67165099, -0.50451637, 2.99982340, 0.51619563],
[7.61044829, -4.06697547, -0.50310368, 3.40313422, 0.51188753],
[8.48589587, -4.50693376, -0.50207774, 3.84963114, 0.50834760],
[9.47178144, -5.00153857, -0.50135204, 4.34933312, 0.50554081],
[10.62908783, -5.58143212, -0.50083601, 4.93304476, 0.50336705],
[12.04507471, -6.29032501, -0.50046912, 5.64457578, 0.50178458],
[13.80713903, -7.17194168, -0.50022214, 6.52775543, 0.50078225],
[15.90569283, -8.22152021, -0.50008523, 7.57806504, 0.50028289],
[18.38258847, -9.46009166, -0.50002563, 8.81691744, 0.50008315],
[21.91502323, -11.22635210, -0.50000448, 10.58327305, 0.50001429],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 5.96e-07
#    Est Max FD Error (based on N/2 run) to R=100.05: 4.26e-7
#    Est Max MOC error for R>100.05: 4.83e-07 
#    Max interpolation error with cubic splines: 5.0523e-07
#    Est Max overall error after cubic interpolation: 1.41e-6
#    

# Run:P8000_G1x35
'P1.35' => [
[-12.49881322, 12.00017800, -0.99999294, -12.50058592, 0.99911287],
[-10.64535281, 10.14676141, -0.99995081, -10.64983703, 0.99775298],
[-9.77048950, 9.27197137, -0.99987443, -9.77744261, 0.99651182],
[-8.79439200, 8.29608166, -0.99966677, -8.80574315, 0.99429422],
[-8.00639012, 7.50847910, -0.99926812, -8.02326544, 0.99149770],
[-7.34555587, 6.84832948, -0.99858578, -7.36910795, 0.98810302],
[-6.77502002, 6.27887962, -0.99750640, -6.80645108, 0.98407940],
[-6.27197361, 5.77745931, -0.99589695, -6.31254043, 0.97939466],
[-5.81985675, 5.32767942, -0.99359798, -5.87091174, 0.97399774],
[-5.40679629, 4.91787674, -0.99042059, -5.46982129, 0.96782382],
[-5.02341869, 4.53894438, -0.98613573, -5.10008452, 0.96078315],
[-4.66191267, 4.18342091, -0.98045939, -4.75416295, 0.95275191],
[-4.31446696, 3.84399445, -0.97300863, -4.42468903, 0.94353336],
[-3.97327027, 3.51359489, -0.96325214, -4.10453440, 0.93283489],
[-3.62582456, 3.18107691, -0.95025764, -3.78258648, 0.92006083],
[-3.24182205, 2.81959400, -0.93166638, -3.43236937, 0.90354250],
[-2.87937868, 2.48577591, -0.90963506, -3.10807622, 0.88555277],
[-2.54724308, 2.18754581, -0.88560625, -2.81697003, 0.86707818],
[-2.23572063, 1.91556944, -0.86006374, -2.54976122, 0.84819337],
[-1.94732819, 1.67120037, -0.83435735, -2.30780708, 0.82960235],
[-1.66768438, 1.44151775, -0.80817236, -2.07842520, 0.81082809],
[-1.38791012, 1.21914940, -0.78141950, -1.85425870, 0.79160535],
[-1.09992665, 0.99807107, -0.75401458, -1.62915875, 0.77168246],
[-0.80599556, 0.78046113, -0.72686059, -1.40530669, 0.75153558],
[-0.50395179, 0.56495638, -0.70040065, -1.18137594, 0.73134496],
[-0.19033627, 0.34935230, -0.67492427, -0.95519706, 0.71121244],
[0.13904068, 0.13111779, -0.65064890, -0.72426815, 0.69121572],
[0.48857883, -0.09222592, -0.62778416, -0.48616082, 0.67145625],
[0.86374543, -0.32365308, -0.60649153, -0.23795623, 0.65202372],
[1.27176294, -0.56699278, -0.58688539, 0.02412407, 0.63298899],
[1.68442611, -0.80567171, -0.57041314, 0.28173367, 0.61588226],
[2.10068813, -1.04017503, -0.55675287, 0.53486682, 0.60066478],
[2.52754985, -1.27531215, -0.54535006, 0.78828431, 0.58699705],
[2.96966828, -1.51424252, -0.53584731, 1.04502098, 0.57469156],
[3.43150506, -1.75982432, -0.52796544, 1.30781058, 0.56360995],
[3.91733458, -2.01468370, -0.52148108, 1.57914465, 0.55365272],
[4.43159154, -2.28144088, -0.51620561, 1.86150722, 0.54474322],
[4.97994998, -2.56328781, -0.51196843, 2.15797852, 0.53680838],
[5.56618895, -2.86239443, -0.50863432, 2.47055930, 0.52981957],
[6.20024001, -3.18402823, -0.50605205, 2.80447853, 0.52369245],
[6.89031635, -3.53252799, -0.50410562, 3.16396436, 0.51839283],
[7.64701952, -3.91341050, -0.50268545, 3.55445149, 0.51388241],
[8.49097028, -4.33719515, -0.50168361, 3.98646376, 0.51009315],
[9.43807454, -4.81199628, -0.50101581, 4.46803107, 0.50700854],
[10.52247598, -5.35504601, -0.50059481, 5.01642070, 0.50456811],
[11.81509235, -6.00193967, -0.50034152, 5.66732551, 0.50269593],
[13.41475628, -6.80218329, -0.50019112, 6.47029163, 0.50136174],
[15.88516542, -8.03770428, -0.50007896, 7.70752952, 0.50044235],
[19.14015409, -9.66533952, -0.50001932, 9.33575338, 0.50009168],
[23.08999916, -11.64029631, -0.50000284, 11.31083487, 0.50001290],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 5.69e-07
#    Est Max FD Error (based on N/2 run) to R=100.04: 4.05e-7
#    Est Max MOC error for R>100.04: 5.7e-07 
#    Max interpolation error with cubic splines: 5.0385e-07
#    Est Max overall error after cubic interpolation: 1.47e-6
#    

# Run:P8000_G1x45
'P1.45' => [
[-12.26556632, 12.00017797, -0.99999273, -12.26736564, 0.99909954],
[-11.58754725, 11.32216656, -0.99998568, -11.59007361, 0.99873525],
[-10.22088578, 9.95555211, -0.99993786, -10.22589523, 0.99748917],
[-9.16417589, 8.89896091, -0.99981819, -9.17268690, 0.99572721],
[-8.26369158, 7.99874218, -0.99955284, -8.27707323, 0.99326762],
[-7.52599434, 7.26153250, -0.99906630, -7.54539795, 0.99021383],
[-6.89984217, 6.63619200, -0.99825797, -6.92646234, 0.98653783],
[-6.35523976, 6.09284825, -0.99700844, -6.39031380, 0.98221245],
[-5.87187778, 5.61134032, -0.99517650, -5.91671037, 0.97719941],
[-5.43475858, 5.17685496, -0.99259170, -5.49076897, 0.97143891],
[-5.03325881, 4.77899517, -0.98905317, -5.10201139, 0.96485979],
[-4.65894230, 4.40961147, -0.98431913, -4.74220584, 0.95736727],
[-4.30368530, 4.06097657, -0.97807223, -4.40356961, 0.94880723],
[-3.96030145, 3.72645791, -0.96989789, -4.07940675, 0.93896730],
[-3.61972967, 3.39788340, -0.95916563, -3.76152614, 0.92747494],
[-3.26697228, 3.06197327, -0.94470155, -3.43674512, 0.91355317],
[-2.87765811, 2.69802125, -0.92416575, -3.08448496, 0.89563698],
[-2.52771279, 2.37848973, -0.90131825, -2.77422101, 0.87721744],
[-2.20416307, 2.09077811, -0.87660578, -2.49340499, 0.85834934],
[-1.90389748, 1.83135090, -0.85099791, -2.23847725, 0.83946763],
[-1.62081671, 1.59408625, -0.82507167, -2.00347587, 0.82071583],
[-1.34348865, 1.36891920, -0.79864831, -1.77848951, 0.80173779],
[-1.06336987, 1.14898750, -0.77162887, -1.55663105, 0.78226986],
[-0.77381309, 0.92956481, -0.74406601, -1.33303488, 0.76216337],
[-0.47864325, 0.71396483, -0.71701137, -1.11105332, 0.74201447],
[-0.17440126, 0.49985913, -0.69077390, -0.88837924, 0.72191727],
[0.14279009, 0.28480633, -0.66560239, -0.66258923, 0.70195207],
[0.47733058, 0.06621045, -0.64170276, -0.43110206, 0.68219450],
[0.83375011, -0.15841754, -0.61928510, -0.19146991, 0.66275217],
[1.21802133, -0.39228828, -0.59849612, 0.05948360, 0.64370580],
[1.63193430, -0.63600122, -0.57969459, 0.32205152, 0.62537437],
[2.04308237, -0.87106703, -0.56426648, 0.57579592, 0.60928575],
[2.46068529, -1.10393537, -0.55143311, 0.82717198, 0.59493184],
[2.89052159, -1.33857826, -0.54073043, 1.08006016, 0.58204058],
[3.33724132, -1.57806906, -0.53182958, 1.33741347, 0.57043617],
[3.80519666, -1.82515039, -0.52447459, 1.60184493, 0.55999239],
[4.29848733, -2.08232008, -0.51845868, 1.87570827, 0.55061930],
[4.82182223, -2.35231624, -0.51359999, 2.16160928, 0.54223953],
[5.37987095, -2.63779861, -0.50974242, 2.46206328, 0.53479715],
[5.97986284, -2.94268716, -0.50673506, 2.78089802, 0.52822345],
[6.62740888, -3.27003764, -0.50445513, 3.12102146, 0.52248783],
[7.33680482, -3.62725675, -0.50276987, 3.48983379, 0.51750714],
[8.11597236, -4.01850009, -0.50158385, 3.89133512, 0.51327340],
[8.99254729, -4.45779447, -0.50078822, 4.33961085, 0.50969768],
[9.98078716, -4.95242809, -0.50030690, 4.84179225, 0.50678513],
[11.11621924, -5.52032907, -0.50005737, 5.41581584, 0.50447904],
[12.48556131, -6.20499981, -0.49996672, 6.10529586, 0.50268943],
[14.17421742, -7.04926771, -0.49997209, 6.95297168, 0.50140472],
[16.81808720, -8.37118336, -0.50000855, 8.27718880, 0.50047307],
[19.88879413, -9.90657244, -0.50000983, 9.81333342, 0.50011700],
[26.98058964, -13.45249433, -0.50000045, 13.35946456, 0.50000357],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 1.92e-07
#    Est Max FD Error (based on N/2 run) to R=30.01: 1.53e-7	
#    Est Max MOC error for R>30.01: 2.18e-07 
#    Max interpolation error with cubic splines: 5.0729e-07
#    Est Max overall error after cubic interpolation: 8.7e-7
#    

# Run:C8000_G1x35
'C1.35' => [
[-6.83587178, 12.00017707, -1.99998588, -6.83720102, 0.99866989],
[-6.45426051, 11.23696227, -1.99996972, -6.45620796, 0.99805069],
[-5.79315930, 9.91480132, -1.99988784, -5.79693472, 0.99621765],
[-5.21223228, 8.75307142, -1.99963926, -5.21899123, 0.99321938],
[-4.76500358, 7.85887510, -1.99911705, -4.77559324, 0.98935878],
[-4.39470180, 7.11875555, -1.99815120, -4.41007032, 0.98452713],
[-4.08184326, 6.49384131, -1.99655226, -4.10290835, 0.97874832],
[-3.80851195, 5.94843024, -1.99406797, -3.83627357, 0.97193336],
[-3.56608577, 5.46542304, -1.99042334, -3.60156658, 0.96405637],
[-3.34683845, 5.02955276, -1.98527619, -3.39115258, 0.95502603],
[-3.14561913, 4.63074074, -1.97823285, -3.19998150, 0.94474868],
[-2.95750090, 4.25943714, -1.96878004, -3.02332324, 0.93304366],
[-2.77901253, 3.90908763, -1.95630502, -2.85794082, 0.91970997],
[-2.60623489, 3.57242595, -1.93995491, -2.70032164, 0.90440103],
[-2.43457700, 3.24118012, -1.91844955, -2.54656745, 0.88653581],
[-2.25645135, 2.90193122, -1.88941179, -2.39053047, 0.86490637],
[-2.06105473, 2.53658357, -1.84847186, -2.22415901, 0.83731970],
[-1.88536574, 2.21570370, -1.80296439, -2.07948851, 0.80903361],
[-1.72300179, 1.92687297, -1.75376156, -1.95044279, 0.78013558],
[-1.56957478, 1.66172942, -1.70170873, -1.83298204, 0.75071833],
[-1.42299932, 1.41619553, -1.64802463, -1.72509962, 0.72111621],
[-1.27902219, 1.18287196, -1.59279308, -1.62342877, 0.69108827],
[-1.13622786, 0.95941674, -1.53687322, -1.52690005, 0.66087957],
[-0.99323561, 0.74366033, -1.48097474, -1.43456055, 0.63070525],
[-0.84876974, 0.53372328, -1.42569497, -1.34561657, 0.60077370],
[-0.70141046, 0.32766187, -1.37144893, -1.25927744, 0.57124440],
[-0.55004297, 0.12410480, -1.31866572, -1.17501819, 0.54232684],
[-0.39336372, -0.07845479, -1.26762947, -1.09227663, 0.51418472],
[-0.23025503, -0.28116208, -1.21862950, -1.01065516, 0.48700866],
[-0.05932394, -0.48539839, -1.17183077, -0.92967379, 0.46093547],
[0.12073208, -0.69232096, -1.12740383, -0.84895441, 0.43611220],
[0.31140936, -0.90321204, -1.08545613, -0.76807916, 0.41265174],
[0.51442029, -1.11948470, -1.04605047, -0.68659020, 0.39063860],
[0.73166661, -1.34264170, -1.00922774, -0.60400608, 0.37013911],
[0.96535873, -1.57439000, -0.97499994, -0.51977936, 0.35119640],
[1.20801198, -1.80718324, -0.94452395, -0.43664563, 0.33447120],
[1.46064584, -2.04228907, -0.91743776, -0.35405770, 0.31976677],
[1.72510822, -2.28164318, -0.89334420, -0.27124907, 0.30686186],
[2.00312943, -2.52695058, -0.87193160, -0.18755267, 0.29557932],
[2.29622539, -2.77964755, -0.85295618, -0.10240358, 0.28577684],
[2.60648560, -3.04160378, -0.83617998, -0.01509700, 0.27731490],
[2.93502379, -3.31382413, -0.82145043, 0.07478180, 0.27009706],
[3.28406824, -3.59822601, -0.80858149, 0.16795292, 0.26400952],
[3.65505662, -3.89606116, -0.79743935, 0.26492101, 0.25896342],
[4.05023314, -4.20923093, -0.78787578, 0.36640869, 0.25486229],
[4.47134493, -4.53924023, -0.77976625, 0.47301557, 0.25161941],
[4.92111222, -4.88836212, -0.77297366, 0.58559543, 0.24914204],
[5.40223876, -5.25885205, -0.76737136, 0.70500083, 0.24734133],
[5.91781495, -5.65326392, -0.76283329, 0.83218443, 0.24612826],
[6.47231917, -6.07520617, -0.75922861, 0.96844188, 0.24541319],
[7.07141658, -6.52917158, -0.75643054, 1.11535775, 0.24510980],
[7.72216023, -7.02068274, -0.75431773, 1.27485426, 0.24513572],
[8.43579322, -7.55839739, -0.75277034, 1.44987933, 0.24541426],
[9.22654876, -8.15318935, -0.75168077, 1.64411783, 0.24587552],
[10.11645260, -8.82175740, -0.75094907, 1.86318101, 0.24645750],
[11.13832414, -9.58886949, -0.75048724, 2.11536494, 0.24710582],
[12.34597812, -10.49501542, -0.75021896, 2.41419926, 0.24777400],
[13.83382617, -11.61110769, -0.75008063, 2.78335504, 0.24842206],
[15.79207920, -13.07988488, -0.75002133, 3.27045135, 0.24901503],
[18.69995298, -15.26081715, -0.75000292, 3.99537286, 0.24951885],
[24.50749830, -19.61648034, -0.75000006, 5.44578171, 0.24988677],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 2.01e-07
#    Est Max FD Error (based on N/2 run) to R=30.01: 1.44e-7
#    Est Max MOC error for R>30.01: 5.19e-08 
#    Max interpolation error with cubic splines: 5.0596e-07
#    Est Max overall error after cubic interpolation: 7.e-7
#    

# Run:C8000_G1x45
'C1.45' => [
[-6.72410370, 12.00017703, -1.99998546, -6.72545290, 0.99864991],
[-5.80725127, 10.16650910, -1.99991096, -5.81062951, 0.99661620],
[-5.17447977, 8.90107868, -1.99968685, -5.18084953, 0.99361092],
[-4.71213442, 7.97662647, -1.99920982, -4.72226608, 0.98982092],
[-4.33175488, 7.21631743, -1.99831154, -4.34660743, 0.98504942],
[-4.01184191, 6.57724717, -1.99680597, -4.03234380, 0.97932008],
[-3.73414356, 6.02303430, -1.99445460, -3.76128113, 0.97256839],
[-3.48847011, 5.53344191, -1.99098554, -3.52326636, 0.96475377],
[-3.26653539, 5.09208090, -1.98605866, -3.31011230, 0.95577650],
[-3.06318853, 4.68886645, -1.97928808, -3.11676202, 0.94554891],
[-2.87383858, 4.31489930, -1.97019369, -2.93878798, 0.93392301],
[-2.69467835, 3.96293809, -1.95817848, -2.77261590, 0.92069685],
[-2.52190072, 3.62590190, -1.94243831, -2.61481216, 0.90555314],
[-2.35113863, 3.29589021, -1.92178767, -2.46164054, 0.88796450],
[-2.17571726, 2.96108934, -1.89413483, -2.30767465, 0.86689028],
[-1.97763324, 2.58972102, -1.85375807, -2.13863276, 0.83917281],
[-1.80037659, 2.26499815, -1.80868746, -1.99235011, 0.81079004],
[-1.63711684, 1.97361221, -1.75977324, -1.86231067, 0.78181947],
[-1.48324193, 1.70675569, -1.70785755, -1.74425124, 0.75235063],
[-1.33668937, 1.46035248, -1.65422563, -1.63614612, 0.72274744],
[-1.19355911, 1.22750004, -1.59915550, -1.53482943, 0.69286085],
[-1.05112673, 1.00371743, -1.54303410, -1.43829155, 0.66266740],
[-0.90878839, 0.78808575, -1.48690376, -1.34611574, 0.63254694],
[-0.76501231, 0.57832133, -1.43129426, -1.25732860, 0.60265454],
[-0.61857173, 0.37274784, -1.37672172, -1.17124745, 0.57318805],
[-0.46816309, 0.16971476, -1.32356678, -1.08722591, 0.54431597],
[-0.31266871, -0.03204900, -1.27219326, -1.00479585, 0.51623312],
[-0.15069563, -0.23405318, -1.22281267, -0.92340891, 0.48907748],
[0.01892380, -0.43740219, -1.17566799, -0.84269623, 0.46302173],
[0.19752415, -0.64330815, -1.13091986, -0.76225607, 0.43820289],
[0.38676213, -0.85324069, -1.08863961, -0.68159829, 0.41471150],
[0.58818348, -1.06842987, -1.04892940, -0.60033674, 0.39265351],
[0.80366665, -1.29036556, -1.01182904, -0.51799471, 0.37209547],
[1.03557838, -1.52091851, -0.97732293, -0.43396575, 0.35306702],
[1.27636939, -1.75245390, -0.94659871, -0.35103092, 0.33624702],
[1.52693514, -1.98612805, -0.91930226, -0.26868591, 0.32144773],
[1.78938197, -2.22411974, -0.89500381, -0.18608151, 0.30843430],
[2.06516368, -2.46788388, -0.87341431, -0.10264041, 0.29704504],
[2.35610922, -2.71913211, -0.85426397, -0.01770662, 0.28712666],
[2.66381690, -2.97931286, -0.83734303, 0.06928024, 0.27855644],
[2.99000361, -3.24993931, -0.82246484, 0.15890154, 0.27122260],
[3.33624868, -3.53238944, -0.80947133, 0.25169617, 0.26502682],
[3.70450844, -3.82833862, -0.79820703, 0.34830521, 0.25987120],
[4.09677739, -4.13948291, -0.78853199, 0.44938120, 0.25566453],
[4.51496712, -4.46745424, -0.78031786, 0.55556207, 0.25231961],
[4.96159790, -4.81436626, -0.77343130, 0.66764766, 0.24974622],
[5.43917476, -5.18232113, -0.76774680, 0.78643891, 0.24785653],
[5.95119089, -5.57418183, -0.76313318, 0.91298630, 0.24656049],
[6.50172550, -5.99324995, -0.75946315, 1.04848492, 0.24577004],
[7.09704600, -6.44447484, -0.75660648, 1.19466610, 0.24539866],
[7.74320605, -6.93262171, -0.75444638, 1.35320582, 0.24536480],
[8.45085032, -7.46590123, -0.75286193, 1.52690491, 0.24559088],
[9.23441347, -8.05534279, -0.75174290, 1.71949647, 0.24600664],
[10.11532220, -8.71719719, -0.75098872, 1.93644287, 0.24654975],
[11.12519662, -9.47533437, -0.75051068, 2.18574167, 0.24716530],
[12.31585202, -10.36874644, -0.75023142, 2.48042318, 0.24780642],
[13.77750039, -11.46519800, -0.75008629, 2.84310845, 0.24843287],
[15.69075275, -12.90022944, -0.75002329, 3.31901435, 0.24900993],
[18.50302430, -15.00946197, -0.75000332, 4.02007078, 0.24950445],
[23.96996114, -19.10966925, -0.75000007, 5.38532407, 0.24987310],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 2.08e-07
#    Est Max FD Error (based on N/2 run) to R=10.00: 1.41e-7
#    Est Max MOC error for R>10.00: 3.31e-07 
#    Max interpolation error with cubic splines: 5.1036e-07
#    Est Max overall error after cubic interpolation: 9.7e-7
#    

# Run:S8000_G1x5
'S1.5' => [
[-4.55347190, 12.00017655, -2.99997789, -4.55455853, 0.99836918],
[-4.18596949, 10.89768266, -2.99993340, -4.18785600, 0.99716762],
[-4.05350291, 10.50029242, -2.99990090, -4.05580457, 0.99654363],
[-3.53125086, 8.93364915, -2.99957337, -3.53629549, 0.99241497],
[-3.21088999, 7.97279652, -2.99888509, -3.21905874, 0.98770108],
[-2.96611234, 7.23886697, -2.99767584, -2.97792424, 0.98219035],
[-2.74944045, 6.58956000, -2.99555672, -2.76582033, 0.97526309],
[-2.56143419, 6.02666149, -2.99221766, -2.58319712, 0.96708012],
[-2.39526308, 5.52982020, -2.98725619, -2.42325063, 0.95759986],
[-2.24502988, 5.08153010, -2.98015298, -2.28017578, 0.94668640],
[-2.10830823, 4.67470204, -2.97040393, -2.15155850, 0.93433451],
[-1.98099250, 4.29731059, -2.95725740, -2.03346780, 0.92030476],
[-1.86028330, 3.94134031, -2.93978238, -1.92331461, 0.90432509],
[-1.74461299, 3.60254965, -2.91695941, -1.81973231, 0.88617123],
[-1.63028637, 3.27070082, -2.88694214, -1.71959114, 0.86513074],
[-1.51325128, 2.93507512, -2.84680022, -1.61977066, 0.84007436],
[-1.38215129, 2.56551342, -2.78858641, -1.51171896, 0.80750618],
[-1.26222591, 2.23495725, -2.72194123, -1.41687750, 0.77351135],
[-1.15201635, 1.93887296, -2.64943770, -1.33350618, 0.73894970],
[-1.04832737, 1.66808369, -2.57236092, -1.25868183, 0.70393585],
[-0.94990122, 1.41876777, -2.49280686, -1.19110532, 0.66897382],
[-0.85392059, 1.18340316, -2.41109682, -1.12857723, 0.63383587],
[-0.75829382, 0.95682102, -2.32760609, -1.06965865, 0.59840576],
[-0.66268995, 0.73829660, -2.24397202, -1.01413721, 0.56316086],
[-0.56627936, 0.52596346, -2.16118195, -0.96152838, 0.52835443],
[-0.46805113, 0.31769822, -2.07986337, -0.91132146, 0.49414166],
[-0.36716153, 0.11189864, -2.00063071, -0.86316978, 0.46071729],
[-0.26283978, -0.09276569, -1.92402766, -0.81681855, 0.42828160],
[-0.15427885, -0.29758936, -1.85046556, -0.77204493, 0.39700452],
[-0.04051211, -0.50404738, -1.78018140, -0.72861269, 0.36699799],
[0.07930129, -0.71326662, -1.71344880, -0.68638481, 0.33839803],
[0.20616912, -0.92656977, -1.65040726, -0.64520602, 0.31128595],
[0.34108278, -1.14515150, -1.59118300, -0.60496892, 0.28573600],
[0.48530676, -1.37055067, -1.53577202, -0.56552662, 0.26176244],
[0.64025190, -1.60441298, -1.48413646, -0.52674498, 0.23935997],
[0.80120318, -1.83947178, -1.43788474, -0.48987875, 0.21923567],
[0.96905391, -2.07725195, -1.39639850, -0.45463714, 0.20112856],
[1.14519394, -2.31983916, -1.35904306, -0.42068751, 0.18476609],
[1.33067731, -2.56871456, -1.32536841, -0.38782498, 0.16995361],
[1.52704710, -2.82590819, -1.29492098, -0.35580666, 0.15649289],
[1.73529796, -3.09263492, -1.26740875, -0.32452387, 0.14425603],
[1.95688371, -3.37064224, -1.24252193, -0.29382609, 0.13310673],
[2.19298433, -3.66127373, -1.22002768, -0.26363059, 0.12294301],
[2.44560688, -3.96683505, -1.19965584, -0.23377749, 0.11364588],
[2.71582301, -4.28844446, -1.18125143, -0.20424667, 0.10514891],
[3.00563062, -4.62830005, -1.16461256, -0.17493132, 0.09736473],
[3.31710309, -4.98863398, -1.14956798, -0.14574682, 0.09022027],
[3.65210404, -5.37139776, -1.13597973, -0.11665093, 0.08365859],
[4.01329045, -5.77941463, -1.12369729, -0.08755433, 0.07761707],
[4.40731377, -6.21990350, -1.11249204, -0.05810889, 0.07199259],
[4.81820019, -6.67494267, -1.10270893, -0.02958500, 0.06697642],
[5.27420201, -7.17564575, -1.09362189, -0.00015873, 0.06221049],
[5.76873701, -7.71438522, -1.08539836, 0.02948574, 0.05779335],
[6.30622590, -8.29571093, -1.07794453, 0.05941745, 0.05369007],
[6.89108409, -8.92411835, -1.07118288, 0.08967334, 0.04987380],
[7.52872273, -9.60512744, -1.06503812, 0.12031142, 0.04631789],
[8.22474951, -10.34441776, -1.05944716, 0.15136545, 0.04300140],
[8.98556948, -11.14846947, -1.05435231, 0.18287292, 0.03990512],
[9.81798577, -12.02413643, -1.04970429, 0.21485519, 0.03701329],
[10.72940023, -12.97885937, -1.04545972, 0.24732605, 0.03431216],
[11.72801384, -14.02087400, -1.04157991, 0.28029735, 0.03178926],
[12.82222710, -15.15858526, -1.03803251, 0.31375879, 0.02943456],
[14.02124025, -16.40119633, -1.03478816, 0.34769920, 0.02723849],
[15.33525337, -17.75891282, -1.03182001, 0.38210984, 0.02519177],
[16.77486651, -19.24232138, -1.02910489, 0.41696736, 0.02328628],
[18.35187966, -20.86321609, -1.02662113, 0.45225435, 0.02151383],
[20.07450421, -22.62968827, -1.02435459, 0.48786109, 0.01987060],
[21.94821060, -24.54703942, -1.02229317, 0.52363275, 0.01835359],
[23.93152640, -26.57270175, -1.02046227, 0.55864290, 0.01698734],
[26.18206273, -28.86726990, -1.01871980, 0.59534983, 0.01566971],
[28.55384534, -31.28157286, -1.01717970, 0.63107890, 0.01449031],
[31.24648644, -34.01843975, -1.01571349, 0.66852348, 0.01335392],
[34.06786600, -36.88228282, -1.01442468, 0.70473637, 0.01234355],
[37.27036438, -40.12895189, -1.01319703, 0.74266510, 0.01137063],
[40.60072005, -43.50142206, -1.01212467, 0.77906540, 0.01051202],
[44.37705921, -47.32156096, -1.01110237, 0.81715822, 0.00968547],
[48.26567624, -51.25157936, -1.01021590, 0.85337778, 0.00896211],
[52.66573626, -55.69467475, -1.00936980, 0.89123735, 0.00826567],
[57.13896187, -60.20814289, -1.00864245, 0.92682677, 0.00766204],
[62.18308763, -65.29405935, -1.00794708, 0.96396919, 0.00708051],
[67.89813408, -71.05257191, -1.00728336, 1.00279183, 0.00652119],
[73.64007289, -76.83466958, -1.00671970, 1.03882718, 0.00604277],
[80.11411922, -83.35042563, -1.00618056, 1.07641714, 0.00558209],
[87.44864276, -90.72833769, -1.00566570, 1.11568964, 0.00513922],
[94.69527613, -98.01442637, -1.00523495, 1.15154648, 0.00476639],
[102.84970671, -106.20981709, -1.00482246, 1.18891096, 0.00440731],
[112.06845379, -115.47115161, -1.00442808, 1.22790555, 0.00406203],
[122.54313624, -125.99018698, -1.00405167, 1.26866839, 0.00373059],
[132.69891545, -136.18550635, -1.00374321, 1.30513352, 0.00345756],
[144.13263510, -147.66028932, -1.00344774, 1.34312382, 0.00319475],
[157.06605372, -160.63641955, -1.00316516, 1.38276439, 0.00294220],
[171.77102629, -175.38589329, -1.00289539, 1.42419650, 0.00269992],
[185.61860804, -189.27204142, -1.00268028, 1.46020680, 0.00250590],
[201.16611766, -204.85957576, -1.00247393, 1.49767680, 0.00231904],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 2.12e-07
#    Est Max FD Error (based on N/2 run) to R=10.00: 1.3e-7
#    Est Max MOC error for R>10.00: -4.59e-07 
#    Max interpolation error with cubic splines: 5.0894e-07
#    Est Max overall error after cubic interpolation: 1.1e-6
#    

# Run:S8000_G1x6
'S1.6' => [
[-4.50135843, 12.00017658, -2.99997732, -4.50245891, 0.99834838],
[-4.13094467, 10.88894903, -2.99993109, -4.13286360, 0.99711889],
[-3.84809547, 10.04042840, -2.99985509, -3.85102997, 0.99559199],
[-3.43277968, 8.79459884, -2.99950364, -3.43825760, 0.99176189],
[-3.12115287, 7.85997402, -2.99873872, -3.12990834, 0.98681450],
[-2.87603930, 7.12508992, -2.99737048, -2.88870732, 0.98089332],
[-2.66373443, 6.48895493, -2.99504138, -2.68118806, 0.97363204],
[-2.47998323, 5.93891456, -2.99142757, -2.50302645, 0.96513074],
[-2.31717138, 5.45227347, -2.98610766, -2.34665769, 0.95531391],
[-2.17011852, 5.01367237, -2.97857703, -2.20697105, 0.94408148],
[-2.03507650, 4.61209339, -2.96823142, -2.08031379, 0.93130380],
[-1.90932842, 4.23966554, -2.95436379, -1.96408551, 0.91683376],
[-1.79018687, 3.88871121, -2.93605299, -1.85580235, 0.90041490],
[-1.67496445, 3.55172852, -2.91203228, -1.75310696, 0.88163763],
[-1.56095269, 3.22143789, -2.88051801, -1.65379790, 0.85989517],
[-1.44304287, 2.88419262, -2.83806176, -1.55390970, 0.83378434],
[-1.31183091, 2.51563614, -2.77713199, -1.44666016, 0.80015737],
[-1.19411900, 2.19260637, -2.70924471, -1.35445082, 0.76591238],
[-1.08550312, 1.90224240, -2.63570731, -1.27312233, 0.73116150],
[-0.98293389, 1.63583322, -2.55776968, -1.19991337, 0.69600873],
[-0.88560956, 1.39075087, -2.47781848, -1.13386413, 0.66108767],
[-0.79020674, 1.15825668, -2.39567317, -1.07246555, 0.62595880],
[-0.69491589, 0.93395713, -2.31188382, -1.01450276, 0.59058924],
[-0.59953598, 0.71745014, -2.22817673, -0.95985075, 0.55549446],
[-0.50309935, 0.50658490, -2.14539254, -0.90795931, 0.52086686],
[-0.40467320, 0.29944959, -2.06420886, -0.85837776, 0.48688264],
[-0.30348676, 0.09461492, -1.98528186, -0.81080470, 0.45375287],
[-0.19862101, -0.10952373, -1.90902858, -0.76492667, 0.42162273],
[-0.08933327, -0.31410187, -1.83590478, -0.72056333, 0.39067856],
[0.02532927, -0.52054499, -1.76614957, -0.67749421, 0.36103121],
[0.14626242, -0.73005665, -1.69999763, -0.63557134, 0.33279809],
[0.27445757, -0.94390715, -1.63759488, -0.59465565, 0.30606170],
[0.41095526, -1.16335152, -1.57904084, -0.55463385, 0.28088368],
[0.55705437, -1.38995726, -1.52432006, -0.51536107, 0.25727283],
[0.71296677, -1.62357894, -1.47375823, -0.47699667, 0.23538139],
[0.87465887, -1.85813044, -1.42858365, -0.44056301, 0.21575249],
[1.04364057, -2.09601936, -1.38800652, -0.40563725, 0.19805296],
[1.22105365, -2.33894591, -1.35147824, -0.37195519, 0.18204987],
[1.40824148, -2.58875965, -1.31850979, -0.33927072, 0.16753326],
[1.60632550, -2.84690808, -1.28873595, -0.30742453, 0.15434564],
[1.81672991, -3.11515449, -1.26180806, -0.27624553, 0.14233558],
[2.04064110, -3.39488738, -1.23746350, -0.24563224, 0.13138929],
[2.27954110, -3.68781449, -1.21544577, -0.21546817, 0.12139543],
[2.53477084, -3.99542188, -1.19554767, -0.18568006, 0.11226493],
[2.80831169, -4.31992090, -1.17754543, -0.15614491, 0.10390106],
[3.10185760, -4.66312722, -1.16127075, -0.12680024, 0.09623266],
[3.41747994, -5.02726091, -1.14655786, -0.09756767, 0.08919006],
[3.75724114, -5.41449234, -1.13326490, -0.06839321, 0.08271513],
[4.12359611, -5.82740438, -1.12125405, -0.03921113, 0.07675193],
[4.51919168, -6.26875693, -1.11040138, -0.00996494, 0.07125142],
[4.94726728, -6.74192892, -1.10058651, 0.01941908, 0.06616615],
[5.41105477, -7.25024633, -1.09170826, 0.04898581, 0.06145830],
[5.91437968, -7.79764612, -1.08366961, 0.07879167, 0.05709179],
[6.46166169, -8.38866409, -1.07638100, 0.10889751, 0.05303384],
[7.05751537, -9.02799905, -1.06976576, 0.13934350, 0.04925783],
[7.70735120, -9.72115794, -1.06375207, 0.17017962, 0.04573875],
[8.41717632, -10.47423359, -1.05827645, 0.20145005, 0.04245498],
[9.19319536, -11.29348117, -1.05328572, 0.23317515, 0.03938954],
[10.04261118, -12.18616553, -1.04873042, 0.26538492, 0.03652607],
[10.97262540, -13.15950685, -1.04457034, 0.29807832, 0.03385233],
[11.99163891, -14.22194204, -1.04076736, 0.33126828, 0.03135573],
[13.10825214, -15.38207659, -1.03728977, 0.36494499, 0.02902621],
[14.33206528, -16.64952018, -1.03410843, 0.39910253, 0.02685391],
[15.67307841, -18.03426316, -1.03119820, 0.43372032, 0.02483027],
[17.14209155, -19.54709292, -1.02853624, 0.46877549, 0.02294711],
[18.73195696, -21.18035577, -1.02612780, 0.50384392, 0.02121537],
[20.48510389, -22.97730324, -1.02390644, 0.53957689, 0.01959339],
[22.35075809, -24.88565708, -1.02192556, 0.57472591, 0.01812594],
[24.32024994, -26.89654770, -1.02016401, 0.60908862, 0.01680342],
[26.54886044, -29.16817087, -1.01848536, 0.64507614, 0.01552703],
[28.89074385, -31.55155644, -1.01699977, 0.68006477, 0.01438364],
[31.54144878, -34.24539101, -1.01558353, 0.71669041, 0.01328102],
[34.31025135, -37.05556739, -1.01433698, 0.75206899, 0.01229981],
[37.44292978, -40.23124254, -1.01314792, 0.78907863, 0.01135408],
[40.68983394, -43.51910171, -1.01210787, 0.82455267, 0.01051863],
[44.35886536, -47.23068354, -1.01111490, 0.86162896, 0.00971355],
[48.52521418, -51.44131079, -1.01016852, 0.90044010, 0.00893907],
[52.82187555, -55.77985071, -1.00934808, 0.93735431, 0.00826177],
[57.69127659, -60.69280963, -1.00856524, 0.97595350, 0.00761014],
[62.64968944, -65.69197905, -1.00789249, 1.01225201, 0.00704577],
[68.25033816, -71.33497936, -1.00724956, 1.05015007, 0.00650244],
[74.60759083, -77.73631198, -1.00663615, 1.08978097, 0.00598028],
[81.00681719, -84.17629306, -1.00611544, 1.12658472, 0.00553398],
[88.23596624, -91.44780111, -1.00561761, 1.16499625, 0.00510458],
[96.44328491, -99.69921797, -1.00514246, 1.20515044, 0.00469212],
[104.56971958, -107.86578484, -1.00474515, 1.24183448, 0.00434519],
[113.73421604, -117.07197726, -1.00436490, 1.28008488, 0.00401133],
[124.11956013, -127.50071110, -1.00400159, 1.32003164, 0.00369059],
[134.16193224, -137.58173345, -1.00370354, 1.35573282, 0.00342614],
[145.43789741, -148.89780522, -1.00341775, 1.39289220, 0.00317139],
[158.15625452, -161.65784209, -1.00314413, 1.43162660, 0.00292635],
[172.57176478, -176.11673682, -1.00288258, 1.47206755, 0.00269104],
[186.10509900, -189.68763380, -1.00267380, 1.50717875, 0.00250241],
[201.25384632, -204.87532757, -1.00247329, 1.54367489, 0.00232057],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 2.02e-07
#    Est Max FD Error (based on N/2 run) to R=10.00: 1.9e-7
#    Est Max MOC error for R>10.00: 3.23e-07 
#    Max interpolation error with cubic splines: 5.0598e-07
#    Est Max overall error after cubic interpolation: 1.e-6

# Run:S8000_G1x35
'S1.35' => [
[-4.65689366, 12.00017659, -2.99997883, -4.65795690, 0.99840429],
[-4.14572715, 10.46670063, -2.99990703, -4.14801744, 0.99656073],
[-3.71898990, 9.18656766, -2.99967672, -3.72333816, 0.99346409],
[-3.39482187, 8.21424313, -2.99913685, -3.40190211, 0.98934486],
[-3.13403177, 7.43221510, -2.99810906, -3.14451767, 0.98419782],
[-2.90739410, 6.75291523, -2.99627482, -2.92215275, 0.97772403],
[-2.71132563, 6.16570246, -2.99331416, -2.73117156, 0.96999715],
[-2.54045157, 5.65456906, -2.98889085, -2.56615236, 0.96108620],
[-2.38581390, 5.19283407, -2.98245909, -2.41830029, 0.95074534],
[-2.24525382, 4.77420546, -2.97352316, -2.28546105, 0.93897778],
[-2.11513871, 4.38804570, -2.96140830, -2.16412689, 0.92561505],
[-1.99274988, 4.02653419, -2.94529473, -2.05174216, 0.91044860],
[-1.87539994, 3.68208577, -2.92408107, -1.94588828, 0.89313865],
[-1.76040147, 3.34734719, -2.89623481, -1.84429664, 0.87317041],
[-1.64394733, 3.01212963, -2.85922008, -1.74395189, 0.84957268],
[-1.51784225, 2.65470391, -2.80728267, -1.63864214, 0.81988299],
[-1.39388577, 2.31057493, -2.74285905, -1.53904468, 0.78638137],
[-1.28081378, 2.00432128, -2.67226658, -1.45202740, 0.75221912],
[-1.17508322, 1.72569585, -2.59679980, -1.37430840, 0.71750682],
[-1.07452954, 1.46849107, -2.51797237, -1.30390699, 0.68248895],
[-0.97744710, 1.22793491, -2.43712101, -1.23934568, 0.64737809],
[-0.88113054, 0.99717771, -2.35424144, -1.17869920, 0.61187811],
[-0.78516784, 0.77525675, -2.27094084, -1.12168338, 0.57645340],
[-0.68866817, 0.56011875, -2.18819534, -1.06775581, 0.54135235],
[-0.59068993, 0.34974051, -2.10672458, -1.01641981, 0.50676986],
[-0.49033417, 0.14234944, -2.02712029, -0.96727701, 0.47289336],
[-0.38686915, -0.06334998, -1.94997576, -0.92007192, 0.43994760],
[-0.27941210, -0.26884171, -1.87566522, -0.87453030, 0.40809014],
[-0.16704401, -0.47554469, -1.80447764, -0.83042081, 0.37745830],
[-0.04902221, -0.68444664, -1.73675847, -0.78762690, 0.34822233],
[0.07571809, -0.89701357, -1.67262050, -0.74595383, 0.32045449],
[0.20807271, -1.11431428, -1.61224605, -0.70530938, 0.29425513],
[0.34925697, -1.33785365, -1.55564868, -0.66554052, 0.26964708],
[0.50064736, -1.56926821, -1.50279800, -0.62650191, 0.24662832],
[0.66029546, -1.80524573, -1.45461309, -0.58884723, 0.22560587],
[0.82628562, -2.04302656, -1.41146837, -0.55300279, 0.20674701],
[1.00008226, -2.28487637, -1.37265649, -0.51858546, 0.18974365],
[1.18292957, -2.53258038, -1.33765775, -0.48533276, 0.17436713],
[1.37602725, -2.78774736, -1.30605228, -0.45304282, 0.16043101],
[1.58068320, -3.05203650, -1.27747487, -0.42153844, 0.14777218],
[1.79814205, -3.32694850, -1.25162865, -0.39068858, 0.13625751],
[2.02966755, -3.61395278, -1.22825709, -0.36038701, 0.12577197],
[2.27696569, -3.91500998, -1.20710162, -0.33049846, 0.11619986],
[2.54132435, -4.23151890, -1.18797760, -0.30096571, 0.10745931],
[2.82464801, -4.56558109, -1.17068107, -0.27168295, 0.09946026],
[3.12881468, -4.91921644, -1.15504257, -0.24257480, 0.09212903],
[3.45589165, -5.29462537, -1.14090742, -0.21357108, 0.08539944],
[3.80813769, -5.69418983, -1.12813468, -0.18460785, 0.07921244],
[4.18800337, -6.12047297, -1.11659596, -0.15562825, 0.07351556],
[4.59813146, -6.57621941, -1.10617402, -0.12658298, 0.06826231],
[5.04215757, -7.06523358, -1.09674606, -0.09737975, 0.06340343],
[5.52330988, -7.59082199, -1.08821815, -0.06798376, 0.05890428],
[6.04561065, -8.15712153, -1.08049604, -0.03833665, 0.05472991],
[6.61367673, -8.76886537, -1.07349332, -0.00837758, 0.05084893],
[7.23252036, -9.43116258, -1.06713437, 0.02194248, 0.04723497],
[7.90755014, -10.14949953, -1.06135287, 0.05266058, 0.04386585],
[8.64497181, -10.93016408, -1.05608800, 0.08381831, 0.04072122],
[9.45138901, -11.77981749, -1.05128769, 0.11544127, 0.03778442],
[10.33400391, -12.70570760, -1.04690613, 0.14754765, 0.03504098],
[11.30081771, -13.71587707, -1.04290245, 0.18015378, 0.03247783],
[12.36023105, -14.81874436, -1.03924193, 0.21325942, 0.03008416],
[13.52104421, -16.02310876, -1.03589459, 0.24685038, 0.02785066],
[14.79325734, -17.33898091, -1.03283230, 0.28092171, 0.02576782],
[16.18727047, -18.77675139, -1.03003086, 0.31545323, 0.02382743],
[17.71468361, -20.34801784, -1.02746777, 0.35043058, 0.02202123],
[19.32603308, -22.00174952, -1.02520246, 0.38457353, 0.02039891],
[21.08888940, -23.80714501, -1.02311954, 0.41916379, 0.01888470],
[22.93851080, -25.69777739, -1.02127707, 0.45280482, 0.01752626],
[25.02855889, -27.83040765, -1.01952125, 0.48802795, 0.01621432],
[27.20782621, -30.05048689, -1.01797621, 0.52205664, 0.01504513],
[29.66876800, -32.55380199, -1.01650268, 0.55765616, 0.01391659],
[32.46099301, -35.39008920, -1.01509969, 0.59495406, 0.01282917],
[35.36335915, -38.33445199, -1.01387461, 0.63075847, 0.01186889],
[38.65271231, -41.66746693, -1.01270706, 0.66823604, 0.01094389],
[42.03797630, -45.09398696, -1.01169493, 0.70387705, 0.01013388],
[45.86616433, -48.96504343, -1.01072911, 0.74113712, 0.00935355],
[49.75549638, -52.89444680, -1.00989899, 0.77616657, 0.00867682],
[54.13819567, -57.31874303, -1.00910551, 0.81272704, 0.00802448],
[59.10002931, -62.32382346, -1.00834824, 0.85094271, 0.00739671],
[64.08120134, -67.34493502, -1.00770518, 0.88641463, 0.00685942],
[69.69275422, -72.99795233, -1.00709014, 0.92341564, 0.00634177],
[76.04422768, -79.39253958, -1.00650282, 0.96207022, 0.00584385],
[82.31356190, -85.70106140, -1.00601143, 0.99735942, 0.00542443],
[89.36128481, -92.78944789, -1.00554087, 1.03412778, 0.00502026],
[97.32021860, -100.79064003, -1.00509094, 1.07249403, 0.00463139],
[106.35263865, -109.86704696, -1.00466146, 1.11259196, 0.00425788],
[115.09942131, -118.65302273, -1.00430948, 1.14845416, 0.00395000],
[124.93470505, -128.52898799, -1.00397226, 1.18580703, 0.00365349],
[136.04503726, -139.68161078, -1.00364969, 1.22477121, 0.00336836],
[148.65854900, -152.33915818, -1.00334167, 1.26548269, 0.00309467],
[160.51910195, -164.23785115, -1.00309600, 1.30085424, 0.00287533],
[173.81612422, -177.57443237, -1.00286028, 1.33764624, 0.00266397],
[188.78931551, -192.58871395, -1.00263443, 1.37597109, 0.00246059],
[205.73037087, -209.57251733, -1.00241841, 1.41595515, 0.00226522],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 2.06e-07
#    Est Max FD Error (based on N/2 run) to R=10.00: 1.5e-7
#    Est Max MOC error for R>10.00: 2.53e-07 
#    Max interpolation error with cubic splines: 5.0813e-07
#    Est Max overall error after cubic interpolation: 9.e-7
#    

# Run:S8000_G1x45
'S1.45' => [
[-4.58381406, 12.00017659, -2.99997819, -4.58489327, 0.99838032],
[-4.22683732, 10.92925884, -2.99993635, -4.22868156, 0.99723113],
[-4.08339717, 10.49894872, -2.99990239, -4.08568464, 0.99656497],
[-3.58611558, 9.00720611, -2.99960363, -3.59094434, 0.99274026],
[-3.26194755, 8.03491869, -2.99895568, -3.26981123, 0.98816193],
[-3.00731609, 7.27142343, -2.99775996, -3.01885639, 0.98260169],
[-2.78784482, 6.61370457, -2.99568243, -2.80391524, 0.97573307],
[-2.59804696, 6.04541320, -2.99239711, -2.61945595, 0.96761901],
[-2.43153993, 5.54752989, -2.98753623, -2.45908554, 0.95827407],
[-2.28086192, 5.09786123, -2.98056132, -2.31547532, 0.94749972],
[-2.14313250, 4.68797065, -2.97092425, -2.18579181, 0.93523714],
[-2.01536887, 4.30917459, -2.95794865, -2.06716168, 0.92134537],
[-1.89465967, 3.95310529, -2.94074554, -1.95687088, 0.90557016],
[-1.77854145, 3.61287695, -2.91816915, -1.85273417, 0.88756794],
[-1.66410286, 3.28054393, -2.88852748, -1.75232273, 0.86674676],
[-1.54717974, 2.94502274, -2.84894088, -1.65239187, 0.84198509],
[-1.41684068, 2.57725961, -2.79181063, -1.54468467, 0.80994408],
[-1.29626636, 2.24449956, -2.72559008, -1.44902752, 0.77608216],
[-1.18549691, 1.94648800, -2.65341208, -1.36494292, 0.74160092],
[-1.08147201, 1.67438911, -2.57668416, -1.28959550, 0.70667204],
[-0.98259796, 1.42350371, -2.49726303, -1.22144132, 0.67169271],
[-0.88639338, 1.18715061, -2.41575749, -1.15850501, 0.63656203],
[-0.79054266, 0.95958525, -2.33237268, -1.09918880, 0.60108929],
[-0.69491589, 0.74054309, -2.24892286, -1.04339781, 0.56582917],
[-0.59841622, 0.52753369, -2.16616546, -0.99048655, 0.53094001],
[-0.50018799, 0.31877730, -2.08486531, -0.94002978, 0.49663735],
[-0.39941037, 0.11270033, -2.00565209, -0.89168420, 0.46312491],
[-0.29520060, -0.09226479, -1.92898253, -0.84513791, 0.43056853],
[-0.18675164, -0.29740740, -1.85527880, -0.80017029, 0.39914458],
[-0.07320885, -0.50399819, -1.78485646, -0.75658827, 0.36900005],
[0.04626862, -0.71318147, -1.71798256, -0.71424730, 0.34026955],
[0.17268855, -0.92629569, -1.65479048, -0.67298512, 0.31303267],
[0.30715431, -1.14472729, -1.59535574, -0.63265620, 0.28734220],
[0.45081841, -1.36983570, -1.53972940, -0.59314582, 0.26323626],
[0.60509171, -1.60327746, -1.48787088, -0.55431475, 0.24070815],
[0.76581904, -1.83857951, -1.44125336, -0.51729776, 0.22040613],
[0.93344582, -2.07657778, -1.39942414, -0.48192069, 0.20213961],
[1.10913794, -2.31905513, -1.36179376, -0.44789081, 0.18565452],
[1.29428539, -2.56796236, -1.32783672, -0.41493476, 0.17072191],
[1.48998333, -2.82473768, -1.29717416, -0.38288423, 0.15717551],
[1.69756234, -3.09104940, -1.26945199, -0.35156916, 0.14485911],
[1.91840334, -3.36855167, -1.24437074, -0.32084929, 0.13364012],
[2.15390500, -3.65885720, -1.22167717, -0.29061328, 0.12340650],
[2.40545742, -3.96351975, -1.20115471, -0.26077605, 0.11406294],
[2.67469841, -4.28435100, -1.18259715, -0.23124602, 0.10551912],
[2.96332588, -4.62319118, -1.16582379, -0.20194885, 0.09769686],
[3.27341426, -4.98228048, -1.15065953, -0.17279607, 0.09052093],
[3.60702829, -5.36380408, -1.13695630, -0.14372552, 0.08392917],
[3.96642459, -5.77013203, -1.12457785, -0.11467997, 0.07786557],
[4.35445215, -6.20426383, -1.11338858, -0.08557921, 0.07227465],
[4.77355143, -6.66870332, -1.10328333, -0.05639822, 0.06711633],
[5.22755639, -7.16746137, -1.09413974, -0.02703931, 0.06234174],
[5.71989367, -7.70404865, -1.08586491, 0.00253609, 0.05791686],
[6.25478412, -8.28280144, -1.07836705, 0.03238761, 0.05380790],
[6.83684336, -8.90843528, -1.07156487, 0.06256562, 0.04998596],
[7.47128268, -9.58625821, -1.06538468, 0.09311963, 0.04642540],
[8.16390987, -10.32216593, -1.05976087, 0.12409467, 0.04310391],
[8.92053007, -11.12200622, -1.05463940, 0.15550465, 0.04000463],
[9.74854649, -11.99327156, -1.04996595, 0.18739865, 0.03710893],
[10.65516100, -12.94319289, -1.04569824, 0.21978340, 0.03440385],
[11.64817463, -13.97959090, -1.04179866, 0.25265927, 0.03187778],
[12.73638791, -15.11129130, -1.03823272, 0.28603157, 0.02951948],
[13.92900106, -16.34749775, -1.03497099, 0.31988970, 0.02731945],
[15.23581418, -17.69800333, -1.03198744, 0.35421515, 0.02526907],
[16.66782732, -19.17380968, -1.02925773, 0.38899643, 0.02335957],
[18.23644047, -20.78629976, -1.02676082, 0.42420736, 0.02158328],
[19.95164562, -22.54539078, -1.02448019, 0.45977497, 0.01993477],
[21.71169951, -24.34675235, -1.02251418, 0.49355937, 0.01849265],
[23.70258008, -26.38053752, -1.02064152, 0.52895241, 0.01709995],
[25.79376998, -28.51311188, -1.01898482, 0.56337249, 0.01585161],
[28.15960190, -30.92193672, -1.01740600, 0.59941130, 0.01464705],
[30.62954645, -33.43310963, -1.01601675, 0.63422971, 0.01357449],
[33.42248466, -36.26887808, -1.01469194, 0.67065889, 0.01254012],
[36.31545923, -39.20263410, -1.01353334, 0.70558051, 0.01162584],
[39.58235106, -42.51187997, -1.01242745, 0.74208223, 0.01074429],
[43.28925580, -46.26284087, -1.01137356, 0.78029358, 0.00989576],
[47.10903241, -50.12426838, -1.01046004, 0.81663753, 0.00915328],
[51.43427990, -54.49281841, -1.00958845, 0.85463932, 0.00843853],
[55.83463501, -58.93367339, -1.00883946, 0.89037381, 0.00781913],
[60.80025303, -63.94135745, -1.00812366, 0.92767946, 0.00722250],
[66.43081180, -69.61567790, -1.00744071, 0.96668573, 0.00664878],
[72.09241035, -75.31771713, -1.00686095, 1.00290353, 0.00615814],
[78.48107849, -81.74839702, -1.00630661, 1.04069630, 0.00568579],
[85.72523110, -89.03626293, -1.00577745, 1.08019448, 0.00523182],
[92.88890317, -96.23969549, -1.00533491, 1.11626990, 0.00484973],
[100.95716131, -104.34924096, -1.00491130, 1.15387508, 0.00448183],
[110.08724177, -113.52226083, -1.00450646, 1.19313524, 0.00412816],
[120.47193861, -123.95169052, -1.00412023, 1.23419173, 0.00378877],
[130.55091541, -134.07055826, -1.00380387, 1.27093376, 0.00350927],
[141.90993724, -145.47102034, -1.00350096, 1.30922719, 0.00324033],
[154.77324419, -158.37744524, -1.00321140, 1.34920038, 0.00298197],
[166.83860501, -170.48012343, -1.00298023, 1.38389866, 0.00277477],
[180.33217550, -184.01237114, -1.00275820, 1.41995822, 0.00257495],
[195.48673270, -199.20707078, -1.00254525, 1.45748438, 0.00238253],
[212.58432877, -216.34639264, -1.00234133, 1.49659527, 0.00219751],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 9.26e-07
#    Est Max FD Error (based on N/2 run) to R=10.00: 1.6e-7
#    Est Max MOC error for R>10.00: -4.7e-07 
#    Max interpolation error with cubic splines: 5.1265e-07
#    Est Max overall error after cubic interpolation: 1.1e-6

# Run:S4000_G5
'S5' => [
[-3.98450494, 12.00033703, -2.99996929, -3.98578557, 0.99807785],
[-3.70145664, 11.15120316, -2.99992821, -3.70341528, 0.99705921],
[-3.52253607, 10.61445509, -2.99987722, -3.52509843, 0.99615165],
[-3.38929735, 10.21475503, -2.99985422, -3.39242745, 0.99529773],
[-2.99764269, 9.03990028, -2.99952665, -3.00328190, 0.99151866],
[-2.69466456, 8.13119986, -2.99882574, -2.70356160, 0.98660028],
[-2.43468111, 7.35171193, -2.99744377, -2.44784644, 0.98013880],
[-2.23182306, 6.74385193, -2.99530554, -2.24970625, 0.97297731],
[-2.04663244, 6.18944071, -2.99184598, -2.07029591, 0.96418158],
[-1.88159557, 5.69606729, -2.98668930, -1.91198037, 0.95393406],
[-1.73783211, 5.26716151, -2.97964881, -1.77562211, 0.94263314],
[-1.60048073, 4.85854908, -2.96957856, -1.64703658, 0.92926203],
[-1.47395993, 4.48363676, -2.95612730, -1.53038577, 0.91424392],
[-1.35561591, 4.13478297, -2.93852713, -1.42315755, 0.89741666],
[-1.24320274, 3.80566694, -2.91584556, -1.32330806, 0.87855930],
[-1.13029967, 3.47808834, -2.88559276, -1.22533277, 0.85644888],
[-1.01631899, 3.15136268, -2.84569759, -1.12915223, 0.83060526],
[-0.89501737, 2.80935319, -2.79106699, -1.03027158, 0.79901212],
[-0.77611021, 2.48132506, -2.72402925, -0.93730766, 0.76397262],
[-0.66717127, 2.18846688, -2.65066830, -0.85598094, 0.72861583],
[-0.56528284, 1.92230298, -2.57252583, -0.78353091, 0.69319794],
[-0.46872671, 1.67777610, -2.49145916, -0.71828334, 0.65809999],
[-0.37557157, 1.44951439, -2.40863275, -0.65858944, 0.62341852],
[-0.28277788, 1.22993995, -2.32365845, -0.60235294, 0.58867424],
[-0.18977793, 1.01782425, -2.23806225, -0.54921277, 0.55423508],
[-0.09595100, 0.81183544, -2.15314439, -0.49880679, 0.52040265],
[-0.00033263, 0.60997526, -2.06969186, -0.45064139, 0.48731311],
[0.09797292, 0.41055061, -1.98838948, -0.40433512, 0.45510406],
[0.19963742, 0.21243991, -1.90995840, -0.35966944, 0.42396496],
[0.30580676, 0.01371952, -1.83464044, -0.31627454, 0.39392719],
[0.41710028, -0.18640800, -1.76298195, -0.27405891, 0.36516911],
[0.53465896, -0.38959570, -1.69510355, -0.23277224, 0.33771888],
[0.65961225, -0.59732396, -1.63113163, -0.19223536, 0.31162011],
[0.79302019, -0.81084246, -1.57121300, -0.15234391, 0.28693428],
[0.93611186, -1.03157941, -1.51539594, -0.11298641, 0.26368894],
[1.08595901, -1.25478695, -1.46500781, -0.07510085, 0.24245839],
[1.24165417, -1.47929220, -1.42005056, -0.03887928, 0.22328022],
[1.40443270, -1.70708425, -1.37979438, -0.00398384, 0.20588091],
[1.57551569, -1.93996999, -1.34365025, 0.02985119, 0.19003950],
[1.75633636, -2.17990783, -1.31110190, 0.06287330, 0.17556000],
[1.94814941, -2.42850356, -1.28176327, 0.09524486, 0.16229869],
[2.15087537, -2.68561143, -1.25546911, 0.12689144, 0.15021009],
[2.36781407, -2.95531638, -1.23166186, 0.15823874, 0.13906512],
[2.59948803, -3.23810414, -1.21020461, 0.18924052, 0.12882519],
[2.84777348, -3.53610689, -1.19083794, 0.22002536, 0.11939300],
[3.11394983, -3.85068957, -1.17338575, 0.25061989, 0.11070943],
[3.40029786, -4.18436440, -1.15763124, 0.28114427, 0.10269320],
[3.70908860, -4.53956770, -1.14339909, 0.31168118, 0.09528107],
[4.04218470, -4.91822668, -1.13055591, 0.34224845, 0.08842973],
[4.40284549, -5.32381565, -1.11893670, 0.37296630, 0.08207682],
[4.79352046, -5.75884212, -1.10843083, 0.40385128, 0.07618702],
[5.21725317, -6.22644573, -1.09892534, 0.43494619, 0.07072187],
[5.67788040, -6.73059610, -1.09030911, 0.46632197, 0.06564148],
[6.18003184, -7.27607095, -1.08247985, 0.49806434, 0.06090819],
[6.72813060, -7.86736599, -1.07535887, 0.53020700, 0.05649588],
[7.32639422, -8.50871777, -1.06888185, 0.56274428, 0.05238519],
[7.98043561, -9.20582223, -1.06298018, 0.59571949, 0.04855179],
[8.69666366, -9.96517066, -1.05759277, 0.62917806, 0.04497363],
[9.48088412, -10.79256965, -1.05267441, 0.66310381, 0.04163678],
[10.33950057, -11.69442676, -1.04818294, 0.69748301, 0.03852739],
[11.28131508, -12.67962670, -1.04407270, 0.73236586, 0.03562699],
[12.31292871, -13.75471011, -1.04031580, 0.76768891, 0.03292762],
[13.44394199, -14.92931841, -1.03687786, 0.80347072, 0.03041498],
[14.68355516, -16.21263756, -1.03373224, 0.83968575, 0.02807878],
[16.04176832, -17.61465328, -1.03085450, 0.87630852, 0.02590896],
[17.52998148, -19.14676634, -1.02822133, 0.91332650, 0.02389505],
[19.12804778, -20.78798642, -1.02585581, 0.95000742, 0.02206155],
[20.88048492, -22.58377347, -1.02368218, 0.98713313, 0.02035570],
[22.68648315, -24.43079128, -1.02179628, 1.02250834, 0.01885854],
[24.73164266, -26.51863266, -1.01999521, 1.05955905, 0.01741329],
[26.81244720, -28.63938799, -1.01844603, 1.09445298, 0.01615757],
[29.16065497, -31.02912218, -1.01696415, 1.13093299, 0.01494508],
[31.82339595, -33.73509511, -1.01554891, 1.16912888, 0.01377636],
[34.50087179, -36.45254951, -1.01434635, 1.20463945, 0.01277468],
[37.52222925, -39.51546525, -1.01319551, 1.24173764, 0.01180837],
[40.94818852, -42.98469350, -1.01209596, 1.28055429, 0.01087779],
[44.33609434, -46.41197895, -1.01117557, 1.31604435, 0.01009312],
[48.15180648, -50.26860541, -1.01029383, 1.35307481, 0.00933629],
[52.46968569, -54.62905935, -1.00945046, 1.39177125, 0.00860753],
[57.38091733, -59.58466806, -1.00864517, 1.43227546, 0.00790706],
[62.14758570, -64.39092998, -1.00798502, 1.46855364, 0.00732933],
[67.51967218, -69.80416682, -1.00735247, 1.50639242, 0.00677267],
[73.60327287, -75.93060361, -1.00674736, 1.54592063, 0.00623720],
[79.30959975, -81.67401958, -1.00626392, 1.58025590, 0.00580725],
[85.69155461, -88.09442985, -1.00579932, 1.61595991, 0.00539217],
[92.85932525, -95.30212807, -1.00535345, 1.65313845, 0.00499202],
[100.94644968, -103.43077096, -1.00492620, 1.69191016, 0.00460688],
[110.11600528, -112.64356966, -1.00451749, 1.73240872, 0.00423680],
[120.56880314, -123.14148675, -1.00412721, 1.77478554, 0.00388186],
[130.02065903, -132.63090389, -1.00382820, 1.81015436, 0.00360883],
[140.60844988, -143.25766620, -1.00354086, 1.84693466, 0.00334555],
[152.52069403, -155.21040349, -1.00326516, 1.88523788, 0.00309205],
[165.98661661, -168.71846618, -1.00300104, 1.92518925, 0.00284836],
[181.28715455, -184.06293266, -1.00274845, 1.96693021, 0.00261450],
[194.17243501, -196.98242901, -1.00256655, 1.99950612, 0.00244558],
[208.46045186, -211.30583359, -1.00239108, 2.03325447, 0.00228221],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 1.38e-06
#    Est Max FD Error (based on N/2 run) to R=100.04: 3.e-7
#    Est Max MOC error for R>100.04: -1e-07 
#    Max interpolation error with cubic splines: 5.0688e-07
#    Est Max overall error after cubic interpolation: 9.e-7
#    

# Run:C4000_G5_rerun
'C5' => [
[-5.76912836, 12.00033752, -1.99997953, -5.77072940, 0.99839770],
[-5.29730654, 11.05670758, -1.99994740, -5.29987409, 0.99742922],
[-4.88101952, 10.22416246, -1.99989731, -4.88491526, 0.99609688],
[-4.25602915, 8.97430879, -1.99964547, -4.26331914, 0.99268480],
[-3.78241589, 8.02736292, -1.99908773, -3.79414622, 0.98820652],
[-3.42569441, 7.31439666, -1.99813557, -3.44249082, 0.98307904],
[-3.15675373, 6.77718033, -1.99680659, -3.17878252, 0.97776570],
[-2.88064728, 6.22614246, -1.99447076, -2.90976488, 0.97054377],
[-2.63521932, 5.73703417, -1.99101148, -2.67255311, 0.96214793],
[-2.41397586, 5.29703982, -1.98610712, -2.46070845, 0.95252306],
[-2.21144538, 4.89543100, -1.97937109, -2.26886420, 0.94157135],
[-2.02334366, 4.52390766, -1.97034208, -2.09288294, 0.92916412],
[-1.84576668, 4.17502135, -1.95843270, -1.92909695, 0.91509975],
[-1.67535548, 3.84254136, -1.94289634, -1.77447912, 0.89910912],
[-1.50844526, 3.51986184, -1.92269035, -1.62590048, 0.88076694],
[-1.33904505, 3.19632447, -1.89596871, -1.47848374, 0.85916914],
[-1.15537991, 2.85138649, -1.85864945, -1.32310702, 0.83215490],
[-0.97825079, 2.52602561, -1.81352819, -1.17828018, 0.80254354],
[-0.81619380, 2.23602426, -1.76422539, -1.05061404, 0.77260348],
[-0.66448411, 1.97229011, -1.71164012, -0.93566463, 0.74248578],
[-0.52027279, 1.72935273, -1.65687177, -0.83074199, 0.71245333],
[-0.38161771, 1.50346082, -1.60105290, -0.73400712, 0.68279302],
[-0.24359003, 1.28641268, -1.54377238, -0.64181763, 0.65302063],
[-0.10535610, 1.07700323, -1.48606379, -0.55359679, 0.62346136],
[0.03403740, 0.87386447, -1.42879538, -0.46872824, 0.59437759],
[0.17600942, 0.67504078, -1.37251098, -0.38638037, 0.56590134],
[0.32170583, 0.47910110, -1.31774903, -0.30596890, 0.53819886],
[0.47263506, 0.28425974, -1.26481937, -0.22678946, 0.51135484],
[0.63005829, 0.08920373, -1.21405959, -0.14835471, 0.48549749],
[0.79525949, -0.10730001, -1.16573994, -0.07022784, 0.46074564],
[0.96992538, -0.30684124, -1.11998028, 0.00815119, 0.43716048],
[1.15549061, -0.51059355, -1.07695529, 0.08716052, 0.41484566],
[1.35406734, -0.72035977, -1.03667840, 0.16740619, 0.39383123],
[1.56745418, -0.93747929, -0.99924779, 0.24929883, 0.37419736],
[1.78992594, -1.15596464, -0.96578870, 0.33053642, 0.35657023],
[2.02070292, -1.37533812, -0.93617799, 0.41097136, 0.34092232],
[2.26194494, -1.59792476, -0.90988625, 0.49149135, 0.32700476],
[2.51509520, -1.82522687, -0.88656467, 0.57266522, 0.31465731],
[2.78194725, -2.05896939, -0.86589310, 0.65513073, 0.30372954],
[3.06425499, -2.30075901, -0.84761342, 0.73947426, 0.29409959],
[3.36365023, -2.55204467, -0.83151646, 0.82622275, 0.28566794],
[3.68178240, -2.81426069, -0.81741944, 0.91589796, 0.27834573],
[4.02062930, -3.08909025, -0.80514856, 1.00910811, 0.27204601],
[4.38195784, -3.37803086, -0.79455729, 1.10640122, 0.26669318],
[4.76848954, -3.68333286, -0.78548806, 1.20858241, 0.26220369],
[5.18172882, -4.00628025, -0.77782548, 1.31613784, 0.25851235],
[5.62517098, -4.34972001, -0.77142224, 1.43008031, 0.25553617],
[6.10071668, -4.71525898, -0.76616512, 1.55101487, 0.25320623],
[6.61326531, -5.10681191, -0.76191263, 1.68031434, 0.25143919],
[7.16729924, -5.52795370, -0.75854348, 1.81924037, 0.25016013],
[7.76889024, -5.98345603, -0.75593601, 1.96945135, 0.24929379],
[8.42769948, -6.48078328, -0.75396857, 2.13349521, 0.24876634],
[9.15297498, -7.02706102, -0.75253621, 2.31381073, 0.24851083],
[9.96095678, -7.63465781, -0.75153244, 2.51457121, 0.24846270],
[10.87687598, -8.32266443, -0.75086111, 2.74218130, 0.24856476],
[11.93655592, -9.11808907, -0.75043978, 3.00568353, 0.24876693],
[13.19921443, -10.06546452, -0.75019672, 3.31995610, 0.24902596],
[14.77066555, -11.24424814, -0.75007230, 3.71151615, 0.24930520],
[16.86632211, -12.81607447, -0.75001927, 4.23427538, 0.24957351],
[20.03340228, -15.19141135, -0.75000276, 5.02510193, 0.24980269],
],

# Run:C8000_G1x6
'C1.6' => [
[-6.59721775, 12.00017708, -1.99998488, -6.59859354, 0.99862328],
[-5.99411966, 10.79399743, -1.99994949, -5.99663570, 0.99748085],
[-5.82022227, 10.44621259, -1.99993230, -5.82321689, 0.99700099],
[-5.18151608, 9.16888781, -1.99975682, -5.18719534, 0.99430530],
[-4.68221894, 8.17050163, -1.99934135, -4.69159214, 0.99058592],
[-4.28448325, 7.37543038, -1.99854254, -4.29846404, 0.98593158],
[-3.95494042, 6.71702405, -1.99718808, -3.97442497, 0.98035286],
[-3.66917986, 6.14658507, -1.99503682, -3.69517918, 0.97372710],
[-3.41757173, 5.64498782, -1.99183139, -3.45110684, 0.96603996],
[-3.19138195, 5.19493635, -1.98724983, -3.23355891, 0.95720486],
[-2.98467584, 4.78477179, -1.98091263, -3.03670406, 0.94712289],
[-2.79286243, 4.40557748, -1.97236313, -2.85609675, 0.93566212],
[-2.61213456, 4.05008426, -1.96104022, -2.68813779, 0.92264236],
[-2.43879704, 3.71138417, -1.94620518, -2.52946092, 0.90778537],
[-2.26865664, 3.38183189, -1.92678530, -2.37643146, 0.89062432],
[-2.09610295, 3.05148448, -1.90102501, -2.22445912, 0.87031679],
[-1.90820867, 2.69755605, -1.86476944, -2.06328573, 0.84462112],
[-1.72669697, 2.36293715, -1.82073787, -1.91250667, 0.81615655],
[-1.56072874, 2.06464647, -1.77258597, -1.77942218, 0.78712117],
[-1.40530725, 1.79306931, -1.72119208, -1.65935260, 0.75762349],
[-1.25763494, 1.54279849, -1.66769411, -1.54964757, 0.72793635],
[-1.11495258, 1.30872849, -1.61290259, -1.44789703, 0.69818268],
[-0.97285613, 1.08352616, -1.55663819, -1.35082810, 0.66801031],
[-0.83130161, 0.86717338, -1.50021593, -1.25840022, 0.63792425],
[-0.68866819, 0.65720455, -1.44418537, -1.16954849, 0.60806488],
[-0.54366038, 0.45180939, -1.38908310, -1.08352291, 0.57861605],
[-0.39493137, 0.24924938, -1.33531288, -0.99963216, 0.54973578],
[-0.24134056, 0.04820244, -1.28326437, -0.91738111, 0.52161256],
[-0.08160699, -0.15272368, -1.23320725, -0.83626395, 0.49440017],
[0.08554899, -0.35479868, -1.18536454, -0.75584109, 0.46824824],
[0.26146194, -0.55924976, -1.13990992, -0.67570326, 0.44329281],
[0.44767659, -0.76743839, -1.09694652, -0.59540093, 0.41963763],
[0.64573868, -0.98061781, -1.05657451, -0.51453763, 0.39738589],
[0.85763865, -1.20041186, -1.01881155, -0.43258662, 0.37659262],
[1.08540730, -1.42836578, -0.98369372, -0.34906222, 0.35731813],
[1.32283905, -1.65809918, -0.95227306, -0.26631445, 0.34016895],
[1.56982160, -1.88975438, -0.92435287, -0.18421812, 0.32505094],
[1.82826335, -2.12534728, -0.89950771, -0.10198227, 0.31173736],
[2.09982101, -2.36653126, -0.87741720, -0.01896263, 0.30005429],
[2.38623051, -2.61494211, -0.85781133, 0.06546768, 0.28985222],
[2.68901128, -2.87196712, -0.84047770, 0.15184524, 0.28101013],
[3.00995311, -3.13918565, -0.82521910, 0.24077015, 0.27341358],
[3.35080543, -3.41811304, -0.81186816, 0.33282197, 0.26696220],
[3.71316011, -3.71012516, -0.80028052, 0.42853869, 0.26156559],
[4.09948777, -4.01729733, -0.79029980, 0.52869261, 0.25712690],
[4.51104500, -4.34073708, -0.78181325, 0.63374621, 0.25356759],
[4.95041288, -4.68260942, -0.77468110, 0.74451282, 0.25079580],
[5.42030073, -5.04517345, -0.76877308, 0.86184058, 0.24872310],
[5.92400454, -5.43113758, -0.76395977, 0.98672703, 0.24726116],
[6.46560698, -5.84380474, -0.76011317, 1.12036528, 0.24632252],
[7.04957736, -6.28676223, -0.75711030, 1.26404307, 0.24582199],
[7.68357462, -6.76599802, -0.75482282, 1.41982997, 0.24567597],
[8.37684531, -7.28866853, -0.75313333, 1.59018221, 0.24580709],
[9.14302655, -7.86520839, -0.75193002, 1.77863541, 0.24614421],
[10.00154701, -8.51037294, -0.75111105, 1.99015714, 0.24662352],
[10.98142848, -9.24608849, -0.75058520, 2.23209760, 0.24718878],
[12.12928747, -10.10745296, -0.75027265, 2.51619001, 0.24779168],
[13.52553678, -11.15488656, -0.75010602, 2.86260535, 0.24839119],
[15.32598756, -12.50533431, -0.75003074, 3.31035922, 0.24895210],
[17.89685356, -14.43352058, -0.75000509, 3.95107495, 0.24944200],
[22.51756936, -17.89906433, -0.75000019, 5.10472023, 0.24982321],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 9.08e-07
#    Est Max FD Error (based on N/2 run) to R=10.00: 2.8e-7
#    Est Max MOC error for R>10.00: -2.04e-07 
#    Max interpolation error with cubic splines: 5.123e-07
#    Est Max overall error after cubic interpolation: 9.8e-7
#    

# Run:S4000_G2x5
'S2.5' => [
[-4.24630291, 12.00033634, -2.99997368, -4.24748848, 0.99822060],
[-3.93056073, 11.05312141, -2.99993213, -3.93246519, 0.99714066],
[-3.68244561, 10.30879676, -2.99987942, -3.68520990, 0.99584799],
[-3.25518598, 9.02712028, -2.99957469, -3.26043934, 0.99210036],
[-2.93026265, 8.05258383, -2.99887597, -2.93882855, 0.98710083],
[-2.70162948, 7.36705423, -2.99776642, -2.71371867, 0.98176982],
[-2.48240141, 6.71006194, -2.99569734, -2.49923153, 0.97457790],
[-2.29318007, 6.14349206, -2.99243434, -2.31558385, 0.96610133],
[-2.12590389, 5.64330309, -2.98756496, -2.15476612, 0.95625934],
[-1.97564644, 5.19488314, -2.98061984, -2.01189525, 0.94499047],
[-1.83860090, 4.78701567, -2.97104804, -1.88323448, 0.93220187],
[-1.71118438, 4.40922944, -2.95812846, -1.76535161, 0.91769186],
[-1.59129959, 4.05556277, -2.94105968, -1.65628825, 0.90129834],
[-1.47642318, 3.71892269, -2.91874272, -1.55378984, 0.88269619],
[-1.36367743, 3.39141631, -2.88955999, -1.45544383, 0.86132434],
[-1.24947281, 3.06352621, -2.85092729, -1.35847591, 0.83621763],
[-1.12537322, 2.71295081, -2.79671047, -1.25660844, 0.80475028],
[-1.00489854, 2.37986639, -2.73053737, -1.16170846, 0.77000951],
[-0.89479345, 2.08310314, -2.65815483, -1.07883241, 0.73488688],
[-0.79156143, 1.81261794, -2.58075763, -1.00477838, 0.69946806],
[-0.69392768, 1.56451656, -2.50054677, -0.93819288, 0.66429657],
[-0.59976817, 1.33289787, -2.41857276, -0.87727995, 0.62942206],
[-0.50589921, 1.10981017, -2.33435509, -0.81984368, 0.59432816],
[-0.41184835, 0.89425670, -2.24948802, -0.76559192, 0.55943177],
[-0.31703928, 0.68499339, -2.16529270, -0.71419009, 0.52506889],
[-0.22045644, 0.47989129, -2.08249166, -0.66511680, 0.49137630],
[-0.12147911, 0.27779697, -2.00196585, -0.61811958, 0.45860058],
[-0.01891888, 0.07652353, -1.92398934, -0.57273659, 0.42678026],
[0.08778074, -0.12471823, -1.84922624, -0.52885644, 0.39614418],
[0.19974606, -0.32770393, -1.77783584, -0.48617426, 0.36673853],
[0.31775749, -0.53344145, -1.71016133, -0.44457838, 0.33870113],
[0.44293471, -0.74343622, -1.64629822, -0.40387918, 0.31207789],
[0.57631004, -0.95892780, -1.58637634, -0.36396749, 0.28693421],
[0.71912307, -1.18139297, -1.53041716, -0.32471557, 0.26329215],
[0.87094803, -1.40974579, -1.47896470, -0.28644252, 0.24139687],
[1.02840882, -1.63891778, -1.43303780, -0.25001944, 0.22170447],
[1.19288702, -1.87114724, -1.39184982, -0.21505339, 0.20390254],
[1.36548779, -2.10810376, -1.35483068, -0.18128647, 0.18776505],
[1.54774828, -2.35191179, -1.32143225, -0.14843671, 0.17307025],
[1.74066120, -2.60384913, -1.29130421, -0.11637323, 0.15967930],
[1.94537582, -2.86533962, -1.26411821, -0.08496677, 0.14746092],
[2.16341547, -3.13821542, -1.23955182, -0.05406392, 0.13628438],
[2.39629697, -3.42422559, -1.21734136, -0.02354881, 0.12604360],
[2.64531667, -3.72479922, -1.19728121, 0.00663983, 0.11665791],
[2.91223433, -4.04188647, -1.17915580, 0.03659817, 0.10804125],
[3.19872911, -4.37729583, -1.16278688, 0.06638791, 0.10012419],
[3.50707476, -4.73348823, -1.14799048, 0.09610750, 0.09283348],
[3.83933527, -5.11263374, -1.13462192, 0.12580701, 0.08611410],
[4.19796800, -5.51731557, -1.12254064, 0.15554925, 0.07991234],
[4.58582280, -5.95051629, -1.11161548, 0.18540258, 0.07417824],
[5.00594188, -6.41538840, -1.10173103, 0.21542185, 0.06886901],
[5.46156019, -6.91525769, -1.09278480, 0.24564938, 0.06394769],
[5.95670617, -7.45427734, -1.08467584, 0.27615160, 0.05937697],
[6.49540190, -8.03654641, -1.07732058, 0.30696323, 0.05512788],
[7.08266418, -8.66719405, -1.07063686, 0.33814529, 0.05117073],
[7.72370511, -9.35150788, -1.06455590, 0.36973437, 0.04748199],
[8.42413308, -10.09515417, -1.05901746, 0.40175538, 0.04404162],
[9.18995362, -10.90417973, -1.05396865, 0.43422194, 0.04083238],
[10.02857017, -11.78606351, -1.04935791, 0.46717479, 0.03783600],
[10.94698474, -12.74781649, -1.04514584, 0.50060570, 0.03504014],
[11.95359840, -13.79787744, -1.04129398, 0.53452923, 0.03243131],
[13.05641169, -14.94422986, -1.03777220, 0.56891805, 0.02999998],
[14.26522486, -16.19669640, -1.03454992, 0.60377612, 0.02773473],
[15.59003800, -17.56527006, -1.03160164, 0.63908466, 0.02562623],
[17.04105115, -19.06012266, -1.02890530, 0.67480778, 0.02366634],
[18.60569923, -20.66803892, -1.02647381, 0.71039523, 0.02187180],
[20.26333122, -22.36771051, -1.02431060, 0.74527315, 0.02025255],
[22.14445713, -24.29257226, -1.02225029, 0.78186007, 0.01868943],
[24.03721906, -26.22575066, -1.02050397, 0.81593273, 0.01734785],
[26.17420629, -28.40472240, -1.01883663, 0.85158354, 0.01605200],
[28.59843222, -30.87262917, -1.01724725, 0.88894080, 0.01480249],
[30.99650919, -33.31043186, -1.01591968, 0.92314167, 0.01374766],
[33.69798952, -36.05315663, -1.01465033, 0.95886897, 0.01272918],
[36.75538377, -39.15343675, -1.01343855, 0.99624535, 0.01174743],
[40.23300215, -42.67572080, -1.01228367, 1.03540914, 0.01080280],
[43.60799498, -46.09053342, -1.01133859, 1.07051903, 0.01002293],
[47.41100230, -49.93489448, -1.01043444, 1.10716836, 0.00927083],
[51.71652443, -54.28342953, -1.00957084, 1.14548231, 0.00854670],
[55.75362470, -58.35774278, -1.00888185, 1.17878357, 0.00796480],
[60.26684051, -62.90951337, -1.00822052, 1.21343032, 0.00740259],
[65.33315655, -68.01582836, -1.00758662, 1.24952488, 0.00686020],
[71.04569799, -73.76992754, -1.00697992, 1.28718183, 0.00633774],
[77.51795873, -80.28543378, -1.00640020, 1.32653007, 0.00583531],
[84.88937891, -87.70193301, -1.00584724, 1.36771535, 0.00535304],
[91.54868865, -94.39873673, -1.00542399, 1.40209678, 0.00498179],
[99.00116684, -101.89008351, -1.00501761, 1.43785496, 0.00462354],
[107.37677058, -110.30603715, -1.00462800, 1.47509659, 0.00427836],
[116.83310856, -119.80432613, -1.00425502, 1.51394138, 0.00394628],
[127.56280101, -130.57770562, -1.00389858, 1.55452427, 0.00362736],
[136.58648115, -139.63538154, -1.00364203, 1.58619326, 0.00339682],
[146.57951439, -149.66354482, -1.00339467, 1.61899761, 0.00317372],
[157.68483902, -160.80521357, -1.00315646, 1.65301713, 0.00295808],
[170.07266779, -173.23068913, -1.00292733, 1.68834039, 0.00274992],
[183.94697080, -187.14404029, -1.00270725, 1.72506606, 0.00254925],
[199.55382838, -202.79145779, -1.00249617, 1.76330453, 0.00235608],
[217.19229654, -220.47212169, -1.00229404, 1.80317983, 0.00217044],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 5.58e-07
#    Est Max FD Error (based on N/2 run) to R=100.02: 4.e-7
#    Est Max MOC error for R>100.02: 8.39e-07 
#    Max interpolation error with cubic splines: 5.0508e-07
#    Est Max overall error after cubic interpolation: 1.7e-6

# Run:P8000_G1x5
'P1.5' => [
[-12.16767816, 12.00017797, -0.99999263, -12.16948986, 0.99909334],
[-11.39154943, 11.22405852, -0.99998398, -11.39422124, 0.99866234],
[-10.01863893, 9.85120033, -0.99993098, -10.02395377, 0.99733572],
[-8.98942474, 8.82211237, -0.99980364, -8.99833168, 0.99552764],
[-8.10081358, 7.93378220, -0.99952299, -8.11473570, 0.99299407],
[-7.37248988, 7.20596927, -0.99901338, -7.39258399, 0.98986275],
[-6.75477389, 6.58909546, -0.99817481, -6.78222671, 0.98611245],
[-6.21610805, 6.05172973, -0.99688446, -6.25217376, 0.98170343],
[-5.73680794, 5.57433996, -0.99499763, -5.78281695, 0.97659293],
[-5.30312570, 5.14336255, -0.99234459, -5.36050976, 0.97072820],
[-4.90443799, 4.74840392, -0.98872203, -4.97478061, 0.96403513],
[-4.53199618, 4.38101300, -0.98387923, -4.61710884, 0.95640758],
[-4.17830144, 4.03408950, -0.97749947, -4.28032722, 0.94769761],
[-3.83585495, 3.70070847, -0.96915380, -3.95745928, 0.93767879],
[-3.49559561, 3.37272594, -0.95819453, -3.64034538, 0.92596686],
[-3.14158842, 3.03603005, -0.94337513, -3.31500285, 0.91172590],
[-2.75571121, 2.67584062, -0.92265140, -2.96658785, 0.89365822],
[-2.40826550, 2.35914667, -0.89964492, -2.65925207, 0.87510316],
[-2.08668142, 2.07374212, -0.87481034, -2.38084092, 0.85612318],
[-1.78829054, 1.81648015, -0.84914633, -2.12817969, 0.83717690],
[-1.50632911, 1.58068557, -0.82315805, -1.89476374, 0.81835793],
[-1.22944396, 1.35641997, -0.79665498, -1.67080079, 0.79930474],
[-0.94932518, 1.13705901, -0.76955443, -1.44963453, 0.77976717],
[-0.65971500, 0.91820402, -0.74194816, -1.22672789, 0.75962660],
[-0.36459856, 0.70326758, -0.71490369, -1.00553383, 0.73949363],
[-0.06004413, 0.48958203, -0.68868514, -0.78339459, 0.71943110],
[0.25777213, 0.27476490, -0.66354810, -0.55794155, 0.69952523],
[0.59311840, 0.05632204, -0.63970975, -0.32669593, 0.67986173],
[0.95091931, -0.16847391, -0.61735598, -0.08695131, 0.66052819],
[1.33706524, -0.40275501, -0.59664622, 0.16439239, 0.64161501],
[1.74972840, -0.64502034, -0.57808940, 0.42536833, 0.62359000],
[2.16036631, -0.87917402, -0.56284591, 0.67811934, 0.60776324],
[2.57785449, -1.11141781, -0.55016508, 0.92884072, 0.59364633],
[3.00791778, -1.34566779, -0.53959236, 1.18135671, 0.58097126],
[3.45525135, -1.58500415, -0.53080375, 1.43862980, 0.56956261],
[3.92389606, -1.83199352, -0.52355242, 1.70308519, 0.55930203],
[4.41815614, -2.08923636, -0.51762993, 1.97718730, 0.55009427],
[4.94257256, -2.35937910, -0.51285788, 2.26344426, 0.54186467],
[5.50230441, -2.64532925, -0.50907773, 2.56463153, 0.53455027],
[6.10566244, -2.95154765, -0.50613707, 2.88513653, 0.52807456],
[6.75457014, -3.27922236, -0.50392886, 3.22591388, 0.52244326],
[7.46739552, -3.63781588, -0.50230572, 3.59650506, 0.51753673],
[8.25102768, -4.03096053, -0.50117756, 4.00035354, 0.51335953],
[9.13289182, -4.47257231, -0.50043722, 4.45143027, 0.50982652],
[10.11942116, -4.96603154, -0.50001038, 4.95289567, 0.50696163],
[11.28683801, -5.54961141, -0.49980804, 5.54327432, 0.50462311],
[12.63929501, -6.22554012, -0.49977089, 6.22447675, 0.50286875],
[14.35052756, -7.08081041, -0.49983391, 7.08376723, 0.50155383],
[16.80494551, -8.30775036, -0.49994030, 8.31348242, 0.50061891],
[19.93826474, -9.87433430, -0.49999697, 9.88124064, 0.50016806],
[23.34046992, -11.57544077, -0.50000185, 11.58263304, 0.50003435],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 5.38e-07
#    Est Max FD Error (based on N/2 run) to R=100.04: 3.8e-7
#    Est Max MOC error for R>100.04: 7.92e-07 
#    Max interpolation error with cubic splines: 5.0482e-07
#    Est Max overall error after cubic interpolation: 1.7e-7
#    

# Run:P8000_G1x6
'P1.6' => [
[-11.99798605, 12.00017796, -0.99999244, -11.99982085, 0.99908177],
[-11.14124491, 11.14344759, -0.99998219, -11.14406225, 0.99858938],
[-9.88456625, 9.88681737, -0.99993239, -9.88985372, 0.99734947],
[-8.76817818, 8.77056851, -0.99979317, -8.77743569, 0.99535086],
[-7.89581448, 7.89849285, -0.99950535, -7.91016798, 0.99277563],
[-7.17280245, 7.17600502, -0.99898239, -7.19346540, 0.98957331],
[-6.55946078, 6.56352374, -0.99812572, -6.58763069, 0.98574572],
[-6.02360700, 6.02899599, -0.99680981, -6.06056511, 0.98124469],
[-5.54680650, 5.55413874, -0.99489076, -5.59389797, 0.97603375],
[-5.11499896, 5.12507905, -0.99219604, -5.17368178, 0.97005435],
[-4.71787351, 4.73173690, -0.98852169, -4.78975645, 0.96323320],
[-4.34668151, 4.36566472, -0.98361424, -4.43360875, 0.95546093],
[-3.99423658, 4.02007332, -0.97715920, -4.09837789, 0.94659378],
[-3.65272743, 3.68773675, -0.96871841, -3.77680098, 0.93639374],
[-3.31309300, 3.36052491, -0.95763460, -3.46074151, 0.92446820],
[-2.95939826, 3.02435359, -0.94264360, -3.13626165, 0.90996568],
[-2.57570821, 2.66650817, -0.92180862, -2.79053854, 0.89168154],
[-2.23009648, 2.35179387, -0.89870815, -2.48554181, 0.87293447],
[-1.90983402, 2.06787871, -0.87377299, -2.20899904, 0.85377291],
[-1.61269295, 1.81201447, -0.84803098, -1.95811634, 0.83468400],
[-1.33192792, 1.57754192, -0.82198574, -1.72640600, 0.81576184],
[-1.05587101, 1.35428128, -0.79540922, -1.50384294, 0.79661855],
[-0.77662512, 1.13595974, -0.76825439, -1.28412376, 0.77703216],
[-0.48770430, 0.91801009, -0.74059072, -1.06254196, 0.75687137],
[-0.19314830, 0.70388959, -0.71349617, -0.84257813, 0.73675583],
[0.11092011, 0.49097918, -0.68724114, -0.62161792, 0.71675346],
[0.42827578, 0.27693513, -0.66208586, -0.39732562, 0.69695422],
[0.76344863, 0.05909907, -0.63822932, -0.16703939, 0.67742842],
[1.12110871, -0.16507865, -0.61587793, 0.07177242, 0.65827877],
[1.50721291, -0.39876665, -0.59518652, 0.32226293, 0.63958999],
[1.91798968, -0.63934856, -0.57674224, 0.58128298, 0.62189955],
[2.32726867, -0.87219244, -0.56157772, 0.83256807, 0.60637554],
[2.74366957, -1.10331677, -0.54896066, 1.08211807, 0.59254114],
[3.17320448, -1.33677176, -0.53843537, 1.33390433, 0.58012006],
[3.61991315, -1.57526832, -0.52969752, 1.59049307, 0.56895469],
[4.08793296, -1.82142268, -0.52249756, 1.85436532, 0.55891995],
[4.58156814, -2.07783251, -0.51662721, 2.12798430, 0.54991748],
[5.10568495, -2.34730864, -0.51190576, 2.41403305, 0.54186548],
[5.66558249, -2.63282441, -0.50817604, 2.71535306, 0.53469974],
[6.26425574, -2.93614814, -0.50530882, 3.03351311, 0.52839598],
[6.91670166, -3.26508223, -0.50314194, 3.37638085, 0.52283065],
[7.62497694, -3.62085329, -0.50158423, 3.74491978, 0.51802635],
[8.41096871, -4.01463559, -0.50050917, 4.15038331, 0.51388158],
[9.28695901, -4.45274843, -0.49983088, 4.59893457, 0.51039159],
[10.28510257, -4.95143933, -0.49945928, 5.10685061, 0.50749275],
[11.45008857, -5.53320048, -0.49932165, 5.69661335, 0.50514118],
[12.82615203, -6.22030590, -0.49935411, 6.39037309, 0.50332236],
[14.51855840, -7.06552941, -0.49949642, 7.24094642, 0.50196532],
[16.94856762, -8.27957844, -0.49970935, 8.45932276, 0.50093713],
[20.62772475, -10.11849922, -0.49990689, 10.30096548, 0.50029814],
],

# Created with: /home/steve/bin/blast_table_fit.pl 
# Wed Mar 21 15:13:16 2018   Inspiron-3668

# Working Directory:
# /home/steve/doc/slh/slh/slh1601/sedov/S4000_G6/moc2_from_W_r1000
#    FD  file(s): ../FD_history_full.30531.txt
#    MOC file(s): shock_history_moc.txt, shock_history_.txt
#    
#    Settings:
#    FD run(s): ../FD_history_full.30531.txt
#    MOC file(s): shock_history_moc.txt, shock_history_.txt
#    MOC1 file(s): moc_restart_G_r1000_J_PR.txt.moc1.1 : from X=18.4207181494312
#    min ln(r): -4.6 
#    max ln(r): 200 
#    max Y : 12
#    Ymin for dYdX Smooththing: 10.5
#    ovp tol: 5e-07     
#    mum tol: 0.001   
#    rmct tol: 1e-06     
#    ovp lin tol: 0.001 
#    Rounding? 0 
#    Decimal places for ln(ovp): 8     
#    Use Cubic Spline for mum?: 0
#    Use Parabolic Smoothing for ovp?: 0
#    rmct re-integration 0=no, 1=trapezoidal, 2=parabolic: 0 
#    Recompute final rmct using cubic spline fit?: 0 
#    Use old table points?: 0
#    Use alternate solution? 0=no 1=yes: 0
#    
#    The initial table of values was constructed using:
#    0 from the high pressure analytic solution
#    30529 from FD calculation 
#    126566 from the MOC calculations
#    0.0366751823062473 alpha = sedov-taylor alpha parameter based on this run
#    0.99999994284261 lambda scale factor implied from average energy Eavg=0.999999725657932
#    0.99999994284261 lambda scale factor actually used to make this table
#    -5.7157391676678e-08 =ln(0.99999994284261) = X and Z shift 
#    Overlap of MOC and MOC1 from 18.4:
#    dY_MOC1_max   = 8.28e-08
#    dz_MOC1_max   = 0.000734
#    ddYdX_MOC1_max= 3.99e-07
#    
#    
#    96 points in the compacted table printed in format '%0.8f'
#    
#    Z end re-integration diff=0 (from 8.23572e+06 to 8.23572e+06)
#    
#    Accurate Interpolation Error Summary:
#    Y err = 5.0796e-07 at X=-1.9138
#    Z err= 1.4272e-06 at X=-3.1387
#    dYdX err = 2.9932e-05 at X=-3.4266
#    
#    Linear Interpolation Error Summary:
#    max_Y_lin_err = 0.00099994 at X=-0.67053
#    max_Z_lin_err = 0.00052013 at X=-0.7844
#    max_dYdX_lin_err = 0.0017176 at X=-0.90509
#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 9.53e-07
#    Est Max FD Error (based on N/2 run) to R=10.00: 1.4e-7
#    Est Max MOC error for R>10.00: -3.57e-07 
#    Max interpolation error with cubic splines: 5.0796e-07
#    Est Max overall error after cubic interpolation: 1.e-6

# Run:S4000_G6
'S6' => [
[-3.92667191, 12.00033721, -2.99996841, -3.92797072, 0.99805055],
[-3.64966973, 11.16934166, -2.99992749, -3.65163824, 0.99704439],
[-3.47254061, 10.63796785, -2.99987665, -3.47510896, 0.99614266],
[-3.34019761, 10.24095487, -2.99985084, -3.34333080, 0.99529308],
[-2.94339255, 9.05065142, -2.99951986, -2.94908123, 0.99144407],
[-2.60884020, 8.04727160, -2.99869478, -2.61825257, 0.98582116],
[-2.37640022, 7.35039247, -2.99737227, -2.38976242, 0.97984034],
[-2.18434978, 6.77492412, -2.99532395, -2.20220680, 0.97301705],
[-1.99132159, 6.19705945, -2.99169198, -2.01523238, 0.96380475],
[-1.83277871, 5.72311385, -2.98669242, -1.86318079, 0.95390754],
[-1.68364091, 5.27819440, -2.97933384, -1.72176223, 0.94212723],
[-1.54808098, 4.87495456, -2.96926905, -1.59491776, 0.92883340],
[-1.42491914, 4.51001823, -2.95610616, -1.48139658, 0.91416476],
[-1.30657512, 4.16116735, -2.93849371, -1.37417862, 0.89732206],
[-1.19277633, 3.82801594, -2.91547092, -1.27312368, 0.87819610],
[-1.08058709, 3.50254650, -2.88529994, -1.17580385, 0.85617631],
[-0.96683033, 3.17649762, -2.84537667, -1.07984271, 0.83034448],
[-0.84418513, 2.83078643, -2.78991734, -0.97992317, 0.79833830],
[-0.72550190, 2.50352527, -2.72266283, -0.88721513, 0.76328847],
[-0.61678690, 2.21142803, -2.64911051, -0.80612914, 0.72795094],
[-0.51512239, 1.94601611, -2.57079966, -0.73390376, 0.69257988],
[-0.41879020, 1.70222961, -2.48958850, -0.66886380, 0.65755321],
[-0.32608292, 1.47523529, -2.40684006, -0.60949960, 0.62304722],
[-0.23315171, 1.25552287, -2.32141692, -0.55321562, 0.58827623],
[-0.14028928, 1.04393782, -2.23563227, -0.50018688, 0.55392755],
[-0.04646235, 0.83819115, -2.15041769, -0.44980737, 0.52014823],
[0.04904075, 0.63684304, -2.06679615, -0.40171939, 0.48716324],
[0.14712236, 0.43815956, -1.98543895, -0.35552617, 0.45509921],
[0.24878686, 0.24035964, -1.90680373, -0.31085710, 0.42403699],
[0.35484755, 0.04218199, -1.83140258, -0.26749303, 0.39411036],
[0.46614107, -0.15757829, -1.75962952, -0.22525247, 0.36543336],
[0.58369975, -0.36036757, -1.69168671, -0.18392995, 0.33806354],
[0.70859834, -0.56757835, -1.62772517, -0.14336214, 0.31205272],
[0.84206099, -0.78072945, -1.56781010, -0.10339244, 0.28742780],
[0.98537658, -1.00132277, -1.51198356, -0.06390041, 0.26421267],
[1.13488899, -1.22353813, -1.46181620, -0.02601321, 0.24308430],
[1.29002317, -1.44676106, -1.41713401, 0.01018350, 0.22401116],
[1.45302563, -1.67439658, -1.37695084, 0.04524612, 0.20661179],
[1.62376942, -1.90635354, -1.34102130, 0.07914305, 0.19081648],
[1.80340380, -2.14427707, -1.30881972, 0.11209698, 0.17643030],
[1.99427911, -2.39124753, -1.27974343, 0.14448119, 0.16321553],
[2.19790078, -2.64908094, -1.25346384, 0.17644630, 0.15105090],
[2.41482057, -2.91834366, -1.22980578, 0.20797098, 0.13988450],
[2.64631723, -3.20050474, -1.20851404, 0.23913665, 0.12962611],
[2.89442315, -3.49789295, -1.18930912, 0.27009511, 0.12017072],
[3.16081992, -3.81234666, -1.17198918, 0.30091688, 0.11144793],
[3.44738419, -4.14589263, -1.15636794, 0.33167031, 0.10339294],
[3.75618692, -4.50074160, -1.14227528, 0.36241922, 0.09594695],
[4.08989371, -4.87973859, -1.12954306, 0.39325689, 0.08904947],
[4.45096174, -5.28544104, -1.11803880, 0.42422584, 0.08265651],
[4.84244217, -5.72103267, -1.10763166, 0.45539261, 0.07672260],
[5.26777892, -6.19008331, -1.09820390, 0.48682197, 0.07120685],
[5.73020853, -6.69589010, -1.08966103, 0.51853358, 0.06607932],
[6.23436127, -7.24322988, -1.08190017, 0.55061188, 0.06130236],
[6.78446062, -7.83638280, -1.07484447, 0.58307773, 0.05685129],
[7.38572456, -8.48065745, -1.06841845, 0.61597919, 0.05269965],
[8.04336604, -9.18130922, -1.06256091, 0.64932954, 0.04882751],
[8.76339406, -9.94439989, -1.05721490, 0.68315146, 0.04521520],
[9.55201449, -10.77615777, -1.05233260, 0.71744517, 0.04184674],
[10.41643092, -11.68382391, -1.04786909, 0.75222399, 0.03870586],
[11.36384540, -12.67459859, -1.04378779, 0.78747079, 0.03577990],
[12.40265902, -13.75690152, -1.04005356, 0.82318564, 0.03305534],
[13.54107230, -14.93890964, -1.03663801, 0.85933459, 0.03052175],
[14.78868547, -16.23022496, -1.03351315, 0.89590483, 0.02816746],
[16.15569862, -17.64103828, -1.03065429, 0.93287487, 0.02598182],
[17.65351178, -19.18274800, -1.02803843, 0.97023011, 0.02395424],
[19.26016582, -20.83250868, -1.02569074, 1.00719515, 0.02211096],
[20.94571046, -22.55956082, -1.02361893, 1.04304049, 0.02046509],
[22.85793597, -24.51500500, -1.02164157, 1.08061390, 0.01887671],
[24.82801071, -26.52597750, -1.01992500, 1.11639414, 0.01748338],
[27.05834001, -28.79886435, -1.01828480, 1.15384789, 0.01613901],
[29.32693731, -31.10729835, -1.01687329, 1.18910274, 0.01497149],
[31.88643075, -33.70820192, -1.01552251, 1.22594183, 0.01384465],
[34.78803962, -36.65292419, -1.01423188, 1.26449485, 0.01275894],
[37.70503097, -39.60978671, -1.01313475, 1.30032116, 0.01182877],
[40.99595405, -42.94215874, -1.01208439, 1.33773339, 0.01093177],
[44.72677569, -46.71613800, -1.01108044, 1.37686226, 0.01006825],
[48.41542100, -50.44406377, -1.01023977, 1.41262398, 0.00934036],
[52.56907558, -54.63853182, -1.00943412, 1.44992457, 0.00863852],
[57.26850138, -59.38042972, -1.00866327, 1.48888953, 0.00796290],
[62.61272292, -64.76892387, -1.00792697, 1.52966060, 0.00731371],
[67.79875065, -69.99445407, -1.00732319, 1.56616610, 0.00677841],
[73.64253617, -75.87929738, -1.00674449, 1.60423082, 0.00626275],
[80.25926478, -82.53876799, -1.00619072, 1.64398355, 0.00576684],
[87.79040925, -90.11448438, -1.00566174, 1.68556989, 0.00529078],
[94.89010298, -97.25283801, -1.00523975, 1.72174456, 0.00490930],
[102.86997634, -105.27286484, -1.00483480, 1.75941522, 0.00454173],
[111.88077153, -114.32542852, -1.00444680, 1.79870326, 0.00418812],
[122.10685692, -124.59503445, -1.00407566, 1.83974561, 0.00384852],
[131.31380246, -133.83816153, -1.00379084, 1.87394620, 0.00358696],
[141.58436475, -144.14621497, -1.00351672, 1.90945776, 0.00333444],
[153.08799289, -155.68874342, -1.00325326, 1.94637981, 0.00309096],
[166.02920309, -168.67037447, -1.00300041, 1.98482370, 0.00285656],
[180.65665271, -183.33989072, -1.00275813, 2.02491460, 0.00263125],
[197.27507723, -200.00216903, -1.00252638, 2.06679386, 0.00241506],
[211.26880472, -214.03005543, -1.00235945, 2.09947264, 0.00225892],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 1.23e-06
#    Est Max FD Error (based on N/2 run) to R=100.08: 8.3e-7
#    Est Max MOC error for R>100.08: 1.92e-07 
#    Max interpolation error with cubic splines: 5.0368e-07
#    Est Max overall error after cubic interpolation: 1.5e-6

# Run:P4000_G6
'P6' => [
[-9.96706457, 12.00046263, -0.99998947, -9.96923005, 0.99891610],
[-9.00486523, 11.03827839, -0.99997245, -9.00837096, 0.99824411],
[-8.37943566, 10.41286963, -0.99995335, -8.38423147, 0.99759647],
[-7.07546810, 9.10902400, -0.99983265, -7.08469259, 0.99536741],
[-6.14888263, 8.18269372, -0.99957792, -6.16358079, 0.99260065],
[-5.40661457, 7.44088999, -0.99911454, -5.42798234, 0.98921340],
[-4.77618655, 6.81123901, -0.99834022, -4.80557475, 0.98511957],
[-4.22635836, 6.26262523, -0.99713371, -4.26519745, 0.98027120],
[-3.73901064, 5.77707291, -0.99535757, -3.78877870, 0.97463843],
[-3.29977288, 5.34038626, -0.99284936, -3.36204714, 0.96816831],
[-2.89677381, 4.94091828, -0.98940744, -2.97331475, 0.96076964],
[-2.52189095, 4.57082048, -0.98479760, -2.61466649, 0.95234704],
[-2.16825241, 4.22357326, -0.97873602, -2.27952196, 0.94277460],
[-1.82773381, 3.89157775, -0.97081864, -1.96030362, 0.93181558],
[-1.49259976, 3.56787645, -0.96048609, -1.65009008, 0.91913369],
[-1.15086474, 3.24188877, -0.94676956, -1.33850077, 0.90406579],
[-0.77071156, 2.88555617, -0.92708474, -0.99842341, 0.88461601],
[-0.41960447, 2.56390849, -0.90436647, -0.69135412, 0.86414181],
[-0.09693333, 2.27598813, -0.87964748, -0.41582451, 0.84337452],
[0.20618523, 2.01325849, -0.85342363, -0.16332224, 0.82245847],
[0.48666883, 1.77753865, -0.82711320, 0.06454267, 0.80223342],
[0.75875134, 1.55611551, -0.80035996, 0.28009320, 0.78217494],
[1.03136452, 1.34164717, -0.77303893, 0.49057338, 0.76201166],
[1.31205279, 1.12859814, -0.74509905, 0.70157890, 0.74155627],
[1.59997725, 0.91809197, -0.71733669, 0.91214753, 0.72125052],
[1.89547672, 0.71014589, -0.69038540, 1.12231972, 0.70143089],
[2.20315254, 0.50178251, -0.66443587, 1.33512704, 0.68213352],
[2.52688966, 0.29075143, -0.63974196, 1.55289035, 0.66346329],
[2.87112390, 0.07461690, -0.61651791, 1.77813183, 0.64551655],
[3.24128020, -0.14949113, -0.59493261, 2.01383358, 0.62837267],
[3.64424609, -0.38511407, -0.57511732, 2.26368662, 0.61209449],
[4.04635569, -0.61297527, -0.55873698, 2.50693052, 0.59809629],
[4.45405567, -0.83789995, -0.54510623, 2.74822811, 0.58592097],
[4.87136772, -1.06293069, -0.53377368, 2.99046005, 0.57527768],
[5.30865944, -1.29418506, -0.52425674, 3.23989559, 0.56580527],
[5.76354834, -1.53081477, -0.51644693, 3.49532744, 0.55748360],
[6.24166468, -1.77614823, -0.51008384, 3.76006047, 0.55013257],
[6.75200174, -2.03508659, -0.50494039, 4.03909000, 0.54357864],
[7.29012761, -2.30567586, -0.50094606, 4.33001102, 0.53783862],
[7.86674773, -2.59360095, -0.49790282, 4.63862697, 0.53275631],
[8.48755617, -2.90196714, -0.49568754, 4.96792549, 0.52826365],
[9.16715901, -3.23827755, -0.49417130, 5.32552583, 0.52425643],
[9.90310592, -3.60159232, -0.49327422, 5.71001292, 0.52074649],
[10.72001722, -4.00435892, -0.49288209, 6.13409133, 0.51761875],
[11.62972237, -4.45272212, -0.49291405, 6.60366227, 0.51484707],
[12.65969852, -4.96058396, -0.49329150, 7.13261615, 0.51237121],
[13.84978902, -5.54802097, -0.49394638, 7.74099869, 0.51013664],
[15.27141558, -6.25084854, -0.49482287, 8.46469529, 0.50808051],
[15.61040736, -6.41862493, -0.49504306, 8.63685927, 0.50766488],
[15.61520987, -6.42100244, -0.49507475, 8.63929732, 0.50765912],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 1.43e-06
#    Est Max FD Error (based on N/2 run) to R=100.04: 2.9e-7
#    Est Max MOC error for R>100.04: -1.11e-07 
#    Max interpolation error with cubic splines: 5.0504e-07
#    Est Max overall error after cubic interpolation: 9.e-7
#    

# Run:C4000_G6
'C6' => [
[-5.67047349, 12.00033772, -1.99997894, -5.67209726, 0.99837493],
[-5.20626531, 11.07193504, -1.99994671, -5.20884955, 0.99741249],
[-4.64419945, 9.94784997, -1.99985764, -4.64873727, 0.99545221],
[-4.10878133, 8.87714593, -1.99959447, -4.11654460, 0.99220824],
[-3.63427235, 7.92844596, -1.99896360, -3.64677705, 0.98742394],
[-3.25627746, 7.17304151, -1.99779903, -3.27457317, 0.98155820],
[-2.99696580, 6.65516687, -1.99630035, -3.02073271, 0.97599714],
[-2.72130720, 6.10520668, -1.99360887, -2.75271358, 0.96820710],
[-2.48371681, 5.63196450, -1.98977696, -2.52367313, 0.95946359],
[-2.26854770, 5.20436454, -1.98439982, -2.31826202, 0.94946666],
[-2.07014517, 4.81133849, -1.97704651, -2.13097753, 0.93807317],
[-1.88562634, 4.44738700, -1.96727051, -1.95903556, 0.92520822],
[-1.71051259, 4.10395788, -1.95441589, -1.79826051, 0.91060915],
[-1.54181093, 3.77558697, -1.93770209, -1.64600235, 0.89400709],
[-1.37487396, 3.45385404, -1.91586824, -1.49832055, 0.87482602],
[-1.20374894, 3.12838120, -1.88683075, -1.35051703, 0.85206860],
[-1.01251218, 2.77137539, -1.84509432, -1.19030712, 0.82276002],
[-0.84082617, 2.45846937, -1.79854934, -1.05155878, 0.79302207],
[-0.68262250, 2.17783037, -1.74810539, -0.92844245, 0.76302493],
[-0.53348471, 1.92104642, -1.69459661, -0.81687483, 0.73288349],
[-0.39169459, 1.68463797, -1.63943934, -0.71506465, 0.70303577],
[-0.25366692, 1.46222064, -1.58302045, -0.62006930, 0.67337673],
[-0.11543298, 1.24738203, -1.52521181, -0.52904616, 0.64359623],
[0.02284086, 1.04048732, -1.46743958, -0.44209180, 0.61422414],
[0.16257357, 0.83944867, -1.41034044, -0.35829025, 0.58540881],
[0.30521738, 0.64229493, -1.35440764, -0.27680976, 0.55726408],
[0.45211586, 0.44736932, -1.30007658, -0.19697921, 0.52991076],
[0.60452529, 0.25326842, -1.24773236, -0.11825619, 0.50347664],
[0.76390920, 0.05845641, -1.19762682, -0.04006672, 0.47805027],
[0.93163323, -0.13834860, -1.15000370, 0.03804040, 0.45373905],
[1.10936090, -0.33865615, -1.10499560, 0.11658822, 0.43061181],
[1.29853253, -0.54360873, -1.06277717, 0.19593885, 0.40877251],
[1.50139464, -0.75511000, -1.02334116, 0.27673294, 0.38824050],
[1.71871599, -0.97344789, -0.98695146, 0.35898359, 0.36918254],
[1.94291258, -1.19102580, -0.95485307, 0.43981061, 0.35228845],
[2.17617951, -1.41035464, -0.92641556, 0.52018891, 0.33726453],
[2.42014861, -1.63321012, -0.90120424, 0.60079721, 0.32391177],
[2.67677310, -1.86152849, -0.87883691, 0.68235542, 0.31205237],
[2.94765590, -2.09682827, -0.85903414, 0.76542098, 0.30155790],
[3.23444225, -2.34060418, -0.84155796, 0.85053705, 0.29231780],
[3.53891465, -2.59442204, -0.82619883, 0.93826828, 0.28423304],
[3.86271833, -2.85970334, -0.81278220, 1.02912892, 0.27721989],
[4.20783529, -3.13812887, -0.80114000, 1.12372409, 0.27119516],
[4.57632961, -3.43143044, -0.79112135, 1.22267861, 0.26608186],
[4.96967443, -3.74087283, -0.78259993, 1.32646530, 0.26181283],
[5.39134303, -4.06929367, -0.77541681, 1.43608927, 0.25830189],
[5.84404294, -4.41891391, -0.76944683, 1.55235110, 0.25547786],
[6.33127618, -4.79256764, -0.76456214, 1.67625954, 0.25326623],
[6.85633250, -5.19292455, -0.76064371, 1.80877227, 0.25159490],
[7.42509423, -5.62462423, -0.75756123, 1.95150115, 0.25038594],
[8.04483036, -6.09333386, -0.75519361, 2.10639823, 0.24956581],
[8.72359753, -6.60529538, -0.75342954, 2.27560665, 0.24906535],
[9.47564338, -7.17139410, -0.75215657, 2.46280637, 0.24881720],
[10.31860470, -7.80503027, -0.75127650, 2.67251409, 0.24876125],
[11.27971037, -8.52678145, -0.75069899, 2.91163119, 0.24884312],
[12.40258243, -9.36949830, -0.75034466, 3.19114325, 0.24901550],
[13.75783756, -10.38625368, -0.75014667, 3.52877503, 0.24923767],
[15.47708928, -11.67584667, -0.75004992, 3.95749075, 0.24947526],
[17.84315119, -13.45045508, -0.75001162, 4.54805116, 0.24969893],
[21.64824685, -16.30429424, -0.75000123, 5.49857799, 0.24988134],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 9.36e-07
#    Est Max FD Error (based on N/2 run) to R=10.00: 1.82e-7
#    Est Max MOC error for R>10.00: -2.26e-07 
#    Max interpolation error with cubic splines: 5.0825e-07
#    Est Max overall error after cubic interpolation: 9.1e-7

# Run:S4000_G4
'S4' => [
[-4.06012486, 12.00033680, -2.99997052, -4.06137960, 0.99811673],
[-3.77304581, 11.13911052, -2.99993025, -3.77497650, 0.99710123],
[-3.50813588, 10.34440347, -2.99987489, -3.51100987, 0.99568300],
[-3.08513094, 9.07549441, -2.99955969, -3.09055818, 0.99183824],
[-2.79827581, 8.21512963, -2.99895525, -2.80663259, 0.98741677],
[-2.52888727, 7.40739589, -2.99766587, -2.54142854, 0.98108469],
[-2.32152752, 6.78598772, -2.99565644, -2.33867790, 0.97409093],
[-2.13822007, 6.23712490, -2.99249414, -2.16084765, 0.96575947],
[-1.96556957, 5.72086622, -2.98747014, -1.99495931, 0.95545244],
[-1.81889501, 5.28314372, -2.98067178, -1.85560762, 0.94427819],
[-1.68050588, 4.87127605, -2.97100652, -1.72580441, 0.93118000],
[-1.55771000, 4.50716656, -2.95860831, -1.61230133, 0.91703323],
[-1.43768318, 4.15301525, -2.94167706, -1.50319698, 0.90048164],
[-1.32292205, 3.81663105, -2.91955875, -1.40090298, 0.88173974],
[-1.21073281, 3.49063081, -2.89072625, -1.30315349, 0.86030134],
[-1.09738189, 3.16502998, -2.85263427, -1.20702615, 0.83520318],
[-0.97556363, 2.82061231, -2.79977954, -1.10713114, 0.80414820],
[-0.85470987, 2.48610007, -2.73369977, -1.01202092, 0.76914089],
[-0.74453600, 2.18880319, -2.66126054, -0.92919903, 0.73383032],
[-0.64164319, 1.91888800, -2.58380067, -0.85549919, 0.69837278],
[-0.54423337, 1.67107746, -2.50316548, -0.78917581, 0.66314924],
[-0.45058834, 1.44049688, -2.42077141, -0.72869935, 0.62836142],
[-0.35765714, 1.21944389, -2.33629758, -0.67192274, 0.59354592],
[-0.26441559, 1.00558891, -2.25087209, -0.61819887, 0.55890895],
[-0.17051991, 0.79823954, -2.16606501, -0.56732649, 0.52486509],
[-0.07490154, 0.59514370, -2.08258549, -0.51874612, 0.49152080],
[0.02318007, 0.39490837, -2.00125956, -0.47214397, 0.45907250],
[0.12462065, 0.19593731, -1.92262967, -0.42718861, 0.42764194],
[0.23031590, -0.00322678, -1.84713152, -0.38361253, 0.39734080],
[0.34118779, -0.20396148, -1.77510253, -0.34119576, 0.36826724],
[0.45807957, -0.40739317, -1.70686255, -0.29979777, 0.34053215],
[0.58213224, -0.61506035, -1.64251067, -0.25922035, 0.31417067],
[0.71442053, -0.82826530, -1.58217242, -0.21934269, 0.28923750],
[0.85639255, -1.04879255, -1.52580717, -0.17998542, 0.26572373],
[1.00623969, -1.27348596, -1.47444548, -0.14182708, 0.24407622],
[1.16176067, -1.49913996, -1.42863115, -0.10542187, 0.22455513],
[1.32406063, -1.72759452, -1.38764221, -0.07044451, 0.20688769],
[1.49467672, -1.96112260, -1.35079491, -0.03655092, 0.19080932],
[1.67460167, -2.20110296, -1.31765083, -0.00357006, 0.17615556],
[1.86520866, -2.44933201, -1.28777574, 0.02869822, 0.16275954],
[2.06772835, -2.70732876, -1.26083395, 0.06038734, 0.15049448],
[2.28344228, -2.97661414, -1.23653369, 0.09160795, 0.13925080],
[2.51360846, -3.25862987, -1.21462501, 0.12244177, 0.12893617],
[2.76012872, -3.55555099, -1.19483919, 0.15302778, 0.11944698],
[3.02453568, -3.86904572, -1.17698179, 0.18342620, 0.11071281],
[3.30870812, -4.20115414, -1.16086142, 0.21371463, 0.10266310],
[3.61471813, -4.55409652, -1.14630509, 0.24396468, 0.09523458],
[3.94483085, -4.93027181, -1.13315716, 0.27424065, 0.08837095],
[4.30170472, -5.33248173, -1.12127153, 0.30461511, 0.08201896],
[4.68819061, -5.76369724, -1.11052009, 0.33514725, 0.07613288],
[5.10953361, -6.22948481, -1.10074304, 0.36603923, 0.07064667],
[5.56096555, -6.72436338, -1.09200319, 0.39676733, 0.06561973],
[6.05632219, -7.26325260, -1.08400211, 0.42807349, 0.06090237],
[6.59622469, -7.84648640, -1.07673390, 0.45973799, 0.05651022],
[7.18529096, -8.47874838, -1.07012494, 0.49178934, 0.05241854],
[7.82893411, -9.16553124, -1.06410616, 0.52426783, 0.04860325],
[8.53296325, -9.91270256, -1.05861812, 0.55719977, 0.04504411],
[9.30358436, -10.72650771, -1.05360937, 0.59059823, 0.04172384],
[10.14780117, -11.61399239, -1.04903312, 0.62447955, 0.03862615],
[11.07281584, -12.58237039, -1.04485019, 0.65883744, 0.03573787],
[12.08642955, -13.63944808, -1.04102553, 0.69366103, 0.03304679],
[13.19744286, -14.79403916, -1.03752668, 0.72894642, 0.03054074],
[14.41525604, -16.05554479, -1.03432521, 0.76468057, 0.02820883],
[15.74946919, -17.43354280, -1.03139677, 0.80083172, 0.02604170],
[17.21088235, -18.93882359, -1.02871817, 0.83737854, 0.02402956],
[18.72257513, -20.49212194, -1.02639287, 0.87233021, 0.02225858],
[20.35280199, -22.16360133, -1.02427615, 0.90725135, 0.02062603],
[22.19989372, -24.05361364, -1.02225644, 0.94385379, 0.01904931],
[24.08012356, -25.97401434, -1.02052059, 0.97834959, 0.01767876],
[26.20299469, -28.13863721, -1.01886146, 1.01443682, 0.01635486],
[28.61136895, -30.59047432, -1.01727817, 1.05224566, 0.01507825],
[31.03416766, -33.05346181, -1.01593380, 1.08741734, 0.01398371],
[33.76931364, -35.83038130, -1.01464814, 1.12418183, 0.01292748],
[36.87207648, -38.97663369, -1.01342060, 1.16267030, 0.01190998],
[39.94162890, -42.08576341, -1.01239376, 1.19787809, 0.01105178],
[43.40012307, -45.58537428, -1.01141063, 1.23463114, 0.01022384],
[47.31530294, -49.54333624, -1.01047085, 1.27305539, 0.00942641],
[51.77026365, -54.04288681, -1.00957405, 1.31329304, 0.00865977],
[56.09570516, -58.40810889, -1.00883931, 1.34934755, 0.00802734],
[60.97223854, -63.32598428, -1.00813567, 1.38696746, 0.00741786],
[66.49663213, -68.89341063, -1.00746290, 1.42628200, 0.00683150],
[71.68016321, -74.11419882, -1.00692568, 1.46044347, 0.00636062],
[77.47923077, -79.95189387, -1.00640963, 1.49597789, 0.00590597],
[83.99442030, -86.50718747, -1.00591460, 1.53299138, 0.00546763],
[91.34766729, -93.90213419, -1.00544048, 1.57160297, 0.00504569],
[99.68792640, -102.28582292, -1.00498713, 1.61194678, 0.00464022],
[109.19867586, -111.84188386, -1.00455444, 1.65417475, 0.00425129],
[117.80145766, -120.48238529, -1.00422308, 1.68942915, 0.00395210],
[127.44099470, -130.16105662, -1.00390478, 1.72609958, 0.00366359],
[138.28969730, -141.05042006, -1.00359949, 1.76429789, 0.00338578],
[150.55727236, -153.36030830, -1.00330714, 1.80414981, 0.00311871],
[164.50082685, -167.34797036, -1.00302767, 1.84579738, 0.00286241],
[176.24686013, -179.12835812, -1.00282648, 1.87830762, 0.00267727],
[189.27511817, -192.19214653, -1.00263247, 1.91199441, 0.00249823],
[203.77825776, -206.73207693, -1.00244561, 1.94694228, 0.00232529],
],

# Created with: /home/steve/bin/blast_table_fit.pl 
# Fri Mar 23 16:52:13 2018   Inspiron-3668

# Working Directory:
# /home/steve/doc/slh/slh/slh1601/sedov/P4000_G5/moc3p_from_AD_r100
#    FD  file(s): ../FD_history_full.26592.txt
#    MOC file(s): shock_history_.txt
#    
#    Settings:
#    FD run(s): ../FD_history_full.26592.txt
#    MOC file(s): shock_history_.txt
#    MOC1 file(s): 
#    min ln(r): -12.3 
#    max ln(r): 20 
#    max Y : 12
#    Ymin for dYdX Smooththing: 10.5
#    ovp tol: 5e-07     
#    mum tol: 0.001   
#    rmct tol: 1e-06     
#    ovp lin tol: 0.001 
#    Rounding? 0 
#    Decimal places for ln(ovp): 8     
#    Use Cubic Spline for mum?: 0
#    Use Parabolic Smoothing for ovp?: 0
#    rmct re-integration 0=no, 1=trapezoidal, 2=parabolic: 0 
#    Recompute final rmct using cubic spline fit?: 0 
#    Use old table points?: 0
#    Use alternate solution? 0=no 1=yes: 0
#    
#    The initial table of values was constructed using:
#    0 from the high pressure analytic solution
#    26591 from FD calculation 
#    65177 from the MOC calculations
#    0.0241337277100575 alpha = sedov-taylor alpha parameter based on this run
#    0.999999656805506 lambda scale factor implied from average energy Eavg=1.00000060224187
#    0.999999656805506 lambda scale factor actually used to make this table
#    -3.43194552936984e-07 =ln(0.999999656805506) = X and Z shift 
#    
#    
#    50 points in the compacted table printed in format '%0.8f'
#    
#    Z end re-integration diff=0 (from 6664.41 to 6664.41)
#    
#    Accurate Interpolation Error Summary:
#    Y err = 5.0376e-07 at X=-5.2405
#    Z err= 3.8537e-06 at X=-7.8541
#    dYdX err = 1.0555e-05 at X=-8.6819
#    
#    Linear Interpolation Error Summary:
#    max_Y_lin_err = 0.00099882 at X=2.5294
#    max_Z_lin_err = 0.00096862 at X=-1.1277
#    max_dYdX_lin_err = 0.00064928 at X=-1.1284
#    
#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 1.31e-06
#    Est Max FD Error (based on N/2 run) to R=100.08: 8.8e-7
#    Est Max MOC error for R>100.08: 2.86e-07 
#    Max interpolation error with cubic splines: 5.0376e-07
#    Est Max overall error after cubic interpolation: 1.7e-6
#    

# Run:P4000_G5
'P5' => [
[-10.18584956, 12.00046251, -0.99998976, -10.18798472, 0.99893130],
[-9.20740530, 11.02203345, -0.99997277, -9.21089015, 0.99825459],
[-8.53886419, 10.35351521, -0.99995019, -8.54373555, 0.99755852],
[-7.19303472, 9.00782055, -0.99981772, -7.20260391, 0.99519355],
[-6.28394378, 8.09899951, -0.99954831, -6.29905881, 0.99238946],
[-5.54292533, 7.35847670, -0.99905415, -5.56488694, 0.98891110],
[-4.93061665, 6.74696646, -0.99825833, -4.96054811, 0.98484174],
[-4.38641170, 6.20401887, -0.99700951, -4.42585805, 0.97995933],
[-3.90343761, 5.72289663, -0.99517837, -3.95387342, 0.97429443],
[-3.46732388, 5.28940894, -0.99259814, -3.53033492, 0.96778828],
[-3.06744883, 4.89315419, -0.98907236, -3.14477286, 0.96036643],
[-2.69506519, 4.52566676, -0.98436095, -2.78866998, 0.95192253],
[-2.34205049, 4.17921158, -0.97814964, -2.45427712, 0.94228944],
[-2.00267915, 3.84856107, -0.97006538, -2.13630612, 0.93128886],
[-1.66715899, 3.52477952, -0.95948580, -1.82592770, 0.91851150],
[-1.32381423, 3.19765964, -0.94540719, -1.51311137, 0.90328260],
[-0.93593190, 2.83476632, -0.92487906, -1.16652190, 0.88330820],
[-0.58818032, 2.51700084, -0.90193915, -0.86282921, 0.86292013],
[-0.26795301, 2.23206620, -0.87706312, -0.58976253, 0.84225710],
[0.03186470, 1.97296860, -0.85088550, -0.34031278, 0.82156552],
[0.31109869, 1.73901273, -0.82455306, -0.11369585, 0.80146461],
[0.58311151, 1.51835170, -0.79774585, 0.10159868, 0.78147136],
[0.85622088, 1.30421372, -0.77038713, 0.31227556, 0.76135352],
[1.13859945, 1.09064385, -0.74236572, 0.52436270, 0.74087888],
[1.42690633, 0.88063024, -0.71472388, 0.73502889, 0.72066622],
[1.72377684, 0.67248257, -0.68786637, 0.94600824, 0.70088519],
[2.03275767, 0.46398987, -0.66207886, 1.15955801, 0.68164312],
[2.35819237, 0.25258997, -0.63757090, 1.37831028, 0.66301396],
[2.70495902, 0.03558773, -0.61452856, 1.60505395, 0.64507350],
[3.07820616, -0.18968629, -0.59314678, 1.84255695, 0.62792209],
[3.48429607, -0.42646152, -0.57358290, 2.09416458, 0.61164614],
[3.88667034, -0.65392805, -0.55755612, 2.33740723, 0.59773689],
[4.29535526, -0.87897379, -0.54421384, 2.57914946, 0.58559964],
[4.71658043, -1.10577836, -0.53306855, 2.82350657, 0.57490510],
[5.15310911, -1.33638393, -0.52382843, 3.07235523, 0.56547800],
[5.61173645, -1.57480049, -0.51618382, 3.32972305, 0.55709938],
[6.09280161, -1.82156129, -0.50998587, 3.59589378, 0.54970561],
[6.60388800, -2.08087440, -0.50500997, 3.87511076, 0.54313438],
[7.14731936, -2.35419976, -0.50112456, 4.16863927, 0.53732536],
[7.73389342, -2.64722316, -0.49815974, 4.48225269, 0.53214896],
[8.36276196, -2.95978095, -0.49602564, 4.81542613, 0.52760070],
[9.03988768, -3.29512531, -0.49459312, 5.17128108, 0.52361300],
[9.78572064, -3.66365043, -0.49373379, 5.56043822, 0.52006697],
[10.61311311, -4.07197515, -0.49336305, 5.98938815, 0.51692528],
[11.54393950, -4.53119955, -0.49340057, 6.46920486, 0.51413525],
[12.59722847, -5.05106619, -0.49377063, 7.00938567, 0.51167161],
[13.81904348, -5.65473735, -0.49440615, 7.63314775, 0.50946796],
[15.28500757, -6.38013705, -0.49525153, 8.37846906, 0.50746187],
[16.11319291, -6.79049184, -0.49572301, 8.79835464, 0.50655002],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 1.31e-06
#    Est Max FD Error (based on N/2 run) to R=100.04: 3.1e-7
#    Est Max MOC error for R>100.04: 6.12e-08 
#    Max interpolation error with cubic splines: 5.0885e-07
#    

# Run:C4000_G4
'C4' => [
[-5.89542521, 12.00033729, -1.99998035, -5.89699387, 0.99843012],
[-5.40703252, 11.02356591, -1.99994780, -5.40959021, 0.99743910],
[-4.90945868, 10.02845615, -1.99987868, -4.91366877, 0.99578131],
[-4.41613952, 9.04191870, -1.99967569, -4.42304347, 0.99307338],
[-4.04643007, 8.30267561, -1.99932049, -4.05643678, 0.98994674],
[-3.60595854, 7.42220787, -1.99837769, -3.62154264, 0.98430770],
[-3.28372237, 6.77846885, -1.99691592, -3.30528733, 0.97823747],
[-3.00380909, 6.21979743, -1.99461975, -3.03242213, 0.97105884],
[-2.75703755, 5.72796673, -1.99122950, -2.79377287, 0.96276076],
[-2.53445051, 5.28524313, -1.98640667, -2.58049506, 0.95322919],
[-2.33089758, 4.88153585, -1.97977299, -2.38752819, 0.94238102],
[-2.14225083, 4.50884751, -1.97088466, -2.21087205, 0.93010621],
[-1.96400206, 4.15853097, -1.95913360, -2.04628738, 0.91616834],
[-1.79336692, 3.82547595, -1.94382325, -1.89126970, 0.90034886],
[-1.62586682, 3.50148225, -1.92384381, -1.74194651, 0.88214792],
[-1.45649336, 3.17777040, -1.89749202, -1.59430017, 0.86077722],
[-1.27320946, 2.83321900, -1.86074188, -1.43892260, 0.83407942],
[-1.09514259, 2.50573583, -1.81591989, -1.29297763, 0.80455812],
[-0.93234502, 2.21399729, -1.76688777, -1.16439475, 0.77467590],
[-0.77991705, 1.94859266, -1.71450357, -1.04858444, 0.74455985],
[-0.63514921, 1.70428890, -1.65992620, -0.94295755, 0.71450715],
[-0.49582234, 1.47687053, -1.60419211, -0.84547384, 0.68475324],
[-0.35743321, 1.25880753, -1.54707921, -0.75277599, 0.65491150],
[-0.21888893, 1.04845896, -1.48952044, -0.66410142, 0.62525515],
[-0.07915623, 0.84433546, -1.43235120, -0.57878497, 0.59603357],
[0.06292444, 0.64484305, -1.37621817, -0.49614711, 0.56743540],
[0.20892721, 0.44794715, -1.32148978, -0.41535604, 0.53954555],
[0.35999795, 0.25235342, -1.26861260, -0.33591218, 0.51252357],
[0.51742118, 0.05669553, -1.21790377, -0.25730692, 0.48649320],
[0.68268198, -0.14051355, -1.16956581, -0.17900360, 0.46154408],
[0.85728826, -0.34065229, -1.12376995, -0.10052856, 0.43776859],
[1.04274038, -0.54497686, -1.08066906, -0.02147237, 0.41526360],
[1.24098236, -0.75511819, -1.04030814, 0.05870386, 0.39407896],
[1.45392134, -0.97254274, -1.00275581, 0.14045981, 0.37428215],
[1.67706490, -1.19243413, -0.96897301, 0.22193413, 0.35641578],
[1.90842726, -1.41306581, -0.93906329, 0.30251373, 0.34057021],
[2.14997962, -1.63660410, -0.91251164, 0.38303335, 0.32650141],
[2.40357774, -1.86493971, -0.88892193, 0.46420571, 0.31402148],
[2.67095090, -2.09973312, -0.86798641, 0.54664691, 0.30298407],
[2.95347810, -2.34226869, -0.84947314, 0.63083418, 0.29327869],
[3.25302849, -2.59420771, -0.83315282, 0.71737307, 0.28479284],
[3.57116066, -2.85691056, -0.81884692, 0.80676415, 0.27743755],
[3.90996769, -3.13215701, -0.80637655, 0.89965233, 0.27112106],
[4.27111217, -3.42135978, -0.79559990, 0.99656139, 0.26576802],
[4.65638119, -3.72603654, -0.78637922, 1.09805676, 0.26130237],
[5.06881203, -4.04868692, -0.77856152, 1.20503541, 0.25763846],
[5.51084599, -4.39133067, -0.77202252, 1.31823850, 0.25470262],
[5.98498693, -4.75604397, -0.76664053, 1.43843166, 0.25242007],
[6.49553204, -5.14628084, -0.76228101, 1.56683856, 0.25070904],
[7.04636324, -5.56516605, -0.75882421, 1.70457682, 0.24949328],
[7.64415337, -6.01793517, -0.75614227, 1.85346024, 0.24869369],
[8.29716230, -6.51099873, -0.75411671, 2.01569108, 0.24823608],
[9.01543850, -7.05209327, -0.75263740, 2.19391118, 0.24805063],
[9.81482138, -7.65328927, -0.75159744, 2.39219555, 0.24807190],
[10.71774145, -8.33157791, -0.75090086, 2.61625372, 0.24824134],
[11.76002209, -9.11397367, -0.75046191, 2.87512695, 0.24850782],
[12.99748087, -10.04246147, -0.75020766, 3.18284550, 0.24882702],
[14.53153167, -11.19320138, -0.75007677, 3.56482604, 0.24916174],
[16.56518703, -12.71852967, -0.75002061, 4.07188164, 0.24948001],
[19.60986427, -15.00206526, -0.75000299, 4.83193201, 0.24975273],
[25.70261672, -19.57163436, -0.75000008, 6.35434386, 0.24994552],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 1.12e-06
#    Est Max FD Error (based on N/2 run) to R=61.00: 3.9e-7
#    Est Max MOC error for R>61.00: -3.83e-08 
#    Max interpolation error with cubic splines: 5.0678e-07
#    Est Max overall error after cubic interpolation: 9.3e-7
#    

# Run:C4000_G2x5
'C2.5' => [
[-6.19744871, 12.00033680, -1.99998245, -6.19893089, 0.99851673],
[-5.64232470, 10.89010467, -1.99994674, -5.64490830, 0.99741312],
[-4.85991111, 9.32536558, -1.99977980, -4.86556920, 0.99432655],
[-4.42906859, 8.46383338, -1.99947121, -4.43778646, 0.99124649],
[-4.01434908, 7.63473576, -1.99878927, -4.02757478, 0.98669508],
[-3.66882413, 6.94428809, -1.99758817, -3.68755381, 0.98111823],
[-3.37323571, 6.35408196, -1.99565643, -3.39847548, 0.97449849],
[-3.11392405, 5.83693166, -1.99273550, -3.14673327, 0.96677657],
[-2.88282765, 5.37686640, -1.98853839, -2.92429754, 0.95791888],
[-2.67210893, 4.95841970, -1.98268235, -2.72347555, 0.94778255],
[-2.47751322, 4.57332359, -1.97474492, -2.54012451, 0.93626791],
[-2.29541938, 4.21463968, -1.96422328, -2.37078621, 0.92323792],
[-2.12164922, 3.87445493, -1.95041624, -2.21160663, 0.90841406],
[-1.95280553, 3.54658713, -1.93242338, -2.05962130, 0.89144683],
[-1.78396184, 3.22221743, -1.90878939, -1.91073525, 0.87165731],
[-1.60675081, 2.88668348, -1.87670086, -1.75834836, 0.84760204],
[-1.42044053, 2.54087861, -1.83380507, -1.60307747, 0.81855624],
[-1.25171212, 2.23534813, -1.78640985, -1.46741360, 0.78902823],
[-1.09514259, 1.95955683, -1.73545585, -1.34618429, 0.75918155],
[-0.94645266, 1.70545088, -1.68169904, -1.23552157, 0.72907714],
[-0.80499728, 1.47140913, -1.62686449, -1.13448420, 0.69932557],
[-0.66560378, 1.24853700, -1.57062565, -1.03908128, 0.66945477],
[-0.52587107, 1.03306401, -1.51342139, -0.94763367, 0.63947438],
[-0.38564847, 0.82485291, -1.45647734, -0.86005049, 0.60984119],
[-0.24336612, 0.62164316, -1.40030354, -0.77536881, 0.58067582],
[-0.09774247, 0.42175951, -1.34540403, -0.69290452, 0.55213869],
[0.05239970, 0.22379682, -1.29221784, -0.61211068, 0.52439509],
[0.20847934, 0.02615929, -1.24101837, -0.53238510, 0.49755698],
[0.37164234, -0.17227492, -1.19211144, -0.45333703, 0.47177977],
[0.54339713, -0.37296332, -1.14564607, -0.37445757, 0.44715629],
[0.72539312, -0.57738964, -1.10172390, -0.29524622, 0.42376734],
[0.91909311, -0.78671247, -1.06048468, -0.21534295, 0.40172203],
[1.12660355, -1.00267963, -1.02193188, -0.13417517, 0.38106037],
[1.34981264, -1.22668735, -0.98613349, -0.05131663, 0.36185724],
[1.58066607, -1.45058241, -0.95440982, 0.03020515, 0.34485595],
[1.82087038, -1.67636219, -0.92624211, 0.11118372, 0.32980433],
[2.07243088, -1.90612902, -0.90118089, 0.19242556, 0.31647957],
[2.33707034, -2.14158429, -0.87889495, 0.27457558, 0.30471657],
[2.61631183, -2.38416694, -0.85912636, 0.35817744, 0.29438495],
[2.91182700, -2.63539032, -0.84164678, 0.44379626, 0.28536665],
[3.22526109, -2.89670379, -0.82626330, 0.53197380, 0.27755945],
[3.55817268, -3.16946298, -0.81281023, 0.62322234, 0.27087266],
[3.91269473, -3.45547787, -0.80111919, 0.71820958, 0.26521229],
[4.29121657, -3.75674132, -0.79103850, 0.81766645, 0.26049119],
[4.69481898, -4.07420841, -0.78246041, 0.92198801, 0.25664045],
[5.12619196, -4.41011969, -0.77524146, 1.03199997, 0.25357243],
[5.58844808, -4.76703334, -0.76924672, 1.14863588, 0.25120234],
[6.08408404, -5.14703065, -0.76435703, 1.27267776, 0.24945056],
[6.61758847, -5.55371984, -0.76043891, 1.40540960, 0.24823140],
[7.19403573, -5.99113807, -0.75736758, 1.54825780, 0.24746302],
[7.82108893, -6.46526573, -0.75501773, 1.70328642, 0.24706541],
[8.50559788, -6.98144652, -0.75327924, 1.87235487, 0.24696388],
[9.26180473, -7.55057496, -0.75203511, 2.05914599, 0.24708711],
[10.10894157, -8.18725898, -0.75118269, 2.26857636, 0.24737202],
[11.07503284, -8.91267634, -0.75063075, 2.50774709, 0.24776280],
[12.20369707, -9.75967678, -0.75029926, 2.78764480, 0.24821079],
[13.57054860, -10.78508359, -0.75011975, 3.12724085, 0.24867490],
[15.32039944, -12.09759451, -0.75003652, 3.56279813, 0.24912058],
[17.78166317, -13.94358579, -0.75000677, 4.17648345, 0.24951689],
[22.02156956, -17.12352490, -0.75000037, 5.23519175, 0.24983134],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 2.06e-07
#    Est Max FD Error (based on N/2 run) to R=50.01: 4.76e-7
#    Est Max MOC error for R>50.01: 2.54e-07 
#    Max interpolation error with cubic splines: 5.0745e-07
#    Est Max overall error after cubic interpolation: 1.2e-6
#    

# Run:C8000_G1x5
'C1.5' => [
[-6.02029767, 10.68580564, -1.99994512, -6.02292042, 0.99737412],
[-5.99163203, 10.62847312, -1.99994188, -5.99433115, 0.99729755],
[-5.94068334, 10.52657770, -1.99993564, -5.94352372, 0.99715589],
[-5.85255891, 10.35033777, -1.99992324, -5.85566133, 0.99689308],
[-5.74427895, 10.13379127, -1.99990469, -5.74773674, 0.99653659],
[-5.62961641, 9.90447995, -1.99988477, -5.63349507, 0.99611422],
[-5.11318707, 8.87169040, -1.99979055, -5.11969597, 0.99347090],
[-4.79058670, 8.22660493, -1.99942931, -4.79958419, 0.99096452],
[-4.63841249, 7.92237110, -1.99899180, -4.64889605, 0.98946563],
[-4.60795525, 7.86148911, -1.99900837, -4.61876464, 0.98913681],
[-4.45286968, 7.55148450, -1.99878190, -4.46550290, 0.98729448],
[-4.09242171, 6.83122674, -1.99751788, -4.11058147, 0.98169757],
[-3.79792710, 6.24323308, -1.99552901, -3.82236964, 0.97531273],
[-3.54016033, 5.72920473, -1.99254226, -3.57188026, 0.96789441],
[-3.30870772, 5.26848712, -1.98823402, -3.34881048, 0.95932878],
[-3.11342308, 4.88071538, -1.98277632, -3.16231724, 0.95033415],
[-2.92082584, 4.49954347, -1.97498764, -2.98029302, 0.93951882],
[-2.73729859, 4.13798851, -1.96449912, -2.80897406, 0.92705915],
[-2.56194552, 3.79465508, -1.95070872, -2.64762226, 0.91284306],
[-2.39107146, 3.46279778, -1.93267579, -2.49300089, 0.89648917],
[-2.21946288, 3.13308659, -1.90887504, -2.34076225, 0.87727661],
[-2.03716736, 2.78796663, -1.87612773, -2.18294902, 0.85353305],
[-1.84994493, 2.44055948, -1.83345725, -2.02572147, 0.82541413],
[-1.67996667, 2.13279661, -1.78644626, -1.88782036, 0.79666574],
[-1.52172472, 1.85401860, -1.73598242, -1.76404076, 0.76739627],
[-1.37125304, 1.59674133, -1.68289549, -1.65078363, 0.73769988],
[-1.22725301, 1.35828031, -1.62860202, -1.54667946, 0.70802768],
[-1.08445868, 1.12969979, -1.57270315, -1.44772354, 0.67788865],
[-0.94247420, 0.91039616, -1.51640517, -1.35361605, 0.64772313],
[-0.79991191, 0.69821991, -1.46038300, -1.26341930, 0.61773156],
[-0.65526293, 0.49099880, -1.40511550, -1.17622267, 0.58806129],
[-0.50731774, 0.28715123, -1.35107863, -1.09139554, 0.55890523],
[-0.35484669, 0.08519311, -1.29864376, -1.00837078, 0.53044016],
[-0.19668076, -0.11615959, -1.24812799, -0.92668257, 0.50284512],
[-0.03151743, -0.31824786, -1.19975042, -0.84585770, 0.47626848],
[0.14202113, -0.52238386, -1.15367738, -0.76544965, 0.45084372],
[0.32532443, -0.72978210, -1.11005382, -0.68506292, 0.42669954],
[0.51993722, -0.94173229, -1.06897692, -0.60428261, 0.40393957],
[0.72767218, -1.15970854, -1.03049255, -0.52263373, 0.38263622],
[0.95061468, -1.38535276, -0.99461801, -0.43958982, 0.36284086],
[1.18643745, -1.61595653, -0.96195748, -0.35619342, 0.34491944],
[1.43123305, -1.84779575, -0.93295199, -0.27374507, 0.32913050],
[1.68684939, -2.08288676, -0.90715251, -0.19144180, 0.31523279],
[1.95501679, -2.32299295, -0.88420747, -0.10859201, 0.30303376],
[2.23736314, -2.56968677, -0.86383626, -0.02458398, 0.29237642],
[2.53575399, -2.82467466, -0.84579080, 0.06123104, 0.28311964],
[2.85152701, -3.08916409, -0.82989113, 0.14933064, 0.27515642],
[3.18679851, -3.36498870, -0.81594286, 0.24040228, 0.26837161],
[3.54299456, -3.65338907, -0.80380548, 0.33493952, 0.26267592],
[3.92180400, -3.95582444, -0.79333717, 0.43351455, 0.25797788],
[4.32594575, -4.27456730, -0.78438376, 0.53696950, 0.25417998],
[4.75706466, -4.61103585, -0.77682896, 0.64587424, 0.25120043],
[5.21787483, -4.96749502, -0.77053958, 0.76107914, 0.24894908],
[5.71127551, -5.34635013, -0.76538844, 0.88348419, 0.24733772],
[6.24095322, -5.75060750, -0.76124690, 1.01418695, 0.24627768],
[6.81158147, -6.18401573, -0.75798745, 1.15452772, 0.24568169],
[7.42882127, -6.65105763, -0.75548652, 1.30608715, 0.24546528],
[8.10252321, -7.15935911, -0.75361808, 1.47147128, 0.24554827],
[8.84432642, -7.71785828, -0.75227046, 1.65372458, 0.24585751],
[9.67106158, -8.33936903, -0.75133902, 1.85717233, 0.24632639],
[10.60875229, -9.04357820, -0.75072826, 2.08841783, 0.24689644],
[11.69701627, -9.86034203, -0.75035497, 2.35745050, 0.24751660],
[13.00286711, -10.84003926, -0.75014778, 2.68109483, 0.24814353],
[14.64971604, -12.07532243, -0.75004811, 3.09026890, 0.24873988],
[16.90737508, -13.76862177, -0.75001010, 3.65249230, 0.24927255],
[20.57286727, -16.51775425, -0.75000077, 4.56711208, 0.24970680],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 1.42e-06
#    Est Max FD Error (based on N/2 run) to R=100.08: 9.7e-7
#    Est Max MOC error for R>100.08: -1.88e-06 
#    Max interpolation error with cubic splines: 5.0263e-07
#    Est Max overall error after cubic interpolation: 3.35e-6

# Run:P4000_G4
'P4' => [
[-10.46652274, 12.00046231, -0.99999017, -10.46861471, 0.99895293],
[-9.46496070, 10.99891565, -0.99997325, -9.46841478, 0.99827002],
[-8.96324247, 10.49721309, -0.99995868, -8.96768356, 0.99777463],
[-7.62241143, 9.15649554, -0.99984643, -7.63111185, 0.99563166],
[-6.67958103, 8.21390586, -0.99960615, -6.69355600, 0.99296696],
[-5.93731297, 7.47207122, -0.99917326, -5.95762653, 0.98975005],
[-5.29251445, 6.82802050, -0.99842816, -5.32065222, 0.98575959],
[-4.73581340, 6.27248835, -0.99726650, -4.77312438, 0.98105766],
[-4.24346724, 5.78187679, -0.99554978, -4.29139214, 0.97559202],
[-3.79985585, 5.34074128, -0.99311408, -3.85994961, 0.96930167],
[-3.39435756, 4.93866946, -0.98977134, -3.46830352, 0.96212320],
[-3.01760028, 4.56656138, -0.98528732, -3.10730606, 0.95395132],
[-2.66208637, 4.21727315, -0.97937400, -2.76976681, 0.94464986],
[-2.32031816, 3.88381035, -0.97164375, -2.44868353, 0.93400517],
[-1.98455985, 3.55918987, -0.96155202, -2.13709588, 0.92169439],
[-1.64341644, 3.23334607, -0.94818004, -1.82509270, 0.90710291],
[-1.26883145, 2.88157243, -0.92923286, -1.48870392, 0.88850059],
[-0.91379190, 2.55551883, -0.90672985, -1.17675622, 0.86835613],
[-0.58880517, 2.26472789, -0.88222180, -0.89782431, 0.84792128],
[-0.28357317, 1.99935946, -0.85612274, -0.64213243, 0.82726496],
[-0.00133824, 1.76139636, -0.82986799, -0.41145965, 0.80723103],
[0.27173595, 1.53840533, -0.80316670, -0.19373400, 0.78734242],
[0.54437357, 1.32313948, -0.77593155, 0.01819873, 0.76734901],
[0.82435580, 1.10980104, -0.74809155, 0.23019052, 0.74704147],
[1.11173281, 0.89883032, -0.72035314, 0.44194460, 0.72678727],
[1.40691254, 0.69022525, -0.69335101, 0.65351826, 0.70691684],
[1.71377993, 0.48150665, -0.66734684, 0.86743566, 0.68751313],
[2.03650646, 0.27020717, -0.64257166, 1.08622780, 0.66865961],
[2.37943570, 0.05393770, -0.61925003, 1.31235535, 0.65045669],
[2.74803732, -0.17021549, -0.59755054, 1.54882859, 0.63298064],
[3.14880906, -0.40558306, -0.57762572, 1.79909078, 0.61630998],
[3.55115014, -0.63454227, -0.56104084, 2.04406887, 0.60180558],
[3.95825818, -0.86004791, -0.54726494, 2.28642818, 0.58915005],
[4.37657984, -1.08649043, -0.53576908, 2.53048763, 0.57798949],
[4.81001777, -1.31656465, -0.52621817, 2.77881307, 0.56811536],
[5.26230362, -1.55271252, -0.51834333, 3.03372962, 0.55936280],
[5.73985138, -1.79863796, -0.51189178, 3.29893503, 0.55155802],
[6.24534956, -2.05602544, -0.50670888, 3.57594217, 0.54462493],
[6.77985687, -2.32572473, -0.50265546, 3.86536360, 0.53850919],
[7.35452254, -2.61363501, -0.49953897, 4.17321095, 0.53306176],
[7.97569436, -2.92317182, -0.49724072, 4.50278153, 0.52822596],
[8.64784576, -3.25681562, -0.49565339, 4.85634810, 0.52396805],
[9.38334371, -3.62096313, -0.49465934, 5.24029449, 0.52021562],
[10.19529334, -4.02236533, -0.49416084, 5.66129281, 0.51691935],
[11.11142182, -4.47500997, -0.49407180, 6.13346550, 0.51400297],
[12.15174979, -4.98911327, -0.49432160, 6.66680890, 0.51144769],
[13.35072629, -5.58208842, -0.49484002, 7.27861835, 0.50921035],
[14.78059894, -6.29016355, -0.49557022, 8.00522870, 0.50721873],
[16.61338544, -7.19929002, -0.49648483, 8.93306066, 0.50536912],
[16.83404501, -7.30885568, -0.49660599, 9.04455478, 0.50518526],
[16.83704578, -7.31034592, -0.49663503, 9.04607072, 0.50518276],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 1.72e-06
#    Est Max FD Error (based on N/2 run) to R=100.08: 1.2e-6
#    Est Max MOC error for R>100.08: 3.71e-07 
#    Max interpolation error with cubic splines: 5.0327e-07
#    Est Max overall error after cubic interpolation: 2.1e-6

# Run:P4000_G2x5
'P2.5' => [
[-11.13471153, 12.00046197, -0.99999123, -11.13668814, 0.99901073],
[-9.45523631, 10.32102358, -0.99995471, -9.45981953, 0.99770326],
[-8.26435842, 9.13025006, -0.99985003, -8.27268662, 0.99581928],
[-7.31652958, 8.18265795, -0.99961379, -7.32993895, 0.99325334],
[-6.55551737, 7.42208596, -0.99917434, -6.57519170, 0.99007544],
[-5.91196846, 6.77928197, -0.99843207, -5.93920062, 0.98622325],
[-5.35526741, 6.22374689, -0.99727341, -5.39137310, 0.98167854],
[-4.86229645, 5.73250887, -0.99555849, -4.90868167, 0.97639007],
[-4.41868506, 5.29136843, -0.99312802, -4.47683980, 0.97031262],
[-4.01256196, 4.88867110, -0.98978676, -4.08413394, 0.96336759],
[-3.63517988, 4.51594009, -0.98530254, -3.72201990, 0.95546018],
[-3.27904117, 4.16603295, -0.97938675, -3.38329837, 0.94645776],
[-2.93664815, 3.83195695, -0.97165157, -3.06095257, 0.93615223],
[-2.60050319, 3.50695985, -0.96156139, -2.74821937, 0.92423798],
[-2.25822434, 3.18003250, -0.94816377, -2.43423497, 0.91007661],
[-1.88021745, 2.82507943, -0.92905916, -2.09357138, 0.89188935],
[-1.52508965, 2.49900388, -0.90657713, -1.78023929, 0.87233477],
[-1.19959956, 2.20780906, -0.88209469, -1.49949258, 0.85244073],
[-0.89379815, 1.94197731, -0.85604896, -1.24186515, 0.83228176],
[-0.61129814, 1.70379672, -0.82990592, -1.00949393, 0.81269858],
[-0.33793119, 1.48053767, -0.80334642, -0.78998598, 0.79319746],
[-0.06444353, 1.26453608, -0.77622864, -0.57575101, 0.77348565],
[0.21680687, 1.05013322, -0.74849617, -0.36104595, 0.75335626],
[0.50541293, 0.83812374, -0.72089749, -0.14654802, 0.73319453],
[0.80248763, 0.62800281, -0.69399937, 0.06828242, 0.71327458],
[1.11110800, 0.41786996, -0.66813875, 0.28536503, 0.69373163],
[1.43565355, 0.20509639, -0.64352535, 0.50737128, 0.67463244],
[1.78063403, -0.01282444, -0.62037056, 0.73685203, 0.65606800],
[2.15191851, -0.23905328, -0.59882327, 0.97704048, 0.63809970],
[2.55625544, -0.47705687, -0.57903672, 1.23147463, 0.62080828],
[2.95984614, -0.70734431, -0.56268262, 1.47891559, 0.60573544],
[3.36811553, -0.93420647, -0.54911178, 1.72344940, 0.59248342],
[3.78795115, -1.16228079, -0.53778724, 1.96966027, 0.58069890],
[4.22317539, -1.39421418, -0.52838200, 2.22004973, 0.57019277],
[4.67851610, -1.63296767, -0.52061564, 2.47748530, 0.56079578],
[5.15710146, -1.88054519, -0.51428472, 2.74381064, 0.55240594],
[5.66370398, -2.13973205, -0.50919383, 3.02170830, 0.54491450],
[6.20716835, -2.41530483, -0.50515339, 3.31596773, 0.53819374],
[6.78345415, -2.70548054, -0.50208034, 3.62436476, 0.53228697],
[7.40409313, -3.01633572, -0.49979622, 3.95304747, 0.52706501],
[8.08940250, -3.35824425, -0.49815455, 4.31259247, 0.52239902],
[8.83749753, -3.73047275, -0.49708519, 4.70182074, 0.51834589],
[9.64021832, -4.12922353, -0.49648847, 5.11648170, 0.51493440],
[10.58748417, -4.59938285, -0.49624033, 5.60273255, 0.51185348],
[11.65724245, -5.13024871, -0.49629981, 6.14884075, 0.50927421],
[12.91603584, -5.75516081, -0.49660205, 6.78845806, 0.50709048],
[14.41868094, -6.50173444, -0.49708715, 7.54898636, 0.50527647],
[16.34782696, -7.46130512, -0.49771760, 8.52213089, 0.50372090],
[18.94944892, -8.75714151, -0.49842674, 9.83076142, 0.50239997],
[18.96705075, -8.76591490, -0.49847801, 9.83960452, 0.50239296],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 2.34e-07
#    Est Max FD Error (based on N/2 run) to R=100.03: 1.2e-7
#    Est Max MOC error for R>100.03: -5.8e-08 
#    Max interpolation error with cubic splines: 5.0665e-07
#    Est Max overall error after cubic interpolation: 6.8e-7
#    

# Run:C8000_G1x8
'C1.8' => [
[-6.47117602, 12.00017710, -1.99998420, -6.47258221, 0.99859284],
[-5.85643251, 10.77070777, -1.99994599, -5.85903433, 0.99739486],
[-5.43562994, 9.92913528, -1.99988444, -5.43959563, 0.99602668],
[-4.85761428, 8.77322971, -1.99963292, -4.86469364, 0.99289688],
[-4.40355509, 7.86538267, -1.99909108, -4.41472411, 0.98877370],
[-4.03672454, 7.13221331, -1.99810955, -4.05287914, 0.98373027],
[-3.72442588, 6.50843495, -1.99647839, -3.74655857, 0.97766170],
[-3.45299815, 5.96684645, -1.99396269, -3.48211446, 0.97054885],
[-3.21180370, 5.48632368, -1.99027476, -3.24897446, 0.96232247],
[-2.99367613, 5.05271735, -1.98507533, -3.04005456, 0.95290076],
[-2.79335261, 4.65572753, -1.97796360, -2.85020331, 0.94218012],
[-2.60646610, 4.28691114, -1.96844726, -2.67522423, 0.93000852],
[-2.42932143, 3.93926333, -1.95591128, -2.51166768, 0.91617630],
[-2.25824295, 3.60598253, -1.93953041, -2.35624576, 0.90036023],
[-2.08882457, 3.27912682, -1.91806863, -2.20522302, 0.88200736],
[-1.91403138, 2.94627089, -1.88927957, -2.05293120, 0.85998898],
[-1.71986649, 2.58327043, -1.84811572, -1.88863632, 0.83164404],
[-1.54550012, 2.26489600, -1.80225363, -1.74610191, 0.80270369],
[-1.38461475, 1.97884445, -1.75258572, -1.61929223, 0.77328861],
[-1.23277646, 1.71666586, -1.69996445, -1.50411675, 0.74350189],
[-1.08852645, 1.47530578, -1.64589223, -1.39899382, 0.71382546],
[-0.94732580, 1.24680124, -1.59038896, -1.30030295, 0.68396425],
[-0.80621918, 1.02637833, -1.53372701, -1.20591664, 0.65383120],
[-0.66504143, 0.81385312, -1.47713702, -1.11572883, 0.62389267],
[-0.52227313, 0.60698129, -1.42116596, -1.02878012, 0.59429591],
[-0.37656745, 0.40393900, -1.36629872, -0.94432314, 0.56520272],
[-0.22667182, 0.20317772, -1.31295114, -0.86175271, 0.53678043],
[-0.07141725, 0.00338323, -1.26148256, -0.78058160, 0.50920177],
[0.09058787, -0.19692420, -1.21211545, -0.70027410, 0.48259378],
[0.26056613, -0.39889210, -1.16508361, -0.62044475, 0.45710909],
[0.43995030, -0.60381783, -1.12052843, -0.54066114, 0.43286273],
[0.63039809, -0.81313773, -1.07852808, -0.46045052, 0.40994123],
[0.83356399, -1.02816957, -1.03916618, -0.37939808, 0.38843525],
[1.05147883, -1.25052343, -1.00246343, -0.29698893, 0.36839890],
[1.28275195, -1.47838800, -0.96891904, -0.21395487, 0.35014186],
[1.52256809, -1.70708718, -0.93915079, -0.13197054, 0.33402540],
[1.77280963, -1.93870009, -0.91267923, -0.05021383, 0.31980214],
[2.03537827, -2.17516193, -0.88912001, 0.03206203, 0.30727054],
[2.31176704, -2.41792980, -0.86819486, 0.11542230, 0.29628206],
[2.60368623, -2.66858573, -0.84965418, 0.20046840, 0.28670062],
[2.91285272, -2.92866222, -0.83328882, 0.28778058, 0.27840944],
[3.24084584, -3.19954369, -0.81892550, 0.37788959, 0.27130847],
[3.58953955, -3.48284016, -0.80639894, 0.47140437, 0.26530019],
[3.96073218, -3.78008741, -0.79556511, 0.56891362, 0.26029610],
[4.35643926, -4.09299338, -0.78628573, 0.67106845, 0.25620891],
[4.77914842, -4.42363576, -0.77842386, 0.77864632, 0.25295086],
[5.23072646, -4.77361025, -0.77186153, 0.89227379, 0.25044061],
[5.71386022, -5.14516066, -0.76646920, 1.01279332, 0.24859065],
[6.23284266, -5.54175776, -0.76210976, 1.14144859, 0.24731102],
[6.79155040, -5.96653946, -0.75866132, 1.27937805, 0.24651745],
[7.39564924, -6.42399417, -0.75599787, 1.42816088, 0.24612564],
[8.05419367, -6.92114954, -0.75399364, 1.59020646, 0.24605522],
[8.77762631, -7.46604765, -0.75253621, 1.76826316, 0.24623231],
[9.58158106, -8.07060985, -0.75151829, 1.96635901, 0.24658914],
[10.48968406, -8.75272969, -0.75084204, 2.19050151, 0.24706572],
[11.53715500, -9.53897014, -0.75042160, 2.44958403, 0.24760940],
[12.78320892, -10.47386685, -0.75018254, 2.75848253, 0.24817564],
[14.33425782, -11.63732836, -0.75006325, 3.14386412, 0.24872649],
[16.41191339, -13.19563970, -0.75001492, 3.66119473, 0.24922885],
[19.61899445, -15.60096943, -0.75000153, 4.46125525, 0.24965085],
[27.63119484, -21.61012181, -0.75000013, 6.46309191, 0.24995347],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 2.17e-07
#    Est Max FD Error (based on N/2 run) to R=10.00: 1.07e-7
#    Est Max MOC error for R>10.00: 1.47e-07 
#    Max interpolation error with cubic splines: 5.1138e-07
#    Est Max overall error after cubic interpolation: 7.5e-7

# Run:S8000_G1x8
'S1.8' => [
[-4.42008682, 12.00017657, -2.99997631, -4.42121161, 0.99831188],
[-4.05840712, 10.91515106, -2.99992988, -4.06034291, 0.99709356],
[-3.88966059, 10.40892477, -2.99990823, -3.89215463, 0.99625438],
[-3.54768851, 9.38306834, -2.99971923, -3.55185750, 0.99373405],
[-3.19933386, 8.33817681, -2.99920223, -3.20637348, 0.98940608],
[-2.92241934, 7.50777740, -2.99817131, -2.93310150, 0.98390050],
[-2.69264636, 6.81906277, -2.99636363, -2.70775317, 0.97719490],
[-2.49579406, 6.22948020, -2.99345611, -2.51613336, 0.96924469],
[-2.32368828, 5.71463374, -2.98908285, -2.35007880, 0.96003076],
[-2.16911074, 5.25304103, -2.98275633, -2.20246829, 0.94940781],
[-2.02846996, 4.83412177, -2.97394749, -2.06976358, 0.93730519],
[-1.89846683, 4.44822852, -2.96200394, -1.94877408, 0.92358004],
[-1.77641393, 4.08762279, -2.94612526, -1.83696865, 0.90803541],
[-1.65973583, 3.74503019, -2.92525582, -1.73202329, 0.89035952],
[-1.54594802, 3.41365151, -2.89796590, -1.63183409, 0.87009593],
[-1.43150943, 3.08398463, -2.86191806, -1.53358580, 0.84636837],
[-1.30947942, 2.73765436, -2.81213916, -1.43204728, 0.81708777],
[-1.18516099, 2.39190136, -2.74787449, -1.33254679, 0.78292386],
[-1.07217807, 2.08532470, -2.67715070, -1.24601842, 0.74823270],
[-0.96692144, 1.80744853, -2.60132968, -1.16908658, 0.71316193],
[-0.86724563, 1.55205622, -2.52206646, -1.09974043, 0.67800465],
[-0.77161886, 1.31472482, -2.44097581, -1.03656870, 0.64306167],
[-0.67668684, 1.08694235, -2.35756434, -0.97719328, 0.60779771],
[-0.58184390, 0.86733893, -2.27333193, -0.92121976, 0.57259755],
[-0.48652703, 0.65465568, -2.18962240, -0.86830580, 0.53782122],
[-0.38966853, 0.44659405, -2.10711880, -0.81788125, 0.50360482],
[-0.29049764, 0.24165972, -2.02658207, -0.76961119, 0.47016541],
[-0.18809534, 0.03817749, -1.94850556, -0.72314850, 0.43764943],
[-0.08171896, -0.16504756, -1.87340696, -0.67828528, 0.40624506],
[0.02956143, -0.36946083, -1.80157906, -0.63478257, 0.37606928],
[0.14657545, -0.57620625, -1.73332157, -0.59249160, 0.34725636],
[0.27041995, -0.78679126, -1.66873772, -0.55121375, 0.31986542],
[0.40199074, -1.00226930, -1.60801723, -0.51086547, 0.29399476],
[0.54260958, -1.22429646, -1.55113436, -0.47127357, 0.26965024],
[0.69348310, -1.45422604, -1.49812806, -0.43234988, 0.24686235],
[0.85103212, -1.68639839, -1.45036325, -0.39512122, 0.22623247],
[1.01519889, -1.92089671, -1.40754989, -0.35954390, 0.20765092],
[1.18718462, -2.15958036, -1.36906039, -0.32531077, 0.19085750],
[1.36830096, -2.40431526, -1.33436076, -0.29215703, 0.17562832],
[1.55974384, -2.65669129, -1.30303612, -0.25989252, 0.16178857],
[1.76273187, -2.91824258, -1.27473613, -0.22836269, 0.14918994],
[1.97860688, -3.19058937, -1.24915058, -0.19742831, 0.13770060],
[2.20877282, -3.47536494, -1.22601241, -0.16697318, 0.12720710],
[2.45465897, -3.77418045, -1.20509209, -0.13690497, 0.11761241],
[2.71786143, -4.08880421, -1.18617924, -0.10713530, 0.10882792],
[3.00002453, -4.42102013, -1.16908943, -0.07759363, 0.10077686],
[3.30322476, -4.77307840, -1.15364024, -0.04818831, 0.09338320],
[3.62932858, -5.14694316, -1.13968623, -0.01887264, 0.08658836],
[3.98099611, -5.54545003, -1.12707019, 0.01044800, 0.08032789],
[4.36047793, -5.97092595, -1.11567337, 0.03980711, 0.07455606],
[4.77101804, -6.42677527, -1.10536494, 0.06929150, 0.06922083],
[5.21545225, -6.91590387, -1.09604476, 0.09893016, 0.06428546],
[5.69760995, -7.44227462, -1.08760669, 0.12879389, 0.05970963],
[6.22171413, -8.01022855, -1.07995716, 0.15894486, 0.05545882],
[6.79198219, -8.62405509, -1.07301820, 0.18941470, 0.05150615],
[7.41362691, -9.28907197, -1.06671311, 0.22025885, 0.04782439],
[8.09225725, -10.01097005, -1.06097580, 0.25151830, 0.04439098],
[8.83387918, -10.79581609, -1.05574898, 0.28322027, 0.04118708],
[9.64529649, -11.65047598, -1.05098062, 0.31539432, 0.03819539],
[10.53371143, -12.58219023, -1.04662637, 0.34805370, 0.03540168],
[11.50692525, -13.59878678, -1.04264724, 0.38120402, 0.03279340],
[12.57333860, -14.70868130, -1.03900899, 0.41484272, 0.03035946],
[13.74195178, -15.92087806, -1.03568144, 0.44895927, 0.02808984],
[15.02256491, -17.24517760, -1.03263744, 0.48354075, 0.02597515],
[16.42557804, -18.69196949, -1.02985297, 0.51856578, 0.02400678],
[17.96279120, -20.27305610, -1.02730533, 0.55402360, 0.02217588],
[19.57153564, -21.92387816, -1.02506997, 0.58835223, 0.02054462],
[21.34396330, -23.73885178, -1.02299868, 0.62337104, 0.01901131],
[23.22348520, -25.65981230, -1.02114827, 0.65776456, 0.01762295],
[25.35099951, -27.83038928, -1.01938491, 0.69379346, 0.01628282],
[27.59961494, -30.12078071, -1.01781653, 0.72901486, 0.01507614],
[30.14736245, -32.71196017, -1.01632175, 0.76590217, 0.01391258],
[32.82830959, -35.43484148, -1.01499879, 0.80176811, 0.01287122],
[35.86728765, -38.51742781, -1.01373759, 0.83931678, 0.01186789],
[39.04665655, -41.73865918, -1.01262754, 0.87559343, 0.01097585],
[42.64987020, -45.38540922, -1.01156879, 0.91355133, 0.01011686],
[46.39103969, -49.16807914, -1.01064284, 0.94994457, 0.00935873],
[50.62643906, -53.44662925, -1.00975901, 0.98799424, 0.00862883],
[54.98075164, -57.84171932, -1.00899169, 1.02413970, 0.00798996],
[59.89973979, -62.80308286, -1.00825846, 1.06188807, 0.00737475],
[65.48403363, -68.43148354, -1.00755896, 1.10137286, 0.00678335],
[71.18060670, -74.16935766, -1.00695797, 1.13851829, 0.00627156],
[77.62704316, -80.65874432, -1.00638379, 1.17731686, 0.00577926],
[84.95932590, -88.03576701, -1.00583615, 1.21790910, 0.00530656],
[92.35561368, -95.47345454, -1.00537146, 1.25562530, 0.00490290],
[100.72936254, -103.89026974, -1.00492737, 1.29501113, 0.00451487],
[108.99738249, -112.19745198, -1.00455556, 1.33095545, 0.00418819],
[118.29522814, -121.53595325, -1.00419929, 1.36839540, 0.00387356],
[128.79962584, -132.08262114, -1.00385843, 1.40745217, 0.00357100],
[140.72675806, -144.05377481, -1.00353287, 1.44826246, 0.00328056],
[152.28288570, -155.64914588, -1.00326592, 1.48475522, 0.00304128],
[165.28451780, -168.69153159, -1.00301004, 1.52275978, 0.00281094],
[179.98119764, -183.43059885, -1.00276518, 1.56239974, 0.00258957],
[196.67856131, -200.17212329, -1.00253126, 1.60381462, 0.00237719],
[212.39108902, -215.92291985, -1.00234463, 1.63979641, 0.00220707],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 5.05e-07
#    Est Max FD Error (based on N/2 run) to R=100.04: 3.6e-7
#    Est Max MOC error for R>100.04: 5.69e-07 
#    Max interpolation error with cubic splines: 5.0459e-07
#    Est Max overall error after cubic interpolation: 1.4e-6
#    

# Run:P8000_G1x8
'P1.8' => [
[-11.72907336, 12.00017795, -0.99999210, -11.73094871, 0.99906145],
[-9.95528800, 10.22643255, -0.99995152, -9.95984656, 0.99771565],
[-8.76516145, 9.03641707, -0.99984025, -8.77344161, 0.99584351],
[-7.82937016, 8.10087291, -0.99959353, -7.84262144, 0.99333346],
[-7.07480056, 7.34676096, -0.99913676, -7.09417907, 0.99022619],
[-6.43771260, 6.71044123, -0.99837138, -6.46444617, 0.98647890],
[-5.88529891, 6.15922641, -0.99718043, -5.92066363, 0.98206102],
[-5.39600036, 5.67170283, -0.99542489, -5.44134420, 0.97693133],
[-4.95450684, 5.23273914, -0.99293854, -5.01128762, 0.97103180],
[-4.54988256, 4.83161853, -0.98952544, -4.61970027, 0.96429218],
[-4.17337889, 4.45987023, -0.98494847, -4.25803899, 0.95661656],
[-3.81749699, 4.11036206, -0.97891203, -3.91910605, 0.94787081],
[-3.47473804, 3.77611792, -0.97101845, -3.59588478, 0.93784576],
[-3.13666587, 3.44951078, -0.96068706, -3.28074160, 0.92620073],
[-2.79046996, 3.11920926, -0.94690134, -2.96244526, 0.91226857],
[-2.39990598, 2.75314402, -0.92679288, -2.60963379, 0.89393776],
[-2.04731256, 2.43023201, -0.90412173, -2.29771390, 0.87497625],
[-1.72267578, 2.14062461, -0.87949679, -2.01676080, 0.85561381],
[-1.42009257, 1.87834702, -0.85367824, -1.76078026, 0.83615056],
[-1.13773323, 1.64094171, -0.82764708, -1.52736605, 0.81703515],
[-0.86254921, 1.41681723, -0.80112741, -1.30516457, 0.79782628],
[-0.58628895, 1.19923411, -0.77406089, -1.08745320, 0.77829151],
[-0.30074626, 0.98217838, -0.74634398, -0.86809513, 0.75817902],
[-0.00880212, 0.76831572, -0.71896065, -0.64970377, 0.73803513],
[0.29173460, 0.55628537, -0.69235823, -0.43092291, 0.71804832],
[0.60467910, 0.34367571, -0.66680129, -0.20933429, 0.69830339],
[0.93404696, 0.12812586, -0.64252899, 0.01742747, 0.67889422],
[1.28457345, -0.09301167, -0.61973445, 0.25201866, 0.65990495],
[1.66224208, -0.32296095, -0.59856267, 0.49768664, 0.64140213],
[2.07339015, -0.56495166, -0.57918717, 0.75763953, 0.62349820],
[2.48036974, -0.79733995, -0.56334055, 1.00815192, 0.60792313],
[2.89302122, -1.02699500, -0.55017781, 1.25609769, 0.59411137],
[3.31755691, -1.25815178, -0.53920397, 1.50563392, 0.58175502],
[3.75832900, -1.49373188, -0.53008822, 1.75955144, 0.57066973],
[4.21947489, -1.73637525, -0.52256951, 2.02035878, 0.56071728],
[4.70556504, -1.98883277, -0.51642848, 2.29068707, 0.55178391],
[5.21999390, -2.25316848, -0.51149312, 2.57242813, 0.54380209],
[5.76910572, -2.53290722, -0.50758775, 2.86902059, 0.53667959],
[6.35651051, -2.83013234, -0.50457917, 3.18235852, 0.53038186],
[6.99260437, -3.15032624, -0.50231706, 3.51790001, 0.52482064],
[7.68433400, -3.49718925, -0.50068868, 3.87919193, 0.51996923],
[8.44529467, -3.87773462, -0.49957852, 4.27320407, 0.51577168],
[9.29605000, -4.30242656, -0.49888510, 4.71039722, 0.51217273],
[10.26571412, -4.78597437, -0.49852346, 5.20547878, 0.50912857],
[11.38514113, -5.34395593, -0.49842254, 5.77392540, 0.50662529],
[12.67857209, -5.98867450, -0.49851409, 6.42785086, 0.50465716],
[14.29093097, -6.79263797, -0.49874768, 7.24017065, 0.50308632],
[16.44721412, -7.86844488, -0.49908086, 8.32351554, 0.50186150],
[20.49847443, -9.89134941, -0.49953018, 10.35418623, 0.50079806],
],

#    Overall Error Estimate:
#    Energy error to P=0.1 P0: 5.68e-07
#    Est Max FD Error (based on N/2 run) to R=4.00: 1.3e07
#    Est Max MOC error for R>4.00: -3.54e-06 
#    Max interpolation error with cubic splines: 5.0842e-07
#    Est Max overall error after cubic interpolation: 4.1e-6
#    

# Run:S4000_G6x5
'S6.5' => [
[-3.90229432, 12.00033728, -2.99996806, -3.90360033, 0.99803973],
[-3.26297240, 10.08242098, -2.99982620, -3.26638332, 0.99487519],
[-2.85832977, 8.86862994, -2.99941976, -2.86459679, 0.99057182],
[-2.58737371, 8.05600579, -2.99868778, -2.59679703, 0.98580459],
[-2.34978332, 7.34369006, -2.99732213, -2.36326546, 0.97965849],
[-2.14135510, 6.71918025, -2.99500551, -2.15982417, 0.97208672],
[-1.95773199, 6.16953206, -2.99136686, -1.98211415, 0.96308659],
[-1.79829339, 5.69297701, -2.98613995, -1.82933788, 0.95292739],
[-1.65049918, 5.25216254, -2.97856709, -1.68934776, 0.94101730],
[-1.51829821, 4.85902004, -2.96844710, -1.56578630, 0.92784080],
[-1.39043382, 4.48030987, -2.95434073, -1.44810754, 0.91234862],
[-1.27343339, 4.13565081, -2.93630184, -1.34232710, 0.89537663],
[-1.16030639, 3.80474683, -2.91268621, -1.24210082, 0.87603781],
[-1.04901287, 3.48221911, -2.88191239, -1.14580711, 0.85386190],
[-0.93525612, 3.15661114, -2.84098499, -1.05012895, 0.82768084],
[-0.81104340, 2.80715010, -2.78345735, -0.94931345, 0.79484772],
[-0.69347982, 2.48378350, -2.71542619, -0.85789456, 0.75974499],
[-0.58576919, 2.19519287, -2.64138808, -0.77793766, 0.72445694],
[-0.48462580, 1.93195046, -2.56256093, -0.70643679, 0.68908203],
[-0.38878349, 1.69020791, -2.48108601, -0.64205991, 0.65412885],
[-0.29598980, 1.46381736, -2.39780018, -0.58296177, 0.61955517],
[-0.20321378, 1.24531757, -2.31225934, -0.52709226, 0.58486962],
[-0.11028258, 1.03443432, -2.22634534, -0.47433846, 0.55058124],
[-0.01645565, 0.82955516, -2.14123975, -0.42426684, 0.51693907],
[0.07904745, 0.62907233, -2.05787566, -0.37647704, 0.48413365],
[0.17735299, 0.43080346, -1.97671531, -0.33046944, 0.45221197],
[0.27915900, 0.23359789, -1.89845556, -0.28602206, 0.42135042],
[0.38552605, 0.03571794, -1.82340079, -0.24280912, 0.39160244],
[0.49726743, -0.16397059, -1.75196856, -0.20066963, 0.36309149],
[0.61549790, -0.36702995, -1.68432273, -0.15937979, 0.33586005],
[0.74106828, -0.57445427, -1.62074366, -0.11886027, 0.31001210],
[0.87520272, -0.78777174, -1.56127519, -0.07895095, 0.28556687],
[1.01945160, -1.00888748, -1.50584580, -0.03945929, 0.26250182],
[1.16937431, -1.23083719, -1.45626454, -0.00170759, 0.24159559],
[1.32495635, -1.45388338, -1.41211898, 0.03437770, 0.22271926],
[1.48817771, -1.68105244, -1.37249043, 0.06929348, 0.20552139],
[1.66097871, -1.91504588, -1.33671177, 0.10341197, 0.18974650],
[1.84061599, -2.15224894, -1.30504953, 0.13619195, 0.17555122],
[2.03321913, -2.40075233, -1.27621559, 0.16870533, 0.16239146],
[2.23812092, -2.65951655, -1.25024976, 0.20071161, 0.15031306],
[2.45563593, -2.92886059, -1.22696365, 0.23217519, 0.13926126],
[2.68868637, -3.21228386, -1.20593011, 0.26341218, 0.12906456],
[2.93855233, -3.51116633, -1.18696243, 0.29445672, 0.11966197],
[3.20671335, -3.82710495, -1.16987334, 0.32535334, 0.11099108],
[3.49524753, -4.16236445, -1.15446288, 0.35619300, 0.10298075],
[3.80622519, -4.51915101, -1.14056311, 0.38703653, 0.09557367],
[4.14251184, -4.90053151, -1.12800081, 0.41799235, 0.08870672],
[4.50656376, -5.30905353, -1.11664785, 0.44909706, 0.08233838],
[4.90143149, -5.74789243, -1.10637663, 0.48041204, 0.07642500],
[5.33115878, -6.22126996, -1.09705989, 0.51204024, 0.07091991],
[5.79678053, -6.73006680, -1.08864662, 0.54384367, 0.06581943],
[6.30572805, -7.28212124, -1.08098286, 0.57609801, 0.06105487],
[6.86122378, -7.88060810, -1.07401374, 0.60874790, 0.05661455],
[7.46868532, -8.53104280, -1.06766361, 0.64184734, 0.05247157],
[8.13332530, -9.23867215, -1.06187348, 0.67540358, 0.04860708],
[8.86115242, -10.00955074, -1.05658804, 0.70943410, 0.04500213],
[9.65857235, -10.85011168, -1.05175954, 0.74394325, 0.04164040],
[10.53238852, -11.76716979, -1.04734638, 0.77892349, 0.03850752],
[11.49060288, -12.76875941, -1.04330897, 0.81438530, 0.03558836],
[12.54081644, -13.86245923, -1.03961641, 0.85029533, 0.03287207],
[13.69202970, -15.05727625, -1.03623797, 0.88664423, 0.03034619],
[14.95344286, -16.36239273, -1.03314756, 0.92340279, 0.02800013],
[16.33605602, -17.78882215, -1.03031923, 0.96056865, 0.02582200],
[17.85046918, -19.34713342, -1.02773216, 0.99810215, 0.02380261],
[19.44271231, -20.98166564, -1.02545248, 1.03452853, 0.02200082],
[21.22607384, -22.80845077, -1.02330944, 1.07219374, 0.02028722],
[23.11723728, -24.74183997, -1.02140111, 1.10906209, 0.01874444],
[25.26759031, -26.93619714, -1.01958053, 1.14772527, 0.01725720],
[27.53886627, -29.25006161, -1.01796779, 1.18536902, 0.01592671],
[29.91247768, -31.66458630, -1.01654498, 1.22173036, 0.01474233],
[32.59956760, -34.39425425, -1.01518488, 1.25976831, 0.01360050],
[35.65702644, -37.49609676, -1.01388697, 1.29962368, 0.01250174],
[38.84483711, -40.72631476, -1.01275137, 1.33789028, 0.01153278],
[42.46967011, -44.39534598, -1.01166719, 1.37795743, 0.01060079],
[46.21339232, -48.18093526, -1.01072583, 1.41607725, 0.00978589],
[50.46290865, -52.47406406, -1.00982634, 1.45595023, 0.00900208],
[54.79745236, -56.84947705, -1.00905234, 1.49346048, 0.00832345],
[59.70226653, -61.79682341, -1.00831180, 1.53263972, 0.00767037],
[65.28091918, -67.41981319, -1.00760449, 1.57362989, 0.00704303],
[70.90830738, -73.08825732, -1.00700350, 1.61171565, 0.00650711],
[77.27965198, -79.50233930, -1.00642845, 1.65148742, 0.00599177],
[84.53094555, -86.79819665, -1.00587919, 1.69309002, 0.00549710],
[91.72911449, -94.03698602, -1.00541962, 1.73112300, 0.00508130],
[99.86820607, -102.21834703, -1.00497958, 1.77080823, 0.00468147],
[109.11840131, -111.51260169, -1.00455897, 1.81228666, 0.00429768],
[118.08962749, -120.52313856, -1.00421382, 1.84938754, 0.00398152],
[128.19321182, -130.66757977, -1.00388281, 1.88803686, 0.00367722],
[139.62681577, -142.14371487, -1.00356587, 1.92836264, 0.00338482],
[152.63307916, -155.19432846, -1.00326292, 1.97050960, 0.00310435],
[164.88763076, -167.48734608, -1.00302112, 2.00714117, 0.00287976],
[178.65444879, -181.29411459, -1.00278895, 2.04525862, 0.00266349],
[194.19165295, -196.87287432, -1.00256639, 2.08498253, 0.00245557],
[211.81476347, -214.53928147, -1.00235338, 2.12644885, 0.00225600],
],
    };
}

1;
