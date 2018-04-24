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

use Blast::IPS::ShockTablesIndex;
use Blast::IPS::ShockTables;
use Blast::IPS::AlphaTable;

my $rtables_info = $Blast::IPS::ShockTablesIndex::rtables_info;
my $rtables      = $Blast::IPS::ShockTables::rtables;
my $ralpha_table = $Blast::IPS::AlphaTable::ralpha_table;
my $rgamma_table;

INIT {

    # Make an index table for searching gamma...

    # Given the hash of info with data of the form:
    #    table_name   => [ symmetry, gamma,  9.8e-7, 4000, ... ],
    # Invert to make a lookup table of sorted gamma values of the form:
    #     $rgamma->[symmetry]->[gamma, table name]
    # This is needed for searching for a specific gamma value.
    my $rtmp = [];
    foreach my $key ( keys %{$rtables_info} ) {
        my $item = $rtables_info->{$key};
        my ( $symmetry, $gamma ) = @{$item};
        if ( $symmetry != 0 && $symmetry != 1 && $symmetry != 2 ) {
            croak "Unexpected symmetry $symmetry in table $key of Blast::Table";
        }
        push @{ $rtmp->[$symmetry] }, [ $gamma, $key ];
    }

    foreach my $sym ( 0, 1, 2 ) {
        my @sorted = sort { $a->[0] <=> $b->[0] } @{ $rtmp->[$sym] };

        # We require a table of unique gamma values.
        # If there are multiple tables for a given gamma, use the
        # table with least error.
        my @unique;
        my ( $gamma_last, $key_last );
        foreach my $item (@sorted) {
            my ( $gamma, $key ) = @{$item};
            if ( !$gamma_last || $gamma != $gamma_last ) {
                push @unique, $item;
                $gamma_last = $gamma;
                $key_last   = $key;
                next;
            }
            my $err_last = $rtables_info->{$key_last}->[2];
            my $err      = $rtables_info->{$key}->[2];
            next if ( $err_last < $err );
            $unique[-1] = $item;
        }

        $rgamma_table->[$sym] = \@unique;
    }
}

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
# Maybe allow lookup on slope dYdX or dZdX, with quadratic interpolation

# Work is needed to handle extrapolation beyond the ranges of the tables for
# some spherical symmetry variables.

# An easy way to test the extrapolation is to have the driver get a
# builtin table, truncate at both ends, reinstall it and look at the errors at
# both missing ends.

sub new {
    my ( $class, @args ) = @_;
    my $self = {};
    bless $self, $class;
    $self->_setup(@args);
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

    # The following four quantities can be specified:
    my $rtable     = $rinput_hash->{rtable};
    my $table_name = $rinput_hash->{table_name};
    my $gamma      = $rinput_hash->{gamma};
    my $symmetry   = $rinput_hash->{symmetry};

    # Allow input of the older 'ASYM' keyword for the symmetry
    if ( !defined($symmetry) ) { $symmetry = $rinput_hash->{ASYM}; }

    # allow alphanumeric abbreviations for symmetry
    # i.e. 'S', 'sph', 'spherical' etc
    if ( defined($symmetry) ) {
        if    ( $symmetry =~ /^S/i ) { $symmetry = 2 }
        elsif ( $symmetry =~ /^C/i ) { $symmetry = 1 }
        elsif ( $symmetry =~ /^P/i ) { $symmetry = 0 }
    }

    # Input Option 1: a user-supplied table is given
    # symmetry and gamma must be given too
    # table name is optional
    if ( defined($rtable) ) {

        if ( !defined($symmetry) || !defined($gamma) ) {
            carp(<<EOM);
You must give a value for symmetry and gamma when you supply a complete table
EOM
            return;
        }
        if ( !defined($table_name) ) {
            $table_name = _make_table_name( $symmetry, $gamma );
            $table_name .= "_custom";
        }
    }

    # Input Option 2: a builtin table name is given
    # symmetry and gamma must match if given; it is safest not to give them. 
    elsif ( defined($table_name) ) {

        $rtable = get_builtin_table($table_name);
        my $item = $rtables_info->{$table_name};
        if ( defined($item) ) {
	    my $old_symmetry=$symmetry;
	    my $old_gamma=$gamma;
            ( $symmetry, $gamma, my $err_est ) = @{$item};
            if ( defined($old_symmetry) && $old_symmetry != $symmetry ) {
                carp(<<EOM);
symmetry of builtin table name: '$table_name' is '$symmetry' which differs from
your value '$symmetry'.
You should either specify a table name or else symmetry and gamma.
EOM
                return;
            }
            if ( defined($old_gamma) && $old_gamma != $gamma ) {
                carp(<<EOM);
gamma of builtin table name: '$table_name' is '$gamma' which differs from
your value '$gamma'. 
You should either specify a table name or else symmetry and gamma.
EOM
                return;
            }
        }
        else {
            carp(<<EOM);
Not a builtin table name: '$table_name'

To create a table for any symmetry and gamma you can use:
   Blast::IPS->new(symmetry='S', gamma=>1.23); 
where symmetry can be 'P', 'C' or 'S'
and gamma can be any value between 1.1 and 6.5

EOM
            return;
        }
    }

    # Input Option 3: gamma and symmetry are given. This is the normal method.
    else {

        # allow default of spherical symmetry, gamma=1.4
        if ( !defined($gamma) )    { $gamma    = 1.4 }
        if ( !defined($symmetry) ) { $symmetry = 2 }

        # Search for this gamma value
        my $result = _gamma_lookup( $symmetry, $gamma );
	my $old_table_name=$table_name;
        $table_name = $result->{table_name};

        # Get builtin table if there is one
        if ( defined($table_name) ) {
            $rtable = get_builtin_table($table_name);
        }

        # Otherwise, make an interpolated table
        else {

            my $igam1 = $result->{igam1};
            my $igam2 = $result->{igam2};
            my $igam3 = $result->{igam3};
            my $igam4 = $result->{igam4};
            $rtable = _make_intermediate_gamma_table( $symmetry, $gamma,
                $igam1, $igam2, $igam3, $igam4 );
            if ( !$table_name ) {
		$table_name = _make_table_name($symmetry, $gamma);
            }
        }
    }

    # Now check the results
    my $error = "";
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

sub _make_table_name {

    # make a table name by combining the symmetry and gamma;
    # i.e. 'S1.37' would be spherical symmetry with gamma=1.37
    my ( $symmetry, $gamma ) = @_;
    my $table_name;
    my $gamma_pr = sprintf "%.4g", $gamma;
    my $char = $symmetry == 0 ? 'P' : $symmetry == 1 ? 'C' : 'S';
    $table_name = $char . $gamma_pr;
    return $table_name;
}

sub _gamma_lookup {
    my ( $sym, $gamma, $eps ) = @_;

    # Given: a 1d symmetry (0, 1, or 2) and an ideal gas gamma
    # Return: the names of the builtin tables with bounding gammas

    # uses $rgamma_table, the global table of gamma values for the builtin
    # tables

    # There are three possible results:

    # 1. returns nothing if gamma is out of bounds or there is an error

    # 2. If the gamma is within a small tolerance of one of these gammas
    # then returns reference to a hash with
    #     table_name => $table_name

    # 3. If the gamma differs from the table entries but is within
    # bounds of the table, returns four surrounding table names
    # for use in cubic interpolation. These have hash keys
    # table_name1, table_name2, table_name3 and table_name4

    return if ( $sym != 0 && $sym != 1 && $sym != 2 );
    $eps = 1.e-5 unless ($eps);

    my $return_hash = {};

    # lookup this gamma in the table of gamma values for this symmetry
    my $rtab  = $rgamma_table->[$sym];
    my $ntab  = @{$rtab};
    my $icol  = 0;
    my $rhash = {
        _jl     => undef,
        _ju     => undef,
        _rtable => $rtab,
    };
    my ( $j2, $j3 ) = Blast::IPS::_locate_2d( $rhash, $gamma, $icol );
    my ( $gamma_min, $key_min ) = @{ $rtab->[0] };
    my ( $gamma_max, $key_max ) = @{ $rtab->[-1] };

    # Check if out of bounds
    if ( $j2 < 0 ) {
        return if ( $gamma + $eps < $gamma_min );
        $return_hash->{'table_name'} = $key_min;
        $return_hash->{'igam'}       = 0;
        return $return_hash;
    }
    if ( $j3 >= $ntab ) {
        return if ( $gamma - $eps > $gamma_max );
        $return_hash->{'table_name'} = $key_max;
        $return_hash->{'igam'}       = $ntab - 1;
        return $return_hash;
    }

    # Check for an exact match:
    my ( $gamma2, $key2 ) = @{ $rtab->[$j2] };
    my ( $gamma3, $key3 ) = @{ $rtab->[$j3] };
    if ( abs( $gamma - $gamma2 ) < $eps ) {
        $return_hash->{'table_name'} = $key2;
        $return_hash->{'igam'}       = $j2;
        return $return_hash;
    }
    if ( abs( $gamma - $gamma3 ) < $eps ) {
        $return_hash->{'table_name'} = $key3;
        $return_hash->{'igam'}       = $j3;
        return $return_hash;
    }

    # Not an exact match: return indexes of four bounding tables in preparation
    # for cubic interpolation.
    my $j1 = $j2 - 1;
    my $j4 = $j3 + 1;
    if ( $j1 < 0 )      { $j1 = $j2 }
    if ( $j4 >= $ntab ) { $j4 = $j3 }

    $return_hash = {
        igam1 => $j1,
        igam2 => $j2,
        igam3 => $j3,
        igam4 => $j4,
    };
    return ($return_hash);
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
    return ( $rtables->{$table_name} );
}

sub get_rgamma_table {
    return $rgamma_table;
}

sub get_info {
    my ($self) = @_;
    if ( !ref($self) ) {
        my $what = ref $self;
        print STDERR "$what\n";
        croak("get_info not called by a Blast::IPS object");
    }

    # Return some information about this case
    my $rinfo      = {};
    my $table_name = $self->{_table_name};
    my $item       = $rtables_info->{$table_name};
    my (
        $symmetry,     $gamma,    $Max_Error, $N,
        $Energy_Error, $R_FD_MOC, $FD_Error,  $MOC_Error,
        $Interp_Error, $rs2,      $zs2
    ) = @{$item};
    $rinfo->{table_name}      = $table_name;
    $rinfo->{alpha}           = $self->{_alpha};
    $rinfo->{symmetry}        = $self->{_symmetry};
    $rinfo->{gamma}           = $self->{_gamma};
    $rinfo->{Max_Error} = $Max_Error;
    $rinfo->{Xs2} = $rs2?log($rs2):undef;
    $rinfo->{Zs2} = $zs2?log($zs2):undef;
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

        #my ( $symmetry, $gamma, $Max_Error, $N ) = @{$item};
        my (
            $symmetry,     $gamma,    $Max_Error, $N,
            $Energy_Error, $R_FD_MOC, $FD_Error,  $MOC_Error,
            $Interp_Error, $rs2,      $zs2
        ) = @{$item};
        $rtable_index->{$key} = {
            symmetry  => $symmetry,
            gamma     => $gamma,
            Max_Error => $Max_Error,
            N         => $N,
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
        my $letter = $symmetry eq 0 ? 'P' : $symmetry eq 1 ? 'C' : 'S';
        my $str = $letter . $gamma;

        #print STDERR "Checking if table $key similar to $str\n";
        if ( $key !~ /^$str/i ) {
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

    # TODO: Generalize to allow any number of lagrange interpolation points

    return if ( $sym != 0 && $sym != 1 && $sym != 2 );

    my $rtab = $ralpha_table->[$sym];
    my $ntab = @{$rtab};
    my ( $jl, $ju );
    my $icol = 0;

    # A small tolerance to avoid interpolations
    my $eps  = 1.e-6;

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

    my $jbase = $jl - 1;
    if ( $jl <= 0 ) {
        $jbase = $jl;
    }
    elsif ( $ju >= $ntab - 1 ) {
        $jbase = $ju - 3;
    }

    my ( $rx, $ry );
    for ( my $jj = $jbase ; $jj <= $jbase + 3 ; $jj += 1 ) {
        my ( $xx, $yy ) = @{ $rtab->[$jj] };
        push @{$rx}, $xx;
        push @{$ry}, $yy;
    }

    my ( $alpha, $dalpha ) = polint( $gamma, $rx, $ry );

    return ($alpha);
}

sub _make_intermediate_gamma_table {

    # Given: a symmetry, gamma, and the indexes of four bounding builtin tables
    #     such that the gamma of interest is between numbers 2 and 3
    #     and such that at least three are different tables
    # Return: a table of X,Y,dYdX,Z,dZdX for an intermediate gamma
    #     by using cubic interpolation between the four builtin tables

    # TODO: Generalize to allow arbitrary number of interpolation points

    my ( $symmetry, $gamma_x, $igam1, $igam2, $igam3, $igam4 ) = @_;
    my $rtable_new;
    my $alpha_x = alpha_interpolate( $symmetry, $gamma_x );

    # Evaluate the interpolation fraction
    my ( $gamma_1, $table_name_1 ) = @{ $rgamma_table->[$symmetry]->[$igam1] };
    my ( $gamma_2, $table_name_2 ) = @{ $rgamma_table->[$symmetry]->[$igam2] };
    my ( $gamma_3, $table_name_3 ) = @{ $rgamma_table->[$symmetry]->[$igam3] };
    my ( $gamma_4, $table_name_4 ) = @{ $rgamma_table->[$symmetry]->[$igam4] };

    my $alpha_1 = alpha_interpolate( $symmetry, $gamma_1 );
    my $alpha_2 = alpha_interpolate( $symmetry, $gamma_2 );
    my $alpha_3 = alpha_interpolate( $symmetry, $gamma_3 );
    my $alpha_4 = alpha_interpolate( $symmetry, $gamma_4 );

    my $rtable_1 = get_builtin_table($table_name_1);
    my $rtable_2 = get_builtin_table($table_name_2);
    my $rtable_3 = get_builtin_table($table_name_3);
    my $rtable_4 = get_builtin_table($table_name_4);

    my $A_1 = -log( ( $gamma_1 + 1 ) * $alpha_1 );
    my $A_2 = -log( ( $gamma_2 + 1 ) * $alpha_2 );
    my $A_x = -log( ( $gamma_x + 1 ) * $alpha_x );
    my $A_3 = -log( ( $gamma_3 + 1 ) * $alpha_3 );
    my $A_4 = -log( ( $gamma_4 + 1 ) * $alpha_4 );

    # This is the fraction that would be used for linear interpolation
    # of X
    my $ff = ( $A_x - $A_2 ) / ( $A_3 - $A_2 );

    # Use the Y values of the closest table so that we converge to it. But clip
    # the Y to the lower bounds of the other table to avoid extrapolation at
    # long range where the theory is not great.
    my $rtable_closest;
    if ( $ff <= 0.5 ) {
        $rtable_closest = $rtable_2;
    }
    else {
        $rtable_closest = $rtable_3;
    }

    my $lookup = sub {

        # This is a specialized version of lookup for gamma interpolation
        # for this routine which avoids extrapolation

        # Given:
        #     $Y = a value of 'Y' = ln(overpressure) to lookup
        #     $rtab = a builtin blast table

        # Returns undef if out of bounds of table; otheerwise
        # Returns
        # 	[X, Y, dY/dX, Z]

        my ( $Y, $rtab ) = @_;
        my $icol   = 1;
        my $interp = 0;    # use cubic interpolation

        my $ntab  = @{$rtab};
        my $rhash = {
            _rtable => $rtab,
            _jl     => undef,
            _ju     => undef,
        };

        # Locate this point in the table
        my ( $il, $iu ) = _locate_2d( $rhash, $Y, $icol );

        # Handle case before start of the table
        if ( $il < 0 ) {
            return;
        }

        # Handle case beyond end of the table
        elsif ( $iu >= $ntab ) {
            return;
        }

        # otherwise interpolate within table
        else {
            return _interpolate_rows( $Y, $icol, $rtab->[$il],
                $rtab->[$iu], $interp );
        }
    };

    # Loop over all Y values in the closest table to make the interpolated
    # table
    my $rtab_x = [];
    foreach my $item_ref ( @{$rtable_closest} ) {
        my $Y      = $item_ref->[1];
        my $item_1 = $lookup->( $Y, $rtable_1 );
        my $item_2 = $lookup->( $Y, $rtable_2 );
        my $item_3 = $lookup->( $Y, $rtable_3 );
        my $item_4 = $lookup->( $Y, $rtable_4 );

        # Can only interpolate if all 4 tables are defined
        next unless ( $item_1 && $item_2 && $item_3 && $item_4 );

        my ( $X_1, $YY_1, $dYdX_1, $Z_1, $dZdX_1 ) = @{$item_1};
        my ( $X_2, $YY_2, $dYdX_2, $Z_2, $dZdX_2 ) = @{$item_2};
        my ( $X_3, $YY_3, $dYdX_3, $Z_3, $dZdX_3 ) = @{$item_3};
        my ( $X_4, $YY_4, $dYdX_4, $Z_4, $dZdX_4 ) = @{$item_4};

        my $lagrange_interp = sub {

            my ( $Q_1, $Q_2, $Q_3, $Q_4 ) = @_;
            my $rx = [ $A_1, $A_2, $A_3, $A_4 ];
            my $ry = [ $Q_1, $Q_2, $Q_3, $Q_4 ];
            my ($Q_x, $dQ) = polint( $A_x, $rx, $ry );
            return $Q_x;
        };

        my $X_x    = $lagrange_interp->( $X_1,    $X_2,    $X_3,    $X_4 );
        my $Z_x    = $lagrange_interp->( $Z_1,    $Z_2,    $Z_3,    $Z_4 );
        my $dYdX_x = $lagrange_interp->( $dYdX_1, $dYdX_2, $dYdX_3, $dYdX_4 );
        my $dZdX_x = $lagrange_interp->( $dZdX_1, $dZdX_2, $dZdX_3, $dZdX_4 );

        push @{$rtab_x}, [ $X_x, $Y, $dYdX_x, $Z_x, $dZdX_x ];
    }
    return $rtab_x;
}

sub polint {

    #  Slightly modified versions of the "polint" routine from
    #  Press, William H., Brian P. Flannery, Saul A. Teukolsky and
    #  William T. Vetterling, 1986, "Numerical Recipes: The Art of
    #  Scientific Computing" (Fortran), Cambrigde University Press,
    #  pp. 80-82.

    # Given:
    # $xx = an x location where y is required
    # ($rx, $ry) = arrays of ($x,$y) lagrange interpolation points 

    # Return:
    #  $yy = the interpolated value at $xx
    #  $dy = the estimated error

    # Example call:
    # my ( $yy, $dy ) = polint( $xx, $rx, $ry );

    my ( $xx, $rx, $ry ) = @_;

    my $n = @{$rx};

    # Return values
    my ( $yy, $dy );

    #..find the index ns of the closest table entry; initialize the c and d
    # tables
    my ( @c, @d );
    my $ns  = 0;
    my $dif = abs( $xx - $rx->[0] );
    for ( my $i = 0 ; $i < $n ; ++$i ) {
        my $dift = abs( $xx - $rx->[$i] );
        if ( $dift < $dif ) {
            $ns  = $i;
            $dif = $dift;
        }
        $c[$i] = $ry->[$i];
        $d[$i] = $ry->[$i];
    }

    #..first guess for y
    $yy = $ry->[$ns];

    #..for each column of the table, loop over the c's and d's and update them
    $ns = $ns - 1;
    for ( my $m = 0 ; $m < $n - 1 ; ++$m ) {
        for ( my $i = 0 ; $i < $n - $m - 1 ; ++$i ) {
            my $ho  = $rx->[$i] - $xx;
            my $hp  = $rx->[ $i + $m + 1 ] - $xx;
            my $w   = $c[ $i + 1 ] - $d[$i];
            my $den = $ho - $hp;
            if ( $den == 0.0 ) {
                print STDERR "polint: 2 rx entries are the same \n";
                return;
            }
            $den   = $w / $den;
            $d[$i] = $hp * $den;
            $c[$i] = $ho * $den;
        }

	# after each column is completed, decide which correction c or d, to
	# add to the accumulating value of y, that is, which path to take in
	# the table by forking up or down. ns is updated as we go to keep track
	# of where we are. the last dy added is the error indicator.
        if ( 2 * ( $ns + 1 ) < $n - $m - 1 ) {
            $dy = $c[ $ns + 1 ];
        }
        else {
            $dy = $d[$ns];
            $ns = $ns - 1;
        }
        $yy += $dy;
    }
    return wantarray ? ( $yy, $dy ) : $yy;
}

1;
