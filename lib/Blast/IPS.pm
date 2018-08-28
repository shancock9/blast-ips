package Blast::IPS;

# This program can evaluate the shock overpressure for a one dimensional point,
# line or sheet source explosion in an ideal homogeneous atmosphere.  This
# problem does not have an analytic solution, so we interpolate tables of
# values.

# The builtin tables cover the three one-dimensional symmetries (plane,
# cylindrical, spherical) and numerous values of the ideal gas gamma (between
# 1.1 and 6.5.  They were prepared with calculations using the finite
# difference method and the method of characteristics.  The estimated relative
# accuracy of interpolated shock overpressures depends on the table and
# whether interpolation is made in gamma, but it is on the order of 1.e-6
# in most cases and well below 1.e-5 in all cases.

# All calculations are done in dimensionless units. The calling program must
# handle conversion to and from physical units.

# TODO list:

# Update for specifying z/x:
# _update_toa to also set
#   [I_ZmX] = $Z-$X
#   [I_dZmXdX] = $dZdX-1
# allow 'Z-X' queries

# - Needed to check the extrapolation method beyond the ranges of the tables for
# some spherical symmetry variables.

# - An easy way to test the extrapolation is to have the driver get a
# builtin table, truncate at both ends, reinstall it and look at the errors at
# both missing ends.

# - allow lookup on slope dYdX or dZdX, with quadratic interpolation
#   and handle case where they exceed limiting values

# - Need clear 'out of bounds' signal

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

use Blast::IPS::AlphaTable;
use Blast::IPS::BlastInfo qw(get_blast_info);
use Blast::IPS::ImpulseTables;
use Blast::IPS::PzeroFit;
use Blast::IPS::PzeroTail;
use Blast::IPS::ShockTables;
use Blast::IPS::ShockTablesIndex;
use Blast::IPS::EnergyTables;
use Blast::IPS::TailShockTables;

my $rshock_tables_info = $Blast::IPS::ShockTablesIndex::rshock_tables_info;
my $rshock_tables      = $Blast::IPS::ShockTables::rshock_tables;
my $ralpha_table       = $Blast::IPS::AlphaTable::ralpha_table;
my $rpzero_fit         = $Blast::IPS::PzeroFit::rpzero_fit;
my $rpzero_tail        = $Blast::IPS::PzeroTail::rpzero_tail;
my $rimpulse_tables    = $Blast::IPS::ImpulseTables::rimpulse_tables;
my $renergy_tables     = $Blast::IPS::EnergyTables::renergy_tables;
my $rtail_shock_tables = $Blast::IPS::TailShockTables::rtail_shock_tables;
my $rgamma_table;

BEGIN {

    # Shock Table variable layout
    my $i = 0;
    use constant {
        I_X      => $i++,
        I_Y      => $i++,
        I_dYdX   => $i++,
        I_Z      => $i++,
        I_dZdX   => $i++,
        I_T      => $i++,
        I_dTdX   => $i++,
        I_ZmX    => $i++,
        I_dZmXdX => $i++,
        I_E1     => $i++,
        I_dE1dX  => $i++,
        I_Er     => $i++,
        I_dErdX  => $i++,
    };
}

INIT {

    # Make an index table for searching gamma...

    # Given the hash of info with data of the form:
    #    table_name   => [ symmetry, gamma,  9.8e-7, 4000, ... ],
    # Invert to make a lookup table of sorted gamma values of the form:
    #     $rgamma->[symmetry]->[gamma, table name]
    # This is needed for searching for a specific gamma value.
    my $rtmp = [];
    foreach my $key ( keys %{$rshock_tables_info} ) {
        my $item = $rshock_tables_info->{$key};
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
            my $err_last = $rshock_tables_info->{$key_last}->[2];
            my $err      = $rshock_tables_info->{$key}->[2];
            next if ( $err_last < $err );
            $unique[-1] = $item;
        }

        $rgamma_table->[$sym] = \@unique;
    }

    sub _merge_energy_tables {

        # Merge the energy tables into the shock tables. The ShockTables and
        # EnergyTables have the same X values in column 1 but they have been
        # stored separately for flexibility.  The disadvantage is that we have
        # to merge them and at the same time verify that the X values are still
        # identical.
        my ($table_name);
        my @errors;
        my $table_error = sub {
            my ($msg) = @_;
            print STDERR "Error merging table $table_name: $msg\n";
            push @errors, $table_name;
        };
        foreach $table_name ( sort keys %{$rshock_tables} ) {
            my $renergy_table = $renergy_tables->{$table_name};
            my $rshock_table  = $rshock_tables->{$table_name};
            my $num0          = @{$rshock_table};
            my $num1          = @{$renergy_table};
            my $rnew_table;

            if ( $num0 != $num1 ) {
                $table_error->("number of rows differ: $num0 != $num1");
                next;
            }
            for ( my $i = 0 ; $i < $num0 ; $i++ ) {
                my $X0 = $rshock_table->[$i]->[0];
                my $X1 = $renergy_table->[$i]->[0];
                if ( $X0 != $X1 ) {
                    $table_error->("X values differ for row $i: $X0 != $X1");
                    next;
                }

                # Combine variables in shock table and energy table,
                # leaving spots for T and dTdX, Z-X and dZ-X/dX
                my @vars = @{ $renergy_table->[$i] };
                shift @vars;    # remove leading X
                push @{$rnew_table},
                  [ @{ $rshock_table->[$i] }, 0, 0, 0, 0, @vars ];
            }

            # Success, install the new table
            $rshock_tables->{$table_name} = $rnew_table;
        }
        if (@errors) {
            print STDERR
              "Tables differ. Run blast_energy_integral_primary.pl\n";
            local $" = ')(';
            print STDERR "(@errors)\n";
            croak "Table Errors; run terminated\n";
        }
        return;
    }
    _merge_energy_tables();

    # TODO: merge_impulse_tables
}

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
    my ( $rimpulse_table, $rtail_shock_table, $loc_gamma_table, $rigam_4,
        $rigam_6 );

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

        ( $rtable, $rimpulse_table, $rtail_shock_table ) =
          get_builtin_tables($table_name);
        my $item = $rshock_tables_info->{$table_name};
        if ( defined($item) ) {
            my $old_symmetry = $symmetry;
            my $old_gamma    = $gamma;
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
        my $old_table_name = $table_name;
        $table_name      = $result->{table_name};
        $loc_gamma_table = $result->{loc_gamma_table};

        # Get builtin table if there is one
        if ( defined($table_name) ) {
            ( $rtable, $rimpulse_table, $rtail_shock_table ) =
              get_builtin_tables($table_name);
        }

        # Otherwise, make an interpolated table
        else {
            $rigam_6 = $result->{rigam_6};
            $rigam_4 = $result->{rigam_4};
            if ( defined($rigam_6) && defined($rigam_4) ) {
                $rtable = _make_intermediate_gamma_table( $symmetry, $gamma,
                    $rigam_6, $rigam_4 );

                # FIXME: make interpolated impulse table too
                if ( !$table_name ) {
                    $table_name = _make_table_name( $symmetry, $gamma );
                }
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
        _update_toa_table($rtable);
        $num_table = @{$rtable};
    }

    my $rztables = _make_z_tables($rimpulse_table);

    $self->{_rtable}                  = $rtable;
    $self->{_rimpulse_table}          = $rimpulse_table;
    $self->{_rztables}                = $rztables;
    $self->{_rtail_shock_table_table} = $rtail_shock_table;
    $self->{_table_name}              = $table_name;
    $self->{_loc_gamma_table}         = $loc_gamma_table;
    $self->{_rigam_4}                 = $rigam_4;
    $self->{_rigam_6}                 = $rigam_6;
    $self->{_gamma}                   = $gamma;
    $self->{_symmetry}                = $symmetry;
    $self->{_jl}                      = -1;
    $self->{_ju}                      = $num_table;
    $self->{_jl2}                     = -1;
    $self->{_ju2}                     = undef;
    $self->{_error}                   = $error;

    if ($error) { carp "$error\n" }
    else {
        $self->_end_model_setup();
    }
    return;
}

sub _end_model_setup {
    my ($self) = @_;

    # Add additional quantities of interest to $self
    my $rtable   = $self->{_rtable};
    my $gamma    = $self->{_gamma};
    my $symmetry = $self->{_symmetry};

    # Set parameters for evaluation beyond the ends of the table
    # Set asymptotic wave parameters, given a table and gamma

    # Near region: we find the value of alpha at the first table point
    # based on R^(symmetry+1)*u**2=constant
    my ( $X_near, $Y_near, $dYdX_near, $Z_near, $dZdX_near ) =
      @{ $rtable->[0] }[ I_X, I_Y, I_dYdX, I_Z, I_dZdX ];
    my $lambda = exp($X_near);
    my $Prat   = exp($Y_near) + 1;
    my $q      = 2 * $gamma / ( ( $gamma + 1 ) * $Prat + ( $gamma - 1 ) );
    my $uovcsq = ( 2 / ( $gamma + 1 ) * ( 1 - $q ) )**2 / $q;
    my $C      = $lambda**( $symmetry + 1 ) * $uovcsq;
    my $alpha =
      ( 4 / ( ( $gamma + 1 ) * ( $symmetry + 3 ) ) )**2 / ( $gamma * $C );
    $self->{_alpha} = $alpha;

    # Distant region
    my ( $X_far, $Y_far, $dYdX_far, $Z_far, $dZdX_far ) =   #@{ $rtable->[-1] };
      @{ $rtable->[-1] }[ I_X, I_Y, I_dYdX, I_Z, I_dZdX ];
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

    $self->_set_global_info();
    return;
}

sub _make_table_name {

    # Make a table name by combining the symmetry and gamma;
    # i.e. make the name 'S1.37' for spherical symmetry with gamma=1.37

    my ( $symmetry, $gamma ) = @_;
    my $table_name;
    my $gamma_pr = sprintf "%.4g", $gamma;
    my $char = $symmetry == 0 ? 'P' : $symmetry == 1 ? 'C' : 'S';
    $table_name = $char . $gamma_pr;
    return $table_name;
}

sub set_interpolation_points {
    my ( $jfloor, $ntab, $NLAG ) = @_;

    # Find the index range for NLAG valid lagrange interpolation points
    # Given:
    #   $jfloor is the first index before the point of interest
    #   $ntab is the number of points in the table
    #   $NLAG is the number of interpolation points desired
    # Return:
    #   a reference to a list of consecutive indexes
    #   the actual number may be less than NLAG for small tables

    return if ( $ntab <= 0 || $NLAG <= 0 );

    # First add points on both sides (will be lopsided if NLAG is odd)
    #my $j_lo = $jfloor - int( $NLAG / 2 ); # ORIGINAL, lopsided
    my $j_lo = $jfloor - int( ( $NLAG - 1 ) / 2 );    # Corrected
    my $j_hi = $j_lo + $NLAG - 1;

    #print STDERR "jfloor=$jfloor, j_lo=$j_lo, j_hi=$j_hi\n";

    # Shift points if too near an edge
    if ( $j_lo < 0 ) {
        my $nshift = 0 - $j_lo;
        $j_lo += $nshift;
        $j_hi += $nshift;
    }
    if ( $j_hi > $ntab - 1 ) {
        my $nshift = $ntab - 1 - $j_hi;
        $j_lo += $nshift;
        $j_hi += $nshift;
    }

    # Be sure points are in bounds
    if ( $j_lo < 0 )         { $j_lo = 0 }
    if ( $j_hi > $ntab - 1 ) { $j_hi = $ntab - 1 }

    #print STDERR "returning j_lo=$j_lo, j_hi=$j_hi\n";
    return [ $j_lo .. $j_hi ];
}

sub _gamma_lookup {
    my ( $symmetry, $gamma ) = @_;

    # Look for a table matching a given symmetry and gamma

    # Given:
    #   $symmetry = a 1d symmetry (0, 1, or 2) and
    #   $gamma = an ideal gas gamma
    # Return:
    #   a return hash, with three possible results as follows;

    # CASE 1. returns nothing if gamma is out of bounds or there is an error

    # CASE 2. If the gamma is within a small tolerance of one of these gammas
    # then returns reference to a hash with
    #     table_name => $table_name

    # CASE 3. If the gamma differs from the table entries but is within
    # bounds of the table, returns index of bounding tables
    # for use in interpolation:
    #    rigam_4 => $rigam_4,
    #    rigam_6 => $rigam_6,
    # where
    # rigam_4 is a list of indexes for 4 point interpolation
    # rigam_6 is a list of indexes for 6 point interpolation

    # for all cases, also returns
    #   loc_gamma_table => has location in the table of gamma values
    # $return_hash->{'loc_gamma_table'} =
    # [index (if exact match), lower index, upper index, $ntab ];

    # uses $rgamma_table, the global table of gamma values for the builtin
    # tables

    return if ( $symmetry != 0 && $symmetry != 1 && $symmetry != 2 );

    # The tolerance for testing gamma is currently fixed. It should be very
    # small because we can construct very accurate interpolated tables.
    my $eps = 1.e-7;

    my $return_hash = {};

    # lookup this gamma in the table of gamma values for this symmetry
    my $rtab  = $rgamma_table->[$symmetry];
    my $ntab  = @{$rtab};
    my $icol  = 0;
    my $rhash = {
        _jl     => undef,
        _ju     => undef,
        _rtable => $rtab,
    };
    my ( $j2, $j3 ) = _locate_2d( $rhash, $gamma, $icol );
    my ( $gamma_min, $key_min ) = @{ $rtab->[0] };
    my ( $gamma_max, $key_max ) = @{ $rtab->[-1] };

    # save the location of this gamma in the gamma table so that we
    # can do further intepolations later

    # Check if out of bounds (CASE 1)
    if ( $j2 < 0 ) {
        return if ( $gamma + $eps < $gamma_min );
        $return_hash->{'table_name'} = $key_min;
        $return_hash->{'loc_gamma_table'} = [ 0, $j2, $j3, $ntab ];
        return $return_hash;
    }
    if ( $j3 >= $ntab ) {
        return if ( $gamma - $eps > $gamma_max );
        $return_hash->{'table_name'} = $key_max;
        $return_hash->{'loc_gamma_table'} = [ $ntab - 1, $j2, $j3, $ntab ];
        return $return_hash;
    }

    # Check for an exact match (CASE 2)
    my ( $gamma2, $key2 ) = @{ $rtab->[$j2] };
    my ( $gamma3, $key3 ) = @{ $rtab->[$j3] };
    if ( abs( $gamma - $gamma2 ) < $eps ) {
        $return_hash->{'table_name'} = $key2;
        $return_hash->{'loc_gamma_table'} = [ $j2, $j2, $j3, $ntab ];
        return $return_hash;
    }
    if ( abs( $gamma - $gamma3 ) < $eps ) {
        $return_hash->{'table_name'} = $key3;
        $return_hash->{'loc_gamma_table'} = [ $j3, $j2, $j3, $ntab ];
        return $return_hash;
    }

    # CASE 3: Not an exact match: return indexes of bounding tables in
    # preparation for interpolation.

    # At the start of the table, in the low gamma range, results are better if
    # we reduce the number of lagrange points to keep the points them from
    # being too lopsided.  Thus, if we are interpolating near the first table
    # point, we will use 4 point interpolation instead of 6 point.
    # At the upper end of the table it is not necessary or beneficial to do
    # this.
    my $six = $j2 == 0 ? 4 : $j2 == 1 ? 5 : 6;

    my $rigam_4 = set_interpolation_points( $j2, $ntab, 4 );
    my $rigam_6 = set_interpolation_points( $j2, $ntab, $six );
    $return_hash->{'rigam_4'}         = $rigam_4;
    $return_hash->{'rigam_6'}         = $rigam_6;
    $return_hash->{'loc_gamma_table'} = [ undef, $j2, $j3, $ntab, ];
    return ($return_hash);
}

sub _check_keys {
    my ( $rtest, $rvalid, $msg, $exact_match ) = @_;

    # Check the keys of a hash for validity:
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
        @missing_keys = sort @missing_keys;
        @unknown_keys = sort @unknown_keys;
        croak(<<EOM);
------------------------------------------------------------------------
Blast::IPS program error detected checking hash keys
Message is: '$msg'
Valid keys are: (@expected_keys)
Keys not seen : (@missing_keys)
Unknown key(s): (@unknown_keys)
------------------------------------------------------------------------
EOM
    }
    return;
}

sub get_builtin_tables {
    my ($table_name) = @_;
    return (
        $rshock_tables->{$table_name},
        $rimpulse_tables->{$table_name},
        $rtail_shock_tables->{$table_name}
    );
}

sub get_rgamma_table {
    return $rgamma_table;
}

sub get_info {
    my ($self) = @_;

    # Return some information about this case
    if ( !ref($self) ) {
        my $what = ref $self;
        print STDERR "$what\n";
        croak("get_info not called by a Blast::IPS object");
    }

    # Return all values in the _rinfo hash
    my $rinfo  = {};
    my $_rinfo = $self->{_rinfo};
    if ( defined($_rinfo) ) {
        foreach my $key ( keys %{$_rinfo} ) {
            $rinfo->{$key} = $_rinfo->{$key};
        }
    }

    # and return a few other key values
    $rinfo->{alpha}    = $self->{_alpha};
    $rinfo->{symmetry} = $self->{_symmetry};
    $rinfo->{gamma}    = $self->{_gamma};
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
    foreach my $key ( sort keys %{$rshock_tables_info} ) {
        my $item = $rshock_tables_info->{$key};

        my (
            $symmetry,     $gamma,    $Max_Error, $N,
            $Energy_Error, $R_FD_MOC, $FD_Error,  $MOC_Error,
            $Interp_Error, $rs2,      $zs2
        ) = @{$item};
        $rtable_index->{$key} = {
            symmetry     => $symmetry,
            gamma        => $gamma,
            Max_Error    => $Max_Error,
            N            => $N,
            r_tail_shock => $rs2,
            z_tail_shock => $zs2,
        };
    }
    return ($rtable_index);
}

sub check_tables {

    # Should be called at program installation to check the tables

    my @info_keys  = keys %{$rshock_tables_info};
    my @table_keys = keys %{$rshock_tables};

    # Check that the keys of the two table hashes are the same
    my @missing_info_keys =
      grep { !exists $rshock_tables_info->{$_} } @table_keys;
    my @missing_table_keys =
      grep { !exists $rshock_tables->{$_} } @info_keys;
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
        my $rtable = $rshock_tables->{$key};

        #print STDERR "Checking table $key\n";

        # Check for a naming error
        my $item = $rshock_tables_info->{$key};
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
    if ( !$is_monotonic->(I_X) > 0 ) {
        $error .= "X is not monotonic increasing\n";
    }

    # Y must be decreasing
    if ( !$is_monotonic->(I_Y) < 0 ) {
        $error .= "Y is not monotonic increasing\n";
    }

    # dYdX must be negative
    foreach my $row ( @{$rtable} ) {
        if ( $row->[I_dYdX] >= 0 ) {
            $error .= "dYdX not everywhere negative\n";
            last;
        }
    }

    # Z must be increasing
    if ( !$is_monotonic->(I_Z) > 0 ) {
        $error .= "Z is not monotonic increasing\n";
    }

    # dZdX must be positive
    foreach my $row ( @{$rtable} ) {
        if ( $row->[I_dZdX] <= 0 ) {
            $error .= "dZdX not everywhere positive\n";
            last;
        }
    }

    return $error;
}

sub _update_toa_row {
    my ($row) = @_;

    # Set the values of T and dTdX, where T=ln(time of arrival)
    # in a single table row of shock values
    # Also set Z-X and its slope.
    # These enable evaluations as function of time and of ln(z/r) or ln(t/r)

    my ( $X, $Y, $dYdX, $Z, $dZdX ) =
      @{$row}[ I_X, I_Y, I_dYdX, I_Z, I_dZdX ];
    my $z    = exp($Z);
    my $dzdX = $z * $dZdX;
    my $x    = exp($X);
    my $T    = log( $x - $z );
    my $dTdX = ( $x - $dzdX ) / ( $x - $z );
    $row->[I_T]      = $T;
    $row->[I_dTdX]   = $dTdX;
    $row->[I_ZmX]    = $Z - $X;
    $row->[I_dZmXdX] = $dZdX - 1;
    return;
}

sub _make_z_tables {
    my ($rimpulse_table) = @_;
    return unless defined($rimpulse_table);

# Create three tables of z values from the impulse table:
# [T_pose, z_pose] = time history of points on the end of the positive phase
# [T_pmin, z_pmin] = time history of points on the minimum overpressure characteristic
# [T_nege, z_nege] = time history of points on the end of the negitive phase

    # These will allow us to define key points on wave profiles

    my ( $rzpose_table, $rzpmin_table, $rznege_table );
    foreach my $rrow ( @{$rimpulse_table} ) {
        my (
            $X,         $Y,         $Z,        $rpint_pos, $rpint_neg,
            $z_pose_rs, $z_nege_rs, $Qint_pos, $rovp_min,  $z_pmin_rs,
            $ke_pos,    $work_pos,  $Disp_pos
        ) = @{$rrow};
        my $rs        = exp($X);
        my $t_pose_rs = $rs - $z_pose_rs;
        my $t_nege_rs = $rs - $z_nege_rs;
        my $t_pmin_rs = $rs - $z_pmin_rs;
        my $T_pose_rs = log($t_pose_rs);
        my $T_pmin_rs = log($t_pmin_rs);
        my $T_nege_rs = log($t_nege_rs);
        push @{$rzpose_table}, [ $T_pose_rs, $z_pose_rs ];
        push @{$rzpmin_table}, [ $T_pmin_rs, $z_pmin_rs ];
        push @{$rznege_table}, [ $T_nege_rs, $z_nege_rs ];
    }
    return [ $rzpose_table, $rzpmin_table, $rznege_table ];
}

sub _update_toa_table {
    my ($rtable) = @_;

    # Update the values of T and dTdX, where T=ln(time of arrival)
    # in a table of shock values
    foreach my $row ( @{$rtable} ) {
        _update_toa_row($row);
    }
    return;
}

sub get_z {
    my ( $self, $T, $rztab, $z_default ) = @_;

    # Given T, lookup the value of z in a table of the form [T,z]
    # return -t before the start of the table
    # return the supplied default beyond the end of the table
    return unless ($rztab);
    my $jl = $self->{_jl2};
    my $ju = $self->{_ju2};

    # Locate this point in the table
    my $ntab  = @{$rztab};
    my $rhash = {
        _jl     => $jl,
        _ju     => $ju,
        _rtable => $rztab,
    };
    ( $jl, $ju ) = _locate_2d( $rhash, $T, 0 );
    my $zz;
    my $tt = exp($T);
    if ( $jl < 0 ) {

        # before start of table;
        # assume pos phase ends at origin (r=0), so
        $zz = -$tt;
    }
    elsif ( $ju >= $ntab ) {

        # off end of table; use supplied default
        $zz = $z_default;
    }
    else {
        my $rrow = _interpolate_table_rows( $T, $rztab, $jl );
        if ($rrow) {
            $zz = $rrow->[1];
        }
    }
    return $zz;
}

sub get_impulse {
    my ( $self, $result ) = @_;
    my $rimpulse_table = $self->{_rimpulse_table};
    return unless ($rimpulse_table);
    my $jl       = $self->{_jl2};
    my $ju       = $self->{_ju2};
    my $symmetry = $self->{_symmetry};
    my $gamma    = $self->{_gamma};
    my $rreturn_hash;
    my ( $X, $Y, $dYdX, $Z, $dZdX ) =
      @{$result}[ I_X, I_Y, I_dYdX, I_Z, I_dZdX ];

    # Locate this point in the table
    my $ntab  = @{$rimpulse_table};
    my $rhash = {
        _jl     => $jl,
        _ju     => $ju,
        _rtable => $rimpulse_table,
    };
    ( $jl, $ju ) = _locate_2d( $rhash, $X, 0 );
    $self->{_jl2} = $jl;
    $self->{_ju2} = $ju;
    return unless ( defined($jl) );

    my (
        $rpint_pos, $rpint_neg, $z_pose_rs, $z_nege_rs, $Qint_pos,
        $rovp_min,  $z_pmin_rs, $ke_pos,    $work_pos,  $Disp_pos,
    );

    # Handle case before start of the table
    if ( $jl < 0 ) {

        my $rrow = $rimpulse_table->[0];
        my (
            $X_b,         $Y_b,         $Z_b,         $rpint_pos_b,
            $rpint_neg_b, $z_pose_rs_b, $z_nege_rs_b, $Qint_pos_b,
            $rovp_min_b,  $z_pmin_rs_b, $ke_pos_b,    $work_pos_b,
            $Disp_pos_b,
        ) = @{$rrow};

        # ke_pos becomes constant (the KE of the similarity solution)
        $ke_pos = $ke_pos_b;

        # positive work equals total work since there is no negative work
        ### my ( $X, $Y, $dYdX, $Z, $dZdX, $T, $dTdX, $ZmX, $dZmXdX, $E1, $dE1dX, $Er, $dErdX ) = @{$result};
        my $Er = $result->[I_Er];
        $work_pos = 1 - $gamma * $Er;

        # ovp_min becomes constant (and equal to min central pressure)
        my $r_b    = exp($X_b);
        my $r      = exp($X);
        my $rpow_b = $r_b**( $symmetry / 2 );
        my $rpow   = $r**( $symmetry / 2 );
        my $rprat  = ( $r / $r_b )**( $symmetry / 2 );
        $rovp_min = $rovp_min_b * $rprat;

        # time of ovp_min becomes a constant (time of central pressure min), so
        # t=r-z constant or z-r=constant or
        $z_pmin_rs = $z_pmin_rs_b + $r - $r_b;

        # pint_neg becomes constant (and equal to min central pressure)
        $rpint_neg = $rpint_neg_b * $rprat;

        # z coordinates of zero overpressure lines are constant at early time
        # (they approach the times of the central pressure)
        $z_pose_rs = $z_pose_rs_b;
        $z_nege_rs = $z_nege_rs_b;

        # We will use log extrapolation for these variables back to the origin:
        $Disp_pos  = $Disp_pos_b;
        $Qint_pos  = $Qint_pos_b;
        $rpint_pos = $rpint_pos_b;

# Here are the theoretical slopes for dQint_pos/dX from the similarity solution.
# The table values are very close to these.
# -0.5 for spherical symmetry
#    0 for cylindrical symmetry
# +0.5 for plane symmetry

      # Positive phase impulse and displacement are harder to do theoretically
      # because they end when absolute pressure hits 1, so the integral does not
      # go to infinity.  So it seems best to just use log slopes of the computed
      # values to continue the solution before the first table point.

        if ( @{$rimpulse_table} > 2 ) {

            # Get the next row so that we can calculate the slope
            my $rrow_a = $rimpulse_table->[1];
            my (
                $X_a,         $Y_a,         $Z_a,         $rpint_pos_a,
                $rpint_neg_a, $z_pose_rs_a, $z_nege_rs_a, $Qint_pos_a,
                $rovp_min_a,  $z_pmin_rs_a, $ke_pos_a,    $work_pos_a,
                $Disp_pos_a,
            ) = @{$rrow_a};

            my $F_a = log($rpint_pos_a);
            my $F_b = log($rpint_pos_b);

            my $dX   = ( $X_a - $X_b );
            my $dQdX = ( $Qint_pos_a - $Qint_pos_b ) / $dX;
            my $dDdX = ( $Disp_pos_a - $Disp_pos_b ) / $dX;
            my $dFdX = ( $F_a - $F_b ) / $dX;

            $Qint_pos = $Qint_pos_b + $dQdX * ( $X - $X_b );
            $Disp_pos = $Disp_pos_b + $dDdX * ( $X - $X_b );
            my $FF = $F_b + $dFdX * ( $X - $X_b );
            $rpint_pos = exp($FF);
        }
    }

    # Handle case beyond end of the table
    elsif ( $ju >= $ntab ) {

        my $rrow = $rimpulse_table->[-1];
        my (
            $X_e,         $Y_e,         $Z_e,         $rpint_pos_e,
            $rpint_neg_e, $z_pose_rs_e, $z_nege_rs_e, $Qint_pos_e,
            $rovp_min_e,  $z_pmin_rs_e, $ke_pos_e,    $work_pos_e,
            $Disp_pos_e,
        ) = @{$rrow};

        # FIXME: rovp_min and z_pmin_rs at long range require some work
        # Set Tentative values for now
        $z_pmin_rs = $z_pmin_rs_e;
        $rovp_min  = $rovp_min_e;

# Qint_pos varis in proportion to Sigma/r^n at long range
# so ln(q_pos) is proportional to ln(p*r^n/2 r^-n) = ln(p)-n/2 ln(r) = Y-(n/2)X
# also note that this implies dln(Q)/dX = dY/dX-n/2 [agrees with tables fairly well]
        $Qint_pos = $Qint_pos_e + $Y - $Y_e - ( $symmetry / 2 ) * ( $X - $X_e );

        # radius^(n/2) x impulse becomes constant at long range
        $rpint_pos = $rpint_pos_e;
        $rpint_neg = $rpint_neg_e;

        # z coordinates of zero overpressure lines become constant at long range
        $z_pose_rs = $z_pose_rs_e;
        $z_nege_rs = $z_nege_rs_e;

# peak displacement times r^(n/2) becomes constant and equals the integral of Sigma
# so ln{disp * (n/2)r} = const = ln(disp)+(n/2)X = Disp_e + n/2 X_e
        $Disp_pos = $Disp_pos_e + ( $symmetry / 2 ) * ( $X_e - $X );

        # KE and Work are proportional to Sigma = ovp*r^(n/2), so
        # ln(W) - ln(We) = ln(ovp)-ln(ovp_e) + (n/2)*(X-X_e)
        # ln(W) = ln(We) + Y-Y + (n/2)*(X-X_e)
        my $dW = $Y - $Y_e + ( $symmetry / 2 ) * ( $X - $X_e );
        if ( $work_pos_e && $ke_pos_e ) {
            $work_pos = exp( log($work_pos_e) + $dW );
            $ke_pos   = exp( log($ke_pos_e) + $dW );
        }
    }

    # otherwise interpolate within table
    else {
        my $rrow = _interpolate_table_rows( $X, $rimpulse_table, $jl );
        if ($rrow) {
            (
                my $X_i,    my $Y_i,    my $Z_i,   $rpint_pos, $rpint_neg,
                $z_pose_rs, $z_nege_rs, $Qint_pos, $rovp_min,  $z_pmin_rs,
                $ke_pos,    $work_pos,  $Disp_pos
            ) = @{$rrow};
        }
    }
    my $disp_pos = defined($Disp_pos) ? exp($Disp_pos) : undef;
    my $qint_pos = defined($Qint_pos) ? exp($Qint_pos) : undef;
    $rreturn_hash->{rpint_pos} = $rpint_pos;
    $rreturn_hash->{rpint_neg} = $rpint_neg;
    $rreturn_hash->{rovp_min}  = $rovp_min;
    $rreturn_hash->{z_pmin_rs} = $z_pmin_rs;
    $rreturn_hash->{z_pose_rs} = $z_pose_rs;
    $rreturn_hash->{z_nege_rs} = $z_nege_rs;
    $rreturn_hash->{KE_pos}    = $ke_pos;
    $rreturn_hash->{W_pos}     = $work_pos;
    $rreturn_hash->{disp_pos}  = $disp_pos;
    $rreturn_hash->{qint_pos}  = $qint_pos;
    return $rreturn_hash;
}

sub _interpolate_table_rows {

    my ( $xx, $rtable, $jl, $icx, $NLAG ) = @_;

    # interpolate all row values of a table to an arbitrary x value

    # $xx = the value of the X variable
    # $rtable   = table
    # $jl = lower bound row index
    # $icx = column number of X variable [default is 0]
    # $NLAG = #Lagrange points [default 4]

    $icx  = 0 unless defined($icx);
    $NLAG = 4 unless defined($NLAG);

    # number of lagrange points cannot exceed number of table points
    my $npts = @{$rtable};
    if ( $NLAG > $npts ) { $NLAG = $npts }
    my $NH = int( $NLAG / 2 );

    my $jpoly_l = -1;

    my $rrow;

    if ( $jl < 0 || $jl + 1 >= $npts ) {

        # CALL ERROR: should not happen with proper call
        return;

    }

    # We want to have about equal numbers around this point
    $jpoly_l = $jl - $NH + 1;

    # But keep the points within the old table
    my $jpoly_u = $jpoly_l + $NLAG - 1;
    my $dj = $jpoly_u - ( $npts - 1 );
    if ( $dj > 0 ) {
        $jpoly_l -= $dj;
    }
    $jpoly_l = 0 if ( $jpoly_l < 0 );

    my $nvars = @{ $rtable->[0] };

    my $rlag_x = [];
    my $jj     = $jpoly_l;
    for ( my $n = 0 ; $n < $NLAG ; $n++ ) {
        last if ( $jj > $npts - 1 );    ## For safety
        push @{$rlag_x}, $rtable->[$jj]->[$icx];
        $jj++;
    }

    # now loop to define all values
    # Note that we are also interpolating the x variable as a control
    for ( my $icy = 0 ; $icy < $nvars ; $icy++ ) {
        my $rlag_y = [];
        my $jj     = $jpoly_l;
        for ( my $n = 0 ; $n < $NLAG ; $n++ ) {
            last if ( $jj > $npts - 1 );    ## For safety
            push @{$rlag_y}, $rtable->[$jj]->[$icy];
            $jj++;
        }
        my $yy = polint( $xx, $rlag_x, $rlag_y );
        $rrow->[$icy] = $yy;
    }
    return $rrow;
}

sub get_phase_lengths {

    # FIXME: rewrite this, combining both get_ routines

    my ( $self, $rs, $zs ) = @_;
    my $symmetry = $self->{_symmetry};
    my $Tneg     = 0;
    my $Lneg     = 0;
    my ( $Tpos, $Lpos ) = $self->get_positive_phase( $rs, $zs );
    if ( $symmetry == 2 ) {
        ( $Tneg, $Lneg ) = $self->get_negative_phase( $rs, $zs );
        $Tneg -= $Tpos if ($Tneg);
        $Lneg -= $Lpos if ($Lneg);

        # Approximation for very long range as an N wave forms:
        # At very long range the negative phase length increases
        # due to the second shock so that the negative phase
        # duration does not get less than the positive phase
        if ( $Tneg < $Tpos ) { $Tneg = $Tpos }
        if ( $Lneg < $Lpos ) { $Lneg = $Lpos }
    }
    return ( $Tpos, $Lpos, $Tneg, $Lneg );
}

sub _set_global_info {
    my ($self) = @_;

    # Given a shock radius and z value, return the positive phase duration and
    # length
    my $symmetry        = $self->{_symmetry};
    my $gamma           = $self->{_gamma};
    my $rhash           = get_blast_info( $symmetry, $gamma );
    my $loc_gamma_table = $self->{_loc_gamma_table};
    my $rigam_4         = $self->{_rigam_4};
    if ( !defined($rhash) && defined($loc_gamma_table) && defined($rigam_4) ) {

        # If this table is an interpolated table then we can interpolate
        # P0 from nearby fits.
        # FIXME: this works now but needs a lot of refinement for better
        # Use best transformation for each variable.
        # Sintegral needs (alpha+2) weighting.
        my ( $jj, $jl, $ju, $num ) = @{$loc_gamma_table};

        #if ( defined($loc_gamma_table) && defined($rigam_4) ) {
        if ( $jl >= 0 && $ju < $num ) {

            my $alpha = alpha_interpolate( $symmetry, $gamma );
            my $rtab;
            foreach my $igam ( @{$rigam_4} ) {
                my ( $gamma_i, $table_name_i ) =
                  @{ $rgamma_table->[$symmetry]->[$igam] };
                my $rhash_i = get_blast_info( $symmetry, $gamma_i );
                my $alpha_i = alpha_interpolate( $symmetry, $gamma_i );
                push @{$rtab}, [ $rhash_i, $gamma_i, $alpha_i, $igam ];
            }

            # loop to interpolate all values
            my @keys = keys %{ $rtab->[0]->[0] };
            foreach my $key (@keys) {
                my ( $rx, $ry );
                my $nogo;
                foreach my $item ( @{$rtab} ) {
                    my ( $rhash_i, $gamma_i, $alpha_i, $igam ) = @{$item};
                    my $val_i = $rhash_i->{$key};
                    if ( !defined($val_i) || $val_i == 0 ) {

                        # no interpolation if missing or zero bounding points
                        $nogo ||= ( $igam eq $jl || $igam eq $ju );
                        last if $nogo;
                        next;
                    }

                    # FIXME: Sintegral should be weighted with (alpha+2)
                    push @{$rx}, $gamma_i;
                    push @{$ry}, $val_i;
                }

                my $val = 0;
                if ( !$nogo ) {
                    $val = polint( $gamma, $rx, $ry );
                }
                $rhash->{$key} = $val;
            }

            # Old Linear interpolation, for reference
            # execute to see comparison with polynomial interpolation
            if (0) {
                my ( $jj, $jl, $ju, $num ) = @{$loc_gamma_table};
                if ( $jl >= 0 && $ju < $num ) {
                    my $rtab    = $rgamma_table->[$symmetry];
                    my $gamma_l = $rtab->[$jl]->[0];
                    my $gamma_u = $rtab->[$ju]->[0];
                    my $rhash_l = get_blast_info( $symmetry, $gamma_l );
                    my $rhash_u = get_blast_info( $symmetry, $gamma_u );

                    if ( defined($rhash_l) && defined($rhash_u) ) {
                        foreach my $key ( keys %{$rhash_l} ) {
                            my $val_l = $rhash_l->{$key};
                            my $val_u = $rhash_u->{$key};

                            if ( defined($val_l) && defined($val_u) ) {
                                my $val = _linear_interpolation(
                                    $gamma,   $gamma_l, $val_l,
                                    $gamma_u, $val_u
                                );
                                my $vsave = $rhash->{$key};
                                $rhash->{$key} = $val;

                                # DEBUG: comparison of methods
                                print STDERR
                                  "key=$key polint=$vsave linear=$val\n";
                            }
                        }
                    }
                }
            }
        }
    }

    # Store this hash of values in self
    $self->{_rinfo} = $rhash;
    return;
}

sub get_positive_phase {
    my ( $self, $rs, $zs ) = @_;

    # Given a shock radius and z value, return the positive phase duration and
    # length
    my $symmetry = $self->{_symmetry};
    my $gamma    = $self->{_gamma};
    my $rpz_fit  = $rpzero_fit->{$symmetry}->{$gamma};
    my ( $Tpos, $Lpos ) = ( 0, 0 );

    if ( defined($rpz_fit) ) {

        # Handle builtin table
        ( $Tpos, $Lpos ) = _positive_phase( $rs, $zs, $symmetry, $rpz_fit );
    }
    else {

        # No P0 fit is defined for this gamma; if this table is an interpolated
        # table then we can interpolate P0 from nearby fits.
        my $loc_gamma_table = $self->{_loc_gamma_table};
        if ( defined($loc_gamma_table) ) {
            my ( $jj, $jl, $ju, $num ) = @{$loc_gamma_table};
            if ( $jl >= 0 && $ju < $num ) {
                my $rtab      = $rgamma_table->[$symmetry];
                my $gamma_l   = $rtab->[$jl]->[0];
                my $gamma_u   = $rtab->[$ju]->[0];
                my $rpz_fit_l = $rpzero_fit->{$symmetry}->{$gamma_l};
                my $rpz_fit_u = $rpzero_fit->{$symmetry}->{$gamma_u};
                if ( defined($rpz_fit_l) && defined($rpz_fit_u) ) {
                    my ( $Tpos_l, $Lpos_l ) =
                      _positive_phase( $rs, $zs, $symmetry, $rpz_fit_l );
                    my ( $Tpos_u, $Lpos_u ) =
                      _positive_phase( $rs, $zs, $symmetry, $rpz_fit_u );
                    $Tpos =
                      _linear_interpolation( $gamma, $gamma_l, $Tpos_l,
                        $gamma_u, $Tpos_u );
                    $Lpos =
                      _linear_interpolation( $gamma, $gamma_l, $Lpos_l,
                        $gamma_u, $Lpos_u );
                }
            }
        }
    }
    return wantarray ? ( $Tpos, $Lpos ) : $Tpos;
}

sub _positive_phase {
    my ( $rs, $zs, $symmetry, $rpz_fit ) = @_;

    # Given:
    # $rs - a shock radius
    # $zs - the z at that radius
    # $symmetry - 0, 1, or 2
    # $rpz_fit - the zero pressure fit parameters

    # return:
    # $Tpos - the positive phase duration (in time)
    # $Lpos - the positive phase length (in space)

    # Note: Returns 0's if no positive phase
    # $Tpos = 0 if no positive phase duration
    # $Lpos = 0 if no positive phase length
    my $Tpos = 0;
    my $Lpos = 0;
    if ( defined($rpz_fit) && @{$rpz_fit} ) {
        my ( $tPc0, $t0, $r0, $tmin_fit, $rmin_fit, $zinf, $AA, $BB ) =
          @{$rpz_fit};

        my $AH = $symmetry / 2;

        # First find the positive phase duration
        my $z = $zs;
        if ( $rs >= $rmin_fit ) {

            # in analytical region
            my $r     = $rs;
            my $rpow  = $r**$AH;
            my $denom = $rpow - $BB + $AA;
            if ( $denom > 0 ) { $z = $zinf * ( 1 - $AA / $denom ); }
        }
        else {

            # before analytical region; use straight line to origin
            my $zc       = 0 - $tPc0;
            my $zmin_fit = $rmin_fit - $tmin_fit;
            $z = $zc + ( $zmin_fit - $zc ) * $rs / $rmin_fit;
            if ( $z > $zs ) { $z = $zs }
        }
        $Tpos = ( $zs - $z );

        my $ts = $rs - $zs;
        if ( $ts >= $tmin_fit ) {

            # in analytical section...
            # We must iterate to find the positive length
            # We start with the previous z
            my $r     = $rs - ( $zs - $z );
            my $rpow  = $r**$AH;
            my $denom = $rpow - $BB + $AA;
            if ( $denom > 0 ) { $z = $zinf * ( 1 - $AA / $denom ); }

            # re-iterate
            $r     = $rs - ( $zs - $z );
            $rpow  = $r**$AH;
            $denom = $rpow - $BB + $AA;
            if ( $denom > 0 ) { $z = $zinf * ( 1 - $AA / $denom ); }
            $Lpos = $zs - $z;
        }
        elsif ( $ts >= $tPc0 ) {

            # before analytical section...connect line to origin
            my $zc       = 0 - $tPc0;
            my $zmin_fit = $rmin_fit - $tmin_fit;
            $z =
              $zc +
              ( $zmin_fit - $zc ) * ( $ts - $tPc0 ) / ( $tmin_fit - $tPc0 );
            if ( $z > $zs ) { $z = $zs }
            $Lpos = $zs - $z;
        }
        else {

            # before tPc0: no positive length
        }
    }
    return wantarray ? ( $Tpos, $Lpos ) : $Tpos;
}

sub get_negative_phase {
    my ( $self, $rs, $zs ) = @_;

    # Given a shock radius and z value, return the negative phase duration and
    # length. This is zero except for spherical symmetry
    my $symmetry = $self->{_symmetry};
    my $gamma    = $self->{_gamma};
    my $rpz_tail = $rpzero_tail->{$symmetry}->{$gamma};
    my ( $Tneg, $Lneg ) = ( 0, 0 );

    if ( defined($rpz_tail) ) {

        # Handle builtin table
        ( $Tneg, $Lneg ) = _negative_phase( $rs, $zs, $symmetry, $rpz_tail );
    }
    else {

        # No P0 fit is defined for this gamma; if this table is an interpolated
        # table then we can interpolate P0 from nearby fits.
        my $loc_gamma_table = $self->{_loc_gamma_table};
        if ( defined($loc_gamma_table) ) {
            my ( $jj, $jl, $ju, $num ) = @{$loc_gamma_table};
            if ( !$jj && $jl >= 0 && $ju < $num ) {
                my $rtab       = $rgamma_table->[$symmetry];
                my $gamma_l    = $rtab->[$jl]->[0];
                my $gamma_u    = $rtab->[$ju]->[0];
                my $rpz_tail_l = $rpzero_tail->{$symmetry}->{$gamma_l};
                my $rpz_tail_u = $rpzero_tail->{$symmetry}->{$gamma_u};
                if ( defined($rpz_tail_l) && defined($rpz_tail_u) ) {
                    my ( $Tneg_l, $Lneg_l ) =
                      _negative_phase( $rs, $zs, $symmetry, $rpz_tail_l );
                    my ( $Tneg_u, $Lneg_u ) =
                      _negative_phase( $rs, $zs, $symmetry, $rpz_tail_u );
                    $Tneg =
                      _linear_interpolation( $gamma, $gamma_l, $Tneg_l,
                        $gamma_u, $Tneg_u );
                    $Lneg =
                      _linear_interpolation( $gamma, $gamma_l, $Lneg_l,
                        $gamma_u, $Lneg_u );
                }
            }
        }
    }
    return wantarray ? ( $Tneg, $Lneg ) : $Tneg;
}

sub _negative_phase {
    my ( $rs, $zs, $symmetry, $rpz_tail ) = @_;

    # Given:
    # $rs - a shock radius
    # $zs - the z at that radius
    # $symmetry - 0, 1, or 2
    # $rpz_tail - the zero pressure fit parameters

    # return:
    # $Tneg - time to arrival of second zero p
    # $Lneg - distance to arrival of second zero p

    # Returns 0's if no negative phase:
    # $Tneg = 0 if no negative phase duration
    # $Lneg = 0 if no negative phase length
    my $Tneg = 0;
    my $Lneg = 0;
    if ( defined($rpz_tail) && @{$rpz_tail} ) {
        my ( $t0, $r0, $ztail ) = @{$rpz_tail};
        if ( $rs >= $r0 ) { $Tneg = ( $zs - $ztail ); }
        my $ts = $rs - $zs;
        if ( $ts >= $t0 ) { $Lneg = $zs - $ztail }
    }
    return wantarray ? ( $Tneg, $Lneg ) : $Tneg;
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

sub wavefront {

    # Lookup wavefront properties at a given range.
    # This is the main entry for external calls.

    # Given:
    #     a hash of the form
    #     id_Q => Q

    # where
    #     Q = a value to lookup

    #     id_Q defines what Q contains:
    #       'X'    if Q = ln(R), where R=scaled range
    #       'Y'    if Q = ln(overpressure ratio)
    #       'dYdX' if Q = dY/dX = d(ln p)/d(ln r)
    #       'Z'    if Q = ln(R-cT), where T=scaled time of arrival
    #       'dZdX' if Q = dZ/dX
    #       'T'    if Q = ln(T), where T=scaled time of arrival
    #       'dTdX' if Q = dT/dX
    #       'Z-X'  if Q = Z-X
    #       'dZ-XdX' if Q = d(Z-X)/dX

    # The hash may also contain (for debugging):
    #      'interpolation_flag' => $interp
    #	     interp = 0 for cubic [This is the default and recommended method]
    #	     interp = 1 for linear [This is only intended for error checking]

    # Returns a hash with numerous quantities of interest

    my ( $self, @args ) = @_;

    # do not proceed unless the table is good
    my $error = $self->{_error};
    if ($error) {
        carp "$error";
        return;
    }
    my $gamma    = $self->{_gamma};
    my $Sint_pos = $self->{_rinfo}->{Sintegral_pos};
    my $Sint_neg = $self->{_rinfo}->{Sintegral_neg};

    my $rtab = $self->{_rtable};
    my $ntab = @{$rtab};

    if ( @args == 0 ) {
        croak "Missing call parameters";
    }

    # default input is a hash
    my %input_hash  = @args;
    my $rinput_hash = \%input_hash;

    # but also allow a hash ref
    if ( @args == 1 ) {
        my $arg0    = $args[0];
        my $reftype = ref($arg0);
        if ( $reftype && $reftype eq 'HASH' ) {
            $rinput_hash = $arg0;
        }
    }

    # FIXME: move to initialization block
##?    @q = qw(
##?      grep
##?      keys
##?      map
##?      reverse
##?      sort
##?      split
##?    );
##?    @is_keyword_returning_list{@q} = (1) x scalar(@q);
    # For interpolation variables, this lists the corresponding table columns.
    # For the interpolation flag it lists the default value
    my %valid_input_keys = (
        X                  => I_X,
        Y                  => I_Y,
        dYdX               => I_dYdX,
        Z                  => I_Z,
        dZdX               => I_dZdX,
        T                  => I_T,
        dTdX               => I_dTdX,
        'Z-X'              => I_ZmX,
        'dZ-XdX'           => I_dZmXdX,
        E1                 => I_E1,       # primary shock residual energy
        dE1dX              => I_dE1dX,
        Er                 => I_Er,       # total shock residual energy
        dErdX              => I_dErdX,
        interpolation_flag => 0,
    );


    # Validate input keys
    _check_keys( $rinput_hash, \%valid_input_keys,
        "Checking for valid input keys" );

    # this interpolation flag is for debugging;
    # the default and recommended interpolation is cubic
    my $interp = $valid_input_keys{'interpolation_flag'};    #0;

    # get table column number and value to search
    my $icol;
    my $Q;
    my $id_Q;
    foreach my $key ( keys %{$rinput_hash} ) {
        my $val = $rinput_hash->{$key};
        if ( $key eq 'interpolation_flag' ) { $interp = $val }
        else {

            if ( defined($Q) ) {
                croak(
"Expecting only one input value but saw both $id_Q and $key\n"
                );
            }
            $Q    = $val;
            $id_Q = $key;
            $icol = $valid_input_keys{$key};
        }
    }

    # Locate this point in the table
    my ( $il, $iu ) = $self->_locate_2d( $Q, $icol );

    my $result;

    # Handle case before start of the table
    if ( $il < 0 ) {
        $result = $self->_short_range_calc( $Q, $icol );
    }

    # Handle case beyond end of the table
    elsif ( $iu >= $ntab ) {
        $result = $self->_long_range_calc( $Q, $icol );
    }

    # otherwise interpolate within table
    else {
        $result =
          _interpolate_rows( $Q, $icol, $rtab->[$il], $rtab->[$iu], $interp );
    }

    my (
        $X,   $Y,      $dYdX, $Z,     $dZdX, $T, $dTdX,
        $ZmX, $dZmXdX, $E1,   $dE1dX, $Er,   $dErdX
      )
      = @{$result}[
      I_X,      I_Y,  I_dYdX,  I_Z,  I_dZdX, I_T, I_dTdX, I_ZmX,
      I_dZmXdX, I_E1, I_dE1dX, I_Er, I_dErdX
      ];

    # PATCH: Avoid trouble in case user adds a table without E variables
    if ( !defined($E1) ) { $E1 = 0; $dE1dX = 0 }
    if ( !defined($Er) ) { $Er = 0; $dErdX = 0 }

    # Fix possible minor problem in which slight differences in interpolation
    # can cause Er to be slightly below E1 at the threshold of a tail shock.
    if ( $Er < $E1 ) { $Er = $E1 }

    my $W_blast = 1 - $gamma * $Er;
    my $W_atm   = ( $gamma - 1 ) * $Er;

    # Get impulse and related values
    my $rimpulse_hash = $self->get_impulse($result);
    my (
        $rpint_pos, $rpint_neg, $rovp_min, $z_pmin_rs, $z_pose_rs,
        $z_nege_rs, $KE_pos,    $W_pos,    $disp_pos,  $qint_pos
    );
    if ($rimpulse_hash) {
        $rpint_pos = $rimpulse_hash->{rpint_pos};
        $rpint_neg = $rimpulse_hash->{rpint_neg};
        $rovp_min  = $rimpulse_hash->{rovp_min};
        $z_pmin_rs = $rimpulse_hash->{z_pmin_rs};
        $z_pose_rs = $rimpulse_hash->{z_pose_rs};
        $z_nege_rs = $rimpulse_hash->{z_nege_rs};
        $KE_pos    = $rimpulse_hash->{KE_pos};
        $W_pos     = $rimpulse_hash->{W_pos};
        $disp_pos  = $rimpulse_hash->{disp_pos};
        $qint_pos  = $rimpulse_hash->{qint_pos};
    }

    # FIXME: evaluate sigma and Sigma = sigma*r**(symmetry/2)

    # find z coordinates along a profile in space which includes this shock
    my ( $z_pose_ts, $z_pmin_ts, $z_nege_ts ) =
      ( $z_pose_rs, $z_pmin_rs, $z_nege_rs );

    my $rztables = $self->{_rztables};
    if ( defined($rztables) ) {
        my ( $rzpose_table, $rzpmin_table, $rznege_table ) = @{$rztables};
        $z_pose_ts = $self->get_z( $T, $rzpose_table, $z_pose_ts );
        $z_pmin_ts = $self->get_z( $T, $rzpmin_table, $z_pmin_ts );
        $z_nege_ts = $self->get_z( $T, $rznege_table, $z_nege_ts );
    }

    # backup analytical fit to zero pressure lines
    my $rs = exp($X);
    my $zs = exp($Z);
    my ( $Tpos, $Lpos, $Tneg, $Lneg ) = $self->get_phase_lengths( $rs, $zs );

    # FUTURE: use more accurate values if possible
    #if ( defined($z_pose_rs) ) { $Tpos = $zs - $z_pose_rs; }
    #if ( defined($z_pose_ts) ) { $Lpos = $zs - $z_pose_ts; }

    # TableLoc  shows which table rows were interpolated
    # TableVars shows the interpolated row values

    my $return_hash = {
        'X'       => $X,
        'Y'       => $Y,
        'dYdX'    => $dYdX,
        'Z'       => $Z,
        'dZdX'    => $dZdX,
        'T'       => $T,
        'dTdX'    => $dTdX,
        'E1'      => $E1,
        'dE1dX'   => $dE1dX,
        'Er'      => $Er,
        'dErdX'   => $dErdX,
        'W_blast' => $W_blast,
        'W_atm'   => $W_atm,
        'KE_pos'  => $KE_pos,
        'W_pos'   => $W_pos,
        'Tpos'    => $Tpos,
        'Lpos'    => $Lpos,

        #'Lneg'        => $Lneg,
        'TableLoc'    => [ $il, $iu, $ntab ],
        'TableVars'   => $result,
        'Ixr_pos_lim' => $Sint_pos * $gamma,
        'Ixr_neg_lim' => $Sint_neg * $gamma,
        'Sint_pos'    => $Sint_pos,
        'Sint_neg'    => $Sint_neg,
        'Ixr_pos'     => $rpint_pos,
        'Ixr_neg'     => $rpint_neg,
        'rovp_min_rs' => $rovp_min,
        'z_pmin_rs'   => $z_pmin_rs,
        'z_pose_rs'   => $z_pose_rs,
        'z_nege_rs'   => $z_nege_rs,
        'disp_pos'    => $disp_pos,
        'qint_pos'    => $qint_pos,
        'z_pmin_ts'   => $z_pmin_ts,
        'z_pose_ts'   => $z_pose_ts,
        'z_nege_ts'   => $z_nege_ts,
    };
    return $return_hash;
}

sub lookup {

    #####################################################################
    # This is the OLD main entry for external calls.
    # Please call 'wavefront' instead; it will be more general.
    # TO BE DELETED!
    #####################################################################

    # Given:
    #     $Q = a value to lookup

    #     $id_Q defines what Q contains:
    #       0 or 'X' => ln(R)    [this is the default if $id_Q not specified]
    #       1 or 'Y' => ln(overpressure)
    #       3 or 'Z' => ln(z=x-ct), where t=time of arrival
    #  OLD  5 or 'W' => ln(T)
    #  NEW  5 or 'T' => ln(T)

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
        elsif ( $id_Q eq 'W' ) { $icol = 5 }    # OLD
        elsif ( $id_Q eq 'T' ) { $icol = 5 }    # NEW
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
        push @{$rtable_gen}, $self->wavefront( 'X' => $X )->{'TableVars'};
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

    # FIXME:
    # for icol=2,4,6  we should just return the first line of the table

    my $ovprat_from_lambda = sub {
        my ($lambda) = @_;

        # FIXME: patch to avoid divide by 0
        if ( $lambda < 1.e-80 ) { $lambda = 1.e-80 }

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
        # theory".  That printed equation is missing a factor of gamma.
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

    my $d2Y_dX2_i = 0;

    # FIXME: This still uses simple theory
    my $lnT_i     = $lnt_0 + $delta * ( $X_i - $X_0 );
    my $T_i       = exp($lnT_i);
    my $lambda_i  = exp($X_i);
    my $z_i       = $lambda_i - $T_i;
    my $dz_dX_i   = $lambda_i - $T_i * $delta;
    my $d2z_dX2_i = $lambda_i - $T_i * ($delta)**2;
    my $Z_i       = $z_i;
    my $dZ_dX_i   = $dz_dX_i;
    $Z_i     = log($z_i);
    $dZ_dX_i = $dz_dX_i / $z_i;

    my $dzdX_i = $z_i * $dZ_dX_i;
    my $dTdX_i =
      $lambda_i != $z_i ? ( $lambda_i - $dzdX_i ) / ( $lambda_i - $z_i ) : 0;

    # calculate $eresidual = residual heat energy from origin to shock
    # multiply by gamma-1 to get work on atmosphere
    # multiply by gamma to get energy not available for blast at this range
    my $symp   = $symmetry + 1;
    my $gammam = $gamma - 1;
    my $pi     = 4 * atan2( 1, 1 );
    my $acon =
        ( $symmetry == 2 ) ? 4 * $pi
      : ( $symmetry == 1 ) ? 2 * $pi
      :                      1;
    my $ovp_i  = exp($Y_i);
    my $C      = $ovp_i * $lambda_i**$symp;
    my $qq     = $symp * $gammam / $gamma;
    my $gammap = $gamma + 1;
    my $vol =
      $acon * $gammam / $gammap * $C**( 1 / $gamma ) * $lambda_i**$qq / $qq;
    my $vol1    = $acon * $lambda_i**$symp / $symp;
    my $E1_i    = ( $vol - $vol1 ) / $gammam;
    my $dE1dX_i = 0;                               # FIXME; use exact derivative
    my $Er_i    = $E1_i;
    my $dErdX_i = $dE1dX_i;

    # The result
    my $result_i = [
        $X_i,     $Y_i,    $dY_dX_i,    $Z_i,         $dZ_dX_i,
        $T_i,     $dTdX_i, $Z_i - $X_i, $dZ_dX_i - 1, $E1_i,
        $dE1dX_i, $Er_i,   $dErdX_i
    ];

    return $result_i;
}

sub _long_range_calc {
    my ( $self, $Q, $icol ) = @_;
    my $symmetry = $self->{_symmetry};

    # FIXME:
    # for icol=2,4,6  we should just return the first line of the table
    my $result;
    if ( $symmetry == 2 ) {
        $result = $self->_long_range_sphere( $Q, $icol );
    }
    else {
        $result = $self->_long_range_non_sphere( $Q, $icol );
    }

    # bring derived variables T and dT/dX up to date
    _update_toa_row($result);
    return $result;
}

sub _add_long_range_energy {
    my ( $self, $result_i ) = @_;
    my $rtab     = $self->{_rtable};
    my $symmetry = $self->{_symmetry};
    my $gamma    = $self->{_gamma};

    # Add the energy variables to the extrapolated result

    # current extrapolated state values
    my ( $X_i, $Y_i, $dYdX_i ) = @{$result_i}[ I_X, I_Y, I_dYdX ];

    # last table point is reference state
##    my (
##        $X_e,    $Y_e,  $dYdX_e,  $Z_e, $dZdX_e, $T_e,
##        $dTdX_e, $ZmX_e, $dZmXdX_e, $E1_e, $dE1dX_e, $Er_e, $dErdX_e
##    ) = @{ $rtab->[-1] };
    my (
        $X_e,   $Y_e,      $dYdX_e, $Z_e,     $dZdX_e, $T_e, $dTdX_e,
        $ZmX_e, $dZmXdX_e, $E1_e,   $dE1dX_e, $Er_e,   $dErdX_e
      )
      = @{ $rtab->[-1] }[
      I_X,    I_Y,   I_dYdX,   I_Z,  I_dZdX, I_T,
      I_dTdX, I_ZmX, I_dZmXdX, I_Er, I_dErdX
      ];

    # Tentative initialize to table end in case work is zero
    my $E1_i    = $E1_e;
    my $Er_i    = $Er_e;
    my $dE1dX_i = $dE1dX_e;
    my $dErdX_i = $dErdX_e;

    return unless defined($Er_i);

    # Work at end state
    my $w_e = 1 - $gamma * $Er_e;

    if ( $w_e > 0 ) {

        # Work and its slope at extrapolated point
        my $w_i =
          exp( log($w_e) + 0.5 * $symmetry * ( $X_i - $X_e ) + $Y_i - $Y_e );
        my $dwdX_i = $w_i * ( 0.5 * $symmetry + $dYdX_i );

        # total residual energy and its slope
        $Er_i    = ( 1 - $w_i ) / $gamma;
        $dErdX_i = -$dwdX_i / $gamma;

        # primary shock energy and its slope
        if ( $dE1dX_e > 0 ) {
            my $kk = $dErdX_e / $dE1dX_e;
            $E1_i = $E1_e + ( $Er_i - $Er_e ) / $kk;
            $dE1dX_i = $dErdX_i / $kk;
        }
    }

    $result_i->[I_E1]    = $E1_i;
    $result_i->[I_dE1dX] = $dE1dX_i;
    $result_i->[I_Er]    = $Er_i;
    $result_i->[I_dErdX] = $dErdX_i;
    return;
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

    # FIXME: use indexes
##?    my (
##?        $X_e,    $Y_e,  $dYdX_e,  $Z_e,  $dZdX_e, $T_e,
##?        $dTdX_e, $ZmX, $dZmXdX, $E1_e, $dE1dX_e, $E_e, $dEdX_e
##?    ) = @{ $rtab->[-1] };
    my (
        $X_e,   $Y_e,      $dYdX_e, $Z_e,     $dZdX_e, $T_e, $dTdX_e,
        $ZmX_e, $dZmXdX_e, $E1_e,   $dE1dX_e, $Er_e,   $dErdX_e
      )
      = @{ $rtab->[-1] }[
      I_X,    I_Y,   I_dYdX,   I_Z,  I_dZdX, I_T,
      I_dTdX, I_ZmX, I_dZmXdX, I_Er, I_dErdX
      ];

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

            # FIXME: not yet programmed to interpolate E1, E at long range
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

    # This is the partial result so far
    my $result_i = [ $X_i, $Y_i, $dY_dX_i, $Z_i, $dZ_dX_i ];

    # now add the long range energy variables
    $self->_add_long_range_energy($result_i);

##    ##############################################################################
##    # OLD: Here are equations to calculate the increment in residual shock
##    # heating from the primary shock, from 'energy_integral_primary.pl'.
##    # Add this to shock heading at end of table to get long range value
##    # This assumes constant dYdX, so it is basically exact for plane and
##    # cylindrical geometry. For spherical geometry dYdX continues to change
##    # but those tables go so far (X=200) that we probably won't be
##    # evaluating for spherical geometry.
##    my $lambda_i  = exp($X_i);
##    my $lambda_e  = exp($X_e);
##    my $ovp_atm_e = exp($Y_e);
##    my $qq        = 3 * $dYdX_e + $symmetry + 1;
##    my $pi        = 4 * atan2( 1, 1 );
##    my $acon =
##        ( $symmetry == 2 ) ? 4 * $pi
##      : ( $symmetry == 1 ) ? 2 * $pi
##      :                      1;
##    my $coef = $acon * ( $gamma + 1 ) / 12 / $gamma**3;
##    my $dint =
##      $coef * $ovp_atm_e**3 /
##      $lambda_e**( 3 * $dYdX_e ) *
##      ( $lambda_i**$qq - $lambda_e**$qq ) /
##      $qq;
##
##?    # OLD: future PATCH for total shock heating: double this for geometries
##?    # with a tail shock
##?    #unless ($symmetry == 0 && $gamma>1.6) {$dint*=2}
##?
##?    # OLD: Now add dint to residual energy at end of table to get final value
##?    ##############################################################################

    return $result_i;
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
    #   dY_dX is about -0.5 for sym=0 and -0.75 for sym=1
    #   dZ_dX is about +0.5 for sym=0 and +0.25 for sym=1
    # Here we will use these limiting values:
    if ( $symmetry == 0 ) {
        $dY_dX_i = -0.5;
        $dZ_dX_i = 0.5;
    }
    elsif ( $symmetry == 1 ) {
        $dY_dX_i = -0.75;
        $dZ_dX_i = 0.25;
    }

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

    # This is the partial result so far
    my $result_i = [ $X_i, $Y_i, $dY_dX_i, $Z_i, $dZ_dX_i ];

    # now add the long range energy variables
    $self->_add_long_range_energy($result_i);

    return $result_i;
}

sub _locate_2d {
    my ( $self, $xx, $icol, $rtab ) = @_;

    # Binary search for two consecutive table row indexes, jl and ju, of a
    # 2D matrix such that the value of x lies between these two table values.
    # $icol is the column of the variable x
    # If x is out of bounds, returns either jl<0 or ju>=N

    $rtab = $self->{_rtable} unless defined($rtab);
    return unless $rtab;
    my $num = @{$rtab};
    return unless $num;
    my $dx_is_positive = $rtab->[-1]->[$icol] > $rtab->[0]->[$icol];

    # Set search bounds to previous location ...
    my $jl = $self->{_jl};
    my $ju = $self->{_ju};

    # ... but reset to beyond end of table if no longer valid
    $jl = -1
      if ( !defined($jl)
        || $jl < 0
        || $jl >= $num
        || ( $xx > $rtab->[$jl]->[$icol] ne $dx_is_positive ) );
    $ju = $num
      if ( !defined($ju)
        || $ju < 0
        || $ju >= $num
        || ( $xx > $rtab->[$ju]->[$icol] eq $dx_is_positive ) );

    # Loop until the requested point lies in a single interval
    while ( $ju - $jl > 1 ) {
        my $jm = int( ( $jl + $ju ) / 2 );
        if ( $xx > $rtab->[$jm]->[$icol] eq $dx_is_positive ) {
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
    # 	[X,Y,Z,dYdX,dZdX,T,ZmX, dZmXdX, dTdX,E1, dE1dX, E, dEdX]

    my (
        $X_b,      $Y_b,     $dY_dX_b, $Z_b,       $dZ_dX_b,
        $T_b,      $dT_dX_b, $ZmX_b,   $dZmX_dX_b, $E1_b,
        $dE1_dX_b, $E_b,     $dE_dX_b
      )
      = @{$row_b}[
      I_X,      I_Y,  I_dYdX,  I_Z,  I_dZdX, I_T, I_dTdX, I_ZmX,
      I_dZmXdX, I_E1, I_dE1dX, I_Er, I_dErdX
      ];
    my (
        $X_e,      $Y_e,     $dY_dX_e, $Z_e,       $dZ_dX_e,
        $T_e,      $dT_dX_e, $ZmX_e,   $dZmX_dX_e, $E1_e,
        $dE1_dX_e, $E_e,     $dE_dX_e
      )
      = @{$row_e}[
      I_X,      I_Y,  I_dYdX,  I_Z,  I_dZdX, I_T, I_dTdX, I_ZmX,
      I_dZmXdX, I_E1, I_dE1dX, I_Er, I_dErdX
      ];

    # FIXME: for slope interpolation, we are currently doing liner interp
    # We should either iterate or do parabolic interpolation using the
    # second derivatives at the segment ends

    # If this is an inverse problem, first find X=X_i
    my $X_i;
    if ( $icol == I_X ) {
        $X_i = $Q;
    }

    # FIXME: these could be collapsed into two calls
    # if icol is odd do a cubic interpolation of icol with slope icol+1
    elsif ( $icol == I_Y ) {
        ( $X_i, my $dX_dY_i, my $d2X_dY2_i ) =
          _interpolate_scalar( $Q, $Y_b, $X_b, 1 / $dY_dX_b,
            $Y_e, $X_e, 1 / $dY_dX_e, $interp );
    }

    # if icol is even do a slope interpolation
    elsif ( $icol == I_dYdX ) {
        ( $X_i, my $slope1, my $slope2 ) =
          _interpolate_scalar( $Q, $dY_dX_b, $X_b, 0, $dY_dX_e, $X_e, 0, 1 );
    }
    elsif ( $icol == I_Z ) {
        ( $X_i, my $dX_dZ_i, my $d2X_dZ2_i ) =
          _interpolate_scalar( $Q, $Z_b, $X_b, 1 / $dZ_dX_b,
            $Z_e, $X_e, 1 / $dZ_dX_e, $interp );
    }
    elsif ( $icol == I_dZdX ) {
        ( $X_i, my $slope1, my $slope2 ) =
          _interpolate_scalar( $Q, $dZ_dX_b, $X_b, 0, $dZ_dX_e, $X_e, 0, 1 );
    }
    elsif ( $icol == I_T ) {
        ( $X_i, my $dX_dT_i, my $d2X_dT2_i ) =
          _interpolate_scalar( $Q, $T_b, $X_b, 1 / $dT_dX_b,
            $T_e, $X_e, 1 / $dT_dX_e, $interp );
    }
    elsif ( $icol == I_dTdX ) {
        ( $X_i, my $slope1, my $slope2 ) =
          _interpolate_scalar( $Q, $dT_dX_b, $X_b, 0, $dT_dX_e, $X_e, 0, 1 );
    }
    elsif ( $icol == I_ZmX ) {
        ( $X_i, my $dX_dZmX_i, my $d2X_dZmX2_i ) =
          _interpolate_scalar( $Q, $ZmX_b, $X_b, 1 / $dZmX_dX_b,
            $ZmX_e, $X_e, 1 / $dZmX_dX_e, $interp );
    }
    elsif ( $icol == I_dZmXdX ) {
        ( $X_i, my $slope1, my $slope2 ) =
          _interpolate_scalar( $Q, $dZmX_dX_b, $X_b, 0, $dZmX_dX_e, $X_e, 0,
            1 );
    }
    elsif ( $icol == I_E1 ) {
        ( $X_i, my $dX_dE1_i, my $d2X_dE12_i ) =
          _interpolate_scalar( $Q, $E1_b, $X_b, 1 / $dE1_dX_b,
            $E1_e, $X_e, 1 / $dE1_dX_e, $interp );
    }
    elsif ( $icol == I_dE1dX ) {
        ( $X_i, my $slope1, my $slope2 ) =
          _interpolate_scalar( $Q, $dE1_dX_b, $X_b, 0, $dE1_dX_e, $X_e, 0, 1 );
    }
    elsif ( $icol == I_Er ) {
        ( $X_i, my $dX_dE_i, my $d2X_dE2_i ) =
          _interpolate_scalar( $Q, $E_b, $X_b, 1 / $dE_dX_b,
            $E_e, $X_e, 1 / $dE_dX_e, $interp );
    }
    elsif ( $icol == I_dErdX ) {
        ( $X_i, my $slope1, my $slope2 ) =
          _interpolate_scalar( $Q, $dE_dX_b, $X_b, 0, $dE_dX_e, $X_e, 0, 1 );
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
    my ( $T_i, $dT_dX_i, $d2T_dX2_i ) =
      _interpolate_scalar( $X_i, $X_b, $T_b, $dT_dX_b, $X_e, $T_e, $dT_dX_e,
        $interp );
    my ( $E1_i, $dE1_dX_i, $d2E1_dX2_i ) =
      _interpolate_scalar( $X_i, $X_b, $E1_b, $dE1_dX_b, $X_e, $E1_e, $dE1_dX_e,
        $interp );
    my ( $E_i, $dE_dX_i, $d2E_dX2_i ) =
      _interpolate_scalar( $X_i, $X_b, $E_b, $dE_dX_b, $X_e, $E_e, $dE_dX_e,
        $interp );

    my @vars;
    @vars[
      I_X,   I_Y,      I_dYdX, I_Z,     I_dZdX, I_T, I_dTdX,
      I_ZmX, I_dZmXdX, I_E1,   I_dE1dX, I_Er,   I_dErdX
      ]
      = (
        $X_i,      $Y_i,     $dY_dX_i,    $Z_i,         $dZ_dX_i,
        $T_i,      $dT_dX_i, $Z_i - $X_i, $dZ_dX_i - 1, $E1_i,
        $dE1_dX_i, $E_i,     $dE_dX_i
      );
    return [@vars];

    #    return [
    #        $X_i,      $Y_i,     $dY_dX_i,    $Z_i,         $dZ_dX_i,
    #        $T_i,      $dT_dX_i, $Z_i - $X_i, $dZ_dX_i - 1, $E1_i,
    #        $dE1_dX_i, $E_i,     $dE_dX_i
    #    ];
}

sub _interpolate_scalar {
    my ( $xx, $x1, $y1, $dydx1, $x2, $y2, $dydx2, $interp ) = @_;
    my ( $yy, $dydx, $d2ydx2, $d3ydx3 );

    # User tables might not define all variables, such as E
    if ( defined($x1) && defined($y1) && defined($x2) && defined($y2) ) {

        # drop down to linear intepolation if slopes not defined
        if ( !defined($dydx1) || !defined($dydx2) ) { $interp = 1 }

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
    }
    return ( $yy, $dydx, $d2ydx2, $d3ydx3 );
}

sub _linear_interpolation {
    my ( $xx, $x1, $y1, $x2, $y2 ) = @_;
    my $dx21 = $x2 - $x1;
    my $dy21 = $y2 - $y1;
    my $dydx = ( $dx21 == 0 ) ? 0 : $dy21 / $dx21;
    my $yy   = $y1 + $dydx * ( $xx - $x1 );
    return wantarray ? ( $yy, $dydx ) : $yy;
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

    # A small tolerance to avoid interpolations
    my $eps = 1.e-6;

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

    # Define N consecutive lagrange interpolation points;
    # Using 4 points gives sufficient accuracy
    my $NLAG = 4;
    my $rj_interp = set_interpolation_points( $jl, $ntab, $NLAG );

    my ( $rx, $ry );

    # alpha varies approximately as 1/{ (gamma-1)*sqrt(gamma+1) },
    # so we can improve accuracy by interpolating the function
    # alpha*(gamma-1)*sqrt(gamma+1)
    foreach my $jj ( @{$rj_interp} ) {
        my ( $xx, $yy ) = @{ $rtab->[$jj] };
        push @{$rx}, $xx;
        push @{$ry}, $yy * ( $xx - 1 ) * sqrt( $xx + 1 );
    }

    my $ff = polint( $gamma, $rx, $ry );
    my $alpha = $ff / ( ( $gamma - 1 ) * sqrt( $gamma + 1 ) );

    return ($alpha);
}

sub is_monotonic_list {

    # returns 1 if list is monotonic increasing
    # returns -1 if list is monotonic decreasing
    # returns 0 if list is non-monotonic

    my ($rx) = @_;
    my $mono = 0;
    my $i_non_mono;    # index of non-monotonicity
    my $num = @{$rx};
    for ( my $i = 1 ; $i < $num ; $i++ ) {
        my $dF = $rx->[$i] - $rx->[ $i - 1 ];
        if ( $i == 1 ) { $mono = $dF > 0 ? 1 : $dF < 0 ? -1 : 0 }
        elsif ( $dF * $mono < 0 ) {
            $i_non_mono = $i;
            $mono       = 0;
            last;
        }
    }
    return wantarray ? ( $mono, $i_non_mono ) : $mono;
}

sub _make_intermediate_gamma_table {

    # Given: a symmetry, gamma, and a list of the indexes of the builtin tables
    #     to be interpolated
    # Return: a table of X,Y,dYdX,Z,dZdX for an intermediate gamma
    #     by using cubic interpolation between the four builtin tables

    my ( $symmetry, $gamma_x, $rilist_6, $rilist_4 ) = @_;
    $rilist_4 = $rilist_6 unless defined($rilist_4);
    my $rtable_new;
    my $alpha_x = alpha_interpolate( $symmetry, $gamma_x );
    my $A_x = -log( ( $gamma_x + 1 ) * $alpha_x );

    my $rtable_closest;
    my $dA_min;
    my $rlag_points_6;
    my $rlag_points_4;

    my $rA_6;
    foreach my $igam ( @{$rilist_6} ) {
        my ( $gamma, $table_name ) = @{ $rgamma_table->[$symmetry]->[$igam] };
        my $alpha = alpha_interpolate( $symmetry, $gamma );
        my ( $rtable, $rimpulse_table, $rtail_shock_table ) =
          get_builtin_tables($table_name);

        my $A = -log( ( $gamma + 1 ) * $alpha );
        my $dA = ( $A_x - $A );
        if ( !defined($dA_min) || abs($dA) < $dA_min ) {
            $rtable_closest = $rtable;
        }
        push @{$rlag_points_6}, [ $rtable, $A, $gamma, $alpha, $table_name ];
        push @{$rA_6}, $A;
    }

    my $rA_4;
    foreach my $igam ( @{$rilist_4} ) {
        my ( $gamma, $table_name ) = @{ $rgamma_table->[$symmetry]->[$igam] };
        my $alpha = alpha_interpolate( $symmetry, $gamma );
        ##my $rtable = get_builtin_table($table_name);
        my ( $rtable, $rimpulse_table, $rtail_shock_table ) =
          get_builtin_tables($table_name);
        my $A = -log( ( $gamma + 1 ) * $alpha );
        push @{$rlag_points_4}, [ $rtable, $A, $gamma, $alpha, $table_name ];
        push @{$rA_4}, $A;
    }

    # check that A values vary monotonically;
    # otherwise interpolation will fail
    # NOTE: assuming that rA_4 is monotonic if rA_6 is
    if ( !is_monotonic_list($rA_6) ) {

        print STDERR <<EOM;
Program error: non-monotonic A
symmetry=$symmetry gamma = $gamma_x alpha=$alpha_x 
indexes of tables = @{$rilist_6}
A values are (@{$rA_6})
EOM
        return;
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

    # find the minimum Ymax and maximum Ymin of the collection
    my ( $Ymin, $Ymax );
    foreach my $rpoint ( @{$rlag_points_6}, @{$rlag_points_4} ) {
        my ( $rtable, $A ) = @{$rpoint};
        my $Ymax_i = $rtable->[0]->[1];
        my $Ymin_i = $rtable->[-1]->[1];
        if ( !defined($Ymax) || $Ymax_i < $Ymax ) { $Ymax = $Ymax_i }
        if ( !defined($Ymin) || $Ymin_i > $Ymin ) { $Ymin = $Ymin_i }
    }

    # add a small tolerance for floating point comparisons near table edges
    my $eps = 1.e-8;
    $Ymax -= $eps;
    $Ymin += $eps;

    # Set the list of Y points equal to the closest table
    my @Ylist = ($Ymax);
    foreach my $item_ref ( @{$rtable_closest} ) {
        my $YY = $item_ref->[1];
        push @Ylist, $YY if ( $YY > $Ymin && $YY < $Ymax );
    }
    push @Ylist, $Ymin;

    # We are interpolating X and Z with 6 interpolation points.
    # We are interpolating dYdX and dZdX with 4 interpolation points.
    my $rtab_x = [];

    foreach my $Y (@Ylist) {
        my ( $rX, $rdYdX, $rZ, $rdZdX, $rE1, $rdE1dX, $rE, $rdEdX );

        my $missing_item;
        foreach my $rpoint ( @{$rlag_points_6} ) {
            my ( $rtable, $A, $gamma ) = @{$rpoint};
            my $item = $lookup->( $Y, $rtable );
            if ( !defined($item) ) { $missing_item = 1; last }
            my (
                $X,   $YY,     $dYdX, $Z,     $dZdX, $T, $dTdX,
                $ZmX, $dZmXdX, $E1,   $dE1dX, $E,    $dEdX
              )
              = @{$item}[
              I_X,    I_Y,   I_dYdX,   I_Z,  I_dZdX,  I_T,
              I_dTdX, I_ZmX, I_dZmXdX, I_E1, I_dE1dX, I_Er,
              I_dErdX
              ];
            push @{$rX},  $X;
            push @{$rZ},  $Z;
            push @{$rE1}, $E1;
            push @{$rE},  $E;
        }

        # All values must exist in order to interpolate
        # (tables all start and end at slightly different values)
        next if ($missing_item);

        foreach my $rpoint ( @{$rlag_points_4} ) {
            my ( $rtable, $A, $gamma ) = @{$rpoint};
            my $item = $lookup->( $Y, $rtable );
            if ( !defined($item) ) { $missing_item = 1; last }

            my (
                $X,   $YY,     $dYdX, $Z,     $dZdX, $T, $dTdX,
                $ZmX, $dZmXdX, $E1,   $dE1dX, $E,    $dEdX
              )
              = @{$item}[
              I_X,    I_Y,   I_dYdX,   I_Z,  I_dZdX,  I_T,
              I_dTdX, I_ZmX, I_dZmXdX, I_E1, I_dE1dX, I_Er,
              I_dErdX
              ];
            push @{$rdYdX},  $dYdX;
            push @{$rdZdX},  $dZdX;
            push @{$rdE1dX}, $dE1dX;
            push @{$rdEdX},  $dEdX;
        }
        next if ($missing_item);

        my $X_x     = polint( $A_x, $rA_6, $rX );
        my $Z_x     = polint( $A_x, $rA_6, $rZ );
        my $E1_x    = polint( $A_x, $rA_6, $rE1 );
        my $E_x     = polint( $A_x, $rA_6, $rE );
        my $dYdX_x  = polint( $A_x, $rA_4, $rdYdX );
        my $dZdX_x  = polint( $A_x, $rA_4, $rdZdX );
        my $dE1dX_x = polint( $A_x, $rA_4, $rdE1dX );
        my $dEdX_x  = polint( $A_x, $rA_4, $rdEdX );

        # leave spaces for T and dTdX, Z-X and slope
        # FIXME: the four unused vars can now be removed

        my @vars;
        @vars[
          I_X,   I_Y,      I_dYdX, I_Z,     I_dZdX, I_T, I_dTdX,
          I_ZmX, I_dZmXdX, I_E1,   I_dE1dX, I_Er,   I_dErdX
          ]
          = (
            $X_x, $Y, $dYdX_x, $Z_x,  $dZdX_x,  0,
            0,    0,  0,       $E1_x, $dE1dX_x, $E_x,
            $dEdX_x
          );

        push @{$rtab_x}, [@vars];
##?          [
##?            $X_x, $Y, $dYdX_x, $Z_x,  $dZdX_x,  0,
##?            0,    0,  0,       $E1_x, $dE1dX_x, $E_x,
##?            $dEdX_x
##?          ];
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
