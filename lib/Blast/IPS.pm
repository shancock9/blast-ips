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
use Storable qw(dclone);

use Blast::IPS::AlphaTable qw(alpha_interpolate);
use Blast::IPS::MathUtils qw(
  locate_2d
  polint
  set_interpolation_points
  table_row_interpolation
);
use Blast::IPS::Medium;
use Blast::IPS::Utils qw( check_keys );
use Blast::IPS::Data;

my $rgamma_table = Blast::IPS::Data->get_index();
my %valid_input_keys;

BEGIN {

    # Shock Table variable layout scheme
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
        I_ImX      => $i++,
        I_dImXdX   => $i++,
    };

    # Allowed keys for queries to sub wavefront
    %valid_input_keys = (
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
        'I-X'              => I_ImX,       # total shock residual energy
        'dI-XdX'           => I_dImXdX,
        interpolation_flag => 0,
    );

    # Impulse table variable layout scheme
    my $j = 0;
    use constant {
        J_X           => $j++,
        J_Y           => $j++,
        J_Z           => $j++,
        J_rpint_pos   => $j++,
        J_rpint_neg   => $j++,
        J_z_pose_rs   => $j++,
        J_z_nege_rs   => $j++,
        J_Qint_pos    => $j++,
        J_rovp_min    => $j++,
        J_z_pmin_rs   => $j++,
        J_ke_pos      => $j++,
        J_work_pos    => $j++,
        J_Disp_pos    => $j++,
    };

    # Tail shock table variable layout scheme
    my $k = 0;
    use constant {
        K_T           => $k++,
        K_z           => $k++,
        K_S1          => $k++,
        K_S2          => $k++,
    };

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
      file
      hide
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
    check_keys( $rinput_hash, \%valid_input_keys,
        "Checking for valid input keys" );

    # The following four quantities can be specified:
    my $table_name = $rinput_hash->{table_name};
    my $gamma      = $rinput_hash->{gamma};
    my $file       = $rinput_hash->{file};
    my $hide       = $rinput_hash->{hide};

    # Convert symmetry to numeric value;
    # allow alphanumeric abbreviations for symmetry
    # i.e. 'S', 'sph', 'spherical' etc
    my $symmetry   = $rinput_hash->{symmetry};
    if ( defined($symmetry) ) {
        if    ( $symmetry =~ /^S/i ) { $symmetry = 2 }
        elsif ( $symmetry =~ /^C/i ) { $symmetry = 1 }
        elsif ( $symmetry =~ /^P/i ) { $symmetry = 0 }
	$rinput_hash->{symmetry} = $symmetry;
    }

    # Look for a table
    my $rTables = get_builtin_tables(%{$rinput_hash});

    if ( !defined($rTables) ) {

        # No builtin table found, so make interpolated tables
        $rTables = _make_interpolated_tables( $rinput_hash );

        if ( !$rTables ) {
            carp(<<EOM);
Unable to create an interpolated table for symmetry=$symmetry and gamma=$gamma
EOM
            return;
        }
    }

    # Now check the results
    my $symmetry_old = $symmetry;
    my $gamma_old    = $gamma;

    $table_name = $rTables->{table_name};
    $gamma      = $rTables->{gamma};
    $symmetry   = $rTables->{symmetry};

    # Be sure requested symmetry and gamma agree with the values in the table
    if ( defined($symmetry_old) && $symmetry_old != $symmetry ) {
                carp(<<EOM);
symmetry of builtin table name: '$table_name' is '$symmetry' which differs from
your value '$symmetry_old'.  You may have entered the wrong value for this table.
EOM
                return;
    }
    if ( defined($gamma_old) && $gamma_old != $gamma ) {
                carp(<<EOM);
gamma of builtin table name: '$table_name' is '$gamma' which differs from
your value '$gamma_old'.  You may have entered the wrong value for this table.
EOM
                return;
    }

    my $error = "";
    if ( $gamma <= 1.0 ) { $error .= "Bad gamma='$gamma'\n" }
    if ( $symmetry != 0 && $symmetry != 1 && $symmetry != 2 ) {
        $error .= "Bad symmetry='$symmetry'\n";
    }

    my $rshock_table = $rTables->{shock_table};
    my $table_error  = _check_table($rshock_table);
    $error .= $table_error;

    my $num_table;
    if ( !$error ) {
        _update_toa_table($rshock_table);
        $num_table = @{$rshock_table};
    }

    my $rimpulse_table = $rTables->{impulse_table};
    my $rztables       = _make_z_tables($rimpulse_table);
    $rTables->{rztables} = $rztables;

    my $medium =
      Blast::IPS::Medium->new( symmetry => $symmetry, gamma => $gamma );

    # some basic parameters
    $self->{_gamma}      = $gamma;
    $self->{_symmetry}   = $symmetry;
    $self->{_rTables}    = $rTables;
    $self->{_table_name} = $table_name;
    $self->{_medium}     = $medium;

    # saved table lookup locations for shock table
    $self->{_jl} = -1;
    $self->{_ju} = $num_table;

    # saved table lookup locations for impulse table
    $self->{_jl2} = -1;
    $self->{_ju2} = undef;

    if ( !$error ) {


        # FIXME: make both available - put the first one in blast_info
        # Compute alpha. We can either get alpha from the shock table:
        # my $alpha = $self->alpha_from_shock_table();

        # or use the pre-computed values:  
        my $alpha = alpha_interpolate( $symmetry, $gamma );

	# The pre-computed values are very slightly more accurate, but the
	# difference is very small
        $self->{_alpha} = $alpha;

        # Compute parameters for extrapolating spherical waves
        my ( $A_far, $B_far, $Z_zero, $msg ) =
          long_range_parameters( $symmetry, $gamma, $rTables );
        $self->{_A_far}  = $A_far;
        $self->{_B_far}  = $B_far;
        $self->{_Z_zero} = $Z_zero;
        $error .= $msg;

    }

    if ($error) {
        $self->{_error} = $error;
        carp "$error\n";
    }
    return;
}

sub alpha_from_shock_table {
    my ($self)       = @_;
    my $rTables      = $self->{_rTables};
    my $rshock_table = $rTables->{shock_table};
    my $symmetry     = $rTables->{symmetry};
    my $gamma        = $rTables->{gamma};

    # we find the value of alpha at the first table point
    # based on R^(symmetry+1)*u**2=constant
    my ( $X_near, $Y_near, $dYdX_near, $Z_near, $dZdX_near ) =
      @{ $rshock_table->[0] }[ I_X, I_Y, I_dYdX, I_Z, I_dZdX ];
    my $lambda = exp($X_near);
    my $Prat   = exp($Y_near) + 1;
    my $q      = 2 * $gamma / ( ( $gamma + 1 ) * $Prat + ( $gamma - 1 ) );
    my $uovcsq = ( 2 / ( $gamma + 1 ) * ( 1 - $q ) )**2 / $q;
    my $C      = $lambda**( $symmetry + 1 ) * $uovcsq;
    my $alpha =
      ( 4 / ( ( $gamma + 1 ) * ( $symmetry + 3 ) ) )**2 / ( $gamma * $C );
    return $alpha;
}

sub long_range_parameters {
    my ( $symmetry, $gamma, $rTables ) = @_;
    my $rtable = $rTables->{shock_table};

    # Set asymptotic wave parameters for the distant region, given a table and gamma
    my ( $X_far, $Y_far, $dYdX_far, $Z_far, $dZdX_far ) = 
      @{ $rtable->[-1] }[ I_X, I_Y, I_dYdX, I_Z, I_dZdX ];
    my ( $A_far, $B_far, $Z_zero, $msg ) = ( 0, 0, 0, "" );
    if ( $symmetry == 2 ) {
        my $mu = -( $dYdX_far + 1 );
        if ( $mu > 0 ) {
            $A_far = exp($X_far) * exp($Y_far) / ( $gamma * sqrt( 2 * $mu ) );
            $B_far = 0.5 / $mu - $X_far;
        }
        else { $msg = "ending dYdX=$dYdX_far is bad\n" }
        $Z_zero =
          exp($Z_far) - 0.5 * ( $gamma + 1 ) * $A_far * sqrt( $X_far + $B_far );
        $Z_zero = log($Z_zero);
    }
    return ( $A_far, $B_far, $Z_zero, $msg );
}

sub get_index {
    return $rgamma_table;
}

sub get_builtin_tables {
    my @args = @_; 
    my $rTables             = Blast::IPS::Data::get(@args);
    return $rTables;
}

sub load_data_tables {
    my ( $rtab, @gammas ) = @_;

    # Create a cache of tables with different gamma values, for interpolation
    my $rTables_loaded;
    foreach my $igam (@gammas) {
        my ( $gamma, $table_name ) = @{ $rtab->[$igam] };
        if ( !defined( $rTables_loaded->{$table_name} ) ) {
            $rTables_loaded->{$table_name} = get_builtin_tables($table_name);
        }
    }
    return $rTables_loaded;
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
    my $_rinfo = $self->{_rTables}->{blast_info};
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

##sub check_tables {
##
##    # Should be called at program installation to check the tables
##
##    my @info_keys  = keys %{$rshock_tables_info};
##    my @table_keys = keys %{$rshock_tables};
##
##    # Check that the keys of the two table hashes are the same
##    my @missing_info_keys =
##      grep { !exists $rshock_tables_info->{$_} } @table_keys;
##    my @missing_table_keys =
##      grep { !exists $rshock_tables->{$_} } @info_keys;
##    my $error = @missing_info_keys || @missing_table_keys;
##    if ($error) {
##        local $" = ')(';
##        $error = <<EOM;
##------------------------------------------------------------------------
##Program error detected checking hash keys
##Have info but not table for: (@missing_table_keys)
##Have table but no info for: (@missing_info_keys)
##------------------------------------------------------------------------
##EOM
##        croak $error;
##    }
##
##    return $error if ($error);
##
##    foreach my $key (@info_keys) {
##        my $rtable = $rshock_tables->{$key};
##
##        # Check for a naming error
##        my $item = $rshock_tables_info->{$key};
##        my ( $symmetry, $gamma, $err_est ) = @{$item};
##        my $letter = $symmetry eq 0 ? 'P' : $symmetry eq 1 ? 'C' : 'S';
##        my $str = $letter . $gamma;
##
##        #print STDERR "Checking if table $key similar to $str\n";
##        if ( $key !~ /^$str/i ) {
##            return "Table $key: expected to be like $str\n";
##        }
##
##        # Check for numerical problems
##        my $error = _check_table($rtable);
##        if ($error) {
##            return "Table $key:\n" . $error;
##        }
##    }
##    return;
##}

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
        $error .= "Y is not monotonic decreasing\n";
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

    # E1 must be increasing if it is defined
    if ( defined($rtable->[0]->[I_E1]) && !$is_monotonic->(I_E1) > 0 ) {
        $error .= "E1 is not monotonic increasing\n";
    }

    # Er must be increasing if it is defined
    if ( defined($rtable->[0]->[I_Er]) && !$is_monotonic->(I_Er) > 0 ) {
        $error .= "Er is not monotonic increasing\n";
    }

    # we must have 0 <= E1 <= Er <=1
    if ( !$error
        && defined( $rtable->[0]->[I_E1] && defined( $rtable->[0]->[I_Er] ) ) )
    {
        foreach my $row ( @{$rtable} ) {
            my $E1 = $row->[ I_E1 ];
            my $Er = $row->[ I_Er ];
            my $tol = 1.e-6;
            if ( $E1 < 0 || ($E1 > $Er + $tol) || $Er < 0 || $Er > 1 ) {
                $error .= "Table violates 0<= E1=$E1 <= Er=$Er <=1";
                last;
            }
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
        my $T_pose_rs = $t_pose_rs > 0 ? log($t_pose_rs) : 0;
        my $T_pmin_rs = $t_pmin_rs > 0 ? log($t_pmin_rs) : 0;
        my $T_nege_rs = $t_nege_rs > 0 ? log($t_nege_rs) : 0;
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
    my $ntab = @{$rztab};
    my $rrow;
    my $icx         = 0;
    my $NLAG        = 4;
    my $no_extrap_l = 1;
    my $no_extrap_u = 1;
    my $zz;

    ( $rrow, $jl, $ju ) =
      table_row_interpolation( $T, $rztab, $icx, $NLAG, $jl, $ju,
        $no_extrap_l, $no_extrap_u, );

    ##($jl, $ju) = locate_2d($T, 0, $rztab, $jl, $ju);
    if ( $jl < 0 ) {

        # before start of table;
        # assume pos phase ends at origin (r=0), so
        my $tt = exp($T);
        $zz = -$tt;
    }
    elsif ( $ju >= $ntab ) {

        # off end of table; use supplied default
        $zz = $z_default;
    }
    else {
        #my $rrow = _interpolate_table_rows( $T, $rztab, $jl );
        #if ($rrow) {
        $zz = $rrow->[1];
        ##}
    }
    return $zz;
}

sub get_impulse {
    my ( $self, $result ) = @_;
    my $rimpulse_table = $self->{_rTables}->{impulse_table};
    return unless ($rimpulse_table);
    my $jl       = $self->{_jl2};
    my $ju       = $self->{_ju2};
    my $symmetry = $self->{_symmetry};
    my $gamma    = $self->{_gamma};
    my $rreturn_hash;
    my ( $X, $Y, $dYdX, $Z, $dZdX ) =
      @{$result}[ I_X, I_Y, I_dYdX, I_Z, I_dZdX ];

    # Locate this point in the table
    my $ntab = @{$rimpulse_table};
    my $rrow;
    my $icx         = 0;
    my $NLAG        = 4;
    my $no_extrap_l = 1;
    my $no_extrap_u = 1;

    ( $rrow, $jl, $ju ) =
      table_row_interpolation( $X, $rimpulse_table, $icx, $NLAG, $jl, $ju,
        $no_extrap_l, $no_extrap_u, );
    $self->{_jl2} = $jl;
    $self->{_ju2} = $ju;
    return unless ( defined($jl) );

    my (
        $rpint_pos, $rpint_neg, $z_pose_rs, $z_nege_rs, $Qint_pos,
        $rovp_min,  $z_pmin_rs, $ke_pos,    $work_pos,  $Disp_pos,
    );

    # Handle case before start of the table
    if ( $jl < 0 ) {

        $rrow = $rimpulse_table->[0]; # NOTE: Should not be needed; already done
        my (
            $X_b,         $Y_b,         $Z_b,         $rpint_pos_b,
            $rpint_neg_b, $z_pose_rs_b, $z_nege_rs_b, $Qint_pos_b,
            $rovp_min_b,  $z_pmin_rs_b, $ke_pos_b,    $work_pos_b,
            $Disp_pos_b,
            ##) = @{$rrow};
          ) = @{$rrow}[
          J_X,         J_Y,         J_Z,         J_rpint_pos,
          J_rpint_neg, J_z_pose_rs, J_z_nege_rs, J_Qint_pos,
          J_rovp_min,  J_z_pmin_rs, J_ke_pos,    J_work_pos,
          J_Disp_pos,
          ];

        # ke_pos becomes constant (the KE of the similarity solution)
        $ke_pos = $ke_pos_b;

        # positive work equals total work since there is no negative work
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

        # We will use log extrapolation for the following variables back to the
        # origin.  Note, in spherical symmetry, the theoretical positive phase
        # impulse varies like 1/sqrt(r) at r=0 (see charts in Sedov p321), but
        # the range over which this applies is vanishingly small so it is not
        # possible to make use of this.
        $Disp_pos  = $Disp_pos_b;
        $Qint_pos  = $Qint_pos_b;
        $rpint_pos = $rpint_pos_b;

        # Here are the theoretical slopes for dQint_pos/dX from the similarity
        # solution.
        # The table values are very close to these.
        # -0.5 for spherical symmetry
        #    0 for cylindrical symmetry
        # +0.5 for plane symmetry

        # Positive phase impulse and displacement are harder to do
        # theoretically because they end when absolute pressure hits 1, so the
        # integral does not go to infinity.  So it seems best to just use log
        # slopes of the computed values to continue the solution before the
        # first table point.

        if ( @{$rimpulse_table} > 2 ) {

            # Get the next row so that we can calculate the slope
            my $rrow_a = $rimpulse_table->[1];

	    # FIXME: only need to pull out a few of these:
            my (
                $X_a,         $Y_a,         $Z_a,         $rpint_pos_a,
                $rpint_neg_a, $z_pose_rs_a, $z_nege_rs_a, $Qint_pos_a,
                $rovp_min_a,  $z_pmin_rs_a, $ke_pos_a,    $work_pos_a,
                $Disp_pos_a,
                ##) = @{$rrow_a};
              ) = @{$rrow_a}[
              J_X,         J_Y,         J_Z,         J_rpint_pos,
              J_rpint_neg, J_z_pose_rs, J_z_nege_rs, J_Qint_pos,
              J_rovp_min,  J_z_pmin_rs, J_ke_pos,    J_work_pos,
              J_Disp_pos,
              ];

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

        $rrow = $rimpulse_table->[-1];
        my (
            $X_e,         $Y_e,         $Z_e,         $rpint_pos_e,
            $rpint_neg_e, $z_pose_rs_e, $z_nege_rs_e, $Qint_pos_e,
            $rovp_min_e,  $z_pmin_rs_e, $ke_pos_e,    $work_pos_e,
            $Disp_pos_e,
            ##) = @{$rrow};
          ) = @{$rrow}[
          J_X,         J_Y,         J_Z,         J_rpint_pos,
          J_rpint_neg, J_z_pose_rs, J_z_nege_rs, J_Qint_pos,
          J_rovp_min,  J_z_pmin_rs, J_ke_pos,    J_work_pos,
          J_Disp_pos,
          ];

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
        #my $rrow = _interpolate_table_rows( $X, $rimpulse_table, $jl );
        #if ($rrow) {
        (
            my $X_i,    my $Y_i,    my $Z_i,   $rpint_pos, $rpint_neg,
            $z_pose_rs, $z_nege_rs, $Qint_pos, $rovp_min,  $z_pmin_rs,
            $ke_pos,    $work_pos,  $Disp_pos
              ##) = @{$rrow};
          ) = @{$rrow}[
          J_X,         J_Y,         J_Z,         J_rpint_pos,
          J_rpint_neg, J_z_pose_rs, J_z_nege_rs, J_Qint_pos,
          J_rovp_min,  J_z_pmin_rs, J_ke_pos,    J_work_pos,
          J_Disp_pos,
          ];

        #}
    }
    my $disp_pos = defined($Disp_pos) ? exp($Disp_pos) : undef;
    my $qint_pos = defined($Qint_pos) ? exp($Qint_pos) : undef;

    # Switching to impulse instead of impulse time r^(j/2)
    my $r        = exp($X);
    my $rpow     = $r**( $symmetry / 2 );
    my $pint_pos = $rpint_pos / $rpow;
    my $pint_neg = $rpint_neg / $rpow;
    $rreturn_hash->{pint_pos}  = $pint_pos;
    $rreturn_hash->{pint_neg}  = $pint_neg;
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

sub get_ASYM       { my $self = shift; return $self->{_symmetry}; }
sub get_symmetry   { my $self = shift; return $self->{_symmetry} }
sub get_gamma      { my $self = shift; return $self->{_gamma} }
sub get_alpha      { my $self = shift; return $self->{_alpha} }
sub get_error      { my $self = shift; return $self->{_error} }
sub get_table      { my $self = shift; return $self->{_rTables}->{shock_table} }
sub get_Tables     { my $self = shift; return $self->{_rTables} }
sub get_table_name { my $self = shift; return $self->{_table_name} }

sub clone {

    my ($self) = @_;
    my $class = ref($self);
    my $newobj = bless { %{$self} }, $class;
    my $rTables = $newobj->{_rTables};
    $newobj->{_rTables} = dclone($rTables);
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
    my $symmetry = $self->{_symmetry};
    my $gamma    = $self->{_gamma};
    my $rTables  = $self->{_rTables};
    my $Sint_pos = $rTables->{blast_info}->{Sintegral_pos};
    my $Sint_neg = $rTables->{blast_info}->{Sintegral_neg};

    my $rtab = $rTables->{shock_table};
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

    # Validate input keys
    check_keys( $rinput_hash, \%valid_input_keys,
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
    #my ( $il, $iu ) = $self->_locate_2d( $Q, $icol );
    my $il = $self->{_jl};
    my $iu = $self->{_ju};
    ( $il, $iu ) = locate_2d( $Q, $icol, $rtab, $il, $iu );
    $self->{_jl} = $il;
    $self->{_ju} = $iu;

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
          _interpolate_shock_table_rows( $Q, $icol, $rtab->[$il], $rtab->[$iu], $interp );
    }

    my (
        $X,   $Y,      $dYdX, $Z,     $dZdX, $T, $dTdX,
        $ZmX, $dZmXdX, $E1,   $dE1dX, $Er,   $dErdX, $ImX, $dImXdX,
      )
      = @{$result}[
      I_X,      I_Y,  I_dYdX,  I_Z,  I_dZdX, I_T, I_dTdX, I_ZmX,
      I_dZmXdX, I_E1, I_dE1dX, I_Er, I_dErdX, I_ImX, I_dImXdX
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
        $z_nege_rs, $KE_pos,    $W_pos,    $disp_pos,  $qint_pos,
	$pint_pos, $pint_neg,
    );
    if ($rimpulse_hash) {
        $rpint_pos = $rimpulse_hash->{rpint_pos};
        $rpint_neg = $rimpulse_hash->{rpint_neg};
        $pint_pos = $rimpulse_hash->{pint_pos};
        $pint_neg = $rimpulse_hash->{pint_neg};
        $rovp_min  = $rimpulse_hash->{rovp_min};
        $z_pmin_rs = $rimpulse_hash->{z_pmin_rs};
        $z_pose_rs = $rimpulse_hash->{z_pose_rs};
        $z_nege_rs = $rimpulse_hash->{z_nege_rs};
        $KE_pos    = $rimpulse_hash->{KE_pos};
        $W_pos     = $rimpulse_hash->{W_pos};
        $disp_pos  = $rimpulse_hash->{disp_pos};
        $qint_pos  = $rimpulse_hash->{qint_pos};
    }
    elsif(defined($ImX)) {

        # get impulse with the I-X method if gamma interpolation makes
	# the impulse hash unavailable
        $pint_pos = exp( $ImX + $X );
        $rpint_pos = $pint_pos*exp($X)**($symmetry/2);
    }

    # FIXME: evaluate sigma and Sigma = sigma*r**(symmetry/2)

    # find z coordinates along a profile in space which includes this shock
    my ( $z_pose_ts, $z_pmin_ts, $z_nege_ts ) =
      ( $z_pose_rs, $z_pmin_rs, $z_nege_rs );

    my $rztables = $self->{_rTables}->{rztables};
    if ( defined($rztables) ) {
        my ( $rzpose_table, $rzpmin_table, $rznege_table ) = @{$rztables};
        $z_pose_ts = $self->get_z( $T, $rzpose_table, $z_pose_ts );
        $z_pmin_ts = $self->get_z( $T, $rzpmin_table, $z_pmin_ts );
        $z_nege_ts = $self->get_z( $T, $rznege_table, $z_nege_ts );
    }

    # backup analytical fit to zero pressure lines
    #my $rs = exp($X);
    #my $zs = exp($Z);
    #my ( $Tpos, $Lpos, $Tneg, $Lneg ) = $self->get_phase_lengths( $rs, $zs );

    # FUTURE: use more accurate values if possible
    #if ( defined($z_pose_rs) ) { $Tpos = $zs - $z_pose_rs; }
    #if ( defined($z_pose_ts) ) { $Lpos = $zs - $z_pose_ts; }

    # TableLoc  shows which table rows were interpolated
    # TableVars shows the interpolated row values

    # For convenience, compute the positive phase time duration and length 
    # since they are frequently needed.
    my $zs = exp($Z);
    my $Tpos = defined($z_pose_rs) ? $zs - $z_pose_rs : 0;
    my $Lpos = defined($z_pose_ts) ? $zs - $z_pose_ts : 0;

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
        'I-X'     => $ImX,
        'dI-XdX'  => $dImXdX,
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
        'pint_pos'    => $pint_pos,
        'pint_neg'    => $pint_neg,
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

    # add some shock front values
    my $medium=$self->{_medium};
    my $rsf_values =
      $medium->profile_slopes_from_dYdX( exp($X), exp($Y), $dYdX );
    foreach my $key(keys %{$rsf_values}) {
	$return_hash->{$key}=$rsf_values->{$key};
    }

    return $return_hash;
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

    my $rtable = $self->{_rTables}->{shock_table};
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
    my $rtab     = $self->{_rTables}->{shock_table};
    my $symmetry = $self->{_symmetry};
    my $gamma    = $self->{_gamma};
    my $alpha    = $self->{_alpha};
    my $delta    = ( 3 + $symmetry ) / 2;
    my $p_amb    = 1;                

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

    ## FIXME: assuming constant impulse at very close range; needs research
    my ( $ImX_0, $dImXdX_0 ) = @{ $rtab->[0] }[ I_ImX, I_dImXdX ];
    my $ImX_i    = $ImX_0;
    my $dImXdX_i = $dImXdX_0;
    ############

    # The result
    my $result_i = [
        $X_i,     $Y_i,    $dY_dX_i,    $Z_i,         $dZ_dX_i,
        $T_i,     $dTdX_i, $Z_i - $X_i, $dZ_dX_i - 1, $E1_i,
        $dE1dX_i, $Er_i,   $dErdX_i, $ImX_i, $dImXdX_i
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
    my $rtab     = $self->{_rTables}->{shock_table};
    my $symmetry = $self->{_symmetry};
    my $gamma    = $self->{_gamma};

    # Add the energy variables to the extrapolated result

    # current extrapolated state values
    my ( $X_i, $Y_i, $dYdX_i ) = @{$result_i}[ I_X, I_Y, I_dYdX ];

    # last table point is reference state
    my (
        $X_e,   $Y_e,      $dYdX_e, $Z_e,     $dZdX_e, $T_e, $dTdX_e,
        $ZmX_e, $dZmXdX_e, $E1_e,   $dE1dX_e, $Er_e,   $dErdX_e, $ImX_e, $dImXdX_e
      )
      = @{ $rtab->[-1] }[
      I_X,    I_Y,   I_dYdX,   I_Z,  I_dZdX, I_T,
      I_dTdX, I_ZmX, I_dZmXdX, I_E1, I_dE1dX, I_Er, I_dErdX, I_ImX, I_dImXdX,
      ];

    # Tentative initialize to table end in case work is zero
    my $E1_i    = $E1_e;
    my $Er_i    = $Er_e;
    my $dE1dX_i = $dE1dX_e;
    my $dErdX_i = $dErdX_e;
    my $dImXdX_i = $dImXdX_e;
    my $ImX_i = $ImX_e; 

    # In case user adds a table (testing only), we do not want to continue
    return unless ( defined($Er_e) && defined($ImX_e));

    $ImX_i = $ImX_e + $dImXdX_e * ( $X_i - $X_e );

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
    $result_i->[I_ImX]     = $ImX_i;
    $result_i->[I_dImXdX]  = $dImXdX_i;
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
    my $rtab     = $self->{_rTables}->{shock_table};

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
    my $rtab     = $self->{_rTables}->{shock_table};
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

        # Given T=ln(TOA). We have to iterate in this case because assuming
        # constant dT/dX is not sufficiently accurate.  Since Z hardly changes
        # with distance, an accurate first guess is made using the last Z in
        # the table. Simple iteration converges in just a couple of steps.
        $Z_i = $Z_e;
        $X_i = $X_e;
        foreach my $it ( 0 .. 5 ) {
            my $X_last = $X_i;
            $X_i = log( exp($Q) + exp($Z_i) );
            $Z_i = $Z_e + $dZ_dX_i * ( $X_i - $X_e );
            my $dX = $X_i - $X_last;
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

sub _interpolate_shock_table_rows {
    my ( $Q, $icol, $row_b, $row_e, $interp ) = @_;

    # Interpolate all of the variables between two table rows. Each row has
    # these values
    # 	[X,Y,Z,dYdX,dZdX,T,ZmX, dZmXdX, dTdX,E1, dE1dX, E, dEdX]

    my (
        $X_b,      $Y_b,     $dY_dX_b, $Z_b,       $dZ_dX_b,
        $T_b,      $dT_dX_b, $ZmX_b,   $dZmX_dX_b, $E1_b,
        $dE1_dX_b, $E_b,     $dE_dX_b, $ImX_b, $dImX_dX_b,
      )
      = @{$row_b}[
      I_X,      I_Y,  I_dYdX,  I_Z,  I_dZdX, I_T, I_dTdX, I_ZmX,
      I_dZmXdX, I_E1, I_dE1dX, I_Er, I_dErdX, I_ImX, I_dImXdX
      ];
    my (
        $X_e,      $Y_e,     $dY_dX_e, $Z_e,       $dZ_dX_e,
        $T_e,      $dT_dX_e, $ZmX_e,   $dZmX_dX_e, $E1_e,
        $dE1_dX_e, $E_e,     $dE_dX_e, $ImX_e, $dImX_dX_e,
      )
      = @{$row_e}[
      I_X,      I_Y,  I_dYdX,  I_Z,  I_dZdX, I_T, I_dTdX, I_ZmX,
      I_dZmXdX, I_E1, I_dE1dX, I_Er, I_dErdX, I_ImX, I_dImXdX,
      ];

    # FIXME: for slope interpolation, we are currently doing liner interp
    # We should either iterate or do parabolic interpolation using the
    # second derivatives at the segment ends

    # ZERO: If this is an inverse problem, first find X=X_i
    my $X_i;
    if ( $icol == I_X ) {
        $X_i = $Q;
    }

    # FIXME: these could be collapsed into two calls:
    # ODD: icol%2==1 : if icol is odd do a cubic interpolation of icol with
    # slope icol+1
    elsif ( $icol == I_Y ) {
        ( $X_i, my $dX_dY_i, my $d2X_dY2_i ) =
          _interpolate_scalar( $Q, $Y_b, $X_b, 1 / $dY_dX_b,
            $Y_e, $X_e, 1 / $dY_dX_e, $interp );
    }

    # EVEN: icol%2==0 : if icol is even do a slope interpolation
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
    elsif ( $icol == I_ImX ) {
        ( $X_i, my $dX_dI_i, my $d2X_dE2_i ) =
          _interpolate_scalar( $Q, $ImX_b, $X_b, 1 / $dImX_dX_b,
            $ImX_e, $X_e, 1 / $dImX_dX_e, $interp );
    }
    elsif ( $icol == I_dImXdX ) {
        ( $X_i, my $slope1, my $slope2 ) =
          _interpolate_scalar( $Q, $dImX_dX_b, $X_b, 0, $dImX_dX_e, $X_e, 0, 1 );
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
    my ( $ImX_i, $dImX_dX_i, $d2I_dX2_i ) =
      _interpolate_scalar( $X_i, $X_b, $ImX_b, $dImX_dX_b, $X_e, $ImX_e, $dImX_dX_e,
        $interp );

    my @vars;
    @vars[
      I_X,   I_Y,      I_dYdX, I_Z,     I_dZdX, I_T, I_dTdX,
      I_ZmX, I_dZmXdX, I_E1,   I_dE1dX, I_Er,   I_dErdX, I_ImX, I_dImXdX
      ]
      = (
        $X_i,      $Y_i,     $dY_dX_i,    $Z_i,         $dZ_dX_i,
        $T_i,      $dT_dX_i, $Z_i - $X_i, $dZ_dX_i - 1, $E1_i,
        $dE1_dX_i, $E_i,     $dE_dX_i, $ImX_i, $dImX_dX_i,
      );
    return [@vars];
}

sub _interpolate_scalar {
    my ( $xx, $x1, $y1, $dydx1, $x2, $y2, $dydx2, $linear ) = @_;
    my ( $yy, $dydx, $d2ydx2, $d3ydx3 );

    # Be sure all variables are defined
    if ( defined($x1) && defined($y1) && defined($x2) && defined($y2) ) {
	
        $yy     = $y1;
        $dydx   = $dydx1;
        $d2ydx2 = $d3ydx3 = 0;
        if ( $x1 != $x2 ) {

            # drop down to linear intepolation if slopes not defined
            if ( !defined($dydx1) || !defined($dydx2) ) { $linear = 1 }

            if ( defined($linear) && $linear == 1 ) {
                ( $yy, $dydx ) =
                  _linear_interpolation( $xx, $x1, $y1, $x2, $y2 );
                ( $dydx, $d2ydx2 ) =
                  _linear_interpolation( $xx, $x1, $dydx1, $x2, $dydx2 );
                $d3ydx3 = 0;
            }
            else {
                ( $yy, $dydx, $d2ydx2, $d3ydx3 ) =
                  _cubic_interpolation( $xx, $x1, $y1, $dydx1, $x2, $y2,
                    $dydx2 );
            }
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

    # Given the value of a function and its slope at two different points,
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

##############################################################################
# Routines to create new tables for intermediate gamma values by interpolation 
##############################################################################

sub _make_interpolated_tables {
    my ( $rinput_hash) =@_;

    my $table_name = $rinput_hash->{table_name};
    my $gamma      = $rinput_hash->{gamma};
    my $symmetry   = $rinput_hash->{symmetry};
    my $file       = $rinput_hash->{file};
    my $hide       = $rinput_hash->{hide};

    return unless(defined($symmetry) && defined($gamma));

    my $rgtab  = $rgamma_table->[$symmetry];

    # Check for the hidden table option, which is important for testing. In
    # this option we make it appear as if any particular table does not exist
    # in the data modules by removing it from the list of modules. In this way
    # we can force an object to be created with interpolation of the remaining
    # tables. This allows us to make a direct comparison of any data set with
    # data created by interpolating from the remaining data. This allows the
    # accuracy of the interpolation schemes to be assessed.  See 'hide.t'.
    if ($hide) {
        my @filtered = grep { $_->[1] ne $hide } @{$rgtab};
        $rgtab = \@filtered;
    }

    my $rTables;
    my $result = _gamma_lookup( $symmetry, $gamma, $rgtab );

    $table_name         = $result->{table_name};
    my $loc_gamma_table = $result->{loc_gamma_table};
    my $rigam_6         = $result->{rigam_6};
    my $rigam_4         = $result->{rigam_4};

    if ( defined($rigam_6) && defined($rigam_4) ) {

        # load all tables needed for interpolation
        my $rTables_loaded =
          load_data_tables( $rgtab, @{$rigam_6}, @{$rigam_4} );

        my $rtable = _make_interpolated_gamma_table( $symmetry, $gamma,
            $rigam_6, $rigam_4, $rTables_loaded, $rgtab );

        my $rblast_info =
          _make_interpolated_blast_info( $symmetry, $gamma, $rigam_4,
            $loc_gamma_table, $rTables_loaded, $rgtab );

        my $rimpulse_table =
          _make_interpolated_impulse_table( $symmetry, $gamma,
            $rigam_6, $rigam_4, $rTables_loaded, $rgtab );

        my $rtail_shock_table =
          _make_interpolated_tail_shock_table( $symmetry, $gamma,
            $rigam_6, $rigam_4, $rTables_loaded, $rgtab );

        if ( !$table_name ) {
            $table_name = Blast::IPS::Data::make_table_name( $symmetry, $gamma );
        }

        $rTables->{shock_table}      = $rtable;
        $rTables->{impulse_table}    = $rimpulse_table;
        $rTables->{tail_shock_table} = $rtail_shock_table;
        $rTables->{blast_info}       = $rblast_info;
        $rTables->{symmetry}         = $symmetry;
        $rTables->{gamma}            = $gamma;
        $rTables->{table_name}       = $table_name;
    }
    return $rTables;
}

sub _gamma_lookup {
    my ( $symmetry, $gamma, $rtab ) = @_;

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

    # if a gamma table is not given,
    # uses $rgamma_table, the global table of gamma values for the builtin
    # tables

    return if ( $symmetry != 0 && $symmetry != 1 && $symmetry != 2 );

    # The tolerance for testing gamma is currently fixed. It should be very
    # small because we can construct very accurate interpolated tables.
    my $eps = 1.e-7;

    my $return_hash = {};

    # allow some older calls to still work (see gamma_interp.t, which should
    # eventually be eliminated)
    if (!defined($rtab)) { $rtab  = $rgamma_table->[$symmetry]; }

    my $ntab  = @{$rtab};
    my $icol  = 0;

    my ($j2, $j3) = locate_2d($gamma, $icol, $rtab, undef,undef);

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

sub _make_interpolated_gamma_table {

    # Given: a symmetry, gamma, and a list of the indexes of the builtin tables
    #     to be interpolated
    # Return: a table of X,Y,dYdX,Z,dZdX for an intermediate gamma
    #     by using cubic interpolation between the four builtin tables

    my ( $symmetry, $gamma_x, $rilist_6, $rilist_4, $rTables_loaded, $rgtab ) = @_;
    $rilist_4 = $rilist_6 unless defined($rilist_4);

    # Allow older software to still work
    $rgtab  = $rgamma_table->[$symmetry] unless defined($rgtab);

    # Create a cache of needed tables if caller has not done so
    if ( !defined($rTables_loaded) ) {
	$rTables_loaded = load_data_tables($rgtab, @{$rilist_6}, @{$rilist_4});
    }

    my $rtable_new;
    my $alpha_x = alpha_interpolate( $symmetry, $gamma_x );
    my $A_x = -log( ( $gamma_x + 1 ) * $alpha_x );

    my $rtable_closest;
    my $dA_min;
    my $rlag_points_6;
    my $rlag_points_4;

    my $rA_6;
    foreach my $igam ( @{$rilist_6} ) {
        my ( $gamma, $table_name ) = @{ $rgtab->[$igam] };
        my $alpha = alpha_interpolate( $symmetry, $gamma );
        my $rTables = $rTables_loaded->{$table_name};
	my $rtable = $rTables->{shock_table};

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
        my ( $gamma, $table_name ) = @{ $rgtab->[$igam] };
        my $alpha = alpha_interpolate( $symmetry, $gamma );
        my $rTables = $rTables_loaded->{$table_name};
	my $rtable = $rTables->{shock_table};
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

        # Locate this point in the table
        my ( $il, $iu ) = locate_2d( $Y, $icol, $rtab, undef, undef );

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
            return _interpolate_shock_table_rows( $Y, $icol, $rtab->[$il],
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

    # We are interpolating vars, like X and Z, with 6 interpolation points.
    # We are interpolating slopes, like dYdX and dZdX, with 4 interpolation points.
    my $rtab_x = [];

    foreach my $Y (@Ylist) {
        my ( $rX, $rdYdX, $rZ, $rdZdX, $rE1, $rdE1dX, $rE, $rdEdX, $rImX, $rdImXdX );

        my $missing_item;
        foreach my $rpoint ( @{$rlag_points_6} ) {
            my ( $rtable, $A, $gamma ) = @{$rpoint};
            my $item = $lookup->( $Y, $rtable );
            if ( !defined($item) ) { $missing_item = 1; last }
            my (
                $X,   $YY,     $dYdX, $Z,     $dZdX, $T, $dTdX,
                $ZmX, $dZmXdX, $E1,   $dE1dX, $E,    $dEdX, $ImX, $dImXdX,
              )
              = @{$item}[
              I_X,    I_Y,   I_dYdX,   I_Z,  I_dZdX,  I_T,
              I_dTdX, I_ZmX, I_dZmXdX, I_E1, I_dE1dX, I_Er,
              I_dErdX, I_ImX, I_dImXdX,
              ];
            push @{$rX},  $X;
            push @{$rZ},  $Z;
            push @{$rE1}, $E1;
            push @{$rE},  $E;
            push @{$rImX},  $ImX;
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
                $ZmX, $dZmXdX, $E1,   $dE1dX, $E,    $dEdX, $ImX, $dImXdX,
              )
              = @{$item}[
              I_X,    I_Y,   I_dYdX,   I_Z,  I_dZdX,  I_T,
              I_dTdX, I_ZmX, I_dZmXdX, I_E1, I_dE1dX, I_Er,
              I_dErdX, I_ImX, I_dImXdX,
              ];
            push @{$rdYdX},  $dYdX;
            push @{$rdZdX},  $dZdX;
            push @{$rdE1dX}, $dE1dX;
            push @{$rdEdX},  $dEdX;
            push @{$rdImXdX},  $dImXdX;
        }
        next if ($missing_item);

        my $X_x      = polint( $A_x, $rA_6, $rX );
        my $Z_x      = polint( $A_x, $rA_6, $rZ );
        my $E1_x     = polint( $A_x, $rA_6, $rE1 );
        my $E_x      = polint( $A_x, $rA_6, $rE );
        my $ImX_x    = polint( $A_x, $rA_6, $rImX );
        my $dYdX_x   = polint( $A_x, $rA_4, $rdYdX );
        my $dZdX_x   = polint( $A_x, $rA_4, $rdZdX );
        my $dE1dX_x  = polint( $A_x, $rA_4, $rdE1dX );
        my $dEdX_x   = polint( $A_x, $rA_4, $rdEdX );
        my $dImXdX_x = polint( $A_x, $rA_4, $rdImXdX );

        # Patch to keep energies in correct order
        if ( $E_x < $E1_x ) { $E_x = $E1_x }

	# Note: do not need to return derived vars T, dTdX, Z-X and its slope
	# since they are added later.
        my @vars;
        @vars[ I_X, I_Y, I_dYdX, I_Z, I_dZdX, I_E1, I_dE1dX, I_Er, I_dErdX, I_ImX, I_dImXdX ] =
          ( $X_x, $Y, $dYdX_x, $Z_x, $dZdX_x, $E1_x, $dE1dX_x, $E_x, $dEdX_x, $ImX_x, $dImXdX_x );

        push @{$rtab_x}, [@vars];
    }
    return $rtab_x;
}

sub _make_interpolated_impulse_table {

    # Given: a symmetry, gamma, and a list of the indexes of the builtin tables
    #     to be interpolated
    # Return: an impulse table for an intermediate gamma made

    my ( $symmetry, $gamma_x, $rilist_6, $rilist_4, $rTables_loaded, $rgtab ) = @_;
    $rilist_4 = $rilist_6 unless defined($rilist_4);

    # patch to allow older software to still work:
    $rgtab  = $rgamma_table->[$symmetry] unless defined($rgtab);

    # Create a cache of needed tables if caller has not done so
    if ( !defined($rTables_loaded) ) {
        $rTables_loaded =
          load_data_tables( $symmetry, @{$rilist_6}, @{$rilist_4}, $rgtab );
    }

    my $rtable_new;
    my $alpha_x = alpha_interpolate( $symmetry, $gamma_x );
    my $A_x = -log( ( $gamma_x + 1 ) * $alpha_x );

    my $rtable_closest;
    my $dA_min;
    my $rlag_points_6;
    my $rlag_points_4;

    # Prepare for either 4 point or 6 point interpolations, even though we
    # might only use 4 points.
    my $rA_6;
    foreach my $igam ( @{$rilist_6} ) {
        my ( $gamma, $table_name ) = @{ $rgtab->[$igam] };
        my $alpha        = alpha_interpolate( $symmetry, $gamma );
        my $rTables      = $rTables_loaded->{$table_name};
        my $rtable       = $rTables->{impulse_table};
        my $rshock_table = $rTables->{shock_table};

        my $A = -log( ( $gamma + 1 ) * $alpha );
        my $dA = ( $A_x - $A );
        if ( !defined($dA_min) || abs($dA) < $dA_min ) {
            $rtable_closest = $rtable;
        }
        push @{$rlag_points_6},
          [ $rtable, $A, $gamma, $alpha, $rshock_table, $table_name ];
        push @{$rA_6}, $A;
    }

    my $rA_4;
    foreach my $igam ( @{$rilist_4} ) {
        my ( $gamma, $table_name ) = @{ $rgtab->[$igam] };
        my $alpha        = alpha_interpolate( $symmetry, $gamma );
        my $rTables      = $rTables_loaded->{$table_name};
        my $rtable       = $rTables->{impulse_table};
        my $rshock_table = $rTables->{shock_table};
        my $A            = -log( ( $gamma + 1 ) * $alpha );
        push @{$rlag_points_4},
          [ $rtable, $A, $gamma, $alpha, $rshock_table, $table_name ];
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

    # We are interpolating vars, like X and Z, with 6 interpolation points.  
    # We are interpolating slopes, like dYdX and dZdX, with 4 interpolation
    # points.
    my $rtab_x = [];

    foreach my $Y (@Ylist) {

        # arrays for holding tables of values
        my (
            $_X,         $_Y,         $_Z,        $_rpint_pos, $_rpint_neg,
            $_z_pose_rs, $_z_nege_rs, $_Qint_pos, $_rovp_min,  $_z_pmin_rs,
            $_ke_pos,    $_work_pos,  $_Disp_pos
        );

        my $missing_item;
        foreach my $rpoint ( @{$rlag_points_6} ) {
            my ( $rtable, $A, $gamma, $alpha, $rshock_table ) = @{$rpoint};
            my $item = table_row_interpolation( $Y, $rtable, 1, 4 );
            if ( !defined($item) ) { $missing_item = 1; last }
            my (
                $X,         $Y,         $Z,        $rpint_pos, $rpint_neg,
                $z_pose_rs, $z_nege_rs, $Qint_pos, $rovp_min,  $z_pmin_rs,
                $ke_pos,    $work_pos,  $Disp_pos
              )
              = @{$item}[
              J_X,         J_Y,         J_Z,         J_rpint_pos,
              J_rpint_neg, J_z_pose_rs, J_z_nege_rs, J_Qint_pos,
              J_rovp_min,  J_z_pmin_rs, J_ke_pos,    J_work_pos,
              J_Disp_pos,
              ];

            # X and Z are more accurately interpolated from the shock table for
            # this Y value.
            my ( $il, $iu ) = locate_2d( $Y, 1, $rshock_table );
            if ( $il >= 0 && $iu < @{$rshock_table} ) {
                my $result =
                  _interpolate_shock_table_rows( $Y, 1, $rshock_table->[$il],
                    $rshock_table->[$iu], 0 );

                my ( $X_old, $Z_old ) = ( $X, $Z );
                $X = $result->[I_X];
                $Z = $result->[I_Z];
            }

            # 6 point interpolation for these
            push @{$_X}, $X;
            push @{$_Z}, $Z;
        }

        # All values must exist in order to interpolate
        # (tables all start and end at slightly different values)
        next if ($missing_item);

        foreach my $rpoint ( @{$rlag_points_4} ) {
            my ( $rtable, $A, $gamma, $alpha, $rshock_table ) = @{$rpoint};
            my $item = table_row_interpolation( $Y, $rtable, 1, 4 );
            if ( !defined($item) ) { $missing_item = 1; last }

            my (
                $X,         $Y,         $Z,        $rpint_pos, $rpint_neg,
                $z_pose_rs, $z_nege_rs, $Qint_pos, $rovp_min,  $z_pmin_rs,
                $ke_pos,    $work_pos,  $Disp_pos
              )
              = @{$item}[
              J_X,         J_Y,         J_Z,         J_rpint_pos,
              J_rpint_neg, J_z_pose_rs, J_z_nege_rs, J_Qint_pos,
              J_rovp_min,  J_z_pmin_rs, J_ke_pos,    J_work_pos,
              J_Disp_pos,
              ];

            my $rs        = exp($X);
            my $t_pose_rs = $rs - $z_pose_rs;
            my $t_nege_rs = $rs - $z_nege_rs;
            my $t_pmin_rs = $rs - $z_pmin_rs;

            # 4 point interpolation variables go here:
            push @{$_rpint_pos}, $rpint_pos;
            push @{$_rpint_neg}, $rpint_neg;
            push @{$_z_pose_rs}, $t_pose_rs;    #$z_pose_rs;
            push @{$_z_nege_rs}, $t_nege_rs;    #$z_nege_rs;
            push @{$_Qint_pos},  $Qint_pos;
            push @{$_rovp_min},  $rovp_min;
            push @{$_z_pmin_rs}, $t_pmin_rs;    #$z_pmin_rs;
            push @{$_ke_pos},    $ke_pos;
            push @{$_work_pos},  $work_pos;
            push @{$_Disp_pos},  $Disp_pos;

        }
        next if ($missing_item);

        my $X_x         = polint( $A_x, $rA_6, $_X );
        my $Z_x         = polint( $A_x, $rA_6, $_Z );
        my $rpint_pos_x = polint( $A_x, $rA_4, $_rpint_pos );
        my $rpint_neg_x = polint( $A_x, $rA_4, $_rpint_neg );
        my $t_pose_rs_x = polint( $A_x, $rA_4, $_z_pose_rs );
        my $t_nege_rs_x = polint( $A_x, $rA_4, $_z_nege_rs );
        my $Qint_pos_x  = polint( $A_x, $rA_4, $_Qint_pos );
        my $rovp_min_x  = polint( $A_x, $rA_4, $_rovp_min );
        my $t_pmin_rs_x = polint( $A_x, $rA_4, $_z_pmin_rs );
        my $ke_pos_x    = polint( $A_x, $rA_4, $_ke_pos );
        my $work_pos_x  = polint( $A_x, $rA_4, $_work_pos );
        my $Disp_pos_x  = polint( $A_x, $rA_4, $_Disp_pos );

        my $rs_x        = exp($X_x);
        my $z_pose_rs_x = $rs_x - $t_pose_rs_x;
        my $z_nege_rs_x = $rs_x - $t_nege_rs_x;
        my $z_pmin_rs_x = $rs_x - $t_pmin_rs_x;

        my @vars;
        @vars[
          J_X,         J_Y,         J_Z,        J_rpint_pos, J_rpint_neg,
          J_z_pose_rs, J_z_nege_rs, J_Qint_pos, J_rovp_min,  J_z_pmin_rs,
          J_ke_pos,    J_work_pos,  J_Disp_pos
          ]
          = (
            $X_x,         $Y,           $Z_x,         $rpint_pos_x,
            $rpint_neg_x, $z_pose_rs_x, $z_nege_rs_x, $Qint_pos_x,
            $rovp_min_x,  $z_pmin_rs_x, $ke_pos_x,    $work_pos_x,
            $Disp_pos_x,
	    );

        push @{$rtab_x}, [@vars];
    }
    return $rtab_x;
}

sub _make_interpolated_tail_shock_table {

    # FIXME: TBD
    my ( $symmetry, $gamma, $rigam_6, $rigam_4, $rTables_loaded, $rgtab ) = @_;
    #  [ T, z, S1, S2 ]
    # where:
    #  T = ln(t) where t=scaled time for this point
    #  z = r-ct for this point
    #  (S1, S2) = values of Sigma on either side of the jump
    my $rtail_shock_table;
    return $rtail_shock_table;
}

sub _make_interpolated_blast_info {
    my ( $symmetry, $gamma, $rigam_4, $loc_gamma_table, $rTables_loaded, $rgtab ) = @_;

    my $rblast_info;
    if( defined($loc_gamma_table) && defined($rigam_4) ) {

        # If this table is an interpolated table then we can interpolate
        # P0 from nearby fits.
        # FIXME: this works now but needs a lot of refinement for better
        # Use best transformation for each variable.
        # Sintegral needs (alpha+2) weighting.
        my ( $jj, $jl, $ju, $num ) = @{$loc_gamma_table};

        if ( $jl >= 0 && $ju < $num ) {

            my $alpha = alpha_interpolate( $symmetry, $gamma );
            my $rtab;
            foreach my $igam ( @{$rigam_4} ) {
                my ( $gamma_i, $table_name_i ) =
                  @{ $rgtab->[$igam] };
                  ##@{ $rgamma_table->[$symmetry]->[$igam] };

                #my $rblast_info_i = get_blast_info( $symmetry, $gamma_i );
                my $rblast_info_i = $rTables_loaded->{$table_name_i}->{blast_info};
                my $alpha_i = alpha_interpolate( $symmetry, $gamma_i );
                push @{$rtab}, [ $rblast_info_i, $gamma_i, $alpha_i, $igam ];
            }

            # loop to interpolate all values
            my @keys = keys %{ $rtab->[0]->[0] };
            foreach my $key (@keys) {
                my ( $rx, $ry );
                my $nogo;
                foreach my $item ( @{$rtab} ) {
                    my ( $rblast_info_i, $gamma_i, $alpha_i, $igam ) = @{$item};
                    my $val_i = $rblast_info_i->{$key};
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
                $rblast_info->{$key} = $val;
            }

            # Old Linear interpolation, for reference
            # execute to see comparison with polynomial interpolation
            if (0) {
                my ( $jj, $jl, $ju, $num ) = @{$loc_gamma_table};
                if ( $jl >= 0 && $ju < $num ) {
                    ##my $rtab          = $rgamma_table->[$symmetry];
                    my $rtab          = $rgtab;

                    my ( $gamma_l, $table_name_l ) = @{$rtab->[$jl]};
                    my ( $gamma_u, $table_name_u ) = @{$rtab->[$ju]};
                    my $rblast_info_l = $rTables_loaded->{$table_name_l}->{blast_info};
                    my $rblast_info_u = $rTables_loaded->{$table_name_u}->{blast_info};

                    if ( defined($rblast_info_l) && defined($rblast_info_u) ) {
                        foreach my $key ( keys %{$rblast_info_l} ) {
                            my $val_l = $rblast_info_l->{$key};
                            my $val_u = $rblast_info_u->{$key};

                            if ( defined($val_l) && defined($val_u) ) {
                                my $val = _linear_interpolation(
                                    $gamma,   $gamma_l, $val_l,
                                    $gamma_u, $val_u
                                );
                                my $vsave = $rblast_info->{$key};
                                $rblast_info->{$key} = $val;

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
    return $rblast_info;
}

1;
