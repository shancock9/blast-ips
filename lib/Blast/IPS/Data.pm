package Blast::IPS::Data;

# This package handles access to individual data modules for
# a given symmetry and gamma value.  It returns a reference to
# the data in the appropriate module in Data/.  The actual
# data returned will differ from the data in the modules due to
# several post-processing steps which are done.

# usage:
#   my $rtables = get('symmetry'=>$symmetry, 'gamma'=>$gamma);
# where 
#   $rtables->{shock_table}   has the shock table 
#   $rtables->{energy_table}  has the energy table 
#   $rtables->{impulse_table} has the impulse table 
#   etc (See Blast::Data::README.txt)

use strict;
use warnings;
use 5.006;
our $VERSION = 0.1.1;
use Carp;
use Blast::IPS::Data::Index;
use Storable qw(dclone);

# This hash reference holds references to all loaded data tables. Each set of
# data has a unique hash key:  
our $rBlastData;

# Local hash lookup table
my $rmodule_names;

# array with BlastInfo names
my @vnames;

my $rindex_table; 

BEGIN {

    $rindex_table = Blast::IPS::Data::Index->get_index();
    sub get_index {
	return $rindex_table;
    }

    # Make a lookup table for module names
    #my $rindex_table = $Blast::IPS::Data::Index::rindex_table;
    #my $rindex_table = Blast::IPS::Data::Index->get_index();

    foreach my $symmetry ( 0, 1, 2 ) {
        next unless defined( $rindex_table->[$symmetry] );
        foreach my $item ( @{ $rindex_table->[$symmetry] } ) {
            my ( $gamma, $key, $module_name ) = @{$item};
            $rmodule_names->{$key} = $module_name;
        }
    }

    # These are the variable names currently in the 'blast_info' array
    @vnames = qw(
      symmetry
      gamma
      t_pcenter_zero_1
      t_pcenter_zero_2
      t_pcenter_min
      pcenter_min
      t_pcenter_max_2
      pcenter_max_2
      impulse_center_pos
      impulse_center_neg
      t_u_plus_c_zero
      dpdt_center_max
      t_dpdt_center_max
      pcenter_dpdt_center_max
      pcenter_initial_ratio
      KE_initial_ratio
      Sintegral_pos
      Sintegral_neg
      r_tail_shock
      z_tail_shock
      r_tail_formed
      z_tail_formed
    );
}

sub make_table_key {

    # Make a standard table key by concatenating the symmetry letter and the gamma;
    # For example, make the key 'S1.37' for spherical symmetry with gamma=1.37
    # This key will be used as the hash key for a module in the global table.
    my ( $sym, $gamma ) = @_;
    $gamma=~s/0+$//;
    my $table_key;
    my $gamma_pr = sprintf "%.4g", $gamma;
    if ( $sym =~ /\d/ ) {
        $sym = $sym == 0 ? 'P' : $sym == 1 ? 'C' : 'S';
    }
    $table_key = $sym . $gamma_pr;
    return $table_key;
}

sub _decrypt_table_key {
    my ($key) = @_;
    my ( $symmetry, $gamma );
    if ( $key =~ /^([SCP])(.*)$/ ) {
        $symmetry = $1;
        $gamma    = $2;
        $symmetry = ( $symmetry eq 'S' ) ? 2 : ( $symmetry eq 'C' ) ? 1 : 0;
    }
    return ( $symmetry, $gamma );
}

sub _merge_shock_and_energy_tables {

    # Merge the shock and energy tables. They have the same X values in 
    # column 1 but they have been stored separately for flexibility and readability.
    # The disadvantage is that we have to merge them and at the same time
    # verify that the X values are still identical.
    my ($key)=@_;
    my $table_error = sub {
        my ($msg) = @_;
        croak "Error merging table $key: $msg\n";
    };
    my $rshock_table  = $rBlastData->{$key}->{shock_table};
    my $renergy_table = $rBlastData->{$key}->{energy_table};
    my $num0          = @{$rshock_table};
    my $num1          = @{$renergy_table};
    my $rnew_table;

    if ( $num0 != $num1 ) {
        $table_error->("number of rows differ: $num0 != $num1");
        return;
    }
    for ( my $i = 0 ; $i < $num0 ; $i++ ) {
        my $X0 = $rshock_table->[$i]->[0];
        my $X1 = $renergy_table->[$i]->[0];
        if ( $X0 != $X1 ) {
            $table_error->("X values differ for row $i: $X0 != $X1");
            return;
        }

        # Combine variables in shock table and energy table,
        # leaving spots for T and dTdX, Z-X and dZ-X/dX
        my @vars = @{ $renergy_table->[$i] };
        shift @vars;    # remove leading X
        push @{$rnew_table}, [ @{ $rshock_table->[$i] }, 0, 0, 0, 0, @vars ];
    }

    # Success, replace shock table with combined shock and energy table
    $rBlastData->{$key}->{shock_table} = $rnew_table;
    return;
}

sub _convert_blast_info_to_hash {
    my ($key)=@_;
    my $rarray = $rBlastData->{$key}->{blast_info};
    return unless defined($rarray);
    my @vals    = @{$rarray};
    my $nvals   = @{$rarray};
    my $nvnames = @vnames;
    if ( $nvals + 2 != $nvnames ) { croak "nvals=$nvals + 2 != $nvnames\n" }

    my $rhash;
    for ( my $i = 0 ; $i < $nvals ; $i++ ) {
        my $name = $vnames[ $i + 2 ];
        $rhash->{$name} = $rarray->[$i];
    }

    # replace the array with a hash
    $rBlastData->{$key}->{blast_info} = $rhash;
}

sub _merge_shock_and_blast_info {
    my ($key) = @_;
    my $rarray = $rBlastData->{$key}->{shock_table_info};
    return unless defined($rarray);
    my (
        $symmetry,     $gamma,    $Max_Error, $N,
        $Energy_Error, $R_FD_MOC, $FD_Error,  $MOC_Error,
        $Interp_Error, $rs2,      $zs2
    ) = @{$rarray};

    # Merge selected values to the blast_table_info:
    $rBlastData->{$key}->{blast_info}->{Max_Error} = $Max_Error;
    $rBlastData->{$key}->{blast_info}->{N}         = $N;
    return;
}


sub get {

    # return reference to the data tables for a specific case, if they exists

    # The call parameter(s) are either the hash key for a module (like 'S1.4')
    # or a hash of call parameter values (either 'name' or 'symmetry' and
    # 'gamma')

    # The data is loaded and cached in $rDataTables->{$module_key},
    # where the key $module_key is in the standard name format (like 'S1.4')

    # Example of calls to load a table for spherical symmetry and gamma=1.4:
    # $rtables = get('gamma'=>1.4, 'symmetry'=>'S')
    # $rtables = get('gamma'=>1.4, 'symmetry'=>2)
    # $rtables = get('S1.4')
    # $rtables = get('table_name'=>'S1.4')

    # We return a reference to a copy of the tables to avoid changing values in
    # the cashed table
    my $rhash = _get_raw(@_);
    return unless $rhash;
    return dclone($rhash);
}

sub _get_raw {

    # Load a data module if necessary and return the hash key for it in 
    # $rDataTables, or return undef if it does not exist

    # a reference consulted:
    # https://stackoverflow.com/questions/1917261/how-can-i-dynamically-include-perl-modules-without-using-eval

    my @args = @_;
    my ( $key, $symmetry, $gamma );
    if ( @args == 1 ) {
        $key = $args[0];
	($symmetry, $gamma) = _decrypt_table_key($key);
    }
    else {
	my $nargs=@args;
        my (%hash) = @args;
        $key      = $hash{'table_name'};
        $symmetry = $hash{'symmetry'};
        $gamma    = $hash{'gamma'};
        if ( !$key ) { $key = make_table_key( $symmetry, $gamma ) }
        else         { ( $symmetry, $gamma ) = _decrypt_table_key($key) }
    }

    my $module_name   = $rmodule_names->{$key};
    return unless ($module_name);
    if ( !defined( $rBlastData->{$key} ) ) {
        eval {
            my $module = "Blast::IPS::Data::$module_name";
            ( my $file = $module ) =~ s|::|/|g;
            require $file . '.pm';
            $module->import();
            1;
        } or do {
            my $error = $@;
            print STDERR $error, "\n";
            return;
        };

	# now do any post-processing steps 
	_merge_shock_and_energy_tables($key);
        _convert_blast_info_to_hash($key); 
        _merge_shock_and_blast_info($key); 

        $rBlastData->{$key}->{symmetry}   = $symmetry;
        $rBlastData->{$key}->{gamma}      = $gamma;
        $rBlastData->{$key}->{table_name} = $key;
    }
    return $rBlastData->{$key}; 
}
1;
