The blast tables for a specific case are stored in a hash with the layout given
below. They are referenced through the index module Index.pm.  You must run
'make' to update the index whenever a data module is added or removed.

The hash keys encode the symmetry and gamma; for example key 'S1.4' is for
Spherical Symmetry and gamma=1.4.

The basic keys are as follows:

$rBlastData->{'S1.4'} = {
    shock_table      => [...],
    energy_table     => [...],
    impulse_table    => [...],
    shock_table_info => [...],
    tail_shock_table => [...],
};

The tables were prepared from calculations using the finite difference method
and method of characteristics.  Note that the first two tables must have the 
same values in the first rows (X values) so that they can later be combined 
into a single table.

The formats of the individual tables is as follows:

***   shock_table      => [...],

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

***   energy_table     => [...],

    # Tables of residual energy due to primary shock for an explosion in an
    # ideal homogeneous atmosphere.  The tables were created by integrating the
    # computed shock histories and interpolating to the same X values as the
    # main shock tables.

    # The table format is:
    #  [ X, E1, dE1/dX, E, dEdX, I-X, d(I-X)dX ]
    # where
    #  X = ln(lambda) where lambda = scaled range
    #  E1 = residual energy from primary shock
    #  dE1dX = derivative of E1 with respect to X
    #  E = residual energy from both primary shock and tail shock
    #  dEdX = derivative of E with respect to X
    #  I-X = positive phase overpressure impulse
    #  d(I-X)/dX = derivative of I with respect to X


***   impulse_table    => [...],


    # [ X, Y, Z, rpint_pos, rpint_neg, z_pose_th, z_nege_th, Qint_pos, rovp_min_th, z_pmin_th,
    #       ke_pos, work_pos, Disp_pos]
    # where
    #  X = ln(lambda) where lambda = scaled range of the shock front
    #  Y = ln(ovp ratio) at shock front
    #  Z = ln(r-ct) at shock front
    #  rpint_pos = r**(n/2) x scaled positive phase overpressure impulse, where n=0,1,2 is the symmetry
    #  rpint_neg = r**(n/2) x scaled negative phase overpressure impulse, where n=0,1,2 is the symmetry
    #  z_pose_th = z value at end of positive phase of the time history, where z=r-ct
    #  z_nege_th = z value at end of negative phase of the time history, if any
    #  Qint_pos = ln(dynamic pressure impulse = integral of rho*u**2)
    #  rovp_min_th = r**(n/2) * minimum overpressure ratio on a time history at this location
    #  z_pmin_th = z value at minimum overpressure ratio in a time history
    #  ke_pos = kinetic energy in positive phase of the wave
    #  work_pos = work of positive phase of the wave
    #  Disp_pos = ln(peak displacement of a Lagrange point )

***   shock_table_info => [...],

      # symmetry = 0,1,2 for Plane, Cylindrical, Spherical
      # gamma = ideal gas gamma
      # Max-Error = the sum of max absolute FD + MOC + interpolation errors;
      # N=number of points from center to shock in Finite Difference run
      # FD-Error = maximum error in FD run by comparing with N/2 run
      # Interp-Error = maximum cubic interpolation error in overpressure
      # MOC-Error = maximum error in MOC run
      # Energy-Error = energy error in the FD run at 0.1 shock overpressure ratio
      # rs2 = radius of formation of second shock
      # zs2 = value of z for formation of second shock
      shock_table_info => [ 2, 1.4, 7.21e-07, 32000, 1.27e-08, 10, 5.4e-08, 1.67e-07, 5e-07, 2090.5, -0.6309, ],

***   tail_shock_table => [...],
    # Tables of values along the trajectory of the second shock, if any
    # Note that information about tail shock work is contained in EnergyTables.pm
    # The table format is:
    #  [ T, z, S1, S2 ]
    # where:
    #  T = ln(t) where t=scaled time for this point
    #  z = r-ct for this point
    #  (S1, S2) = values of Sigma on either side of the jump
    tail_shock_table => [
        [ 7.6454603, -0.6309,     -0.039668803, -0.039668803 ],
};

