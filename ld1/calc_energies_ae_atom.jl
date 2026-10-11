function calc_energies_ae_atom(
    zed, grid, rho, V_h, V_xc, epsxc, excgga, Nwf,
    Nspin, Enl, oc
)

    Nrmesh = grid.Nrmesh

    f1 = zeros(Float64, Nrmesh)
    f2 = zeros(Float64, Nrmesh)
    f3 = zeros(Float64, Nrmesh)
    f4 = zeros(Float64, Nrmesh)
    f5 = zeros(Float64, Nrmesh)
    for ir in 1:Nrmesh
        rhotot = rho[ir,1]
        if Nspin == 2
            rhotot = rhotot + rho[ir,2]
        end
        #
        # The integral for the energy due to the interaction with nuclei
        f1[ir] = -zed/grid.r[ir]*rhotot
        #
        # The integral for the Hartree energy
        f2[ir] = V_h[ir]*rhotot
        #
        # The integral for the exchange and correlation energy
        f3[ir] = epsxc[ir]*rhotot + excgga[ir]
        #
        # The integral to be subtracted to the sum of the eigenvalues to
        # get the kinetic energy
        #
        f5[ir] = -V_xc[ir,1]*rho[ir,1] - f1[ir] - f2[ir] - f4[ir]
        if Nspin == 2
            f5[ir] -= V_xc[ir,2]*rho[ir,2]
        end
    end
    #
    # Now compute the integrals
    E_nucl = integ_0_inf_dr(f1, grid, Nrmesh, 1)
    E_h = 0.5*integ_0_inf_dr(f2, grid, Nrmesh, 2)
    E_xc = integ_0_inf_dr(f3, grid, Nrmesh, 2)
    #
    # The kinetic energy is the sum of the eigenvalues plus the f5 integral
    E_kin = integ_0_inf_dr(f5, grid, Nrmesh, 1)
    for iwf in 1:Nwf
        if oc[iwf] > 0.0
            E_kin += oc[iwf]*Enl[iwf]
        end
    end
    E_tot = E_kin + E_nucl + E_h + E_xc
    # Return as a NamedTuple
    return (;
        E_tot, E_kin, E_nucl, E_h, E_xc
    )

end

