# In this function loop over atoms are done directly
function potxcmt!(
    atoms, mt_vars, 
    rhomt, epsxcmt, vxcmt
)

    # FIXME: make this an input argument
    xc_calc = LibxcXCCalculator(x_id=1, c_id=12)

    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    for isp in 1:atoms.Nspecies
        #
        N = mt_vars.npmt[isp]
        rho = zeros(Float64, N)
        epsxc = zeros(Float64, N) # we combine x and c contrib to energy
        vxc = zeros(Float64, N)
        #
        for ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end

            # Convert from Ylm to "real" space
            backward_SHT!(mt_vars, isp, rhomt[ia], rho)

            calc_epsxc_Vxc_LDA!(xc_calc, rho, epsxc, vxc)

            # convert from "real" to spherical harmonics (Ylm)
            forward_SHT!(mt_vars, isp, epsxc, epsxcmt[ia])
            forward_SHT!(mt_vars, isp, vxc, vxcmt[ia])
        end
    end

    return
end