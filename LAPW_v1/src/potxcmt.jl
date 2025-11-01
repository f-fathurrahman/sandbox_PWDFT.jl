# In this function loop over atoms are done directly
function potxcmt!(
    atoms, mt_vars, 
    rhomt, epsxcmt, vxcmt
)
    #
    # FIXME: Pass xc_calc as an input
    xc_calc = LibxcXCCalculator(x_id=1, c_id=12)
    #
    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    #
    for isp in 1:Nspecies
        # Allocate species dependent variables
        N = mt_vars.npmt[isp]
        rho = zeros(Float64, N)
        epsxc = zeros(Float64, N) # we combine x and c contrib to energy
        vxc = zeros(Float64, N)
        #
        for ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end
            #
            # Convert from Ylm to "real" space
            backward_SHT!(mt_vars, isp, rhomt[ia], rho)
            #
            # Only LDA is supported for the moment
            calc_epsxc_Vxc_LDA!(xc_calc, rho, epsxc, vxc)
            #
            # convert from "real" to spherical harmonics (Ylm)
            forward_SHT!(mt_vars, isp, epsxc, epsxcmt[ia])
            forward_SHT!(mt_vars, isp, vxc, vxcmt[ia])
        end
    end
    return
end


# In this function loop over atoms are done directly
# Spin-polarized
function potxcmt!(
    atoms, mt_vars, 
    rhomt, magmt, epsxcmt, vxcmt, bxcmt
)
    ndmag = size(magmt, 2)
    @assert ndmag == 1

    #
    # FIXME: Pass xc_calc as an input
    xc_calc = LibxcXCCalculator(x_id=1, c_id=12)
    #
    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    #
    for isp in 1:Nspecies
        # Allocate species dependent variables
        N = mt_vars.npmt[isp]
        mag = zeros(Float64, N)
        rho = zeros(Float64, N)
        #
        rhoupdn = zeros(Float64, N, 2)
        epsxc = zeros(Float64, N) # we combine x and c contrib to energy
        vxcupdn = zeros(Float64, N, 2)
        #
        vxc = zeros(Float64, N)
        bxc = zeros(Float64, N)
        #
        for ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end
            #
            # Convert from Ylm to "real" space
            backward_SHT!(mt_vars, isp, rhomt[ia], rho)
            backward_SHT!(mt_vars, isp, magmt[ia], mag)
            #
            # prepare rhoupdn
            @views rhoupdn[:,1] = 0.5*( rho[:] + mag[:] )
            @views rhoupdn[:,2] = 0.5*( rho[:] - mag[:] )
            #
            # Only LDA is supported for the moment
            calc_epsxc_Vxc_LDA!(xc_calc, rhoupdn, epsxc, vxcupdn)
            #
            @views vxc[:] = 0.5*( vxcupdn[:,1] + vxcupdn[:,2] )
            @views bxc[:] = 0.5*( vxcupdn[:,1] - vxcupdn[:,2] )
            #
            # convert from "real" to spherical harmonics (Ylm)
            forward_SHT!(mt_vars, isp, epsxc, epsxcmt[ia])
            forward_SHT!(mt_vars, isp, vxc, vxcmt[ia])
            @views forward_SHT!(mt_vars, isp, bxc, bxcmt[ia][:,1]) # only for spinpol ndmag=1
        end
    end
    return
end