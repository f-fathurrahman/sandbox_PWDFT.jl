function potxcir!(rhoir, epsxcir, vxcir)
    #
    # FIXME: Pass xc_calc as an input
    xc_calc = LibxcXCCalculator(x_id=1, c_id=12)
    # Only LDA is supported for the moment
    calc_epsxc_Vxc_LDA!(xc_calc, rhoir, epsxcir, vxcir)
    return
end

function potxcir!(rhoir, magir, epsxcir, vxcir, bxcir)
    xc_calc = LibxcXCCalculator(x_id=1, c_id=12)
    N = size(rhoir, 1)
    ndmag = size(magir, 2)
    @assert ndmag == 1

    rhoupdn = zeros(Float64, N, 2)
    vxcupdn = zeros(Float64, N, 2)
    @views rhoupdn[:,1] = 0.5*(rhoir[:] + magir[:,1])
    @views rhoupdn[:,2] = 0.5*(rhoir[:] - magir[:,1])

    calc_epsxc_Vxc_LDA!(xc_calc, rhoupdn, epsxcir, vxcupdn)

    # Set output for vxcir and bxcir
    @views vxcir[:] = 0.5*( vxcupdn[:,1] + vxcupdn[:,2] )
    @views bxcir[:,1] = 0.5*( vxcupdn[:,1] - vxcupdn[:,2] )

    return
end