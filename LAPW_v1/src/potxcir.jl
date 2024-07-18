function potxcir!(rhoir, epsxcir, vxcir)
    #
    # FIXME: Pass xc_calc as an input
    xc_calc = LibxcXCCalculator(x_id=1, c_id=12)
    # Only LDA is supported for the moment
    calc_epsxc_Vxc_LDA!(xc_calc, rhoir, epsxcir, vxcir)
    return
end