function _read_paw_data(xroot, Nr::Int64, lmax::Int64)
    #
    # Read some additional data that are present in case of PAW
    #
    pp_full_wfc = LightXML.get_elements_by_tagname(xroot, "PP_FULL_WFC")
    Nwfc = parse(Int64, LightXML.attributes_dict(pp_full_wfc[1])["number_of_wfc"])
    println("Nwfc = ", Nwfc)
    # XXX Compare Nwfc with number of projectors ?

    aewfc = zeros(Float64, Nr, Nwfc)
    pswfc = zeros(Float64, Nr, Nwfc) # different from chi
    for iwf in 1:Nwfc
        tagname = "PP_AEWFC."*string(iwf)
        @views _read_xml_str_vector!(pp_full_wfc[1], tagname, Nr, aewfc[:,iwf])
        #
        tagname = "PP_PSWFC."*string(iwf)
        @views _read_xml_str_vector!(pp_full_wfc[1], tagname, Nr, pswfc[:,iwf])
    end

    #
    # These are data under PP_PAW
    #
    pp_paw = LightXML.get_elements_by_tagname(xroot, "PP_PAW")
    core_energy = parse(Float64, LightXML.attributes_dict(pp_paw[1])["core_energy"])
    println("core_energy = ", core_energy)

    pp_occ = LightXML.get_elements_by_tagname(pp_paw[1], "PP_OCCUPATIONS")
    Nocc = parse(Int64, LightXML.attributes_dict(pp_occ[1])["size"])
    println("Nocc = ", Nocc)
    # Nocc should be the same as Nwf?
    # This is paw.oc
    paw_pp_occ = zeros(Float64, Nocc)
    _read_xml_str_vector!(pp_paw[1], "PP_OCCUPATIONS", Nocc, paw_pp_occ)
    println("paw_pp_occ = ", paw_pp_occ)

    # This will be paw.ae_rho_atc
    pp_ae_nlcc = zeros(Float64, Nr)
    _read_xml_str_vector!(pp_paw[1], "PP_AE_NLCC", Nr, pp_ae_nlcc)
    println("Done reading pp_ae_nlcc")

    # paw.ae_vloc
    pp_ae_vloc = zeros(Float64, Nr)
    _read_xml_str_vector!(pp_paw[1], "PP_AE_VLOC", Nr, pp_ae_vloc)
    println("Done reading pp_ae_vloc")

    # Other parameters from PP_AUGMENTATION
    # CALL xmlr_readtag('shape', upf%paw%augshape )
    # CALL xmlr_readtag('cutoff_r', upf%paw%raug )
    # CALL xmlr_readtag('cutoff_r_index', upf%paw%iraug )
    # CALL xmlr_readtag('augmentation_epsilon', upf%qqq_eps )
    # CALL xmlr_readtag('l_max_aug', upf%paw%lmax_aug )

    # paw.core_energy

    # Prepare paw.pfunc

    # Prepare paw.ptfunc

    println("Pass here ...")

    return
end

