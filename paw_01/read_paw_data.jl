# a helper function to read xml strings into array
function _read_xml_str_vector!(root_elem, tagname, Nr, f)
    # Get the element by its tagname
    pp_f = LightXML.get_elements_by_tagname(root_elem, tagname)[1]
    # Get the content (as string)
    pp_f_str = LightXML.content(pp_f)
    # Replace new line with space (to be parsed)
    pp_f_str = replace(pp_f_str, "\n" => " ")
    # Parse into an array
    spl_str = split(pp_f_str, keepempty=false)
    for i in 1:Nr
        f[i] = parse(eltype(f), spl_str[i])
    end
    return
end



function _read_paw_data(xroot, Nr::Int64, lmax::Int64, Nproj::Int64)
    #
    # Read some additional data that are present in case of PAW
    #
    pp_full_wfc = LightXML.get_elements_by_tagname(xroot, "PP_FULL_WFC")
    Nwfc = parse(Int64, LightXML.attributes_dict(pp_full_wfc[1])["number_of_wfc"])
    println("Nwfc = ", Nwfc)
    # XXX Compare Nwfc with number of projectors ?
    println("Nproj = ", Nproj)

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
    oc = paw_pp_occ

    # This will be paw.ae_rho_atc
    pp_ae_nlcc = zeros(Float64, Nr)
    _read_xml_str_vector!(pp_paw[1], "PP_AE_NLCC", Nr, pp_ae_nlcc)
    println("Done reading pp_ae_nlcc")
    ae_rho_atc = pp_ae_nlcc

    # paw.ae_vloc
    pp_ae_vloc = zeros(Float64, Nr)
    _read_xml_str_vector!(pp_paw[1], "PP_AE_VLOC", Nr, pp_ae_vloc)
    println("Done reading pp_ae_vloc")
    ae_vloc = pp_ae_vloc

    # Other parameters from PP_AUGMENTATION
    pp_nonlocal = LightXML.get_elements_by_tagname(xroot, "PP_NONLOCAL")
    pp_aug = LightXML.get_elements_by_tagname(pp_nonlocal[1], "PP_AUGMENTATION")

    # These attributes should be present in case of PAW
    #
    augshape = LightXML.attributes_dict(pp_aug[1])["shape"]
    println("augshape = ", augshape)

    raug = parse(Float64, LightXML.attributes_dict(pp_aug[1])["cutoff_r"])
    iraug = parse(Int64, LightXML.attributes_dict(pp_aug[1])["cutoff_r_index"])
    println("raug = ", raug, " iraug = ", iraug)

    qqq_eps = parse(Float64, LightXML.attributes_dict(pp_aug[1])["augmentation_epsilon"])
    println("qqq_eps = ", qqq_eps)

    lmax_aug = parse(Int64, LightXML.attributes_dict(pp_aug[1])["l_max_aug"])
    println("lmax_aug = ", lmax_aug)

    # CALL xmlr_readtag('shape', upf%paw%augshape )
    # CALL xmlr_readtag('cutoff_r', upf%paw%raug )
    # CALL xmlr_readtag('cutoff_r_index', upf%paw%iraug )
    # CALL xmlr_readtag('augmentation_epsilon', upf%qqq_eps )
    # CALL xmlr_readtag('l_max_aug', upf%paw%lmax_aug )


    # Prepare paw.pfunc
    pfunc = zeros(Float64, Nr, Nproj, Nproj)
    # NOTE: Nproj should be equal to Nwfc
    # FIXME: if psp has_so is true also prepare pfunc_rel
    for nb in 1:Nproj, mb in 1:nb
        pfunc[1:Nr, nb, mb] .= aewfc[1:Nr,nb] .* aewfc[1:Nr,mb]
        pfunc[(iraug+1):Nr,nb,mb] .= 0.0 # force to zero
        pfunc[1:Nr,mb,nb] .= pfunc[1:Nr,nb,mb] # symmetric w.r.t nb and mb
    end

    # Prepare paw.ptfunc
    # Pseudo wavefunctions (not only the ones for oc > 0)
    # All-electron wavefunctions
    ptfunc = zeros(Float64, Nr, Nproj, Nproj)
    for nb in 1:Nproj, mb in 1:nb
        ptfunc[1:Nr,nb,mb] .= pswfc[1:Nr,nb] .* pswfc[1:Nr,mb]
        ptfunc[(iraug+1):Nr,nb,mb] .= 0.0
        ptfunc[1:Nr,mb,nb] .= ptfunc[1:Nr,nb,mb]
    end

    paw_data = PAWData_UPF(
        ae_rho_atc, pfunc, ptfunc,
        ae_vloc, oc, raug, iraug, lmax_aug, core_energy, augshape
    )

    println("Pass here ...")

    return
end

