struct PsPot_UPF
    pspfile::String
    atsymb::String
    zval::Float64
    #
    is_nlcc::Bool
    is_ultrasoft::Bool
    is_paw::Bool
    # Radial vars
    Nr::Int64
    r::Array{Float64,1}
    rab::Array{Float64,1}
    #
    V_local::Array{Float64,1}
    # Projectors
    Nproj::Int64
    proj_l::Array{Int64,1}
    rcut_l::Array{Float64,1}
    proj_func::Array{Float64,2}
    Dij::Array{Float64,2}
    # Augmentation stuffs
    nqf::Int64
    nqlc::Int64
    qqq::Array{Float64,2}
    q_with_l::Bool
    qfuncl::Array{Float64,3}
    #
    Nwfc::Int64
    chi::Array{Float64,2}
    rhoatom::Array{Float64,1}
end

function PsPot_UPF( upf_file::String )
    
    xdoc = LightXML.parse_file(upf_file)

    # get the root element
    xroot = LightXML.root(xdoc)  # an instance of XMLElement
    
    #
    # Read some information from header
    #
    pp_header = LightXML.get_elements_by_tagname(xroot, "PP_HEADER")
    atsymb = LightXML.attributes_dict(pp_header[1])["element"]
    zval = Int64( parse( Float64, LightXML.attributes_dict(pp_header[1])["z_valence"] ) )
    lmax = parse( Int64, LightXML.attributes_dict(pp_header[1])["l_max"] )
    Nr = parse(Int64,LightXML.attributes_dict(pp_header[1])["mesh_size"])
    # Data generated by atompaw seems to have wrong mesh_size information in the header
    # I solved this problem by modifying the UPF file manually.

    str1 = LightXML.attributes_dict(pp_header[1])["core_correction"]
    if str1 == "F"  # XXX F is not parsed as False
        is_nlcc = false
    elseif str1 == "T"
        is_nlcc = true
    else
        is_nlcc = parse(Bool, str1)
    end

    str1 = LightXML.attributes_dict(pp_header[1])["is_ultrasoft"]
    if str1 == "F"  # XXX F is not parsed as False
        is_ultrasoft = false
    elseif str1 == "T"
        is_ultrasoft = true
    else
        is_ultrasoft = parse(Bool, str1)
    end

    str1 = LightXML.attributes_dict(pp_header[1])["is_paw"]
    if str1 == "F"  # XXX F is not parsed as False
        is_paw = false
    elseif str1 == "T"
        is_paw = true
    else
        is_paw = parse(Bool, str1)
    end

    #
    # Read radial mesh information: r and rab
    #
    pp_mesh = LightXML.get_elements_by_tagname(xroot, "PP_MESH")
    pp_r = LightXML.get_elements_by_tagname(pp_mesh[1], "PP_R")
    pp_r_str = LightXML.content(pp_r[1])
    pp_r_str = replace(pp_r_str, "\n" => " ")
    spl_str = split(pp_r_str, keepempty=false)

    @assert(length(spl_str) == Nr)

    r = zeros(Float64, Nr)
    for i = 1:Nr
        r[i] = parse(Float64, spl_str[i])
    end

    pp_rab = LightXML.get_elements_by_tagname(pp_mesh[1], "PP_RAB")
    pp_rab_str = LightXML.content(pp_rab[1])
    pp_rab_str = replace(pp_rab_str, "\n" => " ")
    spl_str = split(pp_rab_str, keepempty=false)

    @assert(length(spl_str) == Nr)

    rab = zeros(Float64, Nr)
    for i in 1:Nr
        rab[i] = parse(Float64, spl_str[i])
    end

    #
    # Core correction
    #
    if is_nlcc
        rho_atc = zeros(Float64,Nr)
        pp_nlcc = LightXML.get_elements_by_tagname(xroot, "PP_NLCC")
        pp_nlcc_str = LightXML.content(pp_nlcc[1])
        pp_nlcc_str = replace(pp_nlcc_str, "\n" => " ")
        spl_str = split(pp_nlcc_str, keepempty=false)
        for i in 1:Nr
            rho_atc[i] = parse(Float64,spl_str[i])
        end
    else
        rho_atc = zeros(Float64,1)
    end

    #
    # Local potential
    #
    pp_local = LightXML.get_elements_by_tagname(xroot, "PP_LOCAL")
    pp_local_str = LightXML.content(pp_local[1])
    pp_local_str = replace(pp_local_str, "\n" => " ")
    spl_str = split(pp_local_str, keepempty=false)

    @assert(length(spl_str) == Nr)

    V_local = zeros(Float64, Nr)
    for i in 1:Nr
        V_local[i] = parse(Float64, spl_str[i])*0.5 # convert to Hartree
    end

    #
    # Nonlocal projector
    #
    Nproj = parse(Int64,LightXML.attributes_dict(pp_header[1])["number_of_proj"])
    pp_nonlocal = LightXML.get_elements_by_tagname(xroot, "PP_NONLOCAL")
    proj_func = zeros(Float64,Nr,Nproj)
    proj_l = zeros(Int64,Nproj)
    rcut_l = zeros(Float64,Nproj)
    kbeta = zeros(Int64,Nproj)
    for iprj in 1:Nproj
        pp_beta = LightXML.get_elements_by_tagname(pp_nonlocal[1], "PP_BETA."*string(iprj))
        #
        proj_l[iprj] = parse( Int64, LightXML.attributes_dict(pp_beta[1])["angular_momentum"] )
        kbeta[iprj] = parse( Int64, LightXML.attributes_dict(pp_beta[1])["cutoff_radius_index"] )
        rcut_l[iprj] = r[kbeta[iprj]]
        #
        pp_beta_str = LightXML.content(pp_beta[1])
        pp_beta_str = replace(pp_beta_str, "\n" => " ")
        spl_str = split(pp_beta_str, keepempty=false)
        for i in 1:Nr
            proj_func[i,iprj] = parse(Float64,spl_str[i])*0.5 # Convert to Hartree
        end
    end
    kkbeta = maximum(kbeta)

    #
    # Dij matrix elements
    #
    Dij = zeros(Nproj,Nproj)
    Dij_temp = zeros(Nproj*Nproj)
    pp_dij = LightXML.get_elements_by_tagname(pp_nonlocal[1], "PP_DIJ")
    pp_dij_str = LightXML.content(pp_dij[1])
    pp_dij_str = replace(pp_dij_str, "\n" => " ")
    spl_str = split(pp_dij_str, keepempty=false)
    for i in 1:Nproj*Nproj
        Dij_temp[i] = parse(Float64,spl_str[i])
    end
    Dij = reshape(Dij_temp,(Nproj,Nproj))*2  # convert to Hartree

    #
    # augmentation stuffs:
    #
    nqf, nqlc, qqq, q_with_l, qfuncl = 
    _read_us_aug(is_ultrasoft, pp_nonlocal, r, Nproj, proj_l, kkbeta)

    #
    # Pseudo wave function
    #
    pp_pswfc = LightXML.get_elements_by_tagname(xroot, "PP_PSWFC")
    Nwfc = parse(Int64, LightXML.attributes_dict(pp_header[1])["number_of_wfc"] )
    chi = zeros(Float64,Nr,Nwfc)
    for iwf in 1:Nwfc
        tagname = "PP_CHI."*string(iwf)
        pp_chi = LightXML.get_elements_by_tagname(pp_pswfc[1], tagname)
        #
        pp_chi_str = LightXML.content(pp_chi[1])
        pp_chi_str = replace(pp_chi_str, "\n" => " ")
        spl_str = split(pp_chi_str, keepempty=false)
        for i in 1:Nr
            chi[i,iwf] = parse(Float64, spl_str[i])
        end
    end

    # rho atom
    rhoatom = zeros(Float64,Nr)
    pp_rhoatom = LightXML.get_elements_by_tagname(xroot, "PP_RHOATOM")
    pp_rhoatom_str = LightXML.content(pp_rhoatom[1])
    pp_rhoatom_str = replace(pp_rhoatom_str, "\n" => " ")
    spl_str = split(pp_rhoatom_str, keepempty=false)
    for i in 1:Nr
        rhoatom[i] = parse(Float64, spl_str[i])
    end

    LightXML.free(xdoc)

    return PsPot_UPF(upf_file, atsymb, zval,
        is_nlcc, is_ultrasoft, is_paw,
        Nr, r, rab, V_local, Nproj, proj_l, rcut_l, proj_func, Dij,
        nqf, nqlc, qqq, q_with_l, qfuncl,
        Nwfc, chi,
        rhoatom
    )

end


function _read_us_aug(
    is_ultrasoft::Bool,
    pp_nonlocal,
    r::Array{Float64,1},
    Nproj::Int64,
    proj_l,
    kkbeta::Int64
)

    if !is_ultrasoft
        # Dummy data
        nqf = -1
        nqlc = -1
        qqq = zeros(1,1)
        q_with_l = false
        qfuncl = zeros(1,1,1)
        return nqf, nqlc, qqq, q_with_l, qfuncl
    end

    Nr = size(r,1)
    pp_aug = LightXML.get_elements_by_tagname(pp_nonlocal[1], "PP_AUGMENTATION")
    
    # number of angular momenta in Q
    nqlc = parse(Int64, LightXML.attributes_dict(pp_aug[1])["nqlc"])
    
    # number of Q coefficients
    nqf = parse(Int64, LightXML.attributes_dict(pp_aug[1])["nqf"])
    
    # XXX Also need to read q_with_l
    # For GBRV this is false
    str1 = LightXML.attributes_dict(pp_aug[1])["q_with_l"]
    if str1 == "F"
        q_with_l = false
    elseif str1 == "T"
        q_with_l = true
    else
        q_with_l = parse(Bool, str1)
    end

    pp_q = LightXML.get_elements_by_tagname(pp_aug[1], "PP_Q")
    pp_q_str = LightXML.content(pp_q[1])
    pp_q_str = replace(pp_q_str, "\n" => " ")
    spl_str = split(pp_q_str, keepempty=false)

    qqq = zeros(Nproj,Nproj)
    qqq_temp = zeros(Nproj*Nproj)
    for i in 1:Nproj*Nproj
        qqq_temp[i] = parse(Float64, spl_str[i])
    end
    qqq = reshape(qqq_temp, (Nproj,Nproj)) # XXX convert to Ha?

    if nqf > 0
        qfcoef_tmp = zeros(Float64, nqf*nqlc*Nproj*Nproj)
        pp_qfcoef = LightXML.get_elements_by_tagname(pp_aug[1], "PP_QFCOEF")
        pp_qfcoef_str = LightXML.content(pp_qfcoef[1])
        pp_qfcoef_str = replace(pp_qfcoef_str, "\n" => " ")
        spl_str = split(pp_qfcoef_str, keepempty=false)
        for i in 1:length(qfcoef_tmp)
            qfcoef_tmp[i] = parse(Float64, spl_str[i])
        end
        qfcoef = reshape(qfcoef_tmp, nqf, nqlc, Nproj, Nproj)
        #
        pp_rinner = LightXML.get_elements_by_tagname(pp_aug[1], "PP_RINNER")
        pp_rinner_str = LightXML.content(pp_rinner[1])
        pp_rinner_str = replace(pp_rinner_str, "\n" => " ")
        spl_str = split(pp_rinner_str, keepempty=false)
        rinner = zeros(Float64, nqlc)
        for i in 1:nqlc
            rinner[i] = parse(Float64, spl_str[i])
        end
    end


    Nq = Int64( Nproj*(Nproj+1)/2 )
    qfunc = zeros(Float64, Nr, Nq)
    qfuncl = zeros(Float64,Nr,Nq,nqlc) # last index is l

    if !q_with_l
        for iprj in 1:Nproj, jprj in iprj:Nproj
            tagname = "PP_QIJ."*string(iprj)*"."*string(jprj)
            pp_qij = LightXML.get_elements_by_tagname(pp_aug[1], tagname)
            #
            first_idx = parse( Int64, LightXML.attributes_dict(pp_qij[1])["first_index"] )
            second_idx = parse( Int64, LightXML.attributes_dict(pp_qij[1])["second_index"] )
            comp_idx = parse( Int64, LightXML.attributes_dict(pp_qij[1])["composite_index"] )
            #
            pp_qij_str = LightXML.content(pp_qij[1])
            pp_qij_str = replace(pp_qij_str, "\n" => " ")
            spl_str = split(pp_qij_str, keepempty=false)
            # FIXME" using comp_idx?
            for i in 1:Nr
                qfunc[i,comp_idx] = parse(Float64,spl_str[i])
            end
        end

        # Prepare qfuncl
        for nb in 1:Nproj, mb in nb:Nproj
            # ijv is the combined (nb,mb) index
            ijv = round(Int64, mb*(mb-1)/2) + nb
            l1 = proj_l[nb]
            l2 = proj_l[mb]
            # copy q(r) to the l-dependent grid 
            for l in range( abs(l1-l2), stop=(l1+l2), step=2)
                #@printf("nb, mb, l = %d %d %d\n", nb, mb, l)
                @views qfuncl[1:Nr,ijv,l+1] = qfunc[1:Nr,ijv] # index l starts from 1
            end
            #
            if nqf > 0
                for l in range(abs(l1-l2), stop=l1+l2, step=2)
                    if rinner[l+1] > 0.0
                        ilast = 0
                        for ir in 1:kkbeta
                            if r[ir] < rinner[l+1]
                                ilast = ir
                            end
                        end
                        @views _setqfnew!( nqf, qfcoef[:,l+1,nb,mb], ilast, r, l, 2, qfuncl[:,ijv,l+1] )
                    end # if
                end # for
            end # if
        end
    else
        # Read QIJL
        for nb in 1:Nproj, mb in nb:Nproj
            ijv = round(Int64, mb*(mb-1)/2) + nb
            l1 = proj_l[nb]
            l2 = proj_l[mb]
            for l in range( abs(l1-l2), stop=(l1+l2), step=2)
                tagname = "PP_QIJL."*string(nb)*"."*string(mb)*"."*string(l)
                pp_qijl = LightXML.get_elements_by_tagname(pp_aug[1], tagname)
                #
                first_idx = parse( Int64, LightXML.attributes_dict(pp_qijl[1])["first_index"] )
                second_idx = parse( Int64, LightXML.attributes_dict(pp_qijl[1])["second_index"] )
                comp_idx = parse( Int64, LightXML.attributes_dict(pp_qijl[1])["composite_index"] )
                am_idx = parse( Int64, LightXML.attributes_dict(pp_qijl[1])["angular_momentum"] )
                if am_idx != l
                    println("WARNING: am_idx != l")
                end
                #
                str1 = LightXML.content(pp_qijl[1])
                str1 = replace(str1, "\n" => " ")
                spl_str = split(str1, keepempty=false)
                # FIXME" using comp_idx?
                for i in 1:Nr
                    qfuncl[i,comp_idx,l+1] = parse(Float64,spl_str[i])
                end
            end
        end
    end

    return nqf, nqlc, qqq, q_with_l, qfuncl
end


#
# Adapted from QE:
# subroutine setqfnew in upflib/upf_to_internal
#
function _setqfnew!(nqf, qfcoef, Nr, r, l, n, rho)
    #
    # FIXME: n is always equal to 2 ?
    #
    # Computes the Q function from its polynomial expansion (r < rinner)
    # On input: nqf = number of polynomial coefficients
    #    qfcoef(1:nqf) = the coefficients defining Q
    #          Nr = number of mesh point
    #        r(1:Nr)= the radial mesh
    #             l = angular momentum
    #             n = additional exponent, result is multiplied by r^n
    # On output:
    #      rho(1:Nr)= r^n * Q(r)
    #
    # INTEGER,  INTENT(in):: nqf, l, mesh, n
    # REAL(dp), INTENT(in) :: r(mesh), qfcoef(nqf)
    # REAL(dp), INTENT(out) :: rho(mesh)
    for ir in 1:Nr
        rr = r[ir]^2
        rho[ir] = qfcoef[1]
        for i in 2:nqf
           rho[ir] = rho[ir] + qfcoef[i] * rr^(i-1)
        end
        rho[ir] = rho[ir]*r[ir]^(l + n)
    end
    return
end
