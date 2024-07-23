# energy step used for numerical calculation of energy derivatives
function get_deapwlo()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_deapwlo, LIBLAPW), Float64 ))
end

# maximum allowable APW order
function get_maxapword()
    return 4 # parameter, hardcoded
end

# maximum of apword over all angular momenta and species
function get_apwordmax()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_apwordmax, LIBLAPW), Int32 )) |> Int64
end

# polynomial order used for APW radial derivatives
function get_npapw()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_npapw, LIBLAPW), Int32 )) |> Int64
end

# maximum nlorb over all species
function get_nlomax()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_nlomax, LIBLAPW), Int32 )) |> Int64
end

# total number of local-orbitals
function get_nlotot()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_nlotot, LIBLAPW), Int32 )) |> Int64
end

# maximum number of local-orbitals
function get_maxlorb()
    return 200 # parameter, hardcoded
end

# maximum allowable local-orbital order
function get_maxlorbord()
    return 5 # parameter, hardcoded
end


# band energy search tolerance
function get_epsband()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_epsband, LIBLAPW), Float64 ))
end

# maximum allowed change in energy during band energy search; enforced only if
# default energy is less than zero
function get_demaxbnd()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_demaxbnd, LIBLAPW), Float64 ))
end

# minimum default linearisation energy over all APWs and local-orbitals
function get_e0min()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_e0min, LIBLAPW), Float64 ))
end

# difference between linearisation and Fermi energies when autolinengy is .true.
function get_dlefe()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_dlefe, LIBLAPW), Float64 ))
end

# maximum lorbord over all species
function get_lorbordmax()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_lorbordmax, LIBLAPW), Int32 )) |> Int64
end

# polynomial order used for local-orbital radial derivatives
function get_nplorb()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_nplorb, LIBLAPW), Int32 )) |> Int64
end

# maximum lorbl over all species
function get_lolmax()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_lolmax, LIBLAPW), Int32 )) |> Int64
end

# (lolmax+1)^2
function get_lolmmax()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_lolmmax, LIBLAPW), Int32 )) |> Int64
end

# if autolinengy is .true. THEN  the fixed linearisation energies are set to the
# Fermi energy minus dlefe
function get_autolinengy()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_autolinengy, LIBLAPW), Bool ))
end

# lorbcnd is .true. if conduction state local-orbitals should be added
function get_lorbcnd()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_lorbcnd, LIBLAPW), Bool ))
end

# conduction state local-orbital order
function get_lorbordc()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_lorbordc, LIBLAPW), Int32 )) |> Int64
end


# excess order of the APW and local-orbital functions
function get_nxoapwlo()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_nxoapwlo, LIBLAPW), Int32 )) |> Int64
end

# excess local orbitals
function get_nxlo()
    return unsafe_load(cglobal( (:__m_apwlo_MOD_nxlo, LIBLAPW), Int32 )) |> Int64
end

# total number of APW coefficients (l, m and order) for each species
function get_lmoapw()
    symbol = :__m_apwlo_MOD_lmoapw
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Int64, (nspecies,))
    # only load nspecies data
end

# number of local-orbitals
function get_nlorb()
    symbol = :__m_apwlo_MOD_nlorb
    nspecies = get_nspecies()
    return _load_automatic_array(symbol, Int64, (nspecies,))
    # only load nspecies data
end

# APW order
function get_apword()
    symbol = :__m_apwlo_MOD_apword
    maxlapw = get_maxlapw()
    maxspecies = get_maxspecies()
    nspecies = get_nspecies()
    apword = _load_automatic_array(symbol, Int64, (maxlapw+1,maxspecies))
    apword = OffsetArray( apword, 0:maxlapw, 1:maxspecies )
    return OffsetArray(apword[0:maxlapw, 1:nspecies], 0:maxlapw, 1:nspecies) # views of OffsetArray ?
end


# APW initial linearization energies
function get_apwe0()
    symbol = :__m_apwlo_MOD_apwe0
    maxlapw = get_maxlapw()
    maxapword = get_maxapword()
    maxspecies = get_maxspecies()
    nspecies = get_nspecies()
    apwe0 = _load_automatic_array(symbol, Float64, (maxapword,maxlapw+1,maxspecies))
    apwe0 = OffsetArray( apwe0, 1:maxapword, 0:maxlapw, 1:maxspecies )
    return OffsetArray( apwe0[1:maxapword, 0:maxlapw, 1:nspecies], 1:maxapword, 0:maxlapw,1:nspecies )
    # return nspecies data only, need to convert back to OffsetArray
end

# APW derivative order
function get_apwdm()
    symbol = :__m_apwlo_MOD_apwdm
    maxlapw = get_maxlapw()
    maxapword = get_maxapword()
    maxspecies = get_maxspecies()
    nspecies = get_nspecies()
    apwdm = _load_automatic_array(symbol, Float64, (maxapword,maxlapw+1,maxspecies))
    apwdm = OffsetArray( apwdm, 1:maxapword, 0:maxlapw, 1:maxspecies )
    return OffsetArray( apwdm[1:maxapword, 0:maxlapw, 1:nspecies], 1:maxapword, 0:maxlapw,1:nspecies )
    # return nspecies data only, need to convert back to OffsetArray
end

# XXX only return lmaxapw instead of maxlapw

# apwve is .true. if the linearization energies are allowed to vary
function get_apwve()
    symbol = :__m_apwlo_MOD_apwve
    maxlapw = get_maxlapw()
    maxapword = get_maxapword()
    maxspecies = get_maxspecies()
    nspecies = get_nspecies()
    apwve = _load_automatic_array(symbol, Int32, (maxapword,maxlapw+1,maxspecies))
    apwve = convert( Array{Bool,3}, apwve )
    apwve = OffsetArray( apwve, 1:maxapword, 0:maxlapw, 1:maxspecies )
    return OffsetArray( apwve[1:maxapword, 0:maxlapw, 1:nspecies], 1:maxapword, 0:maxlapw,1:nspecies )
    # return nspecies data only, need to convert back to OffsetArray
end

# lorbve is .true. if the linearisation energies are allowed to vary
function get_lorbve()
    symbol = :__m_apwlo_MOD_lorbve
    maxlorbord = get_maxlorbord()
    maxlorb = get_maxlorb()
    maxspecies = get_maxspecies()
    lorbve = _load_automatic_array(symbol, Int32, (maxlorbord,maxlorb,maxspecies))
    return convert( Array{Bool,3}, lorbve ) # read as Int32, then convert to Bool
end

# local-orbital derivative order
function get_lorbdm()
    symbol = :__m_apwlo_MOD_lorbdm
    maxlorbord = get_maxlorbord()
    maxlorb = get_maxlorb()
    maxspecies = get_maxspecies()
    return _load_automatic_array(symbol, Int64, (maxlorbord,maxlorb,maxspecies))
end

# local-orbital initial energies
function get_lorbe0()
    symbol = :__m_apwlo_MOD_lorbe0
    maxlorbord = get_maxlorbord()
    maxlorb = get_maxlorb()
    maxspecies = get_maxspecies()
    return _load_automatic_array(symbol, Float64, (maxlorbord,maxlorb,maxspecies))
end

# local-orbital order
function get_lorbord()
    symbol = :__m_apwlo_MOD_lorbord
    maxlorb = get_maxlorb()
    maxspecies = get_maxspecies()
    return _load_automatic_array(symbol, Int64, (maxlorb,maxspecies))
end



# local-orbital angular momentum
function get_lorbl()
    symbol = :__m_apwlo_MOD_lorbl
    maxlorb = get_maxlorb()
    maxspecies = get_maxspecies()
    return _load_automatic_array(symbol, Int64, (maxlorb,maxspecies))
end



# APW linearisation energies
function get_apwe()
    symbol = :__m_apwlo_MOD_apwe
    lmaxapw = get_lmaxapw()
    apwordmax = get_apwordmax()
    natmtot = get_natmtot()
    apwe = _load_allocatable_array(symbol, Float64, (apwordmax,lmaxapw+1,natmtot))
    return OffsetArray( apwe, 1:apwordmax, 0:lmaxapw, 1:natmtot )
end

# APW radial functions
function get_apwfr()
    symbol = :__m_apwlo_MOD_apwfr
    nrmtmax = get_nrmtmax()
    apwordmax = get_apwordmax()
    lmaxapw = get_lmaxapw()
    natmtot = get_natmtot()
    apwfr = _load_allocatable_array(symbol, Float64, (nrmtmax,2,apwordmax,lmaxapw+1,natmtot))
    return OffsetArray(apwfr, 1:nrmtmax, 1:2, 1:apwordmax, 0:lmaxapw, 1:natmtot)
end

# derivate of radial functions at the muffin-tin surface
function get_apwdfr()
    symbol = :__m_apwlo_MOD_apwdfr
    apwordmax = get_apwordmax()
    lmaxapw = get_lmaxapw()
    natmtot = get_natmtot()
    apwdfr = _load_allocatable_array(symbol, Float64, (apwordmax,lmaxapw+1,natmtot))
    return OffsetArray(apwdfr, 1:apwordmax, 0:lmaxapw, 1:natmtot)
end

# local-orbital energies
#ALLOCATE( lorbe(lorbordmax,maxlorb,natmtot) )
function get_lorbe()
    symbol = :__m_apwlo_MOD_lorbe
    lorbordmax = get_lorbordmax()
    maxlorb = get_maxlorb()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (lorbordmax,maxlorb,natmtot))
end

#  ALLOCATE( lofr(nrmtmax,2,nlomax,natmtot) )
function get_lofr()
    symbol = :__m_apwlo_MOD_lofr
    nrmtmax = get_nrmtmax()
    nlomax = get_nlomax()
    natmtot = get_natmtot()
    return _load_allocatable_array(symbol, Float64, (nrmtmax,2,nlomax,natmtot))
end

