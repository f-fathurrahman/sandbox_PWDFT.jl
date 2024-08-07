mutable struct SymmetryVars
    # type of symmetry allowed for the crystal
    #  0 : only the identity element is used
    #  1 : full symmetry group is used
    #  2 : only symmorphic symmetries are allowed
    symtype::Int64
    #
    # number of Bravais lattice point group symmetries
    nsymlat::Int64
    #
    # Bravais lattice point group symmetries
    symlat::Vector{Matrix{Int64}} #(3,3,48)
    #
    # determinants of lattice symmetry matrices (1 or -1)
    symlatd::Vector{Int64} #(48)
    #
    # index to inverses of the lattice symmetries
    isymlat::Vector{Int64} #(48)
    #
    # lattice point group symmetries in Cartesian coordinates
    symlatc::Vector{Matrix{Float64}} #(3,3,48)
    #
    # tshift is .true. if atomic basis is allowed to be shifted
    tshift::Bool
    #
    # tsyminv is .true. if the crystal has inversion symmetry
    tsyminv::Bool 
    #
    # maximum of symmetries allowed
    maxsymcrys::Int64 # 192
    #
    # number of crystal symmetries
    nsymcrys::Int64
    #
    # crystal symmetry translation vector in lattice and Cartesian coordinates
    vtlsymc::Array{Float64,2} # (3,maxsymcrys)
    vtcsymc::Array{Float64,2} # (3,maxsymcrys)
    #
    # tv0symc is .true. if the translation vector is zero
    tv0symc::Vector{Bool} #(maxsymcrys)
    #
    # spatial rotation element in lattice point group for each crystal symmetry
    lsplsymc::Vector{Int64} #(maxsymcrys)
    #
    # global spin rotation element in lattice point group for each crystal symmetry
    lspnsymc::Vector{Int64} #(maxsymcrys)
    #
    # equivalent atom index for each crystal symmetry
    ieqatom::Array{Int64,2}
    #
    # eqatoms(ia,ja,is) is .true. if atoms ia and ja are equivalent
    eqatoms::Array{Bool,2}
    #
    # number of site symmetries
    nsymsite::Vector{Int64} #(:)
    #
    # site symmetry spatial rotation element in lattice point group
    lsplsyms::Array{Int64,2} #(:,:)
    # XXX: used in genlmirep which will be used in DOS calculation
    #
    # site symmetry global spin rotation element in lattice point group
    lspnsyms::Array{Int64,2} #(:,:)
end

# FIXME: some variables need to be declared as vector of matrix
function SymmetryVars()
    symtype = 1
    nsymlat = 1
    
    #symlat = zeros(Int64,3,3,48)
    symlat = Vector{Matrix{Int64}}(undef,48)
    for i in 1:48
        symlat[i] = zeros(Int64, 3, 3)
    end
    
    symlatd = zeros(Int64, 48)
    isymlat = zeros(Int64, 48)
    
    #symlatc = zeros(Float64,3,3,48)
    symlatc = Vector{Matrix{Float64}}(undef,48)
    for i in 1:48
        symlatc[i] = zeros(Float64, 3, 3)
    end

    tshift = false
    tsyminv = false
    maxsymcrys = 192
    nsymcrys = 1 # must be updated properly
    vtlsymc = zeros(Float64,3,maxsymcrys)
    vtcsymc = zeros(Float64,3,maxsymcrys)
    tv0symc = zeros(Bool,maxsymcrys)
    lsplsymc = zeros(Bool,maxsymcrys)
    lspnsymc = zeros(Int64,maxsymcrys)
    
    # must be properly initialized later
    ieqatom = zeros(Int64,1,1)
    eqatoms = zeros(Bool,1,1)

    # must be properly initialized later
    nsymsite = zeros(Int64,1)
    lsplsyms = zeros(Int64,1,1)
    lspnsyms = zeros(Int64,1,1)

    return SymmetryVars( 
        symtype, nsymlat, symlat, symlatd, isymlat, symlatc,
        tshift, tsyminv, maxsymcrys, nsymcrys, vtlsymc, vtcsymc,
        tv0symc, lsplsymc, lspnsymc, ieqatom, eqatoms,
        nsymsite, lsplsyms, lspnsyms
    )
end
