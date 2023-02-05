# To mimic type paw_in_upf
mutable struct PAWData_UPF

    # REAL(DP),ALLOCATABLE :: ae_rho_atc(:) AE core charge (pseudo ccharge is already included in upf)
    ae_rho_atc::Vector{Float64}

    # REAL(DP),ALLOCATABLE :: pfunc(:,:,:),&! Psi_i(r)*Psi_j(r)
    pfunc::Array{Float64,3}    
    
    #pfunc_rel(:,:,:) ! Psi_i(r)*Psi_j(r) small component
    # NOT IMPLEMENTED

    ptfunc::Array{Float64,3}  # as above, but for pseudo

    #REAL(DP),ALLOCATABLE :: ae_vloc(:)    ! AE local potential (pseudo vlocis already included in upf)
    ae_vloc::Vector{Float64}
    
    #aewfc_rel(:,:) ! as above, but for pseudo
    # NOT IMPLEMENTED

    #REAL(DP),ALLOCATABLE :: oc(:) ! starting occupation used to init becsum
    # they differ from US ones because they
    # are indexed on BETA functions, non on WFC
    oc::Vector{Float64}

    # REAL(DP) :: raug          ! augfunction max radius
    raug::Float64

    # INTEGER  :: iraug         ! index on rgrid closer to, and >, raug
    iraug::Int64

    # INTEGER  :: lmax_aug      ! max angmom of augmentation functions, it is ==
    # ! to 2* max{l of pseudized wavefunctions}
    # ! note that nqlc of upf also includes the angmom of
    # ! empty virtual channel used to generate local potential
    lmax_aug::Int64

    # REAL(DP)         :: core_energy   ! constant to add in order to get all-electron energy
    core_energy::Float64

    # CHARACTER(len=12):: augshape      ! shape of augmentation charge
    augshape::String
end