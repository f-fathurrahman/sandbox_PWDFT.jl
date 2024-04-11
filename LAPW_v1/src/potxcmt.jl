function potxcmt!(
    ia,
    atoms, atsp_vars,
    mt_vars, 
    rhomt, vxcmt
)

    isp = atoms.atm2species[ia]
    n = mt_vars.npmt[isp]

    # allocate local arrays
    rho = zeros(Float64, n)
    ex = zeros(Float64, n)
    ec = zeros(Float64, n)
    vxc = zeros(Float16, n)

    # Only for LDA, no spin
    
    nr = nrmt(is)
    nri = nrmti(is)

    backward_SHT!(mt_vars, isp, rhomt, rho)

    # LDA, spin-unpolarized case only
    CALL xcifc(xctype_, n=n, tempa=swidth, rho=rho, ex=ex, ec=ec, vx=vx, vc=vc)

    ! convert exchange and correlation energy densities to spherical harmonics
    CALL my_rfsht(nr, nri, ex, exmt_(:,ias))
    CALL my_rfsht(nr, nri, ec, ecmt_(:,ias))
    ! convert exchange-correlation potential to spherical harmonics
    CALL my_rfsht(nr, nri, vxc, vxcmt_(:,ias))
    write(*,*) 'my_potxcmt: shape(exmt)  = ', shape(exmt_)
    write(*,*) 'my_potxcmt: shape(ecmt)  = ', shape(ecmt_)
    write(*,*) 'my_potxcmt: shape(vxcmt) = ', shape(vxcmt_)
  ELSE 
    exmt_(1:n,ias) = ex(1:n)
    ecmt_(1:n,ias) = ec(1:n)
    vxcmt_(1:n,ias) = vxc(1:n)
  ENDIF 