ENERGY_VARS_LIST = [
    :evalsum,
    :engykn,
    :engykncr,
    :engynn,
    :engyen,
    :engyhar,
    :engycl,
    :engyvcl,
    :engymad,
    :engyvxc,
    :engybxc,
    :engybext,
    :engyx,
    :engyc,
    :entrpy,
    :engyts,
    :engytot
]

for v in ENERGY_VARS_LIST
    str_funcs_energy = """
    function get_$(v)()
        return unsafe_load(cglobal((:__m_energy_MOD_$(v), LIBLAPW), Float64))
    end
    """
    eval(Meta.parse(str_funcs_energy))
end



#=

!--------------------------!
!     energy variables     !
!--------------------------!

! eigenvalue sum
real(8) evalsum

! electron kinetic energy
real(8) engykn

! core electron kinetic energy
real(8) engykncr

! nuclear-nuclear energy
real(8) engynn

! electron-nuclear energy
real(8) engyen

! Hartree energy
real(8) engyhar

! Coulomb energy (E_nn + E_en + E_H)
real(8) engycl

! electronic Coulomb potential energy
real(8) engyvcl

! Madelung term
real(8) engymad

! exchange-correlation potential energy
real(8) engyvxc

! exchange-correlation effective field energy
real(8) engybxc

! energy of external global magnetic field
real(8) engybext

! exchange energy
real(8) engyx

! correlation energy
real(8) engyc

! electronic entropy
real(8) entrpy

! entropic contribution to free energy
real(8) engyts

! total energy
real(8) engytot
=#