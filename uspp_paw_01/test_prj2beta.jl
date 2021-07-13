using Printf
using SpecialFunctions: sphericalbesselj
import LightXML

using PWDFT
const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("PsPot_UPF.jl")
include("integ_simpson.jl")

function _build_prj_interp_table( psp::PsPot_UPF, pw::PWGrid )

    ecutwfc = pw.ecutwfc
    CellVolume = pw.CellVolume

    cell_factor = 1.0 # XXX HARDCODED
    dq = 0.01 # XXX HARDCODED

    ndm = psp.kkbeta
    Nproj = psp.Nproj

    nqx = floor( Int64, (sqrt(2*ecutwfc)/dq + 4)*cell_factor )

    psp.prj_interp_table = zeros(Float64,nqx,Nproj)

    aux = zeros(Float64, ndm)
    pref = 4*pi/sqrt(CellVolume)

    for ibeta in 1:Nproj
        l = psp.proj_l[ibeta]
        for iq in 1:nqx
            qi = (iq - 1) * dq
            for ir in 1:psp.kkbeta
                jlqr = sphericalbesselj(l, qi*psp.r[ir])
                aux[ir] = psp.proj_func[ir,ibeta] * psp.r[ir] * jlqr
            end
            vqint = integ_simpson( psp.kkbeta, aux, psp.rab )
            psp.prj_interp_table[iq, ibeta] = vqint * pref
        end
    end

    return
end


function eval_proj_G(psp::PsPot_UPF, iprjl::Int64, Gm::Float64)
    #
    dq = 0.01 # HARDCODED
    tab = psp.prj_interp_table
    #
    # Interpolation procedure
    px = Gm/dq - floor(Int64, Gm/dq)
    ux = 1.0 - px
    vx = 2.0 - px
    wx = 3.0 - px
    i0 = floor(Int64, Gm/dq) + 1
    i1 = i0 + 1
    i2 = i0 + 2
    i3 = i0 + 3
    Vq = tab[i0,iprjl] * ux * vx * wx / 6.0 +
         tab[i1,iprjl] * px * vx * wx / 2.0 -
         tab[i2,iprjl] * px * ux * wx / 2.0 +
         tab[i3,iprjl] * px * ux * vx / 6.0
    return Vq
end

function test_prj2beta(
    atoms::Atoms, pspots
)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    atpos = atoms.positions

    # 4: indexed from 0:3
    # 0, 1, 2, 3  -> l indexed
    # 1, 2, 3, 4  -> l + 1

    # -3:3
    # -3, -2, -1, 0, 1, 2, 3  -> m
    #  1,  2,  3, 4, 5, 6, 7  -> 4 + m, lmax = 3 + 1

    # -2, -1, 0, 1, 2  -> m
    #  1,  2, 3, 4, 5  -> 3 + m, lmax = 2 + 1

    prj2beta = Array{Int64}(undef,3,Natoms,4,7)
    prj2beta[:,:,:,:] .= -1   # set to invalid index

    NbetaNL = 0

    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        for l in 0:psp.lmax, m in -l:l
            #println("------")
            for iprj in 1:psp.Nproj_l[l+1]
                NbetaNL = NbetaNL + 1
                prj2beta[iprj,ia,l+1,m+psp.lmax+1] = NbetaNL
                @printf("ibeta=%3d ia=%3d l=%3d m=%3d\n", NbetaNL, ia, l, m)
            end
        end
    end

    println("NbetaNL = ", NbetaNL)

    Gm = 1.0
    for ia in 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        iprjl = 0
        for l in 0:psp.lmax
            for iprj in 1:psp.Nproj_l[l+1]
                iprjl = iprjl + 1
                for m in -l:l
                    ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                    pg = eval_proj_G(psp, iprjl, Gm)
                    @printf("ibeta=%3d iprjl=%3d ia=%3d l,m=(%3d,%3d) Vg=%18.10f\n",
                        ibeta, iprjl, ia, l, m, Vg)
                end
            end
        end
    end

    return
end

function main_ZnO()
    atoms = Atoms( xyz_string_frac=
        """
        4

        Zn      0.3333333   0.6666667   0.0000000
        Zn      0.6666667   0.3333333   0.5000000
        O       0.3333333   0.6666667   0.3450000
        O       0.6666667   0.3333333   0.8450000
        """, in_bohr=true,
        LatVecs = gen_lattice_hexagonal( 3.2495*ANG2BOHR, 5.2069*ANG2BOHR ) )

    #pspfiles = [joinpath(DIR_PSP, "Zn-q2.gth"),
    #            joinpath(DIR_PSP, "O-q6.gth")]

    pw = PWGrid(15.0, atoms.LatVecs)

    pspfiles = ["/home/efefer/pseudo/HGH/Zn.pz-hgh.UPF",
                "/home/efefer/pseudo/HGH/O.pz-hgh.UPF"]


    Nspecies = atoms.Nspecies
    @assert length(pspfiles) == Nspecies
    pspots = Array{PsPot_UPF,1}(undef,Nspecies)
    #pspots = Array{PsPot_GTH,1}(undef,Nspecies)
    for isp in 1:Nspecies
        #pspots[isp] = PsPot_GTH(pspfiles[isp])
        pspots[isp] = PsPot_UPF(pspfiles[isp])
        _build_prj_interp_table(pspots[isp], pw)
        #println(pspots[isp])
    end

    test_prj2beta(atoms, pspots)

end

main_ZnO()