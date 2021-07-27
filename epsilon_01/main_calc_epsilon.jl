using Printf
using PWDFT
import Serialization

function calc_dipole_matrix!( Ham, psiks, ik, M; metal_like=false )

    @assert Ham.electrons.Nspin == 1

    fill!( M, 0.0 + im*0.0 )
    
    psi = psiks[ik]
    Npw = Ham.pw.gvecw.Ngw[ik]
    Nstates = Ham.electrons.Nstates
    idx_gw2g = Ham.pw.gvecw.idx_gw2g[ik]
    Focc = Ham.electrons.Focc
    G = Ham.pw.gvec.G
    k = Ham.pw.gvecw.kpoints.k[:,ik]

    FULL_OCC = 2.0 # for non-spinpol case
    for jst in 1:Nstates
        if Focc[jst,ik] < FULL_OCC # empty states
            for ist in 1:Nstates
                if ist == jst
                    continue
                end
                if Focc[ist,ik] >= 0.5e-4*FULL_OCC # occupied states
                    for igw in 1:Npw
                        ig = idx_gw2g[igw]
                        caux = conj(psi[igw,ist])*psi[igw,jst]
                        for i in 1:3
                            M[i,ist,jst] = M[i,ist,jst] + ( G[i,ig] + k[i] ) * caux
                        end
                    end
                end
            end
        end # if
    end
  
    # The diagonal terms are taken into account only if the system is treated like a metal, not
    # in the intraband therm. Because of this we can recalculate the diagonal component of the dipole
    # tensor directly as we need it for the intraband term, without interference with interband one.
    if metal_like
        for ist in 1:Nstates
            for igw in 1:Npw
                ig = idx_gw2g[igw]
                caux = conj(psi[igw,ist])*psi[igw,ist]
                for i in 1:3
                    M[i,ist,ist] = M[i,ist,ist] + ( G[i,ig] + k[i] ) * caux
                end
            end
        end
    end

    return
end


function main()
    
    print("Read data ...")
    Ham = Serialization.deserialize("Ham_nscf.data");
    psiks = Serialization.deserialize("psiks.data");
    evals = Serialization.deserialize("evals.data");
    println(" done")

    Nstates = Ham.electrons.Nstates
    M_aux = zeros(ComplexF64,3,Nstates,Nstates)
    M = zeros(Float64,3,Nstates,Nstates)

    # Grid
    Nw = 500
    wmin = 0.0
    wmax = 30.0
    #
    wgrid = zeros(Float64,Nw)
    dw = (wmax - wmin)/(Nw-1)
    for iw in 1:Nw
        wgrid[iw] = wmin + (iw-1)*dw
    end

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Focc = Ham.electrons.Focc
    FULL_OCC = 2.0 # for non-spinpol case

    εi = zeros(Float64,3,Nw)
    εr = zeros(Float64,3,Nw)

    intersmear = 0.2
    shift = 0.0
    CellVolume = Ham.pw.CellVolume

    for ik in 1:Nkpt
        println("ik = ", ik)
        calc_dipole_matrix!( Ham, psiks, ik, M_aux; metal_like=false )
        for i in 1:length(M)
            M[i] = real( M_aux[i] * conj(M_aux[i]) )
        end
        #
        for jst in 1:Nstates    
            if Focc[jst,ik] < FULL_OCC
                for ist in 1:Nstates
                    if ist == jst
                        continue
                    end
                    if Focc[ist,ik] >= 0.5e-4*FULL_OCC
                        #
                        #@printf("%d %d\n", ist, jst)
                        #
                        if abs(Focc[jst,ik]-Focc[ist,ik]) < 1.0e-3*FULL_OCC
                            continue
                        end
                        # transition energy
                        #
                        ΔE = ( evals[jst,ik] - evals[ist,ik] ) * Ha2eV + shift
                        #
                        # loop over frequencies
                        #
                        for iw in 1:Nw
                            w = wgrid[iw]
                            #
                            ddw = (ΔE^2 - w^2)
                            ddw2 = ddw^2
                            denum = ( ddw2 + intersmear^2 * w^2 )* ΔE
                for i in 1:3
                    εi[i,iw] = εi[i,iw] + M[i,ist,jst]*Ha2eV^3*Focc[ist,ik]*intersmear*w/denum   
                    εr[i,iw] = εr[i,iw] + M[i,ist,jst]*Ha2eV^3*Focc[ist,ik]*ddw/denum
                end
                        end
                    end # if
                end
            end # if
        end
    end

    #C =  64.0*pi/(CellVolume*Nkpt)
    #C =  32.0*pi/(CellVolume*Nkpt)
    C = 8*pi/(CellVolume*Nkpt)
    εr[:,:] = 1.0 .+ εr[:,:]*C
    εi[:,:] =        εi[:,:]*C

    Serialization.serialize("wgrid.data", wgrid)
    Serialization.serialize("epsr.data", εr)
    Serialization.serialize("epsi.data", εi)

    return
end

main()