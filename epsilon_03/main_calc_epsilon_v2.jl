using Printf
using PWDFT
import Serialization
using LinearAlgebra: det

function calc_epsilon(atoms, electrons, kpoints, evals, Mk; metal_like=false)

    Nstates = electrons.Nstates
    Nspin = electrons.Nspin

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

    Nkpt = kpoints.Nkpt
    Focc = repeat(electrons.Focc,1,Nkpt)
    FULL_OCC = 2.0 # for non-spinpol case

    εi = zeros(Float64,3,Nw)
    εr = zeros(Float64,3,Nw)

    Γ = 0.2
    #shift = 0.0
    shift = 0.55 # in eV for Si
    CellVolume = det(atoms.LatVecs)

    for ik in 1:Nkpt
        #
        println("ik = ", ik)
        M = Mk[ik]
        for jst in 1:Nstates    
            if Focc[jst,ik] >= FULL_OCC
                continue
            end
                
            for ist in 1:Nstates

                if ist == jst
                    continue
                end

                if Focc[ist,ik] < 0.5e-4*FULL_OCC
                    continue # skip the ist
                end

                if abs(Focc[jst,ik]-Focc[ist,ik]) < 1.0e-3*FULL_OCC
                    continue
                end
                #
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
                    denum = ( ddw2 + Γ^2 * w^2 )* ΔE
                    ff = Ha2eV^3*Focc[ist,ik]
                    for i in 1:3
                        εi[i,iw] = εi[i,iw] + M[i,ist,jst]*ff*Γ*w/denum   
                        εr[i,iw] = εr[i,iw] + M[i,ist,jst]*ff*ddw/denum
                    end
                end
            end
        end # states

        if metal_like
            for ist in 1:Nstates
                for iw in 1:Nw
                    w = wgrid[iw]
                    #
                    #df = w0gauss( (evals[ist,ik] - E_fermi)/degauss, ngauss)
                    #denum = (( w^4 + intrasmear^2 * w^2 )*degauss )
                    #mmf = M[i,ist,ist] * 0.5 * FULL_OCC * Ha2eV^2
                    denum = 1
                    Γ = 1 # FIXME !
                    for i in 1:3
                        εi[i,iw] = εi[i,iw] + mmf * df * Γ * w / denum
                        εr[i,iw] = εr[i,iw] - mmf * df * w^2 / denum
                    end
                end # iw
            end
        end

    end # over k-points

    C = 8*pi/(CellVolume*Nkpt)

    εr[:,:] = 1.0 .+ εr[:,:]*C
    εi[:,:] =        εi[:,:]*C

    Serialization.serialize("wgrid.data", wgrid)
    Serialization.serialize("epsr.data", εr)
    Serialization.serialize("epsi.data", εi)

    return
end

function load_data()
    print("Read data ...")
    atoms = Serialization.deserialize("atoms.data")
    evals = Serialization.deserialize("evals_nscf.data")
    kpoints = Serialization.deserialize("kpoints_nscf.data")
    electrons = Serialization.deserialize("electrons_nscf.data")
    Mk = Serialization.deserialize("Mk.data")
    println(" done")
    return atoms, electrons, kpoints, evals, Mk
end

function main()
    atoms, electrons, kpoints, evals, Mk = load_data()
    @time calc_epsilon(atoms, electrons, kpoints, evals, Mk)
    #@time calc_epsilon(atoms, electrons, kpoints, evals)
end

main()