

Dielectric tensor:
$$
\epsilon_{\alpha\beta}=1+\frac{4\pi e^{2}}{\Omega N_{\mathbf{k}}m^{2}}\sum_{n,n'}\sum_{\mathbf{k}}
\frac{\hat{\mathbf{M}}^{nn'}_{\alpha\beta}}{\left(E_{\mathbf{k}n'}-E_{\mathbf{k}n}\right)^{2}}\left(\frac{f(E_{\mathbf{k}n})}{E_{\mathbf{k}n'}-E_{\mathbf{k}n}+\hbar\omega+\imath\hbar\Gamma}+\frac{f(E_{\mathbf{k}n})}{E_{\mathbf{k}n'}-E_{\mathbf{k}n}-\hbar\omega+\imath\hbar\Gamma}\right)
$$

Dipole matrix elements:
$$
\begin{align*}
\hat{\mathbf{M}}_{\alpha\beta} & =\left\langle u_{\mathbf{k}n'}\left|\hat{\mathbf{p}}_{\alpha}\right|u_{\mathbf{k}n}\right\rangle \left\langle u_{\mathbf{k}n}\left|\hat{\mathbf{p}}_{\beta}^{*}\right|u_{\mathbf{k}n'}\right\rangle \\
 & =\left(\sum_{\mathbf{G}}a_{\mathbf{k}n\mathbf{G}}^{*}a_{\mathbf{k}n'\mathbf{G}}(G_{\alpha}+k_{\alpha})\right)\left(\sum_{\mathbf{G}}a_{\mathbf{k}n\mathbf{G}}^{*}a_{\mathbf{k}n'\mathbf{G}}(G_{\beta}+k_{\beta})\right)
\end{align*}
$$
Calculate dipole matrix elements (when squared will result in $\hat{\mathbf{M}}$):

```julia
# ... snipped
fill!( M, 0.0 + im*0.0 )
FULL_OCC = 2.0 # for non-spinpol case
for jst in 1:Nstates
    if Focc[jst,ik] >= FULL_OCC # skip if states are occupied
        continue
    end
    for ist in 1:Nstates
        #
        if ist == jst
            continue
        end
        #
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
end
  
```

Diagonal:

```julia
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
```


Calculating epsilon (diagonal)

```julia

# Do the following for all k-points

calc_dipole_matrix!( Ham, psiks, ik, M_aux; metal_like=metal_like )
for i in 1:length(M)
    M[i] = real( M_aux[i] * conj(M_aux[i]) )
end
        
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
        
        # transition energy
        ΔE = ( evals[jst,ik] - evals[ist,ik] ) * Ha2eV + shift
        
        # loop over frequencies
        for iw in 1:Nw
            w = wgrid[iw]
            #
            ddw = (ΔE^2 - w^2)
            ddw2 = ddw^2
            denum = ( ddw2 + intersmear^2 * w^2 )* ΔE
            ff = Ha2eV^3*Focc[ist,ik]
            for i in 1:3
                εi[i,iw] = εi[i,iw] + M[i,ist,jst]*ff*intersmear*w/denum   
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
            intersmear = 1
            for i in 1:3
                εi[i,iw] = εi[i,iw] + mmf * df * intrasmear * w / denum
                εr[i,iw] = εr[i,iw] - mmf * df * w^2 / denum
            end
        end # iw
    end
end

```