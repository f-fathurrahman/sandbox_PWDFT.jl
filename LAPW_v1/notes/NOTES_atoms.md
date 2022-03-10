rdiracint

Integrates the radial Dirac equation from $r=0$ outwards. This involves using the predictor-corrector method to solve the coupled first-order equations (in atomic units)

$$
\begin{align*}\left(\frac{\mathrm{d}}{\mathrm{d}r}+\frac{\kappa}{r}\right)G_{\kappa} & =\frac{1}{c}\left(2E_{0}+E-V\right)F_{\kappa}\\
\left(\frac{\mathrm{d}}{\mathrm{d}r}-\frac{\kappa}{r}\right)F_{\kappa} & =-\frac{1}{c}\left(E-V\right)G_{\kappa}
\end{align*}
$$
where $G_{\kappa}=rg_{\kappa}$ and $F_{\kappa}=rf_{\kappa}$ are the major and minor components multiplied by $r$, respectively; $V$ is the external potential; $E_{0}$ is the electron rest energy; $E$ is the eigen energy (excluding $E_{0}$); and $\kappa=l$ for $j=l-\frac{1}{2}$ or $\kappa=-(l+1)$ for $j=l+\frac{1}{2}$.

rschrodint

Integrates the scalar relativistic radial Schrodinger equation from $r=0$ outwards. This involves using the predictor-corrector method to solve the coupled first-order equations (in atomic units)
$$
\begin{align*}\frac{\mathrm{d}}{\mathrm{d}r}P_{l} & =2MQ_{l}+\frac{1}{r}P_{l}\\
\frac{\mathrm{d}}{\mathrm{d}r}Q_{l} & =-\frac{1}{r}Q_{l}+\left[\frac{l(l+1)}{2Mr^{2}}+(V-E)\right]P_{l}
\end{align*}
$$
where $V$ is the external potential, $E$ is the eigen energy and $M=1+\dfrac{(E-V)}{2c^{2}}$.
Following the convention of Koelling and Harmon, J. Phys. C: Solid State Phys. 10 3107 (1977),
the functions $P_{l}$ and $Q_{l}$ are defined by
$$
\begin{align*}P_{l} & =rg_{l}\\
Q_{l} & =\frac{r}{2M}\frac{\mathrm{d}g_{l}}{\mathrm{d}r}
\end{align*}
$$
where $g_{l}$ is the major component of the Dirac equation (see the routine rdiracint).

