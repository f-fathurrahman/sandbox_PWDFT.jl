using PWDFT

Rhoe = [1.1 ;; 2.1]

xc_calc = LibxcXCCalculator(x_id=1, c_id=12)
epsxc = calc_epsxc_VWN(xc_calc, Rhoe)
Vxc = calc_Vxc_VWN(xc_calc, Rhoe)

println()
println("From Libxc")
println("epsxc = ", epsxc)
println("Vxc = ", Vxc)

ρ = Rhoe[1,1] + Rhoe[1,2]
ζ = ( Rhoe[1,1] - Rhoe[1,2] ) / ρ
x_res = PWDFT.XC_x_slater_spin(ρ, ζ)
c_res = PWDFT.XC_c_pw_spin(ρ, ζ)
xc_res = x_res .+ c_res # the results are epsxc, Vxc_up, Vxc_dn

println()
println("Internal result")
println(xc_res)

println()
println("Diff epsxc = ", abs(xc_res[1] - epsxc[1]))
println("Diff Vxc up = ", abs(xc_res[2] - Vxc[1]))
println("Diff Vxc dn = ", abs(xc_res[3] - Vxc[2]))


