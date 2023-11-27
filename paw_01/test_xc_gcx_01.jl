using PWDFT

arho = 5.0
gradx = [2.0, 3.0, 4.0]
grhoe2 = gradx[1]^2 + gradx[2]^2 + gradx[3]^2

sx, v1x, v2x = PWDFT.XC_x_pbe( arho, grhoe2 )
sc, v1c, v2c = PWDFT.XC_c_pbe( arho, grhoe2 )

