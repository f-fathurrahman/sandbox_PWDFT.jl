CSC_PARAMS_SOFT_PARA = [
    7.40, 1.120, 1.890, 0.0964,  0.0250,   2.0, 3.0, 2.431, 0.0142, 2.922
]

CSC_PARAMS_SOFT_FERRO = [
    5.24, 0.0, 1.568, 0.12856, 0.003201, 2.0, 3.0, 0.0538, 1.56e-5, 2.958
]

function f_aux(a, rs):
    return -(rs + a[5]*rs^2)*log(1 + a[8]*rs + a[9]*rs^a[10]) /
            (2*(a[1] + a[2]*rs + a[3]*rs^a[6] + a[4]*rs^a[7]))
end

function calc_exc_1d_csc(rs, z)
    return f_aux(CSC_PARAMS_SOFT_PARA, rs) + 
         ( f_aux(CSC_PARAMS_SOFT_FERRO, rs) - f_aux(CSC_PARAMS_SOFT_PARA, rs) )*z^2
end

