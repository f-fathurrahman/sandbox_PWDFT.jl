using DelimitedFiles

function read_complex_matrix(fname_re, fname_im)
    dat_re = readdlm(fname_re)
    dat_im = readdlm(fname_im)
    Ndata = length(dat_re)
    @assert Ndata == length(dat_im)
    _N = round(Int64, sqrt(Ndata))
    mat = reshape(dat_re, (_N,_N)) + im*reshape(dat_im, (_N,_N))
    return mat
end

Hc_ref = read_complex_matrix("fort.100", "fort.101")
Sc_ref = read_complex_matrix("fort.102", "fort.103")
vc_ref = read_complex_matrix("fort.104", "fort.105")

Nstates = 3
#ew_ref, vc_ref = eigen(Hc[1:Nstates,1:Nstates], Sc[1:Nstates,1:Nstates]);