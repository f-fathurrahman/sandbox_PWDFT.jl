# used in _gen_sht_matrices
# convert (theta, phi) angles to vectors
# tp: theta, phi arrays
function sctovec!(Ndata, tp, v)
    # XXX Check data Ndata <= size(tp,2)
    # size of size tp is (2,Ndata)
    for i in 1:Ndata
        t1     = sin(tp[1,i])
        v[1,i] = t1*cos(tp[2,i])
        v[2,i] = t1*sin(tp[2,i])
        v[3,i] = cos(tp[1,i])
    end
    return
end