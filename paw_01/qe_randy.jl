mutable struct QERandom
    ir::Vector{Int64}
    iy::Int64
    idum::Int64
    first::Bool
end


function QERandom()
    ntab = 97 # hardcoded
    ir = zeros(Int64, ntab)
    iy = 0
    idum = 0
    first = true
    return QERandom(ir, iy, idum, first)
end



#
# x=randy(n): reseed with initial seed idum=n ( 0 <= n <= ic, see below)
#             if randy is not explicitly initialized, it will be
#             initialized with seed idum=0 the first time it is called
# x=randy() : generate uniform real(DP) numbers x in [0,1]
#
# From randy function in random_numbers.f90 of QE
function qe_randy!( qe_random::QERandom; irand::Union{Nothing,Int64}=nothing)
    
    m    = 714025
    ia   = 1366
    ic   = 150889
    ntab = 97

    rm = 1.0/m
    
    if irand != nothing
       qe_random.idum = min(abs(irand), ic) 
       qe_random.first = true
    end

    if qe_random.first
        qe_random.first = false
        qe_random.idum = mod( ic - qe_random.idum, m )
        for j in 1:ntab
            qe_random.idum = mod(ia*qe_random.idum + ic, m)
            qe_random.ir[j] = qe_random.idum
        end
        qe_random.idum = mod(ia*qe_random.idum + ic, m)
        qe_random.iy = qe_random.idum
    end
    j = floor(Int64, 1 + (ntab*qe_random.iy)/m)
    
    if (j > ntab) || (j <  1)
        error("j is out of range in qe_randy" * string(abs(j)+1))
    end

    qe_random.iy = qe_random.ir[j]
    res = qe_random.iy*rm
    qe_random.idum = mod(ia*qe_random.idum + ic, m)
    qe_random.ir[j] = qe_random.idum

    return res

end

qe_random = QERandom()
qe_randy!(qe_random, irand=2)