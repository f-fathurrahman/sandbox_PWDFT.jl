function findband!(l, r, Vr, E;
    epsband=1e-12, demax=2.5, de0=0.001, maxstp=250
)

#=
  IMPLICIT NONE
  ! arguments
  REAL(8), INTENT(in) :: sol
  INTEGER, INTENT(in) :: l,nr
  REAL(8), INTENT(in) :: r(nr),vr(nr)
  REAL(8), INTENT(in) :: eps,demax
  REAL(8), INTENT(inout) :: e
  LOGICAL, INTENT(out) :: fnd
  ! local variables
  LOGICAL :: ft,fb
  ! maximum number of steps
  INTEGER, PARAMETER :: maxstp=250
  INTEGER :: ip,ie,nn
  ! initial step size
  REAL(8), PARAMETER :: de0=0.001d0
  REAL(8) :: de,et,eb,t,tp
  ! automatic arrays
  REAL(8) :: p0(nr),p1(nr),q0(nr),q1(nr)
=#

    ft = false
    fb = false
    fnd = false
    et = E
    eb = E

    nr = size(Vr, 1)
    p0 = zeros(Float64, nr)
    p1 = zeros(Float64, nr)
    q0 = zeros(Float64, nr)
    q1 = zeros(Float64, nr)
    # two-pass loop
    for ip in 1:2
        # find the top of the band
        tp = 0.0
        de = de0
        for ie in 1:maxstp
            et += de # increment de
            if E < 0.0
                if et > E + demax
                    break
                end
            end
            nn, et = rschrodint!(l, et, r, Vr, p0, p1, q0, q1)
            t = p0[nr]
            if ie > 1
                if t*tp <= 0.0
                    if abs(de) < epsband
                        if fb
                            @goto LABEL10
                        end
                        ft = true
                        eb = et + 5.0*abs(de0)
                        break
                    end
                    de = -0.5*de
                else
                    de = 1.5*de
                end
            end
            tp = t
        end # maxstp
        #
        if fb
            return fnd, E
        end
        # find the bottom of the band
        tp = 0.0
        de = -de0
        for ie in 1:maxstp
            eb += de
            if eb < E-demax
                return fnd, E
            end
            nn, eb = rschrodint!(l, eb, r, Vr, p0, p1, q0, q1)
            t = p1[nr]
            if ie > 1
                if t*tp <= 0.0
                    if abs(de) < epsband
                        if ft
                            @goto LABEL10
                        end
                        fb = true
                        et = eb - 5.0*abs(de0)
                        break 
                    end
                    de = -0.5*de
                else 
                    de = 1.5*de
                end
            end
            tp = t
        end # do maxstp
    end # do ip
    return fnd, E
    #
    @label LABEL10  # continue
    # set the band energy halfway between top and bottom
    E = (et + eb)/2.0
    fnd = true
    #
    return fnd, E
end

