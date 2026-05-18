"""
Serialize some Elk's global variables using Julia's `serialize`.
"""
function serialize_variables()

    var_list = [
        :atposl, :atposc,
        :rhomt, :rhoir, :rhosp,
        :vclmt, :vclir,
        :exmt, :ecmt, :vxcmt,
        :exir, :ecir, :vxcir,
        :apwfr, :apwdfr,
        :haa, :hloa, :hlolo,
        :oalo, :ololo
    ]

    for v in var_list
        println("Write variable $(v)")
        str_prog = """
        serialize(joinpath(ELK_DATA_DIR, "$(v).jldat"), get_$(v)())
        """
        eval(Meta.parse(str_prog))
    end

    return
end

# This is a separate function because apwalm is calculated on the fly
# for each k-points
function serialize_apwalm()
    nkpt = get_nkpt()
    for ik in 1:nkpt
        apwalm = get_apwalm(ik)
        println("Writing apwalm for ik = ", ik)
        serialize(joinpath(ELK_DATA_DIR, "apwalm_$(ik).jldat"), apwalm) 
    end
    return
end

function serialize_evalfv_evecfv()
    nkpt = get_nkpt()
    for ik in 1:nkpt
        evalfv = get_evalfv_from_file(ik)
        println("Writing evalfv for ik = ", ik)
        serialize(joinpath(ELK_DATA_DIR, "evalfv_$(ik).jldat"), evalfv) 
        #
        evecfv = get_evecfv_from_file(ik)
        println("Writing evecfv for ik = ", ik)
        serialize(joinpath(ELK_DATA_DIR, "evecfv_$(ik).jldat"), evecfv) 
    end
    return
end
