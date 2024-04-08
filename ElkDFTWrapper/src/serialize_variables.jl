"""
Serialize some Elk's global variables using Julia's `serialize`.
"""
function serialize_variables()

    var_list = [
        :rhomt, :rhoir, :rhosp,
        :vclmt, :vclir
    ]

    for v in var_list
        println("Write variable $(v)")
        str_prog = """
        serialize(joinpath(ELK_DATA_DIR, "$(v).dat"), get_$(v)())
        """
        eval(Meta.parse(str_prog))
    end

    return
end

