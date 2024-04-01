function write_elk_variables()
    rhomt = elk_get_rhomt()
    serialize(joinpath(ELK_DATA_DIR, "rhomt.dat"), rhomt)
    return
end

