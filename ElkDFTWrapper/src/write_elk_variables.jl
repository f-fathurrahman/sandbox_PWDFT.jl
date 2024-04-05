"""
Serialize some Elk's global variables using Julia's `serialize`.
"""
function write_elk_variables()
    # Muffin tin charge density
    rhomt = elk_get_rhomt()
    serialize(joinpath(ELK_DATA_DIR, "rhomt.dat"), rhomt)
    #
    rhoir = elk_get_rhoir()
    serialize(joinpath(ELK_DATA_DIR, "rhoir.dat"), rhoir)
    return
end

