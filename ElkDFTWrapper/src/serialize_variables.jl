"""
Serialize some Elk's global variables using Julia's `serialize`.
"""
function serialize_variables()
    # Muffin tin charge density
    rhomt = get_rhomt()
    serialize(joinpath(ELK_DATA_DIR, "rhomt.dat"), rhomt)
    #
    rhoir = get_rhoir()
    serialize(joinpath(ELK_DATA_DIR, "rhoir.dat"), rhoir)
    return
end

