import LightXML
using PWDFT: PsPot_UPF

include("PAWData_UPF.jl")
include("read_paw_data.jl")

function test_paw(upf_file)
    xdoc = LightXML.parse_file(upf_file)
    # get the root element
    xroot = LightXML.root(xdoc)  # an instance of XMLElement
    #
    # Read some information from header
    #
    pp_header = LightXML.get_elements_by_tagname(xroot, "PP_HEADER")
    #atsymb = LightXML.attributes_dict(pp_header[1])["element"]
    #zval = Int64( parse( Float64, LightXML.attributes_dict(pp_header[1])["z_valence"] ) )
    lmax = parse( Int64, LightXML.attributes_dict(pp_header[1])["l_max"] )
    Nr = parse(Int64,LightXML.attributes_dict(pp_header[1])["mesh_size"])
    Nproj = parse(Int64,LightXML.attributes_dict(pp_header[1])["number_of_proj"])
    #
    _read_paw_data(xroot, Nr, lmax, Nproj)
    return
end

function test_main()
    prefix_dir = "/home/efefer/pseudo/PSLIB/"
    files = readlines("PSLIB_PAW")
    for f in files
        filepath = joinpath(prefix_dir, f)
        println("\nfilepath = ", filepath)
        test_paw(filepath)
        #
        #psp = PsPot_UPF(joinpath(prefix_dir, f))
        #println("psp.nqf = ", psp.nqf)
        #println("psp.nqlc = ", psp.nqlc)
        #println("psp.q_with_l = ", psp.q_with_l)
        #println("psp.is_paw = ", psp.is_paw)
    end
end

test_main()
