import LightXML
using PWDFT: PsPot_UPF

function _read_xml_str_vector!(root_elem, tagname, Nr, f)
    # Get the element by its tagname
    pp_f = LightXML.get_elements_by_tagname(root_elem, tagname)[1]
    # Get the content (as string)
    pp_f_str = LightXML.content(pp_f)
    # Replace new line with space (to be parsed)
    pp_f_str = replace(pp_f_str, "\n" => " ")
    # Parse into an array
    spl_str = split(pp_f_str, keepempty=false)
    for i in 1:Nr
        f[i] = parse(eltype(f), spl_str[i])
    end
    return
end

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
    #
    _read_paw_data(xroot, Nr, lmax)
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
