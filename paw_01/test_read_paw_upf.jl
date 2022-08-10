import LightXML

function _read_xml_str_vector!(root_elem, tagname, Nr, f)
    pp_f = LightXML.get_elements_by_tagname(root_elem, tagname)
    pp_f_str = LightXML.content(pp_f[1])
    pp_f_str = replace(pp_f_str, "\n" => " ")
    spl_str = split(pp_f_str, keepempty=false)
    for i in 1:Nr
        f[i] = parse(eltype(f), spl_str[i])
    end
    return
end

upf_file = "/home/efefer/pseudo/PSLIB/Ni.pbe-n-kjpaw_psl.0.1.UPF"

xdoc = LightXML.parse_file(upf_file)

# get the root element
xroot = LightXML.root(xdoc)  # an instance of XMLElement

#
# Read some information from header
#
pp_header = LightXML.get_elements_by_tagname(xroot, "PP_HEADER")
atsymb = LightXML.attributes_dict(pp_header[1])["element"]
zval = Int64( parse( Float64, LightXML.attributes_dict(pp_header[1])["z_valence"] ) )
lmax = parse( Int64, LightXML.attributes_dict(pp_header[1])["l_max"] )
Nr = parse(Int64,LightXML.attributes_dict(pp_header[1])["mesh_size"])


pp_full_wfc = LightXML.get_elements_by_tagname(xroot, "PP_FULL_WFC")
Nwfc = parse(Int64, LightXML.attributes_dict(pp_full_wfc[1])["number_of_wfc"])
println("Nwfc = ", Nwfc)
# XXX Compare Nwfc with number of projectors ?

aewfc = zeros(Float64, Nr, Nwfc)
pswfc = zeros(Float64, Nr, Nwfc) # different from chi
for iwf in 1:Nwfc
    tagname = "PP_AEWFC."*string(iwf)
    @views _read_xml_str_vector!(pp_full_wfc[1], tagname, Nr, aewfc[:,iwf])
    #
    tagname = "PP_PSWFC."*string(iwf)
    @views _read_xml_str_vector!(pp_full_wfc[1], tagname, Nr, pswfc[:,iwf])
end

pp_paw = LightXML.get_elements_by_tagname(xroot, "PP_PAW")
core_energy = parse(Float64, LightXML.attributes_dict(pp_paw[1])["core_energy"])
println("core_energy = ", core_energy)

pp_occ = LightXML.get_elements_by_tagname(pp_paw[1], "PP_OCCUPATIONS")
Nocc = parse(Int64, LightXML.attributes_dict(pp_occ[1])["size"])
println("Nocc = ", Nocc)
# Nocc should be the same as Nwf?
paw_pp_occ = zeros(Float64, Nocc)
_read_xml_str_vector!(pp_paw[1], "PP_OCCUPATIONS", Nocc, paw_pp_occ)
println("paw_pp_occ = ", paw_pp_occ)

pp_ae_nlcc = zeros(Float64, Nr)
_read_xml_str_vector!(pp_paw[1], "PP_AE_NLCC", Nr, pp_ae_nlcc)
println("Done reading pp_ae_nlcc")

pp_ae_vloc = zeros(Float64, Nr)
_read_xml_str_vector!(pp_paw[1], "PP_AE_VLOC", Nr, pp_ae_vloc)
println("Done reading pp_ae_vloc")

