const _SHT_arrays = [
    :rbshti, :rbshto, :rfshti, :rfshto,
    :zbshti, :zbshto, :zfshti, :zfshto
]

const _SHT_arrays_type = [
    Float64, Float64, Float64, Float64,
    ComplexF64, ComplexF64, ComplexF64, ComplexF64
]

const _SHT_arrays_dims = [
    :lmmaxi, :lmmaxo, :lmmaxi, :lmmaxo,
    :lmmaxi, :lmmaxo, :lmmaxi, :lmmaxo
]

for (s,t,d) in zip(_SHT_arrays, _SHT_arrays_type, _SHT_arrays_dims)
    str_prog = """
    function get_$s()
        symbol = :__m_sht_MOD_$s
        $d = get_$d()
        return _load_allocatable_array(symbol, $t, ($d,$d))
    end
    """
    eval(Meta.parse(str_prog))
end
