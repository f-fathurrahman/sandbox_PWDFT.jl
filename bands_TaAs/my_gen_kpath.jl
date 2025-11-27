function my_gen_kpath(
    atoms::Atoms, path_str_dash::String; Δk = 0.02
)
    # manually edited from ASE
    dict_spec_kpts = Dict(
        "G" => [0., 0., 0.],
        "N" => [0. , 0.5, 0. ],
        "P" => [-0.5 ,  0.5 ,  0.25],
        "S" => [ 0.        ,  0.54370934, -0.27185467],
        "S1" => [0.        , 0.45629066, 0.27185467],
        "X" => [-0.5,  0.5,  0. ],
        "Y" => [-0.45629066,  0.54370934, -0.04370934],
        "Y1" => [-0.45629066,  0.45629066,  0.5       ],
        "Z" => [0. , 0. , 0.5]
    )

    path_str = split(path_str_dash,"-",keepempty=false)
    println("path_str = ", path_str)
    Nkpt_spec = length(path_str)
    kpt_spec = zeros(3,Nkpt_spec)
    kpt_spec_labels = path_str
    for ik in 1:Nkpt_spec
        label = path_str[ik]
        kvec = dict_spec_kpts[label]
        kpt_spec[:,ik] = kvec
        println("kpt_spec[:,ik] = ", kpt_spec[:,ik])
    end

    # distance between two adjacent special points in path_str
    d = zeros(Nkpt_spec-1)
    for ik = 1:Nkpt_spec-1
        d[ik] = norm( kpt_spec[:,ik+1] - kpt_spec[:,ik] )
    end

    # estimate number of kpoints between two adjacent special points in kpath
    # based on Δk
    Nk = zeros(Int64,Nkpt_spec-1)
    for ik = 1:Nkpt_spec-1
        Nk[ik] = round(Int64, d[ik]/Δk)
        println("Nk[ik] = ", Nk[ik])
    end

    ipk = 0
    Nkpt_on_path = sum(Nk .- 1) + 1 # total no. of kpt on the path
    kpt = zeros(3,Nkpt_on_path)
    for ik = 1:Nkpt_spec-1
        dvec = kpt_spec[:,ik+1] - kpt_spec[:,ik]
        dk = dvec./(Nk[ik]-1)
        for iik = 1:Nk[ik]-1
            kvec = kpt_spec[:,ik] + (iik-1).*dk
            ipk = ipk + 1
            kpt[:,ipk] = kvec
        end
    end
    println("Nkpt_on_path = ", Nkpt_on_path)

    # The last point
    kpt[:,Nkpt_on_path] = kpt_spec[:,Nkpt_spec]

    # convert to Cartesian form
    RecVecs = 2*pi*inv(atoms.LatVecs')
    kpt_cart = RecVecs*kpt
    kpt_spec_cart = RecVecs*kpt_spec
    #
    wk = ones(Nkpt_on_path) # not used for non-scf calculations
    #
    return KPoints(Nkpt_on_path, kpt_cart, wk, RecVecs), kpt_spec_cart, kpt_spec_labels

end
