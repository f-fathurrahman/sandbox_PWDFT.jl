using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("INC_gamma_only.jl")

function do_calc(;gamma_only=true)

    Random.seed!(1234)

    # In Bohr
    A = 18.65
    B = 18.65
    C = 18.65

    atoms = Atoms(xyz_string_frac="""
    96

    O      0.3342   0.3858   0.1702
    O      0.6681   0.0772   0.0996
    O      0.4894   0.2752   0.9664
    O      0.8374   0.0165   0.8885
    O      0.4164   0.1406   0.5374
    O      0.9200   0.2746   0.7479
    O      0.9298   0.6759   0.0546
    O      0.5197   0.5773   0.2470
    O      0.8453   0.2308   0.2531
    O      0.7080   0.4512   0.1102
    O      0.6297   0.6739   0.4697
    O      0.6889   0.4428   0.8100
    O      0.9273   0.8859   0.2748
    O      0.7260   0.9572   0.4514
    O      0.9539   0.6297   0.6304
    O      0.0634   0.4197   0.2665
    O      0.0088   0.1409   0.5073
    O      0.8042   0.4375   0.4942
    O      0.6981   0.1893   0.5833
    O      0.2908   0.5198   0.5234
    O      0.1835   0.3971   0.7808
    O      0.2450   0.2024   0.3474
    O      0.3062   0.1574   0.7995
    O      0.4465   0.8389   0.2161
    O      0.1931   0.5738   0.0049
    O      0.1879   0.9665   0.1899
    O      0.0876   0.6667   0.4085
    O      0.2637   0.7476   0.7879
    O      0.5203   0.6329   0.7214
    O      0.3714   0.8763   0.4899
    O      0.7915   0.8290   0.7054
    O      0.1436   0.9899   0.9297
    H      0.3742   0.3360   0.0929
    H      0.3150   0.3226   0.2472
    H      0.9405   0.1804   0.4516
    H      0.6062   0.0086   0.1372
    H      0.7231   0.0354   0.0264
    H      0.5460   0.2002   0.9884
    H      0.5502   0.3269   0.9029
    H      0.8537   0.1049   0.8480
    H      0.4251   0.0432   0.5210
    H      0.3722   0.1543   0.6275
    H      0.7041   0.0486   0.4909
    H      0.9378   0.2191   0.6620
    H      0.9083   0.7339   0.1314
    H      0.5940   0.5403   0.1919
    H      0.4413   0.5140   0.2303
    H      0.7975   0.1623   0.1939
    H      0.9351   0.2490   0.2169
    H      0.7495   0.3794   0.1764
    H      0.7769   0.5261   0.0949
    H      0.7242   0.6488   0.4514
    H      0.5812   0.6156   0.4033
    H      0.7610   0.3991   0.7586
    H      0.7253   0.4549   0.9012
    H      0.9805   0.8118   0.3158
    H      0.8774   0.9234   0.3508
    H      0.8266   0.9432   0.8185
    H      0.6439   0.9054   0.4658
    H      0.9099   0.7082   0.6739
    H      0.1146   0.3475   0.3129
    H      0.1177   0.4592   0.1936
    H      0.0984   0.1864   0.4941
    H      0.8732   0.4903   0.5432
    H      0.8466   0.4047   0.4105
    H      0.5983   0.1849   0.5758
    H      0.7255   0.2866   0.5619
    H      0.2468   0.6010   0.4948
    H      0.3655   0.5512   0.5766
    H      0.2218   0.4478   0.6984
    H      0.2272   0.3027   0.7736
    H      0.0089   0.3172   0.7749
    H      0.3206   0.1756   0.4126
    H      0.2308   0.1250   0.2800
    H      0.2420   0.1027   0.8550
    H      0.3752   0.1963   0.8703
    H      0.3477   0.8559   0.1901
    H      0.4738   0.7418   0.2103
    H      0.0297   0.6626   0.0605
    H      0.2668   0.5360   0.0599
    H      0.1691   0.9883   0.0958
    H      0.0952   0.9378   0.2237
    H      0.0161   0.6488   0.4768
    H      0.0940   0.5739   0.3729
    H      0.2728   0.6995   0.8762
    H      0.3551   0.7479   0.7517
    H      0.5812   0.5576   0.7579
    H      0.5605   0.6617   0.6349
    H      0.1741   0.5139   0.9284
    H      0.3026   0.8058   0.5084
    H      0.4067   0.8573   0.3940
    H      0.7040   0.7965   0.7365
    H      0.7716   0.8608   0.6097
    H      0.0148   0.6043   0.7064
    H      0.1834   0.9112   0.8662
    H      0.0492   0.9894   0.9207
    """, in_bohr=true, LatVecs=diagm([A,B,C]))
    pspfiles = get_default_psp(atoms)

    write_xsf("TEMP_32H2O.xsf", atoms)
    
    ecutwfc = 15.0
    
    Ham = HamiltonianGamma(atoms, pspfiles, ecutwfc)
    psis = randn_BlochWavefuncGamma(Ham)
    
    if gamma_only
        @time KS_solve_Emin_PCG_dot!( Ham, psis, NiterMax=200 )
    else
        Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
        println(Ham_)
        psiks = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis )
        @time KS_solve_Emin_PCG_dot!( Ham_, psiks, startingrhoe=:random,
            skip_initial_diag=true, NiterMax=200 )
    end

end


do_calc(gamma_only=true)