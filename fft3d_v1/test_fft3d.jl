using LinearAlgebra
using Random
using Printf
using BenchmarkTools

using PWDFT

function test_multicolumn( ecutwfc::Float64, Nstates::Int64 )

    pw = PWGrid( ecutwfc, gen_lattice_sc(16.0) )
    Ns = pw.Ns
    Npoints = prod(Ns)

    in1 = rand(ComplexF64, Npoints, Nstates )

    println("\n\nMulti column version")
    println("Data size: ", Ns)
    println("Nstates  : ", Nstates)

    println("\nNaive version")
    @time out2 = R_to_G( Ns, in1 )
    @time in2 = G_to_R( Ns, out2 )

    println("\nUsing plan")
    @time out3 = R_to_G( pw, in1 )
    @time in3 = G_to_R( pw, out3 )

    println("\nCheck")
    println("diff out = ", sum(abs.(out3 - out2))/Npoints)
    println("diff in  = ", sum(abs.(in2 - in1))/Npoints)
    println("diff in  = ", sum(abs.(in3 - in1))/Npoints)
end

function test_multicolumn_btime( ecutwfc::Float64, Nstates::Int64 )
    pw = PWGrid( ecutwfc, gen_lattice_sc(16.0) )
    Ns = pw.Ns
    Npoints = prod(Ns)

    in1 = rand(ComplexF64, Npoints, Nstates )

    println("\n\nMulti column version")
    println("Data size: ", Ns)
    println("Nstates  : ", Nstates)

    println("\nNaive version")
    @btime R_to_G( $Ns, $in1 )
    @btime G_to_R( $Ns, $in1 )

    println("\nUsing plan")
    @btime R_to_G( $pw, $in1 )
    @btime G_to_R( $pw, $in1 )
end


function test_singlecolumn_btime( ecutwfc::Float64 )
    pw = PWGrid( ecutwfc, gen_lattice_sc(16.0) )
    Ns = pw.Ns
    Npoints = prod(Ns)

    in1 = rand(ComplexF64, Npoints)
    in1r = reshape(in1,Ns)

    println("\n\nSingle column version")
    println("Data size: ", Ns)

    println("\nUsing plan with reshape previously")
    @btime R_to_G( $pw, $in1r )
    @btime G_to_R( $pw, $in1r )

    println("\nUsing plan with reshape previously (in-place)")
    @btime R_to_G!( $pw, $in1r )
    @btime G_to_R!( $pw, $in1r )

    println("\nNaive version")
    @btime R_to_G( $Ns, $in1 )
    @btime G_to_R( $Ns, $in1 )

    println("\nUsing plan")
    @btime R_to_G( $pw, $in1 )
    @btime G_to_R( $pw, $in1 )

    println("\nUsing plan (in place)")
    @btime R_to_G!( $pw, $in1 )
    @btime G_to_R!( $pw, $in1 )

end

function test_singlecolumn( ecutwfc::Float64 )

    pw = PWGrid( ecutwfc, gen_lattice_sc(16.0) )
    Ns = pw.Ns
    Npoints = prod(Ns)

    in1 = rand(ComplexF64, Npoints)

    println("\n\nSingle column version")
    println("Data size:", Ns)

    println("\nNaive version")
    @time out2 = R_to_G( Ns, in1 )
    @time in2 = G_to_R( Ns, out2 )

    println("\nUsing plan")
    @time out3 = R_to_G( pw, in1 )
    @time in3 = G_to_R( pw, out3 )

    println("\nCheck")
    println("diff out = ", sum(abs.(out3 - out2))/Npoints) 
    println("diff in  = ", sum(abs.(in2 - in1))/Npoints)
    println("diff in  = ", sum(abs.(in3 - in1))/Npoints)

end

#test_singlecolumn( 15.0 )
#test_singlecolumn( 20.0 )
#test_singlecolumn( 25.0 )

#test_multicolumn( 15.0, 10 )
#test_multicolumn( 20.0, 10 )
#test_multicolumn( 25.0, 10 )

#test_multicolumn_btime(25.0, 10)

test_singlecolumn_btime(25.0)