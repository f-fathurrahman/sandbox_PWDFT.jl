#import "@preview/jlyfish:0.1.0": *

#read-julia-output(json("test02-jlyfish.json"))

Check package status

#jl(```julia
using Pkg
Pkg.activate("v1.12", shared=true)
Pkg.status()
```)

This is some code julia code:

#jl(```julia
using InteractiveUtils
println("Hello from stdout!")
versioninfo()

using Pkg
Pkg.status()
```)

This is another code that define some variables:

#jl(```julia
a = 13.1
b = 2.1
x = rand(5)
println(x)
```)

Code again:

#{
let code1 = [
```julia
println("sin(x) = ", sin.(x))
cos.(x) |> println
```
]

code1
}

//#set text(size: 8pt, font: "JuliaMono")
//#jl(code1)
