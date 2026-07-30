import Pkg

Pkg.activate(@__DIR__)
Pkg.resolve()
Pkg.instantiate()

using Ephemerides

include("benchmarks.jl")
include("adbenchmarks.jl")

run_ad_benchmarks()
