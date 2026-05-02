# Performance Benchmarks

`Ephemerides.jl` now ships runnable benchmarking and profiling entrypoints under `benchmark\` and `profiling\`, so changes on this fork can be measured against the same DE440 and PA440 workloads over time.

## Running the benchmark suite

From the repository root:

```julia
julia --startup-file=no benchmark\runbenchmarks.jl
```

For CPU profiling of broader sweeps:

```julia
julia --startup-file=no profiling\runprofiles.jl
```

The benchmark runner activates the benchmark workspace, instantiates it from the local checkout, and downloads the official DE440 and Moon PA440 kernels to `benchmark\data\` when `EPHEMERIDES_DE440_PATH` and `EPHEMERIDES_PA440_PATH` are not already set.

Useful environment variables:

1. `EPHEMERIDES_DE440_PATH`: reuse an existing DE440 kernel instead of downloading one.
2. `EPHEMERIDES_PA440_PATH`: reuse an existing Moon PA440 kernel instead of downloading one.
3. `EPHEMERIDES_BENCH_SECONDS`: change the BenchmarkTools trial time per benchmark.
4. `EPHEMERIDES_BENCH_RESULTS`: choose where the JSON results file is written.
5. `EPHEMERIDES_PROFILE_RESULTS_DIR`: choose where the profiling text reports are written.

The default suite covers:

1. Provider construction from DE440 alone and from the combined DE440 + PA440 kernel set.
2. Representative `ephem_vector3`, `ephem_vector6`, `ephem_vector12`, `ephem_rotation3`, and `ephem_rotation6` queries.
3. Sequential and shuffled Earth-track sweeps to expose cache-friendly and cache-unfriendly access patterns.
4. A broader multi-body DE440 position sweep.
5. Sequential, shuffled, and derivative PA440 orientation sweeps using `moon_pa_de440_200625.bpc` (axes `31008`).

Benchmark results are saved to `benchmark\results\latest.json`, and flat CPU profile reports are written under `profiling\results\`.

## Historical package benchmarks

The historical comparisons below were measured against CALCEPH and SPICE. They remain useful context for the package's original performance goals, but the new `benchmark\` workflow is the recommended way to benchmark this fork.

```@raw html
<p align="center">
<img src="https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl/assets/85893254/e54ac790-7421-47ff-9b68-d35bdea74de5" width="512"/>
</p>
```

```@raw html
<p align="center">
<img src="https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl/assets/85893254/ec2df247-9dde-44e9-82f4-4e8372317b8e" width="512"/>
</p>
```
