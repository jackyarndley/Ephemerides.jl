using BenchmarkTools

const DEFAULT_RESULTS_PATH = joinpath(@__DIR__, "results", "latest.json")
include("common.jl")

function build_benchmark_suite()
    eph, de440_path, pa440_path = build_reference_provider()
    t0 = 0.0

    position_times = benchmark_position_times()
    random_position_times = shuffled_copy(position_times)
    rotation_times = benchmark_rotation_times()
    random_rotation_times = shuffled_copy(rotation_times)
    earth_route = prepare_ephemeris(
        eph, 3, 399; timespan=extrema(position_times)
    )
    third_body_routes = map(
        pair -> prepare_ephemeris(eph, pair...; timespan=extrema(position_times)),
        POSITION_PAIRS,
    )
    pa440_route = prepare_orientation(
        eph, 1, PA440_AXES_ID; timespan=extrema(rotation_times)
    )

    suite = BenchmarkGroup()
    suite["load"] = BenchmarkGroup()
    suite["queries"] = BenchmarkGroup()
    suite["prepared"] = BenchmarkGroup()
    suite["positions"] = BenchmarkGroup()
    suite["rotations"] = BenchmarkGroup()

    suite["load"]["EphemerisProvider(de440)"] = @benchmarkable EphemerisProvider($de440_path)
    suite["load"]["EphemerisProvider(de440 + pa440)"] =
        @benchmarkable EphemerisProvider([$de440_path, $pa440_path])

    suite["queries"]["ephem_vector3(3, 399, t0)"] =
        @benchmarkable ephem_vector3($eph, 3, 399, $t0)
    suite["queries"]["ephem_vector6(3, 301, t0)"] =
        @benchmarkable ephem_vector6($eph, 3, 301, $t0)
    suite["queries"]["ephem_vector12(2, 299, t0)"] =
        @benchmarkable ephem_vector12($eph, 2, 299, $t0)
    suite["queries"]["ephem_rotation3(1, PA440, t0)"] =
        @benchmarkable ephem_rotation3($eph, 1, $PA440_AXES_ID, $t0)
    suite["queries"]["ephem_rotation6(1, PA440, t0)"] =
        @benchmarkable ephem_rotation6($eph, 1, $PA440_AXES_ID, $t0)

    suite["prepared"]["ephem_vector3(earth_route, t0)"] =
        @benchmarkable ephem_vector3($earth_route, $t0)
    suite["prepared"]["ephem_rotation3(pa440_route, t0)"] =
        @benchmarkable ephem_rotation3($pa440_route, $t0)
    suite["prepared"]["third-body-position-sweep"] = @benchmarkable begin
        acc = 0.0
        for t in $position_times
            acc += sum(state -> state[1], ephem_vector3($third_body_routes, t))
        end
        acc
    end

    suite["positions"]["sequential-earth-track"] = @benchmarkable begin
        acc = 0.0
        for t in $position_times
            acc += ephem_vector3($eph, 3, 399, t)[1]
        end
        acc
    end

    suite["positions"]["random-earth-track"] = @benchmarkable begin
        acc = 0.0
        for t in $random_position_times
            acc += ephem_vector3($eph, 3, 399, t)[1]
        end
        acc
    end

    suite["positions"]["multi-body-position-sweep"] = @benchmarkable begin
        acc = 0.0
        for t in $position_times
            for (from, to) in $POSITION_PAIRS
                acc += ephem_vector3($eph, from, to, t)[1]
            end
        end
        acc
    end

    suite["rotations"]["sequential-pa440-axes"] = @benchmarkable begin
        acc = 0.0
        for t in $rotation_times
            acc += ephem_rotation3($eph, 1, $PA440_AXES_ID, t)[1]
        end
        acc
    end

    suite["rotations"]["random-pa440-axes"] = @benchmarkable begin
        acc = 0.0
        for t in $random_rotation_times
            acc += ephem_rotation3($eph, 1, $PA440_AXES_ID, t)[1]
        end
        acc
    end

    suite["rotations"]["pa440-orientation-sweep"] = @benchmarkable begin
        acc = 0.0
        for t in $rotation_times
            acc += ephem_rotation6($eph, 1, $PA440_AXES_ID, t)[4]
        end
        acc
    end

    return suite, de440_path, pa440_path
end

function print_benchmark_summary(group::BenchmarkGroup, prefix::String = "")
    for key in sort!(collect(keys(group)); by = string)
        label = isempty(prefix) ? string(key) : prefix * " / " * string(key)
        value = group[key]

        if value isa BenchmarkGroup
            print_benchmark_summary(value, label)
            continue
        end

        estimate = BenchmarkTools.median(value)
        println(
            rpad(label, 44),
            BenchmarkTools.prettytime(estimate.time),
            " | ",
            BenchmarkTools.prettymemory(estimate.memory),
        )
    end
end

function run_benchmarks(; results_path::AbstractString = get(ENV, "EPHEMERIDES_BENCH_RESULTS", DEFAULT_RESULTS_PATH))
    suite, de440_path, pa440_path = build_benchmark_suite()
    tune!(suite)

    trial_seconds = parse(Float64, get(ENV, "EPHEMERIDES_BENCH_SECONDS", "1.0"))
    results = run(suite; seconds = trial_seconds)

    mkpath(dirname(results_path))
    BenchmarkTools.save(results_path, results)

    println("DE440 kernel: $(de440_path)")
    println("PA440 kernel: $(pa440_path)")
    println("Saved results to $(results_path)")
    print_benchmark_summary(results)

    return results
end
