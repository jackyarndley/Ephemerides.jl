using BenchmarkTools
using ChainRulesCore: NoTangent, rrule
using DifferentiationInterface:
    AutoFiniteDiff,
    AutoForwardDiff,
    AutoZygote,
    derivative
import FiniteDiff
import ForwardDiff
import Zygote

function build_ad_benchmark_suite()
    provider, _, _ = build_reference_provider()
    time = 86400.0
    timespan = extrema(benchmark_position_times())
    route = prepare_ephemeris(provider, 3, 301; timespan)
    routes = map(
        pair -> prepare_ephemeris(provider, pair...; timespan),
        POSITION_PAIRS,
    )

    position = epoch -> ephem_vector3(route, epoch)
    state6 = epoch -> ephem_vector6(route, epoch)
    state9 = epoch -> ephem_vector9(route, epoch)
    scalar_batch = epoch -> sum(sum, ephem_vector3(routes, epoch))

    forward_backend = AutoForwardDiff()
    finite_backend = AutoFiniteDiff(fdtype=Val(:central))
    reverse_backend = AutoZygote()

    _, pullback = rrule(ephem_vector3, route, time)
    output_tangent = Ref((1.0, 1.0, 1.0))

    suite = BenchmarkGroup()
    suite["analytic"] = BenchmarkGroup()
    suite["forward"] = BenchmarkGroup()
    suite["finite"] = BenchmarkGroup()
    suite["reverse"] = BenchmarkGroup()
    suite["chainrules"] = BenchmarkGroup()
    suite["third-body"] = BenchmarkGroup()

    suite["analytic"]["vector6 state"] =
        @benchmarkable ephem_vector6($route, $time)
    suite["analytic"]["vector9 state"] =
        @benchmarkable ephem_vector9($route, $time)
    suite["forward"]["vector3"] =
        @benchmarkable derivative($position, $forward_backend, $time)
    suite["forward"]["vector6"] =
        @benchmarkable derivative($state6, $forward_backend, $time)
    suite["forward"]["vector9"] =
        @benchmarkable derivative($state9, $forward_backend, $time)
    suite["finite"]["vector3"] =
        @benchmarkable derivative($position, $finite_backend, $time)
    suite["reverse"]["vector3"] =
        @benchmarkable derivative($position, $reverse_backend, $time)
    suite["reverse"]["vector6"] =
        @benchmarkable derivative($state6, $reverse_backend, $time)
    suite["chainrules"]["rrule construction"] =
        @benchmarkable rrule(ephem_vector3, $route, $time)
    suite["chainrules"]["reused pullback"] =
        @benchmarkable $pullback($output_tangent[])
    suite["third-body"]["analytic vector6 batch"] =
        @benchmarkable ephem_vector6($routes, $time)
    suite["third-body"]["forward scalar derivative"] =
        @benchmarkable derivative($scalar_batch, $forward_backend, $time)
    suite["third-body"]["reverse scalar derivative"] =
        @benchmarkable derivative($scalar_batch, $reverse_backend, $time)

    return suite
end

function run_ad_benchmarks()
    suite = build_ad_benchmark_suite()
    tune!(suite)
    trial_seconds = parse(
        Float64, get(ENV, "EPHEMERIDES_BENCH_SECONDS", "1.0")
    )
    results = run(suite; seconds=trial_seconds)
    print_benchmark_summary(results)
    return results
end
