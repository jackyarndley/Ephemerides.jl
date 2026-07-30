using DifferentiationInterface
using ChainRulesCore: NoTangent, frule
using FiniteDiff: FiniteDiff
using ForwardDiff: ForwardDiff
using Zygote: Zygote

@testset "DifferentiationInterface compatibility" begin
    test_data = artifact"testdata"
    point_provider = EphemerisProvider(joinpath(test_data, "spk2_ex1.bsp"))
    orientation_provider = EphemerisProvider(joinpath(test_data, "pa421.bpc"))
    time = 86400.0

    point_route = prepare_ephemeris(
        point_provider, 3, 301; timespan=(time - 1, time + 1)
    )
    orientation_route = prepare_orientation(
        orientation_provider, 1, 31006; timespan=(time - 1, time + 1)
    )

    forward_backend = AutoForwardDiff()
    finite_backend = AutoFiniteDiff(fdtype=Val(:central))
    reverse_backend = AutoZygote()

    position = epoch -> ephem_vector3(point_route, epoch)
    state6 = epoch -> ephem_vector6(point_route, epoch)
    state9 = epoch -> ephem_vector9(point_route, epoch)
    angles = epoch -> ephem_rotation3(orientation_route, epoch)

    expected_velocity = ephem_vector6(point_route, time)[4:6]
    expected_acceleration = ephem_vector9(point_route, time)[7:9]
    expected_jerk = ephem_vector12(point_route, time)[10:12]
    expected_angle_rate = ephem_rotation6(orientation_route, time)[4:6]

    rule_value, rule_tangent = frule(
        (NoTangent(), NoTangent(), 1.0), ephem_vector3, point_route, time
    )
    @test rule_value == ephem_vector3(point_route, time)
    @test rule_tangent ≈ expected_velocity atol=1e-13 rtol=1e-13

    for backend in (forward_backend, reverse_backend)
        @test derivative(position, backend, time) ≈ expected_velocity atol=1e-13 rtol=1e-13
        @test derivative(state6, backend, time) ≈
              ephem_vector9(point_route, time)[4:9] atol=1e-13 rtol=1e-13
        @test derivative(state9, backend, time) ≈
              ephem_vector12(point_route, time)[4:12] atol=1e-13 rtol=1e-13
        @test derivative(angles, backend, time) ≈ expected_angle_rate atol=1e-13 rtol=1e-13
    end

    @test derivative(position, finite_backend, time) ≈
          expected_velocity atol=1e-7 rtol=1e-7
    @test derivative(state6, finite_backend, time)[4:6] ≈
          expected_acceleration atol=1e-10 rtol=1e-6
    @test derivative(state9, finite_backend, time)[7:9] ≈
          expected_jerk atol=1e-12 rtol=1e-5
    @test derivative(angles, finite_backend, time) ≈
          expected_angle_rate atol=1e-10 rtol=1e-6

    provider_position = epoch -> ephem_vector3(point_provider, 3, 301, epoch)
    @test derivative(provider_position, reverse_backend, time) ≈
          expected_velocity atol=1e-13 rtol=1e-13

    duplicate_provider = EphemerisProvider([
        joinpath(test_data, "spk2_ex1.bsp"),
        joinpath(test_data, "spk2_ex1.bsp"),
    ])
    composite_route = prepare_ephemeris(duplicate_provider, 3, 301)
    @test derivative(
        epoch -> ephem_vector3(composite_route, epoch), reverse_backend, time
    ) ≈ expected_velocity atol=1e-13 rtol=1e-13

    barycenter_route = prepare_ephemeris(
        point_provider, 0, 10; timespan=(time - 1, time + 1)
    )
    routes = (point_route, barycenter_route)
    scalar_batch = epoch -> sum(sum, ephem_vector3(routes, epoch))
    expected_batch_rate = sum(
        sum(state[4:6]) for state in ephem_vector6(routes, time)
    )
    @test derivative(scalar_batch, reverse_backend, time) ≈
          expected_batch_rate atol=1e-13 rtol=1e-13
end
