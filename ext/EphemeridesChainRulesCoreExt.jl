module EphemeridesChainRulesCoreExt

using ChainRulesCore:
    AbstractZero,
    NoTangent,
    ZeroTangent,
    unthunk
import ChainRulesCore: frule, rrule
import Ephemerides
using Ephemerides:
    AbstractPreparedEphemerisRoute,
    AbstractPreparedOrientationRoute,
    EphemerisProvider,
    ephem_rotation3,
    ephem_rotation6,
    ephem_rotation9,
    ephem_rotation12,
    ephem_vector3,
    ephem_vector6,
    ephem_vector9,
    ephem_vector12
using StaticArraysCore: SVector

@inline function split_state(::Val{N}, state) where {N}
    value = SVector{N}(ntuple(index -> state[index], Val(N)))
    derivative = SVector{N}(ntuple(index -> state[index + 3], Val(N)))
    return value, derivative
end

@inline contract_output(::AbstractZero, derivative) = ZeroTangent()
@inline contract_output(cotangent, derivative::SVector) =
    sum(index -> cotangent[index]*derivative[index], eachindex(derivative))

@inline scale_output(::AbstractZero, derivative) = ZeroTangent()
@inline scale_output(time_tangent, derivative) = time_tangent*derivative

@inline function value_and_time_derivative(higher_order, width::Val, arguments...)
    return split_state(width, higher_order(arguments...))
end

function ephemeris_rrule(higher_order, width::Val, arguments...)
    value, time_derivative = value_and_time_derivative(
        higher_order, width, arguments...
    )
    constant_tangents = ntuple(_ -> NoTangent(), Val(length(arguments) - 1))

    function ephemeris_pullback(output_tangent)
        time_tangent = contract_output(unthunk(output_tangent), time_derivative)
        return NoTangent(), constant_tangents..., time_tangent
    end

    return value, ephemeris_pullback
end

@inline function ephemeris_frule(
    higher_order, width::Val, argument_tangents, arguments...
)
    value, time_derivative = value_and_time_derivative(
        higher_order, width, arguments...
    )
    time_tangent = unthunk(last(argument_tangents))
    return value, scale_output(time_tangent, time_derivative)
end

for (function_name, higher_order, width, route_type) in (
    (:ephem_vector3, :ephem_vector6, 3, AbstractPreparedEphemerisRoute),
    (:ephem_vector6, :ephem_vector9, 6, AbstractPreparedEphemerisRoute),
    (:ephem_vector9, :ephem_vector12, 9, AbstractPreparedEphemerisRoute),
    (:ephem_rotation3, :ephem_rotation6, 3, AbstractPreparedOrientationRoute),
    (:ephem_rotation6, :ephem_rotation9, 6, AbstractPreparedOrientationRoute),
    (:ephem_rotation9, :ephem_rotation12, 9, AbstractPreparedOrientationRoute),
)
    function_object = getfield(Ephemerides, function_name)
    higher_order_object = getfield(Ephemerides, higher_order)

    @eval begin
        rrule(
            ::typeof($function_object),
            eph::EphemerisProvider,
            from::Int,
            to::Int,
            time::Number,
        ) = ephemeris_rrule($higher_order_object, Val($width), eph, from, to, time)

        frule(
            argument_tangents,
            ::typeof($function_object),
            eph::EphemerisProvider,
            from::Int,
            to::Int,
            time::Number,
        ) = ephemeris_frule(
            $higher_order_object,
            Val($width),
            argument_tangents,
            eph,
            from,
            to,
            time,
        )

        rrule(
            ::typeof($function_object),
            route::$route_type,
            time::Number,
        ) = ephemeris_rrule($higher_order_object, Val($width), route, time)

        frule(
            argument_tangents,
            ::typeof($function_object),
            route::$route_type,
            time::Number,
        ) = ephemeris_frule(
            $higher_order_object, Val($width), argument_tangents, route, time
        )
    end
end

end
