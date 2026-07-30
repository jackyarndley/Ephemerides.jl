# Reading Ephemeris Data

This tutorials will walk you through the basic features and interfaces that allow you to compute translation and orientation data from binary ephemeris kernels.

## Computing state vectors

`Ephemerides.jl` allows the computation of a relative position between two points and its higher order derivatives up to order 3 (i.e., velocity, acceleration and jerk). All these computations are natively
thread-safe and compatible with Automatic Differentiation (AD) with respect to time through [DifferentiationInterface.jl](https://juliadiff.org/DifferentiationInterface.jl/stable/).

In particular, the following methods are available to compute translation data:

```
ephem_vector3(eph, from, to, time)
ephem_vector6(eph, from, to, time)
ephem_vector9(eph, from, to, time)
ephem_vector12(eph, from, to, time)
```

They all share the same interface, requiring an `EphemerisProvider` object as the first input. `from` and `to` are integer numbers representing the ID of the center and target points that we desired. The `time` argument is expressed in TDB seconds since J2000.0.

!!! note
    Differently, from traditional ephemerides readers, `Ephemerides.jl` is only meant to read the data stored in the binary kernels and it does not perform any concatenation of state vectors. This means that if data from point 399 is expressed with respect to point 3, we will only be able to compute the relative position of 339 with respect to 3 or viceversa, but not of 399 with respect to another point. The reason behind this is that `Ephemerides.jl` is meant to be used in combination with [FrameTransformations.jl](https://juliaspacemissiondesign.github.io/FrameTransformations.jl/stable/), which already enables tranformations between different user-defined point and axes.

An example to compute the position of the Moon (399) with respect to the Earth-Moon Barycenter (3) at J2000 (time = 0), is the following:

```julia
using Ephemerides

# Load an ephemeris kernel containing the requested data
kernel = "path_to_de440"
eph = EphemerisProvider(kernel)

# Compute the position
pos = ephem_vector3(eph, 3, 399, 0)
```

Instead, if one desires the whole state vector, up to the jerk components, the functions become:

```julia
# Compute position and velocity
pv = ephem_vector6(eph, 3, 399, 0)

# Compute position, velocity and acceleration
pva = ephem_vector9(eph, 3, 399, 0)

# Compute position, velocity, acceleration and jerk
pvaj = ephem_vector12(eph, 3, 399, 0)
```

In all these examples, the returned data is always in the form of a `StaticArray` in order to minimise memory allocations.

## Preparing routes for repeated calls

When the same point pair is evaluated repeatedly—such as inside a third-body
perturbation model—resolve the route once:

```julia
route = prepare_ephemeris(
    eph, 3, 399; timespan=(0.0, 365.25 * 86400.0)
)

for time in integration_times
    position = ephem_vector3(route, time)
    # use position ...
end
```

The `timespan` keyword is optional. Supplying it lets `prepare_ephemeris` select one
concrete segment when that segment covers the complete integration interval. Calls
through the prepared route then bypass provider lookup and segment dispatch while
remaining thread-safe, allocation-free, and compatible with automatic differentiation.
If multiple priority-ordered kernels are required to cover the interval, the prepared
route retains them and uses a small per-thread link cache.

For several repeatedly queried bodies, keep the prepared routes in a tuple and evaluate
the tuple as a static batch. This preserves concrete dispatch for every segment and
returns a tuple of static vectors without allocating:

```julia
third_body_routes = (
    prepare_ephemeris(eph, 0, 10; timespan=propagation_span),
    prepare_ephemeris(eph, 0, 5; timespan=propagation_span),
    prepare_ephemeris(eph, 0, 6; timespan=propagation_span),
)

sun_state, jupiter_state, saturn_state = ephem_vector6(third_body_routes, time)
```

Orientations have the equivalent interface:

```julia
route = prepare_orientation(eph, 1, 31008; timespan=(0.0, 86400.0))
angles = ephem_rotation3(route, 0.0)
```

## Automatic differentiation backends

Ephemerides does not select or depend on one AD backend. Use
[DifferentiationInterface.jl](https://juliadiff.org/DifferentiationInterface.jl/stable/)
to select the backend at the call site:

```julia
using DifferentiationInterface
using ForwardDiff: ForwardDiff
using Zygote: Zygote

point_route = prepare_ephemeris(eph, 3, 399; timespan=(0.0, 86400.0))
position(time) = ephem_vector3(point_route, time)

velocity_forward = derivative(position, AutoForwardDiff(), 0.0)
velocity_reverse = derivative(position, AutoZygote(), 0.0)
```

Forward-mode backends propagate their number types through the interpolation routines.
When ChainRulesCore is loaded, Ephemerides also loads an extension with analytic `frule`
and `rrule` definitions for `ephem_vector3`, `ephem_vector6`, `ephem_vector9`, and their
orientation equivalents. These rules apply to provider and prepared-route calls and avoid
tracing mutation of internal kernel caches.

The test suite validates `AutoForwardDiff`, `AutoZygote`, and central
`AutoFiniteDiff`. Other DifferentiationInterface backends may work, but backends that do
not consume ChainRules rules must support the package's cache mutation and data-dependent
segment selection. For higher time derivatives, use `ephem_vector9` or
`ephem_vector12` directly, or a forward-mode backend; nested reverse-mode differentiation
is not part of the supported contract.

!!! warning
    SPK segments of types 1, 5, 15, 17, and 21 do not natively support acceleration and jerk computations. However, these values can be computed by Automatic Differentiation (AD) of the position and/or velocity components.


## Computing orientation angles

Similarly to position components, `Ephemerides.jl` also allows the computation of orientation angles and their derivatives up to order 3. All these computations are natively thread-safe and compatible with Automatic Differentiation (AD) with respect to time through DifferentiationInterface.

In particular, the following methods are available to compute orientation data:

```
ephem_rotation3(eph, from, to, time)
ephem_rotation6(eph, from, to, time)
ephem_rotation9(eph, from, to, time)
ephem_rotation12(eph, from, to, time)
```

Again, they all share the same interface, requiring an [`EphemerisProvider`](@ref) object as the first input. `from` and `to` are integer numbers representing the ID of the reference and target axes that we desired. The `time` argument is expressed in TDB seconds since J2000.0.

An example to compute the Euler angles of the PA440 axes (31008) with respect to the ICRF (1) at J2000 (time = 0), is the following:

```julia
using Ephemerides

# Load an ephemeris kernel containing the requested data
kernel = "path_to_pa440"
eph = EphemerisProvider(kernel)

# Compute the position
ang = ephem_rotation3(eph, 1, 31008, 0)
```

Instead, if one desires the whole vector, up to the 3rd order derivative, the functions become:

```julia
# Compute angles and derivatives
pv = ephem_rotation6(eph, 1, 31008, 0)
pva = ephem_rotation9(eph, 1, 31008, 0)
pvaj = ephem_rotation12(eph, 1, 31008, 0)
```

The returned orientation data is always in the form of a `StaticArray` in order to minimise memory allocations.


!!! note
    Differently from the translational data contained in SPK kernels, the orientation angles can only be computed in one direction, i.e., if the orientation of the Moon's Principal Axes (PA) is defined with respect to the ICRF, it is not possible to compute the rotation from the PA to the ICRF with this routine.
