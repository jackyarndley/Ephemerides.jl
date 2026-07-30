

"""
    TwoBodyUniversalCache

A container to store precomputed quantities required for the two-body propagation based on
universal variables. Only the quantities that depend on the initial state are computed.
"""
mutable struct TwoBodyUniversalCache

    pos::Vector{Float64}
    vel::Vector{Float64}

    F::Float64
    qr::Float64

    bq::Float64
    bv::Float64
    br0::Float64

    bound::Float64

    gm::Float64
    p::Vector{Float64}

end

function TwoBodyUniversalCache(gm::Number, deg::Int=11)

    # Pre-compute the pairs required for the Stumpff functions series evaluation.
    # The truncation degree is set to 11, similarly to SPICE.

    p = zeros(2(deg-1))
    @inbounds for j in eachindex(p)
        p[j] = 1/(j*(j+1))
    end

    TwoBodyUniversalCache(zeros(3), zeros(3), 0, 0, 0, 0, 0, 0, gm, p)

end

function TwoBodyUniversalCache(gm::Number, pos, vel, deg::Int=11)
    c = TwoBodyUniversalCache(gm, deg)
    c.pos[:] = pos
    c.vel[:] = vel
    update_cache!(c)
    return c
end


"""
    update_cache!(c::TwoBodyUniversalCache)

Update the precomputed values in the cache using the position and velocity in `c`.
"""
function update_cache!(c::TwoBodyUniversalCache)

    # Compute the initial position vector magnitude and radial velocity
    r0 = vnorm(c.pos)
    vr = vdot(c.pos, c.vel)

    # Compute the angular momentum
    h0 = vcross(c.pos, c.vel)
    h2 = vdot(h0, h0)

    # Compute the eccentricity from the eccentricity vector
    @inbounds eccentricity_vector = vcross(c.vel, h0)/c.gm - SVector{3}(c.pos)/r0
    ecc = vnorm(eccentricity_vector)

    # Compute Q = a * (1-e)
    q = h2/(c.gm*(1+ecc))
    b = sqrt(q/c.gm)

    # Precompute the parameters that do not depend on time_offset for Kepler's Equation
    c.F = 1 - ecc
    c.qr = q/r0

    c.bq  = b*q
    c.bv  = b*b*vr
    c.br0 = b*r0

    # Compute the upper and lower boundaries on the value of x that guarantee there
    # will not be any overflow problems in the Universal Variables Kepler's equation.

    maxc = max(1, abs(c.br0), abs(c.bv), abs(c.bq), abs(q/c.bq))
    dpmax = 1.7976931348623157e308

    if c.F < 0
        f = log(dpmax / 2) - log(maxc)
        rF = sqrt(-c.F)
        c.bound = min(f/rF, (f + 3/2*log(-c.F))/rF)
    else
        c.bound = exp((log(3/2) + log(dpmax) - log(maxc))/3)
    end

    nothing
end


"""
    propagate_twobody(cache::TwoBodyUniversalCache, time_offset::Number)

Propagate the state vector in `cache` of `time_offset` using the universal variables formulation for
Kepler's Equation and the Lagrange coefficients f and g.

!!! note
    This routine is valid for any type of orbit and uses a bisection method to find the
    root of the universal variables Kepler's equation. The algorithm has been directly
    taken from the SPICE toolkit `prob2b.f`.

### References
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)
"""
function propagate_twobody(cache::TwoBodyUniversalCache, time_offset)

    # There is no need to propagate
    if time_offset == 0
        T = promote_type(typeof(time_offset), eltype(cache.pos))
        return SVector{3, T}(cache.pos), SVector{3, T}(cache.vel)
    end

    # Retrieve the usefull parameters
    F = cache.F
    qr = cache.qr

    bq = cache.bq
    bv = cache.bv
    br0 = cache.br0

    x = min(cache.bound, max(-cache.bound, time_offset/bq))
    fx_squared = F*x*x

    # Compute the Stumpff functions
    c0, c1, c2, c3 = stumpff(fx_squared, cache.p)

    # Compute Kepler's equation
    kfun = x*(br0*c1 + x*(bv*c2 + x*bq*c3))

    # Compute the bracketing interval lower and upper boundaries. Since due to some
    # rounding errors it may happend that the root has not been bracketed, in those cases
    # we keep doubling the endpoint farthest from zero until the root is bracketed.

    if time_offset > 0
        lb, ub = zero(x), x

        while kfun < time_offset
            lb, ub = ub, 2ub

            # Make sure that x remains within the bounds we have previously defined.
            xold, x = x, min(cache.bound, max(-cache.bound, ub))

            if x == xold
                throw(
                    ErrorException(
                        "The input time_offset is beyond the range of time for which we can"*
                        " reliably propagate the states."
                    )
                )
            end

            fx_squared = F*x*x
            c0, c1, c2, c3 = stumpff(fx_squared, cache.p)
            kfun = x*(br0*c1 + x*(bv*c2 + x*bq*c3))
        end

    else
        lb, ub = x, zero(x)

        while kfun > time_offset
            lb, ub = 2lb, lb

            # Make sure that x remains within the bounds we have previously defined
            xold, x = x, min(cache.bound, max(-cache.bound, lb))

            if x == xold
                throw(
                    ErrorException(
                        "The input time_offset is beyond the range of time for which we can"*
                        " reliably propagate the states."
                    )
                )
            end

            fx_squared = F*x*x
            c0, c1, c2, c3 = stumpff(fx_squared, cache.p)
            kfun = x*(br0*c1 + x*(bv*c2 + x*bq*c3))
        end
    end

    # Apply the bisection method to find the solution of the equation
    x = min(ub, max(lb, (ub+lb)/2))

    fx_squared = F*x*x
    c0, c1, c2, c3 = stumpff(fx_squared, cache.p)

    tol, maxiter = 1e-14, 100
    iter, err  = 0, tol + 1
    while iter < maxiter && err > tol && x > lb && x < ub

        # Compute Kepler's equation
        fun = x*(br0*c1 + x*(bv*c2 + x*bq*c3)) - time_offset

        if fun > 0
            ub = x
        elseif fun < 0
            lb = x
        else
            ub, lb = x, x
        end

        xold, x = x, min(ub, max(lb, (ub+lb)/2))

        fx_squared = F*x*x
        c0, c1, c2, c3 = stumpff(fx_squared, cache.p)

        iter += 1
        err = abs(x - xold)

    end

    # Now that x has been found, we can compute the Lagrange functions f, g, df, dg
    # and then propagate our state.

    x2 = x*x
    x3 = x*x2

    br = br0*c0 + x*(bv*c1 + x*(bq*c2))

    f = 1 - qr*x2*c2
    g = time_offset - bq*x3*c3
    df = - qr/br*x*c1
    dg = 1 - bq/br*x2*c2

    # Compute the position and velocity at t + time_offset
    pos = SVector{3}(cache.pos)
    vel = SVector{3}(cache.vel)

    pnew = f*pos + g*vel
    vnew = df*pos + dg*vel

    return pnew, vnew

end

"""
    stumpff(x::Number, p::AbstractVector)

Compute Stumpff's functions from C0 up to C3 at `x`.

!!! note
    This routine uses the trigonometrical expressions of the functions when the absolute
    value of `x` is greater or equal to 1. If that is not the case, the C2 and C3 functions
    are computed from a truncated expression of the Maclaurin series at order 11, which
    guarantees a higher precision and avoid overflow errors when `x` is null.

### References
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)
"""
function stumpff(x, p)

    if x >= 1
        sqrt_x = sqrt(x)

        c0 = cos(sqrt_x)
        c1 = sin(sqrt_x)/sqrt_x
        c2 = (1 - c0)/x
        c3 = (1 - c1)/x

        return c0, c1, c2, c3

    elseif x <= -1
        sqrt_x = sqrt(-x)

        c0 = cosh(sqrt_x)
        c1 = sinh(sqrt_x)/sqrt_x
        c2 = (1 - c0)/x
        c3 = (1 - c1)/x

        return c0, c1, c2, c3

    end

    # Manually compute the Maclaurin series
    n = length(p)

    # Compute the series expansion for C3:
    c3 = one(x)
    @inbounds for j = n:-2:4
        c3 = 1 - x*p[j]*c3
    end
    c3 *= p[2]

    # Compute the series expansion for C2:
    c2 = one(x)
    @inbounds for j = n-1:-2:3
        c2 = 1 - x*p[j]*c2
    end
    c2 *= p[1]

    # Compute C0 and C1 from the recursion formula:
    c1 = 1 - x*c3
    c0 = 1 - x*c2

    return c0, c1, c2, c3

end
