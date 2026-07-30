
"""
    chebyshev(coefficients, t::Number, idx::Int, N::Int)

Evaluate a sum of Chebyshev polynomials of the first kind at `t` using a
recursive algorithm. It simultaneously evaluates the three state components. `idx` is the
index of the starting row (in 0-based notation) in the matrix of coefficients `coefficients`.

!!! note
    `t` is the normalized abscissa in `[-1, 1]`.

### See Also
See also [`chebyshev_derivative1`](@ref), [`chebyshev_derivative2`](@ref) and [`chebyshev_derivative3`](@ref)

"""
function chebyshev(coefficients, t::Number, idx::Int, N::Int)
    ix, iy, iz = idx + 1, idx + 2, idx + 3

    # Only the two preceding polynomial values are needed by the recurrence.
    # Keeping them in scalar rolling state avoids work-buffer memory traffic.
    poly_prev2 = one(t)
    poly_prev1 = t
    poly = 2t*t - one(t)

    @inbounds begin
        x = coefficients[ix, 1] + t*coefficients[ix, 2] + poly*coefficients[ix, 3]
        y = coefficients[iy, 1] + t*coefficients[iy, 2] + poly*coefficients[iy, 3]
        z = coefficients[iz, 1] + t*coefficients[iz, 2] + poly*coefficients[iz, 3]

        for j = 4:N
            poly_prev2, poly_prev1 = poly_prev1, poly
            poly = 2t*poly_prev1 - poly_prev2
            x += poly*coefficients[ix, j]
            y += poly*coefficients[iy, j]
            z += poly*coefficients[iz, j]
        end
    end

    return x, y, z
end


"""
    chebyshev_derivative1(coefficients, t::Number, idx::Int, N::Int, time_scale)

Evaluate a sum of Chebyshev polynomials of the first kind and its derivative at `t`
using a recursive algorithm. It simultaneously evaluates the three state components. `idx`
is the index of the starting row (in 0-based notation) in the matrix of coefficients `coefficients`.

!!! note
    `t` is the normalized abscissa in `[-1, 1]`.

### See Also
See also [`chebyshev`](@ref), [`chebyshev_derivative2`](@ref) and [`chebyshev_derivative3`](@ref)

"""
function chebyshev_derivative1(
    coefficients, t::Number, idx::Int, N::Int, time_scale::Number
)
    ix, iy, iz = idx + 1, idx + 2, idx + 3

    poly_prev2, poly_prev1, poly = one(t), t, 2t*t - one(t)
    poly_d1_prev2, poly_d1_prev1, poly_d1 = zero(t), one(t), 4t

    @inbounds begin
        x = coefficients[ix, 1] + t*coefficients[ix, 2] + poly*coefficients[ix, 3]
        y = coefficients[iy, 1] + t*coefficients[iy, 2] + poly*coefficients[iy, 3]
        z = coefficients[iz, 1] + t*coefficients[iz, 2] + poly*coefficients[iz, 3]
        vx = coefficients[ix, 2] + poly_d1*coefficients[ix, 3]
        vy = coefficients[iy, 2] + poly_d1*coefficients[iy, 3]
        vz = coefficients[iz, 2] + poly_d1*coefficients[iz, 3]

        for j = 4:N
            poly_prev2, poly_prev1 = poly_prev1, poly
            poly_d1_prev2, poly_d1_prev1 = poly_d1_prev1, poly_d1
            poly = 2t*poly_prev1 - poly_prev2
            poly_d1 = 2t*poly_d1_prev1 + 2poly_prev1 - poly_d1_prev2
            x += poly*coefficients[ix, j]
            y += poly*coefficients[iy, j]
            z += poly*coefficients[iz, j]
            vx += poly_d1*coefficients[ix, j]
            vy += poly_d1*coefficients[iy, j]
            vz += poly_d1*coefficients[iz, j]
        end
    end

    return x, y, z, time_scale*vx, time_scale*vy, time_scale*vz
end


"""
    chebyshev_derivative2(coefficients, t::Number, idx::Int, N::Int, time_scale)

Evaluate a sum of Chebyshev polynomials of the first kind and its two derivatives at `t`
using a recursive algorithm. It simultaneously evaluates the three state components. `idx`
is the index of the starting row (in 0-based notation) in the matrix of coefficients `coefficients`.

!!! note
    `t` is the normalized abscissa in `[-1, 1]`.

### See Also
See also [`chebyshev`](@ref), [`chebyshev_derivative1`](@ref) and [`chebyshev_derivative3`](@ref)

"""
function chebyshev_derivative2(
    coefficients, t::Number, idx::Int, N::Int, time_scale
)
    ix, iy, iz = idx + 1, idx + 2, idx + 3

    poly_prev2, poly_prev1, poly = one(t), t, 2t*t - one(t)
    poly_d1_prev2, poly_d1_prev1, poly_d1 = zero(t), one(t), 4t
    poly_d2_prev2, poly_d2_prev1, poly_d2 = zero(t), zero(t), 4one(t)

    @inbounds begin
        x = coefficients[ix, 1] + t*coefficients[ix, 2] + poly*coefficients[ix, 3]
        y = coefficients[iy, 1] + t*coefficients[iy, 2] + poly*coefficients[iy, 3]
        z = coefficients[iz, 1] + t*coefficients[iz, 2] + poly*coefficients[iz, 3]
        vx = coefficients[ix, 2] + poly_d1*coefficients[ix, 3]
        vy = coefficients[iy, 2] + poly_d1*coefficients[iy, 3]
        vz = coefficients[iz, 2] + poly_d1*coefficients[iz, 3]
        ax = poly_d2*coefficients[ix, 3]
        ay = poly_d2*coefficients[iy, 3]
        az = poly_d2*coefficients[iz, 3]

        for j = 4:N
            poly_prev2, poly_prev1 = poly_prev1, poly
            poly_d1_prev2, poly_d1_prev1 = poly_d1_prev1, poly_d1
            poly_d2_prev2, poly_d2_prev1 = poly_d2_prev1, poly_d2
            poly = 2t*poly_prev1 - poly_prev2
            poly_d1 = 2t*poly_d1_prev1 + 2poly_prev1 - poly_d1_prev2
            poly_d2 = 2t*poly_d2_prev1 + 4poly_d1_prev1 - poly_d2_prev2
            x += poly*coefficients[ix, j]
            y += poly*coefficients[iy, j]
            z += poly*coefficients[iz, j]
            vx += poly_d1*coefficients[ix, j]
            vy += poly_d1*coefficients[iy, j]
            vz += poly_d1*coefficients[iz, j]
            ax += poly_d2*coefficients[ix, j]
            ay += poly_d2*coefficients[iy, j]
            az += poly_d2*coefficients[iz, j]
        end
    end

    time_scale_squared = time_scale*time_scale
    return (
        x, y, z,
        time_scale*vx, time_scale*vy, time_scale*vz,
        time_scale_squared*ax, time_scale_squared*ay, time_scale_squared*az,
    )
end

"""
    chebyshev_derivative3(coefficients, t::Number, idx::Int, N::Int, time_scale)

Evaluate a sum of Chebyshev polynomials of the first kind and its three derivatives at `t`
using a recursive algorithm. It simultaneously evaluates the three state components. `idx`
is the index of the starting row (in 0-based notation) in the matrix of coefficients `coefficients`.

!!! note
    `t` is the normalized abscissa in `[-1, 1]`.

### See Also
See also [`chebyshev`](@ref), [`chebyshev_derivative1`](@ref) and [`chebyshev_derivative2`](@ref)

"""
function chebyshev_derivative3(
    coefficients, t::Number, idx::Int, N::Int, time_scale
)
    ix, iy, iz = idx + 1, idx + 2, idx + 3

    poly_prev2, poly_prev1, poly = one(t), t, 2t*t - one(t)
    poly_d1_prev2, poly_d1_prev1, poly_d1 = zero(t), one(t), 4t
    poly_d2_prev2, poly_d2_prev1, poly_d2 = zero(t), zero(t), 4one(t)
    poly_d3_prev2, poly_d3_prev1, poly_d3 = zero(t), zero(t), zero(t)

    @inbounds begin
        x = coefficients[ix, 1] + t*coefficients[ix, 2] + poly*coefficients[ix, 3]
        y = coefficients[iy, 1] + t*coefficients[iy, 2] + poly*coefficients[iy, 3]
        z = coefficients[iz, 1] + t*coefficients[iz, 2] + poly*coefficients[iz, 3]
        vx = coefficients[ix, 2] + poly_d1*coefficients[ix, 3]
        vy = coefficients[iy, 2] + poly_d1*coefficients[iy, 3]
        vz = coefficients[iz, 2] + poly_d1*coefficients[iz, 3]
        ax = poly_d2*coefficients[ix, 3]
        ay = poly_d2*coefficients[iy, 3]
        az = poly_d2*coefficients[iz, 3]
        jx = zero(t)
        jy = zero(t)
        jz = zero(t)

        for j = 4:N
            poly_prev2, poly_prev1 = poly_prev1, poly
            poly_d1_prev2, poly_d1_prev1 = poly_d1_prev1, poly_d1
            poly_d2_prev2, poly_d2_prev1 = poly_d2_prev1, poly_d2
            poly_d3_prev2, poly_d3_prev1 = poly_d3_prev1, poly_d3
            poly = 2t*poly_prev1 - poly_prev2
            poly_d1 = 2t*poly_d1_prev1 + 2poly_prev1 - poly_d1_prev2
            poly_d2 = 2t*poly_d2_prev1 + 4poly_d1_prev1 - poly_d2_prev2
            poly_d3 = 2t*poly_d3_prev1 + 6poly_d2_prev1 - poly_d3_prev2
            x += poly*coefficients[ix, j]
            y += poly*coefficients[iy, j]
            z += poly*coefficients[iz, j]
            vx += poly_d1*coefficients[ix, j]
            vy += poly_d1*coefficients[iy, j]
            vz += poly_d1*coefficients[iz, j]
            ax += poly_d2*coefficients[ix, j]
            ay += poly_d2*coefficients[iy, j]
            az += poly_d2*coefficients[iz, j]
            jx += poly_d3*coefficients[ix, j]
            jy += poly_d3*coefficients[iy, j]
            jz += poly_d3*coefficients[iz, j]
        end
    end

    time_scale_squared = time_scale*time_scale
    time_scale_cubed = time_scale_squared*time_scale
    return (
        x, y, z,
        time_scale*vx, time_scale*vy, time_scale*vz,
        time_scale_squared*ax, time_scale_squared*ay, time_scale_squared*az,
        time_scale_cubed*jx, time_scale_cubed*jy, time_scale_cubed*jz,
    )
end

# Compatibility methods for callers using the former work-buffer API. Scalar rolling
# recurrence no longer reads or mutates these buffers.
@inline chebyshev(::InterpCache, coefficients, t::Number, idx::Int, N::Int, ibuff=1) =
    chebyshev(coefficients, t, idx, N)
@inline chebyshev_derivative1(
    ::InterpCache, coefficients, t::Number, idx::Int, N::Int, time_scale::Number, ibuff=1
) = chebyshev_derivative1(coefficients, t, idx, N, time_scale)
@inline chebyshev_derivative2(
    ::InterpCache, coefficients, t::Number, idx::Int, N::Int, time_scale, ibuff=1
) = chebyshev_derivative2(coefficients, t, idx, N, time_scale)
@inline chebyshev_derivative3(
    ::InterpCache, coefficients, t::Number, idx::Int, N::Int, time_scale, ibuff=1
) = chebyshev_derivative3(coefficients, t, idx, N, time_scale)


"""
    chebyshev_integral(coefficients, t::Number, N::Int, time_scale, tlen, midpoint_position)

Evaluate the integral of a sum of Chebyshev polynomials of the first kind using a recursive
algorithm. It simultaneously evaluates the three state components. `tlen` is the size of the
record interval, `time_scale` is the timescale factor, and `midpoint_position` contains the position coefficients
at the midpoint (i.e., when the integral is evaluated at `t = 0`).

!!! note
    `t` is the normalized abscissa in `[-1, 1]`.
"""
function chebyshev_integral(coefficients, t::Number, N::Int, time_scale, tlen, midpoint_position)
    poly_prev2 = t
    poly_prev1 = 2t*t - one(t)
    integral_t2 = t*t/2

    @inbounds begin
        x = t*coefficients[1, 1] + integral_t2*coefficients[1, 2]
        y = t*coefficients[2, 1] + integral_t2*coefficients[2, 2]
        z = t*coefficients[3, 1] + integral_t2*coefficients[3, 2]

        for j in 4:(N + 1)
            poly = 2t*poly_prev1 - poly_prev2
            ipoly = (poly/(j - 1) - poly_prev2/(j - 3))/2

            # Starting from the integral of T5, the above formula misses a constant of
            # integration on all the odd polynomials (T5, T7, T9, T11, etc...). The sign of
            # this constant also alternatively changes between 1 and -1.
            # Here we are including that constant to ensure that the integral expression is
            # null when evaluated at t = 0.

            if isodd(j)
                d = (1/(j - 1) + 1/(j - 3))/2
                ipoly += ((j - 5) % 4 == 0) ? -d : d
            end

            x += ipoly*coefficients[1, j - 1]
            y += ipoly*coefficients[2, j - 1]
            z += ipoly*coefficients[3, j - 1]
            poly_prev2, poly_prev1 = poly_prev1, poly
        end

        # Multiplies by tlen/2 because the t coordinate we are using is expressed
        # as function of the original epoch as t = 2(x - mid)/tlen, so that
        # dt/dx = 2/tlen
        x = midpoint_position[1] + x*tlen/time_scale/2
        y = midpoint_position[2] + y*tlen/time_scale/2
        z = midpoint_position[3] + z*tlen/time_scale/2

    end

    return x, y, z

end

# Compatibility with the former work-buffer methods.
@inline chebyshev_integral(
    ::InterpCache, coefficients, t::Number, N::Int, time_scale, tlen, midpoint_position
) = chebyshev_integral(coefficients, t, N, time_scale, tlen, midpoint_position)
@inline chebyshev_integral(
    ::InterpCache, coefficients, t::Number, N::Int, time_scale, tlen, midpoint_position, ibuff::Int
) = chebyshev_integral(coefficients, t, N, time_scale, tlen, midpoint_position)
