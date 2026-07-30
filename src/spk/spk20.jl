
"""
    SPKSegmentHeader20(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 20
"""
function SPKSegmentHeader20(daf::DAF, desc::DAFSegmentDescriptor)

    i0 = 8*(final_address(desc)-7)

    # Length and time scales in km and seconds
    dscale = get_float(array(daf), i0, endian(daf))
    tscale = get_float(array(daf), i0 + 8, endian(daf))

    # Everything is now transformed in seconds past J2000 and kms
    DJ2000, D2S = 2451545, 86400

    # Integer and fractional part of the initial juliad date
    initjd = get_float(array(daf), i0 + 16, endian(daf))
    initfr = get_float(array(daf), i0 + 24, endian(daf))

    # Start epoch and interval length (in seconds)
    tstart = ((initjd - DJ2000) + initfr)*D2S
    tlen = get_float(array(daf), i0 + 32, endian(daf))*D2S

    # Number of elements in each record
    rsize = Int(get_float(array(daf), i0 + 40, endian(daf)))

    # Byte size of each logical record
    recsize = 8*rsize

    # Polynomial degree
    order = div(rsize - 3, 3) - 1

    # Polynomial group size (number of coefficients required for the interpolation)
    N = order + 1

    # Number of records
    n = Int(get_float(array(daf), i0 + 48, endian(daf)))

    # Initial segment address
    iaa = initial_address(desc)

    SPKSegmentHeader20(dscale, tscale, tstart, tlen, recsize, order, N, n, iaa)
end

"""
    SPKSegmentCache20(head::SPKSegmentHeader20)

Initialise the cache for an SPK segment of type 20.
"""
function SPKSegmentCache20(head::SPKSegmentHeader20)
    SPKSegmentCache20(
        -1,
        MVector(0.0, 0.0, 0.0),
        zeros(3, max(3, head.N))
        )
end

"""
    SPKSegmentType20(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 20.
"""
function SPKSegmentType20(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader20(daf, desc)
    caches = build_thread_cache(() -> SPKSegmentCache20(header))

    SPKSegmentType20(header, caches)

end

@inline spk_field(::SPKSegmentType20) = SPK_SEGMENTLIST_MAPPING[20]

function spk_vector3(daf::DAF, seg::SPKSegmentType20, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve length and time scales
    length_scale = head.dscale
    time_scale = head.tscale

    # Find the logical record containing the Chebyshev coefficients at `time`
    index = find_logical_record(head, time)
    get_coefficients!(daf, head, data, index)

    # Normalise the time argument between [-1, 1]
    t = normalise_time(head, time, index)

    # Compute the position
    x, y, z = chebyshev_integral(data.A, t, head.N, time_scale, head.tlen, data.p)

    return SVector{3}(length_scale*x, length_scale*y, length_scale*z)

end


function spk_vector6(daf::DAF, seg::SPKSegmentType20, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve length and time scales
    length_scale = head.dscale
    time_scale = head.tscale

    velocity_scale = length_scale/time_scale

    # Find the logical record containing the Chebyshev coefficients at `time`
    index = find_logical_record(head, time)
    get_coefficients!(daf, head, data, index)

    # Normalise the time argument between [-1, 1]
    t = normalise_time(head, time, index)

    # Compute the velocity
    vx, vy, vz = chebyshev(data.A, t, 0, head.N)

    # Compute the position
    x, y, z = chebyshev_integral(data.A, t, head.N, time_scale, head.tlen, data.p)

    return SVector{6}(length_scale*x, length_scale*y, length_scale*z, velocity_scale*vx, velocity_scale*vy, velocity_scale*vz)

end

function spk_vector9(daf::DAF, seg::SPKSegmentType20, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve length and time scales
    length_scale = head.dscale
    time_scale = head.tscale

    velocity_scale = length_scale/time_scale

    # Find the logical record containing the Chebyshev coefficients at `time`
    index = find_logical_record(head, time)
    get_coefficients!(daf, head, data, index)

    # Normalise the time argument between [-1, 1]
    t = normalise_time(head, time, index)

    # Compute the velocity and acceleration
    vx, vy, vz, ax, ay, az =
        chebyshev_derivative1(data.A, t, 0, head.N, 2/head.tlen)

    # Compute the position
    x, y, z = chebyshev_integral(data.A, t, head.N, time_scale, head.tlen, data.p)

    return SVector{9}(
        length_scale*x, length_scale*y, length_scale*z,
        velocity_scale*vx, velocity_scale*vy, velocity_scale*vz,
        velocity_scale*ax, velocity_scale*ay, velocity_scale*az
    )

end

function spk_vector12(daf::DAF, seg::SPKSegmentType20, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve length and time scales
    length_scale = head.dscale
    time_scale = head.tscale

    velocity_scale = length_scale/time_scale

    # Find the logical record containing the Chebyshev coefficients at `time`
    index = find_logical_record(head, time)
    get_coefficients!(daf, head, data, index)

    # Normalise the time argument between [-1, 1]
    t = normalise_time(head, time, index)

    # Compute the velocity and acceleration and jerk
    vx, vy, vz, ax, ay, az, jx, jy, jz = chebyshev_derivative2(
        data.A, t, 0, head.N, 2/head.tlen
    )

    # Compute the position
    x, y, z = chebyshev_integral(data.A, t, head.N, time_scale, head.tlen, data.p)

    return SVector{12}(
        length_scale*x, length_scale*y, length_scale*z,
        velocity_scale*vx, velocity_scale*vy, velocity_scale*vz,
        velocity_scale*ax, velocity_scale*ay, velocity_scale*az,
        velocity_scale*jx, velocity_scale*jy, velocity_scale*jz
    )

end


"""
    find_logical_record(head::SPKSegmentHeader20, time::Number)
"""
function find_logical_record(head::SPKSegmentHeader20, time::Number)

    # The index is returned in 0-based notation (i.e., the first record has index 0)
    index = floor(Int, (time - head.tstart)/head.tlen)

    if index == head.n
        # This happens only when the time equals the final segment time
        index -= 1
    end

    return index
end

"""
    get_coefficients!(daf::DAF, head::SPKSegmentHeader20, cache::SPKSegmentCache20, index::Int)
"""
function get_coefficients!(
    daf::DAF, head::SPKSegmentHeader20, cache::SPKSegmentCache20, index::Int
    )

    # Check whether the coefficients for this record are already loaded
    index == cache.id && return nothing
    cache.id = index

    # Address of desired logical record
    k = 8*(head.iaa-1) + head.recsize*index

    # For type 20 the position constants follow each component's coefficients.
    bytes = array(daf)
    little_endian = endian(daf)
    stride = head.N + 1
    GC.@preserve bytes begin
        source = Ptr{Float64}(pointer(bytes, k + 1))
        @inbounds for component in 1:3
            offset = (component - 1)*stride
            for coefficient in 1:head.N
                cache.A[component, coefficient] = get_num(
                    unsafe_load(source, offset + coefficient), little_endian
                )
            end
            cache.p[component] =
                get_num(unsafe_load(source, offset + stride), little_endian)
        end
    end

    nothing

end

"""
    normalise_time(head::SPKSegmentHeader20, time::Number, index::Int)
"""
function normalise_time(head::SPKSegmentHeader20, time::Number, index::Int)
    tbeg = head.tstart + head.tlen*index
    hlen = head.tlen/2
    return (time - tbeg)/hlen - 1
end
