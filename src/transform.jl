export ephem_vector3, ephem_vector6, ephem_vector9, ephem_vector12,
       ephem_rotation3, ephem_rotation6, ephem_rotation9, ephem_rotation12

@inline _covers(link::SPKLink, time::Number) =
    initial_time(link) <= time <= final_time(link)

@inline function _cached_route(
    cache::RouteCache, from::Int, to::Int, time::Number
)
    slot = route_cache_slot(from, to)
    link = @inbounds cache.links[slot]

    if @inbounds(cache.from[slot] == from && cache.to[slot] == to) &&
            link !== nothing && _covers(link, time)
        return link, @inbounds(cache.priority[slot])
    end

    return nothing, 0
end

@inline function _cache_link!(
    cache::RouteCache, from::Int, to::Int, link::SPKLink, priority::Int
)
    slot = route_cache_slot(from, to)
    @inbounds begin
        cache.from[slot] = from
        cache.to[slot] = to
        cache.priority[slot] = priority
        cache.links[slot] = link
    end
    return link
end

function _find_link(
    table::FlatLinkTable,
    cache::RouteCache,
    from::Int,
    to::Int,
    time::Number,
    object::AbstractString,
)
    cached, priority = _cached_route(cache, from, to, time)
    cached !== nothing && priority == 1 && return cached

    candidates = get(table, (from, to), nothing)
    candidates === nothing && throw(jEph.EphemerisError(
        "ephemeris data for $object with NAIFId $to with respect to $object $from is unavailable."
    ))

    return _find_candidate(candidates, cache, from, to, time, object)
end

function _find_candidate(
    candidates::Vector{SPKLink},
    cache::RouteCache,
    from::Int,
    to::Int,
    time::Number,
    object::AbstractString,
)
    cached, _ = _cached_route(cache, from, to, time)

    for (priority, link) in pairs(candidates)
        _covers(link, time) || continue
        link === cached && return cached
        return _cache_link!(cache, from, to, link, priority)
    end

    throw(jEph.EphemerisError(
        "ephemeris data for $object with NAIFId $to with respect to $object $from " *
        "is not available at $time seconds since J2000."
    ))
end

@inline function _find_spk_link(
    eph::EphemerisProvider, from::Int, to::Int, time::Number
)
    return _find_link(
        spk_routes(eph), thread_cache(eph.spkcache), from, to, time, "point"
    )
end

@inline function _find_pck_link(
    eph::EphemerisProvider, from::Int, to::Int, time::Number
)
    return _find_link(
        pck_routes(eph), thread_cache(eph.pckcache), from, to, time, "axes"
    )
end

@inline function _segment_from_link(daf::DAF, link::SPKLink)
    lid = list_id(link)
    eid = element_id(link)
    segments = segment_list(daf)

    lid == 1 && return @inbounds segments.spk2[eid]
    lid == 2 && return @inbounds segments.spk9[eid]
    lid == 3 && return @inbounds segments.spk1[eid]
    lid == 4 && return @inbounds segments.spk14[eid]
    lid == 5 && return @inbounds segments.spk15[eid]
    lid == 6 && return @inbounds segments.spk8[eid]
    lid == 7 && return @inbounds segments.spk19[eid]
    lid == 8 && return @inbounds segments.spk20[eid]
    lid == 9 && return @inbounds segments.spk5[eid]
    lid == 10 && return @inbounds segments.spk17[eid]

    throw(jEph.EphemerisError("invalid internal segment-list id $lid."))
end

@inline _evaluate_segment(::Val{3}, daf::DAF, segment, time::Number) =
    spk_vector3(daf, segment, time)
@inline _evaluate_segment(::Val{6}, daf::DAF, segment, time::Number) =
    spk_vector6(daf, segment, time)
@inline _evaluate_segment(::Val{9}, daf::DAF, segment, time::Number) =
    spk_vector9(daf, segment, time)
@inline _evaluate_segment(::Val{12}, daf::DAF, segment, time::Number) =
    spk_vector12(daf, segment, time)

# Preserve the established low-level link evaluators for downstream users. Direction
# handling remains the caller's responsibility, as in the original methods.
@inline spk_vector3(daf::DAF, link::SPKLink, time::Number) =
    _evaluate_segment(Val(3), daf, _segment_from_link(daf, link), time)
@inline spk_vector6(daf::DAF, link::SPKLink, time::Number) =
    _evaluate_segment(Val(6), daf, _segment_from_link(daf, link), time)
@inline spk_vector9(daf::DAF, link::SPKLink, time::Number) =
    _evaluate_segment(Val(9), daf, _segment_from_link(daf, link), time)
@inline spk_vector12(daf::DAF, link::SPKLink, time::Number) =
    _evaluate_segment(Val(12), daf, _segment_from_link(daf, link), time)

@inline pck_vector3(daf::DAF, link::SPKLink, time::Number) =
    _evaluate_segment(Val(3), daf, _segment_from_link(daf, link), time)
@inline pck_vector6(daf::DAF, link::SPKLink, time::Number) =
    _evaluate_segment(Val(6), daf, _segment_from_link(daf, link), time)
@inline pck_vector9(daf::DAF, link::SPKLink, time::Number) =
    _evaluate_segment(Val(9), daf, _segment_from_link(daf, link), time)
@inline pck_vector12(daf::DAF, link::SPKLink, time::Number) =
    _evaluate_segment(Val(12), daf, _segment_from_link(daf, link), time)

@inline function _evaluate_link(
    ::Val{N}, eph::EphemerisProvider, link::SPKLink, time::Number
) where {N}
    daf = @inbounds get_daf(eph, file_id(link))
    segment = _segment_from_link(daf, link)
    return _evaluate_segment(Val(N), daf, segment, time)
end

@inline function _evaluate_point(
    ::Val{N}, eph::EphemerisProvider, from::Int, to::Int, time::Number
) where {N}
    link = _find_spk_link(eph, from, to, time)
    return factor(link) * _evaluate_link(Val(N), eph, link, time)
end

@inline function _evaluate_orientation(
    ::Val{N}, eph::EphemerisProvider, from::Int, to::Int, time::Number
) where {N}
    link = _find_pck_link(eph, from, to, time)
    return _evaluate_link(Val(N), eph, link, time)
end

"""
    ephem_vector3(eph, from, to, time)

Compute the position of `to` relative to `from` at `time`, expressed in the kernel
timescale as seconds since J2000.
"""
@inline ephem_vector3(eph::EphemerisProvider, from::Int, to::Int, time::Number) =
    _evaluate_point(Val(3), eph, from, to, time)

"""
    ephem_vector6(eph, from, to, time)

Compute position and velocity of `to` relative to `from`.
"""
@inline ephem_vector6(eph::EphemerisProvider, from::Int, to::Int, time::Number) =
    _evaluate_point(Val(6), eph, from, to, time)

"""
    ephem_vector9(eph, from, to, time)

Compute position, velocity, and acceleration of `to` relative to `from`.
"""
@inline ephem_vector9(eph::EphemerisProvider, from::Int, to::Int, time::Number) =
    _evaluate_point(Val(9), eph, from, to, time)

"""
    ephem_vector12(eph, from, to, time)

Compute position through jerk of `to` relative to `from`.
"""
@inline ephem_vector12(eph::EphemerisProvider, from::Int, to::Int, time::Number) =
    _evaluate_point(Val(12), eph, from, to, time)

"""
    ephem_rotation3(eph, from, to, time)

Compute the orientation angles of axes `to` relative to `from`.
"""
@inline ephem_rotation3(eph::EphemerisProvider, from::Int, to::Int, time::Number) =
    _evaluate_orientation(Val(3), eph, from, to, time)

"""
    ephem_rotation6(eph, from, to, time)

Compute orientation angles and their first derivatives.
"""
@inline ephem_rotation6(eph::EphemerisProvider, from::Int, to::Int, time::Number) =
    _evaluate_orientation(Val(6), eph, from, to, time)

"""
    ephem_rotation9(eph, from, to, time)

Compute orientation angles and their first two derivatives.
"""
@inline ephem_rotation9(eph::EphemerisProvider, from::Int, to::Int, time::Number) =
    _evaluate_orientation(Val(9), eph, from, to, time)

"""
    ephem_rotation12(eph, from, to, time)

Compute orientation angles and their first three derivatives.
"""
@inline ephem_rotation12(eph::EphemerisProvider, from::Int, to::Int, time::Number) =
    _evaluate_orientation(Val(12), eph, from, to, time)

function _timespan(timespan)
    timespan === nothing && return nothing
    length(timespan) == 2 || throw(ArgumentError("timespan must contain exactly two epochs."))
    first_time, last_time = timespan
    first_time <= last_time || throw(ArgumentError("timespan must be ordered."))
    return Float64(first_time), Float64(last_time)
end

function _check_composite_coverage(links::Vector{SPKLink}, timespan)
    timespan === nothing && return nothing
    first_time, last_time = timespan
    intervals = sort!(
        [(initial_time(link), final_time(link)) for link in links];
        by=first,
    )

    covered_to = first_time
    for (start_time, stop_time) in intervals
        stop_time < covered_to && continue
        start_time <= covered_to || throw(jEph.EphemerisError(
            "prepared ephemeris route has a coverage gap at $covered_to seconds since J2000."
        ))
        covered_to = max(covered_to, stop_time)
        covered_to >= last_time && return nothing
    end

    throw(jEph.EphemerisError(
        "prepared ephemeris route does not cover the requested final epoch $last_time."
    ))
end

function _direct_ephemeris_route(eph::EphemerisProvider, link::SPKLink)
    daf = @inbounds get_daf(eph, file_id(link))
    return DirectEphemerisRoute(daf, _segment_from_link(daf, link), link)
end

function _direct_orientation_route(eph::EphemerisProvider, link::SPKLink)
    daf = @inbounds get_daf(eph, file_id(link))
    return DirectOrientationRoute(daf, _segment_from_link(daf, link), link)
end

"""
    prepare_ephemeris(eph, from, to; timespan=nothing)

Resolve a point-ephemeris route once for repeated calls. When one segment covers
`timespan`, the returned route holds its concrete evaluator and bypasses provider routing.
"""
function prepare_ephemeris(
    eph::EphemerisProvider, from::Int, to::Int; timespan=nothing
)
    links = get(spk_routes(eph), (from, to), nothing)
    links === nothing && throw(jEph.EphemerisError(
        "ephemeris data for point with NAIFId $to with respect to point $from is unavailable."
    ))

    span = _timespan(timespan)
    if span !== nothing
        first_time, last_time = span
        for link in links
            initial_time(link) <= first_time && last_time <= final_time(link) &&
                return _direct_ephemeris_route(eph, link)
        end
        _check_composite_coverage(links, span)
    elseif length(links) == 1
        return _direct_ephemeris_route(eph, only(links))
    end

    return CompositeEphemerisRoute(
        eph, from, to, links, build_thread_cache(RouteCache)
    )
end

"""
    prepare_orientation(eph, from, to; timespan=nothing)

Resolve a PCK orientation route once for repeated calls.
"""
function prepare_orientation(
    eph::EphemerisProvider, from::Int, to::Int; timespan=nothing
)
    links = get(pck_routes(eph), (from, to), nothing)
    links === nothing && throw(jEph.EphemerisError(
        "ephemeris data for axes with NAIFId $to with respect to axes $from is unavailable."
    ))

    span = _timespan(timespan)
    if span !== nothing
        first_time, last_time = span
        for link in links
            initial_time(link) <= first_time && last_time <= final_time(link) &&
                return _direct_orientation_route(eph, link)
        end
        _check_composite_coverage(links, span)
    elseif length(links) == 1
        return _direct_orientation_route(eph, only(links))
    end

    return CompositeOrientationRoute(
        eph, from, to, links, build_thread_cache(RouteCache)
    )
end

@inline function _evaluate_prepared(
    ::Val{N}, route::DirectEphemerisRoute{S, F}, time::Number
) where {N, S, F}
    _covers(route.link, time) || throw(jEph.EphemerisError(
        "prepared ephemeris route is not available at $time seconds since J2000."
    ))
    value = _evaluate_segment(Val(N), route.daf, route.segment, time)
    return F == 1 ? value : -value
end

@inline function _evaluate_prepared(
    ::Val{N}, route::CompositeEphemerisRoute, time::Number
) where {N}
    link = _find_candidate(
        route.links,
        thread_cache(route.cache),
        route.from,
        route.to,
        time,
        "point",
    )
    return factor(link) * _evaluate_link(Val(N), route.eph, link, time)
end

@inline function _evaluate_prepared(
    ::Val{N}, route::DirectOrientationRoute, time::Number
) where {N}
    _covers(route.link, time) || throw(jEph.EphemerisError(
        "prepared orientation route is not available at $time seconds since J2000."
    ))
    return _evaluate_segment(Val(N), route.daf, route.segment, time)
end

@inline function _evaluate_prepared(
    ::Val{N}, route::CompositeOrientationRoute, time::Number
) where {N}
    link = _find_candidate(
        route.links,
        thread_cache(route.cache),
        route.from,
        route.to,
        time,
        "axes",
    )
    return _evaluate_link(Val(N), route.eph, link, time)
end

@inline ephem_vector3(route::AbstractPreparedEphemerisRoute, time::Number) =
    _evaluate_prepared(Val(3), route, time)
@inline ephem_vector6(route::AbstractPreparedEphemerisRoute, time::Number) =
    _evaluate_prepared(Val(6), route, time)
@inline ephem_vector9(route::AbstractPreparedEphemerisRoute, time::Number) =
    _evaluate_prepared(Val(9), route, time)
@inline ephem_vector12(route::AbstractPreparedEphemerisRoute, time::Number) =
    _evaluate_prepared(Val(12), route, time)

@inline ephem_rotation3(route::AbstractPreparedOrientationRoute, time::Number) =
    _evaluate_prepared(Val(3), route, time)
@inline ephem_rotation6(route::AbstractPreparedOrientationRoute, time::Number) =
    _evaluate_prepared(Val(6), route, time)
@inline ephem_rotation9(route::AbstractPreparedOrientationRoute, time::Number) =
    _evaluate_prepared(Val(9), route, time)
@inline ephem_rotation12(route::AbstractPreparedOrientationRoute, time::Number) =
    _evaluate_prepared(Val(12), route, time)

"""
    ephem_vector3(routes::Tuple, time)

Evaluate a static tuple of prepared ephemeris routes at one epoch. The tuple is mapped at
compile time, preserving each route's concrete segment evaluator and avoiding allocations.
The same batch form is available for `ephem_vector6`, `ephem_vector9`, and
`ephem_vector12`.
"""
@inline ephem_vector3(
    routes::Tuple{Vararg{AbstractPreparedEphemerisRoute}}, time::Number
) = map(route -> ephem_vector3(route, time), routes)
@inline ephem_vector6(
    routes::Tuple{Vararg{AbstractPreparedEphemerisRoute}}, time::Number
) = map(route -> ephem_vector6(route, time), routes)
@inline ephem_vector9(
    routes::Tuple{Vararg{AbstractPreparedEphemerisRoute}}, time::Number
) = map(route -> ephem_vector9(route, time), routes)
@inline ephem_vector12(
    routes::Tuple{Vararg{AbstractPreparedEphemerisRoute}}, time::Number
) = map(route -> ephem_vector12(route, time), routes)

@inline ephem_rotation3(
    routes::Tuple{Vararg{AbstractPreparedOrientationRoute}}, time::Number
) = map(route -> ephem_rotation3(route, time), routes)
@inline ephem_rotation6(
    routes::Tuple{Vararg{AbstractPreparedOrientationRoute}}, time::Number
) = map(route -> ephem_rotation6(route, time), routes)
@inline ephem_rotation9(
    routes::Tuple{Vararg{AbstractPreparedOrientationRoute}}, time::Number
) = map(route -> ephem_rotation9(route, time), routes)
@inline ephem_rotation12(
    routes::Tuple{Vararg{AbstractPreparedOrientationRoute}}, time::Number
) = map(route -> ephem_rotation12(route, time), routes)
