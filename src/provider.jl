
export EphemerisProvider, prepare_ephemeris, prepare_orientation

"""
    EphemerisProvider(file::String)
    EphemerisProvider(files::Vector{String})

Create an `EphemerisProvider` instance by loading a single or multiple binary ephemeris
kernel files specified by `files`. Currently, only NAIF Double precision Array File (DAF)
kernels (i.e., SPK and PCK) are accepted.

### Example
```julia-repl
julia> eph = EphemerisProvider("PATH_TO_KERNEL")
EphemerisProvider([...])

julia> eph = EphemerisProvider(["PATH_TO_KERNEL_1", "PATH_TO_KERNEL_2"])
EphemerisProvider([])
```
"""
struct EphemerisProvider <: jEph.AbstractEphemerisProvider
    files::Vector{DAF}
    spklinks::SPKLinkTable
    pcklinks::SPKLinkTable
    spkroutes::FlatLinkTable
    pckroutes::FlatLinkTable
    spkcache::ThreadCache{RouteCache}
    pckcache::ThreadCache{RouteCache}
end

EphemerisProvider(files::AbstractString) = EphemerisProvider([files])
function EphemerisProvider(files::Vector{<:AbstractString})

    # Initial parsing of each DAF file
    ndafs = length(files)

    dafs = Vector{DAF}(undef, ndafs)
    @inbounds for fid = ndafs:-1:1
        dafs[fid] = DAF(files[fid])
        initialise_segments!(dafs[fid])
    end

    spklinks, pcklinks = create_linktables(dafs)
    return EphemerisProvider(
        dafs,
        spklinks,
        pcklinks,
        flatten_linktable(spklinks),
        flatten_linktable(pcklinks),
        build_thread_cache(RouteCache),
        build_thread_cache(RouteCache),
    )

end

abstract type AbstractPreparedEphemerisRoute end
abstract type AbstractPreparedOrientationRoute end

"""
    DirectEphemerisRoute

A prepared point-ephemeris route backed by one concrete SPK segment. Construct instances
with [`prepare_ephemeris`](@ref).
"""
struct DirectEphemerisRoute{S <: AbstractSPKSegment, F} <: AbstractPreparedEphemerisRoute
    daf::DAF
    segment::S
    link::SPKLink
end

function DirectEphemerisRoute(daf::DAF, segment::S, link::SPKLink) where {
    S <: AbstractSPKSegment
}
    return DirectEphemerisRoute{S, factor(link)}(daf, segment, link)
end

"""
    CompositeEphemerisRoute

A prepared point-ephemeris route that retains multiple priority-ordered segments because
no single segment covers the requested interval.
"""
struct CompositeEphemerisRoute <: AbstractPreparedEphemerisRoute
    eph::EphemerisProvider
    from::Int
    to::Int
    links::Vector{SPKLink}
    cache::ThreadCache{RouteCache}
end

"""
    DirectOrientationRoute

A prepared orientation route backed by one concrete PCK segment. Construct instances with
[`prepare_orientation`](@ref).
"""
struct DirectOrientationRoute{S <: AbstractSPKSegment} <: AbstractPreparedOrientationRoute
    daf::DAF
    segment::S
    link::SPKLink
end

"""
    CompositeOrientationRoute

A prepared orientation route that retains multiple priority-ordered PCK segments.
"""
struct CompositeOrientationRoute <: AbstractPreparedOrientationRoute
    eph::EphemerisProvider
    from::Int
    to::Int
    links::Vector{SPKLink}
    cache::ThreadCache{RouteCache}
end

function Base.show(io::IO, eph::EphemerisProvider)
    print(io, "$(length(get_daf(eph)))-kernel EphemerisProvider")
end

function Base.show(io::IO, ::MIME"text/plain", eph::EphemerisProvider)
    println(io, eph, ":")
    for daf in get_daf(eph)
        println(io, " \"", filepath(daf), "\"")
    end
end

"""
    daf(eph::EphemerisProvider)

Return the [`DAF`](@ref) files stored in the ephemeris provider.
"""
@inline get_daf(eph::EphemerisProvider) = eph.files

"""
    daf(eph::EphemerisProvider, id::Int)

Return the [`DAF`](@ref) file in the ephemeris provider at index `id`.
"""
@inline get_daf(eph::EphemerisProvider, id::Int) = get_daf(eph)[id]

"""
    spk_links(eph::EphemerisProvider)

Return the [`SPKLinkTable`] for the SPK segments.
"""
@inline spk_links(eph::EphemerisProvider) = eph.spklinks

"""
    pck_links(eph::EphemerisProvider)

Return the [`SPKLinkTable`](@ref) for the PCK segments.
"""
@inline pck_links(eph::EphemerisProvider) = eph.pcklinks

@inline spk_routes(eph::EphemerisProvider) = eph.spkroutes
@inline pck_routes(eph::EphemerisProvider) = eph.pckroutes
