
"""
    SPKLink

A link object to create a mapping between [`DAFSegmentDescriptor`](@ref) and its actual
location within an [`EphemerisProvider`](@ref) object.

### Fields
- `desc` -- `DAFSegmentDescriptor` for the segment associated to this link
- `fid` -- `Int` index of the DAF containg the link data.
- `lid` -- `Int` field number in the [`SPKSegmentList`](@ref) for this segment type.
- `eid` -- `Int` index of the inner segment list that stores this SPK segment.
- `fct` -- `Int` 1 or -1 depending on whether the (from, to) directions must be reversed.

### See Also
See also [`SPKLinkTable`](@ref), [`SPKSegmentList`](@ref) and [`add_spklinks!`](@ref).
"""
struct SPKLink
    desc::DAFSegmentDescriptor
    fid::Int
    lid::Int
    eid::Int
    fct::Int
end

"""
    descriptor(link::SPKLink)

Return the SPK/PCK segment descriptor associated to this link.
"""
@inline descriptor(link::SPKLink) = link.desc

"""
    file_id(link::SPKLink)

Return the DAF file index.
"""
@inline file_id(link::SPKLink) = link.fid

"""
    list_id(link::SPKLink)

Return the index of the list containing the segments of the given SPK/PCK type.
"""
@inline list_id(link::SPKLink) = link.lid

"""
    element_id(link::SPKLink)

Return the segment index in the inner SPK/PCK segment list.
"""
@inline element_id(link::SPKLink) = link.eid

"""
    factor(link::SPKLink)

Return the direction multiplicative factor.
"""
@inline factor(link::SPKLink) = link.fct

"""
    initial_time(link::SPKLink)

Return the initial epoch of the interval for which ephemeris data are contained in the
segment associated to this link, in seconds since J2000.0
"""
@inline initial_time(link::SPKLink) = initial_time(descriptor(link))

"""
    final_time(link::SPKLink)

Return the final epoch of the interval for which ephemeris data are contained in the
segment associated to this link, in seconds since J2000.0
"""
@inline final_time(link::SPKLink) = final_time(descriptor(link))

"""
    reverse_link(link::SPKLink)

Reverse the sign, i.e. change the sign of the multiplicative factor, of the link.
"""
function reverse_link(link::SPKLink)
    SPKLink(
        descriptor(link), file_id(link), list_id(link), element_id(link), -factor(link)
    )
end


"""
    SPKLinkTable

Dictionary object providing all the [`SPKLink`](@ref) available between a set of (from, to)
objects
"""
const SPKLinkTable = Dict{Int, Dict{Int, Vector{SPKLink}}}

"""
    FlatLinkTable

Internal routing table keyed by `(from, to)`. Unlike [`SPKLinkTable`](@ref), this layout
requires only one dictionary lookup on an uncached query.
"""
const FlatLinkTable = Dict{Tuple{Int, Int}, Vector{SPKLink}}

const ROUTE_CACHE_SIZE = 8

"""
    RouteCache

Small direct-mapped cache for the last active links used by a Julia thread. A cache miss
falls back to `FlatLinkTable`; collisions only affect performance, never precedence.
"""
mutable struct RouteCache
    from::Vector{Int}
    to::Vector{Int}
    priority::Vector{Int}
    links::Vector{Union{Nothing, SPKLink}}
end

RouteCache() = RouteCache(
    fill(typemin(Int), ROUTE_CACHE_SIZE),
    fill(typemin(Int), ROUTE_CACHE_SIZE),
    zeros(Int, ROUTE_CACHE_SIZE),
    Union{Nothing, SPKLink}[nothing for _ in 1:ROUTE_CACHE_SIZE],
)

@inline function route_cache_slot(from::Int, to::Int)
    from_bits = reinterpret(UInt, from)
    to_bits = reinterpret(UInt, to)
    # Fold higher target-ID bits into the low bits used by the power-of-two slot mask.
    # Rotating small positive NAIF IDs by a large amount leaves those low bits empty and
    # makes every route with the same center collide.
    mixed = xor(xor(from_bits * UInt(0x9e3779b97f4a7c15), to_bits), to_bits >> 3)
    return Int(mixed & UInt(ROUTE_CACHE_SIZE - 1)) + 1
end

"""
    flatten_linktable(table::SPKLinkTable)

Create the single-lookup `(from, to)` representation of a legacy link table while
retaining each candidate vector's NAIF priority order.
"""
function flatten_linktable(table::SPKLinkTable)
    flat = FlatLinkTable()
    sizehint!(flat, sum(length, values(table); init=0))

    for (to, from_map) in table
        for (from, candidates) in from_map
            flat[(from, to)] = candidates
        end
    end

    return flat
end

"""
    create_linktables(dafs::Vector{DAF})

Create the SPK and PCK [`SPKLinkTable`](@ref) for all the segments stores in the input DAFs.
"""
function create_linktables(dafs::Vector{DAF})

    spklinks, pcklinks = SPKLinkTable(), SPKLinkTable()
    for j in reverse(eachindex(dafs))
        add_spklinks!(is_spk(dafs[j]) ? spklinks : pcklinks, dafs[j], j)
    end

    return spklinks, pcklinks

end

"""
    add_spklinks!(table::SPKLinkTable, daf::DAF, fid::Int)

Insert in the input [`SPKLinkTable`](@ref) all the SPK or PCK links associated to
the segment descriptors of the input DAF.
"""
function add_spklinks!(table::SPKLinkTable, daf::DAF, fid::Int)

    # Initialise the number of elements contained in each list
    nfields = fieldcount(SPKSegmentList)
    counter = zeros(nfields)

    for desc in descriptors(daf)

        segtype = segment_type(desc)

        # Get the field index of the list for this segment type
        lid = SPK_SEGMENTLIST_MAPPING[segtype]

        counter[lid] += 1

        # Create the forward and backward SPKLink if not already available
        # The forward link always goes from target to center!
        f_map = get!(table, target(desc), Dict{Int, Vector{SPKLink}}())

        f_link = SPKLink(desc, fid, lid, counter[lid], 1)
        push!(get!(f_map, center(desc), SPKLink[]), f_link)

        # Populate with both the forward and backward links, initialising the
        # SPKLink key if a link between the two bodies was yet to be found.
        # NOTE: the backward link is created only for SPK files
        if is_spk(daf)
            b_map = get!(table, center(desc), Dict{Int, Vector{SPKLink}}())
            push!(get!(b_map, target(desc), SPKLink[]), reverse_link(f_link))
        end
    end

end
