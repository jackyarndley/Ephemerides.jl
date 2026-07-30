"""
    DAF_RECORD_LENGTH

DAF record length, in bytes.

### References
- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)
"""
const DAF_RECORD_LENGTH = 1024;


"""
    get_record(array::Vector{UInt8}, index::Integer)

Retrieve a whole DAF record at position `index`.
"""
@inline @views function get_record(array, index::Integer)
    @inbounds return array[1+DAF_RECORD_LENGTH*(index-1):DAF_RECORD_LENGTH*index]
end


# Parsing Utils
# =====================

"""
    is_little_endian(array::Vector{UInt8})

Return true if the array corresponds to the string indicating a little-endian format.
"""
function is_little_endian(array::Vector{UInt8})
    # TODO: add support for VAX-GFLT and VAX-DFLT
    endian = get_string(array, 88, 8);
    if endian == "LTL-IEEE"
        return true
    elseif endian == "BIG-IEEE"
        return false
    else
        throw(ErrorException("The endiannes could not be recognised!"))
    end
end


# Read from array the string at address with given length
function get_string(array, address::Integer, bytes::Integer)
    # address is in 0-index notation!
    @inbounds rstrip(String(@view(array[address+1:address+bytes])))
end

@inline get_num(x::Number, lend::Bool) = lend ? htol(x) : hton(x)

@inline function get_int(array, address::Integer, lend::Bool)
    # address is in 0-index notation!
    GC.@preserve array begin
        value = unsafe_load(Ptr{Int32}(pointer(array, address + 1)))
        return get_num(value, lend)
    end
end

@inline function get_float(array, address::Integer, lend::Bool)
    # address is in 0-index notation!
    GC.@preserve array begin
        value = unsafe_load(Ptr{Float64}(pointer(array, address + 1)))
        return get_num(value, lend)
    end
end

"""
    read_doubles_rowmajor!(dest, bytes, address, nrows, ncols, little_endian)

Read a row-major block of IEEE `Float64` values from a mapped DAF byte array into the
column-major matrix `dest`. `address` uses the DAF parser's zero-based byte convention.
"""
@inline function read_doubles_rowmajor!(
    dest::Matrix{Float64},
    bytes::Vector{UInt8},
    address::Int,
    nrows::Int,
    ncols::Int,
    little_endian::Bool,
)
    GC.@preserve bytes begin
        source = Ptr{Float64}(pointer(bytes, address + 1))
        @inbounds for row in 1:nrows
            offset = (row - 1)*ncols
            for column in 1:ncols
                dest[row, column] =
                    get_num(unsafe_load(source, offset + column), little_endian)
            end
        end
    end
    return dest
end


# Vector Utils
# =====================

vhat(u) = u/vnorm(u)
vnorm(u) = sqrt(vdot(u, u))
vsep(u, v) = acos(min(1, max(-1, vdot(vhat(u), vhat(v)))))

@inbounds vdot(u, v) = u[1]*v[1] + u[2]*v[2] + u[3]*v[3]
function vcross(u, v)
    return @inbounds SVector{3}(
        u[2]*v[3]-u[3]*v[2], u[3]*v[1]-u[1]*v[3], u[1]*v[2]-u[2]*v[1]
    )
end

function vrot(v1, k, angle)

    s, c = sincos(angle)
    u = vcross(k, v1)

    # Apply Rodrigues rotation formula to rotate `v1` around `k`.
    return v1 + (1-c)*vcross(k, u) + s*u

end
