using Downloads
using Random

const DE440_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"
const PA440_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc"

const DEFAULT_DE440_PATH = joinpath(@__DIR__, "data", "de440.bsp")
const DEFAULT_PA440_PATH = joinpath(@__DIR__, "data", "moon_pa_de440_200625.bpc")
const PA440_AXES_ID = 31008
const BENCHMARK_SEED = 440

const POSITION_PAIRS = (
    (3, 399),
    (3, 301),
    (0, 10),
    (0, 5),
    (0, 6),
)

function resolve_kernel_path(env_var::AbstractString, default_path::AbstractString, url::AbstractString)
    if haskey(ENV, env_var)
        path = abspath(ENV[env_var])
        isfile(path) || error("$(env_var) does not exist: $(path)")
        return path
    end

    mkpath(dirname(default_path))
    if !isfile(default_path)
        println("Downloading $(basename(default_path)) to $(default_path)")
        Downloads.download(url, default_path)
    end

    return default_path
end

resolve_de440_path() = resolve_kernel_path("EPHEMERIDES_DE440_PATH", DEFAULT_DE440_PATH, DE440_URL)
resolve_pa440_path() = resolve_kernel_path("EPHEMERIDES_PA440_PATH", DEFAULT_PA440_PATH, PA440_URL)

function build_reference_provider()
    de440_path = resolve_de440_path()
    pa440_path = resolve_pa440_path()
    eph = EphemerisProvider([de440_path, pa440_path])

    return eph, de440_path, pa440_path
end

benchmark_position_times(n::Int = 1024) = collect(range(0.0, step = 86400.0, length = n))
benchmark_rotation_times(n::Int = 1024) = collect(range(0.0, step = 21600.0, length = n))

function shuffled_copy(values)
    shuffled = copy(values)
    Random.seed!(BENCHMARK_SEED)
    shuffle!(shuffled)
    return shuffled
end
