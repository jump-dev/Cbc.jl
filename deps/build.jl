
const DEPS_FILE = joinpath(@__DIR__, "deps.jl")

if isfile(DEPS_FILE)
    rm(DEPS_FILE)
end

using Libdl

function get_paths(path)
    extension = if Sys.isapple()
        ".dylib"
    elseif Sys.iswindows()
        Libdl.dlext
    elseif Sys.isunix()
        ".so"
    else
        error("System not supported.")
    end
    libcbcsolver = escape_string(joinpath(path, "libCbcSolver" * extension))
    libcbc       = escape_string(joinpath(path, "libCbc" * extension))
    return libcbcsolver, libcbc
end

function write_depsfile(libcbcsolver, libcbc)
    open(DEPS_FILE, "w") do io
        println(io, "const libcbcsolver = \"$(libcbcsolver)\"")
        println(io, "const libCbc = \"$(libcbc)\"")
    end
end

path = get(ENV, "JULIA_CBC_LIBRARY_PATH", nothing)
@show path
if path === nothing
    if VERSION >= v"1.3"
        # Use Cbc_jll instead.
        exit()
    else
        # Use version from BinaryProvider.jl
        include(joinpath(@__DIR__, "build_binary_provider.jl"))
    end
else
    # use the version from JULIA_CBC_LIBRARY_PATH!
    libcbcsolver, libcbc = get_paths(path)
    d = Libdl.dlopen_e(libcbcsolver)
    if d != C_NULL
        write_depsfile(libcbcsolver, libcbc)
    else
        error("""
        Unable to build Cbc.jl because your environment variable `JULIA_CBC_LIBRARY_PATH`
        points to an invalid location. It is currently:

        $(ENV["JULIA_CBC_LIBRARY_PATH"])

        Set it to point to the directory containing `libCbcSolver`, or to install the
        default binaries, run:

        delete!(ENV, "JULIA_CBC_LIBRARY_PATH")
        import Pkg
        Pkg.build("Cbc")

        then restart Julia for the changes to take effect.
        """)
    end
end

