module Cbc

if haskey(ENV,"JULIA_CBC_LIBRARY_PATH") || VERSION < v"1.3"
    deps_file = joinpath(dirname(@__DIR__), "deps", "deps.jl")
    if isfile(deps_file)
        using Libdl
        include(deps_file)
    else
        error("Cbc not properly installed. Please run import `Pkg; Pkg.build(\"Cbc\")`.")
    end
else
    import Cbc_jll: libcbcsolver
end

using CEnum

include("gen/ctypes.jl")
include("gen/libcbc_common.jl")
include("gen/libcbc_api.jl")

const _CBC_VERSION = VersionNumber(unsafe_string(Cbc_getVersion()))

if !(v"2.10.0" <= _CBC_VERSION <= v"2.10.5")
    error(
        "You have installed version $_CBC_VERSION of Cbc, which is not " *
        "supported by Cbc.jl. If the version change was breaking, changes " *
        "will need to be made to the Julia code. Please open an issue at " *
        "https://github.com/JuliaOpt/Cbc.jl."
    )
end

include("MOI_wrapper.jl")

# TODO(odow): remove at Cbc.jl v1.0.0.
function CbcSolver(args...; kwargs...)
    error(
        "`CbcSolver` is no longer supported. If you are using JuMP, upgrade " *
        "to the latest version and use `Cbc.Optimizer` instead. If you are " *
        "using MathProgBase (e.g., via `lingprog`), you will need to upgrade " *
        "to MathOptInterface (https://github.com/JuliaOpt/MathOptInterface.jl)."
    )
end
export CbcSolver

end
