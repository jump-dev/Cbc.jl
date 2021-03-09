import Clang
import Cbc_jll
import CoinUtils_jll

FILES = [
    (CoinUtils_jll.artifact_dir, "Coin_C_defines.h"),
    (Cbc_jll.artifact_dir, "Cbc_C_Interface.h"),
]

for (dir, file) in FILES
    cp(
        joinpath(dir, "include", "coin", file),
        joinpath(@__DIR__, file);
        force = true,
    )
end

GEN_DIR = joinpath(@__DIR__, "..", "src", "gen")

wc = Clang.init(
    headers = map(i -> joinpath(@__DIR__, i[2]), FILES),
    output_file = joinpath(GEN_DIR, "libcbc_api.jl"),
    common_file = joinpath(GEN_DIR, "libcbc_common.jl"),
    header_wrapped = (root, current) -> root == current,
    header_library = x -> "libcbcsolver",
    clang_diagnostics = true,
)

run(wc)

rm(joinpath(GEN_DIR, "LibTemplate.jl"))
