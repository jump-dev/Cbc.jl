# TODO(odow):
#
# This script can be used to build the C interface to Cbc. However, it requires
# you to manually do the following steps first:
#
# 1) Copy Cbc_C_Interface.h from cbc into this /scripts directory
# 2) Copy Coin_C_defines.h from coinutils into this /scripts directory
#
# It should be possible to build the wrapper using the jll's, but I couldn't
# figure out how to do that, and I didn't have enough time to spend on it.

import Clang

LIBCLP_HEADERS = [
    joinpath(@__DIR__, "Coin_C_defines.h"),
    joinpath(@__DIR__, "Cbc_C_Interface.h"),
]

wc = Clang.init(
    headers = LIBCLP_HEADERS,
    output_file = joinpath(@__DIR__, "..", "src", "libcbc_api.jl"),
    common_file = joinpath(@__DIR__, "..", "src", "libcbc_common.jl"),
    header_wrapped = (root, current) -> root == current,
    header_library = x -> "libcbcsolver",
    clang_diagnostics = true,
)

run(wc)
