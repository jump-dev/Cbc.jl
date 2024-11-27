# Copyright (c) 2013: Cbc.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module Cbc

import Cbc_jll: libcbcsolver
import MathOptInterface as MOI
import SparseArrays

function __init__()
    version_str = unsafe_string(Cbc_getVersion())
    version = if version_str == "devel"
        # Support un-released versions of Cbc. These may differ in C API
        # compatibility! Use at your own peril.
        v"2.10.5"
    else
        VersionNumber(version_str)
    end
    if !(v"2.10.0" <= version < v"2.11")
        error("""
        You have installed version $version of Cbc, which is not supported by
        Cbc.jl We require Cbc version 2.10. After installing Cbc 2.10, run:

            import Pkg
            Pkg.rm("Cbc")
            Pkg.add("Cbc")

        If you have a newer version of Cbc installed, changes may need to be made
        to the Julia code. Please open an issue at
        https://github.com/jump-dev/Cbc.jl.
        """)
    end
    return
end

include("gen/libcbc_common.jl")
include("gen/libcbc_api.jl")
include("MOI_wrapper/MOI_wrapper.jl")

# Cbc exports all `Cbc_xxx` symbols. If you don't want all of these symbols in
# your environment, then use `import Cbc` instead of `using Cbc`.

for sym in filter(s -> startswith("$s", "Cbc_"), names(@__MODULE__, all = true))
    @eval export $sym
end

end
