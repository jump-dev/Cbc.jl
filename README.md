# COIN-OR Branch and Cut Interface (Cbc.jl)

[![](https://www.coin-or.org/wordpress/wp-content/uploads/2014/08/COINOR.png)](https://www.coin-or.org)

`Cbc.jl` is an interface to the **[COIN-OR Branch and Cut](https://projects.coin-or.org/Cbc)**
solver. It provides a complete interface to the low-level C API, as well as an
implementation of the solver-independent `MathOptInterface`
API's

*Note: This wrapper is maintained by the JuMP community and is not a COIN-OR
project.*

[![Build Status](https://github.com/jump-dev/Cbc.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/Cbc.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/Cbc.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/Cbc.jl)

## Installation

The package can be installed with `Pkg.add`.

```
julia> import Pkg; Pkg.add("Cbc")
```

Cbc.jl will use [BinaryProvider.jl](https://github.com/JuliaPackaging/BinaryProvider.jl) to automatically install the Cbc binaries. This should work for both the official Julia binaries from `https://julialang.org/downloads/` and source-builds.

## Custom Installation

To install custom built Clp binaries set the environmental variable `JULIA_CBC_LIBRARY_PATH` and call `import Pkg; Pkg.build("Cbc")`. For instance, if the libraries are installed in `/opt/lib` just call
```julia
ENV["JULIA_CBC_LIBRARY_PATH"] = "/opt/lib"
import Pkg; Pkg.build("Cbc")
```
If you do not want BinaryProvider to download the default binaries on install set  `JULIA_CBC_LIBRARY_PATH`  before calling `import Pkg; Pkg.add("Cbc")`.

To switch back to the default binaries clear `JULIA_CBC_LIBRARY_PATH` and call `import Pkg; Pkg.build("Cbc")`.

### Using with **[JuMP](https://github.com/jump-dev/JuMP.jl)**

Use `Cbc.Optimizer` to use Cbc with JuMP:
```julia
using Cbc
using JuMP
model = Model(Cbc.Optimizer)
set_optimizer_attribute(model, "logLevel", 1)
```

Options are solver-dependent, and unfortunately not well documented.

The following options are likely to be the most useful:

* ``seconds`` -- Solution timeout limit. (Must be a ``Float64``)
* ``logLevel`` -- Set to 1 to enable solution output.
* ``maxSolutions`` -- Terminate after this many feasible solutions have been found.
* ``maxNodes`` -- Terminate after this many branch-and-bound nodes have been evaluated.
* ``allowableGap`` -- Terminate after optimality gap is less than this value (on an absolute scale).
* ``ratioGap`` -- Terminate after optimality gap is smaller than this relative fraction.
* ``threads`` -- Set the number of threads to use for parallel branch & bound

The complete list of parameters can be found by running the ``cbc`` executable and typing ``?`` at the prompt.

[Cbc]: https://projects.coin-or.org/Cbc
