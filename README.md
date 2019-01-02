# COIN-OR Branch and Cut Interface (Cbc.jl)

[![](https://www.coin-or.org/wordpress/wp-content/uploads/2014/08/COINOR.png)](https://www.coin-or.org)

`Cbc.jl` is an interface to the **[COIN-OR Branch and Cut](https://projects.coin-or.org/Cbc)**
solver. It provides a complete interface to the low-level C API, as well as an
implementation of the solver-independent `MathProgBase` and `MathOptInterface`
API's.   

*Note: This wrapper is maintained by the JuliaOpt community and is not a COIN-OR
project.*


[![Build Status](https://travis-ci.org/JuliaOpt/Cbc.jl.svg?branch=master)](https://travis-ci.org/JuliaOpt/Cbc.jl)

## Installation

The package is registered in `METADATA.jl` and so can be installed with `Pkg.add`.

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

### Using with **[MathProgBase]**

Cbc provides a solver object that can be passed to ``mixintprog`` in MathProgBase (and used to create instances of the solver-independent ``AbstractMathProgModel`` type):

```julia
using Cbc
using MathProgBase
mixintprog(..., CbcSolver(Option1=value1,Option2=value2,...))
```

see the MathProgBase documentation for further information.

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

In addition, we provide the julia-specific option ``check_warmstart`` which, if set to ``false``, will tell the wrapper to pass along the warmstart solution regardless of if it satisfies the constraints of the problem. The default value is ``true``.

### Using the C interface

The low-level C interface is available in the ``CbcCInterface`` submodule:
```julia
using Cbc.CbcCInterface
```

Using this interface is not recommended.

[Cbc]: https://projects.coin-or.org/Cbc
[MathProgBase]: https://github.com/JuliaOpt/MathProgBase.jl
