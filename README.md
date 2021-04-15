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

```julia
import Pkg
Pkg.add("Cbc")
```

In addition to installing the Cbc.jl package, this will also download and
install the Cbc binaries. (You do not need to install Cbc separately.) If you
require a custom build of Cbc, see the **Custom Installation** instructions
below.

## Using with JuMP

Use `Cbc.Optimizer` to use Cbc with [JuMP](https://github.com/jump-dev/JuMP.jl):
```julia
using Cbc
using JuMP
model = Model(Cbc.Optimizer)
set_optimizer_attribute(model, "logLevel", 1)
```

## Options

Options are, unfortunately, not well documented.

The following options are likely to be the most useful:

* `seconds` -- Solution timeout limit.

    For example, `set_optimizer_attribute(model, "seconds", 60.0)`.

* `logLevel` -- Set to 1 to enable solution output.

    For example, `set_optimizer_attribute(model, "logLevel", 1)`.

* `maxSolutions` -- Terminate after this many feasible solutions have been found.

    For example, `set_optimizer_attribute(model, "maxSolutions", 1)`.

* `maxNodes` -- Terminate after this many branch-and-bound nodes have been evaluated.

    For example, `set_optimizer_attribute(model, "maxNodes", 1)`.

* `allowableGap` -- Terminate after optimality gap is less than this value (on an absolute scale).

    For example, `set_optimizer_attribute(model, "allowableGap", 0.05)`.

* `ratioGap` -- Terminate after optimality gap is smaller than this relative fraction.

    For example, `set_optimizer_attribute(model, "allowableGap", 0.05)`.

* `threads` -- Set the number of threads to use for parallel branch & bound.

    For example, `set_optimizer_attribute(model, "threads", 2)`.

The complete list of parameters can be found by running the `cbc` executable and
typing `?` at the prompt.

On Julia 1.3 and above, you can start the `cbc` executable from Julia as follows:
```julia
using Cbc_jll
Cbc_jll.cbc() do exe
    run(`$(exe)`)
end
```

## Custom Installation

To install custom built Cbc binaries, use the environmental variable
`JULIA_CBC_LIBRARY_PATH` to point to the path at which you installed Cbc (the
folder containing `libCbcSolver`). For example, on Mac, after installing Cbc
with `brew install cbc`, use:
```julia
ENV["JULIA_CBC_LIBRARY_PATH"] = "/usr/local/Cellar/cbc/2.10.5/lib"
import Pkg
Pkg.add("Cbc")
Pkg.build("Cbc")
```
Replace `"/usr/local/Cellar/cbc/2.10.5/lib"` with a different path as
appropriate.

**You must have `JULIA_CBC_LIBRARY_PATH` set _every_ time you run `using Cbc`,
not just when you install it.**

Switch back to the default binaries as follows:
```julia
delete!(ENV, "JULIA_CBC_LIBRARY_PATH")
import Pkg
Pkg.build("Cbc")
```
