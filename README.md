[![](https://www.coin-or.org/wordpress/wp-content/uploads/2014/08/COINOR.png)](https://www.coin-or.org)

# Cbc.jl

[![Build Status](https://github.com/jump-dev/Cbc.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/Cbc.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/Cbc.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/Cbc.jl)

`Cbc.jl` is a wrapper for the [COIN-OR Branch and Cut (Cbc)](https://projects.coin-or.org/Cbc)
solver.

The wrapper has two components:
 * a thin wrapper around the complete C API
 * an interface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface)

*Note: This wrapper is maintained by the JuMP community and is not a COIN-OR
project.*

## Installation

Install Cbc.jl using `Pkg.add`:
```julia
import Pkg; Pkg.add("Cbc")
```

In addition to installing the Cbc.jl package, this will also download and
install the Cbc binaries. (You do not need to install Cbc separately.)

To use a custom binary, read the [Custom solver binaries](https://jump.dev/JuMP.jl/stable/developers/custom_solver_binaries/)
section of the JuMP documentation.

## Use with JuMP

To use Cbc with [JuMP](https://github.com/jump-dev/JuMP.jl), use `Cbc.Optimizer`:
```julia
using JuMP, Cbc
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

    For example, `set_optimizer_attribute(model, "ratioGap", 0.05)`.

* `threads` -- Set the number of threads to use for parallel branch & bound.

    For example, `set_optimizer_attribute(model, "threads", 2)`.

The complete list of parameters can be found by running the `cbc` executable and
typing `?` at the prompt.

You can start the `cbc` executable from Julia as follows:
```julia
using Cbc_jll
Cbc_jll.cbc() do exe
    run(`$(exe)`)
end
```
