[![](https://www.coin-or.org/wordpress/wp-content/uploads/2014/08/COINOR.png)](https://www.coin-or.org)

# Cbc.jl

[![Build Status](https://github.com/jump-dev/Cbc.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/Cbc.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/Cbc.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/Cbc.jl)

[Cbc.jl](https://github.com/jump-dev/Cbc.jl) is a wrapper for the [COIN-OR Branch and Cut (Cbc)](https://projects.coin-or.org/Cbc)
solver.

The wrapper has two components:

 * a thin wrapper around the complete C API
 * an interface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)

## Affiliation

This wrapper is maintained by the JuMP community and is not a COIN-OR project.

## License

`Cbc.jl` is licensed under the [MIT License](https://github.com/jump-dev/Cbc.jl/blob/master/LICENSE.md).

The underlying solver, [coin-or/Cbc](https://github.com/coin-or/Cbc), is
licensed under the [Eclipse public license](https://github.com/coin-or/Cbc/blob/master/LICENSE).

## Installation

Install Cbc using `Pkg.add`:
```julia
import Pkg
Pkg.add("Cbc")
```

In addition to installing the Cbc.jl package, this will also download and
install the Cbc binaries. You do not need to install Cbc separately.

To use a custom binary, read the [Custom solver binaries](https://jump.dev/JuMP.jl/stable/developers/custom_solver_binaries/)
section of the JuMP documentation.

## Use with JuMP

To use Cbc with JuMP, use `Cbc.Optimizer`:
```julia
using JuMP, Cbc
model = Model(Cbc.Optimizer)
set_attribute(model, "logLevel", 1)
```

## MathOptInterface API

The COIN Branch-and-Cut (Cbc) optimizer supports the following constraints and attributes.

List of supported objective functions:

 * [`MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}`](@ref)

List of supported variable types:

 * [`MOI.Reals`](@ref)

List of supported constraint types:

 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.EqualTo{Float64}`](@ref)
 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.GreaterThan{Float64}`](@ref)
 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.Interval{Float64}`](@ref)
 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.LessThan{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.EqualTo{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.GreaterThan{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.Integer`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.Interval{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.LessThan{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.ZeroOne`](@ref)
 * [`MOI.VectorOfVariables`](@ref) in [`MOI.SOS1{Float64}`](@ref)
 * [`MOI.VectorOfVariables`](@ref) in [`MOI.SOS2{Float64}`](@ref)

List of supported model attributes:

 * [`MOI.ObjectiveSense()`](@ref)

## Options

Options are, unfortunately, not well documented.

The following options are likely to be the most useful:

| Parameter      | Example | Explanation                                       |
| -------------- | ------- | ------------------------------------------------- |
| `seconds`      | `60.0`  | Solution timeout limit                            |
| `logLevel`     | `2`     | Set to 0 to disable solution output               |
| `maxSolutions` | `1`     | Terminate after this many feasible solutions have been found |
| `maxNodes`     | `1`     | Terminate after this many branch-and-bound nodes have been evaluated |
| `allowableGap` | `0.05`  | Terminate after optimality gap is less than this value (on an absolute scale) |
| `ratioGap`     | `0.05`  | Terminate after optimality gap is smaller than this relative fraction |
| `threads`      | `1`     | Set the number of threads to use for parallel branch & bound |

The complete list of parameters can be found by running the `cbc` executable and
typing `?` at the prompt.

Start the `cbc` executable from Julia as follows:
```julia
using Cbc_jll
Cbc_jll.cbc() do exe
    run(`$(exe)`)
end
```
