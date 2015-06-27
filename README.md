Cbc.jl
=========

[![Build Status](https://travis-ci.org/JuliaOpt/Cbc.jl.svg?branch=master)](https://travis-ci.org/JuliaOpt/Cbc.jl)

This is a Julia interface to the mixed-integer linear programming solver **[Cbc]**. Cbc is a high-performance open-source solver. A basic high-level ``mixintprog`` interface is provided through the **[MathProgBase]** package.

The supported platforms are Linux, OS X, and Windows. Binaries are provided for Windows and OS X, and will be installed by default. On Linux, Cbc will be automatically compiled from source. Ensure that a C++ compiler is installed first; on Debian-based systems, install the ``build-essential`` package.

[Cbc]: https://projects.coin-or.org/Cbc


### Using with **[MathProgBase]**


Cbc provides a solver object that can be passed to ``mixintprog`` in MathProgBase (and used to create instances of the solver-independent ``AbstractMathProgModel`` type):

```julia
    using Cbc
    using MathProgBase
    mixintprog(..., CbcSolver(Option1=value1,Option2=value2,...))
```

see the MathProgBase documentation for further information.

[MathProgBase]: https://github.com/JuliaOpt/MathProgBase.jl

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

### Using the C interface

The low-level C interface is available in the ``CbcCInterface`` submodule:
```
    using Cbc.CbcCInterface
```

Using this interface is not recommended.

