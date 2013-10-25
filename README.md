Cbc.jl
=========

This is a somewhat complete Julia interface to the mixed-integer linear programming solver **[Cbc]** via the **[CoinMP]** C library. Cbc is a high-performance open-source solver. A basic high-level ``mixintprog`` interface is provided through the **[MathProgBase]** package.

The supported platforms are Linux, OS X, and Windows. Binaries are provided for Windows (Vista-7-8, 32-bit only) and OS X, and will be installed by default. Windows users must additionally install the Visual Studio **[redistributable]**. On Linux, CoinMP will be automatically compiled from source. Ensure that a C++ compiler is installed first; on Debian-based systems, install the ``build-essential`` package. 

[Cbc]: https://projects.coin-or.org/Cbc
[CoinMP]: https://projects.coin-or.org/CoinMP
[redistributable]: http://www.microsoft.com/en-us/download/details.aspx?id=30679


### Using with **[MathProgBase]**


Cbc provides a solver object that can be passed to ``mixintprog`` in MathProgBase (and used to create instances of the solver-independent ``AbstractMathProgModel`` type):

```julia
    using Clp
    using MathProgBase
    mixintprog(..., ClpSolver(Option1=value1,Option2=value2,...))
```

see the MathProgBase documentation for further information.

[MathProgBase]: https://github.com/mlubin/MathProgBase.jl

Options are solver-dependent, and unfortunately not well documented.
The following options are likely to be the most useful:

* ``MipMaxSeconds`` -- Solution timeout limit. (Must be a ``Float64``)
* ``LogLevel`` -- Set to 1 to enable solution output.
* ``MipMaxSolutions`` -- Terminate after this many feasible solutions have been found.
* ``MipMaxNodes`` -- Terminate after this many branch-and-bound nodes have been evaluated.
* ``MipAllowableGap`` -- Terminate after optimality gap is less than this value (on an absolute scale).
* ``MipFractionalGap`` -- Terminate after optimality gap is smaller than this relative fraction.

### Using the C interface

The low-level C interface is available in the ``CoinMPInterface`` submodule:
```
    using Cbc.CoinMPInterface
```

Using this interface is not recommended.

