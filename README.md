Cbc.jl
=========

This is a somewhat complete Julia interface to the mixed-integer linear programming solver **[Cbc]** via the **[CoinMP]** C library. Cbc is a high-performance open-source solver. A basic high-level ``mixintprog`` interface is provided through the **[MathProgBase]** package. Solve options and callbacks are available but unforunately undocumented due to the absence of documentation for CoinMP. 

This package should work on Linux, OS X, and Windows, but it has not been widely tested. Binaries are provided for Windows (Vista-7-8). Windows users must install the Visual Studio **[redistributable]** package. OS X users will need a proper build environment. Please report any issues.

[Cbc]: https://projects.coin-or.org/Cbc
[CoinMP]: https://projects.coin-or.org/CoinMP
[redistributable]: http://www.microsoft.com/en-us/download/details.aspx?id=30679
[MathProgBase]: https://github.com/mlubin/MathProgBase.jl
