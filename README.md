Cbc.jl
=========

This is a somewhat complete Julia interface to the mixed-integer linear programming solver **[Cbc]** via the **[CoinMP]** C library. Cbc is a high-performance open-source solver. A basic high-level ``mixintprog`` interface is provided through the **[MathProgBase]** package. Solve options and callbacks are available but unforunately undocumented due to the absence of documentation for CoinMP. 

The supported platforms are Linux, OS X, and Windows. Binaries are provided for Windows (Vista-7-8, 32-bit only) and OS X, and will be installed by default. Windows users must additionally install the Visual Studio **[redistributable]**. OS X users will be prompted to install the ``Homebrew`` package. On Linux, CoinMP will be automatically compiled from source. If you are on a Debian-like system, ensure a C++ compiler is installed first, e.g. ``sudo apt-get install build-essential``.

[Cbc]: https://projects.coin-or.org/Cbc
[CoinMP]: https://projects.coin-or.org/CoinMP
[redistributable]: http://www.microsoft.com/en-us/download/details.aspx?id=30679
[MathProgBase]: https://github.com/mlubin/MathProgBase.jl
