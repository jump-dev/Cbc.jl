CoinMP.jl
=========

This is a somewhat complete Julia interface to the mixed-integer linear programming solver **[Cbc]** via the **[CoinMP]** C library. Cbc is a high-performance open-source solver. We currently provide a basic high-level ``mixintprog`` interface (see in-line comments). Solve options and callbacks are available but unforunately undocumented due to the absence of documentation for CoinMP. This package *should* work on Linux, OS X, and Windows (!), but it has not been widely tested.

[Cbc]: https://projects.coin-or.org/Cbc
[CoinMP]: https://projects.coin-or.org/CoinMP
