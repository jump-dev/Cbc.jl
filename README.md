CoinMP.jl
=========

This is a (still incomplete) interface to the Coin-OR mixed-integer linear programming solver **[Cbc]** via the **[CoinMP]** C library. We currently provide a basic high-level ``mixintprog`` interface (see in-line comments). Solve options are not yet available (patches are welcome). This package *should* work on Linux, OS X, and Windows (!), but it has not been widely tested.

[Cbc]: https://projects.coin-or.org/Cbc
[CoinMP]: https://projects.coin-or.org/CoinMP
