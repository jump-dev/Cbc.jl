
using CoinMP

function knapsack()
    f = Float64[3, 2, 3, 1, 2]
    A = [7. 2. 9. 3. 1.]
    types = fill(INTEGER, 5)
    capacity = 10.
    (z, x, flag) = mixintprog(-f, A, [], [capacity], [], ones(5), types; LogLevel=0)
    println("Solution status: $flag")
    println("Optimal value: $z")
    println("Solution vector: $x") # should be -7

end

knapsack()

