using MathProgBase
using Cbc
using Test

function knapsack()
    f = Float64[3, 2, 3, 1, 2]
    A = [7. 2. 9. 3. 1.]
    capacity = 10.
    solution = mixintprog(-f, A, '<', capacity, :Int, 0, 1, CbcSolver(logLevel=0))
    println("Solution status: ", solution.status)
    println("Optimal value: ", solution.objval)
    print("Solution vector: ") 
    show(solution.sol)
    println()
    @test solution.objval == -7
end

knapsack()

