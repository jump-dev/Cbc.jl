
using MathProgBase

function knapsack()
    f = Float64[3, 2, 3, 1, 2]
    A = [7. 2. 9. 3. 1.]
    capacity = 10.
    solution = mixintprog(-f, A, '<', capacity, 'I', 0, 1, MIPSolver(:Cbc,LogLevel=0))
    println("Solution status: ", solution.status)
    println("Optimal value: ", solution.objval) # should be -7
    print("Solution vector: ") 
    show(solution.sol)
    println()

end

knapsack()

