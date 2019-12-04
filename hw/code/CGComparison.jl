### Question 34
using LinearAlgebra, Statistics
include("GradientDescent.jl")
include("genEigenvalueMat.jl")
"""
Function to compare convergence of Conjugant Gradient under different matrices
"""
function CGComparison(; maxiter :: Integer = 50, ϵ :: Float64 = 1e-14, dim :: Integer = 100)
    #Generate matrix to number of iterations
    for i in (1.01, 1.25, 1.5, 2, 4, 8, 15, 25, 45, 75)
        #Generate a linear system
        A = genEigenvalueMat([1,i])
        x_sol = randn(2)
        b = A*x_sol

        #Solve the system using gradient descent
        x, t, r = GradientDescent(A,b,method = 1)
        printstyled("Error for $i : 1 Eigenvalue Ratio \n")
        printstyled("$r\n\n", color=:red)
    end
end



