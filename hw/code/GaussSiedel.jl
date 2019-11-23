### Added Question 13
using LinearAlgebra, Statistics

"""
- `b :: Vector in R^n`
- `A :: A matrix in R^(nxn)`
- `method :: f for forward, b for backward`
- `iter :: Iterations, an integer`
"""
function GaussSeidel(A :: Matrix, b :: Vector; method :: Char = 'b', iter :: Integer = 10, ϵ :: Float64 = 1e-14)
    #Error checking
    n, m = size(A)
    n == length(b) || @error "Dimension mismatch between coefficient matrix and constant vector"
    m == n || @error "A must be a square matrix"
    sum(broadcast(abs,diag(A)).>0) == n || @error "All diagonal entries must be non-zero"
    uppercase(method) in ('F','B') || @error "Method must be either f for forward or b for backward"

    #If forwards, start at 1 if backwards, start at n
    if uppercase(method) == 'F'
        start = 1; finish = n; steps = 1;
    else
        start = n; finish = 1; steps = -1;
    end

    #Start with vector of zeros as x
    x = zeros(n)
    t = 0
    while t <= iter
        for i in start:steps:finish
            x[i] = x[i] + (b[i] - (A*x)[i])/A[i,i]
        end
        t += 1
    end
    return (x)
end


function backGSExamples()
    dimens = [100, 500, 1000]

    for dim in dimens
        printstyled("Generate a $dim x $dim system.\n",color=:yellow)
        T = randn(dim,dim)
        A = T + I.*maximum(sum(abs.(T),dims=2))
        x_sol = randn(dim)
        b = A*x_sol

        printstyled("Solving the linear system using ")
        printstyled("Backwards Gauss-Seidel\n", color=:blue)
        x = GaussSeidel(A, b)

        printstyled("Error: ")
        printstyled("$(norm(x - x_sol))\n\n", color=:red)
    end
end

function forwardGSExamples()
    dimens = [100, 500, 1000]

    for dim in dimens
        printstyled("Generate a $dim x $dim system.\n",color=:yellow)
        T = randn(dim,dim)
        A = T + I.*maximum(sum(abs.(T),dims=2))
        x_sol = randn(dim)
        b = A*x_sol

        printstyled("Solving the linear system using ")
        printstyled("Forwards Gauss-Seidel\n", color=:blue)
        x = GaussSeidel(A, b, method='f')

        printstyled("Error: ")
        printstyled("$(norm(x - x_sol))\n\n", color=:red)
    end
end

