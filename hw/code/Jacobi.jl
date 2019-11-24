### Question 31
using LinearAlgebra, Statistics

"""
- `b :: Vector in R^n`
- `A :: A matrix in R^(nxn)`
- `x :: Initial guess at x`
- `iter :: Iterations, an integer`
"""
function Jacobi(A :: Matrix, b :: Vector; x = 'f', iter :: Integer = 10, ϵ :: Float64 = 1e-14)
    #Error checking
    n, m = size(A)
    n == length(b) || @error "Dimension mismatch between coefficient matrix and constant vector"
    m == n || @error "A must be a square matrix"
    sum(broadcast(abs,diag(A)).>0) == n || @error "All diagonal entries must be non-zero"

    #Start with vector of zeros as x
    if x == 'f'
        x = zeros(n)
    elseif length(x) != n
        @error "Initial guess at x must be same size as b"
    end
    t = 0
    while t < iter
        for i = 1:n
            sigma = 0
            for j = 1:n
                if j != i
                    sigma += A[i,j]*x[j]
                end
            end
            x[i] = (b[i] - sigma)/A[i,i]
        end
        t += 1
    end
    return (x)
end


function JacobiExamples()
    dimens = [100, 500, 1000]

    for dim in dimens
        printstyled("Generate a $dim x $dim system.\n",color=:yellow)
        T = randn(dim,dim)
        A = T + I.*maximum(sum(abs.(T),dims=2))
        x_sol = randn(dim)
        b = A*x_sol

        printstyled("Solving the linear system using ")
        printstyled("Jacobi Method\n", color=:blue)
        x = Jacobi(A, b)

        printstyled("Error: ")
        printstyled("$(norm(x - x_sol))\n\n", color=:red)
    end
end

