### Added Question 10
using LinearAlgebra, Statistics, SparseArrays

"""
Uses the power Method to return the largest eigen value of a diagonalizable matrix
- `A :: A full rank column matrix in R^(nxn).
- `iter :: Number of iterations`
"""
function powerMethod(A :: Matrix; iter :: Integer = 500)
    #Get size of the matrix A
    n, m = size(A)
    m == n || @error "Matrix A is not square"

    p_k = randn(n)

    for i=1:iter
        #Estimate the primary eigenvector
        p_k = (A*p_k)/norm(A*p_k)
    end

    #Calculate estimate of largest eigenvalue
    lambda = p_k'*A*p_k

    return (lambda, p_k)
end


function powerMethodExamples()
    dims = [25, 100, 500]
    ϵ = 1e-14
    for dim in dims
        printstyled("\nGenerate a $dim-dimensional square matrix.\n", color=:yellow)
        A = randn(dim,dim)

        printstyled("Estimate largest eigen value/vector using", color=:yellow)
        printstyled(" Power Method\n", color=:cyan)
        y = powerMethod(A)

        printstyled("Compare to Julia's internal eigenvalue solver")
    end
    nothing
end
