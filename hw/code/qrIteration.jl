### Added Question 10
using LinearAlgebra, Statistics, SparseArrays
include("qrHouseholder.jl")

"""
Uses the QR iteration to alue of a diagonalizable matrix
- `A :: A full rank column matrix in R^(nxn).
- `iter :: Number of iterations`
"""
function qrIteration(A :: Matrix; iter :: Integer = 500)
    #Get size of the matrix A
    n, m = size(A)
    m == n || @error "Matrix A is not square"

    QR = qrHouseholder(A)
    gamma = QR.Q'*A*QR.Q

    for i=1:iter
        QR = qrHouseholder(gamma)
        gamma = QR.R*QR.Q
    end

    return gamma
end


function qrIterationExamples()
    dims = [10, 25, 50]
    ϵ = 1e-14
    for dim in dims
        printstyled("\nGenerate a $dim-dimensional square matrix.\n", color=:yellow)
        A = randn(dim,dim)

        printstyled("Estimate eigenvalues", color=:yellow)
        printstyled(" QR Iteration\n", color=:cyan)
        y = sort(diag(qrIteration(A)), rev=true)

        printstyled("Compare to Julia's internal eigenvalue solver (eigvals)\n", color=:yellow)
        y_sol = eigvals(A)

        printstyled("Sum of Eigenvalues from qrIteration: ", color=:cyan)
        printstyled("$(sum(y))\n")
        printstyled("Sum of Eigenvalues from eigvals: ", color=:cyan)
        printstyled("$(sum(y_sol))\n")
        #printstyled("Sum of Aboslute Difference between methods: ", color=:cyan)
        #printstyled("$(sum(abs.(y-y_sol)))")

    end
    nothing
end
