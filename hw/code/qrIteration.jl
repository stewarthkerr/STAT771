### Added Question 10
### Look at qr_iteration.jl code
using LinearAlgebra, Statistics, SparseArrays, UnicodePlots
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
    data = zeros(Float64, iter)


    for i=1:iter
        QR = qrHouseholder(gamma)
        gamma = QR.R*QR.Q
        data[i] = sum(abs.(gamma)) - sum(abs.(diag(gamma)))

    end

    return gamma, data
end


function qrIterationExamples()
    dims = [10, 50]
    for dim in dims
        printstyled("\nGenerate a $dim-dimensional symmetric matrix.\n",color=:yellow)
        A = randn(dim,dim)
        A = (A'+A)

        printstyled("Computing eigenvalues with ")
        printstyled("QR Iteration\n", color=:cyan)
        gamma, data = qrIteration(A, iter=dim*60)

        printstyled("Plot of off-diagonal element sums: \n")
        println(scatterplot(log.(10,data)))

        printstyled("\nError between QR Iteration and Julia's Solver: ")
        printstyled("$(norm(sort(diag(gamma)) - sort(eigvals(A))))\n", color=:green)
    end
    nothing
end
