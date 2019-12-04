using LinearAlgebra, UnicodePlots

"""
Given a symmetrix `A :: Matrix` computes the eigenvalue decomposition.
"""
function qrIter( A :: Matrix; iterations :: Int64 = 100)
    issymmetric(A) || @error "A must be a symmetric matrix."

    n = size(A,1)

    #Compute Initial Gamma
    Γ = A

    #Storing Off-Diagonal Elements
    dta = zeros(Float64, iterations)

    #Begin iterating
    for i = 1:iterations
        U, S = qr(Γ)
        Γ = S*U
        dta[i] = sum(abs.(Γ)) - sum( abs.(diag(Γ)))
    end

    return Γ, dta
end

function spectrumExamples()
    dims = [10, 50]
    for dim in dims
        printstyled("\nGenerate a $dim-dimensional symmetric matrix.\n",color=:yellow)
        A = randn(dim,dim)
        A = (A'+A)

        printstyled("Computing eigenvalues with ")
        printstyled("qrIter\n", color=:cyan)
        Γ, data = qrIter(A, iterations=dim*60)

        printstyled("Plot of off-diagonal element sums: \n")
        println(scatterplot(log.(10,data)))

        printstyled("\nError between QR Iteration and Julia's Solver: ")
        printstyled("$(norm(sort(diag(Γ)) - sort(eigvals(A))))\n", color=:green)
    end
    nothing

end
