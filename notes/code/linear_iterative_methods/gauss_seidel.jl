using LinearAlgebra

"""
Implements Gauss Seidel method
"""
function gaussSeidel(A,b, x; iter = 100)

    #Check diagonals of A
    D = Diagonal(A)
    minimum(abs.(diag(D))) < 1e-14 && @error "Diagonal elements are nearly zero."

    #Get D+E
    L = LowerTriangular(A)

    #get F
    F = UpperTriangular(A) - D

    #Compute G and f
    f = L\b
    G = - (L \ F)

    #Iterate
    for i = 1:iter
        x = G*x + f
    end

    return x
end

function gsExamples()
    dimens = [100, 500, 1000]

    for dim in dimens
        printstyled("Generate a $dim x $dim system.\n",color=:yellow)
        T = randn(dim,dim)
        A = T + I.*maximum(sum(abs.(T),dims=2))
        x_sol = randn(dim)
        b = A*x_sol

        printstyled("Solving the linear system using ")
        printstyled("gaussSeidel\n", color=:blue)
        x = gaussSeidel(A, b, zeros(dim))

        printstyled("Error: ")
        printstyled("$(norm(x - x_sol))\n\n", color=:red)
    end
end
