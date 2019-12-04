using LinearAlgebra

"""
Implements Jacobi method
"""
function jacobi(A,b, x; iter = 100)

    #Check diagonals of A
    D = diag(A)
    minimum(abs.(D)) < 1e-14 && @error "Diagonal elements are nearly zero."

    #Compute G and f
    f = b./D
    G = -(A ./ D' - I)

    #Iterate
    for i = 1:iter
        x = G*x + f
    end

    return x
end

function jacobiExamples()
    dimens = [100, 500, 1000]

    for dim in dimens
        printstyled("Generate a $dim x $dim system.\n",color=:yellow)
        T = randn(dim,dim)
        A = T + I.*maximum(sum(abs.(T),dims=2))
        x_sol = randn(dim)
        b = A*x_sol

        printstyled("Solving the linear system using ")
        printstyled("jacobi\n", color=:blue)
        x = jacobi(A, b, zeros(dim))

        printstyled("Error: ")
        printstyled("$(norm(x - x_sol))\n\n", color=:red)
    end
end
