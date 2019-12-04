using LinearAlgebra

"""
Given arguments `A :: Matrix` and `b :: Vector`,
computes the minimum 2-norm solution that solves
|Ay - b|₂
"""
function min2LinReg(A :: Matrix, b :: Vector)
    n, m = size(A)
    n == length(b) || @error "Dimension mismatch between coefficient matrix and constant vector"

    #Get (thin) SVD of A
    U, S, V = svd(A)

    #Get c
    c = U'*b

    #Get z, note S is a vector
    z = c./S

    #Transform solution back
    return V*z
end

function p1Examples()
    dims = [(100,20), (250,250), (500,1000)]
    for dim in dims
        printstyled("\nGenerate a $dim-dimensional system.\n", color=:yellow)
        A = randn(dim)
        b = randn(dim[1])

        printstyled("Solving system with ")
        printstyled("min2LinReg\n", color=:cyan)
        y = min2LinReg(A, b)

        printstyled("Error: ")
        printstyled("$(norm(A\b-y))\n", color=:red)
    end
    nothing
end
