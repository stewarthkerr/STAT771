using LinearAlgebra

"""
Given arguments `A :: Matrix` and `b :: Vector`, either
- Solves the linear system if it has a unique solution
- Returns an error describing why the system cannot be solved
"""
function linSys(A :: Matrix, b :: Vector; ϵ :: Float64 = 1e-14)
    n, m = size(A)
    n == length(b) || @error "Dimension mismatch between coefficient matrix and constant vector"
    m <= n || @error "Underdetermined system"

    #Get QR Decomposition with pivoting
    F = qr(A, Val(true))

    #Check R against ϵ
    (sum(abs.(diag(F.R)) .< ϵ) > 0) && @error "A has effective rank less than $m"

    #Otherwise solve with (your) triangular system solver
    c = (F.Q'*b)[1:m]
    x = F.P*(F.R[1:m,1:m] \ c)

    return x
end


function p1Examples()
    dims = [(100,20), (500,100), (1000,500)]
    for dim in dims
        printstyled("\nGenerate a $dim-dimensional linear system.\n", color=:yellow)
        A = randn(dim)
        x_sol = randn(dim[2])
        b = A*x_sol

        printstyled("Solving linear system with ")
        printstyled("linSys\n", color=:cyan)
        y = linSys(A, b)

        printstyled("Error: ")
        printstyled("$(norm(x_sol-y))\n", color=:red)
    end
    nothing
end
