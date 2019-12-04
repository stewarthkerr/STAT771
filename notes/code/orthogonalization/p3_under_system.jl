using LinearAlgebra

"""
Given arguments `A :: Matrix` and `b :: Vector`, either
- Solves the underdetermined linear system for minimum norm solution
- Returns an error describing why the system cannot be solved
"""
function linSysUnder(A :: Matrix, b :: Vector; ϵ :: Float64 = 1e-14)
    n, m = size(A)
    (n >= m) && @error "System is not underdetermined"
    (n != length(b)) && @error "System dimensions do not align"

    #Get QR Decomposition
    F = qr(A, Val(true))

    #Get numerical rank of A
    r = sum(abs.(diag(F.R)) .> ϵ)

    #Get c and verify second compontent is zero
    d = F.Q'b
    norm(d[r+1:end]) > ϵ && @error "System is inconsistent"
    c = d[1:r]

    #Compute R^{inv}S
    P = F.R[1:r,1:r] \ F.R[1:r, r+1:m]

    #Compute R^{inv}c
    e = F.R[1:r,1:r]\c

    #Solve for z₂
    z₂ = (P'*P + I)\(P'*e)

    #Recover z₁
    z₁ = -P*z₂ + e

    #Recover the solution
    return F.P*vcat(z₁, z₂)
end


function p3Examples()
    dims = [(100,300), (500,1500), (1000,3000)]
    for dim in dims
        printstyled("\nGenerate a $dim-dimensional linear system.\n",color=:yellow)
        A = randn(dim)
        x_sol = vcat(randn(dim[1]), zeros(dim[2]-dim[1]))
        b = A*x_sol

        printstyled("Norm of the generating solution: ")
        printstyled("$(norm(x_sol))\n", color=:red)

        printstyled("Solving linear system with ")
        printstyled("linSysUnder\n", color=:cyan)
        y = linSysUnder(A,b)

        printstyled("Norm of solution: ")
        printstyled("$(norm(y))\n", color=:red)

        printstyled("Residual Error: ")
        printstyled("$(norm(A*y - b))\n", color=:green)

        printstyled("Error to A\\b: ")
        printstyled("$(norm(y - A\b))\n\n", color=:green)
    end
    nothing
end
