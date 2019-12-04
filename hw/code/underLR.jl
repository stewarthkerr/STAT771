### Question 14
using LinearAlgebra, Statistics
#Include upperTriSolve.jl from question 12
include("upperTriSolve.jl")

"""
- `b :: Vector in R^n`
- `A :: A matrix in R^(nxm) with n >= m, and full column rank`
Determines if A is full rank. If not, throws an error.
    If A is full rank, then solves min(Ax-b)2norm.
"""
function underLR(A :: Matrix, b :: Vector; ϵ :: Float64 = 1e-14)
    #Error checking
    n, m = size(A)
    (n >= m) && @error "System is not underdetermined"
    (n != length(b)) && @error "System dimensions do not align"

    #Get QR Decomposition
    F = qr(A, Val(true))

    #Get numerical rank of A
    r = sum(abs.(diag(F.R)) .> ϵ)
    printstyled("Rank of A: ")
    printstyled("$(r)\n", color=:red)

    #Get c
    d = F.Q'b
    c1 = d[1:r]
    c2 = d[r+1:n]

    #print norm of c2 -- error
    #printstyled("Norm of the c2 error: ")
    #printstyled("$(norm(c2))\n", color=:red)

    #Need to solve Rz1 + Sz2 - c1 = 0 --> Rz1 = c1 - Sz2

    # 1. First, using inverse
    #Compute R^{inv}S
    P = F.R[1:r,1:r] \ F.R[1:r, r+1:m]
    #Compute R^{inv}c1
    e = F.R[1:r,1:r]\c1
    #Solve for z₂
    z₂ = (P'*P + I)\(P'*e)
    #Recover z₁
    z₁ = -P*z₂ + e
    #Recover the solution
    return F.P*vcat(z₁, z₂)

    #Second, solve for piTx using back substitution
    #piTx = upperTriSolve(F.R,c)

    #Have to multiply by F.P to return original x (not permuted x)
    #return (x, F.P*piTx)
end


function underLRExamples()
    dims = [(100,300), (500,1500), (1000,3000)]
    for dim in dims
        printstyled("\nGenerate a $dim-dimensional linear system.\n",color=:yellow)
        A = randn(dim)
        x_sol = vcat(randn(dim[1]), zeros(dim[2]-dim[1]))
        b = A*x_sol+randn(dim[1])

        printstyled("Norm of the generating solution: ")
        printstyled("$(norm(x_sol))\n", color=:red)

        printstyled("Solving linear system with ")
        printstyled("underLR\n", color=:cyan)
        y = underLR(A,b)

        printstyled("Norm of solution: ")
        printstyled("$(norm(y))\n", color=:red)

        printstyled("Residual Error: ")
        printstyled("$(norm(A*y - b))\n", color=:green)

        printstyled("Error to A\\b: ")
        printstyled("$(norm(y - A\b))\n\n", color=:green)
    end
    nothing
end
