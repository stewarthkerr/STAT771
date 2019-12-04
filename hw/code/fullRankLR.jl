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
function fullRankLR(A :: Matrix, b :: Vector; ϵ :: Float64 = 1e-14)
    #Error checking
    n, m = size(A)
    n == length(b) || @error "Dimension mismatch between coefficient matrix and constant vector"
    m <= n || @error "Dimensions of A are off, m > n"

    #Get QR Decomposition with pivoting
    F = qr(A, Val(true))

    #Check R against ϵ
    (sum(abs.(diag(F.R)) .< ϵ) > 0) && @error "A has effective rank less than $m"

    #Otherwise, solve the system first using inverse
    c = (F.Q'*b)[1:m] 
    x = F.P*(inv(F.R[1:m,1:m]) * c)

    #Second, solve for piTx using back substitution
    piTx = upperTriSolve(F.R,c)

    #Have to multiply by F.P to return original x (not permuted x)
    return (x, F.P*piTx)
end


function fullRankLRExamples()
    dims = [(100,20), (500,100), (1000,500)]
    for dim in dims
        printstyled("\nGenerate a $dim-dimensional linear system.\n", color=:yellow)
        A = randn(dim)
        x_sol = randn(dim[2])
        b = A*x_sol+randn(dim[1])

        printstyled("Solving linear system with ")
        printstyled("fullRankLR\n", color=:cyan)
        y = fullRankLR(A, b)

        printstyled("Julia Inverse Error: ")
        printstyled("$(norm(x_sol-y[1]))\n", color=:red)

        printstyled("upperTriSolve Backsub Error: ")
        printstyled("$(norm(x_sol-y[2]))\n", color=:red)
    end
    nothing
end
