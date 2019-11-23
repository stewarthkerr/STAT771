### Question 17
# http://www.seas.ucla.edu/~vandenbe/133A/lectures/cls.pdf
using LinearAlgebra, Statistics
#Include upperTriSolve.jl from question 12
include("upperTriSolve.jl")
include("lowerTriSolve.jl")
include("mgs.jl")

"""
Solves the contrained LR problem
- `b :: Vector in R^n`
- `A :: A matrix in R^(nxm)`
- `C :: A matrix in R^(pxm)`
- `d :: A vector in R^p`
"""
function constrainedLR(A :: Matrix, b :: Vector, C :: Matrix, d :: Vector ; ϵ :: Float64 = 1e-14)
    #Error checking
    n, m = size(A)
    p, m2 = size(C)
    (m == m2) || @error "Dimensions of A do not match dimensions of C"
    (n != length(b)) && @error "A dimension does not meet b dimension"
    (p == length(d)) || @error "C dimension does not meet c dimension"
    #Assumption 2: p <= m <= n+p
    (p <= m <= (n+p)) || @error "Dimensions of A,C are incompatible"

    #vcat A and C and perform QR (with pivoting) on new matrix
    AC = vcat(A,C)
    #Assumption 1: AC has linearly independent columns
    rank(AC) == m || @error "AC does not have linearly independent columns"

    #Step 1 - Get QR Decompositions
    F = mgs(AC)
    R = F.R
    Q1 = F.Q[1:n,:]
    Q2 = F.Q[(n+1):(p+n),:]
    Q2t = Q2'

    F2 = mgs(Q2t)
    Q3 = F2.Q
    R2 = F2.R
    R2t = R2'

    #Step 2 - Solve R2'u=d by forward substitution and compute c = Q3'Q1'b-u
    u = lowerTriSolve(R2t,d)
    c = Q3'*Q1'*b-u

    #Step3 solve R2w=c by back subtituion and compute y = Q1'b-Q2'w
    w = upperTriSolve(R2,c)
    y = Q1'*b-Q2'*w

    #Step 4 - Solve Rx=y
    x = upperTriSolve(R,y)

    #https://stanford.edu/class/ee103/julia_slides/julia_least_squares_slides.pdf
    #kkt_sol = [2*A'*A C'; C zeros(p,p)] \ [2*A'*b; d]
    #x2 = kkt_sol[1:m]

    return x
end


function constrainedLRExamples()
    dims = [(300,100,50), (500,250,250), (600,600,600)]
    for dim in dims
        printstyled("\nGenerate a $((dim[1],dim[2]))-dimensional constrained linear system.\n",color=:yellow)
        A = randn(dim[1],dim[2])
        C = randn(dim[3],dim[2])
        y = randn(dim[2])
        d = C*y
        b = randn(dim[1])

        printstyled("Norm of the generating solution: ")
        printstyled("$(norm(y))\n", color=:red)

        printstyled("Solving linear system with ")
        printstyled("Constrained LR\n", color=:cyan)
        y_sol = constrainedLR(A,b,C,d)

        printstyled("Norm of computed solution: ")
        printstyled("$(norm(y_sol))\n", color=:red)

        printstyled("Residual Error: ")
        printstyled("$(norm(A*y_sol - b))\n", color=:green)

        printstyled("Error to A\\b: ")
        printstyled("$(norm(y_sol - A\b))\n\n", color=:green)
    end
    nothing
end



