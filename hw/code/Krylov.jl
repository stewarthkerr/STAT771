### Added Question 22
# http://www.sam.math.ethz.ch/~mhg/pub/biksm.pdf
using LinearAlgebra, Statistics, Random
include("Arnoldi.jl")
include("qrGivens.jl")

"""
Uses the Arnoldi's Method to generate an orthonormal basis of the Krylov space of A wrt vector v
and then solves the linear system
- `A :: A square matrix in R^(nxn) 
- `x :: Vector in R^n, initial guess at solution`
- `b :: Vector in R^n`
- `iter:: Integer number of iterations to perform`
"""
function Krylov(A :: Matrix, b :: Vector; x = 'f', iter :: Integer = 10, ϵ :: Float64 = 1e-14)
    #Error checking
    m, n = size(A)
    m == n || @error "A is not a square matrix"
    m == length(b) || @error "Dimension mismatch between A and b"

    #Create container to hold results
    # results[1] = x_n
    # results[2] = r_n
    results = fill(Float64[], iter)

    #Start with vector of zeros as x
    if x == 'f'
        x = zeros(n)
    elseif length(x) != n
        @error "Initial guess at x must be same size as b"
    end

    #Calculate initial error and populate first result
    r0 = A*x-b
    results[1] = x
    
    #Get the Krylov space of A wrt r0
    krylov_space = Arnoldi(A,r0)

    #Let Krylov_Space[1] = A, and then solve the resulting system Ax=b
    #using Givens rotations and back substitution
    QR = qrGivens(krylov_space[1])

    #Calculate iterates of Krylov solver
    for i in 2:n
        x = x + krylov_space[1][i-1]

        #Solve for Ax = b
        results[i] = x
    end
        
    return results
end

function KrylovExamples()
    dimens = [25, 50, 100]

    #=
    for dim in dimens
        printstyled("Generate a $dim system.\n",color=:yellow)
        A = randn(dim,dim)
        x_sol = randn(dim[2])
        b = A*x_sol

        printstyled("Solving the linear system using ")
        printstyled("Cyclical Kaczmarz\n", color=:blue)
        x = Kaczmarz(A, b, iter = 20, method = "c")

        printstyled("Error: ")
        printstyled("$(norm(x - x_sol))\n\n", color=:red)

        printstyled("Solving the linear system using ")
        printstyled("Randomized Cyclical Kaczmarz (does not account for row norms of A)\n", color=:blue)
        x = Kaczmarz(A, b, iter = 20, method = "rc")

        printstyled("Error: ")
        printstyled("$(norm(x - x_sol))\n\n", color=:red)

        printstyled("Solving the linear system using ")
        printstyled("Randomized Kaczmarz\n", color=:blue)
        printstyled("NOTE: Randomized Kaczmarz does not currently cycle through every row in each iteration\n", color=:yellow)
        x = Kaczmarz(A, b, iter = 20, method = "r")

        printstyled("Error: ")
        printstyled("$(norm(x - x_sol))\n\n", color=:red)

    end
    =#
end
