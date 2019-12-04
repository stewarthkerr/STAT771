### Added Question 22
#https://en.wikipedia.org/wiki/Arnoldi_iteration
using LinearAlgebra, Statistics, Random

"""
Uses the Arnoldi's Method to generate an orthonormal basis of the Krylov space of A wrt vector v
- `v :: Vector in R^n`
- `A :: A square matrix in R^(nxn) 
"""
function Arnoldi(A :: Matrix, v :: Vector; ϵ :: Float64 = 1e-14)
    #Error checking
    m, n = size(A)
    m == n || @error "A is not a square matrix"
    m == length(v) || @error "Dimension mismatch between A and v"

    h = zeros(n,n)
    Q = fill(Float64[], n)

    #Initialize first vector
    q = v/norm(v)
    Q[1] = q

    for i in 2:n
        qn = A*q
        for j = 1:(i-1)
            h[j,i-1] = Q[j]'*qn
            qn = qn - h[j,i-1]*Q[j]
        end
        h[i, i-1] = norm(qn)
        if h[i, i-1] > ϵ
            q = qn / h[i, i-1]
            Q[i] = q
        #If the norm of the new vector is zero, then we've hit the grade
        else 
            grade = i
            i = n+1
        end
    end
        
    return Q, h
end

function ArnoldiExamples()
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
