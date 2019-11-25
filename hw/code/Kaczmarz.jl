### Question 32
using LinearAlgebra, Statistics, Random

"""
- `b :: Vector in R^m`
- `A :: A full rank matrix in R^(mxn) (n <= m)
- `x :: Initial guess at x`
- `w :: Relaxation parameter, must be > 0`
- `method :: c for cyclical, r for randomized, rc for randomized cyclical`
- `iter :: Iterations, an integer`
"""
function Kaczmarz(A :: Matrix, b :: Vector; x = 'f', w :: Float64 = 1.0, method :: String = "r", iter :: Integer = 50, ϵ :: Float64 = 1e-14)
    #Error checking
    m, n = size(A)
    m == length(b) || @error "Dimension mismatch between coefficient matrix and constant vector"
    rank(A) == n || @error "A must be full column rank"
    0 < w  || @error "Relaxation constant w must be greater than 0"
    uppercase(method) in ("C","R","RC") || @error "Method must be either c for cyclical, r for randomized, or rc for randomized cyclical"

    #If cyclical start at 1
    if uppercase(method)=="C"
        rows = collect(1:m)
    elseif uppercase(method)=="R" #Randomly choose row to cycle through according to row_norms
        row_norms = zeros(m)
        row_cdf = zeros(m)
        for i = 1:m
            row_norms[i] = norm(A[i])^2
        end
        row_probs = row_norms./sum(row_norms)
        row_cdf[1] = row_probs[1]
        for i = 2:m
            row_cdf[i]  = row_probs[i] + sum(row_cdf[i-1]) 
        end
    end

    #Start with vector of zeros as x
    if x == 'f'
        x = zeros(n)
    elseif length(x) != n
        @error "Initial guess at x must be same size as b"
    end
    t = 0
    while t < iter 
        if uppercase(method)=="RC" #Randomly permute rows to cycle through
            rows = randperm(m)
        elseif uppercase(method)=="R" #Randomly choose row to cycle through according to row_norms
            rows = collect(m+1:2*m) #Just to initiliaze rows
            runi = rand(m)
            for i = 1:m
                rows[i] = argmax(runi[i] .< row_cdf)
                #while argmax(runi[i] .< row_cdf) in rows
                #    runi[i] = rand(1)[1]
                #    rows[i] = argmax(runi[i] .< row_cdf)
                #end
            end
        end
        for i in rows
            x = x + w*(A[i,:]*(b[i] - (A[i,:]'*x))/(norm(A[i,:])^2))
        end
        t += 1
    end

    return (x)
end


function KaczmarzExamples()
    dimens = [(100,20), (500,100)]

    for dim in dimens
        printstyled("Generate a $dim system.\n",color=:yellow)
        A = randn(dim)
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
end
