#Do not change the name of this file.

function solverOne(arg1, arg2) #Do not change the name of this function
A = arg1
b = arg2

x = Jacobi(A,b, iter = 100)

return x

end

function solverTwo(arg1, arg2) #Do not change the name of this function
A = arg1
b = arg2

#Default Gauss Seidel solution method is backwards
x = GaussSeidel(A,b, iter = 100)

return x

end

#Feel free to add any helper functions below this line,
#but do not call any of them ``solverOne'' or ``solverTwo"


using LinearAlgebra

"""
- `b :: Vector in R^n`
- `A :: A matrix in R^(nxn)`
- `x :: Initial guess at x`
- `iter :: Iterations, an integer`
"""
function Jacobi(A :: Matrix, b :: Vector; x = 'f', iter :: Integer = 10, ϵ :: Float64 = 1e-14)
    #Error checking
    n, m = size(A)
    n == length(b) || @error "Dimension mismatch between coefficient matrix and constant vector"
    m == n || @error "A must be a square matrix"
    sum(broadcast(abs,diag(A)).>0) == n || @error "All diagonal entries must be non-zero"

    #Start with vector of zeros as x
    if x == 'f'
        x = zeros(n)
    elseif length(x) != n
        @error "Initial guess at x must be same size as b"
    end
    t = 0
    while t < iter
        for i = 1:n
            sigma = 0
            for j = 1:n
                if j != i
                    sigma += A[i,j]*x[j]
                end
            end
            x[i] = (b[i] - sigma)/A[i,i]
        end
        t += 1
    end
    return (x)
end

function GaussSeidel(A :: Matrix, b :: Vector; x = 'f', method :: Char = 'b', iter :: Integer = 10, ϵ :: Float64 = 1e-14)
    #Error checking
    n, m = size(A)
    n == length(b) || @error "Dimension mismatch between coefficient matrix and constant vector"
    m == n || @error "A must be a square matrix"
    sum(broadcast(abs,diag(A)).>0) == n || @error "All diagonal entries must be non-zero"
    uppercase(method) in ('F','B') || @error "Method must be either f for forward or b for backward"

    #If forwards, start at 1 else backwards, start at n
    if uppercase(method) == 'F'
        start = 1; finish = n; steps = 1;
    else
        start = n; finish = 1; steps = -1;
    end

    #If guess given, guess all zeros
    if x == 'f'
        x = zeros(n)
    elseif length(x) != n
        @error "Initial guess at x must be same size as b"
    end

    #Peform the Gauss Seidel algorithm
    t = 0
    while t < iter
        for i in start:steps:finish
            x[i] = x[i] + (b[i] - (A*x)[i])/A[i,i]
        end
        t += 1
    end
    return x
end


function solverOneExamples()
    dimens = [100, 500, 1000]

    for dim in dimens
        printstyled("Generate a $dim x $dim diagonally dominant system.\n",color=:yellow)
        T = randn(dim,dim)
        A = T + I.*maximum(sum(abs.(T),dims=2))
        x_sol = randn(dim)
        b = A*x_sol

        printstyled("Solving the linear system using ")
        printstyled("Jacobi Method\n", color=:blue)
        x = solverOne(A, b)

        printstyled("Error: ")
        printstyled("$(norm(x - x_sol))\n\n", color=:red)
    end
end


function solverTwoExamples()
    dimens = [100, 500, 1000]

    for dim in dimens
        printstyled("Generate a $dim x $dim diagonally dominant system.\n",color=:yellow)
        T = randn(dim,dim)
        A = T + I.*maximum(sum(abs.(T),dims=2))
        x_sol = randn(dim)
        b = A*x_sol

        printstyled("Solving the linear system using ")
        printstyled("Backwards Gauss-Seidel\n", color=:blue)
        x = solverTwo(A, b)

        printstyled("Error: ")
        printstyled("$(norm(x - x_sol))\n\n", color=:red)
    end
end


