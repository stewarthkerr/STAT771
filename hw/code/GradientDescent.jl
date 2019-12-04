### Question 34
using LinearAlgebra, Statistics

"""
- `b :: Vector in R^n`
- `A :: A matrix in R^(nxn) (we assume invertible)
- `x :: Initial guess at x`
- `method :: one of the 4 strategies for choosing the step size alpha (see notes - methods 1 = 2, and method 3 = 4)
- `iter :: Iterations, an integer`
"""
function GradientDescent(A :: Matrix, b :: Vector; x = 'f', method :: Int = 2, iter :: Integer = 500, ϵ :: Float64 = 1e-14)
    #Error checking
    n, m = size(A)
    n == length(b) || @error "Dimension mismatch between coefficient matrix and constant vector"
    m == n || @error "A must be an invertible square matrix"
    sum(broadcast(abs,diag(A)).>0) == n || @error "All diagonal entries must be non-zero"
    method in (1,2,3,4) || @error "Method must be one of (1,2,3,4) for choosing step size"

    #Start with vector of zeros as x
    if x == 'f'
        x = zeros(n)
    elseif length(x) != n
        @error "Initial guess at x must be same size as b"
    end

    #Calculate initial error
    r = A*x-b  

    #If method is three or four, alpha does not depend on iteration
    if method in (3,4)
        alpha = 2/(svd(A).S[1]^2 + svd(A).S[n]^2)
    elseif method == 1 #Calculate initial alpha, must update every iteration
        alpha = (norm(A'r)^2)/(norm(A*A'*r)^2)
    elseif method == 2
        alpha = (norm(r)^2)/(norm(r'*A*r)^2)
    end

    conv_iter = t = 0
    while t < iter 
        x = x - alpha*(A'*r)
        r = A*x-b
        if method == 1
            alpha = (norm(A'r)^2)/(norm(A*A'*r)^2)
        elseif method == 2
            alpha = (norm(r)^2)/(norm(r'*A*r)^2)
        end

        #If the error is less than precision, we are done
        if norm(r) < ϵ
            conv_iter = t
            t = iter+1
        else
            conv_iter = t
        end

        t += 1
    end

    return x, conv_iter, norm(r)
end


function GradientDescentExamples()
    dimens = [100, 500]

    for dim in dimens
        printstyled("Generate a $dim x $dim system.\n",color=:yellow)
        T = randn(dim,dim)
        A = T + I.*maximum(sum(abs.(T),dims=2))
        A = (A'+A)
        x_sol = randn(dim)
        b = A*x_sol

        printstyled("Solving the linear system using ")
        printstyled("Gradient Descent Strategies 1,2, & (3/4)\n", color=:blue)
        x1, t1 = GradientDescent(A,b,method = 1)
        x2, t2 = GradientDescent(A,b,method = 2)
        x3, t3 = GradientDescent(A,b,method = 3)

        printstyled("Strategy 1 Error: ")
        printstyled("$(norm(x1 - x_sol))\n\n", color=:red)

        printstyled("Strategy 2 Error: ")
        printstyled("$(norm(x2 - x_sol))\n\n", color=:red)

        printstyled("Strategy (3/4) Error: ")
        printstyled("$(norm(x3 - x_sol))\n\n", color=:red)

    end
    #=
    printstyled("Generate a 50 x 50 system.\n",color=:yellow)
    A = randn(50,50)
    x_sol = randn(50)
    b = A*x_sol

    printstyled("Solving the linear system using ")
    printstyled("Gradient Descent Strategies (1/2) & (3/4)\n", color=:blue)
    x1, t1, e1 = GradientDescent(A,b,method = 1, iter = 10000)
    #x2, t2 = GradientDescent(A,b,method = 2, iter = Inf)
    x3, t3 = GradientDescent(A,b,method = 3, iter = 10000)

    printstyled("Strategy 1 Iterations to Convergence: ")
    printstyled("$t1\n\n", color=:red)

    #printstyled("Strategy 2 Error: ")
    #printstyled("$t2\n\n", color=:red)

    printstyled("Strategy (3/4) Iterations to Convergence: ")
    printstyled("$t3\n\n", color=:red)

    #printstyled("Strategy 1 Final Residual Norm: ")
    #printstyled("$e1\n\n", color=:red)
    =#
end
