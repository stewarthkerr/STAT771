### Added Question 14
using LinearAlgebra, Statistics

"""
- `b :: Vector in R^n`
- `A :: A matrix in R^(nxn)
- `x :: Initial guess at x`
- `w :: Relaxation constant, must be > 0`
- `method :: f for forward, b for backward, s for symmetric`
- `iter :: Iterations, an integer`
"""
function SOR(A :: Matrix, b :: Vector; x = 'f', w :: Float64 = 0.50, method :: Char = 'f', iter :: Integer = 50, ϵ :: Float64 = 1e-14)
    #Error checking
    n, m = size(A)
    n == length(b) || @error "Dimension mismatch between coefficient matrix and constant vector"
    m == n || @error "A must be a square matrix"
    sum(broadcast(abs,diag(A)).>0) == n || @error "All diagonal entries must be non-zero"
    0 < w < 2 || @error "Relaxation constant w must be between 0 and 2"
    uppercase(method) in ('F','B','S') || @error "Method must be either f for forward, b for backward, or s for symmetric"

    #If forwards, start at 1 if backwards, start at n
    if uppercase(method) in ('F','S')
        start = 1; finish = n; steps = 1;
    else
        start = n; finish = 1; steps = -1;
    end

    #Start with vector of zeros as x
    if x == 'f'
        x = zeros(n)
    elseif length(x) != n
        @error "Initial guess at x must be same size as b"
    end
    t = 0
    while t < iter 
        for i in start:steps:finish
            sigma = 0
            for j in 1:n
                if j != i
                    sigma += A[i,j] * x[j]
                end
            end
            x[i] = x[i] + w*(-x[i] + (b[i] - sigma)/A[i,i])
            if uppercase(method) == 'S'
                for i in finish:(-1*steps):start
                    sigma = 0
                    for j in 1:n
                        if j != i
                            sigma += A[i,j] * x[j]
                        end
                    end
                    x[i] = x[i] + w*(-x[i] + (b[i] - sigma)/A[i,i])
                end
            end        
        end
        t += 1
    end

    return (x)
end


function SORExamples()
    dimens = [100, 500]

    for dim in dimens
        printstyled("Generate a $dim x $dim system.\n",color=:yellow)
        T = randn(dim,dim)
        A = T + I.*maximum(sum(abs.(T),dims=2))
        x_sol = randn(dim)
        b = A*x_sol

        printstyled("Solving the linear system using ")
        printstyled("Forward Successive Over Relaxation\n", color=:blue)
        x = SOR(A, b)

        printstyled("Error: ")
        printstyled("$(norm(x - x_sol))\n\n", color=:red)

        printstyled("Solving the linear system using ")
        printstyled("Backward Successive Over Relaxation\n", color=:blue)
        x = SOR(A, b, method = 'b')

        printstyled("Error: ")
        printstyled("$(norm(x - x_sol))\n\n", color=:red)

        printstyled("Solving the linear system using ")
        printstyled("Symmetric Successive Over Relaxation\n", color=:blue)
        x = SOR(A, b, method = 's')

        printstyled("Error: ")
        printstyled("$(norm(x - x_sol))\n\n", color=:red)

    end
end
