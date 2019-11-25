### Question 31
using LinearAlgebra, Statistics
include("SOR.jl")
include("GaussSeidel.jl")
include("Jacobi.jl")

"""
Function to compare convergence of different iterative splitting methods (Jacobi, GaussSeidel, SOR, SSOR)
"""
function iterativeComparison(; maxiter :: Integer = 50, ϵ :: Float64 = 1e-14, dim :: Integer = 100)
    #Generate matrix to hold results
    results = fill(Float64[], 4, maxiter)
    residual_norm = zeros(4,maxiter)
    abs_error = zeros(4,maxiter)
    iterations = repeat([maxiter], 4)

    #Generate linear system
    printstyled("Generate a $dim x $dim system.\n",color=:yellow)
    T = randn(dim,dim)
    A = T + I.*maximum(sum(abs.(T),dims=2))
    x_sol = randn(dim)
    b = A*x_sol       
    x_guess = zeros(dim)
    spectral_radius = maximum(abs.(eigvals(A)))

    printstyled("Error from initial guess of zero vector\n")
    residual_norm0 = norm(x_sol-x_guess)
    abs_error0 = sum(abs.(x_sol-x_guess))
    printstyled("Residual Norm: "); printstyled("$residual_norm0\n", color=:red)
    printstyled("Aboslute Error: "); printstyled("$abs_error0\n", color=:red)

    #Perform first pass of each solver
    results[1,1] = Jacobi(A,b,iter = 1)
    results[2,1] = GaussSeidel(A,b,iter = 1, method = 'f')
    results[3,1] = SOR(A,b,iter = 1, method = 'f')
    results[4,1] = SOR(A,b,iter = 1, method = 's')

    for j in 1:4
        residual_norm[j,1] = norm(x_sol-results[j,1])
        abs_error[j,1] = maximum(abs.(x_sol-results[j,1]))
        #Converged after 1 iteration! Impressive
        if (residual_norm[j,1] < ϵ) && (abs_error[j,1] < ϵ)
            iterations[j] = 1
        end
    end

    #Iterate through solvers
    i = 2
    while i <= maxiter
        #Jacobi 
        results[1,i] = Jacobi(A,b,x = results[1,i-1], iter = 1)
        #GaussSeidel
        results[2,i] = GaussSeidel(A,b,x = results[2,i-1], iter = 1, method = 'f')
        #SOR
        results[3,i] = SOR(A,b,x = results[3,i-1], iter = 1, method = 'f')
        #SSOR
        results[4,i] = SOR(A,b,x = results[4,i-1], iter = 1, method = 's')

        #Calculate errors
        for j in 1:4
            residual_norm[j,i] = norm(x_sol-results[j,i])
            abs_error[j,i] = maximum(abs.(x_sol-results[j,i]))
            if (residual_norm[j,i-1] < ϵ) && (abs_error[j,i-1] < ϵ) && (iterations[j] > i)
                iterations[j] = i
            end
        end

        #If each method has converged, we can stop
        if iterations[1] <= i && iterations[2] <= i && iterations[3] <= i && iterations[4] <= i
            i = maxiter+1
        end
        i += 1
    end

    return (spectral_radius, results, residual_norm, abs_error, iterations)
end



#=
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
=#
