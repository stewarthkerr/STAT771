using LinearAlgebra

"""
Modified Gram-Schmidt with arbitrary inner product matrix

ARGUMENTS
- v :: Vector, to be conjugated
- A :: Matrix, defines inner product
- S :: Vector, a vector of A conjugated vectors

OUTPUT
- c :: Vector, conjugated vector
- S :: Vector, appended with v
"""
function mgs(v :: Vector, A :: Matrix,  S :: Vector)
    for s in S
        v = v - s*( v'*A*s)
    end

    nrm = v'*A*v
    c = nrm > 1e-16 ? v/sqrt( v'*A*v) : zeros(length(v))
    push!(S,c)
    return c, S
end


"""
Conjugated directions method with standard basis vectors

ARGUMENTS
- A :: Matrix, symmetric coefficient matrix
- b :: Vector, constant vector
- x :: Vector, initial estimate

OUTPUT
- x :: Vector, solution
- S :: Vector, set of A conjugated vectors
"""
function conjugateStandard(A :: Matrix, b :: Vector, x :: Vector)
    n = length(x)
    issymmetric(A) || @error "Coefficient matrix must be symmetric"

    normA(z) = sqrt( dot(z, A*z))

    E = diagm(0 => ones(n))

    r₀= A*x - b #residuals
    r = deepcopy(r₀)
    S = []                          #A conjugated vectors
    counter = 1
    while (norm(r) > 1e-14 && counter < n+1)
        s, S = mgs(E[:,counter], A, S)  #Conjugate random vector
        α = dot(r₀, s)              #Step size
        x = x - α*s                    #Update Iterate
        r = A*x - b                 #Update Residual
        counter += 1
    end
    return x, S
end

"""
Conjugated directions method with Gaussian vectors

ARGUMENTS
- A :: Matrix, symmetric coefficient matrix
- b :: Vector, constant vector
- x :: Vector, initial estimate

OUTPUT
- x :: Vector, solution
- S :: Vector, set of A conjugated vectors
"""
function conjugateGaussian(A :: Matrix, b :: Vector, x :: Vector)
    n = length(x)
    issymmetric(A) || @error "Coefficient matrix must be symmetric"

    normA(z) = sqrt( dot(z, A*z))

    r₀= A*x - b #residuals
    r = deepcopy(r₀)
    S = []                          #A conjugated vectors
    counter = 1
    while (norm(r) > 1e-14 && counter < n+1)
        s, S = mgs(randn(n), A, S)  #Conjugate random vector
        α = dot(r₀, s)              #Step size
        x = x - α*s                    #Update Iterate
        r = A*x - b                 #Update Residual
        counter += 1
    end
    return x, S
end

#Experiments
using MatrixDepot
names = ["cauchy","hilb","lehmer"] #Symmetric, Positive Definite Matrices

#View Matrices
for entry in names
    printstyled("Matrix: $entry\n", color=:yellow)
    display( matrixdepot(entry,10))
    print("\n")
end

#Generate and Solve Systems
n = 100
for entry in names
    printstyled("\nSolving (Unperturbed): $entry\n", color=:yellow)
    A = Matrix( matrixdepot(entry,n))
    x_sol = 10*randn(n)
    b = A*x_sol

    printstyled("Initial residual norm: ")
    printstyled("$(norm(b))\n", color=:blue)

    z, S = conjugateStandard(A, b, zeros(n))

    printstyled("(Standard) Final residual norm: ")
    printstyled("$(norm(A*z-b))\n", color=:blue)

    printstyled("(Standard) Absolute Error: ")
    printstyled("$(norm(z-x_sol))\n", color=:blue)

    z, S = conjugateGaussian(A, b, zeros(n))

    printstyled("(Gaussian) Final residual norm: ")
    printstyled("$(norm(A*z-b))\n", color=:blue)

    printstyled("(Gaussian) Absolute Error: ")
    printstyled("$(norm(z-x_sol))\n\n", color=:blue)

end

n = 100
for entry in names
    printstyled("\nSolving (Perturbed): $entry\n", color=:green)
    A = Matrix( matrixdepot(entry,n))
    x_sol = 10*randn(n)
    b = A*x_sol

    printstyled("Initial residual norm: ")
    printstyled("$(norm(b))\n", color=:blue)

    z, S = conjugateStandard(A + I*1e-13, b, zeros(n))

    printstyled("(Standard) Final residual norm: ")
    printstyled("$(norm(A*z-b))\n", color=:blue)

    printstyled("(Standard) Absolute Error: ")
    printstyled("$(norm(z-x_sol))\n", color=:blue)

    z, S = conjugateGaussian(A + I*1e-13, b, zeros(n))

    printstyled("(Gaussian) Final residual norm: ")
    printstyled("$(norm(A*z-b))\n", color=:blue)

    printstyled("(Gaussian) Absolute Error: ")
    printstyled("$(norm(z-x_sol))\n\n", color=:blue)
end
