using LinearAlgebra, UnicodePlots

"""
Implements conjugate gradient method and tracks A-conjugacy and orthogonality

ARGUMENTS:
- `A :: Matrix`, symmetric positive definite matrix
- `b :: Vector`, constant vector
- `x :: Vector`, initial guess

OUTPUTS:
- `x :: Vector`, solution vector
- `conj :: Vector`, search directions
- `resid :: Vector`, residuals
"""
function conjGrad( A :: Matrix, b :: Vector, x :: Vector; track = true)
    issymmetric(A) || error("Coefficient matrix is not symmetric.")
    n = length(b)

    #Initialization (Ignore these quantities for CG question)
    conj = []
    resid = []

    #Begin Algorithm
    r = b - A*x
    p = r
    rsold = dot(r,r)

    push!(conj, p)
    push!(resid, r)

    for i = 1:n
        Ap = A*p
        alpha = rsold/ dot(p, Ap)
        x = x + alpha*p

        r = r - alpha* Ap;
        rsnew = dot(r,r)
        rsnew < 1e-30 && break

        p = r + (rsnew/rsold)*p;
        rsold = rsnew

        #Ignore this update for CG question
        push!(conj, p)
        push!(resid, r)
    end

    return x, conj, resid
end

#Experiments
using MatrixDepot
names = ["cauchy","hilb","lehmer"] #Symmetric, Positive Definite Matrices


#Generate and Solve Systems
n = 100
for entry in names
    printstyled("\nSolving: $entry\n", color=:yellow)
    A = Matrix( matrixdepot(entry,n))
    x_sol = 10*randn(n)
    b = A*x_sol

    printstyled("Initial residual norm: ")
    printstyled("$(norm(b))\n", color=:blue)

    z, conj, resid = conjGrad(A, b, zeros(n))

    printstyled("Final residual norm: ")
    printstyled("$(norm(A*z-b))\n", color=:blue)

    printstyled("Absolute Error: ")
    printstyled("$(norm(z-x_sol))\n\n", color=:blue)

    #Compute Conjugacy Errors
    normA(x) = sqrt(dot(x,A*x))
    C = hcat((conj./normA.(conj))...)
    conj_err = C'*A*C

    println(lineplot( (abs.(conj_err[1,2:end])) , title="Control of Conjugacy", ylabel="cos"))

    #Compute Orthogonalization Error
    R = hcat((resid./norm.(resid))...)
    resid_err = R'*R

    println(lineplot( (abs.(resid_err[1,2:end])), title="Control of Orthogonality",ylabel="cos"))


end
