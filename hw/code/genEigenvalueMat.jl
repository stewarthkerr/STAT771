### Question 35
using LinearAlgebra, Statistics

#=
"""
Generates a dense, symmetric, square 2x2 matrix with user-specified eigenvalues
- `eigenvalues :: List of eigenvalues`
"""
function genEigenvalueMat(eigenvalues)
    #Error checking
    n = length(eigenvalues)
    0 in eigenvalues &&                   @error "Eigenvalues must be nonzero"
    length(unique(eigenvalues)) == n ||   @error "Eigenvalues must be unique"

    #Generate dense orthgonal matrix of size nxn
    x = 1/sqrt(n)
    for i = 1:n
        q = ones(n)*x
        if i >= 2
            q[i] = q[i] * -1
            Q = hcat(Q,q)
        else
            Q = q
        end
    end

    Q'*Q == I || @error "Creating orthogonal matrix failed"

    #Generate the dense symmetric square matrix using user specified eigenvalues
    A = Q'*diagonal(eigenvalues)*Q

    return A
end
=#

"""
Generates a dense, symmetric, square 2x2 matrix with user-specified eigenvalues
- `eigenvalues :: List of eigenvalues`
"""
function genEigenvalueMat(eigenvalues)
    #Error checking
    n = length(eigenvalues)
    @assert(n == 2,"Only supports 2x2 matrices at this time")
    @assert(!(0 in eigenvalues),"Eigenvalues must be nonzero")
    @assert(length(unique(eigenvalues)) == n,"Eigenvalues must be unique")

    #Generate dense orthgonal matrix of size nxn
    x = 1/sqrt(n)
    Q = ones(n,n)*x
    for i = 2:n
        Q[i,i] = x*-1
    end

    #Q'*Q == I || @error "Creating orthogonal matrix failed"

    #Generate the dense symmetric square matrix using user specified eigenvalues
    A = Q'*Diagonal(eigenvalues)*Q

    return A
end


