### Question 19
using LinearAlgebra, Statistics

#Define QR object
struct QR
    Q
    R
end

"""
- `A :: A matrix in R^(nxm), and full column rank`
Determines if A is full rank. If not, throws an error.
"""
function gs(A :: Matrix)
    #Error checking
    n, m = size(A)
    m == rank(A) || @error "A is not full column rank"

    #Compute orthonormal basis
    Q = zeros(n,m)
    Q[:,1] = A[:,1]/sqrt(A[:,1]'*A[:,1])
    for i = 2:m
        Q[:,i] = A[:,i]
        x = zeros(n,1)
        for j = 1:i-1
            x += ((Q[:,i]'*Q[:,j])/(Q[:,j]'*Q[:,j])*Q[:,j])
        end
        Q[:,i] = Q[:,i] - x
        Q[:,i] = Q[:,i]/sqrt(Q[:,i]'*Q[:,i])
    end

    #Note that A = Q*R --> Q^T*A = R
    R = Q'*A

    return QR(Q,R)
end


function gsExamples()
    dims = [(100,20), (500,100), (1000,500)]
    ϵ = 1e-14
    for dim in dims
        printstyled("\nGenerate a $dim-dimensional matrix.\n", color=:yellow)
        A = randn(dim)

        printstyled("Compute QR factorization using", color=:yellow)
        printstyled(" Gram-Schmidt\n", color=:cyan)
        y = gs(A)

        printstyled("Recover A from Q,R \n", color=:yellow)
        AR = y.Q*y.R

        printstyled("Compare original and recovered A\n", color=:yellow)
        error = sum((A-AR).^2)
        printstyled(" sum((A - AR) .<= 1e-14) = $(error .<= ϵ)\n", color=:cyan)

        printstyled("Total error from QR factorization\n", color=:yellow)
        printstyled("$(error)", color=:red)
    end
    nothing
end
