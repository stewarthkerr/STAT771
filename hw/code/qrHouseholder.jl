### Question 23
### http://www.seas.ucla.edu/~vandenbe/133A/lectures/qr.pdf
using LinearAlgebra, Statistics

#Define QR object
struct QR
    Q
    R
end

"""
- `A :: A full rank column matrix in R^(nxm).
"""
function qrHouseholder(A :: Matrix)
    #Get size of the matrix A
    n, m = size(A)
    rank(A) < m && @error "Matrix A is not full column rank"

    #Initialize:
    # Orthogonalized matrix is I_nxn (because skinny QR) 
    # Upper triangular matrix begins as A
    Q = In = Matrix{Float64}(I, n, n)
    R = deepcopy(A)

    #QR via Householder
    for k = 1:m
        l = (n-k+1)

        #Get reflector vector
        a = R[k:n,k]
        Ik = Matrix{Float64}(I, l, l)
        w = a + sign(a[1])*norm(a)*Ik[1:l,1]
        v = w/norm(w)

        #Calculate R
        R[k:n,k:m] = R[k:n,k:m] - 2*v*v'*R[k:n,k:m]

        #Calculate Q
        vkn = zeros(n)
        vkn[k:n] = v
        Qk = In - 2*vkn*vkn'
        Qk = 
        Q = Q*Qk'
    end

    return QR(Q,R)
end


function qrHRExamples()
    dims = [(100,20), (500,100), (1000,500)]
    ϵ = 1e-14
    for dim in dims
        printstyled("\nGenerate a $dim-dimensional matrix.\n", color=:yellow)
        A = randn(dim)

        printstyled("Compute QR factorization using", color=:yellow)
        printstyled(" Householder Reflections\n", color=:cyan)
        y = qrHouseholder(A)

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
