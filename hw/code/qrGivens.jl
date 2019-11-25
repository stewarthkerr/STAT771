### Question 23
using LinearAlgebra, Statistics, SparseArrays

#Define QR object
struct QR
    Q
    R
end

#Define givens rotation function
"""
Returns the givens rotation matrix that when given
- `a :: Vector`
- `i` an index
- `j` an index
will eliminate a[j] using a[i] if a[i] is nonzero.
"""
function givRot(a,i,j)
    n = length(a)
    i > n && @error "Index i is out of bounds"
    j > n && @error "Index j is out of bounds"
    a[i] == 0 && @error "Entry a[i] is zero"

    len = sqrt( a[i]^2 + a[j]^2)
    λ = a[i]/len
    σ = a[j]/len

    I_sparse = sparse(I,n,n)
    G = I_sparse + sparse([i,i,j,j], [i,j,j,i], [λ-1, σ, λ-1, -σ],n,n)
    #G[i,i] = λ
    #G[i,j] = σ
    #G[j,j] = λ
    #G[j,i] = -σ

    return G
end


"""
- `A :: A full rank column matrix in R^(nxm).
"""
function qrGivens(A :: Matrix)
    #Get size of the matrix A
    n, m = size(A)
    final = min(n,m) #This gets the smaller of the dimensions
    rank(A) < m && @error "Matrix A is not full column rank"

    #Initialize:
    # Orthogonalized matrix is I_nxn (because skinny QR) 
    # Upper triangular matrix begins as A
    Q = In = Matrix{Float64}(I, n, n)
    R = deepcopy(A)

    #QR via Givens Rotations
    for j = 1:m
        for i = j+1:n
            #Calculate the Givens rotation to eliminate the next entry
            if R[i,j] >= 1e-14
                G = givRot(R[:,j],j,i)
            else
                G = In
            end

            #Calculate R
            R = G'*R

            #Calculate Q
            Q = Q*G
        end
    end

    return QR(Q,R)
end


function qrGivensExamples()
    dims = [(100,20), (250,75), (500,100)]
    ϵ = 1e-14
    for dim in dims
        printstyled("\nGenerate a $dim-dimensional matrix.\n", color=:yellow)
        A = randn(dim)

        printstyled("Compute QR factorization using", color=:yellow)
        printstyled(" Givens Rotations\n", color=:cyan)
        y = qrGivens(A)

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
