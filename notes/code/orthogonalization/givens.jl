using LinearAlgebra, SparseArrays

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

function grExample()
    dims = [10,100,1000]
    for dim in dims
        printstyled("\nGivens rotation for random $dim-vector\n",color=:yellow)
        a = randn(dim)
        b = deepcopy(a)
        for i = 2:dim
            G = givRot(b, 1, i)
            b = G*b
        end
        printstyled("Verifying that first entry is norm of a: ")
        printstyled("$(abs(b[1] - norm(a)) < 1e-12)\n", color=:red)
        printstyled("Verifying that other entries are zero: ")
        printstyled("$(norm(b[2:end]) < 1e-14)\n\n",color=:red)
    end
end
