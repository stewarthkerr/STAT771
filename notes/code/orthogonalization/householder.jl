using LinearAlgebra

"""
Implementation of Householder Reflections that takes `a :: Vector` and
- returns the `v :: Vector` that produces a rotation such that
- the reflected `a` is a scaling of the first standard basis vector
"""
function householder(a :: Vector; ϵ :: Float64 = 1e-14)
    #Compute the norm
    α = norm(a)
    α < 1e-14 && @error "Vector is too close to zero."

    #Compute v
    v = deepcopy(a)
    v[1] -= α
    v = v/norm(v)

    #return results
    return v, α
end

function hrExample()
    dims = [10, 100, 1000]
    for dim in dims
        printstyled("\nHouseholder reflection for random $dim-vector\n",color=:yellow)
        a = randn(dim)
        v, α = householder(a)
        refl = a - v*(2*dot(v,a))
        printstyled("Verifying that first entry is norm of a: ")
        printstyled("$(abs(refl[1] - α) < 1e-14)\n", color=:red)
        printstyled("Verifying that other entries are zero: ")
        printstyled("$(norm(refl[2:end]) < 1e-14)\n\n",color=:red)
    end
end
