using LinearAlgebra, Distributions

"""
Inexact rank-k approximation to a n by m `A :: Matrix`, given
- `k :: Int64` the approximation level
- `c :: Int64` an integer between k and m
- `p :: Vector` is an m-dimensional vector non-negative sums to 1

Returns
- `H :: Matrix` which approximates the top k left singular vectors
- `S :: Vector` which approximates the top k singular values
"""
function inexactRankK(A :: Matrix, k :: Int64, c :: Int64, p :: Vector)
    #Get distribution over the columns
    dist = Categorical(p)
    cols_ind = rand(dist, c)

    #Divide sampled columns by sqrt(c*probability)
    C = A[:,cols_ind] ./ sqrt.( c*p[cols_ind]')

    #Compute Eigenvalue Decomposition of C'*C (c by c matrix)
    F = eigen(C'*C)

    #Compute top k singular values
    S = sqrt.( F.values[end-k+1:end])

    #Find corresponding top k singular vectors of C
    Z = F.vectors[:,end-k+1:end]

    #Compute left singular vector approximation
    H = (C*Z)./(S')

    return H, S
end

function inexactExamples()
    n, m, r = 30, 500, 15
    #Generate random matrix
    printstyled("\nGenerate a ($n,$m)-dimensional matrix with rank $r.\n", color=:yellow)
    U, _ = qr( randn(n,n))
    V, _ = qr( randn(m,m))
    Σ = [diagm(0 => 10*rand(r)) zeros(r,m-r);
         zeros(n-r,r) zeros(n-r,m-r)]
    A = U*Σ*V'

    #Sample from k, c, p
    p = ones(m)/m
    kc = [(5, 30), (5, 90), (10, 30), (10, 90), (15, 30), (15, 90), (20, 30), (20, 90)]

    for (k,c) in kc
        printstyled("Sampling ")
        printstyled("$c",color=:blue)
        printstyled(" columns for rank-")
        printstyled("$k",color=:blue)
        printstyled(" approximation\n")
        H, S = inexactRankK(A, k, c, p)
        printstyled("Column Space Error: ")
        printstyled("$(sum( (A - H*H'*A).^2))\n\n", color=:red)
    end
    nothing
end
