### Added Question 2e
using LinearAlgebra, Statistics

#Define QR object 
struct QRPI
    Q
    R
    PI
end

"""
Computes the full QR decomposition of an arbitrary matrix using Modified Gram-Schmidt and Pivoting
- `A :: An arbitrary matrix in R^(nxm).
"""
function fullQR_mgs(A :: Matrix)
    n, m = size(A)

    #Initialize Q to A and pi to zeros
    Q = deepcopy(A)
    R = zeros(m,m)
    PI = Matrix{Float64}(I, m, m)
    colnorm = zeros(m)

    for k = 1:m
        #Step 0 - Column Pivoting 
        col_norm = working_norm = working_col = 0
        for r = k:m
            col_norm = norm(Q[:,r])
            if col_norm > working_norm
                working_col = r
                working_norm = col_norm
            end
        end
        if working_norm == 0
            @error "A is not full column rank, it's only rank $k"
        end
        #Pivot in Q
        switch_col = Q[:,k]
        Q[:,k] = Q[:,working_col] 
        Q[:,working_col] = switch_col
        #Pivot in PI
        switch_pi = PI[:,k]
        PI[:,k] = PI[:,working_col]
        PI[:,working_col] = switch_pi
        #Pivot in R
        for c = 1:k-1
            switch_entry = R[c,k]
            R[c,k] = R[c,working_col]
            R[c,working_col] = switch_entry
        end

        #Step 1 - Reorthogonalization of ak wrt q1...qk-1
        for i = 1:k-1
            alpha_i = Q[:,i]'*Q[:,k]
            R[i,k] = R[i,k] + alpha_i
            Q[:,k] = Q[:,k] - (alpha_i*Q[:,i])
        end

        #Step 2 - Normalization
        R[k,k] = norm(Q[:,k])
        Q[:,k] = Q[:,k] / R[k,k]

        #Step 3 - Orthogonalization of ak wrt qk
        for j = k+1:m
            R[k,j] = Q[:,k]'*Q[:,j]
            Q[:,j] = Q[:,j] - R[k,j]*Q[:,k]
        end
    end

    return QRPI(Q,R,PI)
end


function fullQRExamples()
    dims = [(100,20,10), (20,100,10), (500,500,250)]
    ϵ = 1e-14
    for dim in dims
        A = randn(dim[1],dim[3])*randn(dim[3],dim[2])
        printstyled("\nGenerate a $((dim[1],dim[2]))-dimensional matrix with rank $(rank(A)).\n", color=:yellow)

        printstyled("Compute Full QR factorization using", color=:yellow)
        printstyled(" Modified Gram-Schmidt with Column Pivoting\n", color=:cyan)
        y = fullQR_mgs(A)

        printstyled("Recover A from Q,R \n", color=:yellow)
        AR = y.Q*y.R*y.PI'

        printstyled("Compare original and recovered A\n", color=:yellow)
        error = sum((A-AR).^2)
        printstyled(" sum((A - AR) .<= 1e-14) = $(error .<= ϵ)\n", color=:cyan)

        printstyled("Total error from QR factorization\n", color=:yellow)
        printstyled("$(error)", color=:red)
    end
    nothing
end
