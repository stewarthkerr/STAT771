### Question 20
### https://core.ac.uk/download/pdf/82066579.pdf
using LinearAlgebra, Statistics

#Define QR object
struct QRPI
    Q
    R
    PI
end

"""
Performs Modified Gram-Schmidt with Column Pivoting to calculate QR decomposition
- `A :: An arbitrary matrix in R^(nxm)`
"""
function mgs_cp(A :: Matrix)
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


function mgscpExamples()
    dims = [(100,20), (500,100), (1000,500)]
    ϵ = 1e-14
    for dim in dims
        printstyled("\nGenerate a $dim-dimensional matrix.\n", color=:yellow)
        A = randn(dim)

        printstyled("Compute QR factorization using", color=:yellow)
        printstyled(" Modified Gram-Schmidt with Column Pivoting\n", color=:cyan)
        y = mgs_cp(A)

        printstyled("Recover A from Q,R \n", color=:yellow)
        AR = y.Q*y.R*y.PI'

        printstyled("Compare original and recovered A\n", color=:yellow)
        error = sum((A-AR).^2)
        printstyled(" sum((A - AR) .<= 1e-14) = $(error .<= ϵ)\n", color=:cyan)

        printstyled("Total error from QR factorization\n", color=:yellow)
        printstyled("$(error)", color=:red)
    end

    A = randn(100,20)
    A[:,5] = A[:,1]
    A[1,5] = A[1,5] + 1e-15
    printstyled("\nGenerate a (nearly) rank-deficient matrix.\n", color=:yellow)
    printstyled("\nNumber of columns of A = 20, Rank(A) = $(rank(A))\n", color=:yellow)


    printstyled("Compute QR factorization using", color=:yellow)
    printstyled(" Modified Gram-Schmidt with Column Pivoting\n", color=:cyan)
    y = mgs_cp(A)

    printstyled("Recover A from Q,R \n", color=:yellow)
    AR = y.Q*y.R*y.PI'

    printstyled("Compare original and recovered A\n", color=:yellow)
    error = sum((A-AR).^2)
    printstyled(" sum((A - AR) .<= 1e-14) = $(error .<= ϵ)\n", color=:cyan)

    printstyled("Total error from QR factorization\n", color=:yellow)
    printstyled("$(error)", color=:red)



    nothing
end

 #= 
    #MGS with Column Pivoting
    for i = 1:m
        #Find the highest norm column and pivot to work on that column next
        cp_norm = working_norm = working_col = 0
        for k = i:m
            cp_norm = norm(Q[:,k])
            if cp_norm > working_norm
                working_col = k
                working_norm = cp_norm
            end
        end
        #print("\n $working_col, $working_norm")
        PI[i, working_col] = 1
        switch = Q[:,i]
        Q[:,i] = Q[:,working_col] 
        Q[:,working_col] = switch
    
        #Perform orthogonalization
        for j = 1:i-1
            Q[:,i] = Q[:,i] - ((Q[:,i]'*Q[:,j])/(norm(Q[:,j]))*Q[:,j])
        end

        #normalize
        Q[:,i] = Q[:,i]/working_norm
    end

    #Note that A = Q*R --> Q^T*A = R
    R = Q'*A
=#
