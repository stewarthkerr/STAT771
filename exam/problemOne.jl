#Do not change the name of this file.

function solver(arg1, arg2) #Do not change the name of this function
A = arg1
b = arg2

x = fullRankLR(A,b)

return x

end

#Feel free to add any helper functions below this line,
#but do not call any of them ``solver''

using LinearAlgebra

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

"""
- `b :: Vector`
- `R :: Upper triangular matrix`
Determines if R is invertible. If R is not invertible, throws error.
    IF R is invertible, then solves Rx=b.
"""
function upperTriSolve(R,b)
    # Check for R not invertible
    n,m = size(R)
    n != m && @error "R is not a square matrix. Not invertible."
    @assert( minimum(abs.(diag(R))) > 1e-8, "R contains an (almost) 0 along the diagonal. Not invertible.")

    # Check to make sure R,b are same size
    size(R,1) != size(b,1) && @error "R and b are not same size. Can't solve this system."

    #Implementation of back-substition
    x = zeros(m)
    for i = m:-1:1
        od = 0.0
        for j = (i+1):m
            od += (R[i,j]*x[j])
        end
        x[i] = (b[i] - od)/R[i,i]
    end

    return x
end

"""
- `b :: Vector in R^n`
- `A :: A matrix in R^(nxm) with n >= m, and full column rank`
Determines if A is full rank. If not, throws an error.
    If A is full rank, then solves min(Ax-b)2norm.
"""
function fullRankLR(A :: Matrix, b :: Vector; ϵ :: Float64 = 1e-14)
    #Error checking
    n, m = size(A)
    n == length(b) || @error "Dimension mismatch between coefficient matrix and constant vector"
    m <= n || @error "Dimensions of A are off, m > n"

    #Get QR Decomposition with pivoting
    F = mgs_cp(A)

    #Check R against ϵ
    (sum(abs.(diag(F.R)) .< ϵ) > 0) && @error "A has effective rank less than $m"

    #Second, solve for piTx using back substitution
    c = (F.Q'*b)[1:m] 
    piTx = upperTriSolve(F.R,c)

    #Have to multiply by F.P to return original x (not permuted x)
    x = F.PI*piTx
    return x
end

function problemOneExamples()
    dims = [(100,20), (500,100), (1000,500)]
    for dim in dims
        printstyled("\nGenerate a $dim-dimensional linear system.\n", color=:yellow)
        A = randn(dim)
        x_sol = randn(dim[2])
        b = A*x_sol+randn(dim[1])

        printstyled("Solving linear system with ")
        printstyled("fullRankLR\n", color=:cyan)
        y = solver(A, b)

        printstyled("Error:")
        printstyled("$(norm(x_sol-y))\n", color=:red)
    end
    nothing
end

