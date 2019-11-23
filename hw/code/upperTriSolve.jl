### Question 12
using LinearAlgebra, Statistics
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
    1e-14<minimum(diag(R)) || @error "R contains an (almost) 0 along the diagonal. Not invertible."

    # Check to make sure R,b are same size
    size(R,1) != size(b,1) && @error "R and b are not same size. Can't solve this system."

    #If R is invertible, solve the system by Rx=b by multiplying on left by inverse(R)
    #x = inv(R)*b

    #Implementation of back-substition
    x = zeros(m)
    for i = m:-1:1
        od = 0.0
        for j = m:-1:(i+1)
            od += (R[i,j]*x[j])
        end
        x[i] = (b[i] - od)/R[i,i]
    end

    return x
end

"""
- `n:: Dimensions`
"""
function upperTriSolveExample(n)
    #Generate random matrix and vector
    R = UpperTriangular(rand(n,n))
    x_real = rand(n)
    b = R*x_real
    x_sol = upperTriSolve(R,b)

    if n <= 100
        printstyled("\n The real solution is \n")
        printstyled(x_real, color=:green)
        printstyled("\n The computed solution is \n")
        printstyled(x_sol, color=:yellow)
    end
    error = sum((x_real-x_sol).^2)
    printstyled("\n The error is \n")
    printstyled(error, color=:red) 
end

