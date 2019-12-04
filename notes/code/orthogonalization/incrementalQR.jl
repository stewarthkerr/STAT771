using LinearAlgebra

#80 Gigs of Data
m = 10
n = 5*10^7
x_sol = randn(m)

#Data streaming mechanism
function data()
    a = randn(m)
    b = dot(a,x_sol) + randn()
    return hcat(a',b)
end

function givRotSmall(a,i,j)
    n = length(a)
    i > n && @error "Index i is out of bounds"
    j > n && @error "Index j is out of bounds"
    a[i] == 0 && @error "Entry a[i] is zero"

    len = sqrt( a[i]^2 + a[j]^2)
    λ = a[i]/len
    σ = a[j]/len

    G = [λ σ; -σ λ]
    #G[i,i] = λ
    #G[i,j] = σ
    #G[j,j] = λ
    #G[j,i] = -σ

    return G
end

#Incremental QR
function incrementalQR(data,m :: Int64,n :: Int64)
    RC = zeros(m+1,m+1)
    RC[1,:] = data()
    for i = 2:n
        ind = min(i,m+1)
        RC[ind,:] = data() #overwrites ρ since ind = m+1
        for j = 1:ind
            G = givRotSmall(RC[:,j],j,ind)
            if ind == m + 1
                RC[j,:] = G[1,1]*RC[j,:] + G[1,2]*RC[ind,:]
            else
                RC[[j,ind],:] = G*RC[[j,ind],:]
            end
        end
    end
    return RC[1:m,1:m], RC[1:m,m+1]
end


printstyled("Gigs of 64-bit Data to Process: $(n*m*8/1e9) Gigs\n",color=:yellow)
printstyled("Computing least squares solutions:\n\n")
printstyled("n\t\tTime\t\tError\n",color=:blue)
printstyled("-\t\t----\t\t-----\n",color=:blue)
for j = 1:8
    ti_el = @elapsed sol = incrementalQR(data, m, 10^j)
    err = norm(sol[1]\sol[2] - x_sol)^2
    printstyled("$j\t\t$(round(ti_el,digits=4))\t\t$(round(err,sigdigits=5))\n")
end
println()
