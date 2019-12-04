using LinearAlgebra, Random, SparseArrays, UnicodePlots

#Makes Random Clustered Graphs
function makeClusteredGraph(N, c_ind)
    reorder = randperm(N)
    clusters = map( λ -> reorder[λ], c_ind)

    #Generate Sparse Matrix
    i = Int64[]
    j = Int64[]
    for cluster in clusters
        L = length(cluster)
        i_loc = zeros(Int64,L)
        j_loc = zeros(Int64,L)
        for t = 1:L-1
            i_loc[t] = cluster[t]
            j_loc[t] = cluster[t+1]
        end
        i_loc[L] = cluster[L]
        j_loc[L] = cluster[1]
        append!(i, vcat(i_loc,j_loc))
        append!(j, vcat(j_loc,i_loc))
    end

    #Create Spares Array
    L = length(i)
    graph = sparse( i, j, ones(L))
    return graph, clusters
end

########################
#Case 1 N = 10
########################

#Number of Nodes
N = 10
c_ind = [ 1:3, 4:7, 8:10]

#Create Random Graph
graph, partition = makeClusteredGraph(N, c_ind)

#Visualize adjacency matrix
println(spy(graph))

#Create Dense Graph (not necessary)
graph = Matrix(graph);

#Get Eigenvalue Decomposition
Λ, P = eigen(graph);

#Take last columsn of P, note P'1 should be nearly zero for unwanted columns
ind = findall( abs.(P'*ones(N)) .>= 1e-16*N)[end-length(c_ind)+1:end]

#Find Partitions
est_partition = map(λ -> findall(abs.(P[:,λ]) .> 1e-16*N), ind)

sort.(partition)
est_partition

#########################
#Case 2: N = 100
#########################

#Number of Nodes
N = 100
c_ind = [ 1:15, 16:32, 33:58, 59:100]

#Create Random Graph
graph, partition = makeClusteredGraph(N, c_ind)

#Visualize adjacency matrix
println(spy(graph))

#Create Dense Graph (not necessary)
graph = Matrix(graph);

#Get Eigenvalue Decomposition
Λ, P = eigen(graph);

#Take last columsn of P, note P'1 should be nearly zero for unwanted columns
ind = findall( abs.(P'*ones(N)) .>= 1e-16*N)[end-length(c_ind)+1:end]

#Find Partitions
est_partition = map(λ -> findall(abs.(P[:,λ]) .> 1e-16*N), ind)

sort.(partition)
est_partition
