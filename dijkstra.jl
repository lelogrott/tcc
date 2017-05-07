function index_of_not_visited_lowest_value_element(array::Array{Float64}, visited::Array{Int64, 1})
    index = 0
    lowest_value = 999999
    for i=1:length(array)
        if (array[i] < lowest_value && visited[i] == 0)
            index = i
            lowest_value = array[i]
        end
    end
    return index
end
function dijkstra(Graph::Array{Float64,2}, src::Int64, dst::Int64, desired_band::Float64)
    n_vertices = size(Graph, 1)
    dist = Array{Float64, 1}()
    visited = Array{Int64, 1}()
    precedents = Array{Int64}[[] for k in 1:n_vertices]
    path = Array{Int64, 1}()
    for i=1:n_vertices
        push!(dist, 999999.0)
        push!(visited, 0.0)
    end
    dist[src] = 0
    for i=1:n_vertices
        u = index_of_not_visited_lowest_value_element(dist, visited)
        if (u==0)
            continue
        end
        visited[u] = 1
        for j=1:n_vertices
            if ((visited[j]==0) && (Graph[u, j] != 0) && (dist[u] != 999999) && (dist[u] + Graph[u,j] <= dist[j]))
                if (Graph[u, j] >= desired_band) 
                    if (dist[u] + Graph[u, j] == dist[j])
                        push!(precedents[j], u)
                    else
                        splice!(precedents[j], 1:length(precedents[j]), u)
                    end
                    dist[j] = dist[u] + Graph[u, j];
                else
                    visited[j]=1
                end
            end
        end
    end
    if (sum([sum(k) for k in precedents]) > 0) #garantindo que ha algum precedente de algum no
        while(dst != src)
            unshift!(path, dst)
            dst = precedents[dst][1]#just take the first one
        end
    end
    return path
    # println("VERTICES   : ");
    # for i = 1:n_vertices
    #     println(i, " ");
    # end
    # println("\n");
    # println("ESTIMATIVAS: ");
    # for i = 1:n_vertices
    #     println(dist[i], " ");
    # end
    # println("\n");
    # println("PRECEDENTES:\n");
    # for i = 1:n_vertices
    #     println("Precedente(s) vertice ", i, ": ");
    #     println(precedents[i])
    #     println("\n");
    # end

    # println("\nDistancias:\n");
    # for i = 1:n_vertices
    #     println("Vertice ", src," ao vertice ", i," : ", dist[i]);
    # end
end

# Graph = readdlm("BandM_fisico")
# path = dijkstra(Graph, 1, 3, 25.0)
# println(path)
