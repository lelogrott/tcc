include("dijkstra.jl")

MAX_BANDWIDTH = 50
MIN_BANDWIDTH = 10
Σ = sum
Π = prod

type VirtualMachine
    id::Int64
    cp::Float64
    storage::Float64
    memory::Float64
    links::Array{Tuple{Float64,Float64}} #links são enlaces virtuais seguindo a ordem (dest,weight)
    LRC_index::Float64
end

type Resource
    top::Float64
    atual::Float64
    residual::Float64
end

type Server
    cp::Resource # adicionar parametro alfa de overbooking
    storage::Resource
    memory::Resource
    residual_band::Float64
    VMs::Array{VirtualMachine}
    custo_beneficio::Float64
end

Base.isless(x::VirtualMachine, y::VirtualMachine) = x.LRC_index > y.LRC_index
Base.isequal(x::VirtualMachine, y::VirtualMachine) = x.LRC_index == y.LRC_index

function set_lrc_index(VNR)
    for vm in VNR
        vm.LRC_index = Σ(link[2] for link in vm.links) + vm.cp
    end
end

function gera_Gs_BandM(n_nodes)
    Gs=Array{Server}(n_nodes)
    BandM = Array{Float64}(zeros(n_nodes,n_nodes))
    for i=1:n_nodes,j=1:i-1
        BandM[j,i] = BandM[i,j] = rand(0:1)==1? rand(MIN_BANDWIDTH:MAX_BANDWIDTH) : 0
        if (i==n_nodes && Σ(BandM[j:i,j])==0)#if para garantir grafo conexo
            pos = rand(j+1:n_nodes)
            BandM[pos,j] = BandM[j,pos] = rand(MIN_BANDWIDTH:MAX_BANDWIDTH)
        end
    end
    for i=1:n_nodes #inicia vetor de servidores com valores default
        Gs[i] = Server(Resource(100, rand(80:100), 0),Resource(100, rand(80:100), 0),Resource(100, rand(80:100), 0),[],0)
    end
    return Gs, BandM
end

function gera_VNR(n_nodes)
    VNR = Array{VirtualMachine}(0)
    CP_MIN = 1
    CP_MAX = 10
    STORAGE_MIN = 1
    STORAGE_MAX = 10
    MEMORY_MIN = 1
    MEMORY_MAX = 10
    BAND_MIN = 5
    BAND_MAX = 20
    for i=1:n_nodes
        links = rand(0:n_nodes)
        VM = VirtualMachine(rand(CP_MIN:CP_MAX), rand(STORAGE_MIN:STORAGE_MAX), rand(MEMORY_MIN:MEMORY_MAX), [], 0.0)
        for j=1:links
            push!(VM.links, (rand(1:n_nodes), rand(BAND_MIN:BAND_MAX)))
        end
        push!(VNR, VM)
    end
    return VNR
end

function MAUT_norm_criteria_max(array_criteria::Array{Float64})
     min = minimum(array_criteria)
     max = maximum(array_criteria)
     ↑ = (min != max) ? (x-> (x - min)/(max-min)) : (x-> 1.0)
     return copy!(array_criteria, pmap(↑, array_criteria))
end

function MAUT_norm_criteria_min(array_criteria::Array{Float64})
    min = minimum(array_criteria)
    max = maximum(array_criteria)
    ↓ = (min != max) ? (x-> 1+(min-x)/(max - min)) : (x-> 1.0)
    return copy!(array_criteria, pmap(↓, array_criteria))
end

function MAUT_marginalUtilty_custo_beneficio(array_criteria::Array{Float64})
    return copy!(array_criteria, pmap(x -> (exp(x^2)-1)/exp(1), array_criteria))
end

function MAUT_marginalUtility_cp_residual(array_criteria::Array{Float64})
    return copy!(array_criteria, pmap(x -> (exp(x^2)-1)/exp(1), array_criteria))
end

function MAUT_globalUtility(array_criterias, pesos)
    copy!(array_criterias["custo_beneficio"], array_criterias["custo_beneficio"]*pesos["custo_beneficio"])
    copy!(array_criterias["cp_residual"], array_criterias["cp_residual"]*pesos["cp_residual"])
    copy!(array_criterias["storage_residual"], array_criterias["storage_residual"]*pesos["storage_residual"])
    copy!(array_criterias["memory_residual"], array_criterias["memory_residual"]*pesos["memory_residual"])

    result = array_criterias["custo_beneficio"] + array_criterias["cp_residual"] + array_criterias["memory_residual"] + array_criterias["storage_residual"]
    return result
end

function parallel_calc_cp_residual(server::Server, vm_cp::Float64)
    server.cp.residual = server.cp.atual - vm_cp
end

function calc_cp_residual(Gs::Array{Server, 1},vm::VirtualMachine)
    return pmap(x->parallel_calc_cp_residual(x,vm.cp),Gs)
end

function parallel_calc_storage_residual(server::Server, vm_storage::Float64)
    server.storage.residual = server.storage.atual - vm_storage
end

function calc_storage_residual(Gs::Array{Server, 1},vm::VirtualMachine)
    return pmap(x->parallel_calc_storage_residual(x,vm.storage),Gs)
end

function parallel_calc_memory_residual(server::Server, vm_memory::Float64)
    server.memory.residual = server.memory.atual - vm_memory
end

function calc_memory_residual(Gs::Array{Server, 1}, vm::VirtualMachine)
    return pmap(x->parallel_calc_memory_residual(x,vm.memory),Gs)
end

function calc_custo_beneficio(Gs::Array{Server, 1}, BandM::Array{Float64, 2}, vnr::VirtualMachine)
    Band_v = Σ([link[2] for link in vnr.links])
    for i=1:length(Gs)
        if((Gs[i].cp.atual == 0) || (Gs[i].residual_band == 0))
            Gs[i].custo_beneficio = 0    
        else
            Gs[i].custo_beneficio = (vnr.cp / Gs[i].cp.atual) + (Band_v / Gs[i].residual_band)
        end
    end
end

function index_of_higher_value_element(array::Array{Float64, 1})
    index = 0
    higher_value = -1
    println("\tanalisando ", array)
    for i=1:length(array)
        if (array[i] > higher_value)
            index = i
            higher_value = array[i]
        end
    end
    println("\tindex de maior valor ", index)
    return index
end

function all_positive_residual_values_in(node::Server, VM::VirtualMachine)
    if (node.cp.residual < 0)
        return false
    end
    if (node.storage.residual < 0)
        return false
    end
    if (node.memory.residual < 0)
        return false
    end
    if (node.residual_band - Σ([link[2] for link in VM.links]) < 0)
        return false
    end
    return true
end

function do_allocate(Gs::Array{Server, 1}, selected_server::Int64, VM::VirtualMachine)
    ret = false
    for i=1:length(Gs)
        if (i == selected_server && all_positive_residual_values_in(Gs[i], VM))
            Gs[i].cp.atual = Gs[i].cp.residual
            if (Gs[i].cp.atual == 0)
                println(">>>>>>>>>> atual e residual ", Gs[i].cp.atual, Gs[i].cp.residual, selected_server)
            end
            Gs[i].storage.atual = Gs[i].storage.residual
            Gs[i].memory.atual = Gs[i].memory.residual
            push!(Gs[i].VMs, VM)
            Gs[i].residual_band -= Σ([link[2] for link in VM.links])
            ret = true
        end
        Gs[i].storage.residual = 0
        Gs[i].cp.residual = 0
        Gs[i].memory.residual = 0
        Gs[i].custo_beneficio = 0
    end
    return ret
end

function deallocate(Gs::Array{Server, 1}, BandM::Array{Float64, 2}, substrate_network_changes::Array{Tuple{Int64, Int64, Float64}, 1}, first_host_key::Int64, first_host::Int64, results::Dict{Int64,Array{Float64, 1}})
    for key in keys(results)
        # if (key == first_host_key)
        #     continue
        # end
        server = (key == first_host_key) ? first_host : index_of_higher_value_element(results[key])
        println("SERVER", server, " -> ", Gs[server])
        vm = pop!(Gs[server].VMs)
        Gs[server].cp.atual += vm.cp
        println("\t>>CP do server ", server, " -> ", Gs[server].cp.atual)
        Gs[server].storage.atual += vm.storage
        Gs[server].memory.atual += vm.memory
        Gs[server].residual_band += Σ([link[2] for link in vm.links])
        println("\tResetando server ", server)
    end
    for (src, dst, wght) in substrate_network_changes
        BandM[src, dst] += wght
    end
end

BandM = readdlm("BandM_fisico")

Gs_data_from_file = readdlm("Gs_fisico")
n_servers = size(Gs_data_from_file,1)

Gs = Array{Server}(n_servers)
for i=1:length(Gs)
    Gs[i] = Server(Resource(0, 0, 0),Resource(0, 0, 0),Resource(0, 0, 0),0,[],0)
end
for i=1:n_servers
    Gs[i].cp.top = Gs_data_from_file[i, 1]
    Gs[i].cp.atual = Gs_data_from_file[i, 2]
    Gs[i].cp.residual = 0
    Gs[i].storage.top = Gs_data_from_file[i, 3]
    Gs[i].storage.atual = Gs_data_from_file[i, 4]
    Gs[i].storage.residual = 0
    Gs[i].memory.top = Gs_data_from_file[i, 5]
    Gs[i].memory.atual = Gs_data_from_file[i, 6]
    Gs[i].memory.residual = 0
    Gs[i].VMs = []
    Gs[i].custo_beneficio = 0
    Gs[i].residual_band = Σ(BandM[i, :])
    # println("server ", i, " -> band = ", Gs[i].residual_band)
end

# println("\n", Gs, "\n", BandM, "\n")

VNR = Array{VirtualMachine}(0)
VNR_data_from_file = readdlm("VNR")
n_linhas = size(VNR_data_from_file, 1)
n_colunas = size(VNR_data_from_file, 2)
for i=1:n_linhas
    node = [Float64(k) for k in VNR_data_from_file[i, 1:n_colunas] if k != ""]
    VM = VirtualMachine(i, node[1], node[2], node[3], [], 0.0)
    deleteat!(node,[1,2,3])
    while length(node) > 0
        push!(VM.links, (node[1], node[2]))
        deleteat!(node, [1,2])
    end
    push!(VNR, VM)
end

set_lrc_index(VNR)

# println("VNR => ", VNR, "\n")

sort!(VNR)

# println("VNR => ", VNR, "\n")

array_criterias = Dict("custo_beneficio" => Array{Float64}(0), "cp_residual" => Array{Float64}(0), "storage_residual" => Array{Float64}(0), "memory_residual" => Array{Float64}(0))
PESOS = Dict("custo_beneficio" => 0.85, "cp_residual" => 0.05, "storage_residual" => 0.05, "memory_residual" => 0.05)

results = Dict{Int64, Array{Float64, 1}}()
substrate_network_changes = Array{Tuple{Int64, Int64, Float64}, 1}()
iteration = 1
i=1
first_host = -1
first_host_key = -1
while i <= length(VNR)
    VM = VNR[i]
    println(VM)
        println("@@ entrou no calc")
        calc_custo_beneficio(Gs, BandM, VM)
        calc_cp_residual(Gs, VM)
        calc_storage_residual(Gs, VM)
        calc_memory_residual(Gs, VM)

    if (!haskey(results, VM.id))
        array_criterias["custo_beneficio"] = [Float64(k.custo_beneficio) for k in Gs]
        MAUT_norm_criteria_max(array_criterias["custo_beneficio"])
        # quando tentando usar normalizacao de maximizacao, a banda de cada servidor fisico se esgota mais rapido
        MAUT_marginalUtilty_custo_beneficio(array_criterias["custo_beneficio"])

        array_criterias["cp_residual"] = [k.cp.residual for k in Gs]
        MAUT_norm_criteria_min(array_criterias["cp_residual"])
        MAUT_marginalUtility_cp_residual(array_criterias["cp_residual"])


        array_criterias["storage_residual"] = [k.storage.residual for k in Gs]
        MAUT_norm_criteria_min(array_criterias["storage_residual"])


        array_criterias["memory_residual"] = [k.memory.residual for k in Gs]
        MAUT_norm_criteria_min(array_criterias["memory_residual"])

        current_result = MAUT_globalUtility(array_criterias, PESOS)
        
        println("-------------------------\n",array_criterias, "\n", current_result, "\n-------------------------")

        push!(results, VM.id => current_result)        
    end
    println("Alocando VM ", VM.id)

    selected_server = index_of_higher_value_element(haskey(results, VM.id) ? results[VM.id] : current_result)
    println("Server selecionado: ", selected_server)
    if i == 1
        first_host = selected_server
        println("primeira iteracao, setando first_host and first_host_key\n.\n.\n.")
        first_host_key = VM.id
        results[VM.id][selected_server] = -1
    end
    println("Resultados para VM ", VM.id, " : ", results[VM.id])
    if !do_allocate(Gs, selected_server, VM)
        if (Σ(results[first_host_key]) == -length(results[first_host_key]))
            println("Fail to allocate VNR")
            deallocate(Gs, BandM, substrate_network_changes, first_host_key, first_host, results)
            return 0
        end
        println("Alocacao falhou, iniciando processo de desalocacao\n")
        deallocate(Gs, BandM, substrate_network_changes, first_host_key, first_host, results)
        delete!(results, VM.id)
        for key in keys(results)
            if(key != first_host_key)
                delete!(results, key)
            end
        end
        deleteat!(substrate_network_changes, 1:length(substrate_network_changes))
        println("Fim do processo de desalocacao")
        # println(results)
        # return
        i = 1
    else
        link_path = Array{Int64, 1}()
        if (i > 1)
            for link in VM.links
                if (haskey(results, link[1])) #if the node dst is already allocated
                    host_src = index_of_higher_value_element(results[VM.id])
                    host_dst = (link[1] == first_host_key) ? first_host : index_of_higher_value_element(results[link[1]])
                    link_path = dijkstra(BandM, host_src, host_dst, link[2])
                    println(">> PATH FROM ", host_src, " TO ", host_dst, " : ", link_path)
                    for node in link_path
                        BandM[host_src, node] -= link[2]
                        push!(substrate_network_changes, (host_src, node, link[2]))
                        host_src = node
                    end
                end
            end
            println("band allocated\n", BandM)
        end
        i += 1
    end
end
for i=1:length(Gs)
    println("\nSERVER ", i, "\n---------------")
    println("cp top: ", Gs[i].cp.top, "\ncp atual: ", Gs[i].cp.atual, "\ncp residual: ", Gs[i].cp.residual)
    println("storage top: ", Gs[i].storage.top, "\nstorage atual: ", Gs[i].storage.atual, "\nstorage residual: ", Gs[i].storage.residual)
    println("memory top: ", Gs[i].memory.top, "\nmemory atual: ", Gs[i].memory.atual, "\nmemory residual: ", Gs[i].memory.residual)
    println("custo beneficio: ", Gs[i].custo_beneficio)
    println("banda residual: ", Gs[i].residual_band)
    println("VMs: ", Gs[i].VMs)
end
