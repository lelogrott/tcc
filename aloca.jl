include("dijkstra.jl")
include("VNR_generator.jl")
include("fattree.jl")
include("types.jl")

K = 4
BANDWIDTH = 10240.0 
MAX_BANDWIDTH = 50
MIN_BANDWIDTH = 10
Σ = sum
Π = prod

array_criterias = Dict("custo_beneficio" => Array{Float64}(0), "cp_residual" => Array{Float64}(0), "storage_residual" => Array{Float64}(0), "memory_residual" => Array{Float64}(0))
PESOS = Dict("custo_beneficio" => 0.85, "cp_residual" => 0.05, "storage_residual" => 0.05, "memory_residual" => 0.05)

Base.isless(x::VirtualMachine, y::VirtualMachine) = x.LRC_index > y.LRC_index
Base.isequal(x::VirtualMachine, y::VirtualMachine) = x.LRC_index == y.LRC_index

function set_lrc_index(VNR)
    for vm in VNR
        vm.LRC_index = Σ([link[2] for link in vm.links if (length(link) > 0)]) + vm.cp
    end
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
    higher_value = -2
    # println("\tanalisando ", array)
    for i=1:length(array)
        if (array[i] > higher_value)
            index = i
            higher_value = array[i]
        end
    end
    println("\tindex de maior valor ", index + (1 + ((K/2)^2) + (K^2)))
    return index + (1 + ((K/2)^2) + (K^2))
end

function all_positive_residual_values_in(node::Server, VM::VirtualMachine)
    if (node.cp.residual < 0)
        println("fail in cp: ", node.cp.residual)
        return false
    end
    if (node.storage.residual < 0)
        println("fail in storage: ", node.storage.residual)
        return false
    end
    if (node.memory.residual < 0)
        println("fail in memory: ", node.storage.residual)
        return false
    end
    residual_band = node.residual_band - Σ([link[2] for link in VM.links])
    residual_band += Σ([(link[2] * 2) for link in VM.links if link[1] in [vm.id for vm in node.VMs]])
    if (residual_band < 0)
        println("fail in residual band: ", node.residual_band, " | band from vm: ", Σ([link[2] for link in VM.links]))
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
                readline(STDIN)
            end
            Gs[i].storage.atual = Gs[i].storage.residual
            Gs[i].memory.atual = Gs[i].memory.residual
            push!(Gs[i].VMs, VM)
            Gs[i].residual_band -= Σ([link[2] for link in VM.links])
            Gs[i].residual_band += Σ([(link[2] * 2) for link in VM.links if link[1] in [vm.id for vm in Gs[i].VMs]])
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

# BandM = readdlm("BandM_fisico")
BandM = fat_tree(K, BANDWIDTH)
f = open("Gs_fisico", "w")

for i=1:size(BandM, 1)#default configs for every host
    if(i >= (1 + ((K/2)^2) + (K^2)))
        @printf(f, "3000 3000 13500 13500 1500 1500\n") #hosts
    else
        @printf(f, "0 0 0 0 0 0\n") #switches
    end
end
close(f)

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

# # println("\n", Gs, "\n", BandM, "\n")
# READ FROM FILE
# VNR = Array{VirtualMachine}(0)
# VNR_data_from_file = readdlm("VNR")
# n_linhas = size(VNR_data_from_file, 1)
# n_colunas = size(VNR_data_from_file, 2)
# for i=1:n_linhas
#     node = [Float64(k) for k in VNR_data_from_file[i, 1:n_colunas] if k != ""]
#     VM = VirtualMachine(i, node[1], node[2], node[3], [], 0.0)
#     deleteat!(node,[1,2,3])
#     while length(node) > 0
#         push!(VM.links, (node[1], node[2]))
#         deleteat!(node, [1,2])
#     end
#     push!(VNR, VM)
# end

events = generate_event_array(20, 5, BANDWIDTH)

function deallocate_expired_vms(Gs::Array{Server}, current_time::Int64, BandM_in_time::Dict{Int64, Array{Tuple{Int64, Int64, Float64}, 1}}, BandM::Array{Float64, 2})
    deleted_vms = Array{VirtualMachine}(0)
    hosts_vms = Dict()
    for i=Int64(1 + ((K/2)^2) + (K^2)):length(Gs)
        new_VMs = Array{VirtualMachine}(0)
        for j=1:length(Gs[i].VMs)
            vm = Gs[i].VMs[j]
            if (vm.expiration_time == current_time)
                hosts_vms[vm.id] = i
                push!(deleted_vms, vm)
                Gs[i].cp.atual += vm.cp
                Gs[i].storage.atual += vm.storage
                Gs[i].memory.atual += vm.memory
            else
                push!(new_VMs, vm)
            end
        end
        Gs[i].VMs = deepcopy(new_VMs)
    end
    vnr_ids = sort(unique([k.vnr_id for k in deleted_vms]))
    for id in vnr_ids
        if(haskey(BandM_in_time, id))
            for (src, dst, wght) in BandM_in_time[id]
                BandM[src, dst] += wght
            end
        end
    end
    return Gs, BandM, deleted_vms
end

BandM_in_time = Dict{Int64, Array{Tuple{Int64, Int64, Float64}, 1}}()
current_time = 0
REJECTED_VNRS = Array{Array{VirtualMachine}}(0)
for event in events
    current_time += 1

    println("\nTEMPO: ", current_time, "\n")
    # readline(STDIN)

    gs, bandm, deleted_vms = deallocate_expired_vms(Gs, current_time, BandM_in_time, BandM)

    if (!isempty(deleted_vms))
        println("\nDELETED VMS AT TIME: ", current_time, "\n")
        for vm in deleted_vms
            println(vm)
        end
    end

    BandM = deepcopy(bandm)
    Gs = deepcopy(gs)

    if (length(event) == 0)
        #do stuff
    else
        for VNR in event
            # println("VNR => ", VNR, "\n")
            
            set_lrc_index(VNR)

            sort!(VNR)

            # println("VNR => ", VNR, "\n")

            results = Dict{Int64, Array{Float64, 1}}()
            substrate_network_changes = Array{Tuple{Int64, Int64, Float64}, 1}(0)
            iteration = 1
            i=1
            first_host = -1
            first_host_key = -1

            # SET STABLE VARIABLES
            STABLE_BandM = deepcopy(BandM)
            STABLE_Gs = deepcopy(Gs)
            while i <= length(VNR)
                VM = VNR[i]
                println(VM)
                only_hosts = deepcopy(Gs[Int64(1 + ((K/2)^2) + (K^2)):length(Gs)])
                calc_custo_beneficio(only_hosts, BandM, VM)
                calc_cp_residual(only_hosts, VM)
                calc_storage_residual(only_hosts, VM)
                calc_memory_residual(only_hosts, VM)
                if (!haskey(results, VM.id))
                    array_criterias["custo_beneficio"] = [Float64(k.custo_beneficio) for k in only_hosts]
                    MAUT_norm_criteria_max(array_criterias["custo_beneficio"])
                    # quando tentando usar normalizacao de maximizacao, a banda de cada servidor fisico se esgota mais rapido
                    MAUT_marginalUtilty_custo_beneficio(array_criterias["custo_beneficio"])

                    array_criterias["cp_residual"] = [k.cp.residual for k in only_hosts]
                    MAUT_norm_criteria_min(array_criterias["cp_residual"])
                    MAUT_marginalUtility_cp_residual(array_criterias["cp_residual"])


                    array_criterias["storage_residual"] = [k.storage.residual for k in only_hosts]
                    MAUT_norm_criteria_min(array_criterias["storage_residual"])


                    array_criterias["memory_residual"] = [k.memory.residual for k in only_hosts]
                    MAUT_norm_criteria_min(array_criterias["memory_residual"])

                    current_result = MAUT_globalUtility(array_criterias, PESOS)
                    
                    # println("-------------------------\n",array_criterias, "\n", current_result, "\n-------------------------")

                    push!(results, VM.id => current_result) #store the results in a array with the VM's id as key
                end
                Gs[Int64(1 + ((K/2)^2) + (K^2)):length(Gs)] = deepcopy(only_hosts)
                # println("Alocando VM ", VM.id)

                selected_server = Int64(index_of_higher_value_element(haskey(results, VM.id) ? results[VM.id] : current_result))
                # STABLE_HOST = deepcopy(Gs[selected_server])
                # println("Server selecionado: ", selected_server)
                if (i == 1)
                    first_host = selected_server
                    # println("primeira iteracao, setando first_host and first_host_key\n.\n.\n.")
                    first_host_key = VM.id
                    #results[VM.id][selected_server] = -1
                end
                # println("Resultados para VM ", VM.id, " : ", results[VM.id])

                if !do_allocate(Gs, selected_server, VM)
                    # Gs[i].residual_band += Σ([link[2] for link in VM.links])
                    # Gs[i] = STABLE_HOST
                    if (Σ(results[first_host_key]) == -length(results[first_host_key]))
                        # println("Fail to allocate VNR")
                        a = copy(STABLE_Gs)            
                        BandM = copy(STABLE_BandM)
                        push!(REJECTED_VNRS, VNR)
                        break
                    end
                    # println("Alocacao falhou, partindo para proxima opcao\n")

                    #if(selected_server > 0)
                    results[VM.id][(selected_server - Int64(1 + ((K/2)^2) + (K^2)))] = -1
                    #end
                    if (Σ(results[VM.id]) == -length(results[VM.id]))
                        # println("Alocacao falhou. partindo para proxima melhor opcao do primeiro server\n")
                        Gs = deepcopy(STABLE_Gs)
                        # Gs[i] = STABLE_HOST
                        BandM = deepcopy(STABLE_BandM)
                        for key in keys(results)
                            if(key != first_host_key)
                                delete!(results, key)
                            else
                                # println("MATCHING KEY")
                                # println(results[key])
                            end
                        end
                        deleteat!(substrate_network_changes, 1:length(substrate_network_changes))
                        # println("Fim do processo de desalocacao")
                        i = 1
                    end
                    # println(results)
                    # return
                else
                    link_path = Array{Int64, 1}()
                    if (i > 1)
                        for link in VM.links
                            if (haskey(results, link[1])) #if the node dst is already allocated
                                host_src = Int64(index_of_higher_value_element(results[VM.id]))
                                host_dst = Int64((link[1] == first_host_key) ? first_host : index_of_higher_value_element(results[link[1]]))
                                link_path = dijkstra(BandM, host_src, host_dst, link[2])
                                # println(">> PATH FROM ", host_src, " TO ", host_dst, " WITH WIEIGHT ", link[2], " : ", link_path)
                                if (isempty(link_path) && host_src != host_dst)
                                    # println("Nao foi possível alocar banda entre VM ", VM.id, " e ", link[1], "servidores ", host_src, " e ", host_dst, "\n")
                                    # println("Restaurando valores de banda da rede fisica e desalocando VM do server \n")
                                    for (src, dst, wght) in substrate_network_changes
                                        BandM[src, dst] += wght
                                    end
                                    results[VM.id][(host_src - Int64(1 + ((K/2)^2) + (K^2)))] = -1
                                    i-=1
                                    break
                                else
                                    for node in link_path
                                        BandM[host_src, node] -= link[2]
                                        push!(substrate_network_changes, (host_src, node, link[2]))
                                        host_src = node
                                    end
                                    BandM_in_time[VM.vnr_id] = substrate_network_changes
                                    # println("band allocated\n", BandM)
                                end
                            end
                        end
                    end
                    i += 1
                end
                #readline(STDIN)
            end
        end
    end
end
for i=Int64(1 + ((K/2)^2) + (K^2)):length(Gs)
# for i=1:length(Gs)
    println("\nSERVER ", i, "\n---------------")
    println("cp top: ", Gs[i].cp.top, "\ncp atual: ", Gs[i].cp.atual, "\ncp residual: ", Gs[i].cp.residual)
    println("storage top: ", Gs[i].storage.top, "\nstorage atual: ", Gs[i].storage.atual, "\nstorage residual: ", Gs[i].storage.residual)
    println("memory top: ", Gs[i].memory.top, "\nmemory atual: ", Gs[i].memory.atual, "\nmemory residual: ", Gs[i].memory.residual)
    println("custo beneficio: ", Gs[i].custo_beneficio)
    println("banda residual: ", Gs[i].residual_band)
    println("VMs: \n")
    for vm in Gs[i].VMs
        println(vm)
    end
end
for vnr in REJECTED_VNRS
    println("VNR\n----------\n")
    for vm in vnr
        println("\t", vm, "\n")
    end
end
# println("FINAL BAND MATRIX\n")
# display(BandM[Int64(1 + ((K/2)^2) + (K^2)):length(Gs), 1:Int64(1 + ((K/2)^2) + (K^2))])


