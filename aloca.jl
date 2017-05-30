include("dijkstra.jl")
include("VNR_generator.jl")
include("fattree.jl")
include("types.jl")

K = 8
BANDWIDTH = 1024.0 
NUMBER_OF_EVENTS = 100
NUMBER_OF_REQUISITIONS = 150
MAX_VMS = 5
Σ = sum
Π = prod

array_criterias = Dict("custo_beneficio" => Array{Float64}(0), "cp_residual" => Array{Float64}(0), "storage_residual" => Array{Float64}(0), "memory_residual" => Array{Float64}(0))
PESOS = Dict("custo_beneficio" => 0.0, "cp_residual" => 1.0, "storage_residual" => 0.0, "memory_residual" => 0.0, "band_residual" => 0.0)

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

function MAUT_marginalUtility_band_residual(array_criteria::Array{Float64})
    return copy!(array_criteria, pmap(x -> (sqrt(x)), array_criteria))
end

function MAUT_globalUtility(array_criterias, pesos)
    copy!(array_criterias["custo_beneficio"], array_criterias["custo_beneficio"]*pesos["custo_beneficio"])
    copy!(array_criterias["cp_residual"], array_criterias["cp_residual"]*pesos["cp_residual"])
    copy!(array_criterias["storage_residual"], array_criterias["storage_residual"]*pesos["storage_residual"])
    copy!(array_criterias["memory_residual"], array_criterias["memory_residual"]*pesos["memory_residual"])
    copy!(array_criterias["band_residual"], array_criterias["band_residual"]*pesos["band_residual"])

    result = array_criterias["custo_beneficio"] + array_criterias["cp_residual"] + array_criterias["memory_residual"] + array_criterias["band_residual"] + array_criterias["storage_residual"]
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

function calc_residual_band(Gs::Array{Server, 1}, BandM::Array{Float64, 2})
    for i=1:length(Gs)
        Gs[i].residual_band = Σ(BandM[i+Int64(((K/2)^2) + (K^2)), :])
    end
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
    # println("\tindex de maior valor ", index + (((K/2)^2) + (K^2)))
    return index + (((K/2)^2) + (K^2))
end

function all_positive_residual_values_in(node::Server, VM::VirtualMachine)
    if (node.cp.residual < 0)
        # println("fail in cp: ", node.cp.residual)
        return false
    end
    if (node.storage.residual < 0)
        # println("fail in storage: ", node.storage.residual)
        return false
    end
    if (node.memory.residual < 0)
        # println("fail in memory: ", node.storage.residual)
        return false
    end
    residual_band = node.residual_band - Σ([link[2] for link in VM.links])
    residual_band += Σ([(link[2] * 2) for link in VM.links if link[1] in [vm.id for vm in node.VMs]])
    if (residual_band < 0)
        # println("fail in residual band: ", node.residual_band, " | band from vm: ", Σ([link[2] for link in VM.links]))
        return false
    end
    return true
end

function do_allocate(Gs::Array{Server, 1}, selected_server::Int64, VM::VirtualMachine)
    ret = false
    for i=1:length(Gs)
        if (i == selected_server && all_positive_residual_values_in(Gs[i], VM))
            Gs[i].cp.atual = Gs[i].cp.residual
            Gs[i].storage.atual = Gs[i].storage.residual
            Gs[i].memory.atual = Gs[i].memory.residual
            push!(Gs[i].VMs, VM)
            ret = true
        end
        Gs[i].storage.residual = 0
        Gs[i].cp.residual = 0
        Gs[i].memory.residual = 0
        Gs[i].custo_beneficio = 0
    end
    return ret
end

function deallocate_expired_vms(Gs::Array{Server}, current_time::Int64, BandM_in_time::Dict{Int64, Array{Tuple{Int64, Int64, Float64}, 1}}, BandM::Array{Float64, 2})
    deleted_vms = Array{VirtualMachine}(0)
    hosts_vms = Dict()
    for i=Int64(1 + ((K/2)^2) + (K^2)):length(Gs)
        new_VMs = Array{VirtualMachine}(0)
        allocated_vms = [vm.id for vm in Gs[i].VMs]
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
    return Gs, BandM, deleted_vms, vnr_ids
end
# BandM = readdlm("BandM_fisico")
for overbooking in [26,28,33]
    # filename = string("out_", overbooking, ".txt")
    println(overbooking)
    fp = open("out.txt", "a")
    results_for_plot = Array{Float64,1}(0)
    for iteration=1:100
        BandM = fat_tree(K)
        f = open("Gs_fisico", "w")

        for i=1:size(BandM, 1)#default configs for every host
            if(i >= (1 + ((K/2)^2) + (K^2)))
                @printf(f, "%d %d 2048 2048 256 256\n", overbooking, overbooking) #hosts #overbooking de 10%!
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



        events = generate_event_array(NUMBER_OF_EVENTS, NUMBER_OF_REQUISITIONS, MAX_VMS, BANDWIDTH)


        BandM_in_time = Dict{Int64, Array{Tuple{Int64, Int64, Float64}, 1}}()
        current_time = 0
        REJECTED_VNRS = Array{Array{VirtualMachine}}(0)
        total_allocated_vnrs = 0
        total_deleted = 0
        # fp = open("deleted_vnrs.data", "w")
        # fpp = open("rejected_vnrs.data", "w")
        for event in events
            current_time += 1

            # println("\nTEMPO: ", current_time, "\n")
            # readline(STDIN)

            gs, bandm, deleted_vms, vnr_ids = deallocate_expired_vms(Gs, current_time, BandM_in_time, BandM)

            if (!isempty(deleted_vms))
                # println("\nDELETED VMS AT TIME: ", current_time, "\n")
                for vm in deleted_vms
                    # println(vm)
                end
            end
            # total_vnrs -= length(vnr_ids)
            total_deleted += length(vnr_ids)
            
            # @printf(fp, "%d\n",  length(vnr_ids))
            # @printf(fpp, "%d\n",  length(REJECTED_VNRS))
            BandM = deepcopy(bandm)
            Gs = deepcopy(gs)

            if (length(event) == 0)
                #do stuff
            else
                for VNR in event
                    did_allocate = true
                    # println("VNR => ", VNR, "\n")
                    # total_vnrs += 1
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
                        # println(VM)
                        only_hosts = deepcopy(Gs[Int64(1 + ((K/2)^2) + (K^2)):length(Gs)])
                        calc_custo_beneficio(only_hosts, BandM, VM)
                        calc_cp_residual(only_hosts, VM)
                        calc_storage_residual(only_hosts, VM)
                        calc_memory_residual(only_hosts, VM)
                        calc_residual_band(only_hosts, BandM)
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

                            array_criterias["band_residual"] = [k.residual_band for k in only_hosts]
                            MAUT_norm_criteria_max(array_criterias["band_residual"])
                            MAUT_marginalUtility_band_residual(array_criterias["band_residual"])

                            current_result = MAUT_globalUtility(array_criterias, PESOS)
                            
                            # println("-------------------------\n",array_criterias, "\n", current_result, "\n-------------------------")

                            push!(results, VM.id => current_result) #store the results in a array with the VM's id as key
                        end
                        Gs[Int64(1 + ((K/2)^2) + (K^2)):length(Gs)] = deepcopy(only_hosts)
                        selected_server = Int64(index_of_higher_value_element(results[VM.id]))
                        if (i == 1)
                            first_host = selected_server
                            # println("primeira iteracao, setando first_host and first_host_key\n.\n.\n.")
                            first_host_key = VM.id
                            results[VM.id][selected_server - Int64(((K/2)^2) + (K^2))] = -1
                        end
                        # println("Resultados para VM ", VM.id, " : ", results[VM.id])

                        if (!do_allocate(Gs, selected_server, VM) || (Σ(results[VM.id]) == -length(results[VM.id])))
                            if (Σ(results[first_host_key]) == -length(results[first_host_key]))
                                did_allocate = false
                                BandM = deepcopy(STABLE_BandM)
                                push!(REJECTED_VNRS, VNR)
                                break
                            end
                            results[VM.id][(selected_server - Int64(((K/2)^2) + (K^2)))] = -1
                            if (Σ(results[VM.id]) == -length(results[VM.id]))
                                # println("Alocacao falhou. partindo para proxima melhor opcao do primeiro server\n")
                                Gs = deepcopy(STABLE_Gs)
                                BandM = deepcopy(STABLE_BandM)
                                for key in keys(results)
                                    if(key != first_host_key)
                                        delete!(results, key)
                                    else
                                        # println("MATCHING KEY")
                                        # println(results[key])
                                    end
                                end
                                i = 1
                            end
                        else
                            link_path = Array{Int64, 1}()
                            if (i > 1)
                                for k=1:length(VM.links)
                                    link = VM.links[k]
                                    if (haskey(results, link[1])) #if the node dst is already allocated
                                        host_src = Int64(index_of_higher_value_element(results[VM.id]))
                                        host_dst = Int64((link[1] == first_host_key) ? first_host : index_of_higher_value_element(results[link[1]]))
                                        link_path = dijkstra(BandM, host_src, host_dst, link[2])
                                        if (isempty(link_path) && host_src != host_dst)
                                            # println("Nao foi possível alocar banda entre VM ", VM.id, " e ", link[1], "servidores ", host_src, " e ", host_dst, "\n")
                                            # println("Restaurando valores de banda da rede fisica e desalocando VM do server \n")
                                            for (src, dst, wght) in substrate_network_changes
                                                BandM[src, dst] += wght
                                            end
                                            results[VM.id][(host_src - Int64(((K/2)^2) + (K^2)))] = -1
                                            i-=1
                                            break
                                        else
                                            for node in link_path
                                                BandM[host_src, node] -= link[2]
                                                BandM[node, host_src] -= link[2]                                                                                                                 
                                                push!(substrate_network_changes, (host_src, node, link[2]))
                                                push!(substrate_network_changes, (node, host_src, link[2]))
                                                
                                                if(!haskey(BandM_in_time, VM.vnr_id))
                                                    BandM_in_time[VM.vnr_id] = Array{Tuple{Int64, Int64, Float64}, 1}(0)
                                                end
                                                push!(BandM_in_time[VM.vnr_id], (host_src, node, link[2]))
                                                push!(BandM_in_time[VM.vnr_id], (node, host_src, link[2]))
                                                
                                                host_src = node
                                            end
                                        end
                                    end
                                end
                                deleteat!(substrate_network_changes, 1:length(substrate_network_changes))
                            end
                            i += 1
                        end
                    end
                    if(did_allocate)
                        total_allocated_vnrs+=1
                    end
                end
            end
        end
        push!(results_for_plot, total_allocated_vnrs)
        println(total_allocated_vnrs)
    end
    @printf(fp, "%d %d %d\n", overbooking, mean(results_for_plot), std(results_for_plot))
    close(fp)
end
# close(fp)
# close(fpp)
# for i=Int64(1 + ((K/2)^2) + (K^2)):length(Gs)
# # for i=1:length(Gs)
#     println("\nSERVER ", i, "\n---------------")
#     println("cp top: ", Gs[i].cp.top, "\ncp atual: ", Gs[i].cp.atual, "\ncp residual: ", Gs[i].cp.residual)
#     println("storage top: ", Gs[i].storage.top, "\nstorage atual: ", Gs[i].storage.atual, "\nstorage residual: ", Gs[i].storage.residual)
#     println("memory top: ", Gs[i].memory.top, "\nmemory atual: ", Gs[i].memory.atual, "\nmemory residual: ", Gs[i].memory.residual)
#     println("custo beneficio: ", Gs[i].custo_beneficio)
#     println("banda residual: ", Gs[i].residual_band)
#     println("VMs: \n")
#     for vm in Gs[i].VMs
#         println(vm)
#     end
# end
# for vnr in REJECTED_VNRS
#     println("VNR\n----------\n")
#     for vm in vnr
#         println("\t", vm, "\n")
#     end
# end
# println("FINAL BAND MATRIX\n")
# display(BandM[Int64(1 + ((K/2)^2) + (K^2)):length(Gs), 1:Int64(1 + ((K/2)^2) + (K^2))])


