include("types.jl")
import Distributions
# type VirtualMachine
#     id::Int64
#     cp::Float64
#     storage::Float64
#     memory::Float64
#     links::Array{Tuple{Float64,Float64}} #links são enlaces virtuais seguindo a ordem (dest,weight)
#     LRC_index::Float64
# end
function display_events(events)
    for event in events
        println("EVENTS\n----------\n")
        if length(event) > 0
            for vnr in event
                println("\tVNR\n\t----------\n")
                for vm in vnr
                    println("\t\t", vm)
                end
                println("\t\n\t----------\n")
            end
        end
        println("\n----------\n")
    end
end

function generate_event_array(number_of_events::Int64, max_vms::Int64, bandwidth::Float64)
    M3 = Dict("medium" => [1, 4, 3.75], "large" => [2, 32, 7.5], "xlarge" => [4, 80, 15], "2xlarge" => [8, 160, 30])
    M3alias = Dict(1 => "medium", 2 => "large", 3 => "xlarge", 4 => "2xlarge")
    dist = map(x->Int64(round(clamp(x, 0, number_of_events))), rand(Distributions.Normal(number_of_events/2, 0.15 * number_of_events), number_of_events))

    println(dist)
    events = []
    vm_count = 0
    vnr_count = 0
    for i=1:number_of_events
        push!(events, [])
        events_for_i = length(filter(x->x==i, dist))
        for j=1:events_for_i
            vnr_count += 1
            VNR = Array{VirtualMachine}(0)
            expiration_time = rand(i+1:number_of_events + 2)
            # expiration_time = rand(i+1:i+2)
            number_of_vms = rand(3:max_vms)
            for vm=1:number_of_vms
                vm_count += 1
                vm_type = rand(1:4)
                push!(VNR, VirtualMachine(vm_count, M3[M3alias[vm_type]][1], M3[M3alias[vm_type]][2], M3[M3alias[vm_type]][3], [], 0, vnr_count, expiration_time))
            end
            for vm in VNR
                n_links = rand(0:number_of_vms)
                for link=1:n_links
                    dst = rand((vm_count - number_of_vms + 1):vm_count)
                    if(!(dst in [k[1] for k in vm.links]) && (dst != vm.id))        
                        band = rand(Distributions.Uniform(0.01*bandwidth, 0.25*bandwidth), 1)[1]
                        push!(vm.links, (dst, band))
                        push!(VNR[dst - (vm_count - number_of_vms)].links, (vm.id, band))
                    end
                end
            end
            push!(events[i], VNR)
        end
    end
    #display_events(events)
    return events
end
