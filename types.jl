type VirtualMachine
    id::Int64
    cp::Float64
    storage::Float64
    memory::Float64
    links::Array{Tuple{Float64,Float64}} #links são enlaces virtuais seguindo a ordem (dest,weight)
    LRC_index::Float64
    vnr_id::Float64
    expiration_time::Int64
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
