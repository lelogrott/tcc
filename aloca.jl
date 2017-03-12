MAX_BANDWIDTH = 50
MIN_BANDWIDTH = 10
Σ = sum
Π = prod

type VirtualMachine
    cp::Float64
    links::Array{Tuple{Float64,Float64}} #links são enlaces virtuais seguindo a ordem (dest,weight)
    LRC_index::Float64
end

type Server
    cp_top::Float64
    cp_atual::Float64
    residual::Float64
    custo_beneficio::Float64
    desempenho::Float64
    VMs::Array{VirtualMachine}
end

Base.isless(x::VirtualMachine, y::VirtualMachine) = x.LRC_index > y.LRC_index
Base.isequal(x::VirtualMachine, y::VirtualMachine) = x.LRC_index == y.LRC_index

function set_lrc_index(VNR)
    for vi in VNR
        vi.LRC_index = Σ(link[2] for link in vi.links) + vi.cp
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
        Gs[i] = Server(100, rand(80:100), 0, 0, 0, [])
    end
    return Gs, BandM
end

function gera_VNR(n_nodes)
    VNR = Array{VirtualMachine}(0)
    CP_MIN = 1
    CP_MAX = 10
    BAND_MIN = 5
    BAND_MAX = 20
    for i=1:n_nodes
        links = rand(0:n_nodes)
        VM = VirtualMachine(rand(CP_MIN:CP_MAX), [], 0.0)
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
     ↑ = (x-> (x - min)/(max-min))
     return copy!(array_criteria, pmap(↑, array_criteria))
end

function MAUT_norm_criteria_min(array_criteria::Array{Float64})
    min = minimum(array_criteria)
    max = maximum(array_criteria)
    ↓ = (x-> 1+(min-x)/(max - min))
    return copy!(array_criteria, pmap(↓, array_criteria))
end

function MAUT_marginalUtilty_custo_beneficio(array_criteria::Array{Float64})
    return copy!(array_criteria, pmap(x -> (exp(x^2)-1)/exp(1), array_criteria))
end

function MAUT_marginalUtility_residual(array_criteria::Array{Float64})
    return copy!(array_criteria, pmap(x -> (exp(x^2)-1)/exp(1), array_criteria))
end

function MAUT_globalUtility(array_criterias, pesos)
    for i in array_criterias["custo_beneficio"]
        i *= pesos["custo_beneficio"]
    end
    for i in array_criterias["residual"]
        i *= pesos["residual"]
    end
    result = array_criterias["custo_beneficio"] + array_criterias["residual"]
    return result
end

function parallel_calc_residual(server, vnr_cp)
    server.residual = server.cp_atual - vnr_cp
end

function calc_residual(Gs::Array{Server},vnr::VirtualMachine)
    return pmap(x->parallel_calc_residual(x,vnr.cp),Gs)
end

function calc_custo_beneficio(Gs::Array{Server}, BandM::Array{Float64},vnr::VirtualMachine)
    Band_v = Σ([link[2] for link in vnr.links])
    for i=1:length(Gs)
        println(vnr.cp, " / ", Gs[i].cp_atual, " + ", Band_v, " / ", Σ([k for k in BandM[i,1:size(BandM,1)]]), " = ")
        Gs[i].custo_beneficio = (vnr.cp / Gs[i].cp_atual) + (Band_v / Σ([k for k in BandM[i,1:size(BandM,1)]]))
        println(" ", Gs[i].custo_beneficio, "\n\n")
    end
end

#C_{r}(p) = c_{p} - c_{v}, sendo:
#C_{r}(p) a capacidade residual de p;
#c_{p} = capacidade original de p (o total);
#c_{v} = capacidade solicitada pelo cliente (virtual)

#####################################
# TESTE DINAMICO
# MATRIZES GERADAS ALEATORIAMENTE
#

# n_nodes =5

# Gs,BandM = gera_Gs_BandM(n_nodes)
# println(Gs, "\n", BandM)

# MAUT_criteria_custo_beneficio = [Float64(k.custo_beneficio) for k in Gs]
# MAUT_criteria_residual = [Float64(k.residual) for k in Gs]

# println("custo: ", MAUT_criteria_custo_beneficio, "\nresidual: ", MAUT_criteria_residual, "")
# VNR  = gera_VNR(3)

# #println("\n", VNR, "\n")

# for node in VNR
#     println("cp: ", node.cp, "\n")
#     for link in node.links
#         println(link, "\n")
#     end
# end

########################################
# TESTE COM INPUT VIA FILE
#
#
#

BandM = readdlm("BandM_fisico")

Gs_data_from_file = readdlm("Gs_fisico")
n_servers = size(Gs_data_from_file,1)

Gs = Array{Server}(n_servers)
for i=1:length(Gs)
    Gs[i] = Server(0, 0, 0, 0, 0, [])
end
for i=1:length(Gs)
    j=1
    for current_field in fieldnames(Server)
        if current_field != :VMs
            setfield!(Gs[i], current_field, Gs_data_from_file[i, j])
            j += 1
        end
    end
end

println("\n", Gs, "\n", BandM, "\n")

VNR = Array{VirtualMachine}(0)
VNR_data_from_file = readdlm("VNR")
n_linhas = size(VNR_data_from_file, 1)
n_colunas = size(VNR_data_from_file, 2)
for i=1:n_linhas
    node = [Float64(k) for k in VNR_data_from_file[i, 1:n_colunas] if k != ""]
    VM = VirtualMachine(node[1], [], 0.0)
    deleteat!(node,1)
    while length(node) > 0
        push!(VM.links, (node[1], node[2]))
        deleteat!(node, [1,2])
    end
    push!(VNR, VM)
end

set_lrc_index(VNR)

println("VNR => ", VNR, "\n")

sort!(VNR)

println("VNR => ", VNR, "\n")

calc_custo_beneficio(Gs, BandM, VNR[1])
println("custo pos calc: ", [Float64(k.custo_beneficio) for k in Gs], "\n")
calc_residual(Gs, VNR[1])
println("residual pos calc: ", [Float64(k.residual) for k in Gs], "\n")

array_criteria_custo_beneficio = [Float64(k.custo_beneficio) for k in Gs]
MAUT_norm_criteria_max(array_criteria_custo_beneficio)
println("custo pos norm: ", array_criteria_custo_beneficio, "\n")


MAUT_marginalUtilty_custo_beneficio(array_criteria_custo_beneficio)
println("custo pos U: ", array_criteria_custo_beneficio, "\n")

array_criteria_residual = [k.residual for k in Gs]

MAUT_norm_criteria_max(array_criteria_residual)
println("residual pos norm: ", array_criteria_residual, "\n")
MAUT_marginalUtility_residual(array_criteria_residual)
println("residual pos U: ", array_criteria_residual, "\n")

PESOS = Dict("custo_beneficio" => 0.65, "residual" => 0.35)
array_criterias = Dict("custo_beneficio" => array_criteria_custo_beneficio, "residual" => array_criteria_residual)
result = MAUT_globalUtility(array_criterias, PESOS)
println(result)
