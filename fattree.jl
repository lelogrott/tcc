function fat_tree(k::Int64, bandwidth::Float64)
    n_nodes = Int64(k^2 + (k^3 + k^2)/4)

    matrix = zeros(Float64, n_nodes, n_nodes)
    aggregation_switches_end = Int64(((k/2)^2) + (k^2)/2 - (k/2)) #subtraindo ultimos switches.
    offsetter = 1
    control = 0
    for i=1:Int64((k/2)^2) #linking core switches to aggregation switches
        for aggregation_switches_pos=Int64(((k/2)^2)+offsetter):Int64(k/2):Int64(aggregation_switches_end + offsetter)
            matrix[i, aggregation_switches_pos] = bandwidth
            matrix[aggregation_switches_pos, i] = bandwidth
            aggregation_switches_pos += k/2
        end
        control += 1
        if control == k/2 # k/2 ports to connect with
            offsetter += 1
            control = 0
        end
    end

    edges_switches_pos = Int64(((k/2)^2) + (k^2)/2 + 1)
    control = 0
    refer_aggregate_pos = edges_switches_pos - Int64((k^2)/2)
    for i=edges_switches_pos:(edges_switches_pos+Int64((k^2)/2) - 1) #linking edges switches to aggregation switches
        if control == k/2
            refer_aggregate_pos = i - Int64((k^2)/2)
            control = 0
        end
        for offsetter=0:Int64(k/2 - 1)
            matrix[i, refer_aggregate_pos + offsetter] = bandwidth
            matrix[refer_aggregate_pos + offsetter, i] = bandwidth
        end
        control += 1
    end

    hosts_pos = edges_switches_pos + Int64((k^2)/2)
    control = 0
    offsetter = 0
    for i=hosts_pos:(hosts_pos + Int64((k^3)/4) - 1)
        if control == k/2
            offsetter += 1
            control = 0
        end
        matrix[i, edges_switches_pos + offsetter] = bandwidth
        matrix[edges_switches_pos + offsetter, i] = bandwidth
        control += 1
    end
    return matrix
end
