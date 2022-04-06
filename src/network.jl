"""
    connect_nodes!(net::MetaDiGraph, arcs...)

Connect nodes identified by `arc` pairs (e.g., `1 => 2`).
"""
function connect_nodes!(net::AbstractGraph, arcs...)
    for (source, sink) in arcs
        add_edge!(net, source, sink) #add edge between the two
    end
end