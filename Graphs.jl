module Graphs

using LightGraphs, MetaGraphs

struct Graph{T}
    _idx_map::Dict{T, Int} # Map vertices to indices
    _pos_map::Dict{T, Tuple{N, N}} where {N <: Real}; # Map vertices to positions
    is_directed::Bool; # whether the graph is directed
    nodetype::DataType; # data type of node
    _metagraph::MetaGraph; # Internal metagraph
    Graph(nodes::Dict{T, Tuple{N, N}}, directed=false) where {T, N <: Real} = new{T}(
        Dict(v=>i for (i, v) in enumerate(keys(nodes))), nodes, directed, T,
        MetaGraph(directed ? SimpleGraph(length(keys(nodes))) : SimpleDiGraph(length(keys(nodes)))));
end

V(g::Graph) = Set(keys(g._idx_map));

nv(g::Graph) = nv(g._metagraph);

function connect!(g::Graph{T}, s::T, d::T, weight::N = 1) where {T, N <: Real}
    mg = g._metagraph;
    success = add_edge!(mg, s, d);
    if !success
        return false
    else
        set_prop!(mg, s, d, :weight, weight);
    end
end

function connect!(g::Graph{T}, e::Tuple{T, T}, weight::N = 1) where {T, N <: Real}
    connect!(g, e..., weight);
end

E(g::Graph) = Set(edges(g._metagraph));

ne(g::Graph) = ne(g._metagraph);

neighbors(g::Graph{T}, v::T) where {T} = neighbors(g, g._idx_map[v]);

end
