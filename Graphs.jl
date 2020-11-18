module Graphs

using LightGraphs, MetaGraphs
import MetaGraphs

struct Graph{T}
    _idx_map::Dict{T, Int} # Map vertices to indices
    _pos_map::Dict{T, Tuple{N, N}} where {N <: Real}; # Map vertices to positions
    is_directed::Bool; # whether the graph is directed
    nodetype::DataType; # data type of node
    _metagraph::MetaGraph; # Internal metagraph

    function Graph(nodes::Dict{T, Tuple{N, N}}, directed::Bool=false) where {T, N <: Real}
        idx_map = Dict(v=>i for (i, v) in enumerate(keys(nodes)));
        pos_map = nodes;
        mg = MetaGraph(directed ? SimpleGraph(length(keys(nodes))) : SimpleDiGraph(length(keys(nodes))));
        for (n, i) in idx_map
            set_prop!(mg, i, :name, n);
        end
        new{T}(idx_map, pos_map, directed, T, mg);
    end
end

V(g::Graph) = Set(keys(g._idx_map));

nv(g::Graph) = MetaGraphs.nv(g._metagraph);

function connect!(g::Graph{T}, s::T, d::T, weight::N = 1) where {T, N <: Real}
    mg = g._metagraph;
    success = add_edge!(mg, g._idx_map[s], g._idx_map[d]);
    if !success
        return false
    else
        set_prop!(mg, g._idx_map[s], g._idx_map[d], :weight, weight);
    end
end

function connect!(g::Graph{T}, e::Tuple{T, T}, weight::N = 1) where {T, N <: Real}
    connect!(g, e..., weight);
end

E(g::Graph) = Set([(get_prop(g._metagraph, src(e), :name), get_prop(g._metagraph, dst(e), :name))
                  for e in edges(g._metagraph)]);

ne(g::Graph) = MetaGraphs.ne(g._metagraph);

function neighbors(g::Graph{T}, v::T) where {T}
    [get_prop(g._metagraph, i, :name) for i in MetaGraphs.neighbors(g._metagraph, g._idx_map[v])];
end

function BFS(g::Graph{T}, s::T, d::T) where {T}
    if !in(s, keys(g._idx_map)) || !in(d, keys(g._idx_map))
        return T[];
    end

    p = Dict{T, Union{T, Nothing}}(v => nothing for v in keys(g._idx_map)); # predecessor vector
    p[s] = s;
    Q = T[s];

    function getpath(p::Dict{T, Union{T, Nothing}})
        curr = d;
        path = T[];
        while curr != s
            pushfirst!(path, curr);
            curr = p[curr];
        end
        pushfirst!(path, s);
        return path;
    end

    while !isempty(Q)
        curr = popfirst!(Q);

        if curr == d
            return getpath(p);
        end

        for n in neighbors(g, curr)
            if p[n] == nothing
                push!(Q, n);
                p[n] = curr;
            end

            if n == d
                return getpath(p);
            end
        end

    end

end

end
