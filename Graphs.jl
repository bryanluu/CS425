module Graphs

using LightGraphs, MetaGraphs, GraphRecipes, Plots
import MetaGraphs

export Graph;
struct Graph{T}
    _idx_map::Dict{T, Int} # Map vertices to indices
    _pos_map::Dict{T, Tuple{N, N}} where {N <: Real}; # Map vertices to positions
    is_directed::Bool; # whether the graph is directed
    nodetype::DataType; # data type of node
    _metagraph::Union{MetaGraph, MetaDiGraph}; # Internal metagraph

    # Construct graph from:
    # nodemap: a map of nodes to positions
    # L (optional): an adjacency list of edges, with weights (L[u][v] is the weight on edge uv)
    # directed: whether the graph is directed or not
    function Graph(nodemap::Dict{T, Tuple{N, N}}; L::Dict{T, Dict{T, M}}=Dict(v => Dict{T, Int}() for v in keys(nodemap)),
        directed::Bool=false) where {T, N <: Real, M <: Number}
        if length(nodemap) <= 1
            error("Invalid-size dictionary given...(|V| must be at least 2).");
        end
        nodes = keys(nodemap);
        idx_map = Dict(v=>i for (i, v) in enumerate(nodes));
        pos_map = nodemap;
        mg = directed ? MetaDiGraph(SimpleDiGraph(length(nodes))) : MetaGraph(SimpleGraph(length(nodes)));
        for (n, i) in idx_map
            set_prop!(mg, i, :name, n);
        end
        weightfield!(mg, :weight); # set weight to correspond to :weight
        G = new{T}(idx_map, pos_map, directed, T, mg);
        skipped = T[];
        for u in keys(L)
            if !in(u, nodes)
                push!(skipped, u);
                continue # skip key if not a valid node
            end
            for v in keys(L[u])
                if !in(v, nodes)
                    push!(skipped, v);
                    continue
                end
                w = L[u][v]; # edge weight
                connect!(G, u, v, w)
            end
        end
        if !isempty(skipped)
            println();
            println("Warning: adjacency list given does not match nodes...only valid edges added.");
            println("Skipped: $(skipped)\n");
        end
        G;
    end
end

export V, nv;

V(g::Graph) = Set(keys(g._idx_map));

nv(g::Graph) = MetaGraphs.nv(g._metagraph);

export connect!;

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

export E, ne;

E(g::Graph) = Set([(get_prop(g._metagraph, src(e), :name), get_prop(g._metagraph, dst(e), :name))
                  for e in edges(g._metagraph)]);

ne(g::Graph) = MetaGraphs.ne(g._metagraph);

export neighbors;

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

export plot;

function plot(g::Graph; args...)
    N = nv(g);

    mg = copy(g._metagraph);
    for v in 2:N
        add_edge!(mg, 1, v);
    end

    ew(s, d, w) = has_edge(g._metagraph, s, d) ? 1 : 0;

    rev_D = Dict(i => v for (v, i) in g._idx_map);
    x = [g._pos_map[rev_D[i]][1] for i in 1:N];
    y = [g._pos_map[rev_D[i]][2] for i in 1:N];
    names = [rev_D[i] for i in 1:N];

    graphplot(mg, x=x, y=y, names=names, ew=ew, args...);
end

end
