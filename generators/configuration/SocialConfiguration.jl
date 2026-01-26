module SocialConfiguration

using Graphs
using Random
using Logging

export social_configuration_generator,
       degree_sequence_is_graphical,
       random_spanning_tree

const INF = typemax(Int) ÷ 4

struct InfeasibleState <: Exception
    reason::String
end

@inline function sum_with_inf(values::Int...)
    total = 0
    for value in values
        if value == INF
            return INF
        end
        total += value
    end
    return total
end

@inline function residual_edge_count(remaining::Vector{Int})
    total_stubs = 0
    for r in remaining
        total_stubs += max(r, 0)
    end
    return total_stubs ÷ 2
end

struct HeuristicBranch
    name::Symbol
    edges::Vector{NTuple{2, Int}}
    weights::Vector{Float64}
    base_weight::Float64
end

"""
    degree_sequence_is_graphical(degrees::Vector{Int})

Check whether a degree sequence is graphical using the Havel-Hakimi algorithm.

# Arguments
- `degrees::Vector{Int}`: The degree sequence to validate

# Returns
- `Bool`: true if the sequence is graphical, false otherwise
"""
function degree_sequence_is_graphical(degrees::Vector{Int})
    seq = sort(degrees; rev = true)
    while !isempty(seq)
        d1 = first(seq)
        if d1 < 0
            return false
        elseif d1 == 0
            return all(==(0), seq)
        end
        seq = seq[2:end]
        if d1 > length(seq)
            return false
        end
        for i in 1:d1
            seq[i] -= 1
            if seq[i] < 0
                return false
            end
        end
        seq = sort(seq; rev = true)
    end
    return true
end

"""Convert a Prüfer sequence into an edge list."""
function _prufer_sequence_to_tree(seq::Vector{Int})
    n = length(seq) + 2
    degree = ones(Int, n)
    for v in seq
        1 <= v <= n || error("Node in Prüfer sequence out of range: $v")
        degree[v] += 1
    end

    leaves = [i for i in 1:n if degree[i] == 1]
    sort!(leaves)
    edges = NTuple{2, Int}[]

    for v in seq
        u = popfirst!(leaves)
        push!(edges, (u, v))
        degree[u] -= 1
        degree[v] -= 1
        if degree[v] == 1
            push!(leaves, v)
            sort!(leaves)
        end
    end

    push!(edges, (leaves[1], leaves[2]))
    return edges
end

function _sample_bounded_prufer_counts(capacities::Vector{Int}, rng::AbstractRNG)
    n = length(capacities)
    remaining = n - 2
    counts = zeros(Int, n)

    total_capacity = sum(capacities)
    if total_capacity < remaining
        throw(InfeasibleState("Degree sequence cannot host a spanning tree: need $remaining Prüfer slots, have $total_capacity"))
    end

    suffix_cap = zeros(Int, n + 1)
    for i in n:-1:1
        suffix_cap[i] = suffix_cap[i + 1] + capacities[i]
    end

    for i in 1:n-1
        remaining_capacity = suffix_cap[i + 1]
        min_i = max(0, remaining - remaining_capacity)
        max_i = min(capacities[i], remaining)
        if min_i > max_i
            throw(InfeasibleState("Cannot allocate Prüfer counts within degree bounds"))
        end
        counts[i] = rand(rng, min_i:max_i)
        remaining -= counts[i]
    end

    if remaining < 0 || remaining > capacities[n]
        throw(InfeasibleState("Cannot allocate Prüfer counts within degree bounds"))
    end
    counts[n] = remaining
    return counts
end

"""
    random_spanning_tree(degrees::Vector{Int}; rng::AbstractRNG = Random.GLOBAL_RNG)

Generate a random spanning tree while respecting degree upper bounds.

# Arguments
- `degrees::Vector{Int}`: Maximum degree for each node
- `rng::AbstractRNG`: Random number generator (default: Random.GLOBAL_RNG)

# Returns
- `Vector{NTuple{2, Int}}`: List of edges forming the spanning tree

# Throws
- `InfeasibleState`: If the degree sequence cannot form a spanning tree
"""
function random_spanning_tree(degrees::Vector{Int}; rng::AbstractRNG = Random.GLOBAL_RNG)
    n = length(degrees)
    if n <= 1
        return NTuple{2, Int}[]
    end

    if any(<(1), degrees)
        throw(InfeasibleState("Spanning tree requires every node to have degree >= 1"))
    end

    total_capacity = sum(degrees) - n
    if total_capacity < n - 2
        throw(InfeasibleState("Degree sequence cannot form a spanning tree: need at least $(n - 2) excess stubs, have $total_capacity"))
    end

    capacities = degrees .- 1
    counts = _sample_bounded_prufer_counts(capacities, rng)

    seq = Vector{Int}()
    sizehint!(seq, n - 2)
    for (node, cnt) in enumerate(counts)
        for _ in 1:cnt
            push!(seq, node)
        end
    end
    shuffle!(rng, seq)

    return _prufer_sequence_to_tree(seq)
end

canonical_pair(u::Int, v::Int) = u < v ? (u, v) : (v, u)

function _has_edge(graph::Vector{Set{Int}}, u::Int, v::Int)
    return v in graph[u]
end

function _add_edge!(graph::Vector{Set{Int}}, u::Int, v::Int)
    push!(graph[u], v)
    push!(graph[v], u)
end

function _bfs_distances(graph::Vector{Set{Int}}, src::Int)
    n = length(graph)
    dist = fill(INF, n)
    queue = Vector{Int}()
    head = 1
    dist[src] = 0
    push!(queue, src)
    while head <= length(queue)
        u = queue[head]
        head += 1
        du = dist[u]
        for v in graph[u]
            if dist[v] == INF
                dist[v] = du + 1
                push!(queue, v)
            end
        end
    end
    return dist
end

function _initialise_distance_cache(graph::Vector{Set{Int}})
    n = length(graph)
    cache = fill(INF, n, n)
    for src in 1:n
        dist = _bfs_distances(graph, src)
        for dst in 1:n
            cache[src, dst] = dist[dst]
        end
    end
    return cache
end

function _update_distance_cache_after_edge!(cache::Matrix{Int}, u::Int, v::Int)
    n = size(cache, 1)
    
    # Must copy distances before modifying cache
    dist_to_u = copy(cache[:, u])
    dist_to_v = copy(cache[:, v])
    dist_from_u = copy(cache[u, :])
    dist_from_v = copy(cache[v, :])
    
    cache[u, v] = 1
    cache[v, u] = 1
    
    # Update all pairs using the new edge as potential shortcut
    for x in 1:n
        d_xu = dist_to_u[x]
        d_xv = dist_to_v[x]
        for y in x+1:n
            current = cache[x, y]
            if current == 1
                continue
            end
            alt1 = sum_with_inf(d_xu, 1, dist_from_v[y])
            alt2 = sum_with_inf(d_xv, 1, dist_from_u[y])
            new_dist = min(current, alt1, alt2)
            if new_dist != current
                cache[x, y] = new_dist
                cache[y, x] = new_dist
            end
        end
    end
end

function _initialise_candidate_pairs(graph::Vector{Set{Int}}, remaining::Vector{Int})
    n = length(graph)
    pairs = Set{NTuple{2, Int}}()
    candidate_adj = Dict{Int, Set{NTuple{2, Int}}}()
    
    for i in 1:n-1
        for j in i+1:n
            if !_has_edge(graph, i, j) && remaining[i] > 0 && remaining[j] > 0
                pair = (i, j)
                push!(pairs, pair)
                
                # Build adjacency index for O(d) pruning
                if !haskey(candidate_adj, i)
                    candidate_adj[i] = Set{NTuple{2, Int}}()
                end
                if !haskey(candidate_adj, j)
                    candidate_adj[j] = Set{NTuple{2, Int}}()
                end
                push!(candidate_adj[i], pair)
                push!(candidate_adj[j], pair)
            end
        end
    end
    return pairs, candidate_adj
end

function _prune_node!(pairs::Set{NTuple{2, Int}}, adj::Dict{Int, Set{NTuple{2, Int}}}, node::Int)
    if haskey(adj, node)
        to_remove = collect(adj[node])
        for pair in to_remove
            delete!(pairs, pair)
        end
        delete!(adj, node)
        for pair in to_remove
            u, v = pair
            other = (u == node) ? v : u
            if haskey(adj, other)
                delete!(adj[other], pair)
            end
        end
    end
end

function _update_candidate_pairs!(pairs::Set{NTuple{2, Int}}, triadic_index::Dict{NTuple{2, Int}, Set{Int}},
    candidate_adj::Dict{Int, Set{NTuple{2, Int}}}, remaining::Vector{Int}, u::Int, v::Int)
    
    if remaining[u] == 0
        _prune_node!(pairs, candidate_adj, u)
    end
    if remaining[v] == 0
        _prune_node!(pairs, candidate_adj, v)
    end
    
    pair_uv = canonical_pair(u, v)
    if pair_uv in pairs
        delete!(pairs, pair_uv)
        if haskey(triadic_index, pair_uv)
            delete!(triadic_index, pair_uv)
        end
        
        if haskey(candidate_adj, u)
            delete!(candidate_adj[u], pair_uv)
        end
        if haskey(candidate_adj, v)
            delete!(candidate_adj[v], pair_uv)
        end
    end
end

function _initialise_triadic_index(graph::Vector{Set{Int}}, candidate_pairs::Set{NTuple{2, Int}})
    index = Dict{NTuple{2, Int}, Set{Int}}()
    for (u, v) in candidate_pairs
        common = intersect(graph[u], graph[v])
        if !isempty(common)
            index[(u, v)] = common
        end
    end
    return index
end

function _update_triadic_index!(index::Dict{NTuple{2, Int}, Set{Int}}, graph::Vector{Set{Int}},
    candidate_pairs::Set{NTuple{2, Int}}, u::Int, v::Int)
    
    for w in union(graph[u], graph[v])
        if w == u || w == v
            continue
        end
        
        for x in graph[w]
            if x == u || x == v
                continue
            end
            pair = canonical_pair(w, x)
            if pair in candidate_pairs
                if !haskey(index, pair)
                    index[pair] = Set{Int}()
                end
                common = intersect(graph[w], graph[x])
                if !isempty(common)
                    index[pair] = common
                else
                    delete!(index, pair)
                end
            end
        end
    end
end

function _build_degree_branch(candidate_pairs::Set{NTuple{2, Int}}, degrees::Vector{Int})
    edges = collect(candidate_pairs)
    weights = [Float64(degrees[u] * degrees[v]) for (u, v) in edges]
    return HeuristicBranch(:degree, edges, weights, 1.0)
end

function _build_spatial_branch(candidate_pairs::Set{NTuple{2, Int}}, distance_cache::Matrix{Int})
    edges = NTuple{2, Int}[]
    weights = Float64[]
    for pair in candidate_pairs
        u, v = pair
        d = distance_cache[u, v]
        if d != INF
            push!(edges, pair)
            push!(weights, 1.0 / Float64(d))
        end
    end
    return HeuristicBranch(:spatial, edges, weights, 1.0)
end

function _build_triadic_branch(triadic_index::Dict{NTuple{2, Int}, Set{Int}})
    edges = collect(keys(triadic_index))
    weights = [Float64(length(triadic_index[pair])) for pair in edges]
    return HeuristicBranch(:triadic, edges, weights, 1.0)
end

function _select_branch_type(alpha::Float64, beta::Float64, gamma::Float64,
    candidate_pairs::Set{NTuple{2, Int}}, triadic_index::Dict{NTuple{2, Int}, Set{Int}},
    rng::AbstractRNG)
    
    active_weights = Float64[]
    active_types = Symbol[]
    
    if !isempty(candidate_pairs) && alpha > 0
        push!(active_weights, alpha)
        push!(active_types, :degree)
    end
    if !isempty(candidate_pairs) && beta > 0
        push!(active_weights, beta)
        push!(active_types, :spatial)
    end
    if !isempty(triadic_index) && gamma > 0
        push!(active_weights, gamma)
        push!(active_types, :triadic)
    end
    
    if isempty(active_types)
        return nothing
    end
    
    total = sum(active_weights)
    r = rand(rng) * total
    cumsum_weight = 0.0
    for (i, w) in enumerate(active_weights)
        cumsum_weight += w
        if r <= cumsum_weight
            return active_types[i]
        end
    end
    return last(active_types)
end

function _sample_edge_from_branch(branch_type::Symbol, candidate_pairs::Set{NTuple{2, Int}},
    triadic_index::Dict{NTuple{2, Int}, Set{Int}}, distance_cache::Matrix{Int},
    degrees::Vector{Int}, rng::AbstractRNG)
    
    branch = if branch_type == :degree
        _build_degree_branch(candidate_pairs, degrees)
    elseif branch_type == :spatial
        _build_spatial_branch(candidate_pairs, distance_cache)
    elseif branch_type == :triadic
        _build_triadic_branch(triadic_index)
    else
        error("Unknown branch type: $branch_type")
    end
    
    if isempty(branch.edges)
        error("Branch $(branch.name) has no edges")
    end
    
    total_weight = sum(branch.weights)
    if total_weight <= 0
        return rand(rng, branch.edges)
    end
    
    r = rand(rng) * total_weight
    cumsum_weight = 0.0
    for (i, w) in enumerate(branch.weights)
        cumsum_weight += w
        if r <= cumsum_weight
            return branch.edges[i]
        end
    end
    return last(branch.edges)
end

function _feasible_pair(graph::Vector{Set{Int}}, pair::NTuple{2, Int}, remaining::Vector{Int})
    u, v = pair
    return !_has_edge(graph, u, v) && remaining[u] > 0 && remaining[v] > 0
end

function _discard_candidate!(pairs::Set{NTuple{2, Int}}, triadic_index::Dict{NTuple{2, Int}, Set{Int}},
    candidate_adj::Dict{Int, Set{NTuple{2, Int}}}, pair::NTuple{2, Int})
    
    delete!(pairs, pair)
    if haskey(triadic_index, pair)
        delete!(triadic_index, pair)
    end
    
    u, v = pair
    if haskey(candidate_adj, u)
        delete!(candidate_adj[u], pair)
    end
    if haskey(candidate_adj, v)
        delete!(candidate_adj[v], pair)
    end
end

function _handle_infeasibility!(remaining::Vector{Int}, tolerance::Int, reason::String)
    unmet = residual_edge_count(remaining)
    if unmet <= tolerance
        @warn "$reason. Terminating with $unmet unrealised edges (within tolerance)."
        return true
    else
        throw(InfeasibleState("$reason. $unmet unrealised edges exceeds tolerance $tolerance."))
    end
end

function _validate_remaining_zero(remaining::Vector{Int}, tolerance::Int)
    unmet = residual_edge_count(remaining)
    if unmet > tolerance
        throw(InfeasibleState("Graph generation completed with $unmet unrealised edges, exceeding tolerance $tolerance"))
    end
end

function _generate_graph_attempt(degrees::Vector{Int}, alpha::Float64, beta::Float64,
    gamma::Float64, rng::AbstractRNG, unrealised_edge_tolerance::Int)
    n = length(degrees)

    tree_edges = random_spanning_tree(degrees; rng = rng)
    graph = [Set{Int}() for _ in 1:n]
    for (u, v) in tree_edges
        _add_edge!(graph, u, v)
    end

    tree_degrees = [length(graph[i]) for i in 1:n]
    remaining = degrees .- tree_degrees
    if any(<(0), remaining)
        throw(InfeasibleState("Spanning tree violates degree targets"))
    end

    candidate_pairs, candidate_adj = _initialise_candidate_pairs(graph, remaining)
    triadic_index = _initialise_triadic_index(graph, candidate_pairs)
    distance_cache = _initialise_distance_cache(graph)

    while maximum(remaining) > 0
        if isempty(candidate_pairs)
            _handle_infeasibility!(remaining, unrealised_edge_tolerance,
                "No candidate pairs available while stubs remain") && break
        end

        branch_type = _select_branch_type(alpha, beta, gamma, candidate_pairs, triadic_index, rng)
        if branch_type === nothing
            _handle_infeasibility!(remaining, unrealised_edge_tolerance,
                "No active heuristic branches with candidates") && break
        end

        pair = _sample_edge_from_branch(branch_type, candidate_pairs,
            triadic_index, distance_cache, degrees, rng)

        if _feasible_pair(graph, pair, remaining)
            u, v = pair
            _add_edge!(graph, u, v)
            remaining[u] -= 1
            remaining[v] -= 1

            _update_candidate_pairs!(candidate_pairs, triadic_index, candidate_adj, remaining, u, v)
            _update_triadic_index!(triadic_index, graph, candidate_pairs, u, v)
            _update_distance_cache_after_edge!(distance_cache, u, v)
        else
            _discard_candidate!(candidate_pairs, triadic_index, candidate_adj, pair)
        end
    end

    _validate_remaining_zero(remaining, unrealised_edge_tolerance)
    return graph
end

function _adjacency_list_to_graph(adjacency_list::Vector{Set{Int}})
    n = length(adjacency_list)
    graph = SimpleGraph(n)
    for u in 1:n
        for v in adjacency_list[u]
            if u < v
                add_edge!(graph, u, v)
            end
        end
    end
    return graph
end

"""
    social_configuration_generator(degrees::Vector{Int}; alpha::Float64, beta::Float64,
                                   rng::AbstractRNG = Random.GLOBAL_RNG,
                                   unrealised_edge_tolerance::Int = 10000,
                                   max_restart_attempts::Int = 5)

Generate a graph from a degree sequence using a hybrid heuristic approach with
degree-based, spatial, and triadic closure mechanisms.

# Arguments
- `degrees::Vector{Int}`: Target degree sequence for the graph
- `alpha::Float64`: Weight for degree-based edge selection (must be in [0,1])
- `beta::Float64`: Weight for spatial (distance-based) edge selection (must be in [0,1])
- `rng::AbstractRNG`: Random number generator (default: Random.GLOBAL_RNG)
- `unrealised_edge_tolerance::Int`: Maximum number of unrealised edges allowed (default: 10000)
- `max_restart_attempts::Int`: Number of restart attempts if generation fails (default: 5)

Note: gamma (triadic closure weight) is computed as 1 - alpha - beta, so alpha + beta must be <= 1.

# Returns
- `SimpleGraph`: The generated graph matching the degree sequence

# Throws
- `ArgumentError`: If the degree sequence is invalid or parameters are out of range
- `InfeasibleState`: If graph generation fails after all restart attempts

# Example
```julia
degrees = [3, 3, 2, 2, 2, 2]
graph = social_configuration_generator(degrees; alpha=0.4, beta=0.3)
```
"""
function social_configuration_generator(degrees::Vector{Int}; alpha::Float64, beta::Float64,
    rng::AbstractRNG = Random.GLOBAL_RNG, unrealised_edge_tolerance::Int = 10000,
    max_restart_attempts::Int = 5)
    n = length(degrees)
    if n > 1
        any(<(1), degrees) && error("Degree sequence must give every node degree >= 1 to admit a spanning tree")
        sum(degrees) >= 2 * (n - 1) || error("Degree sequence cannot form a connected graph: need sum >= $(2 * (n - 1)) for $n nodes")
    end
    sum(degrees) % 2 == 0 || error("Degree sequence must sum to an even number")
    degree_sequence_is_graphical(degrees) || error("Degree sequence is not graphical")
    gamma = 1.0 - alpha - beta
    gamma >= 0 || error("alpha + beta must be <= 1")
    unrealised_edge_tolerance >= 0 || error("unrealised_edge_tolerance must be non-negative")
    max_restart_attempts > 0 || error("max_restart_attempts must be positive")

    attempt = 1
    while attempt <= max_restart_attempts
        try
            adjacency_list = _generate_graph_attempt(degrees, alpha, beta, gamma, rng,
                unrealised_edge_tolerance)
            return _adjacency_list_to_graph(adjacency_list)
        catch e
            if e isa InfeasibleState
                if attempt == max_restart_attempts
                    rethrow()
                else
                    @warn "Attempt $attempt failed: $(e.reason). Restarting with a new spanning tree."
                end
            else
                rethrow()
            end
        end
        attempt += 1
    end
    error("Failed to generate graph after $max_restart_attempts attempts")
end

end # module SocialConfiguration
