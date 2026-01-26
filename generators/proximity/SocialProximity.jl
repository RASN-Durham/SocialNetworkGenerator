module SocialProximity

using Graphs
using Random

export social_proximity_generator,  # Entry point for experiments
       synthetic_demographic_similarity, uniform_random_similarity, diag_random_similarity


mutable struct BitGraph
    """
    BitGraph

    A temporary simplified graph structure for the network generator, containing:
    - adjacency matrix (BitMatrix)
    - degree vector
    - number of edges
    - neighbour lists
    - maximum degree
    - number of vertices
    - degree^γ vector
    - γ value
    """
    adj_mat::BitMatrix
    degree::Vector{Int}
    num_edges::Int
    neighbours::Vector{Vector{Int}}
    max_degree::Int
    num_vertices::Int
    degree_gamma::Vector{Float64}
    gamma::Float64
end

# A helper function to initialise a graph and necessary structures
function _initialise_graph(n::Int, gamma::Float64)::BitGraph
    """
    _initialise_graph(n::Int, gamma::Float64)

    Initialise a BitGraph with `n` nodes and parameter `gamma`.

    Parameters:
    -----------
    n :: Int
        Number of nodes to initialise.
    gamma :: Float64
        Exponent used for degree-based weighting.

    Returns:
    --------
    BitGraph
        A new BitGraph instance with an adjacency matrix, degree vector,
        neighbour lists, and gamma settings.
    """
    # Using a Boolean adjacency matrix and degree vector for fast updates
    adj = falses(n, n)
    degree = zeros(Int, n)
    neighbours = Vector{Vector{Int}}(undef, n)

    @inbounds for i in 1:n
        neighbours[i] = Int[]
    end

    degree_gamma = zeros(Float64, n)

    return BitGraph(adj, degree, 0, neighbours, 0, n, degree_gamma, gamma)
end

# Fast edge addition
@inline function _add_edge_bitgraph!(g::BitGraph, u::Int, v::Int)
    """
    _add_edge_bitgraph!(g::BitGraph, u::Int, v::Int)

    Add an edge between nodes `u` and `v` in the given BitGraph, if not already present.
    Updates degrees, neighbour lists, and edge count.

    Parameters:
    -----------
    g :: BitGraph
        The graph to which the edge will be added.
    u :: Int
        One endpoint of the edge.
    v :: Int
        The other endpoint of the edge.

    Returns:
    --------
    None. Modifies the graph `g` in place.
    """
    if !g.adj_mat[u,v]
        @inbounds begin
            g.adj_mat[u,v] = g.adj_mat[v,u] = true
            g.degree[u] += 1
            g.degree[v] += 1
            g.degree_gamma[u] = g.degree[u]^g.gamma
            g.degree_gamma[v] = g.degree[v]^g.gamma
            push!(g.neighbours[u], v)
            push!(g.neighbours[v], u)
            g.num_edges += 1
            if g.degree[u] > g.max_degree
                g.max_degree = g.degree[u]
            end
            if g.degree[v] > g.max_degree
                g.max_degree = g.degree[v]
            end
        end
    else
        error("Edge already exists between $u and $v !")
    end
end

function _bitgraph_to_simplegraph(bg::BitGraph)::SimpleGraph
    """
    _bitgraph_to_simplegraph(bg::BitGraph)

    Convert a BitGraph instance into a Graphs.jl SimpleGraph. Iterates over all
    neighbours to extract undirected edges.

    Parameters:
    -----------
    bg :: BitGraph
        The BitGraph to convert.

    Returns:
    --------
    SimpleGraph
        A Graphs.jl SimpleGraph containing the same nodes and edges as `bg`.
    """
    g = SimpleGraph(bg.num_vertices)
    for u in 1:bg.num_vertices
        sort!(bg.neighbours[u])
    end
    g.fadjlist = bg.neighbours
    g.ne = bg.num_edges

    return g
end

####################
# Similarity Functions
####################

function uniform_random_similarity(n::Int, seed::Int)::Matrix{Float64}
    """
    uniform_random_similarity(n::Int, seed::Int)

    Generate a symmetric similarity matrix for `n` items, based on random positions
    in [0,1). The similarity is computed by reflecting distances around 0.5.

    Parameters:
    -----------
    n :: Int
        The number of items for which to compute similarity.
    seed :: Int
        Seed for random number generation.

    Returns:
    --------
    Matrix{Float64}
        An n × n matrix of similarity values in [0,1].
    """
    Random.seed!(seed)
    loc = rand(n)
    
    sims = abs.(2 .* abs.(loc .- loc') .- 1)

    # Optional checks
    @assert maximum(sims) == 1.0
    @assert minimum(sims) < 0.05
    @assert minimum(sims) > 0.0

    return sims
end

function diag_random_similarity(n::Int, seed=123)::Matrix{Float64}
    """
    diag_random_similarity(n::Int, seed=123)

    Generate a similarity matrix for `n` items such that the diagonal is 1.0
    and off-diagonal values are based on the distance between random points in [0,1).

    Parameters:
    -----------
    n :: Int
        The number of items for which to compute similarity.
    seed :: Int (optional, default=123)
        Seed for random number generation.

    Returns:
    --------
    Matrix{Float64}
        An n × n matrix of similarity values in [0,1].
    """
    Random.seed!(seed)
    loc = rand(n)
    sims = zeros(n, n)
    @inbounds for i in 1:n
        li = loc[i]
        sims[i, i] = 1.0
        @inbounds for j in i+1:n
            ij_dist = abs(li - loc[j])
            ij_perturb_sim = 1 - ij_dist
            sims[i, j] = sims[j, i] = ij_perturb_sim
        end
    end

    @assert maximum(sims) == 1
    @assert minimum(sims) < 0.05
    @assert minimum(sims) > 0
    return sims
end

function synthetic_demographic_similarity(n::Int, num_schools::Int, num_majors::Int, 
    num_genders::Int, num_years::Int, delta::Float64; seed=123)
    """
    synthetic_demographic_similarity(n::Int, num_schools::Int, num_majors::Int,
                                     num_genders::Int, num_years::Int, delta::Float64; seed=123)

    Generate a demographic similarity matrix for `n` items, with various demographic attributes.
    A perturbation matrix (scaled by `delta`) is added to demographic similarity to introduce noise.

    Parameters:
    -----------
    n :: Int
        Number of items (nodes).
    num_schools :: Int
        Number of possible schools.
    num_majors :: Int
        Number of possible majors.
    num_genders :: Int
        Number of possible genders.
    num_years :: Int
        Number of possible years (e.g., class years).
    delta :: Float64
        Scaling factor for the perturbation noise.
    seed :: Int (keyword, default=123)
        Seed for random number generation.

    Returns:
    --------
    Matrix{Float64}
        An n × n matrix representing demographic-based similarity, with added noise.
    """
    Random.seed!(seed)
    school = rand(1:num_schools, n)
    major = rand(1:num_majors, n)
    gender = rand(1:num_genders, n)
    years = rand(1:num_years, n)
    perturbation_loc = rand(n)

    # Demographic similarity
    demo_matrix = Float64.(school .== school') .+ Float64.(major .== major') .+
                  Float64.(gender .== gender') .+ ((num_years - 1 .- abs.(years .- years'))/(num_years - 1))
    demo_matrix ./= maximum(demo_matrix)

    # Perturbation matrix
    diff_matrix = abs.(perturbation_loc .- perturbation_loc')
    perturbation_matrix = abs.(2 .* diff_matrix .- 1)

    demo_matrix .+= delta .* perturbation_matrix
    demo_matrix ./= maximum(demo_matrix)

    return demo_matrix
end

####################
# Structural Similarity
####################

function _structural_similarity(g::BitGraph, u::Int, v::Int)::Float64
    """
    _structural_similarity(g::BitGraph, u::Int, v::Int)

    Compute the structural similarity between two nodes `u` and `v` based on
    preferential attachment (normalized by max degree) and their common neighbours.

    Parameters:
    -----------
    g :: BitGraph
        The graph in which the nodes reside.
    u :: Int
        Index of one node.
    v :: Int
        Index of another node.

    Returns:
    --------
    Float64
        The structural similarity value in [0,1].
    """
    if u == v
        return 0.0
    end
    adj, degree = g.adj_mat, g.degree
    degree_u, degree_v = degree[u], degree[v]
    max_degree = g.max_degree
    pref_attach = degree_v / max_degree
    fof = 0.0
    if degree_u > 0 && degree_v > 0
        cnt = count(adj[:, u] .& adj[:, v])
        fof = cnt / min(degree_u, degree_v)
    end
    return (pref_attach + fof) / 2
end

####################
# Pasta Generator Steps
####################

function _select_node_by_similarity(demographic_similarity::Matrix{Float64}, sim_threshold::Float64, alpha::Float64, g::BitGraph, new_node::Int, prob_cumsum::Vector{Float64})::Int
    """
    _select_node_by_similarity(demographic_similarity::Matrix{Float64}, sim_threshold::Float64,
                           alpha::Float64, g::BitGraph, new_node::Int, prob_cumsum::Vector{Float64})

    For the new node, select an existing node to connect based on a mixture of
    structural and demographic similarity. The final similarity must exceed `sim_threshold`.

    Parameters:
    -----------
    demographic_similarity :: Matrix{Float64}
        Matrix of demographic similarity values.
    sim_threshold :: Float64
        Minimum similarity threshold for a valid edge.
    alpha :: Float64
        Mixture parameter: higher alpha emphasizes structural similarity, lower alpha emphasizes demographic similarity.
    g :: BitGraph
        The existing graph.
    new_node :: Int
        The index of the newly added node.
    prob_cumsum :: Vector{Float64}
        A cumulative sum vector used for sampling.

    Returns:
    --------
    Int
        The index of the selected node, or 0 if no suitable node is found.
    """
    demographic_sim_threshold = alpha < 1 ? ( (sim_threshold - alpha)/(1 - alpha) ) : 1

    @inbounds for i in 1:g.num_vertices
        demographic_sim = demographic_similarity[new_node, i]
        if g.adj_mat[new_node, i] || i == new_node || demographic_sim <= demographic_sim_threshold
            prob_cumsum[i] = i==1 ? 0 : prob_cumsum[i-1]
        else
            val = alpha * _structural_similarity(g, new_node, i) + (1 - alpha) * demographic_sim
            val = val > sim_threshold ? val : 0
            prob_cumsum[i] = i==1 ? val : (prob_cumsum[i-1] + val)
        end
    end

    # Pick a node j to connect using cumulative sampling
    total = prob_cumsum[end]
    if total > 0
        r_j = rand() * total
        j = searchsortedfirst(prob_cumsum, r_j)
        return j
    end
    # cannot find a suitable node
    return 0
end

function _form_triad(g::BitGraph, new_node::Int, selected_node::Int, prob_cumsum::Vector{Float64})::Bool
    """
    _form_triad(g::BitGraph, new_node::Int, selected_node::Int, prob_cumsum::Vector{Float64})

    Attempt to form a triad by connecting the new node with a neighbour of the selected node,
    weighted by the neighbour's degree^gamma. If successful, an edge is added.

    Parameters:
    -----------
    g :: BitGraph
        The graph being constructed.
    new_node :: Int
        The newly added node.
    selected_node :: Int
        The node that was connected during the similarity selection step.
    prob_cumsum :: Vector{Float64}
        A cumulative sum vector used for sampling neighbours.

    Returns:
    --------
    Bool
        True if a triad was successfully formed, False otherwise.
    """
    N_j = g.neighbours[selected_node]
    len_N_j = length(N_j)
    if len_N_j > 2
        # Weighted by degree^gamma
        @inbounds for c_ in 1:len_N_j
            c = N_j[c_]
            if c == new_node || g.adj_mat[new_node, c]
                prob_cumsum[c_] = c_ == 1 ? 0 : prob_cumsum[c_-1]
            else
                val = g.degree_gamma[c]
                prob_cumsum[c_] = c_ == 1 ? val : prob_cumsum[c_-1] + val
            end
        end
        
        @inbounds if prob_cumsum[len_N_j] > 0
            r_k = rand() * prob_cumsum[len_N_j]
            k = searchsortedfirst(prob_cumsum[begin:len_N_j], r_k)
            _add_edge_bitgraph!(g, new_node, N_j[k])
            return true
        end
    end
    return false
end

function _add_triadic_closure_edges(link_count::Int, g::BitGraph)::Int
    """
    _add_triadic_closure_edges(link_count::Int, g::BitGraph)

    Add triadic closure edges among existing nodes, up to `link_count` new edges,
    by trying to complete triangles among neighbours of each node.

    Parameters:
    -----------
    link_count :: Int
        The maximum number of edges to add via triadic closure.
    g :: BitGraph
        The graph to update with triadic closure edges.

    Returns:
    --------
    Int
        The number of triadic closure edges actually added.
    """
    added_tri_links = 0
    tri_link_target = g.num_edges + link_count
    while g.num_edges < tri_link_target
        if !_add_single_triadic_edge(g)
            break
        else
            added_tri_links += 1
        end
    end
    return added_tri_links
end

function _add_single_triadic_edge(g::BitGraph)::Bool
    """
    _add_single_triadic_edge(g::BitGraph)

    Attempt to add an edge by finding two unconnected neighbours of some node `u`,
    thus forming a triangle. The first such opportunity is taken.

    Parameters:
    -----------
    g :: BitGraph
        The graph to modify.

    Returns:
    --------
    Bool
        True if a triadic edge was successfully added, False otherwise.
    """
    @inbounds begin
        node_seq = randperm(g.num_vertices)
        for u in node_seq
            if g.degree[u] < 2
                continue
            end
            N_u = g.neighbours[u]
            shuffle!(N_u)
            len_N_u = length(N_u)
            for w_ in 1:len_N_u -1
                for v_ in (w_+1):len_N_u
                    w = N_u[w_]
                    v = N_u[v_]
                    if !g.adj_mat[v, w]
                        _add_edge_bitgraph!(g, v, w)
                        return true
                    end
                end
            end
        end
    end
    return false
end

####################
# Social Proximity Generator
####################

function social_proximity_generator(n::Int, m_0::Int, m_f::Int, gamma::Float64, alpha::Float64, triad_formation_prob::Float64, 
    triadic_closure_prob::Float64, link_count::Int, sim_threshold::Float64, demo_params; seed::Int=123)
    """
    social_proximity_generator(n::Int, m_0::Int, m_f::Int, gamma::Float64, alpha::Float64,
                    triad_formation_prob::Float64, triadic_closure_prob::Float64,
                    link_count::Int, sim_threshold::Float64,
                    demo_params, seed)

    High-level entry point to generate a network using the social proximity approach.
    The function builds a demographic similarity matrix (either random or
    synthetic) and calls an internal generator to construct the graph.

    Parameters:
    -----------
    n :: Int
        Total number of nodes to generate.
    m_0 :: Int
        Lower bound for the number of edges formed per new node.
    m_f :: Int
        Upper bound for the number of edges formed per new node.
    gamma :: Float64
        Exponent used in degree-based weighting.
    alpha :: Float64
        Mixture parameter for structural vs. demographic similarity.
    triad_formation_prob :: Float64
        Probability to attempt triad formation.
    triadic_closure_prob :: Float64
        Probability to attempt triadic closure edge additions.
    link_count :: Int
        Maximum number of edges to add during triadic closure steps.
    sim_threshold :: Float64
        Minimum combined similarity threshold for establishing new edges.
    demo_params
        Either a string ("diag" or "uniform") or a tuple describing synthetic demographics.
    seed
        Seed for random number generation.

    Returns:
    --------
    SimpleGraph
        The resulting social network graph.
    """
    demographic_similarity = if demo_params === "diag"
        diag_random_similarity(n, seed)
    elseif demo_params === "uniform"
        uniform_random_similarity(n, seed)
    else
        # demo_params = (num_schools, num_majors, num_genders, num_years, perturb)
        num_schools, num_majors, num_genders, num_years, perturb = demo_params
        delta = perturb ? 1/16 : 0.0
        synthetic_demographic_similarity(n, num_schools, num_majors, num_genders, num_years, delta; seed=seed)
    end

    return _social_proximity_generator(n, m_0, m_f, gamma, alpha, triad_formation_prob, triadic_closure_prob, link_count, demographic_similarity, sim_threshold; seed=seed)
end

function _social_proximity_generator(n::Int, m_0::Int, m_f::Int, gamma::Float64, alpha::Float64, triad_formation_prob::Float64, 
    triadic_closure_prob::Float64, link_count::Int, demographic_similarity::Matrix{Float64}, sim_threshold::Float64; seed::Int=1234)
    """
    _social_proximity_generator(n::Int, m_0::Int, m_f::Int, gamma::Float64, alpha::Float64,
                     triad_formation_prob::Float64, triadic_closure_prob::Float64,
                     link_count::Int, demographic_similarity::Matrix{Float64},
                     sim_threshold::Float64; seed=1234)

    Internal routine to generate a network with social proximity mechanisms. Assumes a precomputed
    demographic similarity matrix. For each new node, attempts to form edges via:
    (a) Node selection by combined structural and demographic similarity,
    (b) Triad formation with existing neighbours,
    (c) Triadic closure among existing nodes.

    Parameters:
    -----------
    n :: Int
        Total number of nodes to generate.
    m_0 :: Int
        Lower bound for edges added per new node.
    m_f :: Int
        Upper bound for edges added per new node.
    gamma :: Float64
        Exponent used in degree-based weighting.
    alpha :: Float64
        Mixture parameter (structural vs. demographic).
    triad_formation_prob :: Float64
        Probability to attempt triad formation.
    triadic_closure_prob :: Float64
        Probability to attempt triadic closure edge additions.
    link_count :: Int
        Maximum number of edges to add during triadic closure steps.
    demographic_similarity :: Matrix{Float64}
        Matrix of demographic similarity values.
    sim_threshold :: Float64
        Minimum similarity threshold to form new edges.
    seed :: Int (keyword, default=1234)
        Seed for random number generation.

    Returns:
    --------
    SimpleGraph
        The resulting social network graph.
    """
    @assert n >= 3

    Random.seed!(seed)
    g = _initialise_graph(n, gamma)
    # Initialise with a triangle (1-2, 2-3, 1-3)
    _add_edge_bitgraph!(g, 1, 2)
    _add_edge_bitgraph!(g, 2, 3)
    _add_edge_bitgraph!(g, 1, 3)

    prob_cumsum = zeros(Float64, n)
    a_count, b_count, c_count = 0, 0, 0
    
    for new_node in 4:n
        m = rand(m_0:m_f)
        target_edge_count = g.num_edges + m
        failed_a = false
        failed_c = (triadic_closure_prob == 0)

        while g.num_edges < target_edge_count && !(failed_a && failed_c)
            # Select node by similarity
            failed_a = true
            selected_node = _select_node_by_similarity(demographic_similarity, sim_threshold, alpha, g, new_node, prob_cumsum)
            
            if selected_node > 0
                _add_edge_bitgraph!(g, new_node, selected_node)
                failed_a = false
                a_count += 1

                # Form triad
                @inbounds if rand() < triad_formation_prob
                    b_count += _form_triad(g, new_node, selected_node, prob_cumsum) ? 1 : 0
                end
            end
            # Triadic closure
            if rand() < triadic_closure_prob
                failed_c = true
                added_tri_links = _add_triadic_closure_edges(link_count, g)
                if added_tri_links > 0
                    failed_c = false
                    c_count += added_tri_links
                end
            end
        end
    end
    return _bitgraph_to_simplegraph(g)
end

function _bin_counter(adj_mat::BitMatrix, demo_matrix::Matrix{Float64}; bin_size::Int=10)
    """
    _bin_counter(adj_mat::BitMatrix, demo_matrix::Matrix{Float64})

    Count the number of edges in `adj_mat` that fall into different bins based on `demo_matrix`.

    Parameters:
    -----------
    adj_mat :: BitMatrix
        The adjacency matrix of the graph.
    demo_matrix :: Matrix{Float64}
        The demographic similarity matrix.

    Returns:
    --------
    Tuple{Vector{Int}, Vector{Int}, Vector{Int}}
        (sim_bins, edge_bins, skip_k3_bins) - counts for each bin.
    """
    edge_bins = zeros(Int, bin_size)
    sim_bins = zeros(Int, bin_size)
    skip_k3_bins = zeros(Int, bin_size)
    @inbounds for i in 1:size(adj_mat, 1)
        @inbounds for j in i+1:size(adj_mat, 1)
            bin_idx = Int(floor(demo_matrix[i, j] * bin_size))
            bin_idx += bin_idx < bin_size ? 1 : 0
            if adj_mat[i, j]
                edge_bins[bin_idx] += 1
                if i > 3 || j > 3
                    skip_k3_bins[bin_idx] += 1
                end 
            end
            sim_bins[bin_idx] += 1
        end
    end
    return sim_bins, edge_bins, skip_k3_bins
end

end # module
