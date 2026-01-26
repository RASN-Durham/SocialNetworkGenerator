module SocialCircle

using Graphs
using Random
using Statistics
using StatsBase
using Distributions

export social_circle_generator,
       squared_euclidean_distance

"""
    squared_euclidean_distance(point1::Tuple{Int, Int}, point2::Tuple{Int, Int}, rows::Int, cols::Int)

Efficiently computes the wrapped (toroidal) squared Euclidean distance between two points on a toroidal grid.

# Arguments
- `point1::Tuple{Int, Int}`: First point coordinates (x, y)
- `point2::Tuple{Int, Int}`: Second point coordinates (x, y)
- `rows::Int`: Number of rows in the toroidal grid
- `cols::Int`: Number of columns in the toroidal grid

# Returns
- `Float64`: The squared Euclidean distance considering toroidal wrap-around
"""
@fastmath @inline function squared_euclidean_distance(point1::Tuple{Int, Int}, point2::Tuple{Int, Int}, rows::Int, cols::Int)::Float64
    dx = abs(point1[1] - point2[1])
    dy = abs(point1[2] - point2[2])
    dx = min(dx, cols - dx)
    dy = min(dy, rows - dy)

    return dx^2 + dy^2
end

@inline function _is_reachable(node1::Tuple{Int, Int}, node2::Tuple{Int, Int}, distance::Float64, rows::Int, cols::Int)::Bool
    return squared_euclidean_distance(node1, node2, rows, cols) ≤ distance^2
end

"""
    social_circle_generator(rows::Int, cols::Int, num_nodes::Int, node_distances::Vector{<:Number})

Creates a social network graph on a toroidal grid where nodes are placed randomly and 
connections depend on spatial distance.

# Arguments
- `rows::Int`: Number of rows in the toroidal grid
- `cols::Int`: Number of columns in the toroidal grid
- `num_nodes::Int`: Total number of nodes in the network
- `node_distances::Vector{Number}`: Social reach distance for each node (must have length num_nodes)

# Returns
- `SimpleGraph`: The resulting social network as a SimpleGraph

# Example
```julia
using SocialCircle
graph = social_circle_generator(100, 100, 500, fill(10.0, 500))
```
"""
function social_circle_generator(rows::Int, cols::Int, num_nodes::Int, node_distances::Vector{<:Number})::SimpleGraph
    @assert length(node_distances) == num_nodes "The length of node_distances must equal num_nodes"
    
    # Initialise graph
    graph = SimpleGraph(num_nodes)

    # Sample unique grid positions for nodes
    grid_size = rows * cols
    selected_indices = sample(1:grid_size, num_nodes; replace=false)
    
    """
    Compute node position from its index in the selected_indices array
    """
    @inline function _get_node_position(node_idx::Int)::Tuple{Int,Int}
        idx = selected_indices[node_idx]
        return (mod1(idx, cols), div(idx - 1, cols) + 1)
    end
    
    # Create node positions
    node_positions = map(_get_node_position, 1:num_nodes)

    # Create edges based on distance and node reach
    @inbounds for node1 in 1:(num_nodes-1)
        pos1 = node_positions[node1]
        reach1 = node_distances[node1]
        
        # Process potential connections for node1
        for node2 in (node1 + 1):num_nodes
            max_reach = min(reach1, node_distances[node2])
            if _is_reachable(pos1, node_positions[node2], max_reach, rows, cols)
                add_edge!(graph, node1, node2)
            end
        end
    end

    return graph
end

"""
    social_circle_generator(rows::Int, cols::Int, num_nodes::Int, 
                          type_ratios::Vector{<:Number}, 
                          type_distances::Union{Vector{<:Number}, Vector{<:UnivariateDistribution}})

Creates a social network graph with nodes of different types, where types have different social reach distances.

# Arguments
- `rows::Int`: Number of rows in the toroidal grid
- `cols::Int`: Number of columns in the toroidal grid
- `num_nodes::Int`: Total number of nodes in the network
- `type_ratios::Vector{Number}`: Proportion of nodes for each type (should sum to 1.0)
- `type_distances::Vector{Number|Distribution}`: Social reach distance for each node type (can be fixed values or distributions)

# Returns
- `SimpleGraph`: The resulting social network as a SimpleGraph

# Example
```julia
using SocialCircle, Distributions
# 30% type 1 (distance 5), 70% type 2 (distance from Normal(10, 2))
graph = social_circle_generator(100, 100, 500, [0.3, 0.7], [5.0, Normal(10, 2)])
```
"""
function social_circle_generator(rows::Int, cols::Int, num_nodes::Int, type_ratios::Vector{<:Number}, type_distances::Union{Vector{<:Number}, Vector{<:UnivariateDistribution}})::SimpleGraph
    @assert length(type_ratios) == length(type_distances) "The number of node types and distances must match"
    @assert isapprox(sum(type_ratios), 1.0) "Type ratios must sum to 1.0"
    
    # Prepare node_distances vector based on type assignments
    node_distances = Vector{Float64}(undef, num_nodes)
    cum_ratios = cumsum(type_ratios)
    
    # Pre-calculate the node indices where type changes
    type_boundaries = [round(Int, r * num_nodes) for r in cum_ratios]
    
    # Only update nodes that are of a different type
    start_idx = 1
    for (type_idx, end_idx) in enumerate(type_boundaries)
        if end_idx <= start_idx
            continue  # Skip if no nodes of this type
        end
        # Check if this type's distance is a distribution and sample from it if needed
        type_dist = type_distances[type_idx]
        if isa(type_dist, UnivariateDistribution)
            # Sample a distance for each node of this type from the distribution
            node_distances[start_idx:end_idx] .= float.(rand(type_dist, end_idx - start_idx + 1))
        else
            # Use the fixed distance value for all nodes of this type
            node_distances[start_idx:end_idx] .= float.(type_distances[type_idx])
        end
        start_idx = end_idx + 1
    end
    
    # Call the main implementation with the calculated node distances
    return social_circle_generator(rows, cols, num_nodes, node_distances)
end

end # module SocialCircle
