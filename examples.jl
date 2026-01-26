# Example Usage of SocialNetworkGenerator

# This example demonstrates the consistent interface across all three generators

using Pkg
# Pkg.develop(path="path/to/SocialNetworkGenerator")

# Load the three generator modules
include("generators/configuration/SocialConfiguration.jl")
include("generators/proximity/SocialProximity.jl")
include("generators/circle/SocialCircle.jl")

using .SocialConfiguration
using .SocialProximity
using .SocialCircle
using Graphs
using Random
using Statistics

println("=== SocialNetworkGenerator Examples ===\n")

## Example 1: SocialConfiguration Generator
println("1. SocialConfiguration Generator")
println("-" ^ 50)

# Define a degree sequence
degrees = [3, 3, 2, 2, 2, 2]

# Check if it's graphical
if degree_sequence_is_graphical(degrees)
    println("Degree sequence is graphical: ", degrees)
    
    # Generate graph with 40% degree-based, 30% spatial, 30% triadic
    graph_config = social_configuration_generator(
        degrees; 
        alpha=0.4,      # degree-based weight
        beta=0.3,       # spatial weight (gamma=0.3 computed automatically)
        rng=MersenneTwister(42)
    )
    
    println("Generated graph: ", nv(graph_config), " nodes, ", ne(graph_config), " edges")
    println("Degree sequence: ", sort(degree(graph_config), rev=true))
    println()
else
    println("Degree sequence is not graphical!")
end

## Example 2: SocialProximity Generator
println("2. SocialProximity Generator")
println("-" ^ 50)

# Generate with random demographic similarity
graph_prox = social_proximity_generator(
    50,         # num_nodes
    2,          # m_0: min edges per node
    4,          # m_f: max edges per node  
    0.5,        # gamma: degree exponent
    0.7,        # alpha: structural vs demographic mix
    0.3,        # triad_formation_prob
    0.2,        # triadic_closure_prob
    3,          # link_count for triadic closure
    0.5,        # sim_threshold
    "diag";     # demo_params
    seed=42
)

println("Generated graph: ", nv(graph_prox), " nodes, ", ne(graph_prox), " edges")
println()

# Alternative: synthetic demographics
demo_params = (
    2,      # num_schools
    3,      # num_majors
    2,      # num_genders
    4,      # num_years
    true    # add perturbation
)

graph_prox2 = social_proximity_generator(
    30, 1, 3, 0.5, 0.6, 0.4, 0.1, 2, 0.6, demo_params; seed=123
)

println("Synthetic demographics graph: ", nv(graph_prox2), " nodes, ", ne(graph_prox2), " edges")
println()

## Example 3: SocialCircle Generator
println("3. SocialCircle Generator")
println("-" ^ 50)

# Simple version: same reach for all nodes
graph_circle = social_circle_generator(
    50,                 # rows
    50,                 # cols
    200,                # num_nodes
    fill(8.0, 200)     # all nodes have reach distance 8
)

println("Generated graph: ", nv(graph_circle), " nodes, ", ne(graph_circle), " edges")
println("Average degree: ", round(mean(degree(graph_circle)), digits=2))
println()

# Version with node types
type_ratios = [0.4, 0.6]                    # 40% type 1, 60% type 2
type_distances = [5.0, 10.0]                # Type 1: distance 5, Type 2: distance 10

graph_circle2 = social_circle_generator(
    50, 50, 200, type_ratios, type_distances
)

println("Typed nodes graph: ", nv(graph_circle2), " nodes, ", ne(graph_circle2), " edges")
println()

## Demonstrating Utility Functions
println("4. Utility Functions")
println("-" ^ 50)

# From SocialConfiguration
test_degrees = [4, 4, 3, 3, 2, 2]
println("Is [4,4,3,3,2,2] graphical? ", degree_sequence_is_graphical(test_degrees))

tree_edges = random_spanning_tree(test_degrees)
println("Random spanning tree has ", length(tree_edges), " edges")

# From SocialProximity
sim_matrix = uniform_random_similarity(10, 42)
println("Similarity matrix shape: ", size(sim_matrix))

# From SocialCircle  
point1 = (10, 15)
point2 = (45, 48)
dist = squared_euclidean_distance(point1, point2, 50, 50)
println("Toroidal distance² between ", point1, " and ", point2, ": ", dist)

println("\n=== All examples completed successfully! ===")
