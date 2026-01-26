# SocialNetworkGenerator.jl

A Julia package for generating synthetic social networks using different algorithmic approaches. Designed for social scientists and network researchers.

## Overview

This package provides three distinct network generation algorithms, each modeling different aspects of social network formation:

1. **SocialConfiguration** - Generates networks from degree sequences using spatial, degree-based, and triadic closure heuristics
2. **SocialProximity** - Creates networks based on demographic and structural similarity with preferential attachment
3. **SocialCircle** - Builds networks on toroidal grids where connections depend on spatial proximity

## Installation

Currently, this package is structured as separate modules. Load them using `include`:

```julia
include("generators/configuration/SocialConfiguration.jl")
include("generators/proximity/SocialProximity.jl")
include("generators/circle/SocialCircle.jl")

using .SocialConfiguration
using .SocialProximity
using .SocialCircle
```

*Note: A unified package structure with proper Project.toml integration is planned for future releases.*

## Quick Start

### SocialConfiguration

Generate a graph from a degree sequence using hybrid heuristics:

```julia
include("generators/configuration/SocialConfiguration.jl")
using .SocialConfiguration
using Random

degrees = [3, 3, 2, 2, 2, 2]
graph = social_configuration_generator(degrees; alpha=0.4, beta=0.3, rng=MersenneTwister(42))
```

**Parameters:**
- `degrees`: Target degree sequence
- `alpha`: Weight for degree-based edge selection [0,1]
- `beta`: Weight for spatial (distance) edge selection [0,1]
- `gamma` (computed as 1-alpha-beta): Weight for triadic closure

**Useful utilities:**
- `degree_sequence_is_graphical(degrees)` - Validate degree sequences
- `random_spanning_tree(degrees)` - Generate random spanning trees

### SocialProximity

Create networks with demographic similarity and triadic mechanisms:

```julia
include("generators/proximity/SocialProximity.jl")
using .SocialProximity

# With random demographic similarity
graph = social_proximity_generator(
    100,        # n: number of nodes
    2,          # m_0: min edges per node
    5,          # m_f: max edges per node
    0.5,        # gamma: degree exponent
    0.7,        # alpha: structural vs demographic mix
    0.3,        # triad_formation_prob
    0.2,        # triadic_closure_prob
    3,          # link_count: max triadic closure edges
    0.5,        # sim_threshold: minimum similarity
    "diag";     # demo_params: "diag", "uniform", or tuple
    seed=42
)

# With synthetic demographics
demo_params = (
    2,      # num_schools
    5,      # num_majors
    2,      # num_genders
    4,      # num_years
    true    # perturb: add noise
)
graph = social_proximity_generator(100, 2, 5, 0.5, 0.7, 0.3, 0.2, 3, 0.5, demo_params)
```

**Useful utilities:**
- `synthetic_demographic_similarity(n, num_schools, num_majors, num_genders, num_years, delta)` 
- `uniform_random_similarity(n, seed)`
- `diag_random_similarity(n, seed)`

### SocialCircle

Build spatial networks on toroidal grids:

```julia
include("generators/circle/SocialCircle.jl")
using .SocialCircle

# Simple version: same distance for all nodes
graph = social_circle_generator(100, 100, 500, fill(10.0, 500))

# With node types having different reach distances
type_ratios = [0.3, 0.7]  # 30% type 1, 70% type 2
type_distances = [5.0, 10.0]  # Type 1: distance 5, Type 2: distance 10
graph = social_circle_generator(100, 100, 500, type_ratios, type_distances)
```

**Useful utilities:**
- `squared_euclidean_distance(point1, point2, rows, cols)` - Toroidal distance calculation

## Package Structure

```
SocialNetworkGenerator/
├── Project.toml
├── LICENSE
├── README.md
├── examples.jl
└── generators/
    ├── configuration/
    │   └── SocialConfiguration.jl
    ├── proximity/
    │   └── SocialProximity.jl
    └── circle/
        └── SocialCircle.jl
```

## Naming Conventions

This package follows Julia best practices:

- **Function names**: `snake_case` with `!` for mutating functions
- **Type names**: `PascalCase`
- **Parameters**: Spelled-out words (`alpha` not `α`, `gamma` not `γ`)
- **Private functions**: Prefixed with `_` (e.g., `_add_edge!`)
- **British spelling**: `initialise` not `initialize`, `neighbour` not `neighbor`
- **Graph terminology**: Standard graph theory terms (node, edge, degree, adjacency)

## Return Types

All three main generator functions return `SimpleGraph` from Graphs.jl for consistency:

- `social_configuration_generator` → `SimpleGraph`
- `social_proximity_generator` → `SimpleGraph`
- `social_circle_generator` → `SimpleGraph`

## Contributing

This package is designed as an example for social scientists. Contributions welcome!

## Acknowledgments

This work is supported by the Medical Research Council (UK) [MR/W02974X/1].

## License

MIT License - see LICENSE file for details.

## Citation

If you use this package in your research, please cite:

```bibtex
@software{SocialNetworkGenerator2025,
  author = {Zhu, Yiran and Badham, Jennifer and Chueca Del Cerro, Cristina and Kutner, David},
  title = {SocialNetworkGenerator.jl: Synthetic Social Network Generation},
  year = {2025},
  url = {https://github.com/RASN-Durham/SocialNetworkGenerator},
  note = {Supported by the Medical Research Council (UK) [MR/W02974X/1]}
}
```
