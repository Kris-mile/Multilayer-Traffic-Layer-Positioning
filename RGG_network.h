// Ensure M_PI is available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Optimized single-layer Random Geometric Graph (RGG) initialization (adapted for int** raw pointers).
 * @param network 2D array representing the adjacency matrix of the network.
 * @param N_node Number of nodes in the layer (explicitly passed to avoid global variable dependency).
 * @param target_avg_degree The desired average degree for the network.
 * @param D Array to store the degree of each node.
 * @param degree_sum Reference to an integer to store the total sum of degrees.
 * @param seed Seed for the random number generator to ensure reproducibility.
 * @param positions Vector to store the 2D spatial coordinates of each node.
 */
inline void initialize_rgg_network_optimized(
    int** network,
    int N_node,
    double target_avg_degree,
    int* D,
    int& degree_sum,
    unsigned int seed,
    vector<pair<double, double>>& positions)
{
    degree_sum = 0;
    std::mt19937 local_rng(seed);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // 1. Clear the network and degree arrays (O(N^2) is necessary for dense matrices, but logic can be optimized)
    // Note: This assumes memory has already been allocated in the main function
    for (int i = 0; i < N_node; i++) {
        // Using memset is faster, but a loop is kept for compatibility
        // Clear to zero for safety
        for (int j = 0; j < N_node; j++) network[i][j] = 0;
        D[i] = 0;
    }

    // 2. Generate node positions
    positions.clear();
    positions.resize(N_node);
    for (int i = 0; i < N_node; i++) {
        positions[i] = { dis(local_rng), dis(local_rng) };
    }

    // 3. Estimate radius (including boundary correction)
    // Correction factor: Nodes near the boundary in a unit square have fewer neighbors, slightly increase radius to compensate
    double rho = N_node;
    double r_est = sqrt(target_avg_degree / (M_PI * rho));
    if (r_est > 0) r_est *= 1.0 + (1.0 / sqrt(N_node)); // Simple boundary compensation
    double r_squared = r_est * r_est;

    // Internal helper function: execute connections and return total degree count
    auto build_connections = [&](double current_r_sq) -> int {
        // Reset connections (to support recalculation)
        for (int i = 0; i < N_node; i++) {
            D[i] = 0;
            // Since recalculation only happens when adjusting the radius, complete clearing here is safe
            for (int j = 0; j < N_node; j++) network[i][j] = 0;
        }

        int current_connections = 0;
        double current_r = sqrt(current_r_sq);

        // Grid division (Cell List algorithm for optimization)
        int grid_size = static_cast<int>(1.0 / current_r);
        if (grid_size < 1) grid_size = 1;
        // Limit maximum number of grids to prevent memory explosion (e.g., when r is extremely small)
        if (grid_size > 200) grid_size = 200;

        vector<vector<int>> grid(grid_size * grid_size);
        double cell_size = 1.0 / grid_size;

        for (int i = 0; i < N_node; i++) {
            int gx = min(static_cast<int>(positions[i].first / cell_size), grid_size - 1);
            int gy = min(static_cast<int>(positions[i].second / cell_size), grid_size - 1);
            grid[gy * grid_size + gx].push_back(i);
        }

        // Iterate through the grid
        for (int gy = 0; gy < grid_size; gy++) {
            for (int gx = 0; gx < grid_size; gx++) {
                int cell_id = gy * grid_size + gx;
                const auto& nodes_in_cell = grid[cell_id];
                if (nodes_in_cell.empty()) continue;

                // Check current and adjacent grids (using symmetry, only check half: current, right, down, bottom-right, bottom-left)
                int neighbor_offsets[5][2] = { {0,0}, {1,0}, {0,1}, {1,1}, {-1,1} };

                for (auto& offset : neighbor_offsets) {
                    int nx = gx + offset[0];
                    int ny = gy + offset[1];

                    if (nx >= 0 && nx < grid_size && ny >= 0 && ny < grid_size) {
                        int n_cell_id = ny * grid_size + nx;
                        const auto& neighbor_nodes = grid[n_cell_id];

                        bool same_cell = (cell_id == n_cell_id);

                        for (size_t i = 0; i < nodes_in_cell.size(); i++) {
                            int u = nodes_in_cell[i];
                            size_t start_j = same_cell ? (i + 1) : 0;

                            for (size_t j = start_j; j < neighbor_nodes.size(); j++) {
                                int v = neighbor_nodes[j];

                                double dx = positions[u].first - positions[v].first;
                                double dy = positions[u].second - positions[v].second;

                                // Simple bounding box pre-check
                                if (abs(dx) > current_r || abs(dy) > current_r) continue;

                                if (dx * dx + dy * dy <= current_r_sq) {
                                    network[u][v] = network[v][u] = 1; // Populate int** array
                                    D[u]++;
                                    D[v]++;
                                    current_connections += 2;
                                }
                            }
                        }
                    }
                }
            }
        }
        return current_connections;
        };

    // 4. Initial connection phase
    degree_sum = build_connections(r_squared);

    // 5. Fast radius adjustment (if deviation > 5%)
    double actual_avg = (double)degree_sum / N_node;
    if (target_avg_degree > 0 && fabs(actual_avg - target_avg_degree) > target_avg_degree * 0.05) {
        double factor = sqrt(target_avg_degree / (actual_avg + 1e-9));
        r_squared *= factor * factor;
        // Reconnect using optimized algorithm instead of brute-force loop
        degree_sum = build_connections(r_squared);
    }
}

/**
 * @brief Multilayer RGG generation entry point (adapted for main function parameters).
 * @param networks 3D array representing the adjacency matrices of all layers.
 * @param degrees 2D array storing the degree of each node in each layer.
 * @param degree_sums Array storing the total sum of degrees for each layer.
 * @param k_values Array specifying the target average degree for each layer.
 * @param layer_count Total number of layers in the multiplex network.
 * @param seed Seed for the random number generator.
 * @param all_positions Pointer to a vector storing the 2D coordinates of nodes across all layers (optional).
 */
inline void generate_multi_layer_RGG(
    int*** networks,
    int** degrees,
    int* degree_sums,
    double* k_values,
    int layer_count,
    unsigned int seed,
    vector<vector<pair<double, double>>>* all_positions = nullptr)
{
    if (all_positions != nullptr) {
        all_positions->resize(layer_count);
    }

#pragma omp parallel for if(layer_count >= 2)
    for (int layer = 0; layer < layer_count; layer++) {
        unsigned int layer_seed = seed + layer * 1000;
        vector<pair<double, double>> positions;

        // Call single-layer generation function
        initialize_rgg_network_optimized(
            networks[layer],
            N,  // Use N defined in main.cpp
            k_values[layer],
            degrees[layer],
            degree_sums[layer],
            layer_seed,
            positions
        );

        if (all_positions != nullptr) {
            (*all_positions)[layer] = positions;
        }
    }
}