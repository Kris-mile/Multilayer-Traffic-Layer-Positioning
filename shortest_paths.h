/**
 * @file shortest_paths.h
 * @brief Calculates and stores the shortest paths for both individual layers and the global coupled network using BFS.
 */

extern int** ABC;
extern int* D_ABC;
extern int** B_ABC;

extern int*** path_ABC;
extern int** path_num_ABC;
extern int distance_ABC[TOTAL_NODES][TOTAL_NODES];

extern int**** path_single;    // 4D Array: [layer][i][j][k]
extern int*** path_num_single; // Number of shortest paths in single layer [layer][i][j]
extern double coupling_prob;

/**
 * @brief Initializes memory structures for storing global coupled network paths.
 */
void initialize_ABC_path_storage() {
    for (int i = 0; i < TOTAL_NODES; i++) {
        path_num_ABC[i] = new int[TOTAL_NODES];
        path_ABC[i] = new int* [TOTAL_NODES];

        for (int j = 0; j < TOTAL_NODES; j++) {
            path_num_ABC[i][j] = 0;
            path_ABC[i][j] = nullptr;
            distance_ABC[i][j] = -1;
        }
    }
}

/**
 * @brief Frees the memory allocated for the global coupled network paths.
 */
void free_ABC_path_storage() {
    for (int i = 0; i < TOTAL_NODES; i++) {
        for (int j = 0; j < TOTAL_NODES; j++) {
            if (path_ABC[i][j] != nullptr) {
                delete[] path_ABC[i][j];
            }
        }
        delete[] path_num_ABC[i];
        delete[] path_ABC[i];
    }
    delete[] path_num_ABC;
    delete[] path_ABC;
}

/**
 * @brief Calculates all-pairs shortest paths for the multiplex coupled network using BFS.
 */
void multiplex_network_shortest_path_ABC() {
    initialize_ABC_path_storage();

    // Initialize distance matrix
    for (int i = 0; i < TOTAL_NODES; i++) {
        for (int j = 0; j < TOTAL_NODES; j++) {
            if (i == j) distance_ABC[i][j] = 0;
            else distance_ABC[i][j] = (ABC[i][j] ? 1 : -1);
        }
    }

    int allow_vertices[TOTAL_NODES]; // Array to mark if a node has been processed

    // Outer loop: iterate through each source node
    for (int node = 0; node < TOTAL_NODES; node++) {
        int next = node;

        // Distance to itself is 0
        distance_ABC[node][node] = 0;

        // Initialize marking array
        for (int i = 0; i < TOTAL_NODES; i++) {
            allow_vertices[i] = 0;
        }

        // Inner loop: gradually expand the shortest path set
        for (int allow_size = 0; allow_size < TOTAL_NODES; allow_size++) {
            int maximum = INT_MAX;
            next = -1;

            // Find the closest unvisited node
            for (int i = 0; i < TOTAL_NODES; i++) {
                // Find node with a valid path that hasn't been added
                if (distance_ABC[node][i] != -1 && allow_vertices[i] == 0) {
                    if (distance_ABC[node][i] < maximum) {
                        maximum = distance_ABC[node][i];
                        next = i; // Update to nearest node
                    }
                }
            }

            if (next == -1) break;
            allow_vertices[next] = 1;

            // Update distances of neighboring nodes
            for (int i = 0; i < D_ABC[next]; i++) {
                int neighbor = B_ABC[next][i];

                // Skip nodes already processed
                if (allow_vertices[neighbor] == 0) {
                    // New distance = distance from source to current node + 1
                    int new_distance = distance_ABC[node][next] + 1;

                    if (distance_ABC[node][neighbor] == -1 || new_distance < distance_ABC[node][neighbor]) {
                        distance_ABC[node][neighbor] = new_distance;
                    }
                }
            }
        }
    }

    // Path recording phase
    for (int i = 0; i < TOTAL_NODES; i++) {
        for (int j = 0; j < TOTAL_NODES; j++) {

            if (i == j || distance_ABC[i][j] == -1) {
                path_num_ABC[i][j] = 0;
                path_ABC[i][j] = nullptr;
                continue;
            }

            // Calculate path count
            path_num_ABC[i][j] = 0;
            for (int k = 0; k < D_ABC[i]; k++) {
                if (distance_ABC[B_ABC[i][k]][j] == distance_ABC[i][j] - 1) {
                    path_num_ABC[i][j]++;
                }
            }

            // Allocate memory and store paths
            if (path_num_ABC[i][j] > 0) {
                // Store only the next-hop node for each path
                path_ABC[i][j] = new int[path_num_ABC[i][j]];

                // Initialize path array
                for (int k = 0; k < path_num_ABC[i][j]; k++)
                    path_ABC[i][j][k] = -1; // Initialize to invalid value

                int count = 0;
                for (int k = 0; k < D_ABC[i]; k++) {
                    int neighbor = B_ABC[i][k];
                    if (distance_ABC[neighbor][j] == distance_ABC[i][j] - 1) {
                        if (count < path_num_ABC[i][j]) {
                            path_ABC[i][j][count] = neighbor;
                            count++;
                        }
                    }
                }

                // Verify the actual filled count
                if (count != path_num_ABC[i][j]) {
                    cerr << "ERROR: Expected " << path_num_ABC[i][j] << " paths, but filled " << count << endl;
                }
            }
            else {
                path_ABC[i][j] = nullptr;
            }
        }
    }

    // Validate paths
    int invalid_paths = 0;
    for (int i = 0; i < TOTAL_NODES; i++) {
        for (int j = 0; j < TOTAL_NODES; j++) {
            if (path_num_ABC[i][j] > 0) {
                if (path_ABC[i][j] == nullptr) {
                    cerr << "Error: Node " << i << " -> " << j
                        << " path count > 0 but array is null!" << endl;
                    invalid_paths++;
                }
                else {
                    for (int k = 0; k < path_num_ABC[i][j]; k++) {
                        int neighbor = path_ABC[i][j][k];
                        if (neighbor < 0 || neighbor >= TOTAL_NODES) {
                            cerr << "Error: Invalid neighbor node " << neighbor
                                << " (Range: 0-" << (TOTAL_NODES - 1) << ")" << endl;
                            invalid_paths++;
                        }
                    }
                }
            }
        }
    }

    int unreachable_count = 0;
    for (int i = 0; i < TOTAL_NODES; i++) {
        for (int j = 0; j < TOTAL_NODES; j++) {
            if (distance_ABC[i][j] == -1) {
                unreachable_count++;
            }
        }
    }
    //cout << "Total unreachable pairs: " << unreachable_count << endl;
    //cout << "Found " << invalid_paths << " invalid paths" << endl;   
}

/**
 * @brief Initializes memory structures for single-layer shortest paths.
 */
void initialize_single_path_storage() {
    path_single = new int*** [LAYER_COUNT];
    path_num_single = new int** [LAYER_COUNT];

    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        path_single[layer] = new int** [N];
        path_num_single[layer] = new int* [N];

        for (int i = 0; i < N; i++) {
            path_single[layer][i] = new int* [N];
            path_num_single[layer][i] = new int[N];

            for (int j = 0; j < N; j++) {
                path_single[layer][i][j] = nullptr; // Initialize to null pointer
                path_num_single[layer][i][j] = 0;   // Initial path count is 0
            }
        }
    }
}

/**
 * @brief Frees the memory allocated for single-layer paths.
 */
void free_single_path_storage() {
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (path_single[layer][i][j] != nullptr) {
                    delete[] path_single[layer][i][j];
                }
            }
            delete[] path_single[layer][i];
            delete[] path_num_single[layer][i];
        }
        delete[] path_single[layer];
        delete[] path_num_single[layer];
    }
    delete[] path_single;
    delete[] path_num_single;
}

/**
 * @brief Calculates shortest paths for a specified single layer network using BFS.
 * @param layer_index The index of the layer to calculate paths for.
 * @param network The adjacency matrix of the specific layer.
 * @param D The degree array of the specific layer.
 * @param B The adjacency list of the specific layer.
 */
void calculate_single_layer_paths(int layer_index, int** network, int* D, int** B) {
    // Create temporary distance matrix for the current layer
    int** distance = new int* [N];
    for (int i = 0; i < N; i++) {
        distance[i] = new int[N];
        for (int j = 0; j < N; j++) {
            if (i == j) distance[i][j] = 0;
            else distance[i][j] = (network[i][j] ? 1 : -1);
        }
    }

    int allow_vertices[TOTAL_NODES]; // Mark if node is processed

    // Outer loop: iterate through each source node
    for (int node = 0; node < N; node++) {
        // Distance to itself is 0
        distance[node][node] = 0;

        // Initialize marker array
        for (int i = 0; i < N; i++) {
            allow_vertices[i] = 0;
        }

        // Inner loop: gradually expand shortest path set
        for (int allow_size = 0; allow_size < N; allow_size++) {
            int maximum = INT_MAX;
            int next = -1;

            // Find closest unprocessed node and update neighbors
            for (int i = 0; i < N; i++) {
                // Find closer, unvisited node
                if (distance[node][i] != -1 && allow_vertices[i] == 0) {
                    if (distance[node][i] < maximum) {
                        maximum = distance[node][i];
                        next = i;
                    }
                }
            }

            if (next == -1) break;
            allow_vertices[next] = 1;

            // Update distances for neighbors
            for (int i = 0; i < D[next]; i++) {
                int neighbor = B[next][i];

                if (allow_vertices[neighbor] == 0) {
                    int new_distance = distance[node][next] + 1;

                    if (distance[node][neighbor] == -1 || new_distance < distance[node][neighbor]) {
                        distance[node][neighbor] = new_distance;
                    }
                }
            }
        }
    }

    // Path recording
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j || distance[i][j] == -1) {
                path_num_single[layer_index][i][j] = 0;
                path_single[layer_index][i][j] = nullptr;
                continue;
            }

            path_num_single[layer_index][i][j] = 0;
            for (int k = 0; k < D[i]; k++) {
                if (distance[B[i][k]][j] == distance[i][j] - 1) {
                    path_num_single[layer_index][i][j]++;
                }
            }

            if (path_num_single[layer_index][i][j] > 0) {
                path_single[layer_index][i][j] = new int[path_num_single[layer_index][i][j]];

                for (int k = 0; k < path_num_single[layer_index][i][j]; k++)
                    path_single[layer_index][i][j][k] = -1;

                int count = 0;
                for (int k = 0; k < D[i]; k++) {
                    int neighbor = B[i][k];
                    if (distance[neighbor][j] == distance[i][j] - 1) {
                        if (count < path_num_single[layer_index][i][j]) {
                            path_single[layer_index][i][j][count] = neighbor;
                            count++;
                        }
                    }
                }
            }
            else {
                path_single[layer_index][i][j] = nullptr;
            }
        }
    }

    // Free temporary distance matrix
    for (int i = 0; i < N; i++) {
        delete[] distance[i];
    }
    delete[] distance;
}