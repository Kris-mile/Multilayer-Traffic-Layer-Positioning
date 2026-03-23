/**
 * @brief Generates a Barab¨¢si-Albert (BA) scale-free network using preferential attachment.
 * @param network 2D array representing the adjacency matrix of the network.
 * @param initial_nodes The number of initial fully connected nodes (m0).
 * @param m The number of edges to attach from a new node to existing nodes.
 * @param D Array to store the degree of each node.
 * @param degree_sum Reference to an integer to store the total sum of degrees.
 * @param seed Seed for the random number generator to ensure reproducibility.
 */
void initialize_ba_network(int** network, int initial_nodes, int m, int D[], int& degree_sum, unsigned int seed) {
    degree_sum = 0;
    std::mt19937 local_rng(seed);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Initialize the network with all zeros
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            network[i][j] = 0;
        }
        D[i] = 0;
    }

    // Step 1: Create an initial fully connected small network
    for (int i = 0; i < initial_nodes; i++) {
        for (int j = i + 1; j < initial_nodes; j++) {
            network[i][j] = network[j][i] = 1;
            D[i]++;
            D[j]++;
            degree_sum += 2;
        }
    }

    // Step 2: Gradually add new nodes
    for (int new_node = initial_nodes; new_node < N; new_node++) {
        // Calculate the degree distribution probability for all current nodes
        vector<double> probabilities(new_node);
        double total_degree = 0.0;

        for (int i = 0; i < new_node; i++) {
            total_degree += D[i];
        }

        // Calculate cumulative probabilities
        double cumulative = 0.0;
        for (int i = 0; i < new_node; i++) {
            cumulative += D[i] / total_degree;
            probabilities[i] = cumulative;
        }

        // Connect the new node to m existing nodes
        int connections_made = 0;
        while (connections_made < m && connections_made < new_node) {
            double rand_val = dis(local_rng);

            // Select a node to connect based on probability
            for (int i = 0; i < new_node; i++) {
                if (rand_val <= probabilities[i]) {
                    // Check if the connection already exists
                    if (network[new_node][i] == 0) {
                        network[new_node][i] = network[i][new_node] = 1;
                        D[new_node]++;
                        D[i]++;
                        degree_sum += 2;
                        connections_made++;
                    }
                    break;
                }
            }
        }
    }
}

/**
 * @brief Calculates the required parameter 'm' to achieve a target average degree in a BA network.
 * @param target_avg_degree The desired average degree for the network.
 * @param initial_nodes The number of initial nodes (default is 5).
 * @return The parameter 'm' (number of edges per new node), ensuring it is at least 1.
 */
int calculate_m_for_target_degree(double target_avg_degree, int initial_nodes = 5) {
    // In a BA network, the average degree is approximately 2m
    int m = static_cast<int>(round(target_avg_degree / 2.0));
    return max(1, m);  // Ensure m is at least 1
}

/**
 * @brief Builds a list of neighbors for each node based on the adjacency matrix.
 * @param network 2D array representing the adjacency matrix.
 * @param neighbor 2D array to store the list of neighbors for each node.
 */
void build_B(int** network, int** neighbor) {
    for (int i = 0; i < N; i++) {
        int count = 0;
        for (int j = 0; j < N; j++) {
            if (network[i][j] != 0) {
                neighbor[i][count++] = j;
            }
        }
    }
}