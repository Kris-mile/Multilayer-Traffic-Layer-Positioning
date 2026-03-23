/**
 * @brief Initializes an ER network using a target average degree instead of connection probability.
 * @param network 2D array representing the adjacency matrix of the network.
 * @param average_degree The target average degree for the network.
 * @param D Array to store the degree of each node.
 * @param degree_sum Reference to an integer to store the total sum of degrees.
 * @param seed Seed for the random number generator to ensure reproducibility.
 */
void initialize_er_network(int** network, double average_degree, int D[], int& degree_sum, unsigned int seed) {
    degree_sum = 0;
    std::mt19937 local_rng(seed);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Calculate the corresponding connection probability
    double p = average_degree / (N - 1);

    // Initialize all matrix elements and degree array to 0
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            network[i][j] = 0;
        }
        D[i] = 0;
    }

    // Randomly connect nodes based on the calculated probability
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            if (dis(local_rng) < p) {
                network[i][j] = network[j][i] = 1;
                D[i]++;
                D[j]++;
                degree_sum += 2;
            }
        }
    }
}

