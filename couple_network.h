extern int** ABC;
extern int* D_ABC;
extern int** B_ABC;
extern string configType;

// --- Function Declarations ---
void generateMultiLayerCoupledNetwork(int*** networks, int** degrees,
    const std::vector<int>& layerOrder,
    double P, int layerCount);

void coupleMultipleNetworks(int** coupledMatrix, int* coupledDegree,
    int*** networks, int** degrees,
    const std::vector<int>& layerOrder, double P,
    int layerCount, int totalNodes);

std::vector<int> determineCouplingOrder(double* k_values, int layerCount, const std::string& configType);
std::vector<int> getLDCConfig(const std::vector<std::pair<double, int>>& layerPairs, int layerCount);
std::vector<int> getMDCConfig(const std::vector<std::pair<double, int>>& layerPairs, int layerCount);
std::vector<int> getHDCConfig(const std::vector<std::pair<double, int>>& layerPairs, int layerCount);
void printCouplingConfig(const std::vector<int>& order, const std::vector<std::pair<double, int>>& layerPairs, const std::string& configType);

/**
 * @brief Sorts nodes by their degree in descending order.
 * @param D Array containing the degrees of nodes.
 * @param N Total number of nodes in a single layer.
 * @return A vector of node indices sorted by degree (highest to lowest).
 */
std::vector<int> sortNodesByDegree(int* D, int N) {
    std::vector<int> indices(N);
    for (int i = 0; i < N; i++) {
        indices[i] = i;
    }
    std::sort(indices.begin(), indices.end(), [&](int a, int b) { return D[a] > D[b]; });
    return indices;
}

/**
 * @brief Implements the assortative coupling between adjacent layers based on the determined order.
 * @param coupledMatrix The adjacency matrix of the fully coupled multilayer network.
 * @param coupledDegree The degree array of the fully coupled multilayer network.
 * @param networks 3D array of intra-layer networks.
 * @param degrees 2D array of intra-layer degrees.
 * @param layerOrder The specific sequence in which layers are coupled.
 * @param P Inter-layer coupling probability (ratio of nodes to be coupled).
 * @param layerCount Total number of layers.
 * @param totalNodes Total number of nodes across all layers.
 */
void coupleMultipleNetworks(int** coupledMatrix, int* coupledDegree,
    int*** networks, int** degrees,
    const std::vector<int>& layerOrder, double P,
    int layerCount, int totalNodes) {

    if (layerOrder.size() < 2) return;

    int total_inter = static_cast<int>(P * N);

    cout << "Chain coupling configuration: ";
    for (int layer : layerOrder) {
        cout << layer << " ";
    }
    cout << endl;

    // Chain coupling: only connect layers that are adjacent in the specified order
    for (int i = 0; i < layerOrder.size() - 1; i++) {
        int layerA = layerOrder[i];
        int layerB = layerOrder[i + 1];

        cout << "Connecting Layer " << layerA << " <-> Layer " << layerB << endl;

        // Assortative coupling requires sorting nodes by degree first
        std::vector<int> sortedA = sortNodesByDegree(degrees[layerA], N);
        std::vector<int> sortedB = sortNodesByDegree(degrees[layerB], N);

        for (int j = 0; j < total_inter; j++) {
            int nodeA = sortedA[j] + layerA * N;
            int nodeB = sortedB[j] + layerB * N;

            if (coupledMatrix[nodeA][nodeB] == 0) {
                coupledMatrix[nodeA][nodeB] = 1;
                coupledMatrix[nodeB][nodeA] = 1;
                coupledDegree[nodeA]++;
                coupledDegree[nodeB]++;
            }
        }
    }    
}

/**
 * @brief Determines the layer coupling order based on the selected configuration strategy.
 * @param k_values Array containing the average degree of each layer.
 * @param layerCount Total number of layers.
 * @param configType Strategy type ("LDC", "MDC", or "HDC").
 * @return A vector representing the sequence of layers for chain coupling.
 */
std::vector<int> determineCouplingOrder(double* k_values, int layerCount, const std::string& configType) {
    std::vector<int> order;

    // Create pairs of (average degree, layer index)
    std::vector<std::pair<double, int>> layerPairs;
    for (int i = 0; i < layerCount; i++) {
        layerPairs.push_back({ k_values[i], i });
    }

    // Sort layers by average degree in ascending order
    std::sort(layerPairs.begin(), layerPairs.end());

    if (configType == "LDC") {
        order = getLDCConfig(layerPairs, layerCount);
    }
    else if (configType == "MDC") {
        order = getMDCConfig(layerPairs, layerCount);
    }
    else if (configType == "HDC") {
        order = getHDCConfig(layerPairs, layerCount);
    }
    else {
        std::cout << "Unknown configuration type, defaulting to LDC configuration." << std::endl;
        order = getLDCConfig(layerPairs, layerCount);
    }

    return order;
}

/**
 * @brief Generates Low-Degree-Core (LDC) configuration.
 * The layer with the minimum degree is placed at the center, and the rest are placed alternating outwards.
 */
std::vector<int> getLDCConfig(const std::vector<std::pair<double, int>>& layerPairs, int layerCount) {
    std::vector<int> order(layerCount, -1);

    int middleIndex = layerCount / 2;
    order[middleIndex] = layerPairs[0].second;  // Minimum degree at the center

    std::vector<int> remaining;
    for (int i = 1; i < layerCount; i++) {
        remaining.push_back(layerPairs[i].second);
    }

    int left = middleIndex - 1;
    int right = middleIndex + 1;

    for (int i = 0; i < remaining.size(); i++) {
        if (i % 2 == 0 && left >= 0) {
            order[left--] = remaining[i];
        }
        else if (right < layerCount) {
            order[right++] = remaining[i];
        }
    }

    return order;
}

/**
 * @brief Generates Medium-Degree-Core (MDC) configuration.
 * Layers are simply arranged in ascending order of their average degrees.
 */
std::vector<int> getMDCConfig(const std::vector<std::pair<double, int>>& layerPairs, int layerCount) {
    std::vector<int> order;

    // Arranged sequentially from lowest to highest average degree
    for (int i = 0; i < layerCount; i++) {
        order.push_back(layerPairs[i].second);
    }

    return order;
}

/**
 * @brief Generates High-Degree-Core (HDC) configuration.
 * The layer with the maximum degree is placed at the center, and the rest are placed alternating outwards in descending order.
 */
std::vector<int> getHDCConfig(const std::vector<std::pair<double, int>>& layerPairs, int layerCount) {
    std::vector<int> order(layerCount, -1);

    int middleIndex = layerCount / 2;
    order[middleIndex] = layerPairs[layerCount - 1].second;  // Maximum degree at the center

    std::vector<int> remaining;
    for (int i = layerCount - 2; i >= 0; i--) {
        remaining.push_back(layerPairs[i].second);
    }

    int left = middleIndex - 1;
    int right = middleIndex + 1;

    for (int i = 0; i < remaining.size(); i++) {
        if (i % 2 == 0 && left >= 0) {
            order[left--] = remaining[i];
        }
        else if (right < layerCount) {
            order[right++] = remaining[i];
        }
    }

    return order;
}

/**
 * @brief Prints the final coupling configuration to the console.
 */
void printCouplingConfig(const std::vector<int>& order, const std::vector<std::pair<double, int>>& layerPairs, const std::string& configType) {
    std::cout << configType << " configuration - Coupling order:";
    for (int layer : order) {
        std::cout << layer << " ";
    }
    std::cout << std::endl;

    std::cout << "Average degree of each layer: ";
    for (const auto& pair : layerPairs) {
        std::cout << "Layer " << pair.second << "(" << pair.first << ") ";
    }
    std::cout << std::endl;
}

/**
 * @brief Master function to assemble the fully coupled multilayer network.
 */
void generateMultiLayerCoupledNetwork(int*** networks, int** degrees,
    const std::vector<int>& layerOrder,
    double P, int layerCount) {

    int totalNodes = layerCount * N;

    // 1. Copy intra-layer connections into the global adjacency matrix
    for (int layer = 0; layer < layerCount; layer++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                int global_i = layer * N + i;
                int global_j = layer * N + j;
                ABC[global_i][global_j] = networks[layer][i][j];
            }
        }
    }

    // 2. Add inter-layer coupling based on the generated layer order
    coupleMultipleNetworks(ABC, D_ABC, networks, degrees, layerOrder, P, layerCount, totalNodes);

    // 3. Build global neighbor lists for the combined network
    for (int i = 0; i < totalNodes; i++) {
        if (B_ABC[i] != nullptr) {
            delete[] B_ABC[i];
            B_ABC[i] = nullptr;
        }

        // Calculate actual degrees for the combined network
        D_ABC[i] = 0;
        for (int j = 0; j < totalNodes; j++) {
            if (ABC[i][j]) D_ABC[i]++;
        }

        // Allocate and populate neighbor array
        B_ABC[i] = new int[D_ABC[i]];
        int count = 0;
        for (int j = 0; j < totalNodes; j++) {
            if (ABC[i][j]) {
                B_ABC[i][count++] = j;
            }
        }
    }
}

/**
 * @brief Outputs detailed debugging information regarding the layer coupling configuration.
 */
void debug_coupling_config(double* k_values, int layerCount, const string& configType) {
    vector<pair<double, int>> layerPairs;
    for (int i = 0; i < layerCount; i++) {
        layerPairs.push_back({ k_values[i], i });
    }
    sort(layerPairs.begin(), layerPairs.end());

    vector<int> order = determineCouplingOrder(k_values, layerCount, configType);

    cout << "=== " << configType << " Configuration Debug Info ===" << endl;
    cout << "Total layers: " << layerCount << endl;
    cout << "Original average degree of each layer: ";
    for (int i = 0; i < layerCount; i++) {
        cout << "Layer " << i << "(k=" << k_values[i] << ") ";
    }
    cout << endl;

    cout << "Sorted by average degree: ";
    for (const auto& p : layerPairs) {
        cout << "Layer " << p.second << "(k=" << p.first << ") ";
    }
    cout << endl;

    cout << "Middle position index: " << (layerCount / 2) << endl;

    cout << "Coupling order (Layer ID): ";
    for (int layer : order) {
        cout << layer << " ";
    }
    cout << endl;

    cout << "Coupling order (Average Degree): ";
    for (int layer : order) {
        for (int i = 0; i < layerCount; i++) {
            if (i == layer) {
                cout << k_values[i] << " ";
                break;
            }
        }
    }
    cout << endl;

    cout << "Position distribution: [";
    for (int i = 0; i < layerCount; i++) {
        cout << order[i];
        if (i < layerCount - 1) cout << ",";
    }
    cout << "]" << endl;

    // Verify if the middle layer is correct based on the strategy
    int middleIndex = layerCount / 2;
    if (configType == "LDC") {
        cout << "Verification: Middle position " << middleIndex << " should be Layer " << layerPairs[0].second
            << " (Minimum degree k=" << layerPairs[0].first << ")" << endl;
    }
    else if (configType == "HDC") {
        cout << "Verification: Middle position " << middleIndex << " should be Layer " << layerPairs[layerCount - 1].second
            << " (Maximum degree k=" << layerPairs[layerCount - 1].first << ")" << endl;
    }
    cout << endl;
}