/**
 * @file layer_betweenness_analysis.h
 * @brief Analyzes effective betweenness centrality to identify network bottlenecks under varying inter-layer traffic probability (beta/p).
 */

extern double* bI_ABC;      // Intra-layer betweenness
extern double* bE_ABC;      // Inter-layer betweenness

/**
 * @brief Structure to store the maximum node effective betweenness statistics for a layer.
 */
struct MaxNodeBetweenness {
    int layer_id;
    double p;               // Current p (inter-layer probability) value

    // Betweenness values of the node with maximum load
    double max_bI;             // Max intra-layer betweenness (raw value)
    double max_bE;             // Max inter-layer betweenness (raw value)
    int max_bI_node_id;        // Node ID with max intra-layer betweenness
    int max_bE_node_id;        // Node ID with max inter-layer betweenness

    // Effective betweenness of the bottleneck node
    double effective_bI_max;   // Effective intra-layer betweenness = (1-p)/(N-1) * max_bI
    double effective_bE_max;   // Effective inter-layer betweenness = p/[N*(L-1)] * max_bE
    double effective_sum_max;  // Total effective betweenness = effective_bI_max + effective_bE_max
};

/**
 * @brief Calculates the maximum node effective betweenness for a specific layer under a given p value.
 * @param layer_id The index of the layer to analyze.
 * @param p_value The inter-layer routing probability.
 * @return MaxNodeBetweenness structure containing the calculated statistics.
 */
MaxNodeBetweenness calculate_max_node_betweenness(int layer_id, double p_value) {
    MaxNodeBetweenness stats;
    stats.layer_id = layer_id;
    stats.p = p_value;

    // Calculate effective coefficients
    double coeff_internal = (1.0 - p_value) / (N - 1);
    double coeff_external = p_value / (N * (LAYER_COUNT - 1));

    // Initialize statistics
    stats.max_bI = 0.0;
    stats.max_bE = 0.0;
    stats.max_bI_node_id = -1;
    stats.max_bE_node_id = -1;

    // Traverse all nodes in the layer to find the bottleneck
    int start_node = layer_id * N;

    for (int i = 0; i < N; i++) {
        int global_node_id = start_node + i;
        double bI = bI_ABC[global_node_id];
        double bE = bE_ABC[global_node_id];

        // Update max intra-layer betweenness
        if (bI > stats.max_bI) {
            stats.max_bI = bI;
            stats.max_bI_node_id = global_node_id;
        }

        // Update max inter-layer betweenness
        if (bE > stats.max_bE) {
            stats.max_bE = bE;
            stats.max_bE_node_id = global_node_id;
        }
    }

    // Calculate effective betweenness for the bottleneck node
    stats.effective_bI_max = coeff_internal * stats.max_bI;
    stats.effective_bE_max = coeff_external * stats.max_bE;
    stats.effective_sum_max = stats.effective_bI_max + stats.effective_bE_max;

    return stats;
}

/**
 * @brief Analyzes how varying the p value impacts the maximum node effective betweenness across all layers.
 * @param p_start Starting p value.
 * @param p_end Ending p value.
 * @param p_step Increment step for p.
 */
void analyze_p_variation_max_node(double p_start, double p_end, double p_step) {

    cout << "\n=== Analyzing the variation of max node effective betweenness with p ===" << endl;
    cout << "p Range: " << p_start << " -> " << p_end
        << ", Step: " << p_step << endl;
    cout << "Formula: Effective Intra-Betweenness = (1-p)/(N-1) * max_bI" << endl;
    cout << "Formula: Effective Inter-Betweenness = p/[N*(L-1)] * max_bE" << endl;
    cout << "Parameters: N = " << N << ", L = " << LAYER_COUNT << endl;

    // Create output files for each layer
    vector<ofstream> layer_files(LAYER_COUNT);
    vector<string> layer_filenames;

    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        string filename = "betweenness_p_" + to_string(LAYER_COUNT) +
            "layers_layer" + to_string(layer) + ".csv";
        layer_filenames.push_back(filename);
        layer_files[layer].open(filename);

        if (!layer_files[layer]) {
            cerr << "Failed to open file: " << filename << endl;
            return;
        }

        // Write CSV headers
        layer_files[layer] << "p,Max_bI,Max_bE,Max_bI_NodeID,Max_bE_NodeID,"
            << "Effective_bI_max,Effective_bE_max,Effective_Sum_max\n";
    }

    // Create summary comparison file
    string summary_filename = "betweenness_p_" + to_string(LAYER_COUNT) + "layers_summary.csv";
    ofstream summary_file(summary_filename);
    if (!summary_file) {
        cerr << "Failed to open summary file: " << summary_filename << endl;
        return;
    }

    // Summary file headers
    summary_file << "p";
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        summary_file << ",L" << layer << "_Eff_bI_max"
            << ",L" << layer << "_Eff_bE_max"
            << ",L" << layer << "_Eff_Sum_max";
    }
    summary_file << "\n";

    // Iterate through p values
    int step_count = 0;
    int total_steps = static_cast<int>((p_end - p_start) / p_step) + 1;

    for (double p_val = p_start; p_val <= p_end + 1e-9; p_val += p_step) {
        step_count++;

        cout << "\nAnalyzing p = " << fixed << setprecision(3) << p_val
            << " (" << step_count << "/" << total_steps << ")" << endl;

        vector<MaxNodeBetweenness> all_layer_stats;

        summary_file << fixed << setprecision(3) << p_val;

        // Calculate and write statistics for each layer
        for (int layer = 0; layer < LAYER_COUNT; layer++) {
            auto stats = calculate_max_node_betweenness(layer, p_val);
            all_layer_stats.push_back(stats);

            // Write detailed data to layer-specific file
            layer_files[layer] << fixed << setprecision(6)
                << p_val << ","
                << stats.max_bI << ","
                << stats.max_bE << ","
                << stats.max_bI_node_id << ","
                << stats.max_bE_node_id << ","
                << stats.effective_bI_max << ","
                << stats.effective_bE_max << ","
                << stats.effective_sum_max << "\n";

            // Write simplified data to summary file
            summary_file << fixed << setprecision(6)
                << "," << stats.effective_bI_max
                << "," << stats.effective_bE_max
                << "," << stats.effective_sum_max;

            // Display key information on console
            cout << "  Layer " << layer << ": "
                << "Max bI = " << setprecision(4) << stats.max_bI
                << " (Node " << stats.max_bI_node_id << "), "
                << "Max bE = " << stats.max_bE
                << " (Node " << stats.max_bE_node_id << "), "
                << "Eff_bI_max = " << stats.effective_bI_max
                << ", Eff_bE_max = " << stats.effective_bE_max
                << ", Total = " << stats.effective_sum_max << endl;
        }
        summary_file << "\n";
    }

    // Close all files
    for (auto& file : layer_files) {
        file.close();
    }
    summary_file.close();

    cout << "\n=== p variation analysis completed ===" << endl;
    cout << "Files generated:" << endl;
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        cout << "  Layer " << layer << " data: " << layer_filenames[layer] << endl;
    }
    cout << "Summary file: " << summary_filename << endl;

    // Explain file structure
    cout << "\nFile Column Structure:" << endl;
    cout << "Single layer file (8 cols): p, Max_bI, Max_bE, Max_bI_NodeID, Max_bE_NodeID, "
        << "Effective_bI_max, Effective_bE_max, Effective_Sum_max" << endl;
    cout << "Summary file (" << (1 + 3 * LAYER_COUNT) << " cols): p, ";
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        cout << "L" << layer << "_Eff_bI_max, "
            << "L" << layer << "_Eff_bE_max, "
            << "L" << layer << "_Eff_Sum_max";
        if (layer < LAYER_COUNT - 1) cout << ", ";
    }
    cout << endl;
}

/**
 * @brief Quickly retrieves maximum node effective betweenness statistics for all layers at a specific beta.
 * @param p_value The chosen inter-layer transmission probability.
 * @return A vector containing the statistics for each layer.
 */
vector<MaxNodeBetweenness> get_all_layers_max_node_stats(double p_value) {
    vector<MaxNodeBetweenness> all_stats;

    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        all_stats.push_back(calculate_max_node_betweenness(layer, p_value));
    }

    return all_stats;
}