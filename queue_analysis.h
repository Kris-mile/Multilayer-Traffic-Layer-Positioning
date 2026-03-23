/**
 * @file queue_analysis.h
 * @brief Provides theoretical analysis for queue lengths and the critical packet generation rate (r_c).
 */

extern double* bI_ABC;  // Intra-layer betweenness
extern double* bE_ABC;  // Inter-layer betweenness
extern double p;        // Inter-layer transmission probability

/**
 * @brief Calculates the theoretical queue length for each node given a packet generation rate.
 * @param r The packet generation rate.
 * @return Array containing the theoretical queue length for every node.
 */
double* calculate_theoretical_queue_length(double r) {
    double* q_i = new double[TOTAL_NODES]();

    for (int i = 0; i < TOTAL_NODES; i++) {
        double q_I = (r * (1 - p) / (N - 1)) * bI_ABC[i];
        double q_E = (r * p / (N * (LAYER_COUNT - 1))) * bE_ABC[i];
        q_i[i] = q_I + q_E;
    }

    return q_i;
}

/**
 * @brief Calculates the theoretical critical packet generation rate (r_c) for the network.
 * @return The critical rate r_c.
 */
double calculate_critical_rate() {
    double min_r_c = numeric_limits<double>::max();
    int bottleneck_node = -1;

    for (int i = 0; i < TOTAL_NODES; i++) {
        double denominator = ((1 - p) / (N - 1)) * bI_ABC[i] +
            (p / (N * (LAYER_COUNT - 1))) * bE_ABC[i];

        if (denominator > 0) {
            double r_c_i = C_ave / denominator;
            if (r_c_i < min_r_c) {
                min_r_c = r_c_i;
                bottleneck_node = i;
            }
        }
    }

    cout << "Theoretical critical packet generation rate r_c = " << min_r_c << endl;
    if (bottleneck_node != -1) {
        cout << "Bottleneck node: " << bottleneck_node
            << " (Layer " << bottleneck_node / N << ")" << endl;
    }

    return min_r_c;
}

/**
 * @brief Outputs the theoretical queue analysis results to a CSV file.
 * @param r The packet generation rate.
 * @param filename Name of the output file.
 */
void output_queue_analysis(double r, const string& filename) {
    ofstream out_file(filename);
    if (!out_file) {
        cerr << "Failed to open file: " << filename << endl;
        return;
    }

    out_file << "Node_ID,Layer,Internal_Betweenness,External_Betweenness,Queue_Length\n";
    double* q_i = calculate_theoretical_queue_length(r);

    for (int i = 0; i < TOTAL_NODES; i++) {
        int layer = i / N;
        out_file << i << "," << layer << ","
            << bI_ABC[i] << "," << bE_ABC[i] << ","
            << q_i[i] << "\n";
    }

    delete[] q_i;
    out_file.close();
    cout << "Queue analysis results saved to: " << filename << endl;
}

/**
 * @brief Outputs the critical rate (r_c) analysis for every node to a CSV file.
 * @param filename Name of the output file.
 */
void output_critical_analysis(const string& filename) {
    ofstream out_file(filename);
    if (!out_file) {
        cerr << "Failed to open file: " << filename << endl;
        return;
    }

    out_file << "Node_ID,Layer,Internal_Betweenness,External_Betweenness,Critical_Rate\n";

    double global_min_r_c = numeric_limits<double>::max();
    int global_bottleneck_node = -1;

    for (int i = 0; i < TOTAL_NODES; i++) {
        int layer = i / N;
        double denominator = ((1 - p) / (N - 1)) * bI_ABC[i] +
            (p / (N * (LAYER_COUNT - 1))) * bE_ABC[i];

        double r_c_i = (denominator > 0) ? C_ave / denominator : 0;

        out_file << i << "," << layer << ","
            << bI_ABC[i] << "," << bE_ABC[i] << ","
            << r_c_i << "\n";

        // Track global minimum
        if (r_c_i > 0 && r_c_i < global_min_r_c) {
            global_min_r_c = r_c_i;
            global_bottleneck_node = i;
        }
    }

    // Output global critical analysis correctly
    out_file << "\nGlobal Critical Analysis:\n";
    out_file << "Layer_Count," << LAYER_COUNT << "\n";
    out_file << "p," << p << "\n";
    out_file << "Global_Critical_Rate," << global_min_r_c << "\n";
    out_file << "Bottleneck_Node_ID," << global_bottleneck_node << "\n";
    out_file << "Bottleneck_Layer," << global_bottleneck_node / N << "\n";
    out_file << "Bottleneck_Internal_Betweenness," << bI_ABC[global_bottleneck_node] << "\n";
    out_file << "Bottleneck_External_Betweenness," << bE_ABC[global_bottleneck_node] << "\n";

    out_file.close();

    cout << "Saved individual node analysis to: " << filename << endl;   
}

/**
 * @brief Calculates and outputs the theoretical critical rates across a spectrum of p values.
 * @param filename Name of the output file.
 */
void calculate_theoretical_critical_rates(const string& filename) {
    ofstream out_file(filename);
    if (!out_file) {
        cerr << "Failed to open file: " << filename << endl;
        return;
    }

    out_file << "p,Critical_Rate_Rc,Bottleneck_Node_ID,Bottleneck_Layer,Effective_Betweenness\n";

    for (int p_index = 0; p_index <= 20; p_index++) {
        double current_p = p_index * 0.05;

        double min_r_c = numeric_limits<double>::max();
        int bottleneck_node = -1;
        int bottleneck_layer = -1;
        double max_effective_betweenness = 0;

        for (int i = 0; i < TOTAL_NODES; i++) {
            double effective_betweenness = ((1 - current_p) / (N - 1)) * bI_ABC[i] +
                (current_p / (N * (LAYER_COUNT - 1))) * bE_ABC[i];

            if (effective_betweenness > 0) {
                double r_c_i = C_ave / effective_betweenness;

                if (r_c_i < min_r_c) {
                    min_r_c = r_c_i;
                    bottleneck_node = i;
                    bottleneck_layer = i / N;
                    max_effective_betweenness = effective_betweenness;
                }
            }
        }

        out_file << current_p << "," << min_r_c << ","
            << bottleneck_node << "," << bottleneck_layer << ","
            << max_effective_betweenness << "\n";

        //cout << "p=" << current_p << ": r_c=" << min_r_c
           // << " (Bottleneck node: " << bottleneck_node << ", Layer: " << bottleneck_layer << ")\n";
    }

    out_file.close();
    cout << "Saved overall system performance to: " << filename << endl;
}