extern double* bI_ABC;  // 层内介数
extern double* bE_ABC;  // 层间介数
extern double beta;      // 跨层传输概率

// 计算每个节点的理论队列长度
double* calculate_theoretical_queue_length(double r) {
    double* q_i = new double[TOTAL_NODES]();  // 修改为 TOTAL_NODES

    for (int i = 0; i < TOTAL_NODES; i++) {
        double q_I = (r * (1 - beta) / (N - 1)) * bI_ABC[i];
        double q_E = (r * beta / (N * (LAYER_COUNT - 1))) * bE_ABC[i];  // 使用 LAYER_COUNT
        q_i[i] = q_I + q_E;
    }

    return q_i;
}

// 计算临界信息包产生率
double calculate_critical_rate() {
    double min_r_c = numeric_limits<double>::max();
    int bottleneck_node = -1;

    for (int i = 0; i < TOTAL_NODES; i++) {  // 修改为 TOTAL_NODES
        double denominator = ((1 - beta) / (N - 1)) * bI_ABC[i] +
            (beta / (N * (LAYER_COUNT - 1))) * bE_ABC[i];  // 使用 LAYER_COUNT

        if (denominator > 0) {
            double r_c_i = C_ave / denominator;
            if (r_c_i < min_r_c) {
                min_r_c = r_c_i;
                bottleneck_node = i;
            }
        }
    }

    cout << "理论临界信息包产生率 r_c = " << min_r_c << endl;
    if (bottleneck_node != -1) {
        cout << "瓶颈节点: " << bottleneck_node
            << " (层 " << bottleneck_node / N << ")" << endl;
    }

    return min_r_c;
}

// 输出队列分析结果
void output_queue_analysis(double r, const string& filename) {
    ofstream out_file(filename);
    if (!out_file) {
        cerr << "无法打开文件: " << filename << endl;
        return;
    }

    out_file << "Node_ID,Layer,Internal_Betweenness,External_Betweenness,Queue_Length\n";

    double* q_i = calculate_theoretical_queue_length(r);

    for (int i = 0; i < TOTAL_NODES; i++) {  // 修改为 TOTAL_NODES
        int layer = i / N;
        out_file << i << "," << layer << ","
            << bI_ABC[i] << "," << bE_ABC[i] << ","
            << q_i[i] << "\n";
    }

    delete[] q_i;
    out_file.close();
    cout << "队列分析结果已保存至: " << filename << endl;
}

// 输出临界值分析
void output_critical_analysis(const string& filename) {
    ofstream out_file(filename);
    if (!out_file) {
        cerr << "无法打开文件: " << filename << endl;
        return;
    }

    out_file << "Node_ID,Layer,Internal_Betweenness,External_Betweenness,Critical_Rate\n";

    double global_min_r_c = numeric_limits<double>::max();
    int global_bottleneck_node = -1;

    for (int i = 0; i < TOTAL_NODES; i++) {  // 修改为 TOTAL_NODES
        int layer = i / N;
        double denominator = ((1 - beta) / (N - 1)) * bI_ABC[i] +
            (beta / (N * (LAYER_COUNT - 1))) * bE_ABC[i];  // 使用 LAYER_COUNT

        double r_c_i = (denominator > 0) ? C_ave / denominator : 0;

        out_file << i << "," << layer << ","
            << bI_ABC[i] << "," << bE_ABC[i] << ","
            << r_c_i << "\n";

        // 跟踪全局最小值
        if (r_c_i > 0 && r_c_i < global_min_r_c) {
            global_min_r_c = r_c_i;
            global_bottleneck_node = i;
        }
    }

    // 正确输出全局临界值
    out_file << "\nGlobal Critical Analysis:\n";
    out_file << "Layer_Count," << LAYER_COUNT << "\n";
    out_file << "Beta," << beta << "\n";
    out_file << "Global_Critical_Rate," << global_min_r_c << "\n";
    out_file << "Bottleneck_Node_ID," << global_bottleneck_node << "\n";
    out_file << "Bottleneck_Layer," << global_bottleneck_node / N << "\n";
    out_file << "Bottleneck_Internal_Betweenness," << bI_ABC[global_bottleneck_node] << "\n";
    out_file << "Bottleneck_External_Betweenness," << bE_ABC[global_bottleneck_node] << "\n";

    out_file.close();

    cout << "临界值分析已保存至: " << filename << endl;
    cout << "理论临界信息包产生率 r_c = " << global_min_r_c << endl;
    cout << "瓶颈节点: " << global_bottleneck_node
        << " (层 " << global_bottleneck_node / N << ")" << endl;
}

// 计算不同beta值的理论临界值
void calculate_theoretical_critical_rates(const string& filename) {
    ofstream out_file(filename);
    if (!out_file) {
        cerr << "无法打开文件: " << filename << endl;
        return;
    }

    out_file << "Beta,Critical_Rate_Rc,Bottleneck_Node_ID,Bottleneck_Layer,Effective_Betweenness\n";

    for (int beta_index = 0; beta_index <= 20; beta_index++) {
        double current_beta = beta_index * 0.05;

        double min_r_c = numeric_limits<double>::max();
        int bottleneck_node = -1;
        int bottleneck_layer = -1;
        double max_effective_betweenness = 0;

        for (int i = 0; i < TOTAL_NODES; i++) {  // 修改为 TOTAL_NODES
            double effective_betweenness = ((1 - current_beta) / (N - 1)) * bI_ABC[i] +
                (current_beta / (N * (LAYER_COUNT - 1))) * bE_ABC[i];  // 使用 LAYER_COUNT

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

        out_file << current_beta << "," << min_r_c << ","
            << bottleneck_node << "," << bottleneck_layer << ","
            << max_effective_betweenness << "\n";

        cout << "β=" << current_beta << ": r_c=" << min_r_c
            << " (瓶颈节点: " << bottleneck_node << ", 层: " << bottleneck_layer << ")\n";
    }

    out_file.close();
    cout << "理论临界值分析已保存至: " << filename << endl;
}

