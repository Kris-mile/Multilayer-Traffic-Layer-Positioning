extern double* bI_ABC;      // 层内介数
extern double* bE_ABC;      // 层间介数

// 层最大节点有效介数统计
struct MaxNodeBetweenness {
    int layer_id;
    double beta;           // 当前β值

    // 最大节点的介数值
    double max_bI;         // 最大层内介数（原始值）
    double max_bE;         // 最大层间介数（原始值）
    int max_bI_node_id;    // 最大层内介数节点ID
    int max_bE_node_id;    // 最大层间介数节点ID

    // 最大节点的有效介数
    double effective_bI_max;   // 最大节点的有效层内介数 = (1-β)/(N-1) * max_bI
    double effective_bE_max;   // 最大节点的有效层间介数 = β/[N*(L-1)] * max_bE
    double effective_sum_max;  // 最大节点的总有效介数 = effective_bI_max + effective_bE_max
};

// 计算单个层在特定β值下的最大节点有效介数
MaxNodeBetweenness calculate_max_node_betweenness(int layer_id, double beta_value) {
    MaxNodeBetweenness stats;
    stats.layer_id = layer_id;
    stats.beta = beta_value;

    // 计算有效系数
    double coeff_internal = (1.0 - beta_value) / (N - 1);
    double coeff_external = beta_value / (N * (LAYER_COUNT - 1));

    // 初始化统计量
    stats.max_bI = 0.0;
    stats.max_bE = 0.0;
    stats.max_bI_node_id = -1;
    stats.max_bE_node_id = -1;

    // 遍历该层所有节点，找到最大介数节点
    int start_node = layer_id * N;

    for (int i = 0; i < N; i++) {
        int global_node_id = start_node + i;
        double bI = bI_ABC[global_node_id];
        double bE = bE_ABC[global_node_id];

        // 更新最大层内介数
        if (bI > stats.max_bI) {
            stats.max_bI = bI;
            stats.max_bI_node_id = global_node_id;
        }

        // 更新最大层间介数
        if (bE > stats.max_bE) {
            stats.max_bE = bE;
            stats.max_bE_node_id = global_node_id;
        }
    }

    // 计算最大节点的有效介数
    stats.effective_bI_max = coeff_internal * stats.max_bI;
    stats.effective_bE_max = coeff_external * stats.max_bE;
    stats.effective_sum_max = stats.effective_bI_max + stats.effective_bE_max;

    return stats;
}

// 分析β值变化对所有层最大节点有效介数的影响
void analyze_beta_variation_max_node(double beta_start, double beta_end, double beta_step) {

    cout << "\n=== 分析最大节点有效介数随β值的变化 ===" << endl;
    cout << "β范围: " << beta_start << " → " << beta_end
        << ", 步长: " << beta_step << endl;
    cout << "公式: 最大节点有效层内介数 = (1-β)/(N-1) * max_bI" << endl;
    cout << "公式: 最大节点有效层间介数 = β/[N*(L-1)] * max_bE" << endl;
    cout << "参数: N = " << N << ", L = " << LAYER_COUNT << endl;

    // 为每层创建输出文件
    vector<ofstream> layer_files(LAYER_COUNT);
    vector<string> layer_filenames;

    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        string filename = "betweenness_beta_" + to_string(LAYER_COUNT) +
            "layers_layer" + to_string(layer) + ".csv";
        layer_filenames.push_back(filename);
        layer_files[layer].open(filename);

        if (!layer_files[layer]) {
            cerr << "无法打开文件: " << filename << endl;
            return;
        }

        // 写入表头（最大节点的有效介数）
        layer_files[layer] << "Beta,Max_bI,Max_bE,Max_bI_NodeID,Max_bE_NodeID,"
            << "Effective_bI_max,Effective_bE_max,Effective_Sum_max\n";
    }

    // 创建汇总比较文件
    string summary_filename = "betweenness_beta_" + to_string(LAYER_COUNT) + "layers_summary.csv";
    ofstream summary_file(summary_filename);
    if (!summary_file) {
        cerr << "无法打开汇总文件: " << summary_filename << endl;
        return;
    }

    // 汇总文件表头（每层三列：有效层内、有效层间、总和）
    summary_file << "Beta";
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        summary_file << ",L" << layer << "_Eff_bI_max"
            << ",L" << layer << "_Eff_bE_max"
            << ",L" << layer << "_Eff_Sum_max";
    }
    summary_file << "\n";

    // 遍历β值
    int step_count = 0;
    int total_steps = static_cast<int>((beta_end - beta_start) / beta_step) + 1;

    for (double beta_val = beta_start; beta_val <= beta_end + 1e-9; beta_val += beta_step) {
        step_count++;

        cout << "\n分析 β = " << fixed << setprecision(3) << beta_val
            << " (" << step_count << "/" << total_steps << ")" << endl;

        // 存储当前β值下所有层的统计
        vector<MaxNodeBetweenness> all_layer_stats;

        // 写入汇总数据
        summary_file << fixed << setprecision(3) << beta_val;

        // 计算每层统计并写入文件
        for (int layer = 0; layer < LAYER_COUNT; layer++) {
            auto stats = calculate_max_node_betweenness(layer, beta_val);
            all_layer_stats.push_back(stats);

            // 写入层文件（详细数据）
            layer_files[layer] << fixed << setprecision(6)
                << beta_val << ","
                << stats.max_bI << ","
                << stats.max_bE << ","
                << stats.max_bI_node_id << ","
                << stats.max_bE_node_id << ","
                << stats.effective_bI_max << ","
                << stats.effective_bE_max << ","
                << stats.effective_sum_max << "\n";

            // 写入汇总文件（简化数据：仅有效介数）
            summary_file << fixed << setprecision(6)
                << "," << stats.effective_bI_max
                << "," << stats.effective_bE_max
                << "," << stats.effective_sum_max;

            // 显示关键信息
            cout << "  层" << layer << ": "
                << "最大bI=" << setprecision(4) << stats.max_bI
                << "(节点" << stats.max_bI_node_id << "), "
                << "最大bE=" << stats.max_bE
                << "(节点" << stats.max_bE_node_id << "), "
                << "有效bI_max=" << stats.effective_bI_max
                << ", 有效bE_max=" << stats.effective_bE_max
                << ", 合计=" << stats.effective_sum_max << endl;
        }
        summary_file << "\n";
    }

    // 关闭所有文件
    for (auto& file : layer_files) {
        file.close();
    }
    summary_file.close();

    cout << "\n=== β值变化分析完成 ===" << endl;
    cout << "生成的4个文件:" << endl;
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        cout << "  层" << layer << "数据: " << layer_filenames[layer] << endl;
    }
    cout << "汇总比较文件: " << summary_filename << endl;

    // 显示文件列信息
    cout << "\n文件列结构说明:" << endl;
    cout << "单层文件 (8列): Beta, Max_bI, Max_bE, Max_bI_NodeID, Max_bE_NodeID, "
        << "Effective_bI_max, Effective_bE_max, Effective_Sum_max" << endl;
    cout << "汇总文件 (" << (1 + 3 * LAYER_COUNT) << "列): Beta, ";
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        cout << "L" << layer << "_Eff_bI_max, "
            << "L" << layer << "_Eff_bE_max, "
            << "L" << layer << "_Eff_Sum_max";
        if (layer < LAYER_COUNT - 1) cout << ", ";
    }
    cout << endl;
}

// 快速获取特定β值下所有层的最大节点有效介数
vector<MaxNodeBetweenness> get_all_layers_max_node_stats(double beta_value) {
    vector<MaxNodeBetweenness> all_stats;

    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        all_stats.push_back(calculate_max_node_betweenness(layer, beta_value));
    }

    return all_stats;
}

