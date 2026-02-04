extern double* bI_ABC;  // 层内介数
extern double* bE_ABC;  // 层间介数


// 初始化介数数组
void initialize_betweenness_arrays() {
    bI_ABC = new double[TOTAL_NODES]();  // 修改为 TOTAL_NODES
    bE_ABC = new double[TOTAL_NODES]();
}

// 递归计算层内介数
void digui_internal(double pro, int next, int end_node, int layer_index) {
    int w;
    int local_next = next % N;
    int local_end = end_node % N;

    // 检查路径是否存在
    if (path_num_single[layer_index][local_next][local_end] <= 0) {
        return;
    }

    pro = pro / double(path_num_single[layer_index][local_next][local_end]);

    for (w = 0; w < path_num_single[layer_index][local_next][local_end]; w++) {
        int temp_next_hop = path_single[layer_index][local_next][local_end][w];
        int global_node = layer_index * N + temp_next_hop;

        bI_ABC[global_node] += pro;

        if (temp_next_hop != local_end) {
            digui_internal(pro, global_node, end_node, layer_index);
        }
    }
}

// 递归计算层间介数
void digui_external(double pro, int next, int end_node) {
    int w;

    // 检查路径是否存在
    if (path_num_ABC[next][end_node] <= 0) {
        return;
    }

    pro = pro / double(path_num_ABC[next][end_node]);

    for (w = 0; w < path_num_ABC[next][end_node]; w++) {
        int temp_next_hop = path_ABC[next][end_node][w];

        bE_ABC[temp_next_hop] += pro;

        if (temp_next_hop != end_node) {
            digui_external(pro, temp_next_hop, end_node);
        }
    }
}

// 计算层内介数
void calculate_internal_betweenness(int layer_index) {
    int offset = layer_index * N;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) continue;

            int global_i = offset + i;
            int global_j = offset + j;

            if (path_num_single[layer_index][i][j] > 0) {
                double pro = 1.0;
                digui_internal(pro, global_i, global_j, layer_index);
            }
        }
    }
}

// 计算层间介数
void calculate_external_betweenness() {
    // 遍历所有不同层的节点对
    for (int s = 0; s < TOTAL_NODES; s++) {  // 修改为 TOTAL_NODES
        int s_layer = s / N;

        for (int t = 0; t < TOTAL_NODES; t++) {  // 修改为 TOTAL_NODES
            int t_layer = t / N;

            // 只考虑不同层的节点对
            if (s_layer == t_layer) continue;

            if (path_num_ABC[s][t] > 0) {
                double pro = 1.0;
                digui_external(pro, s, t);
            }
        }
    }
}

// 计算增强介数
void cal_betweenness_enhanced() {
    initialize_betweenness_arrays();

    // 计算各层的层内介数
    for (int layer = 0; layer < LAYER_COUNT; layer++) {  // 修改为动态层数
        calculate_internal_betweenness(layer);
    }

    // 计算层间介数
    calculate_external_betweenness();

    cout << "介数计算完成！" << endl;
}

// 输出介数结果
void output_betweenness_results(const string& filename) {
    ofstream out_file(filename);
    if (!out_file) {
        cerr << "无法打开文件: " << filename << endl;
        return;
    }

    out_file << "Node_ID,Layer,Internal_Betweenness,External_Betweenness\n";

    for (int i = 0; i < TOTAL_NODES; i++) {  // 修改为 TOTAL_NODES
        int layer = i / N;
        int node_in_layer = i % N;

        out_file << i << "," << layer << ","
            << bI_ABC[i] << "," << bE_ABC[i] << "\n";
    }

    out_file.close();
    cout << "介数计算结果已保存至: " << filename << endl;
}
