// 优化的RGG网络生成函数
void initialize_rgg_network_optimized(int** network, double target_avg_degree,
    int D[], int& degree_sum,
    unsigned int seed,
    vector<pair<double, double>>& positions) {
    degree_sum = 0;
    std::mt19937 local_rng(seed);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // 清空网络
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            network[i][j] = 0;
        }
        D[i] = 0;
    }

    positions.clear();
    positions.reserve(N);

    // 1. 快速生成节点位置
    for (int i = 0; i < N; i++) {
        positions.push_back({ dis(local_rng), dis(local_rng) });
    }

    // 2. 估算半径（避免平方根计算）
    double r = sqrt(target_avg_degree / (M_PI * (N - 1)));
    double r_squared = r * r;  // 使用平方比较，避免sqrt

    // 3. 网格划分优化：空间分割技术
    int grid_size = max(1, static_cast<int>(1.0 / r));
    if (grid_size > 100) grid_size = 100;  // 限制网格数量

    vector<vector<int>> grid(grid_size * grid_size);

    // 将节点分配到网格
    for (int i = 0; i < N; i++) {
        int gx = static_cast<int>(positions[i].first * grid_size);
        int gy = static_cast<int>(positions[i].second * grid_size);
        gx = min(gx, grid_size - 1);
        gy = min(gy, grid_size - 1);
        int cell_id = gy * grid_size + gx;
        grid[cell_id].push_back(i);
    }

    // 4. 快速连接：只检查相邻网格
    int connections = 0;
    for (int cell = 0; cell < grid.size(); cell++) {
        const vector<int>& nodes_in_cell = grid[cell];

        // 检查当前网格和相邻网格
        int gx = cell % grid_size;
        int gy = cell / grid_size;

        // 只检查当前网格和右侧/下方的相邻网格（避免重复）
        vector<int> neighbor_cells;
        neighbor_cells.push_back(cell);  // 当前网格

        // 右侧网格
        if (gx + 1 < grid_size) {
            neighbor_cells.push_back(cell + 1);
            // 右下网格
            if (gy + 1 < grid_size) neighbor_cells.push_back(cell + grid_size + 1);
        }
        // 下方网格
        if (gy + 1 < grid_size) {
            neighbor_cells.push_back(cell + grid_size);
            // 左下网格（如果可能）
            if (gx - 1 >= 0) neighbor_cells.push_back(cell + grid_size - 1);
        }

        // 连接节点
        for (int i_idx = 0; i_idx < nodes_in_cell.size(); i_idx++) {
            int i = nodes_in_cell[i_idx];
            double xi = positions[i].first;
            double yi = positions[i].second;

            // 检查当前网格内的其他节点
            for (int j_idx = i_idx + 1; j_idx < nodes_in_cell.size(); j_idx++) {
                int j = nodes_in_cell[j_idx];
                double dx = xi - positions[j].first;
                double dy = yi - positions[j].second;
                double dist_sq = dx * dx + dy * dy;

                if (dist_sq <= r_squared) {
                    network[i][j] = network[j][i] = 1;
                    D[i]++; D[j]++;
                    degree_sum += 2;
                    connections++;
                }
            }

            // 检查相邻网格中的节点
            for (int neighbor_cell : neighbor_cells) {
                if (neighbor_cell == cell) continue;  // 当前网格已处理

                const vector<int>& neighbor_nodes = grid[neighbor_cell];
                for (int j : neighbor_nodes) {
                    double dx = xi - positions[j].first;
                    double dy = yi - positions[j].second;
                    double dist_sq = dx * dx + dy * dy;

                    if (dist_sq <= r_squared) {
                        network[i][j] = network[j][i] = 1;
                        D[i]++; D[j]++;
                        degree_sum += 2;
                        connections++;
                    }
                }
            }
        }
    }

    // 5. 快速调整半径（如果必要）
    double actual_avg_degree = static_cast<double>(degree_sum) / N;

    if (fabs(actual_avg_degree - target_avg_degree) > target_avg_degree * 0.2) {
        // 使用快速估计算法调整
        double factor = sqrt(target_avg_degree / actual_avg_degree);
        r_squared *= factor * factor;

        // 重新连接（使用新半径）
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) network[i][j] = 0;
            D[i] = 0;
        }
        degree_sum = 0;
        connections = 0;

        // 简化连接：不重新网格化
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double dx = positions[i].first - positions[j].first;
                double dy = positions[i].second - positions[j].second;
                double dist_sq = dx * dx + dy * dy;

                if (dist_sq <= r_squared) {
                    network[i][j] = network[j][i] = 1;
                    D[i]++; D[j]++;
                    degree_sum += 2;
                    connections++;
                }
            }
        }
    }

    actual_avg_degree = static_cast<double>(degree_sum) / N;
    cout << "快速RGG生成: 平均度=" << actual_avg_degree
        << " (目标=" << target_avg_degree << "), 连接数=" << connections << endl;

    // 6. 确保连通性（快速检查）
    vector<bool> visited(N, false);
    queue<int> q;
    q.push(0);
    visited[0] = true;
    int visited_count = 1;

    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v = 0; v < N; v++) {
            if (network[u][v] && !visited[v]) {
                visited[v] = true;
                q.push(v);
                visited_count++;
            }
        }
    }

    if (visited_count < N) {
        cout << "警告: 网络有 " << (N - visited_count) << " 个孤立组件" << endl;
        // 快速连接主要组件（可选）
    }
}

// 快速多层RGG生成
void generate_multi_layer_RGG(int*** networks, int** degrees, int* degree_sums,
    double* k_values, int layer_count,
    unsigned int seed,
    vector<vector<pair<double, double>>>* all_positions = nullptr) {

    cout << "\n=== 快速生成多层RGG网络 ===" << endl;
    auto start_time = chrono::high_resolution_clock::now();

    if (all_positions != nullptr) {
        all_positions->resize(layer_count);
    }

    // 并行生成各层（如果支持OpenMP）
#pragma omp parallel for if(layer_count >= 3)
    for (int layer = 0; layer < layer_count; layer++) {
        unsigned int layer_seed = seed + layer * 1000;
        vector<pair<double, double>> positions;

        initialize_rgg_network_optimized(networks[layer], k_values[layer],
            degrees[layer], degree_sums[layer],
            layer_seed, positions);

        if (all_positions != nullptr) {
            (*all_positions)[layer] = positions;
        }
    }

    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);

    cout << "多层RGG生成完成，耗时: " << duration.count() << " ms" << endl;
    cout << "总节点数: " << layer_count * N << endl;
}

