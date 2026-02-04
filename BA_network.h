// BA网络生成函数
void initialize_ba_network(int** network, int initial_nodes, int m, int D[], int& degree_sum, unsigned int seed) {
    degree_sum = 0;
    std::mt19937 local_rng(seed);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // 初始化网络全为0
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            network[i][j] = 0;
        }
        D[i] = 0;
    }

    // 步骤1: 创建初始完全连接的小网络
    for (int i = 0; i < initial_nodes; i++) {
        for (int j = i + 1; j < initial_nodes; j++) {
            network[i][j] = network[j][i] = 1;
            D[i]++;
            D[j]++;
            degree_sum += 2;
        }
    }

    // 步骤2: 逐步添加新节点
    for (int new_node = initial_nodes; new_node < N; new_node++) {
        // 计算当前所有节点的度分布概率
        vector<double> probabilities(new_node);
        double total_degree = 0.0;

        for (int i = 0; i < new_node; i++) {
            total_degree += D[i];
        }

        // 计算累积概率
        double cumulative = 0.0;
        for (int i = 0; i < new_node; i++) {
            cumulative += D[i] / total_degree;
            probabilities[i] = cumulative;
        }

        // 新节点连接m个已有节点
        int connections_made = 0;
        while (connections_made < m && connections_made < new_node) {
            double rand_val = dis(local_rng);

            // 根据概率选择连接节点
            for (int i = 0; i < new_node; i++) {
                if (rand_val <= probabilities[i]) {
                    // 检查是否已经连接
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

// 计算BA网络的目标平均度对应的m值
int calculate_m_for_target_degree(double target_avg_degree, int initial_nodes = 5) {
    // BA网络中，平均度 ≈ 2m
    int m = static_cast<int>(round(target_avg_degree / 2.0));
    return max(1, m);  // 确保m至少为1
}

// 构建邻居列表（与ER网络相同）
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