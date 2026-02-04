// 初始化ER网络（使用平均度而不是连接概率）
void initialize_er_network(int** network, double average_degree, int D[], int& degree_sum, unsigned int seed) {
    degree_sum = 0;
    std::mt19937 local_rng(seed);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // 计算对应的连接概率
    double p = average_degree / (N - 1);

    // 初始化所有元素为0
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            network[i][j] = 0;
        }
        D[i] = 0;
    }

    // 随机连接节点
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

// 构建邻居列表
/*void build_B(int** network, int** neighbor) {
    for (int i = 0; i < N; i++) {
        int count = 0;
        for (int j = 0; j < N; j++) {
            if (network[i][j] != 0) {
                neighbor[i][count++] = j;
            }
        }
    }
}*/