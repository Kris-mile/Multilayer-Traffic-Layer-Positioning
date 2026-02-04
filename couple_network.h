extern int** ABC;
extern int* D_ABC;
extern int** B_ABC;
extern string configType;

// 函数声明
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

// 对节点按度数排序
std::vector<int> sortNodesByDegree(int* D, int N) {
    std::vector<int> indices(N);
    for (int i = 0; i < N; i++) {
        indices[i] = i;
    }
    std::sort(indices.begin(), indices.end(), [&](int a, int b) { return D[a] > D[b]; });
    return indices;
}

// 实现coupleMultipleNetworks函数
void coupleMultipleNetworks(int** coupledMatrix, int* coupledDegree,
    int*** networks, int** degrees,
    const std::vector<int>& layerOrder, double P,
    int layerCount, int totalNodes) {

    if (layerOrder.size() < 2) return;

    int total_inter = static_cast<int>(P * N);

    cout << "链式耦合配置: ";
    for (int layer : layerOrder) {
        cout << layer << " ";
    }
    cout << endl;

    // 链式耦合：只连接排序中相邻的层
    for (int i = 0; i < layerOrder.size() - 1; i++) {
        int layerA = layerOrder[i];
        int layerB = layerOrder[i + 1];

        cout << "连接层" << layerA << " <-> 层" << layerB << endl;

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

    cout << "链式耦合完成" << endl;
}

// 耦合配置函数实现
std::vector<int> determineCouplingOrder(double* k_values, int layerCount, const std::string& configType) {
    std::vector<int> order;

    // 创建层索引和平均度的配对
    std::vector<std::pair<double, int>> layerPairs;
    for (int i = 0; i < layerCount; i++) {
        layerPairs.push_back({ k_values[i], i });
    }

    // 按平均度排序（从小到大）
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
        std::cout << "未知配置类型，使用默认LDC配置" << std::endl;
        order = getLDCConfig(layerPairs, layerCount);
    }

    return order;
}

// LDC配置：最小度在中间，剩下的按度从小到大，从中间向两边交替放置
std::vector<int> getLDCConfig(const std::vector<std::pair<double, int>>& layerPairs, int layerCount) {
    std::vector<int> order(layerCount, -1);

    int middleIndex = layerCount / 2;
    order[middleIndex] = layerPairs[0].second;  // 最小度在中间

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

// MDC配置：中间度在中间，按平均度顺序排列
std::vector<int> getMDCConfig(const std::vector<std::pair<double, int>>& layerPairs, int layerCount) {
    std::vector<int> order;

    // 直接按平均度从小到大顺序排列
    for (int i = 0; i < layerCount; i++) {
        order.push_back(layerPairs[i].second);
    }

    return order;
}

// HDC配置：最大度在中间，剩下的按度从大到小，从中间向两边交替放置
std::vector<int> getHDCConfig(const std::vector<std::pair<double, int>>& layerPairs, int layerCount) {
    std::vector<int> order(layerCount, -1);

    int middleIndex = layerCount / 2;
    order[middleIndex] = layerPairs[layerCount - 1].second;  // 最大度在中间

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


void printCouplingConfig(const std::vector<int>& order, const std::vector<std::pair<double, int>>& layerPairs, const std::string& configType) {
    std::cout << configType << "配置 - 耦合顺序: ";
    for (int layer : order) {
        std::cout << layer << " ";
    }
    std::cout << std::endl;

    std::cout << "各层平均度: ";
    for (const auto& pair : layerPairs) {
        std::cout << "层" << pair.second << "(" << pair.first << ") ";
    }
    std::cout << std::endl;
}

void generateMultiLayerCoupledNetwork(int*** networks, int** degrees,
    const std::vector<int>& layerOrder,
    double P, int layerCount) {

    int totalNodes = layerCount * N;

    // 1. 复制各层内部连接
    for (int layer = 0; layer < layerCount; layer++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                int global_i = layer * N + i;
                int global_j = layer * N + j;
                ABC[global_i][global_j] = networks[layer][i][j];
            }
        }
    }

    // 2. 添加层间耦合
    coupleMultipleNetworks(ABC, D_ABC, networks, degrees, layerOrder, P, layerCount, totalNodes);

    // 3. 构建邻居列表
    for (int i = 0; i < totalNodes; i++) {
        if (B_ABC[i] != nullptr) {
            delete[] B_ABC[i];
            B_ABC[i] = nullptr;
        }

        // 计算实际度数
        D_ABC[i] = 0;
        for (int j = 0; j < totalNodes; j++) {
            if (ABC[i][j]) D_ABC[i]++;
        }

        // 分配邻居数组
        B_ABC[i] = new int[D_ABC[i]];
        int count = 0;
        for (int j = 0; j < totalNodes; j++) {
            if (ABC[i][j]) {
                B_ABC[i][count++] = j;
            }
        }
    }
}

// 耦合配置调试函数
void debug_coupling_config(double* k_values, int layerCount, const string& configType) {
    vector<pair<double, int>> layerPairs;
    for (int i = 0; i < layerCount; i++) {
        layerPairs.push_back({ k_values[i], i });
    }
    sort(layerPairs.begin(), layerPairs.end());

    vector<int> order = determineCouplingOrder(k_values, layerCount, configType);

    cout << "=== " << configType << "配置调试信息 ===" << endl;
    cout << "总层数: " << layerCount << endl;
    cout << "各层原始平均度: ";
    for (int i = 0; i < layerCount; i++) {
        cout << "层" << i << "(k=" << k_values[i] << ") ";
    }
    cout << endl;

    cout << "按平均度排序后: ";
    for (const auto& p : layerPairs) {
        cout << "层" << p.second << "(k=" << p.first << ") ";
    }
    cout << endl;

    cout << "中间位置索引: " << (layerCount / 2) << endl;

    cout << "耦合顺序(层编号): ";
    for (int layer : order) {
        cout << layer << " ";
    }
    cout << endl;

    cout << "耦合顺序(平均度): ";
    for (int layer : order) {
        for (int i = 0; i < layerCount; i++) {
            if (i == layer) {
                cout << k_values[i] << " ";
                break;
            }
        }
    }
    cout << endl;

    cout << "位置分布: [";
    for (int i = 0; i < layerCount; i++) {
        cout << order[i];
        if (i < layerCount - 1) cout << ",";
    }
    cout << "]" << endl;

    // 验证中间层是否正确
    int middleIndex = layerCount / 2;
    if (configType == "LDC") {
        cout << "验证: 中间位置" << middleIndex << "应该是层" << layerPairs[0].second
            << "(最小度k=" << layerPairs[0].first << ")" << endl;
    }
    else if (configType == "HDC") {
        cout << "验证: 中间位置" << middleIndex << "应该是层" << layerPairs[layerCount - 1].second
            << "(最大度k=" << layerPairs[layerCount - 1].first << ")" << endl;
    }
    cout << endl;
}
