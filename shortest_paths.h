extern int** ABC;
extern int* D_ABC;
extern int** B_ABC;

extern int*** path_ABC;
extern int** path_num_ABC;
extern int distance_ABC[TOTAL_NODES][TOTAL_NODES];

extern int**** path_single;    // 四维数组: [layer][i][j][k]
extern int*** path_num_single; // 单层网络的最短路径数量 [layer][i][j]
extern double P;

// 初始化路径存储结构
void initialize_ABC_path_storage() {
    for (int i = 0; i < TOTAL_NODES; i++) {
        path_num_ABC[i] = new int[TOTAL_NODES];
        path_ABC[i] = new int* [TOTAL_NODES];

        for (int j = 0; j < TOTAL_NODES; j++) {
            path_num_ABC[i][j] = 0;
            path_ABC[i][j] = nullptr;
            distance_ABC[i][j] = -1;
        }
    }
}

void free_ABC_path_storage() {
    for (int i = 0; i < TOTAL_NODES; i++) {
        for (int j = 0; j < TOTAL_NODES; j++) {
            if (path_ABC[i][j] != nullptr) {
                delete[] path_ABC[i][j];
            }
        }
        delete[] path_num_ABC[i];
        delete[] path_ABC[i];
    }
    delete[] path_num_ABC;
    delete[] path_ABC;
}

void multiplex_network_shortest_path_ABC() {
    initialize_ABC_path_storage();

    // 初始化距离矩阵
    for (int i = 0; i < TOTAL_NODES; i++) {
        for (int j = 0; j < TOTAL_NODES; j++) {
            if (i == j) distance_ABC[i][j] = 0;
            else distance_ABC[i][j] = (ABC[i][j] ? 1 : -1);
        }
    }

    int allow_vertices[TOTAL_NODES];//标记节点是否已处理

    // 外层循环：遍历每个源节点
    for (int node = 0; node < TOTAL_NODES; node++) {
        int next = node;

        // 节点到自身的距离为0
        distance_ABC[node][node] = 0;

        // 初始化标记数组
        for (int i = 0; i < TOTAL_NODES; i++) {
            allow_vertices[i] = 0;
        }

        // 内层循环：逐步扩展最短路径集合
        for (int allow_size = 0; allow_size < TOTAL_NODES; allow_size++) {
            int maximum = INT_MAX;
            int next = -1;

            // 在未处理节点中找出当前距离最小的那个，并更新它的邻居
            for (int i = 0; i < TOTAL_NODES; i++) {

                // 找到更小距离且未加入的节点
                if (distance_ABC[node][i] != -1 && allow_vertices[i] == 0) {
                    if (distance_ABC[node][i] < maximum) {
                        maximum = distance_ABC[node][i];
                        next = i;// 更新为当前最近节点
                    }
                }
            }

            if (next == -1) break;
            allow_vertices[next] = 1;

            // 更新邻居节点的距离
            for (int i = 0; i < D_ABC[next]; i++) {
                int neighbor = B_ABC[next][i];

                // 跳过已加入集合的节点
                if (allow_vertices[neighbor] == 0) {

                    // 计算新距离 = 源点到当前节点距离 + 1
                    int new_distance = distance_ABC[node][next] + 1;

                    if (distance_ABC[node][neighbor] == -1 || new_distance < distance_ABC[node][neighbor]) {
                        distance_ABC[node][neighbor] = new_distance;
                    }
                }
            }
        }
    }

    // 路径记录部分
    for (int i = 0; i < TOTAL_NODES; i++) {
        for (int j = 0; j < TOTAL_NODES; j++) {

            if (i == j || distance_ABC[i][j] == -1) {
                path_num_ABC[i][j] = 0;
                path_ABC[i][j] = nullptr;
                continue;
            }

            // 计算路径数量
            path_num_ABC[i][j] = 0;
            for (int k = 0; k < D_ABC[i]; k++) {
                if (distance_ABC[B_ABC[i][k]][j] == distance_ABC[i][j] - 1) {
                    path_num_ABC[i][j]++;
                }
            }

            // 分配内存并存储路径
            if (path_num_ABC[i][j] > 0) {
                // 每条路径只存储下一跳节点
                path_ABC[i][j] = new int[path_num_ABC[i][j]];

                // 初始化路径数组
                for (int k = 0; k < path_num_ABC[i][j]; k++)
                    path_ABC[i][j][k] = -1; // 初始化为无效值

                int count = 0;
                for (int k = 0; k < D_ABC[i]; k++) {
                    int neighbor = B_ABC[i][k];
                    if (distance_ABC[neighbor][j] == distance_ABC[i][j] - 1) {
                        if (count < path_num_ABC[i][j]) {
                            path_ABC[i][j][count] = neighbor;
                            count++;
                        }
                    }
                }

                // 验证实际填充的数量
                if (count != path_num_ABC[i][j]) {
                    cerr << "ERROR: Expected " << path_num_ABC[i][j] << " paths, but filled " << count << endl;
                }
            }
            else {
                path_ABC[i][j] = nullptr;
            }
        }
    }

    // 验证路径
    int invalid_paths = 0;
    for (int i = 0; i < TOTAL_NODES; i++) {
        for (int j = 0; j < TOTAL_NODES; j++) {
            if (path_num_ABC[i][j] > 0) {
                if (path_ABC[i][j] == nullptr) {
                    cerr << "错误：节点 " << i << " -> " << j
                        << " 路径数>0但数组为空！" << endl;
                    invalid_paths++;
                }
                else {
                    for (int k = 0; k < path_num_ABC[i][j]; k++) {
                        int neighbor = path_ABC[i][j][k];
                        if (neighbor < 0 || neighbor >= TOTAL_NODES) {
                            cerr << "错误：无效邻居节点 " << neighbor
                                << " (范围:0-" << (TOTAL_NODES - 1) << ")" << endl;
                            invalid_paths++;
                        }
                    }
                }
            }
        }
    }

    int unreachable_count = 0;
    for (int i = 0; i < TOTAL_NODES; i++) {
        for (int j = 0; j < TOTAL_NODES; j++) {
            if (distance_ABC[i][j] == -1) {
                unreachable_count++;
            }
        }
    }
    cout << "Total unreachable pairs: " << unreachable_count << endl;
    cout << "Found " << invalid_paths << " invalid paths" << endl;
    cout << "三层耦合网络最短路径计算完毕！" << endl;
}

// 初始化单层路径存储结构
void initialize_single_path_storage() {
    path_single = new int*** [LAYER_COUNT];
    path_num_single = new int** [LAYER_COUNT];

    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        path_single[layer] = new int** [N];
        path_num_single[layer] = new int* [N];

        for (int i = 0; i < N; i++) {
            path_single[layer][i] = new int* [N];
            path_num_single[layer][i] = new int[N];

            for (int j = 0; j < N; j++) {
                path_single[layer][i][j] = nullptr; // 初始化为空指针
                path_num_single[layer][i][j] = 0;   // 初始路径数为0
            }
        }
    }
}

void free_single_path_storage() {
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (path_single[layer][i][j] != nullptr) {
                    delete[] path_single[layer][i][j];
                }
            }
            delete[] path_single[layer][i];
            delete[] path_num_single[layer][i];
        }
        delete[] path_single[layer];
        delete[] path_num_single[layer];
    }
    delete[] path_single;
    delete[] path_num_single;
}

// 计算单层网络的最短路径
void calculate_single_layer_paths(int layer_index, int** network, int* D, int** B) {
    // 为当前层创建临时距离矩阵
    int** distance = new int* [N];
    for (int i = 0; i < N; i++) {
        distance[i] = new int[N];
        for (int j = 0; j < N; j++) {
            if (i == j) distance[i][j] = 0;
            else distance[i][j] = (network[i][j] ? 1 : -1);
        }
    }

    int allow_vertices[TOTAL_NODES]; // 标记节点是否已处理

    // 外层循环：遍历每个源节点
    for (int node = 0; node < N; node++) {
        // 节点到自身的距离为0
        distance[node][node] = 0;

        // 初始化标记数组
        for (int i = 0; i < N; i++) {
            allow_vertices[i] = 0;
        }

        // 内层循环：逐步扩展最短路径集合
        for (int allow_size = 0; allow_size < N; allow_size++) {
            int maximum = INT_MAX;
            int next = -1;

            // 在未处理节点中找出当前距离最小的那个，并更新它的邻居
            for (int i = 0; i < N; i++) {
                // 找到更小距离且未加入的节点
                if (distance[node][i] != -1 && allow_vertices[i] == 0) {
                    if (distance[node][i] < maximum) {
                        maximum = distance[node][i];
                        next = i; // 更新为当前最近节点
                    }
                }
            }

            if (next == -1) break;
            allow_vertices[next] = 1;

            // 更新邻居节点的距离
            for (int i = 0; i < D[next]; i++) {
                int neighbor = B[next][i];

                // 跳过已加入集合的节点
                if (allow_vertices[neighbor] == 0) {
                    // 计算新距离 = 源点到当前节点距离 + 1
                    int new_distance = distance[node][next] + 1;

                    if (distance[node][neighbor] == -1 || new_distance < distance[node][neighbor]) {
                        distance[node][neighbor] = new_distance;
                    }
                }
            }
        }
    }

    // 路径记录部分
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j || distance[i][j] == -1) {
                path_num_single[layer_index][i][j] = 0;
                path_single[layer_index][i][j] = nullptr;
                continue;
            }

            // 计算路径数量
            path_num_single[layer_index][i][j] = 0;
            for (int k = 0; k < D[i]; k++) {
                if (distance[B[i][k]][j] == distance[i][j] - 1) {
                    path_num_single[layer_index][i][j]++;
                }
            }

            // 分配内存并存储路径
            if (path_num_single[layer_index][i][j] > 0) {
                // 每条路径只存储下一跳节点
                path_single[layer_index][i][j] = new int[path_num_single[layer_index][i][j]];

                // 初始化路径数组
                for (int k = 0; k < path_num_single[layer_index][i][j]; k++)
                    path_single[layer_index][i][j][k] = -1; // 初始化为无效值

                int count = 0;
                for (int k = 0; k < D[i]; k++) {
                    int neighbor = B[i][k];
                    if (distance[neighbor][j] == distance[i][j] - 1) {
                        if (count < path_num_single[layer_index][i][j]) {
                            path_single[layer_index][i][j][count] = neighbor;
                            count++;
                        }
                    }
                }
            }
            else {
                path_single[layer_index][i][j] = nullptr;
            }
        }
    }    

    // 释放临时距离矩阵
    for (int i = 0; i < N; i++) {
        delete[] distance[i];
    }
    delete[] distance;
}

