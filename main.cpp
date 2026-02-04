#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>  
#include <ctime>   
#include <vector>
#include <algorithm>
#include <queue>
#include <stack>
#include <random>
#include <cmath>
#include <numeric>  
#include <string>  
#include <iomanip>
#include <chrono>
#include <limits>
#include <unordered_map>
#include <corecrt_math_defines.h>


std::mt19937 global_rng;

using namespace std;

constexpr int N = 1000;      // 网络规模N
constexpr int LAYER_COUNT = 3;   // 定义网络层数，例如5层
constexpr int TOTAL_NODES = LAYER_COUNT * N;  // 总节点数

#include "BA_network.h"
#include "couple_network.h"
#include "shortest_paths.h"
#include "dynamics.h"
#include "calculate_betweenness.h"
#include "queue_analysis.h"
#include "layer_betweenness_analysis.h"
#include "ER_network.h"
#include "RGG_network.h"


// 动态分配内存

int*** networks = new int** [LAYER_COUNT];  // 多层网络 [layer][i][j]

// 耦合网络相关（维度改为TOTAL_NODES）：
int** ABC = new int* [TOTAL_NODES];         // 耦合邻接矩阵
int* D_ABC = new int[TOTAL_NODES];          // 耦合度数组  
int** B_ABC = new int* [TOTAL_NODES];       // 耦合邻居数组

// 新增各层参数数组：
double* k_values = new double[LAYER_COUNT]; // 各层平均度
int* degree_sums = new int[LAYER_COUNT];    // 各层总度数
double* betweenness_layers = new double[LAYER_COUNT * N]; // 各层介数

double* C_ABC = new double[TOTAL_NODES];

double P=1.0;             //层间耦合概率
const int C_ave = 1;
double r_packet_sta, r_packet_end, r_packet_inter;

int** path_num_ABC = new int* [TOTAL_NODES]; // 存储最短路径数量的数组
int*** path_ABC = new int** [TOTAL_NODES];// 三维数组，存储最短路径
int distance_ABC[TOTAL_NODES][TOTAL_NODES];// 存储节点间最短距离
int** degrees = new int* [LAYER_COUNT];  // 各层的度数组

double* bI_ABC = new double[TOTAL_NODES]();  // 层内介数
double* bE_ABC = new double[TOTAL_NODES]();  // 层间介数

int**** path_single;    
int*** path_num_single; 

double beta; 
int current_timestep = 0;
int circuit;
double* eita;
string configType;


int main() {

    cout << "=== 程序开始 ===" << endl;       
    
    // 用户输入参数
    cout << "=== 请输入系统参数 ===" << endl;

    // 耦合比例
    //cout << "请输入层间耦合比例 P (0.0-1.0): ";
    //cin >> P;

    // 跨层传输概率
    //cout << "请输入跨层传输概率 beta (0.0-1.0): ";
    //cin >> beta;

    //string configType;
    cout << "请输入耦合配置类型 (LDC/MDC/HDC): ";
    cin >> configType;

    string network_type = "RGG";  // 新增的网络类型选项
    bool use_rgg = false;         // 控制是否使用RGG

    // 信息包产生率参数
    /*cout << "请输入信息包起始产生率 r_packet_sta (0.0-0.1): ";
    cin >> r_packet_sta;
    cout << "请输入信息包结束产生率 r_packet_end (>起始率): ";
    cin >> r_packet_end;
    cout << "请输入产生率间隔 r_packet_inter (0.001-0.01): ";
    cin >> r_packet_inter;

    // 输入模拟循环次数
    cout << "请输入模拟循环次数 circuit (例如 100): ";
    cin >> circuit;*/

    cout << "请选择网络类型 (1:ER, 2:BA, 3:RGG): ";
    int network_choice;
    cin >> network_choice;

    switch (network_choice) {
    case 1:
        network_type = "ER";
        break;
    case 2:
        network_type = "BA";
        break;
    case 3:
        network_type = "RGG";
        use_rgg = true;
        break;
    default:
        cout << "无效选择，默认使用ER网络" << endl;
        network_type = "ER";
    }

    // 初始化多层网络参数
    k_values[0] = 10.0;  // 层0平均度
    k_values[1] = 20.0;  // 层1平均度  
    k_values[2] = 30.0;  // 层2平均度
    //k_values[3] = 40.0;  // 层3平均度
    //k_values[4] = 50.0;  // 层4平均度
    //k_values[5] = 60.0;  // 层5平均度
    //k_values[6] = 70.0;  // 层6平均度

    // 初始化多层网络数组
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        networks[layer] = new int* [N];
        for (int i = 0; i < N; i++) {
            networks[layer][i] = new int[N];
            for (int j = 0; j < N; j++) {
                networks[layer][i][j] = 0;
            }
        }
    }

    // 初始化耦合网络
    for (int i = 0; i < TOTAL_NODES; i++) {
        ABC[i] = new int[TOTAL_NODES];
        for (int j = 0; j < TOTAL_NODES; j++) {
            ABC[i][j] = 0;
        }
        D_ABC[i] = 0;
        B_ABC[i] = nullptr;
    }

    // 使用当前时间作为全局随机数生成器的种子
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    global_rng.seed(seed);

    if (use_rgg) {
        cout << "\n=== 构建多层随机几何图(RGG)网络 ===" << endl;

        // 初始化degrees数组
        for (int layer = 0; layer < LAYER_COUNT; layer++) {
            degrees[layer] = new int[N];
            for (int i = 0; i < N; i++) {
                degrees[layer][i] = 0;
            }
        }

        // 初始化 degree_sums 数组
        //degree_sums = new int[LAYER_COUNT];
        for (int layer = 0; layer < LAYER_COUNT; layer++) {
            degree_sums[layer] = 0;
        }

        // 存储节点位置（用于分析和可视化）
        vector<vector<pair<double, double>>> rgg_positions;

        // 生成多层RGG - 这里的k_values已在主函数设置
        generate_multi_layer_RGG(networks, degrees, degree_sums,
            k_values, LAYER_COUNT, 42, &rgg_positions);       

    }
    else if (network_type == "BA") {
        cout << "\n=== 构建多层BA网络 ===" << endl;
        cout << "构建多层BA网络..." << endl;

     // 初始化degrees数组
     for (int layer = 0; layer < LAYER_COUNT; layer++) {
         degrees[layer] = new int[N];
         for (int i = 0; i < N; i++) {
             degrees[layer][i] = 0;
         }
     }

     // 构建各层BA网络
     for (int layer = 0; layer < LAYER_COUNT; layer++) {
         // 计算目标平均度对应的m参数
         int m = calculate_m_for_target_degree(k_values[layer]);
         cout << "层 " << layer << ": 目标平均度=" << k_values[layer] << ", m=" << m << endl;

         // 构建BA网络
         initialize_ba_network(networks[layer], 5, m, degrees[layer], degree_sums[layer], 42 + layer);

         // 验证实际平均度
         double actual_avg_degree = static_cast<double>(degree_sums[layer]) / N;
         cout << "层 " << layer << " 实际平均度: " << actual_avg_degree << endl;
     }

    }
    else { // ER
        cout << "\n=== 构建多层ER网络 ===" << endl;
        cout << "构建多层ER网络..." << endl;
        // 初始化degrees数组
        for (int layer = 0; layer < LAYER_COUNT; layer++) {
            degrees[layer] = new int[N];
            for (int i = 0; i < N; i++) {
                degrees[layer][i] = 0;
            }
        }

        for (int layer = 0; layer < LAYER_COUNT; layer++) {
            initialize_er_network(networks[layer], k_values[layer], degrees[layer], degree_sums[layer], 42 + layer);
        }

        cout << "多层网络构建完成，共" << LAYER_COUNT << "层" << endl;
    }
  
    // 构建各层的邻居列表用于最短路径计算
    cout << "计算单层网络最短路径..." << endl;
    initialize_single_path_storage();

    int*** B_layers = new int** [LAYER_COUNT];
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        B_layers[layer] = new int* [N];

        for (int i = 0; i < N; i++) {
            // 重新计算度（确保准确）
            degrees[layer][i] = 0;
            for (int j = 0; j < N; j++) {
                if (networks[layer][i][j]) {
                    degrees[layer][i]++;
                }
            }

            // 分配邻居数组
            B_layers[layer][i] = new int[degrees[layer][i]];
            int count = 0;
            for (int j = 0; j < N; j++) {
                if (networks[layer][i][j]) {
                    B_layers[layer][i][count++] = j;
                }
            }
        }

        calculate_single_layer_paths(layer, networks[layer], degrees[layer], B_layers[layer]);
    }
    cout << "单层网络最短路径计算完成" << endl;


    for (int i = 0; i < TOTAL_NODES; i++) {
        C_ABC[i] = C_ave; // 所有节点统一处理能力
    }

    // 在耦合网络构建之前
    cout << "=== 耦合配置调试输出 ===" << endl;
    debug_coupling_config(k_values, LAYER_COUNT, "LDC");
    debug_coupling_config(k_values, LAYER_COUNT, "MDC");
    debug_coupling_config(k_values, LAYER_COUNT, "HDC");

    // 根据平均度确定耦合配置
    std::vector<int> layerOrder = determineCouplingOrder(k_values, LAYER_COUNT, configType);

    // 创建层索引和平均度的配对用于调试输出
    std::vector<std::pair<double, int>> layerPairs;
    for (int i = 0; i < LAYER_COUNT; i++) {
        layerPairs.push_back({ k_values[i], i });
    }
    std::sort(layerPairs.begin(), layerPairs.end());

    // 打印配置信息
    printCouplingConfig(layerOrder, layerPairs, configType);
   
    // 构建耦合网络
    cout << "耦合多层网络..." << endl;
    generateMultiLayerCoupledNetwork(networks, degrees, layerOrder, P, LAYER_COUNT);
    cout << "网络耦合完成，总节点数：" << TOTAL_NODES << endl;   

    // 在生成耦合网络后调用
    cout << "计算耦合网络最短路径..." << endl;
    multiplex_network_shortest_path_ABC();
    cout << "最短路径计算完成" << endl;

    // 计算介数
    cout << "计算介数..." << endl;
    cal_betweenness_enhanced();
    output_betweenness_results("betweenness_results_" + to_string(LAYER_COUNT) + "layers.csv");

    //analyze_beta_variation_max_node(0.0, 1.0, 0.05);

    // 理论临界值分析
    cout << "进行理论临界值分析..." << endl;
    double theoretical_r_c = calculate_critical_rate();
    cout << "理论临界信息包产生率: " << theoretical_r_c << endl;    

    // 输出临界值分析
    output_critical_analysis("critical_analysis_" + to_string(LAYER_COUNT) + "layers.csv");
     
    // 计算不同beta值的临界值
    calculate_theoretical_critical_rates("theoretical_critical_rates_" + to_string(LAYER_COUNT) + "layers.csv");
    

    eita = new double[circuit + 1];
    
    // 运行模拟
    /*cout << "开始信息包传输模拟..." << endl;
    multi_dyn_spr();
    cout << "模拟完成" << endl;*/


    // 在程序结束前释放多层网络内存
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        for (int i = 0; i < N; i++) {
            delete[] networks[layer][i];
        }
        delete[] networks[layer];
    }
    delete[] networks;

    // 释放degrees数组
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        delete[] degrees[layer];
    }
    delete[] degrees;

    // 释放B_layers数组
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        for (int i = 0; i < N; i++) {
            delete[] B_layers[layer][i];
        }
        delete[] B_layers[layer];
    }
    delete[] B_layers;

    // 释放耦合网络内存
    for (int i = 0; i < TOTAL_NODES; i++) {
        delete[] ABC[i];
        if (B_ABC[i] != nullptr) {
            delete[] B_ABC[i];
        }
    }
    delete[] ABC;
    delete[] B_ABC;
    delete[] D_ABC;

    // 释放其他数组
    delete[] k_values;
    delete[] degree_sums;
    delete[] betweenness_layers;
    delete[] C_ABC;
    delete[] eita;
    // 释放介数数组内存
    delete[] bI_ABC;
    delete[] bE_ABC;
   
    free_ABC_path_storage();
    free_single_path_storage();
    
    // 输出简单的成功信息
    cout << "All results have been exported to files." << endl;
    cout << "\n=== 程序正常结束 ===" << endl;

    return 0;
}

