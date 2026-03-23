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
#include <climits>
#include <omp.h>

std::mt19937 global_rng;

using namespace std;

constexpr int N = 1000;          // Network size N (Number of nodes per layer)
constexpr int LAYER_COUNT = 3;   // Number of layers in the multiplex network
constexpr int TOTAL_NODES = LAYER_COUNT * N;  // Total number of nodes in the system

// Custom Header Inclusions
#include "BA_network.h"
#include "couple_network.h"
#include "shortest_paths.h"
#include "dynamics.h"
#include "calculate_betweenness.h"
#include "queue_analysis.h"
#include "layer_betweenness_analysis.h"
#include "ER_network.h"
#include "RGG_network.h"

// --- Dynamic Memory Allocation ---

int*** networks = new int** [LAYER_COUNT];  // Multilayer network adjacency matrix: [layer][i][j]

// Coupled network variables (Dimensions scaled to TOTAL_NODES):
int** ABC = new int* [TOTAL_NODES];         // Coupled adjacency matrix
int* D_ABC = new int[TOTAL_NODES];          // Coupled degree array
int** B_ABC = new int* [TOTAL_NODES];       // Coupled neighbor array

// Layer-specific parameter arrays:
double* k_values = new double[LAYER_COUNT]; // Average degree of each layer
int* degree_sums = new int[LAYER_COUNT];    // Sum of degrees in each layer
double* betweenness_layers = new double[LAYER_COUNT * N]; // Betweenness centrality for nodes in each layer

double* C_ABC = new double[TOTAL_NODES];    // Node processing capacity

double coupling_prob = 1.0;             // Inter-layer coupling probability
const int C_ave = 1;        // Average node processing capacity
double r_packet_sta, r_packet_end, r_packet_inter; // Packet generation rate parameters

int** path_num_ABC = new int* [TOTAL_NODES]; // Array to store the number of shortest paths
int*** path_ABC = new int** [TOTAL_NODES];   // 3D array to store the actual shortest paths
int distance_ABC[TOTAL_NODES][TOTAL_NODES];  // Matrix to store shortest path distances between nodes
int** degrees = new int* [LAYER_COUNT];      // Degree array for each distinct layer

double* bI_ABC = new double[TOTAL_NODES]();  // Intra-layer betweenness centrality
double* bE_ABC = new double[TOTAL_NODES]();  // Inter-layer betweenness centrality

int**** path_single;
int*** path_num_single;

double p;
int current_timestep = 0;
int circuit;
double* eita;
string configType;
string network_type = "ER";  // Default network type
bool use_rgg = false;         // Flag to control RGG usage

int main() {

    cout << "=== Program Started ===" << endl;

    // User Input Parameters
    cout << "=== Please enter system parameters ===" << endl;    
    
    cout << "Cross-layer prob p (0.0-1.0): ";
    cin >> p;

    cout << "Coupling type (LDC/MDC/HDC): ";
    cin >> configType;

    
    // Packet generation rate inputs  
    cout << "Packet rate start(0.0-0.2): ";
    cin >> r_packet_sta;
    cout << "Packet rate end(> start rate): ";
    cin >> r_packet_end;
    cout << "Packet rate step(0.001-0.01): ";
    cin >> r_packet_inter;
    cout << "Simulations count(e.g., 50-100): ";
    cin >> circuit;
    
    cout << "Select network topology (1: ER, 2: BA, 3: RGG): ";
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
        cout << "Invalid choice. Defaulting to ER network." << endl;
        network_type = "ER";
    }

    // Initialize multilayer network parameters
    k_values[0] = 10.0;  // Average degree for Layer 0
    k_values[1] = 20.0;  // Average degree for Layer 1  
    k_values[2] = 30.0;  // Average degree for Layer 2

    // Allocate memory for multilayer network arrays
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        networks[layer] = new int* [N];
        for (int i = 0; i < N; i++) {
            networks[layer][i] = new int[N];
            for (int j = 0; j < N; j++) {
                networks[layer][i][j] = 0;
            }
        }
    }

    // Allocate memory for coupled network arrays
    for (int i = 0; i < TOTAL_NODES; i++) {
        ABC[i] = new int[TOTAL_NODES];
        for (int j = 0; j < TOTAL_NODES; j++) {
            ABC[i][j] = 0;
        }
        D_ABC[i] = 0;
        B_ABC[i] = nullptr;
    }

    // Seed the global random number generator using the current system time
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    global_rng.seed(seed);

    if (use_rgg) {
        cout << "\n=== Constructing Multilayer Random Geometric Graph (RGG) ===" << endl;

        // Initialize degrees array
        for (int layer = 0; layer < LAYER_COUNT; layer++) {
            degrees[layer] = new int[N];
            for (int i = 0; i < N; i++) {
                degrees[layer][i] = 0;
            }
        }

        // Initialize degree_sums array
        for (int layer = 0; layer < LAYER_COUNT; layer++) {
            degree_sums[layer] = 0;
        }

        // Store node positions for potential analysis and visualization
        vector<vector<pair<double, double>>> rgg_positions;

        // Generate Multilayer RGG 
        generate_multi_layer_RGG(networks, degrees, degree_sums,
            k_values, LAYER_COUNT, 42, &rgg_positions);

    }
    else if (network_type == "BA") {
        cout << "\n=== Constructing Multilayer Scale-Free (BA) Network ===" << endl;

        // Initialize degrees array
        for (int layer = 0; layer < LAYER_COUNT; layer++) {
            degrees[layer] = new int[N];
            for (int i = 0; i < N; i++) {
                degrees[layer][i] = 0;
            }
        }

        // Construct BA network for each layer
        for (int layer = 0; layer < LAYER_COUNT; layer++) {
            // Calculate parameter 'm' based on target average degree
            int m = calculate_m_for_target_degree(k_values[layer]);
            cout << "Layer " << layer << ": Target Avg Degree=" << k_values[layer] << ", m=" << m << endl;

            // Generate BA network
            initialize_ba_network(networks[layer], 5, m, degrees[layer], degree_sums[layer], 42 + layer);

            // Verify actual average degree
            double actual_avg_degree = static_cast<double>(degree_sums[layer]) / N;
            cout << "Layer " << layer << " Actual Avg Degree: " << actual_avg_degree << endl;
        }

    }
    else { // ER
        cout << "\n=== Constructing Multilayer ER Network ===" << endl;

        // Initialize degrees array
        for (int layer = 0; layer < LAYER_COUNT; layer++) {
            degrees[layer] = new int[N];
            for (int i = 0; i < N; i++) {
                degrees[layer][i] = 0;
            }
        }

        // Construct ER network for each layer
        for (int layer = 0; layer < LAYER_COUNT; layer++) {
            initialize_er_network(networks[layer], k_values[layer], degrees[layer], degree_sums[layer], 42 + layer);
        }
    }
   
    cout << "Layer Network: " <<  LAYER_COUNT << " Total nodes: " << TOTAL_NODES << endl;
    // Construct neighbor lists for shortest path calculation 
    cout << "Calculating single-layer shortest paths..." << endl;  
    initialize_single_path_storage();

    int*** B_layers = new int** [LAYER_COUNT];
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        B_layers[layer] = new int* [N];

        for (int i = 0; i < N; i++) {
            // Recalculate degrees to ensure accuracy
            degrees[layer][i] = 0;
            for (int j = 0; j < N; j++) {
                if (networks[layer][i][j]) {
                    degrees[layer][i]++;
                }
            }

            // Allocate and populate neighbor array
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

    for (int i = 0; i < TOTAL_NODES; i++) {
        C_ABC[i] = C_ave; // Uniform processing capacity for all nodes
    }   

    // Determine the coupling order based on average degrees and strategy
    std::vector<int> layerOrder = determineCouplingOrder(k_values, LAYER_COUNT, configType);

    // Create pairs of (average degree, layer index) for debugging output
    std::vector<std::pair<double, int>> layerPairs;
    for (int i = 0; i < LAYER_COUNT; i++) {
        layerPairs.push_back({ k_values[i], i });
    }
    std::sort(layerPairs.begin(), layerPairs.end());
    
    printCouplingConfig(layerOrder, layerPairs, configType);
   
    cout << "Coupling the multilayer network..." << endl;
    generateMultiLayerCoupledNetwork(networks, degrees, layerOrder, coupling_prob, LAYER_COUNT);

    cout << "Calculating shortest paths for the coupled network..." << endl;
    multiplex_network_shortest_path_ABC();

    // Calculate betweenness centrality
    cout << "Calculating betweenness centrality...";
    cal_betweenness_enhanced();
    output_betweenness_results("betweenness_results_" + to_string(LAYER_COUNT) + "layers.csv");

    // Theoretical critical value analysis
    cout << "\n===Performing theoretical critical rate analysis...===" << endl;
    double theoretical_r_c = calculate_critical_rate();

    // Output critical analysis results
    output_critical_analysis("criticality_by_node_" + to_string(LAYER_COUNT) + "layers.csv");

    // Calculate critical rates across varying beta values
    calculate_theoretical_critical_rates("criticality_by_p_" + to_string(LAYER_COUNT) + "layers.csv");

    // Commented out simulation execution
    eita = new double[circuit + 1];
    multi_dyn_spr();
    

    // --- Memory Deallocation ---

    // Free multilayer network memory
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        for (int i = 0; i < N; i++) {
            delete[] networks[layer][i];
        }
        delete[] networks[layer];
    }
    delete[] networks;

    // Free degrees array
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        delete[] degrees[layer];
    }
    delete[] degrees;

    // Free B_layers array
    for (int layer = 0; layer < LAYER_COUNT; layer++) {
        for (int i = 0; i < N; i++) {
            delete[] B_layers[layer][i];
        }
        delete[] B_layers[layer];
    }
    delete[] B_layers;

    // Free coupled network memory
    for (int i = 0; i < TOTAL_NODES; i++) {
        delete[] ABC[i];
        if (B_ABC[i] != nullptr) {
            delete[] B_ABC[i];
        }
    }
    delete[] ABC;
    delete[] B_ABC;
    delete[] D_ABC;

    // Free remaining arrays
    delete[] k_values;
    delete[] degree_sums;
    delete[] betweenness_layers;
    delete[] C_ABC;
    delete[] eita;

    // Free betweenness arrays
    delete[] bI_ABC;
    delete[] bE_ABC;

    free_ABC_path_storage();
    free_single_path_storage();

    // Final success message
    cout << "All results have been exported to files." << endl;
    cout << "\n=== Program Terminated Successfully ===" << endl;

    return 0;
}