/**
 * @file calculate_betweenness.h
 * @brief Calculates the intra-layer and inter-layer effective betweenness centrality using an optimized Brandes algorithm.
 */

 extern double* bI_ABC;  // Intra-layer betweenness
 extern double* bE_ABC;  // Inter-layer betweenness

 // Initialize betweenness arrays
 void initialize_betweenness_arrays() {
     bI_ABC = new double[TOTAL_NODES]();  // Modified to TOTAL_NODES
     bE_ABC = new double[TOTAL_NODES]();
 }

 // Recursively calculate intra-layer betweenness
 void digui_internal(double pro, int next, int end_node, int layer_index) {
     int w;
     int local_next = next % N;
     int local_end = end_node % N;

     // Check if path exists
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

 // Recursively calculate inter-layer betweenness
 void digui_external(double pro, int next, int end_node) {
     int w;

     // Check if path exists
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

 // Calculate intra-layer betweenness
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

 // Calculate inter-layer betweenness
 void calculate_external_betweenness() {
     // Iterate through all pairs of nodes in different layers
     for (int s = 0; s < TOTAL_NODES; s++) {  // Modified to TOTAL_NODES
         int s_layer = s / N;

         for (int t = 0; t < TOTAL_NODES; t++) {  // Modified to TOTAL_NODES
             int t_layer = t / N;

             // Only consider pairs from different layers
             if (s_layer == t_layer) continue;

             if (path_num_ABC[s][t] > 0) {
                 double pro = 1.0;
                 digui_external(pro, s, t);
             }
         }
     }
 }

 // Calculate enhanced betweenness
 void cal_betweenness_enhanced() {
     initialize_betweenness_arrays();

     // Calculate intra-layer betweenness for each layer
     for (int layer = 0; layer < LAYER_COUNT; layer++) {  // Modified to dynamic layer count
         calculate_internal_betweenness(layer);
     }

     // Calculate inter-layer betweenness
     calculate_external_betweenness();

     cout << "Betweenness calculation completed!" << endl;
 }

 // Output betweenness results
 void output_betweenness_results(const string& filename) {
     ofstream out_file(filename);
     if (!out_file) {
         cerr << "Failed to open file: " << filename << endl;
         return;
     }

     out_file << "Node_ID,Layer,Internal_Betweenness,External_Betweenness\n";

     for (int i = 0; i < TOTAL_NODES; i++) {  // Modified to TOTAL_NODES
         int layer = i / N;
         int node_in_layer = i % N;

         out_file << i << "," << layer << ","
             << bI_ABC[i] << "," << bE_ABC[i] << "\n";
     }

     out_file.close();
     cout << "Betweenness calculation results saved to: " << filename << endl;
 }
 
 // Global betweenness arrays (externally defined in main.cpp)
 /*extern double* bI_ABC;  // Intra-layer betweenness
 extern double* bE_ABC;  // Inter-layer betweenness
 extern int*** networks; // Single-layer networks
 extern int** ABC;       // Coupled network
 extern int** degrees;
 extern int* D_ABC;      // Coupled degrees
 */
 /**
  * @brief Initializes the betweenness arrays to zero.
  */
  /*void initialize_betweenness_arrays() {
      std::fill(bI_ABC, bI_ABC + TOTAL_NODES, 0.0);
      std::fill(bE_ABC, bE_ABC + TOTAL_NODES, 0.0);
  }*/

  // ============================================================
  // Core Algorithm: Brandes Algorithm (Single-layer version - for computing bI)
  // ============================================================

  /**
   * @brief Computes the intra-layer betweenness (bI) using Brandes' algorithm for a specific layer.
   * @param layer_idx The index of the layer to process.
   */
   /*void run_brandes_single_layer(int layer_idx) {
       int offset = layer_idx * N;

       // Parallel computation: each thread processes a subset of source nodes s
   #pragma omp parallel for schedule(dynamic)
       for (int s_local = 0; s_local < N; s_local++) {
           // Thread-local variables to avoid frequent memory allocation
           vector<int> sigma(N, 0);      // Shortest path count
           vector<int> d(N, -1);         // Distances
           vector<double> delta(N, 0.0); // Dependency values
           vector<vector<int>> P(N);     // Predecessor node list
           queue<int> Q;                 // BFS queue
           stack<int> S;                 // Visitation order stack

           int s = s_local;

           // 1. BFS Phase
           sigma[s] = 1;
           d[s] = 0;
           Q.push(s);

           while (!Q.empty()) {
               int v = Q.front();
               Q.pop();
               S.push(v);

               // Traverse neighbors (using networks[layer_idx][v] adjacency matrix)
               // Optimization: adjacency list would be better, traversing matrix here
               for (int w = 0; w < N; w++) {
                   if (networks[layer_idx][v][w] == 0) continue; // No edge

                   // Discover new node
                   if (d[w] < 0) {
                       Q.push(w);
                       d[w] = d[v] + 1;
                   }

                   // Accumulate shortest path count
                   if (d[w] == d[v] + 1) {
                       sigma[w] += sigma[v];
                       P[w].push_back(v);
                   }
               }
           }

           // 2. Dependency accumulation (backtracking) phase
           while (!S.empty()) {
               int w = S.top();
               S.pop();

               for (int v : P[w]) {
                   // Brandes formula: delta[v] += (sigma[v]/sigma[w]) * (1 + delta[w])
                   // Note: Calculating "Internal Betweenness" here, assuming this is an isolated layer
                   // w contributes to the path only when w != s
                   delta[v] += ((double)sigma[v] / sigma[w]) * (1.0 + delta[w]);
               }

               // Accumulate to global betweenness (excluding source node)
               if (w != s) {
                   // Use atomic operation to prevent multi-threading conflicts
   #pragma omp atomic
                   bI_ABC[offset + w] += delta[w];
               }
           }
       }
   }
   */
   // ============================================================
   // Core Algorithm: Brandes Algorithm (Coupled version - for computing bE)
   // ============================================================

   /**
    * @brief Computes the inter-layer betweenness (bE) using Brandes' algorithm on the fully coupled network.
    */
    /*void run_brandes_coupled_network() {
        // Preprocessing: Convert adjacency matrix ABC to adjacency list to accelerate BFS
        // This step is crucial for RGG because matrix traversal is too slow
        vector<vector<int>> adj_list(TOTAL_NODES);
        for (int i = 0; i < TOTAL_NODES; i++) {
            // If D_ABC is accurate, we can reserve
            adj_list[i].reserve(20);
            for (int j = 0; j < TOTAL_NODES; j++) {
                if (ABC[i][j]) {
                    adj_list[i].push_back(j);
                }
            }
        }

        // Parallel computation
    #pragma omp parallel for schedule(dynamic)
        for (int s = 0; s < TOTAL_NODES; s++) {
            // Thread-local variables
            vector<long long> sigma(TOTAL_NODES, 0); // Note: Total graph path count might overflow int, using long long
            vector<int> d(TOTAL_NODES, -1);
            vector<double> delta(TOTAL_NODES, 0.0);
            vector<vector<int>> P(TOTAL_NODES);
            queue<int> Q;
            stack<int> S;

            int s_layer = s / N;

            // 1. BFS
            sigma[s] = 1;
            d[s] = 0;
            Q.push(s);

            while (!Q.empty()) {
                int v = Q.front();
                Q.pop();
                S.push(v);

                // Traverse using adjacency list
                for (int w : adj_list[v]) {
                    if (d[w] < 0) {
                        Q.push(w);
                        d[w] = d[v] + 1;
                    }
                    if (d[w] == d[v] + 1) {
                        sigma[w] += sigma[v];
                        P[w].push_back(v);
                    }
                }
            }

            // 2. Backtracking
            while (!S.empty()) {
                int w = S.top();
                S.pop();

                // Key modification: To match the "External Betweenness" definition
                // bE counts the number of times paths (where source s and destination t are in different layers) pass through
                // (1 + delta[w]) in the original Brandes formula means w provides 1 unit of flow as a target node
                // We modify it to: if layer(s) != layer(w), w provides 1 flow as a valid target, else 0

                int w_layer = w / N;
                double flow_contribution = (s_layer != w_layer) ? 1.0 : 0.0;

                // Add dependencies from subsequent nodes
                double total_dependency = flow_contribution + delta[w];

                for (int v : P[w]) {
                    delta[v] += ((double)sigma[v] / sigma[w]) * total_dependency;
                }

                // Accumulate to global
                if (w != s) {
    #pragma omp atomic
                    bE_ABC[w] += delta[w];
                }
            }

        }
        printf("\n");
    }
    */
    /**
     * @brief Main execution function to calculate the enhanced betweenness centrality.
     */
     /*void cal_betweenness_enhanced() {

         // 1. Clear/Initialize
         initialize_betweenness_arrays();

         // 2. Calculate intra-layer betweenness (layers are independent)
         for (int layer = 0; layer < LAYER_COUNT; layer++) {
             run_brandes_single_layer(layer);
         }

         // 3. Calculate inter-layer betweenness (bE - cross-layer traffic only)
         run_brandes_coupled_network();

     }
     */
     /**
      * @brief Outputs the calculated betweenness results to a CSV file.
      * @param filename The name of the output file.
      */
      /*void output_betweenness_results(const string& filename) {
          ofstream out_file(filename);
          if (!out_file) {
              cerr << "Failed to open file: " << filename << endl;
              return;
          }

          out_file << "Node_ID,Layer,Internal_Betweenness,External_Betweenness\n";

          for (int i = 0; i < TOTAL_NODES; i++) {
              int layer = i / N;
              // The betweenness calculated by Brandes algorithm is usually doubled (undirected graph s-t and t-s are counted twice),
              // and your original recursive algorithm calculated 1.0 flow for each pair (s,t).
              // Brandes' result is completely equivalent to traversing all paths.
              // If normalization is needed, you can divide by 2 or (N-1)(N-2) here.
              // Keeping the original cumulative value here.

              out_file << i << "," << layer << ","
                  << bI_ABC[i] << "," << bE_ABC[i] << "\n";
          }

          out_file.close();
          cout << "Betweenness calculation results saved to: " << filename << endl;
      }*/