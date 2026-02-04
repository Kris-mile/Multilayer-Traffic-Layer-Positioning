Multilayer Network Traffic Dynamics Simulator
A comprehensive simulation system for studying information transmission in multilayer complex networks with ER, BA, and RGG topologies.

Overview
This project implements a multilayer coupled network system to investigate information transmission performance under different network topologies, coupling configurations, and parameter settings. The system focuses on betweenness centrality and critical packet generation rate analysis.

Key Features
- Three Network Topologies: ER (Erdős–Rényi), BA (Barabási–Albert), RGG (Random Geometric Graph)
- Three Coupling Configurations: LDC (Low-Degree Central), MDC (Medium-Degree Central), HDC (High-Degree Central)
- Comprehensive Analysis: Betweenness centrality, queue length analysis, critical rate calculation
- Performance Optimized: O(E log V) shortest path algorithms, parallel RGG generation

Quick Start

Prerequisites
- C++11 or higher
- OpenMP (optional, for parallel acceleration)
- CMake (optional)

Compilation
```bash
 Clone the repository
git clone https://github.com/Kris-mile/Multilayer-Traffic-Layer-Positioning.git
cd multilayer-network-traffic

 Compile with optimization
g++ -std=c++11 -O3 -fopenmp main.cpp -o network_simulator
