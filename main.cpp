// Last update: 2024/10/17
// Author: Ricard Arbat Carandell

// Master in Aerospace Engineering - Computational Engineering
// Universitat Politècnica de Catalunya (UPC) - BarcelonaTech
// Overview: Solution of the Smth-Hutton problem using the convection-diffusuion equations.

// Libraries
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "compute_potential.cpp"
#include <chrono>
using namespace std;
using namespace std::chrono;

/**
 * Main function of the program. It initializes the mesh, computes the stream
 * function and calculates the forces around the cilinder.
 */
int main(void)
{

    double Pe = 100; // Peclet number

    // Main mesh
    vector<vector<vector<node>>> mesh(time_steps, vector<vector<node>>(N, vector<node>(M)));
    buildMesh(mesh); // creating the mesh

    // Compute stream function
    auto start = high_resolution_clock::now();
    computeStream(mesh); // stream solver
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    calculateVelocity(mesh); // calculating velocity
    calculateCp(mesh);       // calculating pressure coeficient distribution

    // Print final results
    cout << "### CONVECTION DIFFUSION FLOW RESULTS ###" << endl;
    cout << "Computation time = " << duration.count() / 1000 << "ms" << endl;
    return 0;
}
