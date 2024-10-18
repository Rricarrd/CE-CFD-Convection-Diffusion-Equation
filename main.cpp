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
#include "compute_diff_conv.cpp"
#include <chrono>
using namespace std;
using namespace std::chrono;

/**
 * Main function of the program. It initializes the mesh, computes the stream
 * function and calculates the forces around the cilinder.
 */
int main(void)
{

    double rho_to_gamma = 1;           // density to gamma ratio
    double gamma = rho / rho_to_gamma; // diffusion coefficient
    string scheme = "PDS";             // scheme to be used

    // Main mesh
    vector<vector<vector<node>>> mesh(time_steps, vector<vector<node>>(N, vector<node>(M)));
    build_mesh(mesh); // creating the mesh
    set_smith_hutton_problem(mesh);
    set_mesh_value(mesh[0], 1, "phi");

    // Compute stream function
    auto start = high_resolution_clock::now();
    compute_diffusive_convective(mesh, gamma, scheme); // stream solver
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    // Export data
    export_data(mesh, file_name(time_steps, "output"));

    // Print final results
    cout << "### CONVECTION DIFFUSION FLOW RESULTS ###" << endl;
    cout << "Computation time = " << duration.count() / 1000 << "ms" << endl;
    return 0;
}
