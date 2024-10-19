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

    double Pe = 1000000;      // peclet number
    double gamma = rho / Pe;  // diffusion coefficient
    string scheme = "PDS";    // scheme to be used
    string type = "diagonal"; // type of problem

    if (type == "smith-hutton")
    {
        N = 2 * M;
    }
    else if (type == "diagonal")
    {
        N = M;
    }

    // Main mesh
    vector<vector<vector<node>>> mesh(time_steps, vector<vector<node>>(N, vector<node>(M)));
    build_mesh(mesh, type); // creating the mesh

    set_mesh_value(mesh[0], initial_phi, "phi");

    // Compute stream function
    auto start = high_resolution_clock::now();
    compute_diffusive_convective(mesh, gamma, scheme, type); // stream solver
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    // Export data
    export_data(mesh, file_name(Pe, scheme, "output", type), "output", "last");

    // Print final results
    cout << "### CONVECTION DIFFUSION FLOW RESULTS ###" << endl;
    cout << "Computation time = " << duration.count() / 1000 << "ms" << endl;
    return 0;
}
