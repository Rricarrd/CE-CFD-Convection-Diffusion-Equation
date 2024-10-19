// Author: Ricard Arbat Carandell

// Master in Aerospace Engineering - Computational Engineering
// Universitat Politècnica de Catalunya (UPC) - BarcelonaTech
// Overview: // TODO

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
    double Pe[] = {10, 1000, 1000000};                        // peclet numbers
    double gamma[] = {rho / Pe[0], rho / Pe[1], rho / Pe[2]}; // diffusion coefficient
    string scheme = "HDS";                                    // scheme to be used
    string type = "smith-hutton";                             // type of problem

    for (int i = 0; i < 3; i++)
    {
        // Main mesh
        vector<vector<vector<node>>> mesh(time_steps, vector<vector<node>>(N, vector<node>(M)));
        build_mesh(mesh, type); // creating the mesh

        set_mesh_value(mesh[0], 1, "phi");

        // Compute stream function
        auto start = high_resolution_clock::now();
        compute_diffusive_convective(mesh, gamma[i], scheme, type); // stream solver
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        // Export data
        export_data(mesh, file_name(Pe[i], scheme, "output", type), "smith_hutton", "last");

        // Print final results
        cout << "### CONVECTION DIFFUSION FLOW RESULTS FOR" << Pe[i] << "###" << endl;
        cout << "Computation time = " << duration.count() / 1000000 << "s" << endl;
    }

    return 0;
}
