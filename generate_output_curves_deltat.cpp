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
    double Pe = 1000;                                // peclet numbers
    string schemes[] = {"UDS", "CDS", "HDS", "PDS"}; // scheme to be used
    string type = "smith-hutton";
    double delta_t[] = {0.0001}; // type of problem
    double gamma = rho / Pe;
    for (int k = 0; k < 4; k++)
    {
        ofstream outfile("deltat/" + file_name(Pe, "COMPARISON_delta", "SCHEMES", type, delta_t[k]));

        for (int i = 0; i < 5; i++)
        {
            // Main mesh
            vector<vector<vector<node>>> mesh(time_steps, vector<vector<node>>(N, vector<node>(M)));
            build_mesh(mesh, type); // creating the mesh

            set_mesh_value(mesh[0], 1, "phi");

            // Compute stream function
            auto start = high_resolution_clock::now();
            compute_diffusive_convective(mesh, gamma, delta_t[k], schemes[i], type); // stream solver
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);

            // Export data    string path = folder + "/" + filename;
            export_data_at_outlet(mesh, schemes[i], outfile);
        }
        outfile.close();
    }
    return 0;
}