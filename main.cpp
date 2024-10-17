// Last update: 2024/10/06
// Author: Ricard Arbat Carandell

// Master in Aerospace Engineering - Computational Engineering
// Universitat Politècnica de Catalunya (UPC) - BarcelonaTech
// Overview: Potential flow solution in a channel around a cilinder using the stream function formulation.

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

    // Main mesh
    Parameters p;
    vector<vector<node>> mesh(N, vector<node>(M));
    buildMesh(mesh, p); // creating the mesh
    setStream(mesh, p); // setting initial values
    setRho(mesh, p);

    // Compute stream function
    auto start = high_resolution_clock::now();
    computeStream(mesh, p); // stream solver
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    calculateVelocity(mesh, p);                                 // calculating velocity
    calculateCp(mesh, p);                                       // calculating pressure coeficient distribution
    calculatePressureTemperature(mesh, p);                      // calculating pressure and temperature
    double circulation = calculateCylinderCirculation(mesh, p); // calculate circulation around cilinder
    struct Coefficients c = cylinderForces(mesh, p);            // calculate forces around cilinder
    string name = p.folder + p.mesh_number + "_output.csv";
    exportData(mesh, name); // export data to file output.dat

    name = p.folder + p.mesh_number + "_results.csv";
    ofstream outfile_r(name);
    outfile_r << "L,H,R,N,M,Cl,Cd,Circ,t" << endl;
    outfile_r << p.L << "," << p.H << "," << p.cylinder_r << "," << N << "," << M << "," << c.C_L << "," << c.C_D << "," << circulation << "," << duration.count() / 1000 << endl;

    // Compute analytic stream function
    vector<vector<node>> analytic_mesh(N, vector<node>(M));
    buildMesh(analytic_mesh, p);             // creating the mesh
    computeAnalyticStream(analytic_mesh, p); // stream solver
    calculateVelocity(mesh, p);              // calculating velocity
    calculateCp(mesh, p);                    // calculating pressure coeficient distribution
    calculatePressureTemperature(mesh, p);   // calculating pressure and temperature
    name = p.folder + p.mesh_number + "_analytic_output.csv";
    exportData(analytic_mesh, name); // export data to file output.dat

    vector<vector<node>> error_mesh = analyticError(mesh, analytic_mesh); // calculating error
    name = p.folder + p.mesh_number + "_error_output.csv";
    exportData(error_mesh, name);

    // Print final results
    cout << "### POTENTIAL FLOW RESULTS ###" << endl;
    cout << "C_L = " << c.C_L << endl;
    cout << "C_D = " << c.C_D << endl;
    cout << "Cylinder circulation = " << circulation << endl;
    cout << "Equivalent rotation speed = " << circulation / (2 * M_PI * p.cylinder_r * p.cylinder_r) << endl;
    cout << "Computation time = " << duration.count() / 1000 << "ms" << endl;
    return 0;
}
