// Last update: 2024/10/17
// Author: Ricard Arbat Carandell

// Master in Aerospace Engineering - Computational Engineering
// Universitat Politècnica de Catalunya (UPC)
// Overview: Convection-diffusion calculations

#define pass (void)0

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "mesh.cpp"
using namespace std;

struct Coefficients
{
    double C_L, C_D;
};

/**
 * Compares two matrices and returns the maximum error between them.
 *
 * @param v1 Matrix (Vector of vectors) to compare
 * @param v2 Matrix (Vector of vectors) to compare
 */
double calculateError(vector<vector<double>> v1, vector<vector<double>> v2)
{
    double max_error = delta_convergence;
    for (int i = 0; i < v1.size(); i++)
    {
        for (int j = 0; j < v1[i].size(); j++)
        {
            double difference = abs(v1[i][j] - v2[i][j]);
            if (difference > max_error)
            {
                max_error = difference;
            }
        }
    }
    return max_error;
}

/**
 * Computes the stream function of the mesh nodes. It uses the discretized stream function
 * equation to solve the stream values of the nodes. The function iterates until the error
 * is below a certain threshold.
 *
 * @param mesh Mesh matrix (vector of vectors)
 */
void computeTimeStep(vector<vector<vector<node>>> &mesh, int time_instant)
{
    // TODO: WORKING ON THIS!

    double an, ae, as, aw, ap, b_p;
    double dPE, dPe, dEe, dPS, dPs, dSs, dPW, dPw, dWw, dPN, dPn, dNn;
    double gauss_seidel;
    double error = 1;
    int cont = 0;

    vector<vector<double>> last_value(N, vector<double>(M, 0));
    vector<vector<double>> next_value(N, vector<double>(M, 0));

    // Convergence loop for the variable phi for each time step
    while (error > delta_convergence)
    {

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {

                // BOUNDARY CONDITIONS
                if ((i == 0) || (i == N - 1) || (j = M - 1)) // Left, top and right BC nodes
                {
                    next_value[i][j] = 1 - tanh(10);
                }

                else if (j == 0) // Bottom BC nodes
                {
                    if (i < N / 2)
                    {
                        next_value[i][j] = 1 + tanh(10 * (2 * mesh[t][i][j].x + 1));
                    }
                    else
                    {
                        next_value[i][j] = -1;
                    }
                }
                else // Internal nodes
                {

                    // Discretization of the convective-diffusive

                    ap = ae + aw + an + as + initial_density * dx * dy * 1 / delta_t - source_term_phi * dx * dy * 1;
                    bp = initial_density * dx * dy * (1 / delta_t) * mesh[timeinstant - 1][i][j].phi + source_term_consti

                                                                                                           // Clculating next value
                                                                                                           next_value[i][j] = (ae * phi_e + aw * phi_w + an * phi_n + as * phi_s + bp) / ap;
                }

                gauss_seidel = (last_value[i + 1][j] * ae + last_value[i - 1][j] * aw + last_value[i][j + 1] * an + last_value[i][j - 1] * as + b_p) / ap;

                next_value[i][j] = last_value[i][j] + relaxation_factor * (gauss_seidel - last_value[i][j]);
            }
        }
    }
    error = calculateError(next_value, last_value);

    cont++;

    if (cont == 100)
    {
        printf("Error = %f\n", error);
        cont = 0;
    }
}

// Update the mesh values with the new value
for (int i = 0; i < N; i++)
{
    for (int j = 0; j < M; j++)
    {
        mesh_t[i][j].phi = next_value[i][j];
    }
}
}
