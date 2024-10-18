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

double cds(double P);
double hds(double P);
double eds(double P);
double pds(double P);
double calculate_error(vector<vector<double>> v1, vector<vector<double>> v2);
void fill_values_at_t(vector<vector<double>> &values, double value);

/**
 * Computes the stream function of the mesh nodes. It uses the discretized stream function
 * equation to solve the stream values of the nodes. The function iterates until the error
 * is below a certain threshold.
 *
 * @param mesh Mesh matrix (vector of vectors)
 */
void compute_diffusive_convective(vector<vector<vector<node>>> &mesh, double gamma, string scheme)
{
    // Initializing variables for each iteration
    double aE, aW, aN, aS, aP, aP0, b;
    double De, Dw, Dn, Ds;
    double Fe, Fw, Fn, Fs;
    double ue, uw, vn, vs;
    double Pe, Pw, Pn, Ps;
    double Ae, Aw, An, As;
    double phi_e, phi_w, phi_n, phi_s;
    double gauss_seidel;
    double error = 10;
    int cont = 0;
    bool as_bc = false;

    vector<vector<double>> last_value(N, vector<double>(M, 0));
    vector<vector<double>> next_value(N, vector<double>(M, 0));

    fill_values_at_t(last_value, initial_value);

    // Time loop
    for (int t = 1; t < time_steps; t++)
    {
        error = 10;
        // Convergence loop for the variable phi for each time step
        cout << "Time step: " << t << endl;
        while (error > delta_convergence)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    // Discretization of the convective-diffusive
                    // Velocities at the faces (According to the Smith-Hutton flow + dx or dy effect)
                    ue = 2 * mesh[t][i][j].y * (1 - pow((mesh[t][i][j].x + 0.5 * dx), 2));
                    uw = 2 * mesh[t][i][j].y * (1 - pow((mesh[t][i][j].x - 0.5 * dx), 2));
                    vn = -2 * mesh[t][i][j].x * (1 - pow((mesh[t][i][j].y + 0.5 * dy), 2));
                    vs = -2 * mesh[t][i][j].x * (1 - pow((mesh[t][i][j].y - 0.5 * dy), 2));

                    // Diffusive terms
                    De = gamma * dy / dx;
                    Dw = gamma * dy / dx;
                    Dn = gamma * dx / dy;
                    Ds = gamma * dx / dy;

                    // Mass flows
                    Fe = rho * ue * dy;
                    Fw = rho * uw * dy;
                    Fn = rho * vn * dx;
                    Fs = rho * vs * dx;

                    // Peclet numbers
                    Pe = Fe / De;
                    Pw = Fw / Dw;
                    Pn = Fn / Dn;
                    Ps = Fs / Ds;

                    // Schemes
                    if (scheme == "UDS")
                    {
                        Ae = Pe;
                        Aw = Pw;
                        An = Pn;
                        As = Ps;
                    }
                    else if (scheme == "CDS")
                    {
                        Ae = cds(Pe);
                        Aw = cds(Pw);
                        An = cds(Pn);
                        As = cds(Ps);
                    }
                    else if (scheme == "HDS")
                    {
                        Ae = hds(Pe);
                        Aw = hds(Pw);
                        An = hds(Pn);
                        As = hds(Ps);
                    }
                    else if (scheme == "PDS")
                    {
                        Ae = pds(Pe);
                        Aw = pds(Pw);
                        An = pds(Pn);
                        As = pds(Ps);
                    }
                    else if (scheme == "EDS")
                    {
                        Ae = eds(Pe);
                        Aw = eds(Pw);
                        An = eds(Pn);
                        As = eds(Ps);
                    }

                    // Coefficients
                    aE = De * Ae + max(-Fe, 0.0);
                    aW = Dw * Aw + max(Fw, 0.0);
                    aN = Dn * An + max(-Fn, 0.0);
                    aS = Ds * As + max(Fs, 0.0);
                    aP0 = rho * dx * dy / delta_t;
                    aP = aE + aW + aN + aS + aP0 - S_phi * dx * dy;
                    b = aP0 * mesh[t - 1][i][j].phi + S_c * dx * dy;

                    // BOUNDARY CONDITIONS
                    if (i == 0) // Left BC
                    {
                        next_value[i][j] = 1 - tanh(10);
                    }
                    else if (i == N - 1) // Top BC
                    {
                        next_value[i][j] = 1 - tanh(10);
                    }
                    else if (j == M - 1) // Right BC
                    {
                        next_value[i][j] = 1 - tanh(10);
                    }
                    else if (j == 0) // Bottom inlet BC nodes
                    {
                        if (i < N / 2)
                        {
                            next_value[i][j] = 1 + tanh(10 * (2 * mesh[t][i][j].x + 1));
                        }
                        else
                        {
                            next_value[i][j] = (last_value[i + 1][j] * aE + last_value[i - 1][j] * aW + last_value[i][j + 1] * aN + last_value[i][j] * aS + b) / aP;
                        }
                    }
                    else // Internal nodes or outlet nodes
                    {

                        next_value[i][j] = (last_value[i + 1][j] * aE + last_value[i - 1][j] * aW + last_value[i][j + 1] * aN + last_value[i][j - 1] * aS + b) / aP;
                    }
                }
            }

            // Print matrices

            // Calculate the error
            error = calculate_error(next_value, last_value);
            printf("Error = %f, timestep %i \n ", error, t);

            // Update the last value with the new value
            last_value = next_value;
        }

        // Update the mesh values with the new value
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                mesh[t][i][j].phi = next_value[i][j];
            }
        }
    }
}

/**
 * Compares two matrices and returns the maximum error between them.
 *
 * @param v1 Matrix (Vector of vectors) to compare
 * @param v2 Matrix (Vector of vectors) to compare
 */
double calculate_error(vector<vector<double>> v1, vector<vector<double>> v2)
{
    double max_error = delta_convergence;
    for (int i = 1; i < v1.size(); i++)
    {
        for (int j = 1; j < v1[i].size(); j++)
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
 * Fills the initial values of the phi temporal mesh
 *
 * @param values Matrix (Vector of vectors) containing the next stream values
 * @param value Value to set for all the nodes
 */
void fill_values_at_t(vector<vector<double>> &values, double value)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            values[i][j] = value;
        }
    }
}

// Schemes
double cds(double P)
{
    return 1 - 0.5 * abs(P);
}

double hds(double P)
{
    return max(1 - 0.5 * abs(P), 0.0);
}

double eds(double P)
{
    return abs(P) / (exp(abs(P)) - 1);
}

double pds(double P)
{
    return max(pow(1 - 0.5 * abs(P), 5), 0.0);
}