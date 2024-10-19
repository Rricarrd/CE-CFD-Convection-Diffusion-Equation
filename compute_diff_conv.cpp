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
void evaluate_time_step_smith_hutton(int t, vector<vector<vector<node>>> &mesh, double gamma, string scheme);
void evaluate_time_step_diagonal(int t, vector<vector<vector<node>>> &mesh, double gamma, string scheme);

/**
 * Computes the stream function of the mesh nodes. It uses the discretized stream function
 * equation to solve the stream values of the nodes. The function iterates until the error
 * is below a certain threshold.
 *
 * @param mesh Mesh matrix (vector of vectors)
 */
void compute_diffusive_convective(vector<vector<vector<node>>> &mesh, double gamma, string scheme, string type) // type of problem)
{

    // Time loop
    for (int t = 1; t < time_steps; t++)
    {
        if (type == "smith-hutton")
        {
            evaluate_time_step_smith_hutton(t, mesh, gamma, scheme);
        }
        else if (type == "diagonal")
        {
            evaluate_time_step_diagonal(t, mesh, gamma, scheme);
        }
    }
}

/**
 * Evaluates the time step of the mesh
 *
 * @param t Time instant of the simulation
 * @param mesh Mesh matrix (vector of vectors)
 * @param convergence_criteria Maximum delta for the error
 * @param Pe Peclet number
 */

void evaluate_time_step_smith_hutton(int t, vector<vector<vector<node>>> &mesh, double gamma, string scheme)
{
    // Initializing variables for each iteration
    double aE, aW, aN, aS, aP, aP0, b;
    double De, Dw, Dn, Ds;
    double Fe, Fw, Fn, Fs;
    double ue, uw, vn, vs;
    double Pe, Pw, Pn, Ps;
    double Ae, Aw, An, As;
    double phiE, phiW, phiN, phiS;
    double error = 10;
    int cont = 0;

    vector<vector<double>> last_phi(N, vector<double>(M, 0));
    vector<vector<double>> next_phi(N, vector<double>(M, 0));

    // Copy values from last timestep
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            last_phi[i][j] = mesh[t - 1][i][j].phi; // copy phi from previous time step
        }
    }

    while (error > delta_convergence)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                ue = 2 * mesh[t][i][j].y * (1 - ((mesh[t][i][j].x + 0.5 * dx) * (mesh[t][i][j].x + 0.5 * dx)));
                uw = 2 * mesh[t][i][j].y * (1 - ((mesh[t][i][j].x - 0.5 * dx) * (mesh[t][i][j].x - 0.5 * dx)));
                vn = -2 * mesh[t][i][j].x * (1 - ((mesh[t][i][j].y + 0.5 * dy) * (mesh[t][i][j].y + 0.5 * dy)));
                vs = -2 * mesh[t][i][j].x * (1 - ((mesh[t][i][j].y - 0.5 * dy) * (mesh[t][i][j].y - 0.5 * dy)));

                De = gamma * (dy / dx);
                Dw = gamma * (dy / dx);
                Dn = gamma * (dx / dy);
                Ds = gamma * (dx / dy);

                Fe = rho * dy * ue;
                Fw = rho * dy * uw;
                Fn = rho * dx * vn;
                Fs = rho * dx * vs;

                // Peclet numbers
                Pe = Fe / De;
                Pw = Fw / Dw;
                Pn = Fn / Dn;
                Ps = Fs / Ds;

                // Schemes
                if (scheme == "UDS")
                {
                    Ae = 1;
                    Aw = 1;
                    An = 1;
                    As = 1;
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

                aE = De * Ae + max(-Fe, 0.0);
                aW = Dw * Aw + max(Fw, 0.0);
                aN = Dn * An + max(-Fn, 0.0);
                aS = Ds * As + max(Fs, 0.0);

                aP = aE + aW + aN + aS + (rho * dx * dy / delta_t) - (S_phi * dx * dy);

                b = ((rho * dx * dy * mesh[t - 1][i][j].phi) / delta_t) + (S_c * dx * dy);

                if (i == 0)
                {
                    next_phi[i][j] = 1 - tanh(10); // bc
                }
                else if (i == N - 1)
                {
                    next_phi[i][j] = 1 - tanh(10); // bc
                }
                else if (j == 0)
                {
                    if (i < N / 2) // inlet
                    {
                        next_phi[i][j] = 1 + tanh(10 * (2 * mesh[t][i][j].x + 1)); // bc
                    }
                    else // outlet
                    {
                        phiE = last_phi[i + 1][j];
                        phiW = last_phi[i - 1][j];
                        phiN = last_phi[i][j + 1];
                        phiS = last_phi[i][j]; // bc
                        next_phi[i][j] = (aE * phiE + aW * phiW + aN * phiN + aS * phiS + b) / aP;
                    }
                }
                else if (j == M - 1)
                {
                    next_phi[i][j] = 1 - tanh(10); // bc
                }
                else
                {
                    phiE = last_phi[i + 1][j];
                    phiW = last_phi[i - 1][j];
                    phiN = last_phi[i][j + 1];
                    phiS = last_phi[i][j - 1];
                    next_phi[i][j] = (aE * phiE + aW * phiW + aN * phiN + aS * phiS + b) / aP;
                }
            }
        }

        error = calculate_error(next_phi, last_phi);

        if (cont % 100 == 0)
            cout << "Time step: " << t << " Error: " << error << endl;
        cont++;

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                last_phi[i][j] = next_phi[i][j];
            }
        }
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            mesh[t][i][j].phi = next_phi[i][j];
        }
    }
}

void evaluate_time_step_diagonal(int t, vector<vector<vector<node>>> &mesh, double gamma, string scheme)
{
    // Initializing variables for each iteration
    double aE, aW, aN, aS, aP, aP0, b;
    double De, Dw, Dn, Ds;
    double Fe, Fw, Fn, Fs;
    double ue, uw, vn, vs;
    double Pe, Pw, Pn, Ps;
    double Ae, Aw, An, As;
    double phiE, phiW, phiN, phiS;
    double error = 10;
    int cont = 0;

    vector<vector<double>> last_phi(N, vector<double>(M, 0));
    vector<vector<double>> next_phi(N, vector<double>(M, 0));

    // Copy values from last timestep
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            last_phi[i][j] = mesh[t - 1][i][j].phi; // copy phi from previous time step
        }
    }

    while (error > delta_convergence)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                ue = u0 * cos(0.785398);
                uw = u0 * cos(0.785398);
                vn = u0 * sin(0.785398);
                vs = u0 * sin(0.785398);

                De = gamma * (dy / dx);
                Dw = gamma * (dy / dx);
                Dn = gamma * (dx / dy);
                Ds = gamma * (dx / dy);

                Fe = rho * dy * ue;
                Fw = rho * dy * uw;
                Fn = rho * dx * vn;
                Fs = rho * dx * vs;

                // Peclet numbers
                Pe = Fe / De;
                Pw = Fw / Dw;
                Pn = Fn / Dn;
                Ps = Fs / Ds;

                // Schemes
                if (scheme == "UDS")
                {
                    Ae = 1;
                    Aw = 1;
                    An = 1;
                    As = 1;
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

                aE = De * Ae + max(-Fe, 0.0);
                aW = Dw * Aw + max(Fw, 0.0);
                aN = Dn * An + max(-Fn, 0.0);
                aS = Ds * As + max(Fs, 0.0);

                aP = aE + aW + aN + aS + (rho * dx * dy / delta_t) - (S_phi * dx * dy);

                b = ((rho * dx * dy * mesh[t - 1][i][j].phi) / delta_t) + (S_c * dx * dy);

                if (i == 0) // Left BC
                {
                    next_phi[i][j] = phi_l;
                }
                else if (i == N - 1) // Right BC
                {
                    next_phi[i][j] = phi_h;
                }
                else if (j == M - 1) // Top BC
                {
                    next_phi[i][j] = phi_l;
                }
                else if (j == 0) // Bottom inlet BC nodes
                {
                    next_phi[i][j] = phi_h;
                }
                else
                {
                    phiE = last_phi[i + 1][j];
                    phiW = last_phi[i - 1][j];
                    phiN = last_phi[i][j + 1];
                    phiS = last_phi[i][j - 1];
                    next_phi[i][j] = (aE * phiE + aW * phiW + aN * phiN + aS * phiS + b) / aP;
                }
            }
        }

        error = calculate_error(next_phi, last_phi);

        if (cont % 100 == 0)
            cout << "Time step: " << t << " Error: " << error << endl;
        cont++;

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                last_phi[i][j] = next_phi[i][j];
            }
        }
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            mesh[t][i][j].phi = next_phi[i][j];
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

double pds(double P)
{
    return max(pow(1 - 0.5 * abs(P), 5), 0.0);
}

double eds(double P)
{
    double val = abs(P) / (exp(abs(P)) - 1) + 0.0000000000001;
    return val;
}