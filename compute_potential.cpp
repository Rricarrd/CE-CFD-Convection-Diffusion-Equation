// Last update: 2024/10/06
// Author: Ricard Arbat Carandell

// Master in Aerospace Engineering - Computational Engineering
// Universitat Politècnica de Catalunya (UPC)
// Overview: Stream computing functions

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
double streamsError(vector<vector<double>> v1, vector<vector<double>> v2)
{
    double maxError = 0.0;
    for (int i = 0; i < v1.size(); i++)
    {
        for (int j = 0; j < v1[i].size(); j++)
        {
            double diff = abs(v1[i][j] - v2[i][j]);
            if (diff > maxError)
            {
                maxError = diff;
            }
        }
    }
    return maxError;
}

/**
 * Fills the initial stream values of the mesh nodes.
 *
 * @param initial_stream_value Matrix (Vector of vectors) containing the next stream values
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
void fillStream(vector<vector<double>> &initial_stream_value, vector<vector<node>> &mesh, Parameters p)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (mesh[i][j].is_solid == true)
            {
                initial_stream_value[i][j] = p.solid_stream;
            }
            else
            {
                initial_stream_value[i][j] = p.start_stream;
            }
        }
    }
}

/**
 * Computes the stream function of the mesh nodes. It uses the discretized stream function
 * equation to solve the stream values of the nodes. The function iterates until the error
 * is below a certain threshold.
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
void computeStream(vector<vector<node>> &mesh, Parameters p)
{
    double an, ae, as, aw, ap, b_p;
    double dPE, dPe, dEe, dPS, dPs, dSs, dPW, dPw, dWw, dPN, dPn, dNn;
    double error = p.initial_error;
    double gauss_seidel;

    vector<vector<double>> next_stream_value(N, vector<double>(M));
    vector<vector<double>> last_stream_value(N, vector<double>(M));

    fillStream(next_stream_value, mesh, p);
    int cont = 0;
    while (error > p.delta)
    {
        last_stream_value = next_stream_value;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {

                if (j == 0) // Bottom channel nodes
                {
                    next_stream_value[i][j] = 0;
                }
                else if (j == M - 1) // Top channel nodes
                {
                    next_stream_value[i][j] = p.v_in * p.H;
                }
                else if (i == 0) // Inlet nodes
                {
                    next_stream_value[i][j] = p.v_in * mesh[i][j].y;
                }
                else if (i == N - 1) // Outlet nodes
                {

                    next_stream_value[i][j] = last_stream_value[i - 1][j];
                }
                else // Internal nodes
                {
                    // East node
                    dPE = p.dx;
                    dPe = p.dx / 2;
                    dEe = p.dx / 2;
                    ae = (p.dx / ((dPe * mesh[i][j].rho + dEe * mesh[i + 1][j].rho) / p.rho_in)) * p.dy / dPE;

                    // South node
                    dPS = p.dy;
                    dPs = p.dy / 2;
                    dSs = p.dy / 2;
                    as = (dPS / ((dPs * mesh[i][j].rho + dSs * mesh[i][j - 1].rho) / p.rho_in)) * p.dx / dPS;

                    // West node
                    dPW = p.dx;
                    dPw = p.dx / 2;
                    dWw = p.dx / 2;
                    aw = (dPW / ((dPw * mesh[i][j].rho + dWw * mesh[i - 1][j].rho) / p.rho_in)) * p.dy / dPW;

                    // North node
                    dPN = p.dy;
                    dPn = p.dy / 2;
                    dNn = p.dy / 2;
                    an = (dPN / ((dPn * mesh[i][j].rho + dNn * mesh[i][j + 1].rho) / p.rho_in)) * p.dx / dPN;

                    // Discretization of the stream function equation
                    ap = an + ae + as + aw;
                    gauss_seidel = (last_stream_value[i + 1][j] * ae + last_stream_value[i - 1][j] * aw + last_stream_value[i][j + 1] * an + last_stream_value[i][j - 1] * as + b_p) / ap;
                    next_stream_value[i][j] = last_stream_value[i][j] + p.relaxation_factor * (gauss_seidel - last_stream_value[i][j]);
                }
            }
        }
        error = streamsError(next_stream_value, last_stream_value);
        cont++;
        if (cont == 100)
        {
            printf("Error = %f\n", error);
            cont = 0;
        }
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            mesh[i][j].stream = next_stream_value[i][j];
        }
    }
}

/**
 * Calculates the velocity of the mesh nodes. It uses the stream function values
 * to calculate the velocity components of the nodes.
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
void calculateVelocity(vector<vector<node>> &mesh, Parameters p)
{
    double vxn, vye, vxs, vyw, dPN, dPn, dNn, dPE, dPe, dEe, dPS, dPs, dSs, dPW, dPw, dWw;

    for (int i = 1; i < N - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {

            dPE = p.dy;
            dPe = p.dy / 2;
            dEe = p.dy / 2;
            vye = -(dPE / ((dPe * mesh[i][j].rho + dEe * mesh[i + 1][j].rho) / p.rho_in)) * ((mesh[i + 1][j].stream - mesh[i][j].stream) / dPE);

            dPS = p.dx;
            dPs = p.dx / 2;
            dSs = p.dx / 2;
            vxs = -(dPS / ((dPs * mesh[i][j].rho + dSs * mesh[i][j - 1].rho) / p.rho_in)) * ((mesh[i][j - 1].stream - mesh[i][j].stream) / dPS);

            dPW = p.dy;
            dPw = p.dy / 2;
            dWw = p.dy / 2;
            vyw = (dPW / ((dPw * mesh[i][j].rho + dWw * mesh[i - 1][j].rho) / p.rho_in)) * ((mesh[i - 1][j].stream - mesh[i][j].stream) / dPW);

            // North node
            dPN = p.dx;
            dPn = p.dx / 2;
            dNn = p.dx / 2;
            vxn = (dPN / ((dPn * mesh[i][j].rho + dNn * mesh[i][j + 1].rho) / p.rho_in)) * ((mesh[i][j + 1].stream - mesh[i][j].stream) / dPN);

            mesh[i][j].u = (vxn + vxs) / 2;
            mesh[i][j].v = (vye + vyw) / 2;
        }
    }
}

/**
 * Calculates the circulation around the cylinder. It uses the velocity values
 * to calculate the circulation of each cell around the cylinder, adding all of them up.
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
double calculateCylinderCirculation(vector<vector<node>> &mesh, Parameters p)
{
    double vxn, vye, vxs, vyw, dPN, dPn, dNn, dPE, dPe, dEe, dPS, dPs, dSs, dPW, dPw, dWw, circ = 0;

    for (int i = 1; i < N - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {
            dPN = p.dy;
            dPn = p.dy / 2;
            dNn = p.dy / 2;
            vxn = (dPN / ((dPn * mesh[i][j].rho + dNn * mesh[i][j + 1].rho / p.rho_in))) * ((mesh[i][j + 1].stream - mesh[i][j].stream) / dPN);

            dPE = p.dx;
            dPe = p.dx / 2;
            dEe = p.dx / 2;
            vye = -(dPE / ((dPe * mesh[i][j].rho + dEe * mesh[i + 1][j].rho / p.rho_in))) * ((mesh[i + 1][j].stream - mesh[i][j].stream) / dPE);

            dPS = p.dy;
            dPs = p.dy / 2;
            dSs = p.dy / 2;
            vxs = -(dPS / ((dPs * mesh[i][j].rho + dSs * mesh[i][j - 1].rho / p.rho_in))) * ((mesh[i][j - 1].stream - mesh[i][j].stream) / dPS);

            dPW = p.dx;
            dPw = p.dx / 2;
            dWw = p.dx / 2;
            vyw = (dPW / ((dPw * mesh[i][j].rho + dWw * mesh[i - 1][j].rho / p.rho_in))) * ((mesh[i - 1][j].stream - mesh[i][j].stream) / dPW);

            if (mesh[i][j].is_solid == true)
            {
                if (mesh[i][j - 1].is_solid == false)
                {
                    circ += vxs * p.dx;
                }
                if (mesh[i][j + 1].is_solid == false)
                {
                    circ += -vxn * p.dx;
                }
                if (mesh[i - 1][j].is_solid == false)
                {
                    circ += -vyw * p.dy;
                }
                if (mesh[i + 1][j].is_solid == false)
                {
                    circ += vye * p.dy;
                }
            }
        }
    }
    return circ;
}

/**
 * Calculates the pressure coefficient of the mesh nodes. It uses the velocity values
 * to calculate the pressure coefficient of the nodes.
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
void calculateCp(vector<vector<node>> &mesh, Parameters p)
{
    double v = 0;
    for (int i = 1; i < N - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {
            v = sqrt(mesh[i][j].u * mesh[i][j].u + mesh[i][j].v * mesh[i][j].v);
            mesh[i][j].cp = 1 - (v / p.v_in) * (v / p.v_in);
        }
    }
}

/**
 * Calculates pressure and temperature of the mesh nodes.
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
void calculatePressureTemperature(vector<vector<node>> &mesh, Parameters p)
{
    double v = 0;
    for (int i = 1; i < N - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {
            v = sqrt((mesh[i][j].u * mesh[i][j].u) + (mesh[i][j].v * mesh[i][j].v));
            mesh[i][j].T = p.t_in + ((p.v_in * p.v_in) - (v * v)) / (2 * p.specific_heat);
            mesh[i][j].p = p.p_in * pow(mesh[i][j].T / p.t_in, p.gamma / (p.gamma - 1));
            mesh[i][j].cp = 1 - (v / p.v_in) * (v / p.v_in);
        }
    }
}

/**
 * Calculates the forces around the cylinder. It uses the pressure coefficient values
 * to calculate the lift and drag coefficients of the cylinder.
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
struct Coefficients cylinderForces(vector<vector<node>> &mesh, Parameters p)
{
    Coefficients c;

    double q = 0.5 * p.rho_in * p.p_in * p.v_in;
    for (int i = 1; i < N - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {
            if (mesh[i][j].is_solid == true)
            {
                if (mesh[i][j - 1].is_solid == false)
                {
                    c.C_L += (mesh[i][j - 1].cp * q + p.p_in) * p.dx;
                }

                if (mesh[i][j + 1].is_solid == false)
                {
                    c.C_L += -(mesh[i][j + 1].cp * q + p.p_in) * p.dx;
                }

                if (mesh[i - 1][j].is_solid == false)
                {
                    c.C_D += (mesh[i - 1][j].cp * q + p.p_in) * p.dy;
                }

                if (mesh[i + 1][j].is_solid == false)
                {
                    c.C_D += -(mesh[i + 1][j].cp * q + p.p_in) * p.dy;
                }
            }
        }
    }

    c.C_L = c.C_L / (q * 2 * p.cylinder_r);
    c.C_D = c.C_D / (q * 2 * p.cylinder_r);

    return c;
}

/**
 * Converts Cartesian coordinates to polar coordinates.
 *
 * @param x The x-coordinate in Cartesian coordinates.
 * @param y The y-coordinate in Cartesian coordinates.
 * @return A pair where the first element is the radius (r) and the second element is the angle (theta) in radians.
 */
pair<double, double> cartesianToPolar(double x, double y)
{
    double r = sqrt(x * x + y * y);
    double theta = atan2(y, x);
    return make_pair(r, theta);
}

/**
 * Computes the analytic stream function of the mesh nodes.
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
void computeAnalyticStream(vector<vector<node>> &mesh, Parameters p)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            pair<double, double> polar;
            pair<double, double> cartesian;

            polar = cartesianToPolar(mesh[i][j].x - p.H / 2, mesh[i][j].y - p.H / 2);

            double input_flow = p.v_in * p.H / 2;
            double cylinder = p.v_in * (polar.first - (p.cylinder_r * p.cylinder_r) / polar.first) * sin(polar.second);
            // double circulation = p.cylinder_r * p.cylinder_r * log(polar.first) * p.v_in; // <missing rotation speed

            mesh[i][j].stream = input_flow + cylinder; // + circulation;
        }
    }
}

/**
 * Compares two matrices and returns the maximum error between them.
 *
 * @param v1 Matrix (Vector of vectors) to compare
 * @param v2 Matrix (Vector of vectors) to compare
 */
vector<vector<node>> analyticError(vector<vector<node>> numerical, vector<vector<node>> analytical)
{
    vector<vector<node>> error_matrix(N, vector<node>(M));
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            error_matrix[i][j].stream = numerical[i][j].stream - analytical[i][j].stream;
        }
    }
    return error_matrix;
}