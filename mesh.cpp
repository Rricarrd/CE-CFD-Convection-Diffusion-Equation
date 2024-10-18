// Last update: 2024/10/17
// Author: Ricard Arbat Carandell

// Master in Aerospace Engineering - Computational Engineering
// Universitat Politècnica de Catalunya (UPC) - BarcelonaTech
// Overview: Mesh definition functions

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include "parameters.cpp"
using namespace std;

// node struct
struct node
{
    double x, y, u, v, phi;
};

/**
 * Fills the mesh with nodes as it defines their positions.
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 */
void build_mesh(vector<vector<vector<node>>> &mesh)
{
    for (int t = 0; t < time_steps; t++)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                mesh[t][i][j].x = (i * dx) + (0.5 * dx);
                mesh[t][i][j].y = (j * dy) + (0.5 * dy);
            }
        }
    }
}

/**
 * sets the initial density of the mesh nodes.
 * also checks if the node is solid and sets the solid density
 *
 * @param mesh mesh matrix
 */
void set_mesh_value(vector<vector<node>> &mesh_t, float variable, string name)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (name == "u")
            {
                mesh_t[i][j].u = variable;
            }
            else if (name == "v")
            {
                mesh_t[i][j].v = variable;
            }
            else if (name == "phi")
            {
                mesh_t[i][j].phi = variable;
            }
        }
    }
}

void set_smith_hutton_problem(vector<vector<vector<node>>> &mesh)
{
    for (int t = 0; t < time_steps; t++)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                mesh[t][i][j].u = 2 * mesh[t][i][j].y * (1 - ((mesh[t][i][j].x - 1) * (mesh[t][i][j].x - 1)));
                mesh[t][i][j].v = -2 * (mesh[t][i][j].x - 1) * (1 - (mesh[t][i][j].y * mesh[t][i][j].y));
            }
        }
    }
}

/**
 * Exports the mesh data to a CSV file
 *
 * @param mesh Mesh matrix at time t
 * @param filename Name of the file to be exported
 **/
void export_data(vector<vector<vector<node>>> &mesh, string filename = "output/output.csv")
{

    ofstream outfile(filename);
    outfile << "s,X,Y,U,V,phi" << endl;
    for (int t = 0; t < time_steps; t++)
    {
        for (int i = 1; i < N - 1; i += 1)
        {
            for (int j = 1; j < M - 1; j += 1)
            {
                outfile << t * delta_t << "," << mesh[t][i][j].x << "," << mesh[t][i][j].y << "," << mesh[t][i][j].u << "," << mesh[t][i][j].v << "," << mesh[t][i][j].phi << endl;
            }
        }
    }
    outfile.close();
}

/**
 * Makes a filename with the time and name
 *
 * @param time Time of the simulation
 * @param name Name of the file
 *
 */
string file_name(int time, string name = "output")
{
    std::ostringstream oss;
    oss << name << "_t" << time << ".csv";
    std::string var = oss.str();

    return var;
}

/**
 * Prints the value of phi in the mesh graphically, like a matrix
 *
 * @param mesh Mesh matrix at time t
 **/
void print_phi_matrix(vector<vector<double>> &values)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            cout << values[i][j] << " ";
        }
        cout << endl;
    }
}