// Last update: 2024/10/19
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
 * @param type Type of problem to solve
 */
void build_mesh(vector<vector<vector<node>>> &mesh, string type)
{
    for (int t = 0; t < time_steps; t++)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                mesh[t][i][j].x = -1 + (i * dx) + (0.5 * dx);
                mesh[t][i][j].y = (j * dy) + (0.5 * dy);
            }
        }
    }
}

/**
 * Sets the initial density of the mesh nodes.
 * also checks if the node is solid and sets the solid density
 *
 * @param mesh_t Mesh matrix at time t
 * @param variable Value to set
 * @param name Name of the variable to set
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

/**
 * Exports the mesh data to a CSV file
 *
 * @param mesh Mesh matrix at time t
 * @param filename Name of the file to be exported
 * @param folder Folder where the file will be saved
 * @param how_many How many time steps to export. Either "all" or "last".
 **/
void export_data(vector<vector<vector<node>>> &mesh, string filename = "output.csv", string folder = "smith_hutton", string how_many = "all")
{

    string header = "s,X,Y,phi,U,V";
    string rows = header + "\n";

    if (how_many == "all")
    {
        for (int t = 0; t < time_steps; t++)
        {
            cout << "Exporting data at time " << t << endl;
            for (int i = 1; i < N - 1; i += 1)
            {
                for (int j = 1; j < M - 1; j += 1)
                {
                    string row = to_string(t * delta_t) + "," + to_string(mesh[t][i][j].x) + "," + to_string(mesh[t][i][j].y) + "," + to_string(mesh[t][i][j].phi) + "," + to_string(mesh[t][i][j].u) + "," + to_string(mesh[t][i][j].v);
                    rows.append(row + "\n");
                }
            }
        }
    }
    else if (how_many == "last")
    {
        int t = time_steps - 1;
        cout << "Exporting data last timestep " << endl;
        for (int i = 1; i < N - 1; i += 1)
        {
            for (int j = 1; j < M - 1; j += 1)
            {
                string row = to_string(t * delta_t) + "," + to_string(mesh[t][i][j].x) + "," + to_string(mesh[t][i][j].y) + "," + to_string(mesh[t][i][j].phi) + "," + to_string(mesh[t][i][j].u) + "," + to_string(mesh[t][i][j].v);
                rows.append(row + "\n");
            }
        }
    }

    string path = folder + "/" + filename;

    ofstream outfile(path);
    outfile << rows << endl;
    outfile.close();
}

/**
 * Exports the mesh data to a CSV file
 *
 * @param mesh Mesh matrix at time t
 * @param scheme Scheme used
 * @param outfile File to write the data
 **/
void export_data_at_outlet(vector<vector<vector<node>>> &mesh, string scheme, ofstream &outfile)
{

    string rows = scheme + ",";

    int j = 2;
    for (int i = N / 2; i < N - 1; i++)
    {
        string row = to_string(mesh[time_steps - 1][i][j].phi) + ",";
        rows.append(row);
        cout << "x" << mesh[time_steps - 1][i][j].x << " y" << mesh[time_steps - 1][i][j].y << " phi" << mesh[time_steps - 1][i][j].phi << endl;
    }
    rows.append("\n");

    outfile << rows << endl;
}

/**
 * Makes a filename with the time and name
 *
 * @param Pe Peclet number
 * @param scheme Scheme used
 * @param name Name of the file
 * @param type Type of problem
 * @param delta_t Time step size
 *
 */
string file_name(double Pe, string scheme, string name, string type, double delta_t = 0.001)
{
    std::ostringstream oss;
    if (Pe > 100000.0)
    {
        string Pe_s = "1000000";
        oss << name << "_Pe_" << Pe_s << "_S_" << scheme << "_M_" << M << "_type_" << type << "_t_" << delta_t << ".csv";
    }
    else
    {
        oss << name << "_Pe_" << Pe << "_S_" << scheme << "_M_" << M << "_type_" << type << "_t_" << delta_t << ".csv";
    }

    std::string var = oss.str();

    return var;
}

/**
 * Prints the value of phi in the mesh graphically, like a matrix
 *
 * @param values Matrix (Vector of vectors) containing the phi values
 */

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