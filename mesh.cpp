// Last update: 2024/10/06
// Author: Ricard Arbat Carandell

// Master in Aerospace Engineering - Computational Engineering
// Universitat Politècnica de Catalunya (UPC) - BarcelonaTech
// Overview: Mesh definition functions

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "parameters.cpp"
using namespace std;

// node struct
struct node
{
    double x, y, u, v, cp, stream, rho, T, p;
    bool is_solid;
};

/**
 * Fills the mesh with nodes as it defines their positions.
 * Also defines the cylinder solid nodes by considering if they are inside a circle
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
void buildMesh(vector<vector<node>> &mesh, Parameters p)
{

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            mesh[i][j].x = (i * p.dx) + 0.5 * p.dx;
            mesh[i][j].y = (j * p.dy) + 0.5 * p.dy;

            double dist = sqrt(pow(mesh[i][j].x - p.cylinder_x, 2) + pow(mesh[i][j].y - p.cylinder_y, 2));
            if (dist < p.cylinder_r)
            {
                mesh[i][j].is_solid = true;
            }
            else
            {
                mesh[i][j].is_solid = false;
            }
        }
    }
}

/**
 * Sets the initial density of the mesh nodes.
 * Also checks if the node is solid and sets the solid density
 *
 * @param mesh Mesh matrix
 * @param p Parameters of the simulation
 */
void setRho(vector<vector<node>> &mesh, Parameters p)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (mesh[i][j].is_solid == true)
            {
                mesh[i][j].rho = p.solid_density;
            }
            else
            {
                mesh[i][j].rho = p.initial_density;
            }
        }
    }
}

/**
 * Sets the initial stream value of the mesh nodes.
 * Also checks if the node is solid and sets the solid density.
 *
 * @param mesh Mesh matrix
 * @param p Parameters of the simulation
 */
void setStream(vector<vector<node>> &mesh, Parameters p)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (mesh[i][j].is_solid)
            {
                mesh[i][j].stream = p.solid_stream;
            }
            else
            {
                mesh[i][j].stream = p.start_stream;
            }
        }
    }
}

/**
 * Exports the mesh data to a CSV file
 *
 * @param mesh Mesh matrix
 * @param filename Name of the file to be exported
 */
void exportData(vector<vector<node>> &mesh, string filename = "output/output.csv")
{
    ofstream outfile(filename);
    outfile << "X,Y,U,V,Stream,Density,Cp,Solid,P,T" << endl;
    for (int i = 1; i < N - 1; i += 1)
    {
        for (int j = 1; j < M - 1; j += 1)
        {
            outfile << mesh[i][j].x << "," << mesh[i][j].y << "," << mesh[i][j].u << "," << mesh[i][j].v << "," << mesh[i][j].stream << "," << mesh[i][j].rho << "," << mesh[i][j].cp << "," << mesh[i][j].is_solid << "," << mesh[i][j].p << "," << mesh[i][j].T << endl;
        }
    }
    outfile.close();
}
