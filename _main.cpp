
#include "graph.h"
#include "generate_graph.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


int main()
{

    int Nx = 20; 
    int Ny = 20;
    double Lx = 10;
    double z = 30.;
    long int seed = 1234567890;
    string topologyName = "topology.dat";
    string r0Name = "r0.dat";

    Graph graph = generateGraph(Nx,Ny,Lx,z,seed,0);


    ofstream top(topologyName);
    graph.write(top);
    top.close();

    ofstream out0(r0Name);
    graph.writePositions(out0);
    out0.close();



    return 0;
}

