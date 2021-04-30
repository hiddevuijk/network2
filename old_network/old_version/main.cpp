

#include "graph.h"
#include "generate_network.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>

using namespace std;

int main()
{
    //Graph g(6);
    //g.addEdge(0,1);
    //g.addEdge(1,2);
    //g.addEdge(2,0);
    //g.addEdge(3,4);
    //g.addEdge(4,5);
    //g.addEdge(5,3);

    //g.addEdge(0,3);
    //g.addEdge(1,3);
    //g.addEdge(1,4);
    //g.addEdge(2,4);
    //g.addEdge(2,5);
    //
    //g.addBend(0,2,1);
    //g.addBend(1,0,2);
    //g.addBend(2,1,0);

    //g.polymerize(0,1);
    //g.polymerize(1,2);
    //g.polymerize(2,0);

    //cout << "done " << endl;

    int Nx = 4;
    int Ny = 4;
    double Lx = 4;
    Graph g = generateNetwork(Nx, Ny, Lx);

    Graph::Bend *b = &g.vertices[0].bends[0];

    cout << b->I << endl;
    cout << b->next << endl;
    cout << b->nextBend << endl;
    cout << b->nextBend->I << endl;
    cout << &g.vertices[1].bends[0] << endl;
       

    return 0;
}
