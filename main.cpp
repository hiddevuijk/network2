
#include "graph.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


int main()
{

    Graph g;

    int n0 = g.addNode();
    int n1 = g.addNode();
    int n2 = g.addNode();

    g.addBond(n0,n1);
    g.addBond(n0,n2);

    cout <<g.Nnodes() << endl;
    g.delNode(0);
    cout <<g.Nnodes() << endl;

    g.delNode(0);
    cout <<g.Nnodes() << endl;

    g.delNode(0);
    cout <<g.Nnodes() << endl;

    return 0;
}

