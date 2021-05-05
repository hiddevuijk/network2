
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
    int n3 = g.addNode();
    int n4 = g.addNode();
    int n5 = g.addNode();
    int n6 = g.addNode();
    int n7 = g.addNode();

    g.addBond(n0,n1);

    g.addBond(n1,n2);
    g.addBond(n2,n3);
    g.addBond(n3,n0);

    g.addBend(n0,n3,n1);
    g.addBend(n1,n0,n2);
    g.addBend(n2,n1,n3);
    g.addBend(n3,n2,n0);

    g.polymerize(n0,n1);
    g.polymerize(n1,n2);
    g.polymerize(n2,n3);
    g.polymerize(n3,n0);


    g.addBond(n4,n5);
    g.addBond(n5,n6);
    g.addBond(n6,n7);

    g.addBend(n5,n4,n6);
    g.addBend(n6,n5,n7);
    g.polymerize(n5,n6);


    //g.showBonds();
    //g.showBends(); 
   
    cout << "-----------------------------\n"; 

    g.delBond(n4,n5);
    
    g.showBonds();
    g.showBends(); 

    return 0;
}

