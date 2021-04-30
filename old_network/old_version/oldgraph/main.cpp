

#include "graph.h"
#include "paths.h"
#include "paths_bfs.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>

using namespace std;

int main()
{


    ifstream in("topology.txt");
    Graph g(in);
    g.showAdj();
    cout << endl;

    g.deleteVertex(2);
    g.showAdj();
    cout << endl;

    return 0;

    int v = 2;
    cout << "All nodes not connected to " << v << " are:\n";
    PathsBFS<Graph> paths(g,v);


    for( int vi=0; vi<g.V(); ++vi) {
        cout <<vi << " : " <<  paths.distance(vi) << endl;
        //if( !paths.hasPathTo(vi) ) cout << vi << endl;
    }
    return 0;
}
