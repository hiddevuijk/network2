
#include "graph.h"
#include "generate_graph.h"
#include "network.h"
#include "ConfigFile.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


int main()
{
  double Lx     = 10.;
  double Ly     = Lx;
  double kappa  = 1.;

  Graph graph;
  graph.addNode(0., 1.);
  graph.addNode(1., 1.);
  graph.addNode(2., 1.);

  graph.addBond(0, 1, 0, 0);
  graph.addBond(1, 2, 0, 0);

  graph.addBend(1,0,2);

  Network network(graph, Lx, Ly, kappa);

  cout << network.getBendAngle(0)/acos(-1) << endl;

  network.moveNode(0, 1., 2.);

  cout << network.getBendAngle(0)/acos(-1) << endl;
  cout << network.getBendDphi(0)/acos(-1) << endl;
  cout << network.getBendDistanceji(0) << endl;
  cout << network.getBendDistancejk(0) << endl;

  return 0;

}


