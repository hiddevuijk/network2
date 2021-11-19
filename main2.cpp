
#include "graph.h"
#include "generate_graph.h"
#include "network.h"
#include "ConfigFile.h"
#include "minimize.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


int main()
{
  double Lx = 10.;
  double Ly = Lx;
  double kappa = 1.;

  
  Graph graph; 
  graph.addNode(1., 1.);
  graph.addNode(1., 9.);
  graph.addNode(3., 9.);

  graph.addBond(0, 1, 0, -1);
  graph.addBond(1, 2, 0, 0);

  graph.addBend(1, 0, 2);

  Network network(graph, Lx, Ly, kappa);

  network.moveNode(2, 1., 8.);
  double d;
  d = network.getBendDistanceji(0);
  cout << d << endl;
  d = network.getBendDistancejk(0);
  cout << d << endl;

  cout << endl;




  double phi;
  phi = network.bends_[0].phi0;
  cout << "phi0:\t" <<  phi/M_PI << endl;

  phi = network.getBendAngleCW(0);
  cout << "phiCW:\t" << phi/M_PI << endl;

  phi = network.getBendAngle(0);
  cout << "phi:\t" << phi/M_PI << endl;

  phi = network.getBendDphi(0);
  cout << "Dphi:\t" << phi/M_PI << endl;


  vector<double> F(6, 0.);
  network.bendForce(F);

  cout << "| ";
  for(unsigned int i=0; i < F.size(); ++i) {
    cout << F[i] << " | ";
  }
  cout << endl;
  return 0;

}


