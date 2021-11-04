
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


  // Read data from the input file
  ConfigFile config("input.txt");

  int Nx = config.read<int>("Nx");
  int Ny = Nx;
  assert(Ny % 2 == 0); // due to periodic boundary conditions

  double Lx     = config.read<double>("Lx");
  double Ly     = Lx*sqrt(3/4.);
  double z      = config.read<double>("z");
  double kappa  = config.read<double>("kappa");
  long int seed = config.read<long int>("seed");
  double s      = config.read<double>("s");

  double gmax   = config.read<double>("gmax");
  double dg     = config.read<double>("dg");
  double gamma  = 0;
  double alpha  = config.read<double>("alpha");


  double e      = config.read<double>("e");
  double emax   = config.read<double>("emax");
  int    Nmin   = config.read<int>("Nmin");
  double dt0    = config.read<double>("dt0");
  double dtmax  = config.read<double>("dtmax");
  double dtmin  = config.read<double>("dtmin");
  double finc   = config.read<double>("finc");
  double fdec   = config.read<double>("fdec");
  double alpha0 = config.read<double>("alpha0");
  double falpha = config.read<double>("falpha");
  double m	    = config.read<double>("m");

  string topologyName = config.read<string>("topologyName");
  string r0Name       = config.read<string>("r0Name");
  string rName        = config.read<string>("rName");
  string gammaEName   = config.read<string>("gammaEName");


  Graph graph = generateGraph(Nx, Ny, Lx, z, seed, s);
  Network network(graph, Lx, Ly, kappa);


  // save the topology of the network
  ofstream top(topologyName + ".dat");
  graph.write(top);
  top.close();

  // save the initial positions on the nodes 
  ofstream out0(r0Name + ".dat");
  network.savePositions(out0);
  out0.close();

  // out stream for the data
  ofstream gEout(gammaEName + ".dat");

  double Hs, Hb; // stretching and bending energy


  vector<double> sigma(4);  // stress tensor: s_xx, s_xy, s_yx, s_yy 


  int i=0;

  // deformation before stretching
  //network.stretchXAffine(-0.05);
  //network.stretchYAffine(-0.05);
  //network.minimize(e, dt0, dtmax, dtmin,finc, fdec, Nmin, alpha0, falpha,  m);

  // save  data
  //ofstream outc("rcompressed.dat");
  //network.savePositions(outc);
  //outc.close();

  int NF = 0; // total number of force calculations

  while (fabs(gamma) < gmax) {

    gamma += dg;

    // deform and minimize the energy

    //network.shearAffine(dg);
    //network.stretchX(dg);
    //network.stretchY(dg);
    network.shearAffine(dg);
    network.minimize(e, emax, dt0, dtmax, dtmin,finc, fdec,
                     Nmin, alpha0, falpha,  m);

    Hs = network.bondEnergy();	
    Hb = network.bendEnergy();	
    sigma = network.stress();
    //cout << network.get_forceNorm() << endl;

    // std::out information 
    if (i%10 == 0) {
      cout << gmax << "\t" 
           << gamma << "\t" 
           << Hs + Hb << "\t"
           << network.minimizer.Fnorm/ network.minimizer.N << "\t"
           << network.minimizer.Fmax() <<  "\t"
           << network.minimizer.NF << endl;
      cout << "______________________________________________________" << endl;
    }

    i++;

    NF += network.minimizer.NF;	

    // save energy and stress
    gEout << gamma << '\t'
          << Hs << '\t'
          << Hb << '\t'
          << -1*sigma[0] << '\t'
          << -1*sigma[1] << '\t'
          << -1*sigma[2] << '\t'
          << -1*sigma[3] << endl;

    ofstream out(rName + "_" + to_string(i) + ".dat");
    network.savePositions(out);
    out.close();

   


    // increase gamma increment s.t. gamma is logarithmic
    dg *= alpha;
  }

  cout << NF << endl;

  // save the final positions of the network
  ofstream out(rName + ".dat");
  network.savePositions(out);
  out.close();



  return 0;

}


