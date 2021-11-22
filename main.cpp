
#include "graph.h"
#include "generate_graph.h"
#include "network.h"
#include "ConfigFile.h"
#include "minimize.h"
#include "minimize_gsl.h"

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
  double alpha  = config.read<double>("alpha");


  double error  = config.read<double>("e");
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


  Minimizer minimizer(network, error, emax, dt0, dtmax, dtmin, finc, fdec,
                      Nmin, alpha0, falpha, m);

  double eLine = 1e-9;
  double dLine = 1e-9;
  double e     = 1e-12;
  int maxIter  = 100;
  MinimizerGSL minimizerGSL(&network, eLine, dLine, e, maxIter);

  ofstream top(topologyName + ".dat");
  graph.write(top);
  top.close();

  ofstream out0(r0Name + ".dat");
  network.savePositions(out0);
  out0.close();

  ofstream gEout(gammaEName + ".dat");

  double Hs, Hb; // stretch and bend energy

  vector<double> sigma;
  double V = Lx*Ly;


  int i = 0;
  while (fabs(network.getGamma()) < gmax) {
    network.shearAffine(dg);

    //minimizerGSL.minimize(network);
    minimizer.minimize(network);

    Hs = network.getBondEnergy();
    Hb = network.getBendEnergy();
    sigma = network.getStress();

    if (i % 10 == 0) {
      cout << network.getGamma() << "\t" << (Hs + Hb)/V << endl;
    }
    ++i;
      
    gEout << network.getGamma() << "\t"
          << Hs/V << "\t"
          << Hb/V << "\t"
          << sigma[0]/V << "\t"
          << sigma[1]/V << "\t"
          << sigma[2]/V << "\t"
          << sigma[3]/V << endl;

    ofstream out(rName + "_" + to_string(i) + ".dat");
    network.savePositions(out);
    out.close();

    dg *= alpha;
  } 


  return 0;

}


