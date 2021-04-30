
#include "graph.h"
#include "generate_graph.h"
#include "network.h"

#include "ConfigFile.h"

#include "boost/random.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main()
{

	ConfigFile config("input.txt");

	int Nx = config.read<int>("Nx");
    int Ny = Nx;
	// assert( Ny % 2 == 0 );
    double Lx = config.read<double>("Lx");
    double Ly = Lx*sqrt(3/4.);
    double z = config.read<double>("z");
	double kappa = config.read<double>("kappa");
	long int seed  = config.read<long int>("seed");
    double s = config.read<double>("s");
	
    double gmax = config.read<double>("gmax");
    double dg = config.read<double>("dg");
    double gamma = 0;
	double alpha = config.read<double>("alpha");


	double e = config.read<double>("e");
	double emax = config.read<double>("emax");
    int Nmin = config.read<int>("Nmin");
	double dt0 = config.read<double>("dt0");
	double dtmax = config.read<double>("dtmax");
	double dtmin = config.read<double>("dtmin");
	double finc = config.read<double>("finc");
	double fdec = config.read<double>("fdec");
    double alpha0 = config.read<double>("alpha0");
    double falpha = config.read<double>("falpha");
    double m = config.read<double>("m");

	string topologyName = config.read<string>("topologyName");
	string r0Name = config.read<string>("r0Name");
	string rName = config.read<string>("rName");
	string gammaEName = config.read<string>("gammaEName");


    Graph graph = generateGraph(Nx,Ny,Lx,z, seed, s);
	Network network(graph,Lx,Ly, kappa);


    ofstream top(topologyName);
    graph.write(top);
    top.close();


    ofstream out0(r0Name);
    network.savePositions(out0);
    out0.close();

	ofstream gEout( gammaEName );

	double Hs, Hb;
	vector<double> sigma(4);
	int i=0;

    //network.stretchXAffine(-0.05);
    //network.stretchYAffine(-0.05);
    //network.minimize(e, dt0, dtmax, dtmin,finc, fdec, Nmin, alpha0, falpha,  m);

    ofstream outc("rcompressed.dat");
    network.savePositions(outc);
    outc.close();

    int NF = 0;
    while( fabs(gamma) < gmax ) {
		
        gamma += dg;

		//network.shearAffine(dg);
		//network.stretchX(dg);
		//network.stretchY(dg);
		network.shearAffine(dg);
		network.minimize(e, emax, dt0, dtmax, dtmin,finc, fdec, Nmin, alpha0, falpha,  m);

		Hs = network.edgeEnergy();	
		Hb = network.bendEnergy();	
        sigma = network.stress();
		//cout << network.get_forceNorm() << endl;
		if(i%10 == 0 ){
			cout << gmax << "\t" << gamma << "\t" << Hs + Hb << "\t" << network.minimizer.Fnorm/ network.minimizer.N  <<
                    "\t" << network.minimizer.Fmax() <<  "\t" << network.minimizer.NF << endl;
            cout << "______________________________________________________________" << endl;
		}
		i++;

        NF += network.minimizer.NF;	

		gEout << gamma << '\t';
		gEout << Hs << '\t';
		gEout << Hb << '\t';
		gEout << sigma[0] << '\t';
		gEout << sigma[1] << '\t';
		gEout << sigma[2] << '\t';
		gEout << sigma[3] << endl;

		dg *= alpha;
    }
    cout << NF << endl;

    ofstream out(rName);
    network.savePositions(out);
    out.close();



    return 0;

}

