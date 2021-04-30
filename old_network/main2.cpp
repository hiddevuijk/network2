
#include "graph.h"
#include "generate_graph.h"
#include "network.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

double min( double a, double b) {
	if( a < b ) return a;
	else return b;
}



int main()
{
	double pi = acos(-1);
    int Nx = 14;
    //int Ny = Nx;
    double Lx = Nx;
    double Ly = Lx*sqrt(3/4.);
    //double z = 3.2;
	double kappa = 1.;
	//long int seed  = 123456789;


    //Graph graph = generateGraph(Nx,Ny,Lx,z, seed);
	Graph graph(3);
	graph.setVertexPosition(0, 1,0);
	graph.setVertexPosition(1, 1,1);
	graph.setVertexPosition(2, 2,1);

	graph.addEdge(0,1);
	graph.addEdge(1,2);
	graph.addBend(1,0,2);

    ofstream top("topology.txt");
    graph.write(top);
    top.close();

	Network network(graph,Lx,Ly, kappa);

	cout << "initial phi :\n";
	cout << network.bends[0].phi0 << endl << endl;

	network.bends[0].phi0 = 0.25*pi;

	cout << "target phi :\n";
	cout << network.bends[0].phi0 << endl << endl;

	gsl_vector *df = gsl_vector_alloc(6);
	dEnergy( network.r, &network, df);

	cout << "Forces:\n";
	cout << gsl_vector_get(df, 0) << "\t";
	cout << gsl_vector_get(df, 1) << "\n";
	cout << gsl_vector_get(df, 2) << "\t";
	cout << gsl_vector_get(df, 3) << "\n";
	cout << gsl_vector_get(df, 4) << "\t";
	cout << gsl_vector_get(df, 5) << "\n";
	cout << endl;


    ofstream out0("r0.dat");
    network.savePositions(out0);
    out0.close();

	network.minimize();

	dEnergy( network.r, &network, df);

	cout << "phi:\t" << network.get_phi(0) << endl << endl;

	cout << "Forces:\n";
	cout << gsl_vector_get(df, 0) << "\t";
	cout << gsl_vector_get(df, 1) << "\n";
	cout << gsl_vector_get(df, 2) << "\t";
	cout << gsl_vector_get(df, 3) << "\n";
	cout << gsl_vector_get(df, 4) << "\t";
	cout << gsl_vector_get(df, 5) << "\n";
	cout << endl;


	cout << "energy: \n";
	cout << network.totalEnergy() << endl;
	cout << network.edgeEnergy() << endl;
	cout << network.bendEnergy() << endl;

    ofstream out("r.dat");
    network.savePositions(out);
    out.close();



    return 0;

}

