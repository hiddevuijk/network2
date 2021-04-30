#ifndef GUARD_NETWORK_H
#define GUARD_NETWORK_H

/*
to do:
    -- write EdEnergy

	-- change bend::energy_dEnergy
*/

#include "vec2.h"
#include "graph.h"
#include "generate_graph.h"

#include "gsl/gsl_multimin.h"

#include <math.h>


double energy( const gsl_vector *v, void *params);
void dEnergy( const gsl_vector *v, void *params, gsl_vector *df);
void EdE( const gsl_vector *x, void *params, double *f, gsl_vector *df);

class Network
{
  private:
    class Edge;
    class Bend;

  public:
    Network(const Graph&, double Lx, double Ly, double kappa);
	~Network();

	void shake( boost::mt19937 &rng, double sigma);

    void minimize(double eLine, double dLine, double e);
    void minimize2(double eLine, double dLine, double e);

	// deform Network
    void shear(double delta_gamma); 
    void shearAffine(double delta_gamma); 
    void stretchX(double delta_epsilonX); 
    void stretchXAffine(double delta_epsilonX); 
    void stretchY(double delta_epsilonY); 
    void stretchYAffine(double delta_epsilonY); 

    int get_Nv() const {return Nv; }
    int get_Nedges() const { return Ne; }
    int get_Nbends() const { return Nb; }

    double edgeEnergy() const; 
    double bendEnergy() const; 
    double totalEnergy() const; 

	std::vector<double> stress() const;
	std::vector<double> stress2() const;

    double get_edgeEnergy( int ei, const gsl_vector *r) const
        { return edges[ei].energy(r, this); }

    void get_edgeDEnergy( int ei, const gsl_vector *r, gsl_vector *df) const
        {  edges[ei].dEnergy(r, df, this); }

    double get_edgeEdE( int ei, const gsl_vector *r, gsl_vector *df) const
        { return edges[ei].energy_dEnergy(r, df, this); }

	double get_bendEnergy( int bi, const gsl_vector *r) const
		{ return bends[bi].energy(r, this); }

	void get_bendDEnergy( int bi, const gsl_vector *r, gsl_vector *df) const
		{ bends[bi].dEnergy(r, df, this); }

	double get_bendEdE( int bi, const gsl_vector *r, gsl_vector *df) const 
		{ return bends[bi].energy_dEnergy(r ,df, this); }

	double get_forceNorm() const;

    std::vector<double> get_positions() const;
    void savePositions(std::ostream& out) const;

    void set_Lx(double Lxx);
    void set_Ly(double Lyy);

  private:
    class Edge{
      public:
        int i,j; 
        int xb, yb;
        double l0; 

		double get_l( const gsl_vector *r, const Network *net) const;
		void set_l0( const gsl_vector *r, const Network *net) { l0 = get_l(r,net); }
		void set_l0( double l00) { l0 = l00; }

        double energy(const gsl_vector *r, const Network *net) const;
        void dEnergy( const gsl_vector *r,  gsl_vector *df, const Network *net) const;
        double energy_dEnergy(const gsl_vector *r, gsl_vector *df, const Network *net) const;
    };

    class Bend{
	  public:
		int i,j,k;
		int xib, yib;
		int xkb, ykb;

		double kappa;
		double phi0;

		double get_lji( const gsl_vector *r, const Network *net) const;
		double get_ljk( const gsl_vector *r, const Network *net) const;
		double get_phi( const gsl_vector *r, const Network *net) const;

		void set_phi0( const gsl_vector *r, const Network *net);
	

		double energy( const gsl_vector *r, const Network *net) const;
		void dEnergy ( const gsl_vector *r, gsl_vector *df, const Network *net) const;
		double energy_dEnergy( const gsl_vector *r, gsl_vector *df, const Network *net) const;
    };


  public: // remove
	void set_kappa();
	void set_phi0();
	void set_l0();
	void set_l0(double l00);

	double get_phi(int bi) const { return bends[bi].get_phi(r,this); };
	
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;

    const gsl_multimin_fdfminimizer_type *T2;
    gsl_multimin_fdfminimizer *s2;

    gsl_vector *r;
    gsl_multimin_function_fdf functions;
    std::vector<Edge> edges;
    std::vector<Bend> bends;
  public:
    double kappa;
    double Lx, Ly;
    int Nv,Ne,Nb;

    double gamma;  // shear deformation
	double epsilonX; // bulk deformation in x
	double epsilonY; // bulk deformation in y

    int maxIter = 1000000;
};

///////////////////////////////
// Network Member Functions  //
///////////////////////////////
Network::~Network()
{
	gsl_vector_free(r);
	gsl_multimin_fdfminimizer_free(s);
	gsl_multimin_fdfminimizer_free(s2);
}

Network::Network(const Graph& g, double Lxx, double Lyy, double kappaa)
{
    
    Nv = g.Nvertices();
	Lx = Lxx;
	Ly = Lyy;

    //T = gsl_multimin_fdfminimizer_conjugate_fr;
    //T = gsl_multimin_fdfminimizer_conjugate_pr;
    //T = gsl_multimin_fdfminimizer_vector_bfgs;
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    //T = gsl_multimin_fdfminimizer_steepest_descent;

    s = gsl_multimin_fdfminimizer_alloc(T,2*Nv);

    //T2 = gsl_multimin_fdfminimizer_conjugate_fr;
    //T2 = gsl_multimin_fdfminimizer_conjugate_pr;
    //T2 = gsl_multimin_fdfminimizer_vector_bfgs;
    T2 = gsl_multimin_fdfminimizer_vector_bfgs2;
    //T2 = gsl_multimin_fdfminimizer_steepest_descent;
    s2 = gsl_multimin_fdfminimizer_alloc(T2,2*Nv);

    r = gsl_vector_alloc(2*Nv);

    functions.n = 2*Nv;
    functions.f = &energy;
    functions.df = &dEnergy;
    functions.fdf = &EdE;
    functions.params = this;

    for(int vi=0; vi<Nv; ++vi ) {
        vec2 p = g.getVertexPosition(vi);
        gsl_vector_set( r, 2*vi  , p.x);   
        gsl_vector_set( r, 2*vi+1, p.y);   
    }

    std::vector<std::vector<int> > e = g.getEdges();
    Ne = e.size();
    edges = std::vector<Edge>(Ne);
    for(int ei=0; ei< Ne; ++ei ) {
        edges[ei].i  = e[ei][0]; 
        edges[ei].j  = e[ei][1]; 
        edges[ei].xb = e[ei][2]; 
        edges[ei].yb = e[ei][3]; 
        edges[ei].l0 = e[ei][4]; 
    }

	std::vector<std::vector<int> > b = g.getBends();
	Nb = b.size();
	bends = std::vector<Bend>(Nb);
	for(int bi=0; bi < Nb; ++bi ){
		bends[bi].j   = b[bi][0];

		bends[bi].i   = b[bi][1];
		bends[bi].xib = b[bi][2];
		bends[bi].yib = b[bi][3];

		bends[bi].k   = b[bi][4];
		bends[bi].xkb = b[bi][5];
		bends[bi].ykb = b[bi][6];
	}

	kappa = kappaa;
	set_kappa();

	gamma = 0;
	epsilonX = 0;
	epsilonY = 0;

	set_phi0();
	set_l0();
	
}

double Network::get_forceNorm() const
{    
	gsl_vector *df;
	df = gsl_vector_alloc(2*Nv);
    gsl_vector_scale( df, 0); // set all elements to 0

    for( int ei=0; ei< Ne; ++ei ) {
        get_edgeDEnergy(ei, r, df);
    } 


	for( int bi=0; bi < Nb; ++bi) {
		get_bendDEnergy(bi, r, df);
	}

	double norm = 0;
	double f;
	for(int vi=0; vi<2*Nv; ++ vi) {
		f = gsl_vector_get(df,vi);
		norm += f*f;
	}
		
	gsl_vector_free(df);

	return norm;
}
void Network::shake( boost::mt19937 &rng, double sigma)
{
	boost::normal_distribution<> nd(0.0, sigma);
	boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > random_normal(rng, nd);

	double xy;
	for( int vi=0; vi<2*Nv; ++vi){
		xy = gsl_vector_get(r, vi);
		gsl_vector_set(r, vi, xy + random_normal() );
	}
	
}


void Network::set_kappa()
{
	double lji, ljk;
	for( int bi=0; bi< Nb; ++bi ) {
		lji = bends[bi].get_lji(r,this);
		ljk = bends[bi].get_ljk(r,this);
		//bends[bi].kappa = kappa*( bends[bi].get_lji(r, this) + bends[bi].get_ljk(r, this) ) / 2;
		bends[bi].kappa = kappa*( lji + ljk ) / 2;
	}
}

void Network::set_phi0()
{
	for(unsigned int bi = 0; bi < bends.size(); ++bi) {
		bends[bi].set_phi0(r,this);
	}
}

void Network::set_l0()
{
	for(unsigned int ei=0; ei<edges.size(); ++ei) {
		edges[ei].set_l0(r,this);
	}
}

void Network::set_l0(double l00)
{
	for(unsigned int ei=0; ei<edges.size(); ++ei) {
		edges[ei].set_l0(l00);
	}
}

void Network::shearAffine( double delta_gamma)
{
    gamma += delta_gamma;
    double x,y;
    //affine deformation
    for( int vi=0; vi<Nv; ++vi ) {
        x = gsl_vector_get(r,2*vi);
        y = gsl_vector_get(r,2*vi+1);
        gsl_vector_set( r, 2*vi, x+delta_gamma*y);
    }

}


void Network::shear( double delta_gamma)
{ gamma += delta_gamma; }

void Network::stretchX( double delta_epsilonX )
{ epsilonX += delta_epsilonX; }

void Network::stretchXAffine( double delta_epsilonX)
{
	// add Affine deformation
	double x;
	for( int vi = 0; vi < Nv ; ++vi) { 
		x = gsl_vector_get(r, 2*vi);
		gsl_vector_set(r, 2*vi, x*(1+delta_epsilonX) );
	}
	epsilonX += delta_epsilonX;
}

void Network::stretchY( double delta_epsilonY)
{ epsilonY += delta_epsilonY; }


void Network::stretchYAffine( double delta_epsilonY)
{

	// add Affine deformation
	double y;
	for(int vi=0; vi < Nv; ++vi ) {
		y = gsl_vector_get(r, 2*vi + 1);
		gsl_vector_set(r, 2*vi + 1, y*(1+delta_epsilonY) );
	}

	epsilonY += delta_epsilonY;
}
void Network::minimize( double eLine, double dLine, double e)
{
    // add params to network
    gsl_multimin_fdfminimizer_set(s, &functions, r, eLine, dLine);

    int iter = 0;
    int status;
    do{
        iter++;    
		status = gsl_multimin_fdfminimizer_iterate(s);
        if( status ) break;
        status = gsl_multimin_test_gradient( s->gradient, e);
        
    } while( status == GSL_CONTINUE && iter < maxIter);
    if( iter >= maxIter ) std::cout << "\t Fuck \n";

	//std::cout << "\t" << iter << std::endl;


    // copy to r
    gsl_vector * x = gsl_multimin_fdfminimizer_x(s);
    for(int i=0;i<2*Nv; ++i) {
        gsl_vector_set(r,i, gsl_vector_get(x,i) );
    }
}

void Network::minimize2( double eLine, double dLine, double e)
{
    // add params to network
    gsl_multimin_fdfminimizer_set(s2, &functions, r, eLine, dLine);

    int iter = 0;
    int status;
    do{
        iter++;    
		status = gsl_multimin_fdfminimizer_iterate(s2);
        if( status ) break;
        status = gsl_multimin_test_gradient( s2->gradient, e);
        
    } while( status == GSL_CONTINUE && iter < maxIter);
    if( iter >= maxIter ) std::cout << "\t Fuck \n";

	//std::cout << "\t" << iter << std::endl;

    // copy to r
    gsl_vector * x = gsl_multimin_fdfminimizer_x(s2);
    for(int i=0;i<2*Nv; ++i) {
        gsl_vector_set(r,i, gsl_vector_get(x,i) );
    }
}


std::vector<double> Network::get_positions() const
{
    std::vector<double> positions(2*Nv);
    for(int i=0; i<2*Nv; ++i) {
        positions[i] = gsl_vector_get(r,i);
    }

    return positions;
}

void Network::savePositions(std::ostream& out) const
{
    for(int i=0; i<Nv; ++i) {
        out << gsl_vector_get(r, 2*i) << '\t';
        out << gsl_vector_get(r, 2*i+1) << '\n';
    }
}


double Network::edgeEnergy() const
{
	double e = 0;
    for(int ei=0; ei<Ne; ++ei) {
        e += get_edgeEnergy(ei,r);
    }
	return e/(Lx*Ly);
}

double Network::bendEnergy() const
{
    double e = 0;
	for(int bi=0; bi<Nb; ++bi) {
		e += get_bendEnergy(bi,r);
	}
    return e/(Lx*Ly);
}

double Network::totalEnergy() const
{
    double e = 0;
    for(int ei=0; ei<Ne; ++ei) {
        e += get_edgeEnergy(ei,r);
    }
	for(int bi=0; bi<Nb; ++bi) {
		e += get_bendEnergy(bi,r);
	}
    return e/(Lx*Ly);
}

std::vector<double> Network::stress2() const
{
	double sigmaXX = 0;
	double sigmaXY = 0;
	double sigmaYX = 0;
	double sigmaYY = 0;

	gsl_vector *dH;
    dH = gsl_vector_alloc(2*Nv);
    gsl_vector_scale( dH, 0); // set all elements to 0
	
	double xi,yi,xj,yj, xk,yk;
    for( int ei=0; ei< get_Nedges(); ++ei ) {

        get_edgeDEnergy(ei, r, dH);
	
		xi = gsl_vector_get(r, 2*edges[ei].i);	
		yi = gsl_vector_get(r, 2*edges[ei].i + 1);	

		xj = gsl_vector_get(r, 2*edges[ei].j);	
		yj = gsl_vector_get(r, 2*edges[ei].j + 1);	
	
		sigmaXX += gsl_vector_get(dH, 2*edges[ei].i) * (xj- xi);
		sigmaXY += gsl_vector_get(dH, 2*edges[ei].i) * (yj - yi);
		sigmaYX += gsl_vector_get(dH, 2*edges[ei].i+1) * (xj - xi);
		sigmaYY += gsl_vector_get(dH, 2*edges[ei].i+1) * (yj - yi);
	
		sigmaXX += gsl_vector_get(dH, 2*edges[ei].j) * (xi - xj);
		sigmaXY += gsl_vector_get(dH, 2*edges[ei].j) * (yi - yj);	
		sigmaYX += gsl_vector_get(dH, 2*edges[ei].j+1) * (xi - xj);	
		sigmaYY += gsl_vector_get(dH, 2*edges[ei].j+1) * (yi - yj);	

		gsl_vector_set(dH, 2*edges[ei].i, 0);
		gsl_vector_set(dH, 2*edges[ei].i + 1, 0);
		gsl_vector_set(dH, 2*edges[ei].j, 0);
		gsl_vector_set(dH, 2*edges[ei].j + 1, 0);
		//gsl_vector_scale( dH, 0); // set all elements to 0
    } 

	for( int bi=0; bi < get_Nbends(); ++bi) {
		get_bendDEnergy(bi, r, dH);

		xi = gsl_vector_get(r, 2*bends[bi].i);	
		yi = gsl_vector_get(r, 2*bends[bi].i + 1);	

		xj = gsl_vector_get(r, 2*bends[bi].j);	
		yj = gsl_vector_get(r, 2*bends[bi].j + 1);	

		xk = gsl_vector_get(r, 2*bends[bi].k);	
		yk = gsl_vector_get(r, 2*bends[bi].k + 1);	

	
		sigmaXX += gsl_vector_get(dH, 2*bends[bi].i) * (xj- xi);	
		sigmaXY += gsl_vector_get(dH, 2*bends[bi].i) * (yj - yi);	
		sigmaYX += gsl_vector_get(dH, 2*bends[bi].i+1) * (xj - xi);	
		sigmaYY += gsl_vector_get(dH, 2*bends[bi].i+1) * (yj - yi);	
	
	
		sigmaXX += gsl_vector_get(dH, 2*bends[bi].j) * (xi- xj)/2;	
		sigmaXY += gsl_vector_get(dH, 2*bends[bi].j) * (yi - yj)/2;	
		sigmaYX += gsl_vector_get(dH, 2*bends[bi].j+1) * (xi - xj)/2;	
		sigmaYY += gsl_vector_get(dH, 2*bends[bi].j+1) * (yi - yj)/2;	
		
		sigmaXX += gsl_vector_get(dH, 2*bends[bi].j) * (xk- xj)/2;	
		sigmaXY += gsl_vector_get(dH, 2*bends[bi].j) * (yk - yj)/2;	
		sigmaYX += gsl_vector_get(dH, 2*bends[bi].j+1) * (xk - xj)/2;	
		sigmaYY += gsl_vector_get(dH, 2*bends[bi].j+1) * (yk - yj)/2;	
		
	
		sigmaXX += gsl_vector_get(dH, 2*bends[bi].k) * (xj- xk);	
		sigmaXY += gsl_vector_get(dH, 2*bends[bi].k) * (yj - yk);	
		sigmaYX += gsl_vector_get(dH, 2*bends[bi].k+1) * (xj - xk);	
		sigmaYY += gsl_vector_get(dH, 2*bends[bi].k+1) * (yj - yk);	
	
		gsl_vector_set(dH, 2*bends[bi].i, 0);
		gsl_vector_set(dH, 2*bends[bi].i + 1, 0);
		gsl_vector_set(dH, 2*bends[bi].j, 0);
		gsl_vector_set(dH, 2*bends[bi].j + 1, 0);
		gsl_vector_set(dH, 2*bends[bi].k, 0);
		gsl_vector_set(dH, 2*bends[bi].k + 1, 0);
		//gsl_vector_scale( dH, 0); // set all elements to 0
	}

	std::vector<double> sigma(4);

	sigma[0] = sigmaXX/(2*Lx*Ly);
	sigma[1] = sigmaXY/(2*Lx*Ly);
	sigma[2] = sigmaYX/(2*Lx*Ly);
	sigma[3] = sigmaYY/(2*Lx*Ly);

	gsl_vector_free(dH);

	return sigma;
}

std::vector<double> Network::stress() const
{
	double sigmaXX = 0;
	double sigmaXY = 0;
	double sigmaYX = 0;
	double sigmaYY = 0;

	gsl_vector *dH;
    dH = gsl_vector_alloc(2*Nv);
    gsl_vector_scale( dH, 0); // set all elements to 0

    for( int ei=0; ei< get_Nedges(); ++ei ) {
        get_edgeDEnergy(ei, r, dH);
    } 

	for( int bi=0; bi < get_Nbends(); ++bi) {
		get_bendDEnergy(bi, r, dH);
	}

	for( int vi=0; vi<  Nv; ++vi) {
		sigmaXX += gsl_vector_get(dH, 2*vi    ) * gsl_vector_get(r,2*vi);
		sigmaXY += gsl_vector_get(dH, 2*vi    ) * gsl_vector_get(r,2*vi + 1);
		sigmaYX += gsl_vector_get(dH, 2*vi + 1) * gsl_vector_get(r,2*vi);
		sigmaYY += gsl_vector_get(dH, 2*vi + 1) * gsl_vector_get(r,2*vi + 1);
	}
	
	std::vector<double> sigma(4);
	sigma[0] = sigmaXX/(Lx*Ly);
	sigma[1] = sigmaXY/(Lx*Ly);
	sigma[2] = sigmaYX/(Lx*Ly);
	sigma[3] = sigmaYY/(Lx*Ly);

	gsl_vector_free(dH);

	return sigma;
}

void Network::set_Lx(double Lxx)
{ Lx = Lxx; }

void Network::set_Ly(double Lyy)
{ Ly = Lyy; }



///////////////////////////
// Edge Member functions //
///////////////////////////

double Network::Edge::get_l( const gsl_vector *r, const Network *net) const
{

    double g = net->gamma*net->Ly;
    double dx = gsl_vector_get(r, 2*i) - gsl_vector_get(r,2*j) - (net->Lx)*(1+net->epsilonX)*xb - yb*g;
    double dy = gsl_vector_get(r, 2*i+1) - gsl_vector_get(r,2*j+1) - (net->Ly)*(1+net->epsilonY)*yb;

    return std::sqrt( dx*dx + dy*dy);
}

double Network::Edge::energy(const gsl_vector *r,const Network *net) const
{
    double g = net->gamma*net->Ly;
    double dx = gsl_vector_get(r, 2*i  ) - gsl_vector_get(r,2*j  ) - (net->Lx)*(1+net->epsilonX)*xb - yb*g;
    double dy = gsl_vector_get(r, 2*i+1) - gsl_vector_get(r,2*j+1) - (net->Ly)*(1+net->epsilonY)*yb;


    double l = std::sqrt( dx*dx + dy*dy);
    l -= l0;
    return 0.5*l*l;
}

void Network::Edge::dEnergy( const gsl_vector *r, gsl_vector *df, const Network *net) const
{
    double g = (net->gamma)*(net->Ly);
    double dx = gsl_vector_get(r, 2*i) - gsl_vector_get(r,2*j) - (net->Lx)*(1+net->epsilonX)*xb - yb*g;
    double dy = gsl_vector_get(r, 2*i+1) - gsl_vector_get(r,2*j+1) - (net->Ly)*(1+net->epsilonY)*yb;
    double l = std::sqrt( dx*dx + dy*dy);
	
    dx *= 1-l0/l;
    dy *= 1-l0/l;

    double temp = gsl_vector_get(df, 2*i) + dx;
    gsl_vector_set(df, 2*i, temp );

    temp = gsl_vector_get(df, 2*j) - dx;
    gsl_vector_set(df, 2*j, temp );

    temp = gsl_vector_get(df, 2*i+1) + dy;
    gsl_vector_set(df, 2*i+1, temp);

    temp = gsl_vector_get(df, 2*j+1) - dy;
    gsl_vector_set(df, 2*j+1, temp);

}

double Network::Edge::energy_dEnergy(const gsl_vector *r, gsl_vector *df, const Network *net) const
{

    double g = (net->gamma)*(net->Ly);
    double dx = gsl_vector_get(r, 2*i) - gsl_vector_get(r,2*j) - (net->Lx)*(1+net->epsilonX)*xb - yb*g;
    double dy = gsl_vector_get(r, 2*i+1) - gsl_vector_get(r,2*j+1) - (net->Ly)*(1+net->epsilonY)*yb;
    double l = std::sqrt( dx*dx + dy*dy);

    dx *= 1-l0/l;
    dy *= 1-l0/l;

    double temp = gsl_vector_get(df, 2*i) + dx;
    gsl_vector_set(df, 2*i, temp );

    temp = gsl_vector_get(df, 2*j) - dx;
    gsl_vector_set(df, 2*j, temp );

    temp = gsl_vector_get(df, 2*i+1) + dy;
    gsl_vector_set(df, 2*i+1, temp);

    temp = gsl_vector_get(df, 2*j+1) - dy;
    gsl_vector_set(df, 2*j+1, temp);

    l -= l0;
    return 0.5*l*l;

}


///////////////////////////
// Bend Member Functions //
///////////////////////////

double Network::Bend::get_lji( const gsl_vector *r, const Network *net) const
{
	double xi = gsl_vector_get(r, 2*i);
	double yi = gsl_vector_get(r, 2*i+1);
	
	double xj = gsl_vector_get(r, 2*j);
	double yj = gsl_vector_get(r, 2*j+1);

	double dxji = xi - xj + (net->Lx)*(1+net->epsilonX)*xib + (net->Ly)*(net->gamma)*yib;
	double dyji = yi - yj + (net->Ly)*(1+net->epsilonY)*yib;

	return std::sqrt( dxji*dxji + dyji*dyji );
}

double Network::Bend::get_ljk( const gsl_vector *r, const Network *net) const
{
	double xk = gsl_vector_get(r, 2*k);
	double yk = gsl_vector_get(r, 2*k+1);
	
	double xj = gsl_vector_get(r, 2*j);
	double yj = gsl_vector_get(r, 2*j+1);

	double dxjk = xk - xj + (net->Lx)*(1+net->epsilonX)*xkb + (net->Ly)*(net->gamma)*ykb;
	double dyjk = yk - yj + (net->Ly)*(1+net->epsilonY)*ykb;

	return std::sqrt( dxjk*dxjk + dyjk*dyjk );
}

void Network::Bend::set_phi0( const gsl_vector *r, const Network *net)
{ phi0 = get_phi(r,net); }


double Network::Bend::get_phi( const gsl_vector *r, const Network *net) const
{

	double xi = gsl_vector_get(r, 2*i);
	double yi = gsl_vector_get(r, 2*i+1);
	
	double xj = gsl_vector_get(r, 2*j);
	double yj = gsl_vector_get(r, 2*j+1);

	double xk = gsl_vector_get(r, 2*k);
	double yk = gsl_vector_get(r, 2*k+1);

	double dxji = xi - xj + (net->Lx)*(1+net->epsilonX)*xib + (net->Ly)*(net->gamma)*yib;	
	double dyji = yi - yj + (net->Ly)*(1+net->epsilonY)*yib;
	double dxjk = xk - xj + (net->Lx)*(1+net->epsilonX)*xkb + (net->Ly)*(net->gamma)*ykb;
	double dyjk = yk - yj + (net->Ly)*(1+net->epsilonY)*ykb;

	double a = dyji*dxjk - dxji*dyjk;
	double b = dxji*dxjk + dyji*dyjk;
	double phi = std::atan2(a,b);
	if( phi < 0) phi += 2*std::acos(-1);
	return phi;
}


double Network::Bend::energy( const gsl_vector *r, const Network *net) const
{
	//if( ykb != 0 or yib != 0 or xkb != 0 or xib!=0) return 0;
	double phi = get_phi(r,net);
	double delta_phi = phi - phi0;
	return kappa*delta_phi*delta_phi/2;	
}

void Network::Bend::dEnergy ( const gsl_vector *r, gsl_vector *df, const Network *net) const
{

	
	//if( ykb != 0 or yib != 0 or xkb != 0 or xib!=0) return;
	double xi = gsl_vector_get(r, 2*i);
	double yi = gsl_vector_get(r, 2*i+1);
	
	double xj = gsl_vector_get(r, 2*j);
	double yj = gsl_vector_get(r, 2*j+1);

	double xk = gsl_vector_get(r, 2*k);
	double yk = gsl_vector_get(r, 2*k+1);

	double dxji = xi - xj + (net->Lx)*(1+net->epsilonX)*xib + (net->Ly)*(net->gamma)*yib;	
	double dyji = yi - yj + (net->Ly)*(1+net->epsilonY)*yib;
	double dxjk = xk - xj + (net->Lx)*(1+net->epsilonX)*xkb + (net->Ly)*(net->gamma)*ykb;
	double dyjk = yk - yj + (net->Ly)*(1+net->epsilonY)*ykb;

	double a = dyji*dxjk - dxji*dyjk;
	double b = dxji*dxjk + dyji*dyjk;
	double phi =  std::atan2(a,b);
	if( phi < 0 ) phi += 2*std::acos(-1);

	double Falpha =  kappa*(phi-phi0);
	double A =    b/(a*a+b*b);
	double B = -1*a/(a*a+b*b);

	double F;
	// set F_xi
	F = gsl_vector_get(df, 2*i);
	gsl_vector_set(df, 2*i, F+Falpha*(-A*dyjk + B*dxjk) );

	// set F_yi
	F = gsl_vector_get(df, 2*i+1);
	gsl_vector_set(df, 2*i+1, F+Falpha*(A*dxjk + B*dyjk) );

	// set F_xj
	F = gsl_vector_get(df, 2*j);
	gsl_vector_set(df, 2*j, F+Falpha*(A*(dyjk-dyji) - B*(dxjk + dxji) ) );

	// set F_yj
	F = gsl_vector_get(df, 2*j+1);
	gsl_vector_set(df, 2*j+1, F+Falpha*(A*(dxji-dxjk)-B*(dyjk+dyji) ) );

	// set F_xk
	F = gsl_vector_get(df, 2*k);
	gsl_vector_set(df, 2*k, F+Falpha*(A*dyji +B*dxji ) );

	// set F_yk
	F = gsl_vector_get(df, 2*k+1);
	gsl_vector_set(df, 2*k+1, F+Falpha*(-A*dxji+B*dyji) );


}

double Network::Bend::energy_dEnergy( const gsl_vector *r, gsl_vector *df, const Network *net) const
{
	dEnergy(r,df,net);
	return energy(r,net);
}


//////////////////////
// Energy Functions //
//////////////////////

double energy( const gsl_vector *v, void *params)
{
    Network *net = (Network *) params;
    double e=0;

	// energy from stretching bonds
    for(int ei=0; ei < net->get_Nedges(); ++ei) {
        e += net->get_edgeEnergy(ei, v);
    }

	// energy from bending
	for(int bi=0; bi < net->get_Nbends(); ++bi) {
		e += net->get_bendEnergy(bi, v);
	}
	
    return e;
}

void dEnergy( const gsl_vector *v, void *params, gsl_vector *df)
{
    
    Network *net = (Network *) params; 
    gsl_vector_scale( df, 0); // set all elements to 0
    for( int ei=0; ei< net->get_Nedges(); ++ei ) {
        net->get_edgeDEnergy(ei, v, df);
    } 

	for( int bi=0; bi < net->get_Nbends(); ++bi) {
		net->get_bendDEnergy(bi, v, df);
	}
}


void EdE( const gsl_vector *v, void *params, double *f, gsl_vector *df)
{
    *f = energy(v,params);
    dEnergy( v, params, df);
}





#endif
