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

#include "fire2.h"

#include <math.h>


class Network
{
  private:
    class Edge;
    class Bend;

  public:
    Network(const Graph&, double Lx, double Ly, double kappa);


    void minimize( double e, double emax, double dt0, double dtmax,double dtmin, double finc, double fdec, int Nmin, double alpha0,  double falpha, double m);

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

	void dE( std::vector<double> &F, const  std::vector<double> &r) const;

    double get_edgeEnergy( int ei, const  std::vector<double> &r) const
        { return edges[ei].energy(r, this); }

    void get_edgeDEnergy( int ei, const std::vector<double> &r,  std::vector<double> &F) const
        {  edges[ei].dEnergy(r, F, this); }

    double get_edgeEdE( int ei, const std::vector<double> &r, std::vector<double> &F) const
        { return edges[ei].energy_dEnergy(r, F, this); }

	double get_bendEnergy( int bi, const  std::vector<double> &r) const
		{ return bends[bi].energy(r, this); }

	void get_bendDEnergy( int bi, const  std::vector<double> &r,  std::vector<double> &F) const
		{ bends[bi].dEnergy(r, F, this); }

	double get_bendEdE( int bi, const  std::vector<double> &r,  std::vector<double> &F) const 
		{ return bends[bi].energy_dEnergy(r , F, this); }

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

		double get_l( const  std::vector<double> &r, const Network *net) const;
		void set_l0( const  std::vector<double> &r, const Network *net) { l0 = get_l(r,net); }
		void set_l0( double l00) { l0 = l00; }

        double energy(const  std::vector<double> &r, const Network *net) const;
        void dEnergy( const  std::vector<double> &r,   std::vector<double> &F, const Network *net) const;
        double energy_dEnergy(const  std::vector<double> &r,  std::vector<double> &F, const Network *net) const;
    };

    class Bend{
	  public:
		int i,j,k;
		int xib, yib;
		int xkb, ykb;

		double kappa;
		double phi0;

		double get_lji( const std::vector<double> &r, const Network *net) const;
		double get_ljk( const std::vector<double> &r, const Network *net) const;
		double get_phi( const std::vector<double> &r, const Network *net) const;

		void set_phi0( const std::vector<double> &r, const Network *net);
	

		double energy( const  std::vector<double> &r, const Network *net) const;
		void dEnergy ( const std::vector<double> &r,  std::vector<double> &F, const Network *net) const;
		double energy_dEnergy( const std::vector<double> &r, std::vector<double> &F, const Network *net) const;
    };


  public: // remove
	void set_kappa();
	void set_phi0();
	void set_l0();
	void set_l0(double l00);

	double get_phi(int bi) const { return bends[bi].get_phi(r,this); };
	
    int Nv,Ne,Nb;
    double Lx, Ly;
	Fire<Network> minimizer;

    std::vector<double> r;
    std::vector<Edge> edges;
    std::vector<Bend> bends;
    double kappa;

    double gamma;  // shear deformation
	double epsilonX; // bulk deformation in x
	double epsilonY; // bulk deformation in y

};

///////////////////////////////
// Network Member Functions  //
///////////////////////////////

Network::Network(const Graph& g, double Lxx, double Lyy, double kappaa)
: Nv(g.Nvertices() ), Lx(Lxx), Ly(Lyy),  minimizer(2*Nv, this), r(2*Nv)
{
    
    for(int vi=0; vi<Nv; ++vi ) {
        vec2 p = g.getVertexPosition(vi);
		r[2*vi]     = p.x;		
		r[2*vi + 1] = p.y;		
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
	std::vector<double> f(2*Nv, 0.0);

    for( int ei=0; ei< Ne; ++ei ) {
        get_edgeDEnergy(ei, r, f);
    } 


	for( int bi=0; bi < Nb; ++bi) {
		get_bendDEnergy(bi, r, f);
	}

	double norm = 0;
	for(int vi=0; vi<2*Nv; ++ vi) {
		norm += f[vi]*f[vi];
	}
		
	return norm;
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
    //affine deformation
    for( int vi=0; vi<Nv; ++vi ) {
		r[2*vi] += delta_gamma*r[2*vi+1];
    }

}


void Network::shear( double delta_gamma)
{ gamma += delta_gamma; }

void Network::stretchX( double delta_epsilonX )
{ epsilonX += delta_epsilonX; }

void Network::stretchXAffine( double delta_epsilonX)
{
	// add Affine deformation
	for( int vi = 0; vi < Nv ; ++vi) { 
		r[2*vi] += delta_epsilonX*r[2*vi];
	}
	epsilonX += delta_epsilonX;
}

void Network::stretchY( double delta_epsilonY)
{ epsilonY += delta_epsilonY; }


void Network::stretchYAffine( double delta_epsilonY)
{

	// add Affine deformation
	for(int vi=0; vi < Nv; ++vi ) {
		r[2*vi+1] += delta_epsilonY*r[2*vi+1];
	}

	epsilonY += delta_epsilonY;
}
void Network::minimize( double e, double emax, double dt0, double dtmax,double dtmin, double finc, double fdec, int Nmin, double alpha0,  double falpha, double m)
{
	minimizer.error = e;
	minimizer.error_max = emax;
	minimizer.dt0 = dt0;
	minimizer.dtmax = dtmax;
    minimizer.dtmin = dtmin;
    minimizer.finc = finc;
    minimizer.fdec = fdec;
    minimizer.Nmin = Nmin;
    minimizer.alpha0 = alpha0;
    minimizer.falpha = falpha;
	minimizer.m = m;


	minimizer.minimizeVV(r);
	r = minimizer.x;
}


std::vector<double> Network::get_positions() const
{
    std::vector<double> positions(2*Nv);
    for(int i=0; i<2*Nv; ++i) {
        positions[i] = r[i];
    }

    return positions;
}

void Network::savePositions(std::ostream& out) const
{
    for(int i=0; i<Nv; ++i) {
        out << r[2*i] << '\t' << r[2*i+1] << '\n';
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

std::vector<double> Network::stress() const
{
    double sigmaXX = 0;
    double sigmaXY = 0;
    double sigmaYX = 0;
    double sigmaYY = 0;

    std::vector<double> dH(2*Nv, 0.0);

    double xi, yi, xj,yj,xk,yk;
    int i, j, k;

    for( int ei=0; ei<get_Nedges(); ++ei) {
        get_edgeDEnergy(ei, r, dH);

        i = 2*edges[ei].i;
        j = 2*edges[ei].j;

        xi = r[i];
        yi = r[i+1];
        xj = r[j];
        yj = r[j+1];

        sigmaXX += dH[i  ] * (xj - xi);
        sigmaXY += dH[i  ] * (yj - yi);
        sigmaYX += dH[i+1] * (xj - xi);
        sigmaYY += dH[i+1] * (yj - yi);
        
        sigmaXX += dH[j  ] * (xi - xj);
        sigmaXY += dH[j  ] * (yi - yj);
        sigmaYX += dH[j+1] * (xi - xj);
        sigmaYY += dH[j+1] * (yi - yj);

        dH[i] = 0;
        dH[i+1] = 0;
        dH[j] = 0;
        dH[j+1] = 0;
    }
   
    for(int bi =0;bi< get_Nbends(); ++bi) {
        get_bendDEnergy(bi, r, dH);

        i = 2*bends[bi].i;
        j = 2*bends[bi].j;
        k = 2*bends[bi].k;

        xi = r[i];
        yi = r[i+1];
        xj = r[j];
        yj = r[j+1];
        xk = r[k];
        yk = r[k+1];

        sigmaXX += dH[i  ] * (xj - xi);
        sigmaXY += dH[i  ] * (yj - yi);
        sigmaYX += dH[i+1] * (xj - xi);
        sigmaYY += dH[i+1] * (yj - yi);
        
        sigmaXX += dH[j  ] * (xi - xj);
        sigmaXY += dH[j  ] * (yi - yj);
        sigmaYX += dH[j+1] * (xi - xj);
        sigmaYY += dH[j+1] * (yi - yj);
 
        sigmaXX += dH[j  ] * (xk - xj);
        sigmaXY += dH[j  ] * (yk - yj);
        sigmaYX += dH[j+1] * (xk - xj);
        sigmaYY += dH[j+1] * (yk - yj);
 
        sigmaXX += dH[k  ] * (xj - xk);
        sigmaXY += dH[k  ] * (yj - yk);
        sigmaYX += dH[k+1] * (xj - xk);
        sigmaYY += dH[k+1] * (yj - yk);



        dH[i] = 0;
        dH[i+1] = 0;
        dH[j] = 0;
        dH[j+1] = 0;
        dH[k] = 0;
        dH[k+1] = 0;
    } 
    

    std::vector<double> sigma(4,0);
    sigma[0] = sigmaXX/(2*Lx*Ly);
    sigma[1] = sigmaXY/(2*Lx*Ly);
    sigma[2] = sigmaYX/(2*Lx*Ly);
    sigma[3] = sigmaYY/(2*Lx*Ly);
    return sigma;
}


void Network::dE( std::vector<double> &F, const  std::vector<double> &r) const
{
    
    
	std::fill(F.begin(), F.end(), 0.0); // set all elements to 0
    for( int ei=0; ei< get_Nedges(); ++ei ) {
        get_edgeDEnergy(ei, r, F);
    } 

	for( int bi=0; bi < get_Nbends(); ++bi) {
		get_bendDEnergy(bi, r, F);
	}
	for(int i=0; i<2*Nv; ++i) F[i] *= -1;
}


void Network::set_Lx(double Lxx)
{ Lx = Lxx; }

void Network::set_Ly(double Lyy)
{ Ly = Lyy; }



///////////////////////////
// Edge Member functions //
///////////////////////////

double Network::Edge::get_l( const std::vector<double> &r, const Network *net) const
{

    double g = net->gamma*net->Ly;
	double dx = r[2*i]   - r[2*j]   - (net->Lx)*(1+net->epsilonX)*xb - yb*g;
	double dy = r[2*i+1] - r[2*j+1] - (net->Ly)*(1+net->epsilonY)*yb;
    return std::sqrt( dx*dx + dy*dy);
}

double Network::Edge::energy(const std::vector<double> &r,const Network *net) const
{
    double g = net->gamma*net->Ly;
	double dx = r[2*i]   - r[2*j]   - (net->Lx)*(1+net->epsilonX)*xb - yb*g;
	double dy = r[2*i+1] - r[2*j+1] - (net->Ly)*(1+net->epsilonY)*yb;


    double l = std::sqrt( dx*dx + dy*dy);
    l -= l0;
    return 0.5*l*l;
}

void Network::Edge::dEnergy( const std::vector<double> &r, std::vector<double> &F, const Network *net) const
{
    double g = (net->gamma)*(net->Ly);
	double dx = r[2*i]   - r[2*j]   - (net->Lx)*(1+net->epsilonX)*xb - yb*g;
	double dy = r[2*i+1] - r[2*j+1] - (net->Ly)*(1+net->epsilonY)*yb;
    double l = std::sqrt( dx*dx + dy*dy);
	
    dx *= 1-l0/l;
    dy *= 1-l0/l;

	F[2*i]   += dx;
	F[2*j]   -= dx;
	F[2*i+1] += dy;
	F[2*j+1] -= dy;
}

double Network::Edge::energy_dEnergy(const std::vector<double> &r, std::vector<double> &F, const Network *net) const
{

    double g = (net->gamma)*(net->Ly);
	double dx = r[2*i]   - r[2*j]   - (net->Lx)*(1+net->epsilonX)*xb - yb*g;
	double dy = r[2*i+1] - r[2*j+1] - (net->Ly)*(1+net->epsilonY)*yb;
    double l = std::sqrt( dx*dx + dy*dy);

    dx *= 1-l0/l;
    dy *= 1-l0/l;

	F[2*i]   += dx;
	F[2*j]   -= dx;
	F[2*i+1] += dy;
	F[2*j+1] -= dy;

    l -= l0;
    return 0.5*l*l;

}


///////////////////////////
// Bend Member Functions //
///////////////////////////

double Network::Bend::get_lji( const std::vector<double> &r, const Network *net) const
{

	double xi = r[2*i];
	double yi = r[2*i+1];

	double xj = r[2*j];
	double yj = r[2*j+1];
		

	double dxji = xi - xj + (net->Lx)*(1+net->epsilonX)*xib + (net->Ly)*(net->gamma)*yib;
	double dyji = yi - yj + (net->Ly)*(1+net->epsilonY)*yib;

	return std::sqrt( dxji*dxji + dyji*dyji );
}

double Network::Bend::get_ljk( const std::vector<double> &r, const Network *net) const
{

	double xk = r[2*k];
	double yk = r[2*k+1];
	double xj = r[2*j];
	double yj = r[2*j+1];

	double dxjk = xk - xj + (net->Lx)*(1+net->epsilonX)*xkb + (net->Ly)*(net->gamma)*ykb;
	double dyjk = yk - yj + (net->Ly)*(1+net->epsilonY)*ykb;

	return std::sqrt( dxjk*dxjk + dyjk*dyjk );
}

void Network::Bend::set_phi0( const std::vector<double>  &r, const Network *net)
{ phi0 = get_phi(r,net); }


double Network::Bend::get_phi( const  std::vector<double> &r, const Network *net) const
{

	double xi = r[2*i];
	double yi = r[2*i+1];

	double xj = r[2*j];
	double yj = r[2*j+1];

	double xk = r[2*k];
	double yk = r[2*k+1];

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


double Network::Bend::energy( const std::vector<double> &r, const Network *net) const
{
	//if( ykb != 0 or yib != 0 or xkb != 0 or xib!=0) return 0;
	double phi = get_phi(r,net);
	double delta_phi = phi - phi0;
	return kappa*delta_phi*delta_phi/2;	
}

void Network::Bend::dEnergy ( const std::vector<double> &r, std::vector<double> &F, const Network *net) const
{

	
	double xi = r[2*i];
	double yi = r[2*i+1];

	double xj = r[2*j];
	double yj = r[2*j+1];

	double xk = r[2*k];
	double yk = r[2*k+1];

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


	// set F
	F[2*i]   += Falpha*(-A*dyjk + B*dxjk );
	F[2*i+1] += Falpha*( A*dxjk + B*dyjk );
	F[2*j]   += Falpha*( A*(dyjk-dyji) - B*(dxjk+dxji) );
	F[2*j+1] += Falpha*( A*(dxji-dxjk) - B*(dyjk+dyji) );
	F[2*k]   += Falpha*( A*dyji + B*dxji );
	F[2*k+1] += Falpha*(-A*dxji + B*dyji );

}

double Network::Bend::energy_dEnergy( const  std::vector<double> &r,  std::vector<double> &F, const Network *net) const
{
	dEnergy(r,F,net);
	return energy(r,net);
}







#endif
