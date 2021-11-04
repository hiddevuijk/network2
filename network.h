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
  class Bond;
  class Bend;

 public:
  Network( Graph&, double Lx, double Ly, double kappa);


  void minimize(double e, double emax, double dt0, double dtmax,
                double dtmin, double finc, double fdec, int Nmin,
                double alpha0,  double falpha, double m);

	// deform Network
  void shear(double delta_gamma); 
  void shearAffine(double delta_gamma); 
  void stretchX(double delta_epsilonX); 
  void stretchXAffine(double delta_epsilonX); 
  void stretchY(double delta_epsilonY); 
  void stretchYAffine(double delta_epsilonY); 

  int get_Nnodes() const {return Nnode; }
  int get_Nbonds() const { return Nbond; }
  int get_Nbends() const { return Nbend; }

  double bondEnergy() const; 
  double bendEnergy() const; 
  double totalEnergy() const; 

	std::vector<double> stress() const;
	std::vector<double> stress2() const;

	void dE(std::vector<double> &F, const  std::vector<double> &r) const;

  double get_bondEnergy(int bi, const  std::vector<double> &r) const
         { return bonds[bi].energy(r, this); }

  void get_bondDEnergy(int bi, const std::vector<double> &r,
                       std::vector<double> &F) const
       { bonds[bi].dEnergy(r, F, this); }

  double get_bondEdE(int bi, const std::vector<double> &r,
                     std::vector<double> &F) const
         { return bonds[bi].energy_dEnergy(r, F, this); }

	double get_bendEnergy(int bi, const  std::vector<double> &r) const
         { return bends[bi].energy(r, this); }

	void get_bendDEnergy(int bi, const std::vector<double> &r,
                       std::vector<double> &F) const
       { bends[bi].dEnergy(r, F, this); }

	double get_bendEdE(int bi, const  std::vector<double> &r,
                     std::vector<double> &F) const 
         { return bends[bi].energy_dEnergy(r , F, this); }

	double get_bend_dxij( int bi, const std::vector<double>& r) const
		     { return bends[bi].get_dxij(r,this); }
	double get_bend_dyij( int bi, const std::vector<double>& r) const
         { return bends[bi].get_dyij(r,this); }
	double get_bend_dxkj( int bi, const std::vector<double>& r) const
         { return bends[bi].get_dxkj(r,this); }
	double get_bend_dykj( int bi, const std::vector<double>& r) const
		     { return bends[bi].get_dykj(r,this); }

	double get_bond_dxij( int bi, const std::vector<double>& r) const
		     { return bonds[bi].get_dxij(r,this); }
	double get_bond_dyij( int bi, const std::vector<double>& r) const
		     { return bonds[bi].get_dyij(r,this); }


	double get_forceNorm() const;

  std::vector<double> get_positions() const;
  void savePositions(std::ostream& out) const;

  void set_Lx(double Lxx);
  void set_Ly(double Lyy);

 private:
  class Bond{
    public:
      int i,j; 
      int xb, yb;
      double l0; 

      double get_l(const  std::vector<double> &r, const Network *net) const;

      // return xi - xj with PBC
      double get_dxij(const std::vector<double> &r, const Network *net) const;
      double get_dyij(const std::vector<double> &r, const Network *net) const;

      void set_l0(const std::vector<double> &r, const Network *net)
           { l0 = get_l(r,net); }
      void set_l0(double l00) { l0 = l00; }

      double energy(const std::vector<double> &r, const Network *net) const;
      void dEnergy(const std::vector<double> &r,
                   std::vector<double> &F, const Network *net) const;
      double energy_dEnergy(const std::vector<double> &r,
                            std::vector<double> &F, const Network *net) const;
  };

  class Bend{
   public:
	  int i,j,k;
		int xib, yib;
		int xkb, ykb;

		double kappa;
		double phi0;

		//  return xi - xj with PBC
		double get_dxij(const std::vector<double> &r, const Network *net) const;
		double get_dyij(const std::vector<double> &r, const Network *net) const;
		double get_dxkj(const std::vector<double> &r, const Network *net) const;
		double get_dykj(const std::vector<double> &r, const Network *net) const;

		double get_lji(const std::vector<double> &r, const Network *net) const;
		double get_ljk(const std::vector<double> &r, const Network *net) const;
		double get_phi(const std::vector<double> &r, const Network *net) const;

		void set_phi0(const std::vector<double> &r, const Network *net);
	

		double energy(const  std::vector<double> &r, const Network *net) const;
		void dEnergy (const std::vector<double> &r,  std::vector<double> &F,
                  const Network *net) const;
		double energy_dEnergy(const std::vector<double> &r, std::vector<double> &F,
                          const Network *net) const;
  };


 public: // remove
  void set_kappa();
	void set_phi0();
	void set_l0();
	void set_l0(double l00);

	double get_phi(int bi) const { return bends[bi].get_phi(r,this); };
	
  int Nnode,Nbond,Nbend;
  double Lx, Ly;

  double gamma;  // shear deformation
	double epsilonX; // bulk deformation in x
	double epsilonY; // bulk deformation in y

	Fire<Network> minimizer;

  std::vector<double> r;
  std::vector<Bond> bonds;
  std::vector<Bend> bends;
  double kappa;
};

///////////////////////////////
// Network Member Functions  //
///////////////////////////////

Network::Network(Graph& g, double Lxx, double Lyy, double kappaa)
: Nnode(g.Nnodes()), Nbond(g.Nbonds()), Nbend(g.Nbends()),
	Lx(Lxx), Ly(Lyy), gamma(0), epsilonX(0), epsilonY(0),
	minimizer(2*Nnode, this), r(2*Nnode)
{
    
  for(int ni=0; ni<Nnode; ++ni ) {
    vec2 p = g.getNodePosition(ni);
    r[2*ni]     = p.x;		
    r[2*ni + 1] = p.y;		
  }

  std::vector<std::vector<int> > bo = g.getBonds();
  Nbond = bo.size();
  bonds = std::vector<Bond>(Nbond);
  for(int bi=0; bi< Nbond; ++bi ) {
    bonds[bi].i  = bo[bi][0]; 
    bonds[bi].j  = bo[bi][1]; 
    bonds[bi].xb = bo[bi][2]; 
    bonds[bi].yb = bo[bi][3]; 
    //bonds[bi].l0 = bo[bi][4]; 
  }

  std::vector<std::vector<int> > be = g.getBends();
	Nbend = be.size();
	bends = std::vector<Bend>(Nbend);
	for(int bi=0; bi < Nbend; ++bi ){
    bends[bi].j   = be[bi][0];
		bends[bi].i   = be[bi][1];
		bends[bi].xib = be[bi][2];
		bends[bi].yib = be[bi][3];

		bends[bi].k   = be[bi][4];
		bends[bi].xkb = be[bi][5];
		bends[bi].ykb = be[bi][6];
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
	std::vector<double> f(2*Nnode, 0.0);

  for( int bi=0; bi< Nbond; ++bi ) {
    get_bondDEnergy(bi, r, f);
  } 


  for( int bi=0; bi < Nbend; ++bi) {
    get_bendDEnergy(bi, r, f);
  }

	double norm = 0;
	for(int ni=0; ni<2*Nnode; ++ ni) {
		norm += f[ni]*f[ni];
	}
		
	return norm;
}

void Network::set_kappa()
{
	double lji, ljk;
	for( int bi=0; bi< Nbend; ++bi ) {
		lji = bends[bi].get_lji(r,this);
		ljk = bends[bi].get_ljk(r,this);
		//bends[bi].kappa = kappa*(
    //       bends[bi].get_lji(r, this) + bends[bi].get_ljk(r, this) ) / 2;
		bends[bi].kappa = kappa*(lji + ljk) / 2;
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
	for(unsigned int bi=0; bi<bonds.size(); ++bi) {
		bonds[bi].set_l0(r,this);
	}
}

void Network::set_l0(double l00)
{
	for(unsigned int bi=0; bi<bonds.size(); ++bi) {
		bonds[bi].set_l0(l00);
	}
}

void Network::shearAffine(double delta_gamma)
{
  gamma += delta_gamma;
  //affine deformation
  for( int ni=0; ni<Nnode; ++ni ) {
    r[2*ni] += delta_gamma*r[2*ni+1];
  }
}


void Network::shear(double delta_gamma)
{ gamma += delta_gamma; }

void Network::stretchX(double delta_epsilonX )
{ epsilonX += delta_epsilonX; }

void Network::stretchXAffine(double delta_epsilonX)
{
	// add Affine deformation
	for( int ni = 0; ni < Nnode ; ++ni) { 
		r[2*ni] += delta_epsilonX*r[2*ni];
	}
	epsilonX += delta_epsilonX;
}

void Network::stretchY(double delta_epsilonY)
{ epsilonY += delta_epsilonY; }


void Network::stretchYAffine(double delta_epsilonY)
{

  // add Affine deformation
	for(int ni=0; ni < Nnode; ++ni ) {
		r[2*ni+1] += delta_epsilonY*r[2*ni+1];
	}

	epsilonY += delta_epsilonY;
}
void Network::minimize(double e, double emax, double dt0, double dtmax,
                       double dtmin, double finc, double fdec, int Nmin,
                       double alpha0, double falpha, double m)
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
  std::vector<double> positions(2*Nnode);
  for(int i=0; i<2*Nnode; ++i) {
    positions[i] = r[i];
  }

    return positions;
}

void Network::savePositions(std::ostream& out) const
{
  for(int i=0; i<Nnode; ++i) {
    out << r[2*i] << '\t' << r[2*i+1] << '\n';
  }
}


double Network::bondEnergy() const
{
	double e = 0;
  for(int bi=0; bi<Nbond; ++bi) {
    e += get_bondEnergy(bi,r);
  }
	return e/(Lx*Ly);
}

double Network::bendEnergy() const
{
  double e = 0;
  for(int bi=0; bi<Nbend; ++bi) {
		e += get_bendEnergy(bi,r);
	}
  return e/(Lx*Ly);
}

double Network::totalEnergy() const
{
  double e = 0;
  for(int bi=0; bi<Nbond; ++bi) {
    e += get_bondEnergy(bi,r);
  }
  for(int bi=0; bi<Nbend; ++bi) {
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

  std::vector<double> dH(2*Nnode, 0.0);

  double dxij, dxji, dyij, dyji, dxkj, dxjk, dykj, dyjk;
  int i, j, k;

  for( int bi=0; bi<get_Nbonds(); ++bi) {
    get_bondDEnergy(bi, r, dH);
    dxij = get_bond_dxij(bi,r);
    dyij = get_bond_dyij(bi,r);
    dxji = - dxij;
    dyji = - dyij;

    i = 2*bonds[bi].i;
    j = 2*bonds[bi].j;

    sigmaXX += dH[i  ] * dxji;
    sigmaXY += dH[i  ] * dyji;
    sigmaYX += dH[i+1] * dxji;
    sigmaYY += dH[i+1] * dyji;

    sigmaXX += dH[j  ] * dxij;
    sigmaXY += dH[j  ] * dyij;
    sigmaYX += dH[j+1] * dxij;
    sigmaYY += dH[j+1] * dyij;

    dH[i  ] = 0;
    dH[i+1] = 0;
    dH[j  ] = 0;
    dH[j+1] = 0;
  }

  for(int bi = 0; bi < get_Nbends(); ++bi) {
    get_bendDEnergy(bi, r, dH);
    dxij = get_bend_dxij(bi,r);
    dyij = get_bend_dyij(bi,r);
    dxkj = get_bend_dxkj(bi,r);
    dykj = get_bend_dykj(bi,r);
    dxji = -dxij;
    dyji = -dyij;
    dxjk = -dxkj;
    dyjk = -dykj;

    i = 2*bends[bi].i;
    j = 2*bends[bi].j;
    k = 2*bends[bi].k;

    sigmaXX += dH[i  ] * dxji;
    sigmaXY += dH[i  ] * dyji;
    sigmaYX += dH[i+1] * dxji;
    sigmaYY += dH[i+1] * dyji;

    sigmaXX += dH[j  ] * dxij;
    sigmaXY += dH[j  ] * dyij;
    sigmaYX += dH[j+1] * dxij;
    sigmaYY += dH[j+1] * dyij;

    sigmaXX += dH[j  ] * dxkj;
    sigmaXY += dH[j  ] * dykj;
    sigmaYX += dH[j+1] * dxkj;
    sigmaYY += dH[j+1] * dykj;

    sigmaXX += dH[k  ] * dxjk;
    sigmaXY += dH[k  ] * dyjk;
    sigmaYX += dH[k+1] * dxjk;
    sigmaYY += dH[k+1] * dyjk;

    dH[i  ] = 0;
    dH[i+1] = 0;
    dH[j  ] = 0;
    dH[j+1] = 0;
    dH[k  ] = 0;
    dH[k+1] = 0;
  } 


  std::vector<double> sigma(4,0);
  sigma[0] = sigmaXX/(2*Lx*Ly);
  sigma[1] = sigmaXY/(2*Lx*Ly);
  sigma[2] = sigmaYX/(2*Lx*Ly);
  sigma[3] = sigmaYY/(2*Lx*Ly);

  return sigma;
}
std::vector<double> Network::stress2() const
{
  double sigmaXX = 0;
  double sigmaXY = 0;
  double sigmaYX = 0;
  double sigmaYY = 0;

  std::vector<double> dH(2*Nnode, 0.0);

  double xi, yi, xj,yj,xk,yk;
  int i, j, k;

  for(int bi = 0; bi < get_Nbonds(); ++bi) {
    get_bondDEnergy(bi, r, dH);

    i = 2*bonds[bi].i;
    j = 2*bonds[bi].j;

    xi = r[i  ];
    yi = r[i+1];
    xj = r[j  ];
    yj = r[j+1];

    sigmaXX += dH[i  ] * (xj - xi);
    sigmaXY += dH[i  ] * (yj - yi);
    sigmaYX += dH[i+1] * (xj - xi);
    sigmaYY += dH[i+1] * (yj - yi);

    sigmaXX += dH[j  ] * (xi - xj);
    sigmaXY += dH[j  ] * (yi - yj);
    sigmaYX += dH[j+1] * (xi - xj);
    sigmaYY += dH[j+1] * (yi - yj);

    dH[i  ] = 0;
    dH[i+1] = 0;
    dH[j  ] = 0;
    dH[j+1] = 0;
  }

  for(int bi =0;bi< get_Nbends(); ++bi) {
    get_bendDEnergy(bi, r, dH);

    i = 2*bends[bi].i;
    j = 2*bends[bi].j;
    k = 2*bends[bi].k;

    xi = r[i  ];
    yi = r[i+1];
    xj = r[j  ];
    yj = r[j+1];
    xk = r[k  ];
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



    dH[i  ] = 0;
    dH[i+1] = 0;
    dH[j  ] = 0;
    dH[j+1] = 0;
    dH[k  ] = 0;
    dH[k+1] = 0;
  } 


  std::vector<double> sigma(4,0);
  sigma[0] = sigmaXX/(2*Lx*Ly);
  sigma[1] = sigmaXY/(2*Lx*Ly);
  sigma[2] = sigmaYX/(2*Lx*Ly);
  sigma[3] = sigmaYY/(2*Lx*Ly);

  return sigma;
}


void Network::dE(std::vector<double> &F, const std::vector<double> &r) const
{
  std::fill(F.begin(), F.end(), 0.0); // set all elements to 0
  for( int bi=0; bi< get_Nbonds(); ++bi ) {
    get_bondDEnergy(bi, r, F);
  } 

  for( int bi=0; bi < get_Nbends(); ++bi) {
    get_bendDEnergy(bi, r, F);
  }
  for(int i=0; i<2*Nnode; ++i) F[i] *= -1;
}


void Network::set_Lx(double Lxx)
{ Lx = Lxx; }

void Network::set_Ly(double Lyy)
{ Ly = Lyy; }



///////////////////////////
// Bond Member functions //
///////////////////////////

double Network::Bond::get_l(const std::vector<double> &r,
                            const Network *net) const
{
  double g = net->gamma*net->Ly;
	double dx = r[2*i  ] - r[2*j  ] - (net->Lx)*(1+net->epsilonX)*xb - yb*g;
	double dy = r[2*i+1] - r[2*j+1] - (net->Ly)*(1+net->epsilonY)*yb;
  return std::sqrt( dx*dx + dy*dy);
}

double Network::Bond::get_dxij(const std::vector<double> &r,
                               const Network *net) const
{
  double g = net->gamma*net->Ly;
	return r[2*i] - r[2*j] - (net->Lx)*(1+net->epsilonX)*xb - yb*g;
}

double Network::Bond::get_dyij(const std::vector<double> &r,
                               const Network *net) const
{ return r[2*i+1] - r[2*j+1] - (net->Ly)*(1+net->epsilonY)*yb; }

double Network::Bond::energy(const std::vector<double> &r,const Network *net) const
{
  double g = net->gamma*net->Ly;
	double dx = r[2*i  ] - r[2*j  ] - (net->Lx)*(1+net->epsilonX)*xb - yb*g;
	double dy = r[2*i+1] - r[2*j+1] - (net->Ly)*(1+net->epsilonY)*yb;

  double l = std::sqrt( dx*dx + dy*dy);
  l -= l0;
  return 0.5*l*l;
}

void Network::Bond::dEnergy(const std::vector<double> &r,
                            std::vector<double> &F, const Network *net) const
{
  double g = (net->gamma)*(net->Ly);
	double dx = r[2*i  ] - r[2*j  ] - (net->Lx)*(1+net->epsilonX)*xb - yb*g;
	double dy = r[2*i+1] - r[2*j+1] - (net->Ly)*(1+net->epsilonY)*yb;
  double l = std::sqrt( dx*dx + dy*dy);
	
  dx *= 1-l0/l;
  dy *= 1-l0/l;

	F[2*i  ] += dx;
	F[2*j  ] -= dx;
	F[2*i+1] += dy;
	F[2*j+1] -= dy;
}

double Network::Bond::energy_dEnergy(const std::vector<double> &r,
                                     std::vector<double> &F,
                                     const Network *net) const
{
  double g = (net->gamma)*(net->Ly);
	double dx = r[2*i  ] - r[2*j  ] - (net->Lx)*(1+net->epsilonX)*xb - yb*g;
	double dy = r[2*i+1] - r[2*j+1] - (net->Ly)*(1+net->epsilonY)*yb;
  double l = std::sqrt(dx*dx + dy*dy);

  dx *= 1-l0/l;
  dy *= 1-l0/l;

	F[2*i  ] += dx;
	F[2*j  ] -= dx;
	F[2*i+1] += dy;
	F[2*j+1] -= dy;

  l -= l0;
  return 0.5*l*l;
}


///////////////////////////
// Bend Member Functions //
///////////////////////////

double Network::Bend::get_dxij(const std::vector<double> &r,
                               const Network *net) const
{
	double xi = r[2*i];
	double xj = r[2*j];
	return xi - xj + (net->Lx)*(1+net->epsilonX)*xib + (net->Ly)*(net->gamma)*yib;
}

double Network::Bend::get_dyij(const std::vector<double> &r,
                               const Network *net) const
{
	double yi = r[2*i+1];
	double yj = r[2*j+1];
	return yi - yj + (net->Ly)*(1+net->epsilonY)*yib;
}


double Network::Bend::get_dxkj(const std::vector<double> &r,
                               const Network *net) const
{
	double xk = r[2*k];
	double xj = r[2*j];
	return xk - xj + (net->Lx)*(1+net->epsilonX)*xkb + (net->Ly)*(net->gamma)*ykb;
}

double Network::Bend::get_dykj(const std::vector<double> &r,
                               const Network *net) const
{
	double yk = r[2*k+1];
	double yj = r[2*j+1];
	return yk - yj + (net->Ly)*(1+net->epsilonY)*ykb;
}


double Network::Bend::get_lji(const std::vector<double> &r,
                              const Network *net) const
{
	double xi = r[2*i];
	double yi = r[2*i+1];
	double xj = r[2*j];
	double yj = r[2*j+1];
		
	double dxji = xi - xj + (net->Lx)*(1+net->epsilonX)*xib +
                (net->Ly)*(net->gamma)*yib;
	double dyji = yi - yj + (net->Ly)*(1+net->epsilonY)*yib;

	return std::sqrt( dxji*dxji + dyji*dyji );
}

double Network::Bend::get_ljk(const std::vector<double> &r,
                              const Network *net) const
{
	double xk = r[2*k  ];
	double yk = r[2*k+1];
	double xj = r[2*j  ];
	double yj = r[2*j+1];

	double dxjk = xk - xj + (net->Lx)*(1+net->epsilonX)*xkb +
                (net->Ly)*(net->gamma)*ykb;
	double dyjk = yk - yj + (net->Ly)*(1+net->epsilonY)*ykb;

	return std::sqrt( dxjk*dxjk + dyjk*dyjk );
}

void Network::Bend::set_phi0(const std::vector<double>  &r, const Network *net)
{ phi0 = get_phi(r,net); }


double Network::Bend::get_phi(const  std::vector<double> &r,
                              const Network *net) const
{
	double xi = r[2*i  ];
	double yi = r[2*i+1];

	double xj = r[2*j  ];
	double yj = r[2*j+1];

	double xk = r[2*k  ];
	double yk = r[2*k+1];

	double dxji = xi - xj + (net->Lx)*(1+net->epsilonX)*xib +
                (net->Ly)*(net->gamma)*yib;	
	double dyji = yi - yj + (net->Ly)*(1+net->epsilonY)*yib;
	double dxjk = xk - xj + (net->Lx)*(1+net->epsilonX)*xkb +
                (net->Ly)*(net->gamma)*ykb;
	double dyjk = yk - yj + (net->Ly)*(1+net->epsilonY)*ykb;

	double a = dyji*dxjk - dxji*dyjk;
	double b = dxji*dxjk + dyji*dyjk;
	double phi = std::atan2(a,b);
	if (phi < 0) phi += 2*std::acos(-1);

	return phi;
}


double Network::Bend::energy(const std::vector<double> &r,
                             const Network *net) const
{
	//if( ykb != 0 or yib != 0 or xkb != 0 or xib!=0) return 0;
	double phi = get_phi(r,net);
	double delta_phi = phi - phi0;
	return kappa*delta_phi*delta_phi/2;	
}

void Network::Bend::dEnergy(const std::vector<double> &r,
                            std::vector<double> &F, const Network *net) const
{
	double xi = r[2*i  ];
	double yi = r[2*i+1];
	double xj = r[2*j  ];
	double yj = r[2*j+1];

	double xk = r[2*k  ];
	double yk = r[2*k+1];

	double dxji = xi - xj + (net->Lx)*(1+net->epsilonX)*xib +
                (net->Ly)*(net->gamma)*yib;	
	double dyji = yi - yj + (net->Ly)*(1+net->epsilonY)*yib;
	double dxjk = xk - xj + (net->Lx)*(1+net->epsilonX)*xkb +
                (net->Ly)*(net->gamma)*ykb;
	double dyjk = yk - yj + (net->Ly)*(1+net->epsilonY)*ykb;

	double a = dyji*dxjk - dxji*dyjk;
	double b = dxji*dxjk + dyji*dyjk;
	double phi =  std::atan2(a,b);
	if( phi < 0 ) phi += 2*std::acos(-1);

	double Falpha =  kappa*(phi-phi0);
	double A =    b/(a*a+b*b);
	double B = -1*a/(a*a+b*b);


	// set F
	F[2*i  ] += Falpha*(-A*dyjk + B*dxjk );
	F[2*i+1] += Falpha*( A*dxjk + B*dyjk );
	F[2*j  ] += Falpha*( A*(dyjk-dyji) - B*(dxjk+dxji) );
	F[2*j+1] += Falpha*( A*(dxji-dxjk) - B*(dyjk+dyji) );
	F[2*k  ] += Falpha*( A*dyji + B*dxji );
	F[2*k+1] += Falpha*(-A*dxji + B*dyji );
}

double Network::Bend::energy_dEnergy(const  std::vector<double> &r,
                                     std::vector<double> &F,
                                     const Network *net) const
{
	dEnergy(r,F,net);
	return energy(r,net);
}





#endif
