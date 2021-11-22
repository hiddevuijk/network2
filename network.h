#ifndef GUARD_NETWORK_H
#define GUARD_NETWORK_H

//
// TO DO:
// - Check definition/calculation of angles and differences of angles.
//

#define _USE_MATH_DEFINES

#include "graph.h"

#include <vector>
#include <math.h>


class Bond {
 public: 
  // Bond between node i and node j.
  int i, j;
  // If i->j goes over the right (left) x boundary, xb=1 (xb=-1).
  int xb, yb;

  // Rest length
  double l0;
};

class Bend {
 public:
  // Bend between nodes i-j-k.
  int i, j, k;

  // If j->i goes over the right (left) x boundary, xib=1 (xib=-1).
  // If j->i does not cross a boundary xib = 0.
  int xib, yib;

  // If j->k goes over the right (left) x boundary, xkb=1 (xkb=-1).
  // If j->k does not cross a boundary xkb = 0.
  int xkb, ykb;

  // Rest angle of the bend. The angle is defined as the clockwise rotation
  // of the vector j->i to the vector j->k.
  double phi0;
  // Bending constant
  double kappa;
};

class Network {
 public:
  Network(Graph &graph, double Lx, double Ly, double kappa);

  int getNumberOfNodes() const { return number_of_nodes_; }
  int getNumberOfBonds() const { return number_of_bonds_; }
  int getNumberOfBends() const { return number_of_bends_; }

  void shear(double delta_gamma) { gamma_ += delta_gamma; }
  // Increments gamma and performs affine shear deformation.
  void shearAffine(double delta_gamma);

  void stretchX(double delta_epsilon_x) { epsilon_x_ += delta_epsilon_x; }
  // Increments epsilon_x_ and perform affine deformation of positions_.
  void stretchXAffine(double delta_epsilon_x);

  void stretchY(double delta_epsilon_y) { epsilon_y_ += delta_epsilon_y; }
  // Increments epsilon_y_ and perform affine deformation of positions_.
  void stretchYAffine(double delta_epsilon_y);

  // Sets rest length such that current positions correspond to the rest length.
  void resetRestLength();
  // Sets rest angle such that the current positions correspond to the rest
  // angle.
  void resetRestAngle();
  // Sets kappa of each bend object to kappa * (lij + lkj)/2, where lij and 
  // lkj are the arm lengths.
  void resetKappa(double kappa);

  // Returns the sum of energies stored in all bonds.
  double getBondEnergy() const;
  double getBondEnergy(const std::vector<double> &positions) const;
  // Returns the sum of energies stored in all bends.
  double getBendEnergy() const;
  double getBendEnergy(const std::vector<double> &positions) const;
  // Returns the sum of the previous two functions.
  double getTotalEnergy() const { return getBondEnergy() + getBendEnergy();}
  double getTotalEnergy(const std::vector<double> &positions) const
          { return getBondEnergy(positions) + getBendEnergy(positions);}

  void moveNode(int i, double x, double y)
  { positions_[2*i] = x; positions_[2*i+1] = y; }

  // Adds bond force of bond i to F.
  // The force in the x (y) directions on node i is F[2*ni] (F[2*ni + 1]).
  void bondForce(int ni, std::vector<double> &F) const;    
  void bondForce(int ni, std::vector<double> &F,
                 const std::vector<double> & positions) const;
  // Adds all bond forces to F.
  // The force in the x (y) directions on node i is F[2*i] (F[2*i + 1]).
  void bondForce(std::vector<double> &F) const;    
  void bondForce(std::vector<double> &F,
                 const std::vector<double> &positions) const;    

  // Returns the Clockwise angle of bend i
  double getBendAngleCW(int i) const;
  double getBendAngleCW(int i, const std::vector<double> &positions) const;
  // Returns the angle of bend i
  double getBendAngle(int i) const;
  double getBendAngle(int i, const std::vector<double> &positions) const;

  void bendForce(int i, std::vector<double> &F) const;
  void bendForce(int i, std::vector<double> &F,
                 const std::vector<double> &positions) const;
  void bendForce(std::vector<double> &F) const;
  void bendForce(std::vector<double> &F,
                 const std::vector<double> &positions) const;

  // Replaces F with the total force on each degree of freedom.
  // The force in the x (y) directions on node i is F[2*i] (F[2*i + 1]).
  void force(std::vector<double> &F) const;
  void force(std::vector<double> &F,
             const std::vector<double> &positions) const;
  // Returns the  total force on each degree of freedom.
  std::vector<double> getForce() const;
  std::vector<double> getForce(const std::vector<double> &positions) const;

  void dE(std::vector<double> &F, const std::vector<double> &positions);

  std::vector<double> getStress() const;
  //std::vector<double> getStress(const std::vector<double> &positions) const;

  std::vector<double> getPositions() const { return positions_; }
  double getPosition(int i) const { return positions_[i]; }
  void setPositions(const std::vector<double> &p) { positions_ = p; }
  void setPosition(int i, double r) { positions_[i] = r; }
  void savePositions(std::ostream& out) const;

 //private:
  // Returns the distance between node i and j of bond bi
  double getBondDistance(int bi) const;
  double getBondDistance(int bi, const std::vector<double> &positions) const;
  double getBondEnergy(int i) const;
  double getBondEnergy(int i, const std::vector<double> &positions) const;

  // Returns the distance between node i and j of bend bi
  double getBendDistanceji(int bi) const;
  double getBendDistanceji(int bi, const std::vector<double> &positions) const;
  // Returns the distance between node k and j of bend bi
  double getBendDistancejk(int i) const;
  double getBendDistancejk(int i, const std::vector<double> &positions) const;
  double getBendDphi(int i) const;
  double getBendDphi(int i, const std::vector<double> &positions) const;
  double getBendEnergy(int i) const;
  double getBendEnergy(int i, const std::vector<double> &positions) const;

  double getGamma() const { return gamma_;}
  double getEpsilonX() const { return epsilon_x_;}
  double getEpsilonY() const { return epsilon_y_;}
  // Strain
  double gamma_;
  double epsilon_x_;
  double epsilon_y_;
 
  // Positions of the nodes.
  // The x and y position of node i are positions_[ 2*i ] and 
  // positions_[ 2*i + 1].
  std::vector<double> positions_;

  // x and y size of the system
  const double length_x_, length_y_;

  const int number_of_nodes_;
  const int number_of_bonds_;
  const int number_of_bends_;

  // bonds_ and bends_ are constant after the constructor is called.
  std::vector<Bond> bonds_;
  std::vector<Bend> bends_;

};


Network::Network(Graph &graph, double length_x, double length_y, double kappa)
  : gamma_(0), epsilon_x_(0), epsilon_y_(0), positions_(2*graph.Nnodes()),
    length_x_(length_x), length_y_(length_y), number_of_nodes_(graph.Nnodes()),
    number_of_bonds_(graph.Nbonds()), number_of_bends_(graph.Nbends()),
    bonds_(number_of_bonds_), bends_(number_of_bends_)
{

  // Store the positions, bonds and bends infromation from the graph
  // in the positions_, bonds_ and bends_ vector.
  // For the return objects, see graph.h.


  for (int ni = 0; ni < number_of_nodes_; ++ni) {
    vec2 r = graph.getNodePosition(ni);
    positions_[2*ni]     = r.x;
    positions_[2*ni + 1] = r.y;
  }

  std::vector<std::vector<int> > graph_bonds = graph.getBonds();
  for (int bi = 0; bi < number_of_bonds_; ++bi) {
    bonds_[bi].i  = graph_bonds[bi][0];
    bonds_[bi].j  = graph_bonds[bi][1];
    bonds_[bi].xb = graph_bonds[bi][2];
    bonds_[bi].yb = graph_bonds[bi][3];
  } 


  // Note:
  // The -1 in xib and yib are because in the graph they are defined 
  // as the crossing of the boundary from i to j, in the network (here)
  // they are defined as the crossing of the boundary from j to i.

  std::vector<std::vector<int> > graph_bends = graph.getBends();
  for (int bi = 0; bi < number_of_bends_; ++bi) {
    bends_[bi].j   = graph_bends[bi][0];
    bends_[bi].i   = graph_bends[bi][1];
    bends_[bi].xib = graph_bends[bi][2];
    bends_[bi].yib = graph_bends[bi][3];
    bends_[bi].k   = graph_bends[bi][4];
    bends_[bi].xkb = graph_bends[bi][5];
    bends_[bi].ykb = graph_bends[bi][6];
    bends_[bi].kappa = kappa;
  }

  resetRestLength();
  resetRestAngle();
  resetKappa(kappa);
}

void Network::shearAffine(double delta_gamma)
{
  gamma_ += delta_gamma;
  for (int i = 0; i < number_of_nodes_; ++i) {
    positions_[2*i] += delta_gamma * positions_[2*i + 1];
  }
}

void Network::stretchXAffine(double delta_epsilon_x)
{
  epsilon_x_ += delta_epsilon_x;
  for (int i = 0; i < number_of_nodes_; ++i) {
    positions_[2*i] += delta_epsilon_x * positions_[2*i];
  }
}

void Network::stretchYAffine(double delta_epsilon_y)
{
  epsilon_y_ += delta_epsilon_y;
  for (int i = 0; i < number_of_nodes_; ++i) {
    positions_[2*i + 1] += delta_epsilon_y * positions_[2*i + 1];
  }
}

 
void Network::resetRestLength()
{
  for (int i = 0; i < number_of_bonds_; ++i) {
    bonds_[i].l0 = getBondDistance(i); 
  }
}

void Network::resetRestAngle()
{
  for (int i = 0; i < number_of_bends_; ++i) {
    bends_[i].phi0 = getBendAngleCW(i);
  }
}

void Network::resetKappa(double kappa)
{
  double lji, ljk;
  for (int i = 0; i < number_of_bends_; ++i) {
    lji = getBendDistanceji(i);
    ljk = getBendDistancejk(i);
    bends_[i].kappa = kappa*(lji + ljk)/2;
  }
}

double Network::getBondEnergy() const
{
  double energy = 0;   
  for (int i = 0; i < number_of_bonds_; ++i) {
    energy += getBondEnergy(i);
  }
  return energy;
}

double Network::getBondEnergy(const std::vector<double> &positions) const
{
  double energy = 0;   
  for (int i = 0; i < number_of_bonds_; ++i) {
    energy += getBondEnergy(i, positions);
  }
  return energy;
}


//////////////////////////////
//
// private member functions
//
/////////////////////////////

double Network::getBondDistance(int bi) const
{
  // xi = positions_[i], yi = positions_[i + 1]
  int i = 2*bonds_[bi].i;
  int j = 2*bonds_[bi].j;

  // dx = xi - xj + "periodic bc"
  double dx = positions_[i] - positions_[j]
              - bonds_[bi].xb * length_x_ * (1 + epsilon_x_)
              - bonds_[bi].yb * length_y_ * gamma_ * (1 + epsilon_y_);

  double dy = positions_[i + 1] - positions_[j + 1]
              - bonds_[bi].yb * length_y_ * (1 + epsilon_y_);

  return sqrt(dx*dx + dy*dy);
}

double Network::getBondDistance(int bi,
                const std::vector<double> &positions) const
{
  // xi = positions_[i], yi = positions_[i + 1]
  int i = 2*bonds_[bi].i;
  int j = 2*bonds_[bi].j;

  // dx = xi - xj + "periodic bc"
  double dx = positions[i] - positions[j]
              - bonds_[bi].xb * length_x_ * (1 + epsilon_x_)
              - bonds_[bi].yb * length_y_ * gamma_ * (1 + epsilon_y_);

  double dy = positions[i + 1] - positions[j + 1]
              - bonds_[bi].yb * length_y_ * (1 + epsilon_y_);

  return sqrt(dx*dx + dy*dy);
}


double Network::getBondEnergy(int i) const
{
  double d = bonds_[i].l0 - getBondDistance(i);
  return d*d/2;
}

double Network::getBondEnergy(int i, const std::vector<double> &positions) const
{
  double d = bonds_[i].l0 - getBondDistance(i, positions);
  return d*d/2;
}

double Network::getBendDistanceji(int i) const
{
  // dxij = xi + "PBC" - xj 
  double dxij = positions_[2*bends_[i].i]
                + bends_[i].xib * length_x_ * (1 + epsilon_x_)
                + bends_[i].yib * length_y_ * gamma_ * (1 + epsilon_y_)
                - positions_[2*bends_[i].j];

  double dyij = positions_[2*bends_[i].i + 1]
                + bends_[i].yib * length_y_ * (1 + epsilon_y_)
                - positions_[2*bends_[i].j + 1];

  return sqrt(dxij * dxij + dyij * dyij);
}

double Network::getBendDistanceji(int i,
                const std::vector<double> &positions) const
{
  // dxij = xi + "PBC" - xj 
  double dxij = positions[2*bends_[i].i]
                + bends_[i].xib * length_x_ * (1 + epsilon_x_)
                + bends_[i].yib * length_y_ * gamma_ * (1 + epsilon_y_)
                - positions[2*bends_[i].j];

  double dyij = positions_[2*bends_[i].i + 1]
                + bends_[i].yib * length_y_ * (1 + epsilon_y_)
                - positions[2*bends_[i].j + 1];

  return sqrt(dxij * dxij + dyij * dyij);
}

double Network::getBendDistancejk(int i) const
{
  // dxkj = xk + "PBC" - xj 
  double dxkj = positions_[2*bends_[i].k]
                + bends_[i].xkb * length_x_ * (1 + epsilon_x_)
                + bends_[i].ykb * length_y_ * gamma_ * (1 + epsilon_y_)
                - positions_[2*bends_[i].j];

  double dykj = positions_[2*bends_[i].k + 1] 
                + bends_[i].ykb * length_y_ * (1 + epsilon_y_)
                - positions_[2*bends_[i].j + 1];

  return sqrt(dxkj * dxkj + dykj * dykj);
}

double Network::getBendDistancejk(int i,
                const std::vector<double> &positions) const
{
  // dxkj = xk + "PBC" - xj 
  double dxkj = positions[2*bends_[i].k]
                + bends_[i].xkb * length_x_ * (1 + epsilon_x_)
                + bends_[i].ykb * length_y_ * gamma_ * (1 + epsilon_y_)
                - positions[2*bends_[i].j];

  double dykj = positions[2*bends_[i].k + 1] 
                + bends_[i].ykb * length_y_ * (1 + epsilon_y_)
                - positions[2*bends_[i].j + 1];

  return sqrt(dxkj * dxkj + dykj * dykj);
}

double Network::getBendDphi(int i) const
{
  return getBendAngle(i) - bends_[i].phi0;
}
 
double Network::getBendDphi(int i, const std::vector<double> &positions) const
{
  return getBendAngle(i, positions) - bends_[i].phi0;
}
 
double Network::getBendEnergy() const
{
  double energy = 0;   
  for (int i = 0; i < number_of_bends_; ++i) {
    energy += getBendEnergy(i);
  }
  return energy;
}
 
double Network::getBendEnergy(const std::vector<double> &positions) const
{
  double energy = 0;   
  for (int i = 0; i < number_of_bends_; ++i) {
    energy += getBendEnergy(i, positions);
  }
  return energy;
}

double Network::getBendEnergy(int i) const
{
  double delta_phi = getBendAngle(i) - bends_[i].phi0;
  return bends_[i].kappa * delta_phi * delta_phi / 2;
}
 
double Network::getBendEnergy(int i, const std::vector<double> &positions) const
{
  double delta_phi = getBendAngle(i, positions) - bends_[i].phi0;
  return bends_[i].kappa * delta_phi * delta_phi / 2;
}

void Network::bondForce(int bi, std::vector<double> &F) const
{
  int i = 2 * bonds_[bi].i;
  int j = 2 * bonds_[bi].j;
  // dxij = xi + "PBC" - xj
  double dxij = positions_[i]
                - bonds_[bi].xb * length_x_ * (1 + epsilon_x_)
                - bonds_[bi].yb * length_y_ * gamma_ * (1 + epsilon_y_)
                - positions_[j];

  double dyij = positions_[i + 1]
                - bonds_[bi].yb * length_y_ * (1 + epsilon_y_)
                - positions_[j + 1];

  // f = l0/l - 1
  double f = bonds_[bi].l0/sqrt( dxij*dxij + dyij*dyij ) - 1;


  // F[xi] = (l0/l - 1 ) * (xi - xj)
  F[i] += f * dxij;
  F[j] -= f * dxij;

  F[i + 1] += f * dyij;
  F[j + 1] -= f * dyij;

}

void Network::bondForce(int bi, std::vector<double> &F,
              const std::vector<double> &positions) const
{
  int i = 2 * bonds_[bi].i;
  int j = 2 * bonds_[bi].j;

  // dxij = xi + "PBC" - xj
  double dxij = positions[i]
                - bonds_[bi].xb * length_x_ * (1 + epsilon_x_)
                - bonds_[bi].yb * length_y_ * gamma_ * (1 + epsilon_y_)
                - positions[j];

  double dyij = positions[i + 1]
                - bonds_[bi].yb * length_y_ * (1 + epsilon_y_)
                - positions[j + 1];

  // f = l0/l - 1
  double f = bonds_[bi].l0/sqrt( dxij*dxij + dyij*dyij ) - 1;


  // F[xi] = (l0/l - 1 ) * (xi - xj)
  F[i] += f * dxij;
  F[j] -= f * dxij;

  F[i + 1] += f * dyij;
  F[j + 1] -= f * dyij;
}

void Network::bondForce(std::vector<double> &F) const
{
  for (int i = 0; i < number_of_bonds_; ++i) {
    bondForce(i, F);
  }
}


void Network::bondForce(std::vector<double> &F,
              const std::vector<double> &positions) const
{
  for (int i = 0; i < number_of_bonds_; ++i) {
    bondForce(i, F, positions);
  }
}

double Network::getBendAngleCW(int bi) const
{
  int i = bends_[bi].i;
  int j = bends_[bi].j;
  int k = bends_[bi].k;

  // xi + "PBC"
  double xi = positions_[2*i]
              + bends_[bi].xib * length_x_ * (1 + epsilon_x_)
              + bends_[bi].yib * length_y_ * gamma_ * (1 + epsilon_y_);
  // yi + "PBC"
  double yi = positions_[2*i + 1]
              + bends_[bi].yib * length_y_ * (1 + epsilon_y_);

  double xj = positions_[2*j];
  double yj = positions_[2*j + 1];


  // xk + "PBC"
  double xk = positions_[2*k]
              + bends_[bi].xkb * length_x_ * (1 + epsilon_x_)
              + bends_[bi].ykb * length_y_ * gamma_ * (1 + epsilon_y_);

  // yk + "PBC"
  double yk = positions_[2*k + 1]
              + bends_[bi].ykb * length_y_ * (1 + epsilon_y_);

  double dxji = xi - xj;
  double dyji = yi - yj;

  double dxjk = xk - xj;
  double dyjk = yk - yj;

  double a = dyji * dxjk - dxji * dyjk;
  double b = dxji * dxjk + dyji * dyjk;
  double phi = std::atan2(a, b);
  if (phi < 0) phi += 2 * M_PI;

  return phi;
}

double Network::getBendAngleCW(int bi,
                const std::vector<double> &positions) const
{
  int i = 2*bends_[bi].i;
  int j = 2*bends_[bi].j;
  int k = 2*bends_[bi].k;

  // xi + "PBC"
  double xi = positions[i]
              + bends_[bi].xib * length_x_ * (1 + epsilon_x_)
              + bends_[bi].yib * length_y_ * gamma_ * (1 + epsilon_y_);
  // yi + "PBC"
  double yi = positions[i + 1]
              + bends_[bi].yib * length_y_ * (1 + epsilon_y_);

  double xj = positions[j];
  double yj = positions[j + 1];


  // xk + "PBC"
  double xk = positions[k]
              + bends_[bi].xkb * length_x_ * (1 + epsilon_x_)
              + bends_[bi].ykb * length_y_ * gamma_ * (1 + epsilon_y_);

  // yk + "PBC"
  double yk = positions[k + 1]
              + bends_[bi].ykb * length_y_ * (1 + epsilon_y_);

  double dxji = xi - xj;
  double dyji = yi - yj;

  double dxjk = xk - xj;
  double dyjk = yk - yj;

  double a = dyji * dxjk - dxji * dyjk;
  double b = dxji * dxjk + dyji * dyjk;
  double phi = std::atan2(a, b);
  if (phi < 0) phi += 2 * M_PI;

  return phi;
}

double Network::getBendAngle(int bi) const
{
  double phi = getBendAngleCW(bi);
  if (std::fabs(phi - bends_[bi].phi0) > M_PI) {
    if(phi > 0) phi -= 2 * M_PI;
  }

  return phi;
}

double Network::getBendAngle(int bi, const std::vector<double> &positions) const
{
  double phi = getBendAngleCW(bi, positions);
  if (std::fabs(phi - bends_[bi].phi0) > M_PI) {
    if(phi > 0) phi -= 2 * M_PI;
  }

  return phi;
}


void Network::bendForce(int bi, std::vector<double> &F) const
{
  int i = 2*bends_[bi].i;
  int j = 2*bends_[bi].j;
  int k = 2*bends_[bi].k;

   // xi + "PBC"
  double xi = positions_[i]
              + bends_[bi].xib * length_x_ * (1 + epsilon_x_)
              + bends_[bi].yib * length_y_ * gamma_ * (1 + epsilon_y_);
  // yi + "PBC"
  double yi = positions_[i + 1]
              + bends_[bi].yib * length_y_ * (1 + epsilon_y_);

  double xj = positions_[j];
  double yj = positions_[j + 1];


  // xk + "PBC"
  double xk = positions_[k]
              + bends_[bi].xkb * length_x_ * (1 + epsilon_x_)
              + bends_[bi].ykb * length_y_ * gamma_ * (1 + epsilon_y_);

  // yk + "PBC"
  double yk = positions_[k + 1]
              + bends_[bi].ykb * length_y_ * (1 + epsilon_y_);

  double dxji = xi - xj;
  double dyji = yi - yj;

  double dxjk = xk - xj;
  double dyjk = yk - yj;


  double a   = dyji * dxjk - dxji * dyjk;
  double b   = dxji * dxjk + dyji * dyjk;
  double phi = std::atan2(a, b);

  if (phi < 0) phi += 2 * M_PI;
  // phi is the clockwise angle between the vector j->i and j->k.

  if (std::fabs(phi - bends_[bi].phi0) > M_PI) {
    if(phi > 0) phi -= 2 * M_PI;
  }

	double Falpha =  bends_[bi].kappa * (phi - bends_[bi].phi0);
	double A =    b/(a*a + b*b);
	double B = -1*a/(a*a + b*b);

	// set F
	F[i  ] -= Falpha*(-A*dyjk + B*dxjk );
	F[i+1] -= Falpha*( A*dxjk + B*dyjk );
	F[j  ] -= Falpha*( A*(dyjk - dyji) - B*(dxjk + dxji) );
	F[j+1] -= Falpha*( A*(dxji - dxjk) - B*(dyjk + dyji) );
	F[k  ] -= Falpha*( A*dyji + B*dxji );
	F[k+1] -= Falpha*(-A*dxji + B*dyji );
}

void Network::bendForce(int bi, std::vector<double> &F,
              const std::vector<double> &positions) const
{
  int i = 2*bends_[bi].i;
  int j = 2*bends_[bi].j;
  int k = 2*bends_[bi].k;

   // xi + "PBC"
  double xi = positions[i]
              + bends_[bi].xib * length_x_ * (1 + epsilon_x_)
              + bends_[bi].yib * length_y_ * gamma_ * (1 + epsilon_y_);
  // yi + "PBC"
  double yi = positions[i + 1]
              + bends_[bi].yib * length_y_ * (1 + epsilon_y_);

  double xj = positions[j];
  double yj = positions[j + 1];


  // xk + "PBC"
  double xk = positions[k]
              + bends_[bi].xkb * length_x_ * (1 + epsilon_x_)
              + bends_[bi].ykb * length_y_ * gamma_ * (1 + epsilon_y_);

  // yk + "PBC"
  double yk = positions[k + 1]
              + bends_[bi].ykb * length_y_ * (1 + epsilon_y_);

  double dxji = xi - xj;
  double dyji = yi - yj;

  double dxjk = xk - xj;
  double dyjk = yk - yj;


  double a   = dyji * dxjk - dxji * dyjk;
  double b   = dxji * dxjk + dyji * dyjk;
  double phi = std::atan2(a, b);

  if (phi < 0) phi += 2 * M_PI;
  // phi is the clockwise angle between the vector j->i and j->k.

  if (std::fabs(phi - bends_[bi].phi0) > M_PI) {
    if(phi > 0) phi -= 2 * M_PI;
  }

	double Falpha =  bends_[bi].kappa * (phi - bends_[bi].phi0);
	double A =    b/(a*a + b*b);
	double B = -1*a/(a*a + b*b);

	// set F
	F[i  ] -= Falpha*(-A*dyjk + B*dxjk );
	F[i+1] -= Falpha*( A*dxjk + B*dyjk );
	F[j  ] -= Falpha*( A*(dyjk - dyji) - B*(dxjk + dxji) );
	F[j+1] -= Falpha*( A*(dxji - dxjk) - B*(dyjk + dyji) );
	F[k  ] -= Falpha*( A*dyji + B*dxji );
	F[k+1] -= Falpha*(-A*dxji + B*dyji );
}

void Network::bendForce(std::vector<double> & F) const
{
  for (int i = 0; i < number_of_bends_; ++i) {
    bendForce(i, F);
  }
}

void Network::bendForce(std::vector<double> & F,
              const std::vector<double> &positions) const
{
  for (int i = 0; i < number_of_bends_; ++i) {
    bendForce(i, F, positions);
  }
}


void Network::force(std::vector<double> &F) const 
{
  std::fill(F.begin(),F.end(), 0.0);
  bondForce(F);
  bendForce(F);
}

void Network::force(std::vector<double> &F, 
              const std::vector<double> &positions) const 
{
  std::fill(F.begin(),F.end(), 0.0);
  bondForce(F, positions);
  bendForce(F, positions);
}

std::vector<double> Network::getForce() const 
{
  std::vector<double> F(2 * number_of_nodes_, 0.0);
  bondForce(F);
  bendForce(F);

  return F;
}

std::vector<double> Network::getForce(
                    const std::vector<double> &positions) const 
{
  std::vector<double> F(2 * number_of_nodes_, 0.0);
  bondForce(F, positions);
  bendForce(F, positions);

  return F;
}

void Network::dE(std::vector<double> &F, const std::vector<double> &positions)
{
  force(F, positions);
  //for (int i = 0; i < 2 * number_of_nodes_; ++i) F[i] *= -1;
}



std::vector<double> Network::getStress() const
{
  double sigmaXX = 0;
  double sigmaXY = 0;
  double sigmaYX = 0;
  double sigmaYY = 0;

  std::vector<double> F(2*number_of_nodes_, 0.0);

  double dxij, dxji, dyij, dyji, dxkj, dxjk, dykj, dyjk;
  int i, j, k;

  for( int bi = 0; bi < number_of_bonds_; ++bi) {
    bondForce(bi, F);

    // xi = positions_[i], yi = positions_[i+1]
    i = 2*bonds_[bi].i;
    j = 2*bonds_[bi].j; 

    // dxij = xi + "PBC" - xj
    dxij = positions_[i]
           - bonds_[bi].xb * length_x_ * (1 + epsilon_x_)
           - bonds_[bi].yb * length_y_ * gamma_ * (1 + epsilon_y_)
           - positions_[j];

    dyij = positions_[i + 1]
           - bonds_[bi].yb * length_y_ * (1 + epsilon_y_)
           - positions_[j + 1];

    dxji = - dxij;
    dyji = - dyij;


    sigmaXX += F[i  ] * dxji;
    sigmaXY += F[i  ] * dyji;
    sigmaYX += F[i+1] * dxji;
    sigmaYY += F[i+1] * dyji;

    sigmaXX += F[j  ] * dxij;
    sigmaXY += F[j  ] * dyij;
    sigmaYX += F[j+1] * dxij;
    sigmaYY += F[j+1] * dyij;

    F[i  ] = 0;
    F[i+1] = 0;
    F[j  ] = 0;
    F[j+1] = 0;
  }

  for(int bi = 0; bi < number_of_bends_; ++bi) {
    bendForce(bi, F);

    i = 2*bends_[bi].i;
    j = 2*bends_[bi].j;
    k = 2*bends_[bi].k;

    // xi + "PBC"
    double xi = positions_[i]
                + bends_[bi].xib * length_x_ * (1 + epsilon_x_)
                + bends_[bi].yib * length_y_ * gamma_ * (1 + epsilon_y_);
    // yi + "PBC"
    double yi = positions_[i + 1]
                + bends_[bi].yib * length_y_ * (1 + epsilon_y_);

    double xj = positions_[j];
    double yj = positions_[j + 1];


    // xk + "PBC"
    double xk = positions_[k]
                + bends_[bi].xkb * length_x_ * (1 + epsilon_x_)
                + bends_[bi].ykb * length_y_ * gamma_ * (1 + epsilon_y_);

    // yk + "PBC"
    double yk = positions_[k + 1]
                + bends_[bi].ykb * length_y_ * (1 + epsilon_y_);

    dxij = xi - xj;
    dyij = yi - yj;
    dxkj = xk - xj;
    dykj = yk - yj;

    dxji = -dxij;
    dyji = -dyij;
    dxjk = -dxkj;
    dyjk = -dykj;


    sigmaXX += F[i  ] * dxji;
    sigmaXY += F[i  ] * dyji;
    sigmaYX += F[i+1] * dxji;
    sigmaYY += F[i+1] * dyji;

    sigmaXX += F[j  ] * dxij;
    sigmaXY += F[j  ] * dyij;
    sigmaYX += F[j+1] * dxij;
    sigmaYY += F[j+1] * dyij;

    sigmaXX += F[j  ] * dxkj;
    sigmaXY += F[j  ] * dykj;
    sigmaYX += F[j+1] * dxkj;
    sigmaYY += F[j+1] * dykj;

    sigmaXX += F[k  ] * dxjk;
    sigmaXY += F[k  ] * dyjk;
    sigmaYX += F[k+1] * dxjk;
    sigmaYY += F[k+1] * dyjk;

    F[i  ] = 0;
    F[i+1] = 0;
    F[j  ] = 0;
    F[j+1] = 0;
    F[k  ] = 0;
    F[k+1] = 0;
  } 


  std::vector<double> sigma(4,0);
  sigma[0] = sigmaXX/2;
  sigma[1] = sigmaXY/2;
  sigma[2] = sigmaYX/2;
  sigma[3] = sigmaYY/2;

  return sigma;
}

void Network::savePositions(std::ostream& out) const
{
  for(int i = 0; i < number_of_nodes_; ++i) {
    out << positions_[2*i] << '\t' << positions_[2*i+1] << '\n';
  }
}


#endif
