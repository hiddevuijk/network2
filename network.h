#ifndef GUARD_NETWORK_H
#define GUARD_NETWORK_H

//
// TO DO:
// - Check definition/calculation of angles and differences of angles.
//


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
  // Returns the sum of energies stored in all bends.
  double getBendEnergy() const;
  // Returns the sum of the previous two functions.
  double getTotalEnergy() const;

  void moveNode(int i, double x, double y)
  { positions_[2*i] = x; positions_[2*i+1] = y; }

  // Adds bond force of bond i to F.
  // The force in the x (y) directions on node i is F[2*i] (F[2*i + 1]).
  void bondForce(int i, std::vector<double> &F) const;    
  // Adds all bond forces to F.
  // The force in the x (y) directions on node i is F[2*i] (F[2*i + 1]).
  void bondForce(std::vector<double> &F) const;    

  // Returns the angle of bend i
  double getBendAngle(int i) const;

  void bendForce(int i, std::vector<double> &F) const;
  void bendForce(std::vector<double> & F) const;

  // Replaces F with the total force on each degree of freedom.
  // The force in the x (y) directions on node i is F[2*i] (F[2*i + 1]).
  void force(std::vector<double> &F) const;
  // Returns the  total force on each degree of freedom.
  std::vector<double> getForce() const;

  std::vector<double> getPositions() const { return positions_; }

 //private:
  // Returns the distance between node i and j of bond bi
  double getBondDistance(int bi) const;
  double getBondEnergy(int i) const;

  // Returns the distance between node i and j of bend bi
  double getBendDistanceji(int bi) const;
  // Returns the distance between node k and j of bend bi
  double getBendDistancejk(int i) const;
  double getBendDphi(int i) const;
  double getBendEnergy(int i) const;

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
    bends_[bi].xib = -1*graph_bends[bi][2];
    bends_[bi].yib = -1*graph_bends[bi][3];
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
    bends_[i].phi0 = getBendAngle(i);
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

//////////////////////////////
//
// private member functions
//
/////////////////////////////

double Network::getBondDistance(int i) const
{
  // dx = xi - xj + "periodic bc"
  double dx = positions_[2*bonds_[i].i] - positions_[2*bonds_[i].j]
              - bonds_[i].xb * length_x_ * (1 + epsilon_x_)
              - bonds_[i].yb * length_y_ * gamma_;

  double dy = positions_[2*bonds_[i].i + 1] - positions_[2*bonds_[i].j + 1]
              - bonds_[i].yb * length_y_ * (1 + epsilon_y_);

  return sqrt(dx*dx + dy*dy);
}

double Network::getBondEnergy(int i) const
{
  double d = bonds_[i].l0 - getBondDistance(i);
  return d*d/2;
}

double Network::getBendDistanceji(int i) const
{
  // dxij = xi + "PBC" - xj 
  double dxij = positions_[2*bends_[i].i]
                + bends_[i].xib * length_x_ * (1 + epsilon_x_)
                + bends_[i].yib * length_y_ * gamma_
                - positions_[2*bends_[i].j];

  double dyij = positions_[2*bends_[i].i + 1]
                + bends_[i].yib * length_y_ * (1 + epsilon_y_)
                - positions_[2*bends_[i].j + 1];

  return sqrt(dxij * dxij + dyij * dyij);
}

double Network::getBendDistancejk(int i) const
{
  // dxkj = xk + "PBC" - xj 
  double dxkj = positions_[2*bends_[i].k]
                + bends_[i].xkb * length_x_ * (1 + epsilon_x_)
                + bends_[i].ykb * length_y_ * gamma_
                - positions_[2*bends_[i].j];

  double dykj = positions_[2*bends_[i].k + 1] 
                + bends_[i].ykb * length_y_ * (1 + epsilon_y_)
                - positions_[2*bends_[i].j + 1];

  return sqrt(dxkj * dxkj + dykj * dykj);
}

double Network::getBendDphi(int i) const
{
  return getBendAngle(i) - bends_[i].phi0;
}
  
double Network::getBendEnergy(int i) const
{
  double delta_phi = getBendAngle(i) - bends_[i].phi0;
  return bends_[i].kappa * delta_phi * delta_phi / 2;
}

void Network::bondForce(int i, std::vector<double> &F) const
{
  // dxij = xi + "PBC" - xj
  double dxij = positions_[2*bonds_[i].i]
                + bonds_[i].xb * length_x_ * (1 + epsilon_x_)
                + bonds_[i].yb * length_y_ * gamma_
                - positions_[2*bonds_[i].j];

  double dyij = positions_[2*bonds_[i].i + 1]
                + bonds_[i].yb * length_y_ * (1 + epsilon_y_)
                - positions_[2*bonds_[i].j + 1];

  // f = l0/l - 1
  double f = bonds_[i].l0/sqrt( dxij*dxij + dyij*dyij ) - 1;


  // F[xi] = (l0/l - 1 ) * (xi - xj)
  F[2*bonds_[i].i] += f * dxij;
  F[2*bonds_[i].j] -= f * dxij;

  F[2*bonds_[i].i + 1] += f * dyij;
  F[2*bonds_[i].j + 1] -= f * dyij;

}
void Network::bondForce(std::vector<double> &F) const
{
  for (int i = 0; i < number_of_bonds_; ++i) {
    bondForce(i, F);
  }
}


double Network::getBendAngle(int i) const
{
  // xi + "PBC"
  double xi = positions_[2*bends_[i].i]
              + bends_[i].xib * length_x_ * (1 + epsilon_x_)
              + bends_[i].yib * length_y_ * gamma_;
  // yi + "PBC"
  double yi = positions_[2*bends_[i].i + 1]
              + bends_[i].yib * length_y_ * (1 + epsilon_y_);

  double xj = positions_[2*bends_[i].j];
  double yj = positions_[2*bends_[i].j + 1];


  // xk + "PBC"
  double xk = positions_[2*bends_[i].k]
              + bends_[i].xkb * length_x_ * (1 + epsilon_x_)
              + bends_[i].ykb * length_y_ * gamma_;

  // yk + "PBC"
  double yk = positions_[2*bends_[i].k + 1]
              + bends_[i].ykb * length_y_ * (1 + epsilon_y_);

  double dxji = xi - xj;
  double dyji = yi - yj;

  double dxjk = xk - xj;
  double dyjk = yk - yj;

  double a = dyji * dxjk - dxji * dyjk;
  double b = dxji * dxjk + dyji * dyjk;
  double phi = std::atan2(a, b);
  if (phi < 0) phi += 2 * std::acos(-1);

  return phi;
}

void Network::bendForce(int i, std::vector<double> &F) const
{

}

void Network::bendForce(std::vector<double> & F) const
{
  for (int i = 0; i < number_of_bonds_; ++i) {
    bendForce(i, F);
  }
}

void Network::force(std::vector<double> &F) const 
{
  std::fill(F.begin(),F.end(), 0.0);
  bondForce(F);
  bendForce(F);
}

std::vector<double> Network::getForce() const 
{
  std::vector<double> F(2 * number_of_nodes_, 0.0);
  bondForce(F);
  bendForce(F);

  return F;
}


#endif
