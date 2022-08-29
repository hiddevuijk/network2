#ifndef GUARD_GET_GAMMA_LIST_H
#define GUARD_GET_GAMMA_LIST_H

#include <vector>

std::vector<double> getGammaList( double dg0, double gmax, double gc, double alpha)
{
  double dg = dg0;
  double g = 0;
  std::vector<double> gamma_list;
  gamma_list.push_back(0);
  while (g < 0.8 * gc) {
    g += dg;
    gamma_list.push_back(g);
    dg *= alpha;
  }
  dg = dg0/2;
  while (g < 1.025 * gc) {
    g += dg;
    gamma_list.push_back(g);
    dg *= alpha;
  }
  if( dg < dg0) dg = dg0;
  while (g < 0.5) {
    g += dg;
    gamma_list.push_back(g);
    dg *= alpha;
  }
  alpha = 1 + 5 * (alpha - 1);
  while (g < 2 * gmax) {
    g += dg;
    gamma_list.push_back(g);
    dg *= alpha;
  }
  int gi = 0;
  while (gamma_list[gi] < gc) ++gi;
  gamma_list.insert( gamma_list.begin() + gi, gc); 
  gamma_list[gi-1] = ( gamma_list[gi-2] + gc)/2.;
  gamma_list[gi+1] = ( gamma_list[gi+2] + gc)/2.;
  return gamma_list;
}

#endif // GUARD_GET_GAMMA_LIST_H
