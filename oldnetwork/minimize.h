#ifndef GUARD_MINIMIZE_H
#define GUARD_MINIMIZE_H

#include "fire2.h"

#include <vector>


class Minimizer {
 public:
  Minimizer(Network &network, double error, double maxError, double dt0,
            double dtmax, double dtmin, double finc, double fdec, double Nmin,
            double alpha0, double falpha, double m)
      : fire(2*network.Nnode, &network) {

    fire.error = error;
    fire.error_max = maxError;
    fire.dt0 = dt0;
    fire.dtmax = dtmax;
    fire.dtmin = dtmin;
    fire.finc = finc;
    fire.fdec = fdec;
    fire.Nmin = Nmin;
    fire.alpha0 = alpha0;
    fire.falpha = falpha;
    fire.m = m;


  }

  void minimize(Network &network);

  Fire<Network> fire;

};

void Minimizer::minimize(Network &network) 
{
  fire.minimizeVV(network.r);  
  network.r = fire.x;
}

#endif
