#ifndef GUARD_MINIMIZE_GSL_H
#define GUARD_MINIMIZE_GSL_H

#include "gsl/gsl_multimin.h"

double my_f(const gsl_vector *v, void *params)
{
  Network *network_ptr = (Network *)params;
  for (int i = 0; i < 2*(network_ptr->Nnode); ++i) {
    network_ptr->r[i] = gsl_vector_get(v, i);
  }
  return network_ptr->totalEnergy();
}

void my_df(const gsl_vector *v, void *params, gsl_vector* df)
{
  Network *network_ptr = (Network *)params;

  for (int i = 0; i < 2*(network_ptr->Nnode); ++i) {
    network_ptr->r[i] = gsl_vector_get(v, i);
  }

  std::vector<double> std_df( 2*(network_ptr->Nnode) );
  network_ptr->dE( std_df, network_ptr->r);

  for (int i=0; i<2*(network_ptr->Nnode); ++i) {
      gsl_vector_set(df, i, std_df[i]);
   }
}

void my_fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df)
{
  *f = my_f(v, params);
  my_df(v, params, df);
}


class MinimizerGSL {
 public:
  MinimizerGSL(Network *network, double eLine, double dLine, double e)
    : N(2*network->Nnode), eLine(eLine), dLine(dLine), e(e) {

      my_func.n = N;
      my_func.f = my_f;
      my_func.df = my_df;
      my_func.fdf = my_fdf;
      my_func.params = network;

  }

  void minimize( Network &network);
  int N;
  double eLine, dLine, e;
  gsl_multimin_function_fdf my_func;
};


void MinimizerGSL::minimize(Network &network)
{

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc(T, N);


  
  gsl_vector *x;
  x = gsl_vector_alloc(N);

  for (int i=0;i<N; ++i) {
    gsl_vector_set(x, i, network.r[i]);
  }

  gsl_multimin_fdfminimizer_set(s, &my_func, x, dLine, eLine);

  int iter = 0;
  int status;
  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);
   
    if (status) break; 
    status = gsl_multimin_test_gradient( s->gradient, e);

  } while (status == GSL_CONTINUE && iter < 10000);

  for (int i=0; i<N; ++i) {
    network.r[i] = gsl_vector_get(x, i);
  } 

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);
}


#endif
