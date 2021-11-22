#ifndef GUARD_MINIMIZE_GSL_H
#define GUARD_MINIMIZE_GSL_H

#include "gsl/gsl_multimin.h"


std::vector<double> gsl2vec(const gsl_vector *v)
{
  std::vector<double> vec(v->size);
  for (unsigned int i = 0; i < v->size; ++i) {
    vec[i] = gsl_vector_get(v, i);
  }
  return vec;
}


double my_f(const gsl_vector *v, void *params)
{
  Network *network_ptr = (Network *)params;
  std::vector<double> positions = gsl2vec(v);
  return network_ptr->getTotalEnergy(positions);
}

void my_df(const gsl_vector *v, void *params, gsl_vector *df)
{
  Network *network_ptr = (Network *)params;

  std::vector<double> positions = gsl2vec(v);

  std::vector<double> std_df(2 * network_ptr->getNumberOfNodes(), 0.0);
  network_ptr->dE(std_df, positions);

  // -1 because dE calculates the forces
  for (unsigned int i = 0; i < std_df.size(); ++i) {
     gsl_vector_set(df, i, -1*std_df[i]);
  }
}

void my_fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df)
{
  *f = my_f(v, params);
  my_df(v, params, df);
}


class MinimizerGSL {
 public:
  MinimizerGSL(Network *network, double  eLine, double dLine, double e, int maxIter)
    : N(2 * network->getNumberOfNodes()), eLine(eLine), dLine(dLine), e(e),
      maxIter(maxIter) {

    my_func.n   = N;
    my_func.f   = my_f;
    my_func.df  = my_df;
    my_func.fdf = my_fdf;
    my_func.params = network;

  }

  void minimize(Network &network);
  int N;
  double eLine, dLine, e;
  int maxIter;
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

  for (int i = 0; i < N; ++i) {
    gsl_vector_set(x, i, network.getPosition(i));
  }

  gsl_multimin_fdfminimizer_set(s, &my_func, x, dLine, eLine);

  int iter = 0;
  int status;

  do {
    ++iter;
    status = gsl_multimin_fdfminimizer_iterate(s);

    if (status) break;
    status = gsl_multimin_test_gradient(s->gradient, e);

  } while (status == GSL_CONTINUE && iter < maxIter);

  // set positions network
  for (int i = 0; i < N; ++i) {
    network.setPosition(i, gsl_vector_get(s->x, i));
  }

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);
}



#endif
