#include <iostream>
#include <vector>
#include <string>


#include "gsl/gsl_multimin.h"

#include "minimize.h"

using namespace std;


double my_f(const gsl_vector *v, void *params)
{
  double x,y;
  double *p = (double *)params;
 
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);

  return p[2] * (x - p[0]) * (x - p[0]) +
         p[3] * (y - p[1]) * (y - p[1]) + p[4];
}



void my_df(const gsl_vector *v, void *params, gsl_vector *df)
{
  double x,y;
  double *p = (double *)params;

  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);

  gsl_vector_set(df, 0, 2.0 * p[2] * (x - p[0]) );
  gsl_vector_set(df, 1, 2.0 * p[3] * (y - p[1]) );
}

void my_fdf(const gsl_vector *x, void *params,
            double *f, gsl_vector *df)
{
  *f = my_f(x, params);
  my_df(x, params, df);
}
int main()
{


  size_t iter = 0; 
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  double par[5] = { 1., 2., 10., 20., 30. }; 

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  Network net(2);

  my_func.n = 2;
  my_func.f = my_f;
  my_func.df = my_df;
  my_func.fdf = my_fdf;
  my_func.params = par;


  x = gsl_vector_alloc(2);
  gsl_vector_set(x, 0, 5.0);
  gsl_vector_set(x, 1, 7.0);

  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc(T,2);
  
  gsl_multimin_fdfminimizer_set(s, &my_func, x, 10., 1e-4);

  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);

    if (status) break;

    status = gsl_multimin_test_gradient( s->gradient, 1e-3);

    if (status == GSL_SUCCESS) {
      cout << "MIN FOUND AT: \n";
    }

    cout << iter << "\t"
         << gsl_vector_get(s->x, 0) << "\t"
         << gsl_vector_get(s->x, 1) << "\t"
         << s->f << endl;


  } while (status == GSL_CONTINUE && iter < 100);


  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  return 0;
}
