#ifndef GUARD_MINIMIZE_H
#define GUARD_MINIMIZE_H

class Network {
 public:

  Network(int N): N(N), r(N), k(N) 
  {
     for (int i = 0; i< N; ++i) {
        r[i] = (i%5)*(i%3) + 0.5 *(i%3)*(i*4); 
        k[i] = 0.5 + (i%2);
      }
  }

  int N;
  std::vector<double> r;
  std::vector<double> k;


  double my_f(const gsl_vector *v, void *params)
  {
    double E = 0;
    for(int i = 0;i<N; ++i) {
      E += gsl_vector_get(v, i) * gsl_vector_get(v,i)*k[i]; 
    }
    return E;
  }

  void my_df(const gsl_vector *v, void *params, gsl_vector *df)
  {
    for(int i = 0;i<N; ++i) {
      gsl_vector_set(df, i, 2*gsl_vector_get(v,i)*k[i]);
    }

  }


};

#endif
