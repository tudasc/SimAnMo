#include "sampler.h"
#include <time.h>

inline double my_rand() {

  // return double(drand48());
  return double((rand()) / RAND_MAX);
}

double cachedExp(double arg, double *old_arg, double *old_val) {

  if (arg == *old_arg) {
    return *old_val;
  }
  if (arg == 0)
    return 1.0;

  *old_arg = arg;
  *old_val = exp(arg);

  return *old_val;
}

int roundI(double x) {
  if (x < 0.0)
    return (int)(x - 0.5);
  else
    return (int)(x + 0.5);
}

int SampleZ(double c, double s_square, double t_, int seed,
            drand48_data *private_status) {

  // To speedup the exp calculation
  double argument, old_argument = 0, old_value = 1;

  // printf("%f \n",c);
  // getchar();

  double aleatorio;
  double *alptx = &aleatorio;

  const double s = sqrt(s_square);
  const long minimum_c = floor(c - s * t_);
  const long maximum_c = ceil(c + s * t_);
  int register x;
  double register rho;
  while (1) {

    drand48_r(private_status, alptx);

    x = minimum_c + roundI((maximum_c - minimum_c) * aleatorio);

    argument = -M_PI * (x - c) * (x - c) / s_square;

    // printf("Argument %f Old %f Old_val %f \n", argument, old_argument,
    // old_value);
    // getchar();

    rho = // exp(- M_PI * (x - c) * (x - c) / s_square);
        // exp(argument);//
        cachedExp(argument, &old_argument, &old_value);
    // cachedExp( argument, old_argument, old_value );
    // old_value = rho;
    // old_argument = argument;

    drand48_r(private_status, alptx);

    if (aleatorio <= rho) {
      return x;
    }
  }

  /*

          unsigned int se = seed;
    double s = sqrt(s_square);
    long minimum_c = floor(c - s * t_);
    long maximum_c = ceil(c + s * t_);
    long x;
    double rho;
    while(true) {
      x = minimum_c +
        round((maximum_c - minimum_c) * double(rand_r(&se)));
      rho = exp( - M_PI * (x - c) * (x - c) / s_square);
      if ((double(rand_r(&se)) <= rho)) {
        return x;
      }
    }
  */
}

void SampleLF(struct NodeC *pLF, int n_, int m_, double t_, double **coef_,
              double *s_prime_square_, double **mu_, long **B_, int seed,
              drand48_data *private_status) {

  // cout << "val unit " << seed << endl;

  int i, j;

  for (i = 0; i < n_; ++i)
    (*coef_)[i] = 0;

  for (i = n_ - 1; i >= 0; --i) {
    (*coef_)[i] =
        SampleZ((*coef_)[i], s_prime_square_[i], t_, seed, private_status);
    for (j = 0; j < i; j++) {
      (*coef_)[j] -= ((*coef_)[i] * mu_[i][j]);
    }
  }

  /* RESET VECTOR TO 0 */

  pLF->norm = 0;
  for (i = 0; i < m_; i++)
    pLF->data[i] = 0;

  /* INSERT VALUES ON THE VECTOR */

  for (i = 0; i < m_; i++) {
    for (j = 0; j < n_; j++) {
      pLF->data[i] += (*coef_)[j] * B_[j][i];
    }
    pLF->norm += pLF->data[i] * pLF->data[i];
  }

//  if (pLF->norm < 0) std::cout << "Sampling error\n";
}
