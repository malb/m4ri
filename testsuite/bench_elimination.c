#include <stdlib.h>

#include "config.h"
#include "m4ri.h"
#include "cpucycles.h"
#include "walltime.h"

int main(int argc, char **argv) {
  if (argc < 3) {
    m4ri_die("Parameters m,n (alg,r) expected.\n");
  }

  rci_t m = atoi(argv[1]);
  rci_t n = atoi(argv[2]);

  rci_t r;
  char const *algorithm;
  if (argc >= 4)
    algorithm = argv[3];
  else
    algorithm = "m4ri";
  if (argc == 5)
    r = atoi(argv[4]);
  else
    r = m;

  mzd_t *A = mzd_init(m, n);
  mzd_t *L = mzd_init(m, m);
  mzd_t *U = mzd_init(m, n);
  mzd_randomize(U);
  mzd_randomize(L);
  for (rci_t i = 0; i < m; ++i) {
    mzd_write_bit(U,i,i, 1);
    for (rci_t j = 0; j < i; ++j)
      mzd_write_bit(U,i,j, 0);
    if (i >= r)
      for (rci_t j = i; j < n; ++j)
        mzd_write_bit(U,i,j, 0);
    for (rci_t j = i + 1; j < m; ++j)
      mzd_write_bit(L,i,j, 0);
    mzd_write_bit(L,i,i, 1);
  }
  mzd_mul(A,L,U,0);
  mzd_free(L);
  mzd_free(U);

  double clockZero = 0.0;
  double wt = walltime(&clockZero);
  unsigned long long t = cpucycles();
  if(strcmp(algorithm,"m4ri") == 0)
    r = mzd_echelonize_m4ri(A, 1, 0);
  else if(strcmp(algorithm,"pluq") == 0)
    r = mzd_echelonize_pluq(A, 1);
  else if(strcmp(algorithm,"mmpf") == 0)
    r = _mzd_pluq_mmpf(A, mzp_init(A->nrows), mzp_init(A->ncols), 0);
  else if(strcmp(algorithm,"naive") == 0)
    r = mzd_echelonize_naive(A, 1);
  printf("m: %5d, n: %5d, r: %5d, cpu cycles: %10llu wall time: %lf\n", m, n, r, cpucycles() - t, walltime(&wt));

  mzd_free(A);
}
