#include <stdlib.h>

#include "config.h"
#include "m4ri.h"
#include "cpucycles.h"
#include "benchmarketing.h"

struct elim_params {
  rci_t m;
  rci_t n;
  rci_t r;
  char const *algorithm;  
};

int run(void *_p, double *wt, unsigned long long *cycles) {
  struct elim_params *p = (struct elim_params *)_p;

  mzd_t *A = mzd_init(p->m, p->n);

  mzd_t *L, *U;
  L = mzd_init(p->m, p->m);
  U = mzd_init(p->m, p->n);
  mzd_randomize(U);
  mzd_randomize(L);
  for (rci_t i = 0; i < p->m; ++i) {
    mzd_write_bit(U,i,i, 1);
    for (rci_t j = 0; j < i; ++j)
      mzd_write_bit(U,i,j, 0);
    if (i >= p->r)
      for (rci_t j = i; j < p->n; ++j)
        mzd_write_bit(U,i,j, 0);
    for (rci_t j = i + 1; j < p->m; ++j)
      mzd_write_bit(L,i,j, 0);
    mzd_write_bit(L,i,i, 1);
  }
  mzd_mul(A,L,U,0);
  mzd_free(L);
  mzd_free(U);

  *wt = walltime(0.0);
  *cycles = cpucycles();
  if(strcmp(p->algorithm, "m4ri") == 0)
    p->r = mzd_echelonize_m4ri(A, 1, 0);
  else if(strcmp(p->algorithm, "pluq") == 0)
    p->r = mzd_echelonize_pluq(A, 1);
  else if(strcmp(p->algorithm, "mmpf") == 0)
    p->r = _mzd_pluq_mmpf(A, mzp_init(A->nrows), mzp_init(A->ncols), 0);
  else if(strcmp(p->algorithm, "naive") == 0)
    p->r = mzd_echelonize_naive(A, 1);
  *cycles = cpucycles() - *cycles;
  *wt = walltime(*wt);
  mzd_free(A);
  return 0;
}

int main(int argc, char **argv) {
  global_options(&argc, &argv);

  if (argc < 3) {
    m4ri_die("Parameters m,n (alg,r) expected.\n");
  }

  struct elim_params params;
  params.m = atoi(argv[1]);
  params.n = atoi(argv[2]);

  if (argc >= 4)
    params.algorithm = argv[3];
  else
    params.algorithm = "m4ri";
  if (argc == 5)
    params.r = atoi(argv[4]);
  else
    params.r = params.m;

  /* put this call in run() to benchmark one particular matrix over
     and over again instead of computing the average of various
     matrices.*/
  srandom(17);
  unsigned long long t;
  double wt;
  run_bench(run, (void*)&params, &wt, &t);

  printf("m: %5d, n: %5d, last r: %5d, cpu cycles: %10llu, wall time: %lf\n", params.m, params.n, params.r, t, wt);
}
