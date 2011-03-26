#include <stdlib.h>

#include "config.h"
#include "m4ri.h"
#include "cpucycles.h"
#include "benchmarketing.h"

struct pluq_params {
  rci_t m;
  rci_t n;
  rci_t r;
};

int run(void *_p, double *wt, unsigned long long *cycles) {
  struct pluq_params *p = (struct pluq_params *)_p;

  mzd_t *A = mzd_init(p->m, p->n);
  mzd_t *U;
  mzd_t *L;

  int halfrank = 1;
  if(halfrank) {
    U = mzd_init(p->m, p->n);
    L = mzd_init(p->m, p->m);
    mzd_randomize(U);
    mzd_randomize(L);
#if 0
    for (rci_t i = 0; i < p->m; ++i) {
      for (rci_t j = i + 1; j < p->m; ++j)
        mzd_write_bit(L,i,j, 0);
      mzd_write_bit(L,i,i, 1);
    }
    for(rci_t i = 0; i < MIN(p->m, p->n); ++i) {
      for (rci_t j = 0; j < i; ++j)
        mzd_write_bit(U,i,j, 0);
      mzd_write_bit(U,i,i, 1);
    }
#endif
    for(rci_t i = 0; i < p->m; ++i) {
      mzd_write_bit(U,i,i, 1);
      for(rci_t j = 0; j < i; ++j)
	mzd_write_bit(U,i,j, 0);
      if ((i % 2))
	for(rci_t j = i; j < p->n; ++j)
	  mzd_write_bit(U,i,j, 0);
      for(rci_t j = i + 1; j < p->m; ++j)
	mzd_write_bit(L,i,j, 0);
      mzd_write_bit(L,i,i, 1);
    }

    mzd_mul(A,L,U,0);

  } else {
    mzd_randomize(A);
  }

  mzp_t *P = mzp_init(p->m);
  mzp_t *Q = mzp_init(p->n);

  *wt = walltime(0.0);
  *cycles = cpucycles();
  p->r = mzd_pluq(A, P, Q, 0);
  *wt = walltime(*wt);
  *cycles = cpucycles() - *cycles;

  mzd_free(A);
  mzp_free(P);
  mzp_free(Q);
  if(halfrank) {
    mzd_free(U);
    mzd_free(L);
  }

  return 0;
}

int main(int argc, char **argv) {
  global_options(&argc, &argv);

  if (argc != 3) {
    m4ri_die("Parameters m, n expected.\n");
  }

  struct pluq_params p;
  p.m = atoi(argv[1]);
  p.n = atoi(argv[2]);

  srandom(17);
  unsigned long long t;
  double wt;
  run_bench(run, (void*)&p, &wt, &t);

  printf("m: %5d, n: %5d, r: %5d, cpu cycles: %12llu, wall time: %6.3lf\n", p.m, p.n, p.r, t, wt);
}

