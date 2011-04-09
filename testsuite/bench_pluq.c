#include <stdlib.h>

#include "config.h"
#include "m4ri.h"
#include "cpucycles.h"
#include "benchmarking.h"

struct pluq_params {
  rci_t m;
  rci_t n;
  rci_t r;
};

int run(void *_p, unsigned long long *data, int *data_len) {
  struct pluq_params *p = (struct pluq_params *)_p;
  *data_len = 2;

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

  data[0] = walltime(0);
  data[1] = cpucycles();
  p->r = mzd_pluq(A, P, Q, 0);
  data[0] = walltime(data[0]);
  data[1]= cpucycles() - data[1];

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
  unsigned long long data[2];
  run_bench(run, (void*)&p, data, 2);

  printf("m: %5d, n: %5d, r: %5d, cpu cycles: %12llu, wall time: %6.3lf\n", p.m, p.n, p.r, data[1], data[0] / 1000000.0);
}

