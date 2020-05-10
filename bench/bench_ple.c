#include <stdlib.h>

#include "benchmarking.h"
#include "cpucycles.h"
#include <m4ri/config.h>
#include <m4ri/m4ri.h>

struct pluq_params {
  rci_t m;
  rci_t n;
  rci_t r;
  char *what;
};

int run(void *_p, unsigned long long *data, int *data_len) {
  struct pluq_params *p = (struct pluq_params *)_p;
  *data_len             = 2;

  mzd_t *A = mzd_init(p->m, p->n);
  mzd_t *U;
  mzd_t *L;

  int halfrank = 0;
  if (halfrank) {
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
    for (rci_t i = 0; i < p->m; ++i) {
      mzd_write_bit(U, i, i, 1);
      for (rci_t j = 0; j < i; ++j) mzd_write_bit(U, i, j, 0);
      if ((i % 2))
        for (rci_t j = i; j < p->n; ++j) mzd_write_bit(U, i, j, 0);
      for (rci_t j = i + 1; j < p->m; ++j) mzd_write_bit(L, i, j, 0);
      mzd_write_bit(L, i, i, 1);
    }

    mzd_mul(A, L, U, 0);

  } else {
    mzd_randomize(A);
  }

  mzp_t *P = mzp_init(p->m);
  mzp_t *Q = mzp_init(p->n);

  data[0] = walltime(0);
  data[1] = cpucycles();
  if (strcmp(p->what, "pluq"))
    p->r = mzd_pluq(A, P, Q, 0);
  else if (strcmp(p->what, "ple"))
    p->r = mzd_ple(A, P, Q, 0);
  else
    m4ri_die("Unknown task '%s'", p->what);
  data[0] = walltime(data[0]);
  data[1] = cpucycles() - data[1];

  mzd_free(A);
  mzp_free(P);
  mzp_free(Q);
  if (halfrank) {
    mzd_free(U);
    mzd_free(L);
  }

  return 0;
}

int main(int argc, char **argv) {
  int opts = global_options(&argc, &argv);

  if (opts < 0) {
    bench_print_global_options(stderr);
    exit(-1);
  }

  if (argc != 4) {
    printf("Parameters m,n,what expected.\n");
    printf(" m    -- integer > 0\n");
    printf(" n    -- integer > 0\n");
    printf(" what -- PLUQ or PLE.\n");
    printf("\n");
    bench_print_global_options(stderr);
    m4ri_die("");
  }

  struct pluq_params p;
  p.m    = atoi(argv[1]);
  p.n    = atoi(argv[2]);
  p.what = argv[3];

  srandom(17);
  unsigned long long data[2];
  run_bench(run, (void *)&p, data, 2);

  printf("m: %5d, n: %5d, what: %s, r: %5d, cpu cycles: %12llu, ", p.m, p.n, p.what, p.r, data[1]);
  print_wall_time(data[0] / 1000000.0);
  printf(", ");
  print_cpu_time(data[1] / (double)cpucycles_persecond());
  printf("\n");
}
