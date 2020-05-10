#include <stdlib.h>

#include "benchmarking.h"
#include "cpucycles.h"
#include <m4ri/config.h>
#include <m4ri/m4ri.h>

struct elim_sparse_params {
  rci_t m;
  rci_t n;
  rci_t r;
  char const *algorithm;

  long density;
  int full;
};

int run(void *_p, unsigned long long *data, int *data_len) {
  struct elim_sparse_params *p = (struct elim_sparse_params *)_p;
  *data_len                    = 2;

  mzd_t *A = mzd_init(p->m, p->n);
  for (rci_t i = 0; i < p->m; ++i) {
    for (rci_t j = 0; j < p->n; ++j) {
      if (random() <= p->density) { mzd_write_bit(A, i, j, 1); }
    }
  }

  data[0] = walltime(0);
  data[1] = cpucycles();
  if (strcmp(p->algorithm, "m4ri") == 0)
    p->r = mzd_echelonize_m4ri(A, p->full, 0);
  else if (strcmp(p->algorithm, "cross") == 0)
    p->r = mzd_echelonize(A, p->full);
  else if (strcmp(p->algorithm, "pluq") == 0)
    p->r = mzd_echelonize_pluq(A, p->full);
  else if (strcmp(p->algorithm, "naive") == 0)
    p->r = mzd_echelonize_naive(A, p->full);
  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);
  mzd_free(A);
  return 0;
}

int main(int argc, char **argv) {
  global_options(&argc, &argv);

  if (argc < 3) { m4ri_die("Parameters m,n, (alg,density,full) expected.\n"); }

  struct elim_sparse_params p;
  p.density = RAND_MAX / 10;  // Use a density of 0.1 by default.
  p.full    = 1;

  if (argc >= 4)
    p.algorithm = argv[3];
  else
    p.algorithm = "m4ri";
  if (argc >= 5) p.density = RAND_MAX * strtod(argv[4], NULL);
  if (argc >= 6) p.full = atoi(argv[5]);

  p.m = atoi(argv[1]);
  p.n = atoi(argv[2]);

  /* put this call in run() to benchmark one particular matrix over
     and over again instead of computing the average of various
     matrices.*/
  srandom(17);
  unsigned long long data[2];
  run_bench(run, (void *)&p, data, 2);

  printf("m: %5d, n: %5d, last r: %5d, density: %7.5f, cpu cycles: %10llu, wall time: %lf\n", p.m,
         p.n, p.r, (double)p.density / RAND_MAX, data[1], data[0] / 1000000.0);
}
