#include <stdlib.h>

#include "benchmarking.h"
#include "cpucycles.h"
#include <m4ri/config.h>
#include <m4ri/m4ri.h>

struct inv_params {
  rci_t n;
  int direction;
  char const *algorithm;
};

int run(void *_p, unsigned long long *data, int *data_len) {
  struct inv_params *p = (struct inv_params *)_p;
  *data_len            = 2;

  mzd_t *A = NULL, *L = NULL, *U = NULL, *B = NULL;

  if (p->direction <= 0) {
    L = mzd_init(p->n, p->n);
    mzd_randomize(L);
    for (rci_t i = 0; i < p->n; ++i) {
      for (rci_t j = i + 1; j < p->n; ++j) mzd_write_bit(L, i, j, 0);
      mzd_write_bit(L, i, i, 1);
    }
  }
  if (p->direction >= 0) {
    U = mzd_init(p->n, p->n);
    mzd_randomize(U);
    for (rci_t i = 0; i < p->n; ++i) {
      for (rci_t j = 0; j < i; ++j) mzd_write_bit(U, i, j, 0);
      mzd_write_bit(U, i, i, 1);
    }
  }
  switch (p->direction) {
  case 0:
    A = mzd_mul(NULL, L, U, 0);
    mzd_free(L);
    mzd_free(U);
    break;
  case -1: A = L; break;
  case 1: A = U; break;
  default: m4ri_die("unknown direction '%d'", p->direction);
  };

  data[0] = walltime(0);
  data[1] = cpucycles();
  switch (p->direction) {
  case 0:
    if (strcmp(p->algorithm, "m4ri") == 0)
      B = mzd_inv_m4ri(NULL, A, 0);
    else
      m4ri_die("unknown algorithm: '%s'", p->algorithm);
    break;
  case 1:
    if (strcmp(p->algorithm, "m4ri") == 0) {
      mzd_trtri_upper_russian(A, 0);
      B = mzd_copy(NULL, A);
    } else if (strcmp(p->algorithm, "mm") == 0) {
      mzd_trtri_upper(A);
      B = mzd_copy(NULL, A);
    } else
      m4ri_die("unknown algorithm: '%s'", p->algorithm);
    break;
  case -1: m4ri_die("not implemented error"); break;
  }
  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);
  mzd_free(A);
  mzd_free(B);
  return 0;
}

int main(int argc, char **argv) {
  int opts = global_options(&argc, &argv);

  if (opts < 0) {
    bench_print_global_options(stderr);
    exit(-1);
  }

  if (argc != 4) {
    printf("Parameters n,direction,alg expected.\n");
    printf(" n         -- integer > 0\n");
    printf(" direction -- lower triangular (-1), full (0) or upper triangular (1).\n");
    printf(" algorithm -- 'm4ri' or 'mm' (for direction 1)\n");
    printf("\n");
    bench_print_global_options(stderr);
    m4ri_die("");
  }

  struct inv_params params;
  params.n         = atoi(argv[1]);
  params.direction = atoi(argv[2]);
  params.algorithm = argv[3];

  /* put this call in run() to benchmark one particular matrix over
     and over again instead of computing the average of various
     matrices.*/
  srandom(17);
  unsigned long long data[2];
  run_bench(run, (void *)&params, data, 2);

  double cc_per_op = ((double)data[1]) / powl((double)params.n, 2.807);

  printf("n: %5d, cpu cycles: %10llu, cc/(n^2.807): %.5lf, wall time: %lf\n", params.n, data[1],
         cc_per_op, data[0] / 1000000.0);
}
