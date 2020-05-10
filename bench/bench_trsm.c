#include <stdlib.h>

#include "benchmarking.h"
#include "cpucycles.h"
#include <m4ri/config.h>
#include <m4ri/m4ri.h>

struct trsm_params {
  rci_t m;
  rci_t n;
  int upper;
  int left;
  char const *algorithm;
};

mzd_t *mzd_random_lower(const rci_t n) {
  mzd_t *A = mzd_init(n, n);
  mzd_randomize(A);
  for (rci_t i = 0; i < n; i++) mzd_write_bit(A, i, i, 1);
  mzd_t *L = mzd_extract_l(NULL, A);
  mzd_free(A);
  return L;
}

mzd_t *mzd_random_upper(const rci_t n) {
  mzd_t *A = mzd_init(n, n);
  mzd_randomize(A);
  for (rci_t i = 0; i < n; i++) mzd_write_bit(A, i, i, 1);
  mzd_t *U = mzd_extract_u(NULL, A);
  mzd_free(A);
  return U;
}

int run(void *_p, unsigned long long *data, int *data_len) {
  struct trsm_params *p = (struct trsm_params *)_p;
  *data_len             = 2;

  mzd_t *T = NULL;
  mzd_t *B = mzd_init(p->m, p->n);
  mzd_randomize(B);

  if (p->upper) {
    T = mzd_random_upper(p->m);
  } else {
    T = mzd_random_lower(p->m);
  }

  data[0] = walltime(0);
  data[1] = cpucycles();
  switch (2 * p->upper + p->left) {
  case 3: mzd_trsm_upper_left(T, B, 0); break;
  case 2: mzd_trsm_upper_right(T, B, 0); break;
  case 1: mzd_trsm_lower_left(T, B, 0); break;
  case 0: mzd_trsm_lower_right(T, B, 0); break;
  default: m4ri_die("Parameters for upper (=%d) or left (=%d) not supported", p->upper, p->left);
  }
  data[0] = walltime(data[0]);
  data[1] = cpucycles() - data[1];

  mzd_free(B);
  mzd_free(T);
  return 0;
}

int main(int argc, char **argv) {
  int opts = global_options(&argc, &argv);

  if (opts < 0) {
    bench_print_global_options(stderr);
    exit(-1);
  }

  if (argc != 5) {
    printf("Parameters m,n,upper,left expected.\n");
    printf(" m         -- integer > 0\n");
    printf(" n         -- integer > 0\n");
    printf(" upper     -- 1 for upper triangular, 0 for lower triangular.\n");
    printf(" left      -- 1 for triangular matrix on left, 0 for right\n");
    printf("\n");
    bench_print_global_options(stderr);
    exit(-1);
  }

  struct trsm_params p;
  p.m     = atoi(argv[1]);
  p.n     = atoi(argv[2]);
  p.upper = atoi(argv[3]);
  p.left  = atoi(argv[4]);

  srandom(17);
  unsigned long long data[2];
  run_bench(run, (void *)&p, data, 2);

  /** this has no meaning if m << n **/
  double cc_per_op = (4 * (double)data[1]) / (p.m * powl((double)p.n, 1.807));

  printf("m: %5d, n: %5d, upper: %d, left: %d, cpu cycles: %llu, cc/(n^2.807): %.5lf, wall time: "
         "%lf\n",
         p.m, p.n, p.upper, p.left, data[1], cc_per_op, data[0] / 1000000.0);
}
