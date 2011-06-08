#include <stdlib.h>

#include "config.h"
#include "cpucycles.h"
#include "m4ri.h"
#include "benchmarking.h"

struct mul_params {
  rci_t n;
  int cutoff;
};

int run(void *_p, unsigned long long *data, int *data_len) {
  struct mul_params *p = (struct mul_params *)_p;
  *data_len = 2;

  mzd_t *A = mzd_init(p->n, p->n);
  mzd_t *B = mzd_init(p->n, p->n);
  mzd_randomize(A);
  mzd_randomize(B);

  data[0] = walltime(0);
  data[1] = cpucycles();
  mzd_t *C = mzd_mul(NULL, A, B, p->cutoff);
  data[0] = walltime(data[0]);
  data[1] = cpucycles() - data[1];

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  return (0);
}

int main(int argc, char **argv) {
  global_options(&argc, &argv);

  if (argc != 3) {
    m4ri_die("Parameters n and cutoff expected.\n");
  }

  struct mul_params p;
  p.n = atoi(argv[1]);
  p.cutoff = atoi(argv[2]);

  if (p.n <= 0) {
    m4ri_die("Parameter n must be > 0\n");
  }

  /* put this call in run() to benchmark one particular matrix over
     and over again instead of computing the average of various
     matrices.*/
  srandom(17);
  unsigned long long data[2];
  run_bench(run, (void*)&p, data, 2);

  double cc_per_op = ((double)data[1])/ powl((double)p.n,2.807);

  printf("n: %5d, cutoff: %5d, cpu cycles: %llu, cc/n^2.807: %.5lf, wall time: %lf\n", p.n, p.cutoff, data[1], cc_per_op, data[0] / 1000000.0);
}
