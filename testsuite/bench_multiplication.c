#include <stdlib.h>

#include "config.h"
#include "cpucycles.h"
#include "m4ri.h"
#include "benchmarketing.h"

struct mul_params {
  rci_t n;
  int cutoff;
};

int run(void *_p, double *wt, unsigned long long *cycles) {
  struct mul_params *p = (struct mul_params *)_p;

  mzd_t *A = mzd_init(p->n, p->n);
  mzd_t *B = mzd_init(p->n, p->n);
  mzd_randomize(A);
  mzd_randomize(B);

  *wt = walltime(0.0);
  *cycles = cpucycles();
  mzd_t *C = mzd_mul(NULL, A, B, p->cutoff);
  *wt = walltime(*wt);
  *cycles = cpucycles() - *cycles;

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
  unsigned long long t;
  double wt;
  run_bench(run, (void*)&p, &wt, &t);

  printf("n: %5d, cutoff: %5d, cpu cycles: %llu, wall time: %lf\n", p.n, p.cutoff, t, wt);
}
