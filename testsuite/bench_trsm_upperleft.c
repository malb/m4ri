#include <stdlib.h>

#include "config.h"
#include "cpucycles.h"
#include "m4ri.h"
#include "benchmarketing.h"

struct trsm_params {
  rci_t m;
  rci_t n;
};

int run(void *_p, double *wt, unsigned long long *cycles) {
  struct trsm_params *p = (struct trsm_params *)_p;

  mzd_t *B = mzd_init(p->m, p->n);
  mzd_t *U = mzd_init(p->m, p->m);
  mzd_randomize(B);
  mzd_randomize(U);
  for (rci_t i = 0; i < p->n; ++i){
    for (rci_t j = 0; j < i; ++j)
      mzd_write_bit(U,i,j, 0);
    mzd_write_bit(U,i,i, 1);
  }

  *wt = walltime(0.0);
  *cycles = cpucycles();
  mzd_trsm_upper_left(U, B, 0);
  *wt = walltime(*wt);
  *cycles = cpucycles() - *cycles;

  mzd_free(B);
  mzd_free(U);
  return 0;
}

int main(int argc, char **argv) {
  global_options(&argc, &argv);

  if (argc != 3) {
    m4ri_die("Parameters m, n expected.\n");
  }

  struct trsm_params p;
  p.m = atoi(argv[1]);
  p.n = atoi(argv[2]);
  
  srandom(17);
  unsigned long long t;
  double wt;
  run_bench(run, (void*)&p, &wt, &t);

  printf("m: %5d, n: %5d, cpu cycles: %llu wall time: %lf\n", p.m, p.n, t, wt);
}
