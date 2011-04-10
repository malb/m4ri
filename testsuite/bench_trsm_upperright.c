#include <stdlib.h>

#include "config.h"
#include "cpucycles.h"
#include "m4ri.h"
#include "benchmarking.h"

struct trsm_params {
  rci_t m;
  rci_t n;
};

int run(void *_p, unsigned long long *data, int *data_len) {
  struct trsm_params *p = (struct trsm_params *)_p;
  *data_len = 2;

  mzd_t *B = mzd_init(p->m, p->n);
  mzd_t *U = mzd_init(p->n, p->n);
  mzd_randomize(B);
  mzd_randomize(U);
  for (rci_t i = 0; i < p->n; ++i){
    for (rci_t j = 0; j < i; ++j)
      mzd_write_bit(U,i,j, 0);
    mzd_write_bit(U,i,i, 1);
  }

  data[0] = walltime(0);
  data[1] = cpucycles();
  mzd_trsm_upper_right(U, B, 2048);
  data[0]= walltime(data[0]);
  data[1] = cpucycles() - data[1];

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
  unsigned long long data[2];
  run_bench(run, (void*)&p, data, 2);

  printf("m: %5d, n: %5d, cpu cycles: %llu wall time: %lf\n", p.m, p.n, data[1], data[0] / 1000000.0);
}
