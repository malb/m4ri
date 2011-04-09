#include <stdlib.h>

#include "config.h"
#include "cpucycles.h"
#include "m4ri.h"
#include "benchmarking.h"

struct smallops_params {
  rci_t n;
};

int run_transpose(void *_p, unsigned long long *data, int *data_len) {
  struct smallops_params *p = (struct smallops_params *)_p;
  *data_len = 2;

  mzd_t *A = mzd_init(p->n, p->n);
  mzd_t *B = mzd_init(p->n, p->n);
  mzd_randomize(A);

  data[0] = walltime(0);
  data[1] = cpucycles();
  for(int i = 0; i < 65536; i++) {
    mzd_transpose(A, B);
  }
  data[0] = walltime(data[0]);
  data[1] = cpucycles() - data[1];

  mzd_free(A);
  mzd_free(B);
  return (0);
}


int run_copy(void *_p, unsigned long long *data, int *data_len) {
  struct smallops_params *p = (struct smallops_params *)_p;
  *data_len = 2;

  mzd_t *A = mzd_init(p->n, p->n);
  mzd_t *B = mzd_init(p->n, p->n);
  mzd_randomize(A);

  data[0] = walltime(0);
  data[1] = cpucycles();
  for(int i = 0; i < 65536; i++) {
    mzd_copy(A, B);
  }
  data[0] = walltime(data[0]);
  data[1] = cpucycles() - data[1];

  mzd_free(A);
  mzd_free(B);
  return (0);
}


int run_add(void *_p, unsigned long long *data, int *data_len) {
  struct smallops_params *p = (struct smallops_params *)_p;
  *data_len = 2;

  mzd_t *A = mzd_init(p->n, p->n);
  mzd_t *B = mzd_init(p->n, p->n);
  mzd_t *C = mzd_init(p->n, p->n);
  mzd_randomize(A);
  mzd_randomize(B);

  data[0] = walltime(0);
  data[1] = cpucycles();
  for(int i = 0; i < 65536; i++) {
    _mzd_add(C, A, B);
  }
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
    m4ri_die("Parameter n and one of 'add','copy' or 'transpose' expected.\n");
    
  }

  struct smallops_params p;
  p.n = atoi(argv[1]);

  if (p.n <= 0) {
    m4ri_die("Parameter n must be > 0\n");
  }

  /* put this call in run() to benchmark one particular matrix over
     and over again instead of computing the average of various
     matrices.*/
  srandom(17);
  bench_count = 65536;
  unsigned long long data[2];
  if(strcmp(argv[2],"transpose") == 0) {

    run_bench(run_transpose, (void*)&p, data, 2);

  } else if(strcmp(argv[2],"add") == 0) {

    run_bench(run_add, (void*)&p, data, 2);

  } else if(strcmp(argv[2],"copy") == 0) {

    run_bench(run_copy, (void*)&p, data, 2);

  } else {

    m4ri_die("parameter not understood");
  }

  printf("what: %s, n: %5d, cpu cycles: %llu, cc/n^2: %lf, wall time: %lf\n", argv[2], p.n, data[1], data[1]/((double)p.n)/((double)p.n), data[0] / 1000000.0);
}
