#include <stdlib.h>

#include "config.h"
#include "cpucycles.h"
#include "m4ri.h"
#include "benchmarketing.h"

struct smallops_params {
  rci_t n;
};

int run_transpose(void *_p, double *wt, unsigned long long *cycles) {
  struct smallops_params *p = (struct smallops_params *)_p;

  mzd_t *A = mzd_init(p->n, p->n);
  mzd_t *B = mzd_init(p->n, p->n);
  mzd_randomize(A);

  *wt = walltime(0.0);
  *cycles = cpucycles();
  for(int i = 0; i < 65536; i++) {
    mzd_transpose(A, B);
  }
  *wt = walltime(*wt)/65536.0;
  *cycles = (cpucycles() - *cycles)/65536;

  mzd_free(A);
  mzd_free(B);
  return (0);
}


int run_copy(void *_p, double *wt, unsigned long long *cycles) {
  struct smallops_params *p = (struct smallops_params *)_p;

  mzd_t *A = mzd_init(p->n, p->n);
  mzd_t *B = mzd_init(p->n, p->n);
  mzd_randomize(A);

  *wt = walltime(0.0);
  *cycles = cpucycles();
  for(int i = 0; i < 65536; i++) {
    mzd_copy(A, B);
  }
  *wt = walltime(*wt)/65536.0;
  *cycles = (cpucycles() - *cycles)/65536;

  mzd_free(A);
  mzd_free(B);
  return (0);
}


int run_add(void *_p, double *wt, unsigned long long *cycles) {
  struct smallops_params *p = (struct smallops_params *)_p;

  mzd_t *A = mzd_init(p->n, p->n);
  mzd_t *B = mzd_init(p->n, p->n);
  mzd_t *C = mzd_init(p->n, p->n);
  mzd_randomize(A);
  mzd_randomize(B);

  *wt = walltime(0.0);
  *cycles = cpucycles();
  for(int i = 0; i < 65536; i++) {
    _mzd_add(C, A, B);
  }
  *wt = walltime(*wt)/65536.0;
  *cycles = (cpucycles() - *cycles)/65536.0;

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
  unsigned long long t;
  double wt;
  if(strcmp(argv[2],"transpose") == 0) {

    run_bench(run_transpose, (void*)&p, &wt, &t);

  } else if(strcmp(argv[2],"add") == 0) {

    run_bench(run_add, (void*)&p, &wt, &t);

  } else if(strcmp(argv[2],"copy") == 0) {

    run_bench(run_copy, (void*)&p, &wt, &t);

  } else {

    m4ri_die("parameter not understood");
  }

  printf("what: %s, n: %5d, cpu cycles: %llu, cc/n^2: %lf, wall time: %lf\n", argv[2], p.n, t, t/((double)p.n)/((double)p.n), wt);
}
