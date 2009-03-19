#include <stdlib.h>

#include "cpucycles.h"
#include "walltime.h"
#include "m4ri/m4ri.h"

int main(int argc, char **argv) {
  int n, cutoff;
  unsigned long long t;
  double wt;
  double clockZero = 0.0;

  if (argc != 3) {
    m4ri_die("Parameters n and cutoff expected.\n");
  }
  n = atoi(argv[1]);
  cutoff = atoi(argv[2]);

  if (n<=0) {
    m4ri_die("Parameter n must be > 0\n");
  }

  if (cutoff<=0) {
    m4ri_die("Parameter cutoff must be > 0\n");
  }

  mzd_t *A = mzd_init(n, n);
  mzd_t *B = mzd_init(n, n);
  mzd_randomize(A);
  mzd_randomize(B);

  wt = walltime(&clockZero);
  t = cpucycles();
  mzd_t *C = mzd_mul(NULL, A, B, cutoff);
  printf("n: %5d, cutoff: %5d, cpu cycles: %llu wall time: %lf\n",n, cutoff, cpucycles() - t, walltime(&wt));

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
}
