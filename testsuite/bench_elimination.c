#include <stdlib.h>

#include "cpucycles.h"
#include "walltime.h"
#include "m4ri.h"

int main(int argc, char **argv) {
  int n;
  unsigned long long t;
  double wt;
  double clockZero = 0.0;

  if (argc != 2) {
    m4ri_die("Parameter n expected.\n");
  }
  n = atoi(argv[1]);
  packedmatrix *A = mzd_init(n, n);
  mzd_randomize(A);
  
  wt = walltime(&clockZero);
  t = cpucycles();
  mzd_reduce_m4ri(A, 1, 0, NULL, NULL);
  printf("n: %5d, cpu cycles: %llu wall time: %lf\n",n, cpucycles() - t, walltime(&wt));

  mzd_free(A);
}
