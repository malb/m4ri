#include <stdlib.h>

#include "m4ri/m4ri.h"
#include "cpucycles.h"
#include "walltime.h"

int main(int argc, char **argv) {
  size_t m, n;
  unsigned long long t;
  double wt;
  double clockZero = 0.0;

  if (argc != 3) {
    m4ri_die("Parameters m,n expected.\n");
  }
  m = atoi(argv[1]);
  n = atoi(argv[2]);
  packedmatrix *A = mzd_init(m, n);
  mzd_randomize(A);
  
  wt = walltime(&clockZero);
  t = cpucycles();
  size_t r = mzd_reduce_m4ri(A, 1, 0, NULL, NULL);
  printf("m: %5d, n: %5d, r: %5d, cpu cycles: %llu wall time: %lf\n",m, n, r, cpucycles() - t, walltime(&wt));

  mzd_free(A);
}
