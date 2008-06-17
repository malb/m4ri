#include <stdlib.h>
#include "cpucycles.h"
#include "m4ri.h"

int main(int argc, char **argv) {
  int n, cutoff;
  unsigned long long t;
  m4ri_build_all_codes();

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

  packedmatrix *A = mzd_init(n, n);
  packedmatrix *B = mzd_init(n, n);
  mzd_randomize(A);
  mzd_randomize(B);

  t = cpucycles();
  packedmatrix *C = mzd_mul(NULL, A, B, cutoff);
  printf("n: %5d, cutoff: %5d, cpu cycles: %llu\n",n, cutoff, cpucycles() - t);

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  m4ri_destroy_all_codes();
}
