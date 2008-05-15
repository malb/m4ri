#include <stdlib.h>
#include "cpucycles.h"
#include "m4ri.h"

int main(int argc, char **argv) {
  int n;
  unsigned long t;
  m4ri_build_all_codes();

  if (argc != 2) {
    m4ri_die("Parameters n expected.\n");
  }
  n = atoi(argv[1]);

  if (n<=0) {
    m4ri_die("Parameter n must be > 0\n");
  }

  packedmatrix *A = mzd_init(n, n);
  packedmatrix *B = mzd_init(n, n);
  mzd_randomize(A);
  mzd_randomize(B);

  t = cpucycles();
  packedmatrix *C = mzd_add(A, B, A);
  printf("n: %5d, cpu cycles: %10u\n",n, cpucycles() - t);

  mzd_free(A);
  mzd_free(B);
  //mzd_free(C);
  m4ri_destroy_all_codes();
}
