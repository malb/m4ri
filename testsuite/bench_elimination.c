#include <stdlib.h>
#include "packedmatrix.h"
#include "brilliantrussian.h"
#include "grayflex.h"
#include "cpucycles.h"

int main(int argc, char **argv) {
  int n;
  unsigned long long t;
  packedmatrix *A;

  //initialise the library
  m4ri_build_all_codes();
  if (argc != 2) {
    m4ri_die("Parameter n expected.\n");
  }
  n = atoi(argv[1]);
  A = mzd_init(n, n);
  mzd_randomize(A);
  
  t = cpucycles();
  // standard parameter
  mzd_reduce_m4ri(A, 1, 0, NULL, NULL);
  printf("n: %5d, cpu cycles: %llu\n",n, cpucycles() - t);

  mzd_free(A);
  m4ri_destroy_all_codes();
}
