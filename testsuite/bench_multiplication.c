#include <stdlib.h>

#include "config.h"
#include "cpucycles.h"
#include "walltime.h"
#include "m4ri.h"

int main(int argc, char **argv) {
  if (argc != 3) {
    m4ri_die("Parameters n and cutoff expected.\n");
  }
  rci_t n = (unsigned int)atoi(argv[1]);
  if (n <= 0) {
    m4ri_die("Parameter n must be > 0\n");
  }
  int cutoff = atoi(argv[2]);
  if (cutoff <= 0) {
    m4ri_die("Parameter cutoff must be > 0\n");
  }

  mzd_t *A = mzd_init(n, n);
  mzd_t *B = mzd_init(n, n);
  mzd_randomize(A);
  mzd_randomize(B);

  double clockZero = 0.0;
  double wt = walltime(&clockZero);
  unsigned long long t = cpucycles();
  mzd_t *C = mzd_mul(NULL, A, B, cutoff);
  printf("n: %5d, cutoff: %5d, cpu cycles: %llu wall time: %lf\n", n.val(), cutoff, cpucycles() - t, walltime(&wt));

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
}
