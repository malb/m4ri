#include <stdlib.h>

#include "config.h"
#include "cpucycles.h"
#include "walltime.h"
#include "m4ri.h"

int main(int argc, char **argv) {
  if (argc != 3) {
    m4ri_die("Parameters m, n expected.\n");
  }
  srandom(17);
  rci_t m = atoi(argv[1]);
  rci_t n = atoi(argv[2]);
  mzd_t *B = mzd_init(m, n);
  mzd_t *U = mzd_init(n, n);
  mzd_randomize(B);
  mzd_randomize(U);
  for (rci_t i = 0; i < n; ++i){
    for (rci_t j = 0; j < i; ++j)
      mzd_write_bit(U,i,j, 0);
    mzd_write_bit(U,i,i, 1);
  }

  double wt = walltime(0.0);
  unsigned long long t = cpucycles();
  mzd_trsm_upper_right(U, B, 2048);
  printf("n: %5d, cpu cycles: %llu wall time: %lf\n", n, cpucycles() - t, walltime(0.0));

  mzd_free(B);
  mzd_free(U);
}
