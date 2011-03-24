#include <stdlib.h>

#include "config.h"
#include "cpucycles.h"
#include "walltime.h"
#include "m4ri.h"

int main(int argc, char **argv) {
  if (argc != 3) {
    m4ri_die("Parameters m, n expected.\n");
  }
  rci_t m = (unsigned int)atoi(argv[1]);
  rci_t n = (unsigned int)atoi(argv[2]);
  mzd_t *B = mzd_init(m, n);
  mzd_t *L = mzd_init(n, n);
  mzd_randomize(B);
  mzd_randomize(L);
  for (rci_t i = 0; i < n; ++i){
    for (rci_t j = i + 1; j < n; ++j)
      mzd_write_bit(L,i,j, 0);
    mzd_write_bit(L,i,i, 1);
  }
  
  double clockZero = 0.0;
  double wt = walltime(&clockZero);
  unsigned long long t = cpucycles();
  mzd_trsm_lower_right(L, B, 2048);
  printf("n: %5d, cpu cycles: %llu wall time: %lf\n", n.val(), cpucycles() - t, walltime(&wt));

  mzd_free(B);
  mzd_free(L);
}
