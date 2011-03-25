#include <stdlib.h>

#include "config.h"
#include "m4ri.h"
#include "cpucycles.h"
#include "walltime.h"

int main(int argc, char **argv) {
  if (argc < 3) {
    m4ri_die("Parameters m,n, (alg,r) expected.\n");
  }

  char const *algorithm;
  long density = ~0;
  if (argc >= 4)
    algorithm = argv[3];
  else
    algorithm = "m4ri";
  if (argc >= 5)
    density = RAND_MAX * atof(argv[4]);

  int full = 1;
  if(argc >= 6)
    full = atoi(argv[5]);

  rci_t m = atoi(argv[1]);
  rci_t n = atoi(argv[2]);
  mzd_t *A = mzd_init(m, n);

  for(rci_t i = 0; i < m; ++i) {
    for(rci_t j = 0; j < n; ++j) {
      if(random() <= density) {
        mzd_write_bit(A, i, j, 1);
      }
    }
  }

  double clockZero = 0.0;
  double wt = walltime(&clockZero);
  unsigned long long t = cpucycles();
  rci_t r = 0;
  if(strcmp(algorithm,"m4ri") == 0)
    r = mzd_echelonize_m4ri(A, full, 0);
  else if(strcmp(algorithm,"cross") == 0)
    r = mzd_echelonize(A, full);
  else if(strcmp(algorithm,"pluq") == 0)
    r = mzd_echelonize_pluq(A, full);
  else if(strcmp(algorithm,"naive") == 0)
    r = mzd_echelonize_naive(A, full);
  else
    m4ri_die("Unknown algorithm (%s); should be one of: m4ri, cross, pluq, naive\n", algorithm);
  printf("m: %5d, n: %5d, r: %5d, cpu cycles: %10llu wall time: %lf\n", m, n, r, cpucycles() - t, walltime(&wt));

  mzd_free(A);
}
