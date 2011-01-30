#include <stdlib.h>

#include "m4ri.h"
#include "cpucycles.h"
#include "walltime.h"

int main(int argc, char **argv) {
  size_t m, n, r;
  long density = ~0;
  const char *algorithm;
  unsigned long long t;
  double wt;
  double clockZero = 0.0;

  if (argc < 3) {
    m4ri_die("Parameters m,n, (alg,r) expected.\n");
  }
  if (argc >= 4)
    algorithm = argv[3];
  else
    algorithm = "m4ri";
  if (argc == 5)
    density = RAND_MAX * atof(argv[4]);

  m = atoi(argv[1]);
  n = atoi(argv[2]);
  mzd_t *A = mzd_init(m, n);


  for(size_t i=0; i<m; i++) {
    for(size_t j=0; j<n; j++) {
      if(random() <= density) {
        mzd_write_bit(A, i, j, 1);
      }
    }
  }

  wt = walltime(&clockZero);
  t = cpucycles();
  if(strcmp(algorithm,"m4ri")==0)
    r = mzd_echelonize_m4ri(A, 1, 0);
  else if(strcmp(algorithm,"pluq")==0)
    r = mzd_echelonize_pluq(A, 1);
  else if(strcmp(algorithm,"mmpf")==0)
    r = _mzd_pluq_mmpf(A, mzp_init(A->nrows),mzp_init(A->ncols),0);
  else if(strcmp(algorithm,"naive")==0)
    r = mzd_echelonize_naive(A, 1);
  printf("m: %5d, n: %5d, r: %5d, cpu cycles: %10llu wall time: %lf\n",m, n, r, cpucycles() - t, walltime(&wt));

  mzd_free(A);
}
