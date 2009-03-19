#include <stdlib.h>

#include "cpucycles.h"
#include "walltime.h"
#include "m4ri/m4ri.h"

int main(int argc, char **argv) {
  int n,m;
  unsigned long long t;
  double wt;
  double clockZero = 0.0;

  if (argc != 3) {
    m4ri_die("Parameters m, n expected.\n");
  }
  m = atoi(argv[1]);
  n = atoi(argv[2]);
  mzd_t *B = mzd_init(m, n);
  mzd_t *L = mzd_init(m, m);
  mzd_randomize(B);
  mzd_randomize(L);
  size_t i,j;
  for (i=0; i<n; ++i){
    for (j=i+1; j<n;++j)
      mzd_write_bit(L,i,j, 0);
    mzd_write_bit(L,i,i, 1);
  }
  
  wt = walltime(&clockZero);
  t = cpucycles();
  mzd_trsm_lower_left(L, B, 2048);
  printf("n: %5d, cpu cycles: %llu wall time: %lf\n",n, cpucycles() - t, walltime(&wt));

  mzd_free(B);
  mzd_free(L);
}
