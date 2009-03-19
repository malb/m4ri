#include <stdlib.h>

#include "m4ri/m4ri.h"
#include "cpucycles.h"
#include "walltime.h"

int main(int argc, char **argv) {
  int halfrank = 0;
  size_t m, n, r;
  const char *algorithm;
  unsigned long long t;
  double wt;
  double clockZero = 0.0;

  if (argc < 3) {
    m4ri_die("Parameters m,n expected.\n");
  }
  if (argc == 4)
    algorithm = argv[3];
  else
    algorithm = "m4ri";
  m = atoi(argv[1]);
  n = atoi(argv[2]);
  mzd_t *A = mzd_init(m, n);
  if (!halfrank) {
    mzd_randomize(A);
  } else {
    mzd_t *L, *U;
    L = mzd_init(m, m);
    U = mzd_init(m, n);
    mzd_randomize(U);
    mzd_randomize(L);
    size_t i,j;
    for (i=0; i<m; ++i){
      mzd_write_bit(U,i,i, 1);
      for (j=0; j<i;++j)
        mzd_write_bit(U,i,j, 0);
      if (i%2)
        for (j=i; j<n;++j)
          mzd_write_bit(U,i,j, 0);
      for (j=i+1; j<m;++j)
        mzd_write_bit(L,i,j, 0);
      mzd_write_bit(L,i,i, 1);
    }
    mzd_mul(A,L,U,0);
    mzd_free(L);
    mzd_free(U);
  }
  wt = walltime(&clockZero);
  t = cpucycles();
  if(strcmp(algorithm,"m4ri")==0)
    r = mzd_echelonize_m4ri(A, 1, 0);
  else if(strcmp(algorithm,"pluq")==0)
    r = mzd_echelonize_pluq(A, 1);
  else if(strcmp(algorithm,"naive")==0)
    r = mzd_echelonize_naive(A, 1);
  printf("m: %5d, n: %5d, r: %5d, cpu cycles: %llu wall time: %lf\n",m, n, r, cpucycles() - t, walltime(&wt));

  mzd_free(A);
}
