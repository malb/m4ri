#include <stdlib.h>

#include "m4ri/m4ri.h"
#include "cpucycles.h"
#include "walltime.h"

int main(int argc, char **argv) {
  int halfrank = 1;
  size_t m, n;
  unsigned long long t;
  double wt;
  double clockZero = 0.0;

  if (argc != 3) {
    m4ri_die("Parameters m,n expected.\n");
  }
  m = atoi(argv[1]);
  n = atoi(argv[2]);
  packedmatrix *A = mzd_init(m, n);
  if (!halfrank) {
    mzd_randomize(A);
  } else {
    packedmatrix *L, *U;
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
  size_t r = mzd_reduce_m4ri(A, 1, 0, NULL, NULL);
  printf("m: %5d, n: %5d, r: %5d, cpu cycles: %llu wall time: %lf\n",m, n, r, cpucycles() - t, walltime(&wt));

  mzd_free(A);
}
