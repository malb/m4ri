#include <stdlib.h>

#include "m4ri/m4ri.h"
#include "cpucycles.h"
#include "walltime.h"

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
  packedmatrix *L = mzd_init(m, m);
  packedmatrix *U = mzd_init(m, n);
  packedmatrix *A = mzd_init(m, n);
  mzd_randomize(U);
  mzd_randomize(L);
  size_t i,j;
  for (i=0; i<m; ++i){
    for (j=i+1; j<m;++j)
      mzd_write_bit(L,i,j, 0);
    mzd_write_bit(L,i,i, 1);
    for (j=0; j<i;++j)
      mzd_write_bit(U,i,j, 0);
    mzd_write_bit(U,i,i, 1);
  }

  mzd_mul(A,L,U,0);

  permutation* P = mzp_init(m);
  permutation* Q = mzp_init(n);

  
  wt = walltime(&clockZero);
  t = cpucycles();
  mzd_lqup (A, P, Q,  2048);
  printf("n: %5d, cpu cycles: %12llu, wall time: %6.3lf\n",n, cpucycles() - t, walltime(&wt));

  mzp_free(P);
  mzp_free(Q);
  mzd_free(U);
  mzd_free(L);
}
