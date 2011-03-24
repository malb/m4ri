#include <stdlib.h>

#include "config.h"
#include "m4ri.h"
#include "cpucycles.h"
#include "walltime.h"

int main(int argc, char **argv) {
  if (argc != 3) {
    m4ri_die("Parameters m, n expected.\n");
  }
  rci_t m = (unsigned int)atoi(argv[1]);
  rci_t n = (unsigned int)atoi(argv[2]);
  mzd_t *A = mzd_init(m, n);
  mzd_t *U;
  mzd_t *L;

  int halfrank = 1;
  if(halfrank) {
    U = mzd_init(m, n);
    L = mzd_init(m, m);
    mzd_randomize(U);
    mzd_randomize(L);
#if 0
    size_t i,j;
    for (i=0; i<m; ++i){
      for (j=i+1; j<m;++j)
        mzd_write_bit(L,i,j, 0);
      mzd_write_bit(L,i,i, 1);
    }
    for(i=0; i<MIN(m,n); ++i) {
      for (j=0; j<i;++j)
        mzd_write_bit(U,i,j, 0);
      mzd_write_bit(U,i,i, 1);
    }
#endif
    for(rci_t i = 0; i < m; ++i){
      mzd_write_bit(U,i,i, 1);
      for(rci_t j = 0; j < i; ++j)
	mzd_write_bit(U,i,j, 0);
      if ((i % 2))
	for(rci_t j = i; j < n; ++j)
	  mzd_write_bit(U,i,j, 0);
      for(rci_t j = i + 1; j < m;++j)
	mzd_write_bit(L,i,j, 0);
      mzd_write_bit(L,i,i, 1);
    }

    mzd_mul(A,L,U,0);

  } else {
    mzd_randomize(A);
  }

  mzp_t* P = mzp_init(m);
  mzp_t* Q = mzp_init(n);

  double clockZero = 0.0;
  double wt = walltime(&clockZero);
  unsigned long long t = cpucycles();
  rci_t r = mzd_pluq(A, P, Q, 0);
  printf("m: %5d, n: %5d, r: %5d, cpu cycles: %12llu, wall time: %6.3lf\n", m.val(), n.val(), r.val(), cpucycles() - t, walltime(&wt));

  mzd_free(A);
  mzp_free(P);
  mzp_free(Q);
  if(halfrank) {
    mzd_free(U);
    mzd_free(L);
  }
}
