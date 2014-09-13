#include <m4ri/config.h>
#include <stdlib.h>
#include <m4ri/m4ri.h>

int check_pluq(mzd_t *A, rci_t *r) {
  mzd_t* Acopy = mzd_copy (NULL,A);

  const rci_t m = A->nrows;
  const rci_t n = A->ncols;
  
  mzd_t* L = mzd_init(m, m);
  mzd_t* U = mzd_init(m, n);
  
  mzp_t* P = mzp_init(m);
  mzp_t* Q = mzp_init(n);
  r[0] = mzd_pluq(A, P, Q, 0);

  for (rci_t i = 0; i < r[0]; ++i){
    for (rci_t j = 0; j < i; ++j)
      mzd_write_bit(L, i, j, mzd_read_bit(A,i,j));
    for (rci_t j = i + 1; j < n; ++j)
      mzd_write_bit(U, i, j, mzd_read_bit(A,i,j));
  }
  for (rci_t i = r[0]; i < m; ++i)
    for (rci_t j = 0; j < r[0]; ++j)
      mzd_write_bit(L, i, j, mzd_read_bit(A,i,j));
  for (rci_t i = 0; i < r[0]; ++i){
    mzd_write_bit(L,i,i, 1);
    mzd_write_bit(U,i,i, 1);
  }

  mzd_apply_p_left(Acopy, P);
  mzd_apply_p_right_trans(Acopy, Q);

  mzd_addmul(Acopy, L, U, 0);

  int status = 0;
  for (rci_t i = 0; i < m; ++i)
    for (rci_t j = 0; j < n; ++j){
      if (mzd_read_bit (Acopy,i,j)){
	status = 1;
        break;
      }
    }  
  mzd_free(U);
  mzd_free(L);
  mzd_free(Acopy);
  mzp_free(P);
  mzp_free(Q);
  
  return status;
}

int test_pluq_full_rank (rci_t m, rci_t n){
  printf("pluq: testing full rank m: %5d, n: %5d", m, n);

  mzd_t* U  = mzd_init(m,n);
  mzd_t* L  = mzd_init(m,m);
  mzd_t* U2 = mzd_init(m,n);
  mzd_t* L2 = mzd_init(m,m);
  mzd_t* A = mzd_init(m,n);
  mzd_randomize (U);
  mzd_randomize (L);

  for (rci_t i = 0; i < m; ++i){
    for (rci_t j = 0; j < i && j < n;++j)
      mzd_write_bit(U,i,j, 0);
    for (rci_t j = i + 1; j < m; ++j)
      mzd_write_bit(L,i,j, 0);
    if(i<n)
      mzd_write_bit(U,i,i, 1);
    mzd_write_bit(L,i,i, 1);
  }
  
  mzd_mul(A, L, U, 2048);

  mzd_t* Acopy = mzd_copy (NULL,A);

  mzp_t* P = mzp_init(m);
  mzp_t* Q = mzp_init(n);
  mzd_pluq(A, P, Q, 2048);

  for (rci_t i = 0; i < m; ++i){
    for (rci_t j = 0; j < i && j < n; ++j)
      mzd_write_bit (L2, i, j, mzd_read_bit(A,i,j));
    for (rci_t j = i + 1; j < n; ++j)
      mzd_write_bit (U2, i, j, mzd_read_bit(A,i,j));
  }
  
  for (rci_t i = 0; i < n && i < m; ++i){
    mzd_write_bit(L2,i,i, 1);
    mzd_write_bit(U2,i,i, 1);
  }
  mzd_addmul(Acopy,L2,U2,0);
  int status = 0;
  for (rci_t i = 0; i < m; ++i)
    for (rci_t j=0; j < n; ++j){
      if (mzd_read_bit (Acopy,i,j)){
	status = 1;
      }
    }
  if (status){
    printf(" ... FAILED\n");
  }  else
    printf (" ... passed\n");
  mzd_free(U);
  mzd_free(L);
  mzd_free(U2);
  mzd_free(L2);
  mzd_free(A);
  mzd_free(Acopy);
  mzp_free(P);
  mzp_free(Q);
  return status;
}

int test_pluq_half_rank(rci_t m, rci_t n) {
  printf("pluq: testing half rank m: %5d, n: %5d", m, n);

  mzd_t* U = mzd_init(m, n);
  mzd_t* L = mzd_init(m, m);
  mzd_t* U2 = mzd_init(m, n);
  mzd_t* L2 = mzd_init(m, m);
  mzd_t* A = mzd_init(m, n);
  mzd_randomize (U);
  mzd_randomize (L);

  for (rci_t i = 0; i < m && i < n; ++i) {
    mzd_write_bit(U,i,i, 1);
    for (rci_t j = 0; j < i;++j)
      mzd_write_bit(U,i,j, 0);
    if (i%2)
      for (rci_t j = i; j < n;++j)
	mzd_write_bit(U,i,j, 0);
    for (rci_t j = i + 1; j < m; ++j)
      mzd_write_bit(L,i,j, 0);
    mzd_write_bit(L,i,i, 1);
  }
  
  mzd_mul(A, L, U, 0);

  mzd_t* Acopy = mzd_copy (NULL,A);

  mzp_t* Pt = mzp_init(m);
  mzp_t* Q = mzp_init(n);
  rci_t r = mzd_pluq(A, Pt, Q, 0);

  for (rci_t i = 0; i < r; ++i) {
    for (rci_t j = 0; j < i; ++j)
      mzd_write_bit (L2, i, j, mzd_read_bit(A,i,j));
    for (rci_t j = i + 1; j < n; ++j)
      mzd_write_bit (U2, i, j, mzd_read_bit(A,i,j));
  }
  for (rci_t i = r; i < m; ++i)
    for (rci_t j = 0; j < r;++j)
      mzd_write_bit (L2, i, j, mzd_read_bit(A,i,j));
  for (rci_t i = 0; i < r; ++i){
    mzd_write_bit(L2,i,i, 1);
    mzd_write_bit(U2,i,i, 1);
  }

  mzd_apply_p_left(Acopy, Pt);
  mzd_apply_p_right_trans(Acopy, Q);

  mzd_addmul(Acopy,L2,U2,0);

  int status = 0;
  for (rci_t i = 0; i < m; ++i) {
    for (rci_t j = 0; j < n; ++j){
      if (mzd_read_bit(Acopy,i,j)){
	status = 1;
      }
    }
    if(status)
      break;
  }
  if (status)
    printf(" ... FAILED\n");
  else
    printf (" ... passed\n");
  mzd_free(U);
  mzd_free(L);
  mzd_free(U2);
  mzd_free(L2);
  mzd_free(A);
  mzd_free(Acopy);
  mzp_free(Pt);
  mzp_free(Q);
  return status;
}

int test_pluq_structured(rci_t m, rci_t n) {

  printf("pluq: testing structured m: %5d, n: %5d", m, n);

  mzd_t* A = mzd_init(m, n);
  mzd_t* L = mzd_init(m, m);
  mzd_t* U = mzd_init(m, n);

  for(rci_t i = 0; i < m; i += 2)
    for (rci_t j = i; j < n; ++j)
      mzd_write_bit(A, i, j, 1);

  rci_t r = 0;
  int status = check_pluq(A, &r);
  printf(", rank: %5d ",r);
  
  if (status) {
    printf(" ... FAILED\n");
  }  else
    printf (" ... passed\n");
  mzd_free(A);
  return status;
}

int test_pluq_random(rci_t m, rci_t n) {
  printf("pluq: testing random m: %5d, n: %5d", m, n);

  mzd_t* A = mzd_init(m, n);
  mzd_randomize(A);

  rci_t r = 0;
  int status = check_pluq(A, &r);
  printf(", rank: %5d ",r);
  
  if (status) {
    printf(" ... FAILED\n");
  }  else
    printf (" ... passed\n");
  mzd_free(A);
  return status;
}

int test_pluq_string(rci_t m, rci_t n, const char *str) {
  printf("pluq: testing string m: %5d, n: %5d", m, n);
  
  
  mzd_t *A = mzd_from_str(m, n, str);

  mzd_t *Acopy = mzd_copy(NULL, A);
  mzp_t *P = mzp_init(A->nrows);
  mzp_t *Q = mzp_init(A->ncols);
  _mzd_ple_russian(Acopy, P, Q, 0);
  mzd_print(Acopy);
  rci_t r = 0;
  int status = check_pluq(A, &r);
  printf(", rank: %5d ",r);
  
  if (status) {
    printf(" ... FAILED\n");
  }  else
    printf (" ... passed\n");
  mzd_free(A);
  return status;
}


int main() {
  int status = 0;

  srandom(17);

  status += test_pluq_string(4, 4, "0101011100010110");
  
  status += test_pluq_structured(37, 37);
  status += test_pluq_structured(63, 63);
  status += test_pluq_structured(64, 64);
  status += test_pluq_structured(65, 65);
  status += test_pluq_structured(128, 128);

  status += test_pluq_structured(37, 137);
  status += test_pluq_structured(65, 5);
  status += test_pluq_structured(128, 18);

  status += test_pluq_full_rank(13, 13);
  status += test_pluq_full_rank(37, 37);
  status += test_pluq_full_rank(63, 63);
  status += test_pluq_full_rank(64, 64);
  status += test_pluq_full_rank(65, 65);
  status += test_pluq_full_rank(97, 97); 
  status += test_pluq_full_rank(128, 128);
  status += test_pluq_full_rank(150, 150);
  status += test_pluq_full_rank(256, 256);
  status += test_pluq_full_rank(1024, 1024);

  status += test_pluq_full_rank(13, 11);
  status += test_pluq_full_rank(37, 39);
  status += test_pluq_full_rank(64, 164);
  status += test_pluq_full_rank(97, 92);
  status += test_pluq_full_rank(128, 121);
  status += test_pluq_full_rank(150, 153);
  status += test_pluq_full_rank(256, 258);
  status += test_pluq_full_rank(1024, 1023);

  status += test_pluq_half_rank(64, 64);
  status += test_pluq_half_rank(65, 65);
  status += test_pluq_half_rank(66, 66);
  status += test_pluq_half_rank(127, 127);
  status += test_pluq_half_rank(129, 129);
  status += test_pluq_half_rank(148, 148);
  status += test_pluq_half_rank(132, 132);
  status += test_pluq_half_rank(256, 256);
  status += test_pluq_half_rank(1024, 1024);

  status += test_pluq_half_rank(129, 127);
  status += test_pluq_half_rank(132, 136);
  status += test_pluq_half_rank(256, 251);
  status += test_pluq_half_rank(1024, 2100);

  status += test_pluq_random(63, 63);
  status += test_pluq_random(64, 64);
  status += test_pluq_random(65, 65);

  status += test_pluq_random(128, 128);
  status += test_pluq_random(128, 131);
  status += test_pluq_random(132, 731);
  status += test_pluq_random(150, 150);
  status += test_pluq_random(252, 24);
  status += test_pluq_random(256, 256);
  status += test_pluq_random(1024, 1022);
  status += test_pluq_random(1024, 1024);

  status += test_pluq_random(128, 1280);
  status += test_pluq_random(128, 130);
  status += test_pluq_random(132, 132);
  status += test_pluq_random(150, 151);
  status += test_pluq_random(252, 2);
  status += test_pluq_random(256, 251);
  status += test_pluq_random(1024, 1025);
  status += test_pluq_random(1024, 1021);

  if (!status) {
    printf("All tests passed.\n");
    return 0;
  } else {
    return -1;
  }
}
