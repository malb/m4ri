#include "config.h"
#include <stdlib.h>
#include "m4ri.h"

int test_pluq_full_rank (rci_t m, rci_t n){
  printf("pluq: testing full rank m: %5zu, n: %5zu", m.val(), n.val());

  mzd_t* U = mzd_init (m,n);
  mzd_t* L = mzd_init (m,m);
  mzd_t* U2 = mzd_init (m,n);
  mzd_t* L2 = mzd_init (m,m);
  mzd_t* A = mzd_init (m,n);
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
  printf("pluq: testing half rank m: %5zd, n: %5zd", m.val(), n.val());

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

  printf("pluq: testing structured m: %5zd, n: %5zd", m.val(), n.val());

  mzd_t* A = mzd_init(m, n);
  mzd_t* L = mzd_init(m, m);
  mzd_t* U = mzd_init(m, n);

  for(rci_t i = 0; i < m; i += 2)
    for (rci_t j = i; j < n; ++j)
      mzd_write_bit(A, i, j, 1);

  mzd_t* Acopy = mzd_copy (NULL,A);

  mzp_t* P = mzp_init(m);
  mzp_t* Q = mzp_init(n);
  rci_t r = mzd_pluq(A, P, Q, 0);
  printf(", rank: %5d ",r.val());

  for (rci_t i = 0; i < r; ++i){
    for (rci_t j = 0; j < i; ++j)
      mzd_write_bit(L, i, j, mzd_read_bit(A,i,j));
    for (rci_t j = i + 1; j < n; ++j)
      mzd_write_bit(U, i, j, mzd_read_bit(A,i,j));
  }
  for (rci_t i = r; i < m; ++i)
    for (rci_t j = 0; j < r; ++j)
      mzd_write_bit(L, i, j, mzd_read_bit(A,i,j));
  for (rci_t i = 0; i < r; ++i){
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

  if (status) {
    printf("\n");
    printf(" ... FAILED\n");
  }  else
    printf (" ... passed\n");
  mzd_free(U);
  mzd_free(L);
  mzd_free(A);
  mzd_free(Acopy);
  mzp_free(P);
  mzp_free(Q);
  return status;
}

int test_pluq_random(rci_t m, rci_t n) {
  printf("pluq: testing random m: %5zd, n: %5zd", m.val(), n.val());

  mzd_t* U = mzd_init(m, n);
  mzd_t* L = mzd_init(m, m);
  mzd_t* A = mzd_init(m, n);
  mzd_randomize(A);

  mzd_t* Acopy = mzd_copy (NULL,A);

  mzp_t* P = mzp_init(m);
  mzp_t* Q = mzp_init(n);
  rci_t r = mzd_pluq(A, P, Q, 0);
  printf(", rank: %5d ", r.val());

  for (rci_t i = 0; i < r; ++i){
    for (rci_t j = 0; j < i; ++j)
      mzd_write_bit(L, i, j, mzd_read_bit(A,i,j));
    for (rci_t j = i + 1; j < n; ++j)
      mzd_write_bit(U, i, j, mzd_read_bit(A,i,j));
  }
  for (rci_t i = r; i < m; ++i)
    for (rci_t j = 0; j < r; ++j)
      mzd_write_bit(L, i, j, mzd_read_bit(A,i,j));
  for (rci_t i = 0; i < r; ++i){
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
  if (status) {
    printf(" ... FAILED\n");
  }  else
    printf (" ... passed\n");
  mzd_free(U);
  mzd_free(L);
  mzd_free(A);
  mzd_free(Acopy);
  mzp_free(P);
  mzp_free(Q);
  return status;
}


int main(int argc, char **argv) {
  int status = 0;

  status += test_pluq_structured(37U, 37U);
  status += test_pluq_structured(63U, 63U);
  status += test_pluq_structured(64U, 64U);
  status += test_pluq_structured(65U, 65U);
  status += test_pluq_structured(128U, 128U);

  status += test_pluq_structured(37U, 137U);
  status += test_pluq_structured(65U, 5U);
  status += test_pluq_structured(128U, 18U);

  status += test_pluq_full_rank(13U,13U);
  status += test_pluq_full_rank(37U,37U);
  status += test_pluq_full_rank(63U,63U);
  status += test_pluq_full_rank(64U,64U);
  status += test_pluq_full_rank(65U,65U);
  status += test_pluq_full_rank(97U,97U); 
  status += test_pluq_full_rank(128U,128U);
  status += test_pluq_full_rank(150U,150U);
  status += test_pluq_full_rank(256U,256U);
  status += test_pluq_full_rank(1024U,1024U);

  status += test_pluq_full_rank(13U,11U);
  status += test_pluq_full_rank(37U,39U);
  status += test_pluq_full_rank(64U,164U);
  status += test_pluq_full_rank(97U,92U);
  status += test_pluq_full_rank(128U,121U);
  status += test_pluq_full_rank(150U,153U);
  status += test_pluq_full_rank(256U,258U);
  status += test_pluq_full_rank(1024U,1023U);

  status += test_pluq_half_rank(64U,64U);
  status += test_pluq_half_rank(65U,65U);
  status += test_pluq_half_rank(66U,66U);
  status += test_pluq_half_rank(127U,127U);
  status += test_pluq_half_rank(129U,129U);
  status += test_pluq_half_rank(148U,148U);
  status += test_pluq_half_rank(132U,132U);
  status += test_pluq_half_rank(256U,256U);
  status += test_pluq_half_rank(1024U,1024U);

  status += test_pluq_half_rank(129U,127U);
  status += test_pluq_half_rank(132U,136U);
  status += test_pluq_half_rank(256U,251U);
  status += test_pluq_half_rank(1024U,2100U);

  status += test_pluq_random(63U,63U);
  status += test_pluq_random(64U,64U);
  status += test_pluq_random(65U,65U);

  status += test_pluq_random(128U,128U);
  status += test_pluq_random(128U, 131U);
  status += test_pluq_random(132U, 731U);
  status += test_pluq_random(150U,150U);
  status += test_pluq_random(252U, 24U);
  status += test_pluq_random(256U,256U);
  status += test_pluq_random(1024U,1022U);
  status += test_pluq_random(1024U,1024U);

  status += test_pluq_random(128U,1280U);
  status += test_pluq_random(128U, 130U);
  status += test_pluq_random(132U, 132U);
  status += test_pluq_random(150U,151U);
  status += test_pluq_random(252U, 2U);
  status += test_pluq_random(256U,251U);
  status += test_pluq_random(1024U,1025U);
  status += test_pluq_random(1024U,1021U);

  if (!status) {
    printf("All tests passed.\n");
    return 0;
  } else {
    return -1;
  }
}
