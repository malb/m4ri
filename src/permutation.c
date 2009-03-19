/******************************************************************************
*
*                 M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2008 Martin Albrecht <malb@informatik.uni-bremen.de> 
*
*  Distributed under the terms of the GNU General Public License (GPL)
*  version 2 or higher.
*
*    This code is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    General Public License for more details.
*
*  The full text of the GPL is available at:
*
*                  http://www.gnu.org/licenses/
******************************************************************************/

#include "permutation.h"
#include "packedmatrix.h"

mzp_t *mzp_init(size_t length) {
  size_t i;
  mzp_t *P = (mzp_t*)m4ri_mm_malloc(sizeof(mzp_t));
  P->values = (size_t*)m4ri_mm_malloc(sizeof(size_t)*length);
  P->length = length;
  for (i=0; i<length; i++) {
    P->values[i] = i;
  }
  return P;
}

void mzp_free(mzp_t *P) {
  m4ri_mm_free(P->values);
  m4ri_mm_free(P);
}

mzp_t *mzp_init_window(mzp_t* P, size_t begin, size_t end){
  mzp_t *window = (mzp_t *)m4ri_mm_malloc(sizeof(mzp_t));
  window->values = P->values + begin;
  window->length = end-begin;
  return window;
}

void mzp_free_window(mzp_t* condemned){
  m4ri_mm_free(condemned);
}

void mzp_set_ui(mzp_t *P, unsigned int value) {
  size_t i;
  for (i=0; i<P->length; i++) {
    P->values[i] = i;
  }
}

void mzd_apply_p_left(mzd_t *A, mzp_t *P) {
  size_t i;
  if(A->ncols == 0)
    return;
  for (i=0; i<P->length; i++) {
    assert(P->values[i] >= i);
    mzd_row_swap(A, i, P->values[i]);
  }
}

void mzd_apply_p_left_trans(mzd_t *A, mzp_t *P) {
  long i;
  if(A->ncols == 0)
    return;
  for (i=P->length-1; i>=0; i--) {
    assert(P->values[i] >= (size_t)i);
    mzd_row_swap(A, i, P->values[i]);
  }
}

void mzd_apply_p_right(mzd_t *A, mzp_t *P) {
  size_t i;
  if(A->nrows == 0)
    return;
  const size_t step_size = MAX((CPU_L1_CACHE>>3)/A->width,1);
  for(size_t j=0; j<A->nrows; j+=step_size) {
    size_t stop_row = MIN(j+step_size, A->nrows);
    for (i=0; i<P->length; ++i) {
      assert(P->values[i] >= i);
      mzd_col_swap_in_rows(A, i, P->values[i], j, stop_row);
    }
  }
/*   for (i=0; i<P->length; i++) { */
/*     assert(P->values[i] >= i); */
/*     mzd_col_swap(A, i, P->values[i]); */
/*   } */
}

void mzd_apply_p_right_trans(mzd_t *A, mzp_t *P) {
  int i;
  if(A->nrows == 0)
    return;
  const size_t step_size = MAX((CPU_L1_CACHE>>3)/A->width,1);
  for(size_t j=0; j<A->nrows; j+=step_size) {
    size_t stop_row = MIN(j+step_size, A->nrows);
    for (i=P->length-1; i>=0; --i) {
      assert(P->values[i] >= (size_t)i);
      mzd_col_swap_in_rows(A, i, P->values[i], j, stop_row);
    }
  }
/*   long i; */
/*   for (i=P->length-1; i>=0; i--) { */
/*     assert(P->values[i] >= i); */
/*     mzd_col_swap(A, i, P->values[i]); */
/*   } */
}

void mzd_col_block_rotate(mzd_t *M, size_t zs, size_t ze, size_t de, int copy) {
  size_t i,j;

  const size_t ds = ze;
  const size_t ld_f = (de - ze)/RADIX;
  const size_t ld_r = (de - ds)%RADIX;

  const size_t lz_f = (ze - zs)/RADIX;
  const size_t lz_r = (ze - zs)%RADIX;

  word *data = (word*)m4ri_mm_calloc(DIV_CEIL(de-ze, RADIX), sizeof(word));
  word *begin = (word*)m4ri_mm_calloc(DIV_CEIL(ze-zs, RADIX), sizeof(word));

  for(i=0; i<M->nrows; i++) {
    
    for(j=0; j < ld_f; j++) /* copy out */
      data[j] = mzd_read_bits(M, i, ds + j*RADIX, RADIX);
    if (ld_r)
      data[ld_f] = mzd_read_bits(M, i, ds + ld_f*RADIX, ld_r);

    for(j=0; j < lz_f; j++) /* copy out */
      begin[j] = mzd_read_bits(M, i, zs + j*RADIX, RADIX);
    if (lz_r)
      begin[lz_f] = mzd_read_bits(M, i, zs + lz_f*RADIX, lz_r);

    /* write */
    for(j=0; j<ld_f; j++) {
      mzd_clear_bits(M, i, zs + j*RADIX, RADIX);
      mzd_xor_bits(M, i, zs + j*RADIX, RADIX, data[j]);
    }
    if(ld_r) {
      mzd_clear_bits(M, i, zs + ld_f*RADIX, ld_r);
      mzd_xor_bits(M, i, zs + ld_f*RADIX, ld_r, data[ld_f]);
    }
    
    if (copy) {
      /* zero rest */
      for(j=0; j<lz_f; j++) {
        mzd_clear_bits(M, i, zs + (de - ds) + j*RADIX, RADIX);
        mzd_xor_bits(M, i, zs + (de - ds) + j*RADIX, RADIX, begin[j]);
      }
      if(lz_r) {
        mzd_clear_bits(M, i, zs + (de - ds) + lz_f*RADIX, lz_r);
        mzd_xor_bits(M, i, zs + (de - ds) + lz_f*RADIX, lz_r, begin[lz_f]);
      }
    }
  }
  
  m4ri_mm_free(data);
  m4ri_mm_free(begin);
}

void mzp_print(mzp_t *P) {
  printf("[ ");
  for(size_t i=0; i<P->length; i++) {
    printf("%zu ",P->values[i]);
  }
  printf("]");
}
