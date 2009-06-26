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
  const size_t length = MIN(P->length, A->nrows);
  for (i=0; i<length; i++) {
    assert(P->values[i] >= i);
    mzd_row_swap(A, i, P->values[i]);
  }
}

void mzd_apply_p_left_trans(mzd_t *A, mzp_t *P) {
  long i;
  if(A->ncols == 0)
    return;
  const size_t length = MIN(P->length, A->nrows);
  for (i=length-1; i>=0; i--) {
    assert(P->values[i] >= (size_t)i);
    mzd_row_swap(A, i, P->values[i]);
  }
}

/* optimised column swap operations */

static inline void mzd_write_col_to_rows_blockd(mzd_t *A, mzd_t *B, size_t *permutation, word *write_mask, const size_t start_row, const size_t stop_row, size_t length) {
  assert(A->offset == 0);
  for(size_t i=0; i<length; i+=RADIX) {
    /* optimisation for identity permutations */
    if (write_mask[i/RADIX] == FFFF)
      continue;
    const size_t todo = MIN(RADIX,length-i);
    const size_t a_word = (A->offset+i)/RADIX;
    size_t words[RADIX];
    size_t bits[RADIX];
    size_t bitmasks[RADIX];

    /* we pre-compute bit access in advance */
    for(size_t k=0;k<todo; k++) {
        const size_t colb = permutation[i+k] + B->offset;
        words[k] = colb/RADIX;
        bits[k] = colb%RADIX;
        bitmasks[k] = (ONE<<(RADIX - (bits[k]) - 1));
    }

    for (size_t r=start_row; r<stop_row; r++) {
      word *Brow = B->rows[r-start_row];
      word *Arow = A->rows[r];
      register word value = 0;
      /* we gather the bits in a register word */
      for(register size_t k=0; k<todo; k++) {
        value |= ((Brow[words[k]] & bitmasks[k]) << bits[k]) >> k;
      }
      /* and write the word once */
      Arow[a_word] |= value;
    }
  }
}
/**
 * Implements both apply_p_right and apply_p_right_trans.
 */
void _mzd_apply_p_right_even(mzd_t *A, mzp_t *P, int notrans) {
  assert(A->offset = 0);
  const size_t length = MIN(P->length,A->ncols);
  const size_t width = A->width;
  size_t step_size = MIN(A->nrows, MAX((CPU_L1_CACHE>>3)/A->width,1));

  /* our temporary where we store the columns we want to swap around */
  mzd_t *B = mzd_init(step_size, A->ncols);
  word *Arow;
  word *Brow;

  /* setup mathematical permutation */
  size_t *permutation = m4ri_mm_calloc(sizeof(size_t),A->ncols);
  for(size_t i=0; i<A->ncols; i++)
    permutation[i] = i;

  if (!notrans) {
    for(size_t i=0; i<length; i++) {
      size_t t = permutation[i];
      permutation[i] = permutation[P->values[i]];
      permutation[P->values[i]] = t;
    }
  } else {
    for(size_t i=0; i<length; i++) {
      size_t t = permutation[length-i-1];
      permutation[length-i-1] = permutation[P->values[length-i-1]];
      permutation[P->values[length-i-1]] = t;
    }
  }

  /* we have a bitmask to encode where to write to */
  word *write_mask = m4ri_mm_calloc(sizeof(size_t), length);
  for(size_t i=0; i<A->ncols; i+=RADIX) {
    const size_t todo = MIN(RADIX,A->ncols-i);
    for(size_t k=0; k<todo; k++) {
      if(permutation[i+k] == i+k) {
        write_mask[i/RADIX] |= ONE<<(RADIX - k - 1);
      }
    }
  }

  for(size_t i=0; i<A->nrows; i+=step_size) {
    step_size = MIN(step_size, A->nrows-i);

    for(size_t k=0; k<step_size; k++) {
      Arow = A->rows[i+k];
      Brow = B->rows[k];

      /*copy row & clear those values which will be overwritten */
      for(size_t j=0; j<width; j++) {
        Brow[j] = Arow[j];
        Arow[j] = Arow[j] & write_mask[j];
      }
    }
    /* here we actually write out the permutation */
    mzd_write_col_to_rows_blockd(A, B, permutation, write_mask, i, i+step_size, length);
  }
  m4ri_mm_free(permutation);
  m4ri_mm_free(write_mask);
  mzd_free(B);
}

void _mzd_apply_p_right_trans(mzd_t *A, mzp_t *P) {
  size_t i;
  if(A->nrows == 0)
    return;
  const size_t length = MIN(P->length, A->ncols);
  const size_t step_size = MAX((CPU_L1_CACHE>>3)/A->width,1);
  for(size_t j=0; j<A->nrows; j+=step_size) {
    size_t stop_row = MIN(j+step_size, A->nrows);
    for (i=0; i<length; ++i) {
      assert(P->values[i] >= i);
      mzd_col_swap_in_rows(A, i, P->values[i], j, stop_row);
    }
  }
/*   for (i=0; i<P->length; i++) { */
/*     assert(P->values[i] >= i); */
/*     mzd_col_swap(A, i, P->values[i]); */
/*   } */
}

void _mzd_apply_p_right(mzd_t *A, mzp_t *P) {
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


void mzd_apply_p_right_trans(mzd_t *A, mzp_t *P) {
  if(!A->nrows)
    return;
  if(1/*A->offset*/) {
    _mzd_apply_p_right_trans(A,P);
    return;
  }
  _mzd_apply_p_right_even(A, P, 0); 
}

void mzd_apply_p_right(mzd_t *A, mzp_t *P) {
  if(!A->nrows)
    return;
  if(1/*A->offset*/) {
    _mzd_apply_p_right(A,P);
    return;
  }
  _mzd_apply_p_right_even(A, P, 1); 
}


void mzd_col_block_rotate(mzd_t *M, size_t zs, size_t ze, size_t de) {
  size_t i,j;
  const size_t ds = ze;
/*   const size_t ld_f = (de - ds)/RADIX; */
/*   const size_t ld_r = (de - ds)%RADIX; */
  
  const size_t lz_f = (ze - zs)/RADIX;
  const size_t lz_r = (ze - zs)%RADIX;

  const size_t le_f = (M->ncols - de)/RADIX;
  const size_t le_r = (M->ncols - de)%RADIX;

  size_t n1 = ze;
  size_t r1 = zs;
  size_t r2 = de - ds;
  
  for(i=0; i<M->nrows; i++) {
    
/*     for(j=0; j < i; j++) /\* copy out *\/ */
/*       data->rows[0][j] = M->rows[i][j]; */

    /* write */
    size_t im = (i+1<r2)?i+1:r2;
    size_t ld_f = im / RADIX;
    size_t ld_r = im % RADIX;
    for(j=0; j<ld_f; j++) {
      mzd_clear_bits(M, i, zs + j*RADIX, RADIX);
      mzd_xor_bits(M, i, zs + j*RADIX, RADIX, mzd_read_bits(M,i,ds+j*RADIX,RADIX));
    }
    if(ld_r) {
      mzd_clear_bits(M, i, zs + ld_f*RADIX, ld_r);
      mzd_xor_bits(M, i, zs + ld_f*RADIX, ld_r, mzd_read_bits(M,i,ds+ld_f*RADIX,ld_r));
    }
    //mzd_write_bit(M,i,i+r1,1);

    /* Placing zeros */
   for (j = r1+im; j<n1+im; ++j)
     mzd_write_bit(M,i,j,0);
  }
  // mzd_free(data);
}

void mzp_print(mzp_t *P) {
  printf("[ ");
  for(size_t i=0; i<P->length; i++) {
    printf("%zu ",P->values[i]);
  }
  printf("]");
}

void  mzd_apply_p_right_tri (mzd_t * A, mzp_t * P){
  /* To be optimized */
  size_t i;
  assert(P->length==A->ncols);
  for (i =0 ; i<P->length; ++i){
      assert(P->values[i] >= i);
      if (P->values[i] > i){
	mzd_col_swap_in_rows(A, i, P->values[i], 0,  i);
      }
  }
}
