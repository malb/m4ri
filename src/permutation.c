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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "permutation.h"
#include "packedmatrix.h"

mzp_t *mzp_init(rci_t length) {
  mzp_t *P = (mzp_t*)m4ri_mm_malloc(sizeof(mzp_t));
  P->values = (rci_t*)m4ri_mm_malloc(sizeof(rci_t) * length.val());
  P->length = length;
  for (rci_t i = 0; i < length; ++i) {
    P->values[i.val()] = i;
  }
  return P;
}

void mzp_free(mzp_t *P) {
  m4ri_mm_free(P->values);
  m4ri_mm_free(P);
}

mzp_t *mzp_init_window(mzp_t *P, rci_t begin, rci_t end){
  mzp_t *window = (mzp_t *)m4ri_mm_malloc(sizeof(mzp_t));
  window->values = P->values + begin.val();
  window->length = end - begin;
  return window;
}

void mzp_free_window(mzp_t *condemned){
  m4ri_mm_free(condemned);
}

void mzp_set_ui(mzp_t *P, unsigned int value) {
  assert(value == 1);
  for (rci_t i = 0; i < P->length; ++i) {
    P->values[i.val()] = i;
  }
}

void mzd_apply_p_left(mzd_t *A, mzp_t *P) {
  if(A->ncols == 0)
    return;
  rci_t const length = MIN(P->length, A->nrows);
  for (rci_t i = 0; i < length; ++i) {
    assert(P->values[i.val()] >= i);
    mzd_row_swap(A, i, P->values[i.val()]);
  }
}

void mzd_apply_p_left_trans(mzd_t *A, mzp_t *P) {
  if(A->ncols == 0)
    return;
  rci_t const length = MIN(P->length, A->nrows);
  for (rci_t i = length - 1; i >= 0; --i) {
    assert(P->values[i.val()] >= i);
    mzd_row_swap(A, i, P->values[i.val()]);
  }
}

/* optimised column swap operations */

static inline void mzd_write_col_to_rows_blockd(mzd_t *A, mzd_t *B, rci_t *permutation, wordPtr write_mask, rci_t const start_row, rci_t const stop_row, rci_t length) {
  assert(A->offset == 0);
  for(rci_t i = 0; i < length; i += RADIX) {
    /* optimisation for identity permutations */
    if (write_mask[i / RADIX] == FFFF)
      continue;
    int const todo = MIN(RADIX.val(), (length - i).val());
    wi_t const a_word = (i + A->offset) / RADIX;
    wi_t words[RADIX.val()];
    int bits[RADIX.val()];
    word bitmasks[RADIX.val()];

    /* we pre-compute bit access in advance */
    for(int k = 0; k < todo; ++k) {
        rci_t const colb = permutation[(i + k).val()] + B->offset;
        words[k] = colb / RADIX;
        bits[k] = colb % RADIX;
        bitmasks[k] = ONE << bits[k];
    }

    for (rci_t r = start_row; r < stop_row; ++r) {
      wordPtr Brow = B->rows[r-start_row];
      wordPtr Arow = A->rows[r];
      register word value = 0;

      /* we gather the bits in a register word */
      switch(todo-1) {
      case 63: value |= ((Brow[words[63]] & bitmasks[63]) >> bits[63]) << 63;
      case 62: value |= ((Brow[words[62]] & bitmasks[62]) >> bits[62]) << 62;
      case 61: value |= ((Brow[words[61]] & bitmasks[61]) >> bits[61]) << 61;
      case 60: value |= ((Brow[words[60]] & bitmasks[60]) >> bits[60]) << 60;
      case 59: value |= ((Brow[words[59]] & bitmasks[59]) >> bits[59]) << 59;
      case 58: value |= ((Brow[words[58]] & bitmasks[58]) >> bits[58]) << 58;
      case 57: value |= ((Brow[words[57]] & bitmasks[57]) >> bits[57]) << 57;
      case 56: value |= ((Brow[words[56]] & bitmasks[56]) >> bits[56]) << 56;
      case 55: value |= ((Brow[words[55]] & bitmasks[55]) >> bits[55]) << 55;
      case 54: value |= ((Brow[words[54]] & bitmasks[54]) >> bits[54]) << 54;
      case 53: value |= ((Brow[words[53]] & bitmasks[53]) >> bits[53]) << 53;
      case 52: value |= ((Brow[words[52]] & bitmasks[52]) >> bits[52]) << 52;
      case 51: value |= ((Brow[words[51]] & bitmasks[51]) >> bits[51]) << 51;
      case 50: value |= ((Brow[words[50]] & bitmasks[50]) >> bits[50]) << 50;
      case 49: value |= ((Brow[words[49]] & bitmasks[49]) >> bits[49]) << 49;
      case 48: value |= ((Brow[words[48]] & bitmasks[48]) >> bits[48]) << 48;
      case 47: value |= ((Brow[words[47]] & bitmasks[47]) >> bits[47]) << 47;
      case 46: value |= ((Brow[words[46]] & bitmasks[46]) >> bits[46]) << 46;
      case 45: value |= ((Brow[words[45]] & bitmasks[45]) >> bits[45]) << 45;
      case 44: value |= ((Brow[words[44]] & bitmasks[44]) >> bits[44]) << 44;
      case 43: value |= ((Brow[words[43]] & bitmasks[43]) >> bits[43]) << 43;
      case 42: value |= ((Brow[words[42]] & bitmasks[42]) >> bits[42]) << 42;
      case 41: value |= ((Brow[words[41]] & bitmasks[41]) >> bits[41]) << 41;
      case 40: value |= ((Brow[words[40]] & bitmasks[40]) >> bits[40]) << 40;
      case 39: value |= ((Brow[words[39]] & bitmasks[39]) >> bits[39]) << 39;
      case 38: value |= ((Brow[words[38]] & bitmasks[38]) >> bits[38]) << 38;
      case 37: value |= ((Brow[words[37]] & bitmasks[37]) >> bits[37]) << 37;
      case 36: value |= ((Brow[words[36]] & bitmasks[36]) >> bits[36]) << 36;
      case 35: value |= ((Brow[words[35]] & bitmasks[35]) >> bits[35]) << 35;
      case 34: value |= ((Brow[words[34]] & bitmasks[34]) >> bits[34]) << 34;
      case 33: value |= ((Brow[words[33]] & bitmasks[33]) >> bits[33]) << 33;
      case 32: value |= ((Brow[words[32]] & bitmasks[32]) >> bits[32]) << 32;
      case 31: value |= ((Brow[words[31]] & bitmasks[31]) >> bits[31]) << 31;
      case 30: value |= ((Brow[words[30]] & bitmasks[30]) >> bits[30]) << 30;
      case 29: value |= ((Brow[words[29]] & bitmasks[29]) >> bits[29]) << 29;
      case 28: value |= ((Brow[words[28]] & bitmasks[28]) >> bits[28]) << 28;
      case 27: value |= ((Brow[words[27]] & bitmasks[27]) >> bits[27]) << 27;
      case 26: value |= ((Brow[words[26]] & bitmasks[26]) >> bits[26]) << 26;
      case 25: value |= ((Brow[words[25]] & bitmasks[25]) >> bits[25]) << 25;
      case 24: value |= ((Brow[words[24]] & bitmasks[24]) >> bits[24]) << 24;
      case 23: value |= ((Brow[words[23]] & bitmasks[23]) >> bits[23]) << 23;
      case 22: value |= ((Brow[words[22]] & bitmasks[22]) >> bits[22]) << 22;
      case 21: value |= ((Brow[words[21]] & bitmasks[21]) >> bits[21]) << 21;
      case 20: value |= ((Brow[words[20]] & bitmasks[20]) >> bits[20]) << 20;
      case 19: value |= ((Brow[words[19]] & bitmasks[19]) >> bits[19]) << 19;
      case 18: value |= ((Brow[words[18]] & bitmasks[18]) >> bits[18]) << 18;
      case 17: value |= ((Brow[words[17]] & bitmasks[17]) >> bits[17]) << 17;
      case 16: value |= ((Brow[words[16]] & bitmasks[16]) >> bits[16]) << 16;
      case 15: value |= ((Brow[words[15]] & bitmasks[15]) >> bits[15]) << 15;
      case 14: value |= ((Brow[words[14]] & bitmasks[14]) >> bits[14]) << 14;
      case 13: value |= ((Brow[words[13]] & bitmasks[13]) >> bits[13]) << 13;
      case 12: value |= ((Brow[words[12]] & bitmasks[12]) >> bits[12]) << 12;
      case 11: value |= ((Brow[words[11]] & bitmasks[11]) >> bits[11]) << 11;
      case 10: value |= ((Brow[words[10]] & bitmasks[10]) >> bits[10]) << 10;
      case  9: value |= ((Brow[words[ 9]] & bitmasks[ 9]) >> bits[ 9]) <<  9;
      case  8: value |= ((Brow[words[ 8]] & bitmasks[ 8]) >> bits[ 8]) <<  8;
      case  7: value |= ((Brow[words[ 7]] & bitmasks[ 7]) >> bits[ 7]) <<  7;
      case  6: value |= ((Brow[words[ 6]] & bitmasks[ 6]) >> bits[ 6]) <<  6;
      case  5: value |= ((Brow[words[ 5]] & bitmasks[ 5]) >> bits[ 5]) <<  5;
      case  4: value |= ((Brow[words[ 4]] & bitmasks[ 4]) >> bits[ 4]) <<  4;
      case  3: value |= ((Brow[words[ 3]] & bitmasks[ 3]) >> bits[ 3]) <<  3;
      case  2: value |= ((Brow[words[ 2]] & bitmasks[ 2]) >> bits[ 2]) <<  2;
      case  1: value |= ((Brow[words[ 1]] & bitmasks[ 1]) >> bits[ 1]) <<  1;
      case  0: value |= ((Brow[words[ 0]] & bitmasks[ 0]) >> bits[ 0]) <<  0;
      default:
        break;
      }
/*       for(int k = 0; k < todo; ++k) { */
/*         value |= ((Brow[words[k]] & bitmasks[k]) << bits[k]) >> k; */
/*       } */
      /* and write the word once */
      Arow[a_word] |= value;
    }
  }
}
/**
 * Implements both apply_p_right and apply_p_right_trans.
 */
void _mzd_apply_p_right_even(mzd_t *A, mzp_t *P, rci_t start_row, rci_t start_col, int notrans) {
  assert(A->offset == 0);
  if(A->nrows - start_row == 0)
    return;
  rci_t const length = MIN(P->length, A->ncols);
  wi_t const width = A->width;
  int step_size = MIN((A->nrows - start_row).val(), MAX((CPU_L1_CACHE >> 3) / A->width.val(), 1));

  /* our temporary where we store the columns we want to swap around */
  mzd_t *B = mzd_init((unsigned int)step_size, A->ncols);	// FIXME
  wordPtr Arow;
  wordPtr Brow;

  /* setup mathematical permutation */
  rci_t *permutation = (rci_t*)m4ri_mm_calloc(sizeof(rci_t), A->ncols.val());
  for(rci_t i = 0; i < A->ncols; ++i)
    permutation[i.val()] = i;

  if (!notrans) {
    for(rci_t i = start_col; i < length; ++i) {
      rci_t t = permutation[i.val()];
      permutation[i.val()] = permutation[P->values[i.val()].val()];
      permutation[P->values[i.val()].val()] = t;
    }
  } else {
    for(rci_t i = start_col; i < length; ++i) {
      rci_t t = permutation[(length-i-1).val()];
      permutation[(length-i-1).val()] = permutation[P->values[(length-i-1).val()].val()];
      permutation[P->values[(length-i-1).val()].val()] = t;
    }
  }

  /* we have a bitmask to encode where to write to */
  wordPtr write_mask = (word*)m4ri_mm_calloc(sizeof(word), width.val());
#ifdef M4RI_WRAPWORD
  word::init_array(write_mask, width);
#endif
  for(rci_t i = 0; i < A->ncols; i += RADIX) {
    int const todo = MIN(RADIX.val(), (A->ncols - i).val());
    for(int k = 0; k < todo; ++k) {
      if(permutation[(i+k).val()] == i+k) {
        write_mask[i/RADIX] |= ONE << k;
      }
    }
  }

  for(rci_t i = start_row; i < A->nrows; i += step_size) {
    step_size = MIN(step_size, (A->nrows - i).val());
    for(int k = 0; k < step_size; ++k) {
      Arow = A->rows[i+k];
      Brow = B->rows[rci_t(0) + k];

      /*copy row & clear those values which will be overwritten */
      for(wi_t j = 0U; j < width; ++j) {
        Brow[j] = Arow[j];
        Arow[j] = Arow[j] & write_mask[j];
      }
    }
    /* here we actually write out the permutation */
    mzd_write_col_to_rows_blockd(A, B, permutation, write_mask, i, i + step_size, length);
  }
  m4ri_mm_free(permutation);
  m4ri_mm_free(write_mask.val());
  mzd_free(B);
}

void _mzd_apply_p_right_trans(mzd_t *A, mzp_t *P) {
  if(A->nrows == 0)
    return;
  rci_t const length = MIN(P->length, A->ncols);
  int const step_size = MAX((CPU_L1_CACHE >> 3) / A->width.val(), 1);
  for(rci_t j = 0; j < A->nrows; j += step_size) {
    rci_t stop_row = MIN(j + step_size, A->nrows);
    for (rci_t i = 0; i < length; ++i) {
      assert(P->values[i.val()] >= i);
      mzd_col_swap_in_rows(A, i, P->values[i.val()], j, stop_row);
    }
  }
/*   for (i=0; i<P->length; i++) { */
/*     assert(P->values[i] >= i); */
/*     mzd_col_swap(A, i, P->values[i]); */
/*   } */
}

void _mzd_apply_p_right(mzd_t *A, mzp_t *P) {
  if(A->nrows == 0)
    return;
  int const step_size = MAX((CPU_L1_CACHE >> 3) / A->width.val(), 1);
  for(rci_t j = 0; j < A->nrows; j += step_size) {
    rci_t stop_row = MIN(j + step_size, A->nrows);
    for (rci_t i = P->length - 1; i >= 0; --i) {
      assert(P->values[i.val()] >= i);
      mzd_col_swap_in_rows(A, i, P->values[i.val()], j, stop_row);
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
  if(A->offset) {
    _mzd_apply_p_right_trans(A,P);
    return;
  }
  _mzd_apply_p_right_even(A, P, 0, 0, 0); 
}

void mzd_apply_p_right(mzd_t *A, mzp_t *P) {
  if(!A->nrows)
    return;
  if(A->offset) {
    _mzd_apply_p_right(A,P);
    return;
  }
  _mzd_apply_p_right_even(A, P, 0, 0, 1); 
}

void mzd_apply_p_right_trans_even_capped(mzd_t *A, mzp_t *P, rci_t start_row, rci_t start_col) {
  assert(A->offset == 0);
  if(!A->nrows)
    return;
  _mzd_apply_p_right_even(A, P, start_row, start_col, 0); 
}

void mzd_apply_p_right_even_capped(mzd_t *A, mzp_t *P, rci_t start_row, rci_t start_col) {
  assert(A->offset == 0);
  if(!A->nrows)
    return;
  _mzd_apply_p_right_even(A, P, start_row, start_col, 1); 
}

void mzp_print(mzp_t *P) {
  printf("[ ");
  for(rci_t i = 0; i < P->length; ++i) {
    printf("%zu ", P->values[i.val()].val());
  }
  printf("]");
}

void mzd_apply_p_right_trans_tri(mzd_t *A, mzp_t *P) {
  assert(P->length == A->ncols);
  int const step_size = MAX((CPU_L1_CACHE >> 2) / A->width.val(), 1);

  for(rci_t r = 0; r < A->nrows; r += step_size) {
    rci_t const row_bound = MIN(r + step_size, A->nrows);
    for (rci_t i =0 ; i < A->ncols; ++i) {
      assert(P->values[i.val()] >= i);
      mzd_col_swap_in_rows(A, i, P->values[i.val()], r, MIN(row_bound, i));
    }
  }
}

void _mzd_compress_l(mzd_t *A, rci_t r1, rci_t n1, rci_t r2) {
  /**
   * We are compressing this matrix
\verbatim
           r1           n1
   ------------------------------------------
   | \ \____|___        | A01               |
   |  \     |   \       |                   |
 r1------------------------------------------ 
   |   |    |           | \  \_____         |
   | L1|    |           |  \       \________|
   |   |    |           | L2|               |
   ------------------------------------------
\endverbatim
  *
  * to this matrix
  *
\verbatim
           r1           n1
   ------------------------------------------
   | \ \____|___        | A01               |
   |  \     |   \       |                   |
 r1------------------------------------------ 
   |    \   |           |    \_____         |
   |     \  |           |          \________|
   |      | |           |                   |
   ------------------------------------------
\endverbatim
  */

  if (r1 == n1)
    return;

#if 0

  mzp_t *shift = mzp_init(A->ncols);
  for (rci_t i=r1,j=n1;i<r1+r2;i++,j++){
    mzd_col_swap_in_rows(A, i, j, i, r1+r2);
    shift->values[i] = j;
  }

  mzd_apply_p_right_trans_even_capped(A, shift, r1+r2, 0);
  mzp_free(shift);

#else

  for (rci_t i = r1, j = n1; i < r1 + r2; ++i, ++j){
    mzd_col_swap_in_rows(A, i, j, i, r1 + r2);
  }
  
  word tmp;
  wi_t block;

  for(rci_t i = r1 + r2; i < A->nrows; ++i) {

    rci_t j = r1;

    /* first we deal with the rest of the current word we need to
       write */
    int const rest = RADIX - ((j + A->offset) % RADIX);

    tmp = mzd_read_bits(A, i, n1, rest);
    mzd_clear_bits(A, i, j, rest);
    mzd_xor_bits(A, i, j, rest, tmp);
    
    j += rest;

    /* now each write is simply a word write */

    block = (n1+j-r1+A->offset) / RADIX;

    if (rest%RADIX == 0) {
      for( ;j+RADIX<=r1+r2; j+=RADIX, block++) {
        tmp = A->rows[i][block];
        A->rows[i][j/RADIX] = tmp;
      }
    } else {
      for(; j+RADIX<=r1+r2; j+=RADIX, block++) {
        tmp = (A->rows[i][block] >> rest) | ( A->rows[i][block + 1U] << (RADIX - rest)); 
        A->rows[i][j/RADIX] = tmp;
      }
    }

    /* we deal with the remaining bits. While we could write past the
       end of r1+r2 here, but we have no guarantee that we can read
       past the end of n1+r2. */

    if (j<r1+r2) {
      tmp = mzd_read_bits(A, i, n1+j-r1, (r1+r2-j).val());
      A->rows[i][j/RADIX] = tmp;
    }

    /* now clear the rest of L2 */
    j = r1+r2;
    mzd_clear_bits(A, i, j, RADIX - ((j+A->offset)%RADIX));

    j += RADIX - ((j+A->offset)%RADIX);

    /* it's okay to write the full word, i.e. past n1+r2, because
       everything is zero there anyway. Thus, we can omit the code
       which deals with last few bits. */

    for(; j<n1+r2; j+=RADIX) {
      A->rows[i][j/RADIX] = 0;
    }
  }
 
#endif
}

