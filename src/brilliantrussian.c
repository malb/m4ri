/*******************************************************************
*
*                 M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2007, 2008 Gregory Bard <bard@fordham.edu>
*    Copyright (C) 2008-2010 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
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
*
********************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>

#include "brilliantrussian.h"
#include "xor.h"
#include "grayflex.h"
#include "echelonform.h"

/**
 * \brief Perform Gaussian reduction to reduced row echelon form on a
 * submatrix.
 * 
 * The submatrix has dimension at most k starting at r x c of A. Checks
 * for pivot rows up to row endrow (exclusive). Terminates as soon as
 * finding a pivot column fails.
 *
 * \param A Matrix.
 * \param r First row.
 * \param c First column.
 * \param k Maximal dimension of identity matrix to produce.
 * \param end_row Maximal row index (exclusive) for rows to consider
 * for inclusion.
 */

static inline int _mzd_gauss_submatrix_full(mzd_t *A, rci_t r, rci_t c, rci_t end_row, int k) {
  assert(k <= RADIX);
  rci_t start_row = r;
  rci_t j;
  for (j = c; j < c + k; ++j) {
    int found = 0;
    for (rci_t i = start_row; i < end_row; ++i) {
      /* first we need to clear the first columns */
      word const tmp = mzd_read_bits(A, i, c, j - c + 1);
      if(tmp) {
        for (int l = 0; l < j - c; ++l)
          if (GET_BIT(tmp, l))
            mzd_row_add_offset(A, i, r+l, c+l);
      
        /* pivot? */
        if (mzd_read_bit(A, i, j)) {
          mzd_row_swap(A, i, start_row);
          /* clear above */
          for (rci_t l = r; l < start_row; ++l) {
            if (mzd_read_bit(A, l, j)) {
              mzd_row_add_offset(A, l, start_row, j);
            }
          }
          ++start_row;
          found = 1;
          break;
        }
      }
    }
    if (found == 0) {
      return j - c;
    }
  }
  return j - c;
}

/**
 * \brief Perform Gaussian reduction to upper triangular matrix on a
 * submatrix.
 * 
 * The submatrix has dimension at most k starting at r x c of A. Checks
 * for pivot rows up to row end_row (exclusive). Terminates as soon as
 * finding a pivot column fails.
 *
 * \param A Matrix.
 * \param r First row.
 * \param c First column.
 * \param k Maximal dimension of identity matrix to produce.
 * \param end_row Maximal row index (exclusive) for rows to consider
 * for inclusion.
 */

static inline int _mzd_gauss_submatrix(mzd_t *A, rci_t r, rci_t c, rci_t end_row, int k) {
  rci_t start_row = r;
  int found;
  rci_t j;
  for (j = c; j < c+k; ++j) {
    found = 0;
    for (rci_t i = start_row; i < end_row; ++i) {
      /* first we need to clear the first columns */
      for (int l = 0; l < j - c; ++l)
        if (mzd_read_bit(A, i, c+l))
          mzd_row_add_offset(A, i, r+l, c+l);
      
      /* pivot? */
      if (mzd_read_bit(A, i, j)) {
        mzd_row_swap(A, i, start_row);
        start_row++;
        found = 1;
        break;
      }
    }
    if (found == 0) {
      return j - c;
    }
  }
  return j - c;
}

/**
 * \brief Given a submatrix in upper triangular form compute the
 * reduced row echelon form.
 * 
 * The submatrix has dimension at most k starting at r x c of A. Checks
 * for pivot rows up to row end_row (exclusive). Terminates as soon as
 * finding a pivot column fails.
 *
 * \param A Matrix.
 * \param r First row.
 * \param c First column.
 * \param k Maximal dimension of identity matrix to produce.
 * \param end_row Maximal row index (exclusive) for rows to consider
 * for inclusion.
 */

static inline int _mzd_gauss_submatrix_top(mzd_t *A, rci_t r, rci_t c, int k) {
  rci_t start_row = r;
  for (rci_t j = c; j < c + k; ++j) {
    for (rci_t l = r; l < start_row; ++l) {
      if (mzd_read_bit(A, l, j)) {
        mzd_row_add_offset(A, l, start_row, j);
      }
    }
    ++start_row;
  }
  return k;
}

static inline void _mzd_copy_back_rows(mzd_t *A, mzd_t *U, rci_t r, rci_t c, int k) {
  wi_t const startblock = c / RADIX;
  wi_t const width = A->width - startblock;
  for (int i = 0; i < k; ++i) {
    word const *const src = U->rows[i] + startblock;
    word *const dst = A->rows[r+i] + startblock;
    for (wi_t j = 0; j < width; ++j) {
      dst[j] = src[j];
    }
  }
}

void mzd_make_table(mzd_t *M, rci_t r, rci_t c, int k, mzd_t *T, rci_t *L)
{
  wi_t const homeblock = (c + M->offset) / RADIX;
  word const mask_end = LEFT_BITMASK((M->ncols + M->offset) % RADIX);
  word const pure_mask_begin = RIGHT_BITMASK(RADIX - ((c + M->offset) % RADIX));
  word const mask_begin = (M->width - homeblock != 1) ? pure_mask_begin : pure_mask_begin & mask_end;
  wi_t const wide = M->width - homeblock;

  int const twokay = TWOPOW(k);
  L[0] = 0;
  for (rci_t i = 1; i < twokay; ++i) {
    word *ti = T->rows[i] + homeblock;
    word *ti1 = T->rows[i-1] + homeblock;   

    rci_t const rowneeded = r + codebook[k]->inc[i - 1];
    int const id = codebook[k]->ord[i];
    L[id] = i;

    if (rowneeded >= M->nrows)
      continue;

    word *m = M->rows[rowneeded] + homeblock;

    *ti++ = (*m++ ^ *ti1++) & mask_begin;
    
    wi_t j;
    for(j = 1; j + 8 <= wide - 1; j += 8) {
      *ti++ = *m++ ^ *ti1++;
      *ti++ = *m++ ^ *ti1++;
      *ti++ = *m++ ^ *ti1++;
      *ti++ = *m++ ^ *ti1++;
      *ti++ = *m++ ^ *ti1++;
      *ti++ = *m++ ^ *ti1++;
      *ti++ = *m++ ^ *ti1++;
      *ti++ = *m++ ^ *ti1++;
    }
    switch(wide - j) {
    case 8:  *ti++ = *m++ ^ *ti1++;
    case 7:  *ti++ = *m++ ^ *ti1++;
    case 6:  *ti++ = *m++ ^ *ti1++;
    case 5:  *ti++ = *m++ ^ *ti1++;
    case 4:  *ti++ = *m++ ^ *ti1++;
    case 3:  *ti++ = *m++ ^ *ti1++;
    case 2:  *ti++ = *m++ ^ *ti1++;
    case 1:  *ti++ = (*m++ ^ *ti1++) & mask_end;
    }
  }
}

void mzd_process_rows(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k, mzd_t *T, rci_t *L) {
  wi_t const block = startcol / RADIX;
  wi_t const wide = M->width - block;
  wi_t const count = (wide + 7) / 8;	/* Unrolled loop count */
  int const entry_point = wide % 8;	/* Unrolled loop entry point */

  if(k == 1) {
    word const bm = ONE << (startcol % RADIX);

    rci_t r;
    for (r = startrow; r + 2 <= stoprow; r += 2) {
      word const b0 = M->rows[r+0][block] & bm;
      word const b1 = M->rows[r+1][block] & bm;

      word *m0 = M->rows[r+0] + block;
      word *m1 = M->rows[r+1] + block;
      word *t = T->rows[1] + block;

      wi_t n = count;
      if((b0 & b1)) {
	switch (entry_point) {
	case 0: do { *m0++ ^= *t; *m1++ ^= *t++;
	  case 7:    *m0++ ^= *t; *m1++ ^= *t++;
	  case 6:    *m0++ ^= *t; *m1++ ^= *t++;
	  case 5:    *m0++ ^= *t; *m1++ ^= *t++;
	  case 4:    *m0++ ^= *t; *m1++ ^= *t++;
	  case 3:    *m0++ ^= *t; *m1++ ^= *t++;
	  case 2:    *m0++ ^= *t; *m1++ ^= *t++;
	  case 1:    *m0++ ^= *t; *m1++ ^= *t++;
	  } while (--n > 0);
	}
      } else if(b0) {
	switch (entry_point) {
	case 0: do { *m0++ ^= *t++;
	  case 7:    *m0++ ^= *t++;
	  case 6:    *m0++ ^= *t++;
	  case 5:    *m0++ ^= *t++;
	  case 4:    *m0++ ^= *t++;
	  case 3:    *m0++ ^= *t++;
	  case 2:    *m0++ ^= *t++;
	  case 1:    *m0++ ^= *t++;
	  } while (--n > 0);
	}
      } else if(b1) {
	switch (entry_point) {
	case 0: do { *m1++ ^= *t++;
	  case 7:    *m1++ ^= *t++;
	  case 6:    *m1++ ^= *t++;
	  case 5:    *m1++ ^= *t++;
	  case 4:    *m1++ ^= *t++;
	  case 3:    *m1++ ^= *t++;
	  case 2:    *m1++ ^= *t++;
	  case 1:    *m1++ ^= *t++;
	  } while (--n > 0);
	}
      }
    }

    /* TODO: this code is a bit silly/overkill, it just takes care of the last row */
    for( ; r < stoprow; ++r) {
      rci_t const x0 = L[ mzd_read_bits_int(M, r, startcol, k) ];

      word *m0 = M->rows[r] + block;
      word *t0 = T->rows[x0] + block;

      wi_t n = count;
      switch (entry_point) {
      case 0: do { *m0++ ^= *t0++;
        case 7:    *m0++ ^= *t0++;
        case 6:    *m0++ ^= *t0++;
        case 5:    *m0++ ^= *t0++;
        case 4:    *m0++ ^= *t0++;
        case 3:    *m0++ ^= *t0++;
        case 2:    *m0++ ^= *t0++;
        case 1:    *m0++ ^= *t0++;
        } while (--n > 0);
      }
    }
    return;
  }

  rci_t r;
  for (r = startrow; r + 2 <= stoprow; r += 2) {
    rci_t const x0 = L[ mzd_read_bits_int(M, r+0, startcol, k) ];
    rci_t const x1 = L[ mzd_read_bits_int(M, r+1, startcol, k) ];
    
    word *m0 = M->rows[r+0] + block;
    word *t0 = T->rows[x0] + block;

    word *m1 = M->rows[r+1] + block;
    word *t1 = T->rows[x1] + block;

    wi_t n = count;
    switch (entry_point) {
    case 0: do { *m0++ ^= *t0++; *m1++ ^= *t1++;
      case 7:    *m0++ ^= *t0++; *m1++ ^= *t1++;
      case 6:    *m0++ ^= *t0++; *m1++ ^= *t1++;
      case 5:    *m0++ ^= *t0++; *m1++ ^= *t1++;
      case 4:    *m0++ ^= *t0++; *m1++ ^= *t1++;
      case 3:    *m0++ ^= *t0++; *m1++ ^= *t1++;
      case 2:    *m0++ ^= *t0++; *m1++ ^= *t1++;
      case 1:    *m0++ ^= *t0++; *m1++ ^= *t1++;
      } while (--n > 0);
    }
  }

  for( ; r < stoprow; ++r) {
    rci_t const x0 = L[ mzd_read_bits_int(M, r, startcol, k) ];

    word *m0 = M->rows[r] + block;
    word *t0 = T->rows[x0] + block;

    wi_t n = count;
    switch (entry_point) {
    case 0: do { *m0++ ^= *t0++;
      case 7:    *m0++ ^= *t0++;
      case 6:    *m0++ ^= *t0++;
      case 5:    *m0++ ^= *t0++;
      case 4:    *m0++ ^= *t0++;
      case 3:    *m0++ ^= *t0++;
      case 2:    *m0++ ^= *t0++;
      case 1:    *m0++ ^= *t0++;
      } while (--n > 0);
    }
  }
}

void mzd_process_rows2(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k, mzd_t *T0, rci_t *L0, mzd_t *T1, rci_t *L1) {
  wi_t const blocknum = startcol / RADIX;
  wi_t const wide = M->width - blocknum;
  wi_t const count = (wide + 7) / 8;	/* Unrolled loop count */
  int const entry_point = wide % 8;	/* Unrolled loop entry point */

  int const ka = k / 2;
  int const kb = k - k / 2;

#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(static,512) // MAX((CPU_L1_CACHE >> 3) / wide,
#endif
  for(rci_t r = startrow; r < stoprow; ++r) {
    rci_t const x0 = L0[ mzd_read_bits_int(M, r, startcol, ka)];
    rci_t const x1 = L1[ mzd_read_bits_int(M, r, startcol+ka, kb)];
    if((x0 | x1) == 0)	// x0 == 0 && x1 == 0
      continue;
    word *m0 = M->rows[r] + blocknum;
    word const *t0 = T0->rows[x0] + blocknum;
    word const *t1 = T1->rows[x1] + blocknum;

    wi_t n = count;
    switch (entry_point) {
    case 0: do { *m0++ ^= *t0++ ^ *t1++;
      case 7:    *m0++ ^= *t0++ ^ *t1++;
      case 6:    *m0++ ^= *t0++ ^ *t1++;
      case 5:    *m0++ ^= *t0++ ^ *t1++;
      case 4:    *m0++ ^= *t0++ ^ *t1++;
      case 3:    *m0++ ^= *t0++ ^ *t1++;
      case 2:    *m0++ ^= *t0++ ^ *t1++;
      case 1:    *m0++ ^= *t0++ ^ *t1++;
      } while (--n > 0);
    }
  }
}

void mzd_process_rows3(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k, mzd_t *T0, rci_t *L0, mzd_t *T1, rci_t *L1, mzd_t *T2, rci_t *L2) {
  wi_t const blocknum = startcol / RADIX;
  wi_t const wide = M->width - blocknum;
  wi_t const count = (wide + 7) / 8;	/* Unrolled loop count */
  int const entry_point = wide % 8;	/* Unrolled loop entry point */

  int rem = k % 3;
  
  int const ka = k / 3 + ((rem >= 2) ? 1 : 0);
  int const kb = k / 3 + ((rem >= 1) ? 1 : 0);
  int const kc = k / 3;

#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(static,512) //if(stoprow-startrow > 128)
#endif
  for(rci_t r = startrow; r < stoprow; ++r) {
    rci_t const x0 = L0[ mzd_read_bits_int(M, r, startcol, ka)];
    rci_t const x1 = L1[ mzd_read_bits_int(M, r, startcol+ka, kb)];
    rci_t const x2 = L2[ mzd_read_bits_int(M, r, startcol+ka+kb, kc)];
    if((x0 | x1 | x2) == 0) // x0 == 0 && x1 == 0 && x2 == 0
      continue;

    word *m0 = M->rows[r] + blocknum;
    word const *t0 = T0->rows[x0] + blocknum;
    word const *t1 = T1->rows[x1] + blocknum;
    word const *t2 = T2->rows[x2] + blocknum;

    wi_t n = count;
    switch (entry_point) {
    case 0: do { *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      case 7:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      case 6:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      case 5:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      case 4:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      case 3:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      case 2:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      case 1:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      } while (--n > 0);
    }
  }
}

void mzd_process_rows4(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k, 
                       mzd_t *T0, rci_t *L0, mzd_t *T1, rci_t *L1, mzd_t *T2, rci_t *L2, mzd_t *T3, rci_t *L3) {
  wi_t const blocknum = startcol / RADIX;
  wi_t const wide = M->width - blocknum;
  wi_t const count = (wide + 7) / 8;	/* Unrolled loop count */
  int const entry_point = wide % 8;	/* Unrolled loop entry point */

  int const rem = k % 4;
  
  int const ka = k / 4 + ((rem >= 3) ? 1 : 0);
  int const kb = k / 4 + ((rem >= 2) ? 1 : 0);
  int const kc = k / 4 + ((rem >= 1) ? 1 : 0);
  int const kd = k / 4;

#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(static,512) //if(stoprow-startrow > 128)
#endif
  for(rci_t r = startrow; r < stoprow; ++r) {
    rci_t const x0 = L0[ mzd_read_bits_int(M, r, startcol, ka)];
    rci_t const x1 = L1[ mzd_read_bits_int(M, r, startcol+ka, kb)];
    rci_t const x2 = L2[ mzd_read_bits_int(M, r, startcol+ka+kb, kc)];
    rci_t const x3 = L3[ mzd_read_bits_int(M, r, startcol+ka+kb+kc, kd)];
    if(((x0 | x1) | (x2 | x3)) == 0) // x0 == 0 && x1 == 0 && x2 == 0 && x3 == 0
      continue;

    word *m0 = M->rows[r] + blocknum;
    word const *t0 = T0->rows[x0] + blocknum;
    word const *t1 = T1->rows[x1] + blocknum;
    word const *t2 = T2->rows[x2] + blocknum;
    word const *t3 = T3->rows[x3] + blocknum;
    
    wi_t n = count;
    switch (entry_point) {
    case 0: do { *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      case 7:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      case 6:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      case 5:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      case 4:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      case 3:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      case 2:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      case 1:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      } while (--n > 0);
    }
  }
}

void mzd_process_rows5(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k, 
                       mzd_t *T0, rci_t *L0, mzd_t *T1, rci_t *L1, mzd_t *T2, rci_t *L2, mzd_t *T3, rci_t *L3,
                       mzd_t *T4, rci_t *L4) {
  wi_t const blocknum = startcol / RADIX;
  wi_t const wide = M->width - blocknum;
  wi_t const count = (wide + 7) / 8;	/* Unrolled loop count */
  int const entry_point = wide % 8;	/* Unrolled loop entry point */
  int rem = k % 5;
  
  int const ka = k / 5 + ((rem >= 4) ? 1 : 0);
  int const kb = k / 5 + ((rem >= 3) ? 1 : 0);
  int const kc = k / 5 + ((rem >= 2) ? 1 : 0);
  int const kd = k / 5 + ((rem >= 1) ? 1 : 0);
  int const ke = k / 5;

#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(static,512) //if(stoprow-startrow > 128)
#endif
  for(rci_t r = startrow; r < stoprow; ++r) {
    
    rci_t const x0 = L0[ mzd_read_bits_int(M, r, startcol, ka)];
    rci_t const x1 = L1[ mzd_read_bits_int(M, r, startcol+ka, kb)];
    rci_t const x2 = L2[ mzd_read_bits_int(M, r, startcol+ka+kb, kc)];
    rci_t const x3 = L3[ mzd_read_bits_int(M, r, startcol+ka+kb+kc, kd)];
    rci_t const x4 = L4[ mzd_read_bits_int(M, r, startcol+ka+kb+kc+kd, ke)];

    if(((x0 | x1 | x2) | (x3 | x4)) == 0) // x0 == 0 && x1 == 0 && x2 == 0 && x3 == 0 && x4 == 0
      continue;

    word *m0 = M->rows[r] + blocknum;
    word const *t0 = T0->rows[x0] + blocknum;
    word const *t1 = T1->rows[x1] + blocknum;
    word const *t2 = T2->rows[x2] + blocknum;
    word const *t3 = T3->rows[x3] + blocknum;
    word const *t4 = T4->rows[x4] + blocknum;
    
    wi_t n = count;
    switch (entry_point) {
    case 0: do { *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      case 7:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      case 6:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      case 5:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      case 4:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      case 3:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      case 2:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      case 1:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      } while (--n > 0);
    }
  }
}

void mzd_process_rows6(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k, 
                       mzd_t *T0, rci_t *L0, mzd_t *T1, rci_t *L1, mzd_t *T2, rci_t *L2, mzd_t *T3, rci_t *L3,
                       mzd_t *T4, rci_t *L4, mzd_t *T5, rci_t *L5) {
  wi_t const blocknum = startcol / RADIX;
  wi_t const wide = M->width - blocknum;
  wi_t const count = (wide + 7) / 8;	/* Unrolled loop count */
  int const entry_point = wide % 8;	/* Unrolled loop entry point */

  int const rem = k % 6;
  
  int const ka = k / 6 + ((rem >= 5) ? 1 : 0);
  int const kb = k / 6 + ((rem >= 4) ? 1 : 0);
  int const kc = k / 6 + ((rem >= 3) ? 1 : 0);
  int const kd = k / 6 + ((rem >= 2) ? 1 : 0);
  int const ke = k / 6 + ((rem >= 1) ? 1 : 0);;
  int const kf = k / 6;

#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(static,512) //if(stoprow-startrow > 128)
#endif
  for(rci_t r = startrow; r < stoprow; ++r) {
    rci_t const x0 = L0[ mzd_read_bits_int(M, r, startcol, ka)];
    rci_t const x1 = L1[ mzd_read_bits_int(M, r, startcol+ka, kb)];
    rci_t const x2 = L2[ mzd_read_bits_int(M, r, startcol+ka+kb, kc)];
    rci_t const x3 = L3[ mzd_read_bits_int(M, r, startcol+ka+kb+kc, kd)];
    rci_t const x4 = L4[ mzd_read_bits_int(M, r, startcol+ka+kb+kc+kd, ke)];
    rci_t const x5 = L5[ mzd_read_bits_int(M, r, startcol+ka+kb+kc+kd+ke, kf)];
    
    /* Waste three clocks on OR-ing (modern CPU can do three in
     * parallel) to avoid possible multiple conditional jumps. */
    if(((x0 | x1) | (x2 | x3) | (x4 | x5)) == 0) // x0 == 0 && x1 == 0 && x2 == 0 && x3 == 0 && x4 == 0 && x5 == 0
      continue;

    word *m0 = M->rows[r] + blocknum;
    word const *t0 = T0->rows[x0] + blocknum;
    word const *t1 = T1->rows[x1] + blocknum;
    word const *t2 = T2->rows[x2] + blocknum;
    word const *t3 = T3->rows[x3] + blocknum;
    word const *t4 = T4->rows[x4] + blocknum;
    word const *t5 = T5->rows[x5] + blocknum;

    wi_t n = count;
    switch (entry_point) {
      case 0: do { *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      case 7:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      case 6:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      case 5:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      case 4:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      case 3:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      case 2:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      case 1:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      } while (--n > 0);
    }
  }
}

rci_t _mzd_echelonize_m4ri(mzd_t *A, int const full, int k, int heuristic, double const threshold) {
  /**
   * \par General algorithm
   * \li Step 1.Denote the first column to be processed in a given
   * iteration as \f$a_i\f$. Then, perform Gaussian elimination on the
   * first \f$3k\f$ rows after and including the \f$i\f$-th row to
   * produce an identity matrix in \f$a_{i,i} ... a_{i+k-1,i+k-1},\f$
   * and zeroes in \f$a_{i+k,i} ... a_{i+3k-1,i+k-1}\f$.
   *
   * \li Step 2. Construct a table consisting of the \f$2^k\f$ binary strings of
   * length k in a Gray code.  Thus with only \f$2^k\f$ vector
   * additions, all possible linear combinations of these k rows
   * have been precomputed.
   *
   * \li Step 3. One can rapidly process the remaining rows from \f$i +
   * 3k\f$ until row \f$m\f$ (the last row) by using the table. For
   * example, suppose the \f$j\f$-th row has entries \f$a_{j,i}
   * ... a_{j,i+k-1}\f$ in the columns being processed. Selecting the
   * row of the table associated with this k-bit string, and adding it
   * to row j will force the k columns to zero, and adjust the
   * remaining columns from \f$ i + k\f$ to n in the appropriate way,
   * as if Gaussian elimination had been performed.
   *
   * \li Step 4. While the above form of the algorithm will reduce a
   * system of boolean linear equations to unit upper triangular form,
   * and thus permit a system to be solved with back substitution, the
   * M4RI algorithm can also be used to invert a matrix, or put the
   * system into reduced row echelon form (RREF). Simply run Step 3
   * on rows \f$0 ... i-1\f$ as well as on rows \f$i + 3k
   * ... m\f$. This only affects the complexity slightly, changing the
   * 2.5 coeffcient to 3.
   *
   * \attention This function implements a variant of the algorithm
   * described above. If heuristic is true, then this algorithm, will
   * switch to PLUQ based echelon form computation once the density
   * reaches the threshold.
   */

  rci_t const ncols = A->ncols; 

  if (k == 0) {
    k = m4ri_opt_k(A->nrows, ncols, 0);
    if (k >= 7)
      k = 7;
    if (0.75 * TWOPOW(k) * ncols > CPU_L2_CACHE / 2.0)
      k -= 1;
  }
  int kk = 6 * k;

  mzd_t *U  = mzd_init(kk, ncols);
  mzd_t *T0 = mzd_init(TWOPOW(k), ncols);
  mzd_t *T1 = mzd_init(TWOPOW(k), ncols);
  mzd_t *T2 = mzd_init(TWOPOW(k), ncols);
  mzd_t *T3 = mzd_init(TWOPOW(k), ncols);
  mzd_t *T4 = mzd_init(TWOPOW(k), ncols);
  mzd_t *T5 = mzd_init(TWOPOW(k), ncols);
  rci_t *L0 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L1 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L2 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L3 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L4 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L5 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));

  rci_t last_check = 0;
  rci_t r = 0;
  rci_t c = 0;

  if (heuristic) {
    if (c < ncols && r < A->nrows && _mzd_density(A, 32, 0, 0) >= threshold) {
      wi_t const tmp = c / RADIX;
      rci_t const tmp2 = tmp * RADIX;
      mzd_t *Abar = mzd_init_window(A, r, tmp2, A->nrows, ncols);
      r += mzd_echelonize_pluq(Abar, full);
      mzd_free(Abar);
      c = ncols;
    }
  }

  while(c < ncols) {
    if (heuristic && c > (last_check + 256)) {
      last_check = c;
      if (c < ncols && r < A->nrows && _mzd_density(A, 32, r, c) >= threshold) {
        mzd_t *Abar = mzd_init_window(A, r, (c / RADIX) * RADIX, A->nrows, ncols);
        if (!full) {
          r += mzd_echelonize_pluq(Abar, full);
        } else {
          rci_t r2 = mzd_echelonize_pluq(Abar, full);
          if (r > 0)
            _mzd_top_echelonize_m4ri(A, 0, r, c, r);
          r += r2;
        }
        mzd_free(Abar);
        c = ncols;
        break;
      }
    }

    if(c + kk > ncols) {
      kk = ncols - c;
    }
    int kbar;
    if (full) {
      kbar = _mzd_gauss_submatrix_full(A, r, c, A->nrows, kk);
    } else {
      kbar = _mzd_gauss_submatrix(A, r, c, A->nrows, kk);
      /* this isn't necessary, adapt make_table */
      U = mzd_submatrix(U, A, r, 0, r + kbar, ncols);
      _mzd_gauss_submatrix_top(A, r, c, kbar);
    }

    if (kbar > 5 * k) {
      int const rem = kbar % 6;
      int const ka = kbar / 6 + ((rem >= 5) ? 1 : 0);
      int const kb = kbar / 6 + ((rem >= 4) ? 1 : 0);
      int const kc = kbar / 6 + ((rem >= 3) ? 1 : 0);
      int const kd = kbar / 6 + ((rem >= 2) ? 1 : 0);
      int const ke = kbar / 6 + ((rem >= 1) ? 1 : 0);;
      int const kf = kbar / 6;

      if(full || kbar == kk) {
        mzd_make_table(A, r, c, ka, T0, L0);
        mzd_make_table(A, r+ka, c, kb, T1, L1);
        mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
        mzd_make_table(A, r+ka+kb+kc, c, kd, T3, L3);
        mzd_make_table(A, r+ka+kb+kc+kd, c, ke, T4, L4);
        mzd_make_table(A, r+ka+kb+kc+kd+ke, c, kf, T5, L5);
      }
      if(kbar == kk)
        mzd_process_rows6(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4, T5, L5);
      if(full)
        mzd_process_rows6(A, 0, r, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4, T5, L5);

  } else if (kbar > 4 * k) { 
      int const rem = kbar % 5;
      int const ka = kbar / 5 + ((rem >= 4) ? 1 : 0);
      int const kb = kbar / 5 + ((rem >= 3) ? 1 : 0);
      int const kc = kbar / 5 + ((rem >= 2) ? 1 : 0);
      int const kd = kbar / 5 + ((rem >= 1) ? 1 : 0);
      int const ke = kbar / 5;
      if(full || kbar == kk) {
        mzd_make_table(A, r, c, ka, T0, L0);
        mzd_make_table(A, r+ka, c, kb, T1, L1);
        mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
        mzd_make_table(A, r+ka+kb+kc, c, kd, T3, L3);
        mzd_make_table(A, r+ka+kb+kc+kd, c, ke, T4, L4);
      }
      if(kbar == kk)
        mzd_process_rows5(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4);
      if(full)
        mzd_process_rows5(A, 0, r, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4);
      
    } else if (kbar > 3 * k) {
      int const rem = kbar % 4;
      int const ka = kbar / 4 + ((rem >= 3) ? 1 : 0);
      int const kb = kbar / 4 + ((rem >= 2) ? 1 : 0);
      int const kc = kbar / 4 + ((rem >= 1) ? 1 : 0);
      int const kd = kbar / 4;
      if(full || kbar == kk) {
        mzd_make_table(A, r, c, ka, T0, L0);
        mzd_make_table(A, r+ka, c, kb, T1, L1);
        mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
        mzd_make_table(A, r+ka+kb+kc, c, kd, T3, L3);
      }
      if(kbar == kk)
        mzd_process_rows4(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3);
      if(full)
        mzd_process_rows4(A, 0, r, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3);
      
    } else if (kbar > 2 * k) {
      int const rem = kbar % 3;
      int const ka = kbar / 3 + ((rem >= 2) ? 1 : 0);
      int const kb = kbar / 3 + ((rem >= 1) ? 1 : 0);
      int const kc = kbar / 3;
      if(full || kbar == kk) {
        mzd_make_table(A, r, c, ka, T0, L0);
        mzd_make_table(A, r+ka, c, kb, T1, L1);
        mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
      }
      if(kbar == kk)
        mzd_process_rows3(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1, T2, L2);
      if(full)
        mzd_process_rows3(A, 0, r, c, kbar, T0, L0, T1, L1, T2, L2);
      
    } else if (kbar > k) {
      int const ka = kbar / 2;
      int const kb = kbar - ka;
      if(full || kbar == kk) {
        mzd_make_table(A, r, c, ka, T0, L0);
        mzd_make_table(A, r+ka, c, kb, T1, L1);
      }
      if(kbar == kk)
        mzd_process_rows2(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1);
      if(full)
        mzd_process_rows2(A, 0, r, c, kbar, T0, L0, T1, L1);
      
    } else if(kbar > 0) {
      if(full || kbar == kk) {
        mzd_make_table(A, r, c, kbar, T0, L0);
      }
      if(kbar == kk)
        mzd_process_rows(A, r+kbar, A->nrows, c, kbar, T0, L0);
      if(full)
        mzd_process_rows(A, 0, r, c, kbar, T0, L0);
    }

    if (!full) {
      _mzd_copy_back_rows(A, U, r, c, kbar);
    }

    r += kbar;
    c += kbar;
    if(kk != kbar) {
      rci_t cbar;
      rci_t rbar;
      if (mzd_find_pivot(A, r, c, &rbar, &cbar)) {
        c = cbar;
        mzd_row_swap(A, r, rbar);
      } else {
        break;
      }
      //c++;
    }
  }

  mzd_free(T0);
  m4ri_mm_free(L0);
  mzd_free(T1);
  m4ri_mm_free(L1);
  mzd_free(T2);
  m4ri_mm_free(L2);
  mzd_free(T3);
  m4ri_mm_free(L3);
  mzd_free(T4);
  m4ri_mm_free(L4);
  mzd_free(T5);
  m4ri_mm_free(L5);
  mzd_free(U);
  return r;
}

rci_t _mzd_top_echelonize_m4ri(mzd_t *A, int k, rci_t r, rci_t c, rci_t max_r) {
  rci_t const ncols = A->ncols; 
  int kbar = 0;

  if (k == 0) {
    k = m4ri_opt_k(max_r, A->ncols, 0);
    if (k >= 7)
      k = 7;
    if (0.75 * TWOPOW(k) *A->ncols > CPU_L2_CACHE / 2.0)
      k -= 1;
  }
  int kk = 6 * k;

  mzd_t *U  = mzd_init(kk, A->ncols);
  mzd_t *T0 = mzd_init(TWOPOW(k), A->ncols);
  mzd_t *T1 = mzd_init(TWOPOW(k), A->ncols);
  mzd_t *T2 = mzd_init(TWOPOW(k), A->ncols);
  mzd_t *T3 = mzd_init(TWOPOW(k), A->ncols);
  mzd_t *T4 = mzd_init(TWOPOW(k), A->ncols);
  mzd_t *T5 = mzd_init(TWOPOW(k), A->ncols);
  rci_t *L0 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L1 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L2 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L3 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L4 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L5 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));

  while(c < ncols) {
    if(c+kk > A->ncols) {
      kk = ncols - c;
    }
    kbar = _mzd_gauss_submatrix_full(A, r, c, MIN(A->nrows,r+kk), kk);

    if (kbar > 5 * k) {
      int const rem = kbar % 6;
      int const ka = kbar / 6 + ((rem >= 5) ? 1 : 0);
      int const kb = kbar / 6 + ((rem >= 4) ? 1 : 0);
      int const kc = kbar / 6 + ((rem >= 3) ? 1 : 0);
      int const kd = kbar / 6 + ((rem >= 2) ? 1 : 0);
      int const ke = kbar / 6 + ((rem >= 1) ? 1 : 0);;
      int const kf = kbar / 6;

      mzd_make_table(A, r, c, ka, T0, L0);
      mzd_make_table(A, r+ka, c, kb, T1, L1);
      mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
      mzd_make_table(A, r+ka+kb+kc, c, kd, T3, L3);
      mzd_make_table(A, r+ka+kb+kc+kd, c, ke, T4, L4);
      mzd_make_table(A, r+ka+kb+kc+kd+ke, c, kf, T5, L5);
      mzd_process_rows6(A, 0, MIN(r, max_r), c, kbar, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4, T5, L5);

  } else if (kbar > 4 * k) { 
      int const rem = kbar % 5;
      int const ka = kbar / 5 + ((rem >= 4) ? 1 : 0);
      int const kb = kbar / 5 + ((rem >= 3) ? 1 : 0);
      int const kc = kbar / 5 + ((rem >= 2) ? 1 : 0);
      int const kd = kbar / 5 + ((rem >= 1) ? 1 : 0);
      int const ke = kbar / 5;

      mzd_make_table(A, r, c, ka, T0, L0);
      mzd_make_table(A, r+ka, c, kb, T1, L1);
      mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
      mzd_make_table(A, r+ka+kb+kc, c, kd, T3, L3);
      mzd_make_table(A, r+ka+kb+kc+kd, c, ke, T4, L4);
      mzd_process_rows5(A, 0, MIN(r, max_r), c, kbar, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4);
      
    } else if (kbar > 3 * k) {
      const int rem = kbar%4;
      const int ka = kbar/4 + ((rem >= 3) ? 1 : 0);
      const int kb = kbar/4 + ((rem >= 2) ? 1 : 0);
      const int kc = kbar/4 + ((rem >= 1) ? 1 : 0);
      const int kd = kbar/4;

      mzd_make_table(A, r, c, ka, T0, L0);
      mzd_make_table(A, r+ka, c, kb, T1, L1);
      mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
      mzd_make_table(A, r+ka+kb+kc, c, kd, T3, L3);
      mzd_process_rows4(A, 0, MIN(r, max_r), c, kbar, T0, L0, T1, L1, T2, L2, T3, L3);
      
    } else if (kbar > 2 * k) {
      const int rem = kbar%3;
      const int ka = kbar/3 + ((rem >= 2) ? 1 : 0);
      const int kb = kbar/3 + ((rem >= 1) ? 1 : 0);
      const int kc = kbar/3;

      mzd_make_table(A, r, c, ka, T0, L0);
      mzd_make_table(A, r+ka, c, kb, T1, L1);
      mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
      mzd_process_rows3(A, 0, MIN(r, max_r), c, kbar, T0, L0, T1, L1, T2, L2);
      
    } else if (kbar > k) {
      const int ka = kbar/2;
      const int kb = kbar - ka;
      mzd_make_table(A, r, c, ka, T0, L0);
      mzd_make_table(A, r+ka, c, kb, T1, L1);
      mzd_process_rows2(A, 0, MIN(r, max_r), c, kbar, T0, L0, T1, L1);
      
    } else if(kbar > 0) {
      mzd_make_table(A, r, c, kbar, T0, L0);
      mzd_process_rows(A, 0, MIN(r, max_r), c, kbar, T0, L0);
    }

    r += kbar;
    c += kbar;
    if(kk != kbar) {
      c++;
    }
  }

  mzd_free(T0);
  m4ri_mm_free(L0);
  mzd_free(T1);
  m4ri_mm_free(L1);
  mzd_free(T2);
  m4ri_mm_free(L2);
  mzd_free(T3);
  m4ri_mm_free(L3);
  mzd_free(T4);
  m4ri_mm_free(L4);
  mzd_free(T5);
  m4ri_mm_free(L5);
  mzd_free(U);
  return r;
}

void mzd_top_echelonize_m4ri(mzd_t *M, int k) {
  _mzd_top_echelonize_m4ri(M,k,0,0,M->nrows);
}


mzd_t *mzd_invert_m4ri(mzd_t *m, mzd_t *I, int k) {
  mzd_t *big = mzd_concat(NULL, m, I);
  rci_t size = m->ncols;
  if (k == 0)
    k = m4ri_opt_k(m->nrows, m->ncols, 0);
  
  mzd_echelonize_m4ri(big, TRUE, k);
  
  mzd_t *answer;
  rci_t i;
  for(i = 0; i < size; ++i) {
    if (!mzd_read_bit(big, i,i )) {
      answer = NULL;
      break;
    }
  }
  if (i == size)
    answer = mzd_submatrix(NULL, big, 0, size, size, 2 * size);
  
  mzd_free(big);
  
  return answer;
}

mzd_t *mzd_mul_m4rm(mzd_t *C, mzd_t *A, mzd_t *B, int k) {
  rci_t a = A->nrows;
  rci_t c = B->ncols;

  if(A->ncols != B->nrows) 
    m4ri_die("mzd_mul_m4rm: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);
  if (C == NULL) {
    C = mzd_init(a, c);
  } else {
    if (C->nrows != a || C->ncols != c)
      m4ri_die("mzd_mul_m4rm: C (%d x %d) has wrong dimensions.\n", C->nrows, C->ncols);
  }
  return _mzd_mul_m4rm(C, A, B, k, TRUE);
}

mzd_t *mzd_addmul_m4rm(mzd_t *C, mzd_t *A, mzd_t *B, int k) {
  rci_t a = A->nrows;
  rci_t c = B->ncols;

  if(C->ncols == 0 || C->nrows == 0)
    return C;

  if(A->ncols != B->nrows) 
    m4ri_die("mzd_mul_m4rm A ncols (%d) need to match B nrows (%d) .\n", A->ncols, B->nrows);
  if (C == NULL) {
    C = mzd_init(a, c);
  } else {
    if (C->nrows != a || C->ncols != c)
      m4ri_die("mzd_mul_m4rm: C has wrong dimensions.\n");
  }
  return _mzd_mul_m4rm(C, A, B, k, FALSE);
}

mzd_t *_mzd_mul_m4rm(mzd_t *C, mzd_t *A, mzd_t *B, int k, int clear) {
  /**
   * The algorithm proceeds as follows:
   * 
   * Step 1. Make a Gray code table of all the \f$2^k\f$ linear combinations
   * of the \f$k\f$ rows of \f$B_i\f$.  Call the \f$x\f$-th row
   * \f$T_x\f$.
   *
   * Step 2. Read the entries 
   *    \f$a_{j,(i-1)k+1}, a_{j,(i-1)k+2} , ... , a_{j,(i-1)k+k}.\f$
   *
   * Let \f$x\f$ be the \f$k\f$ bit binary number formed by the
   * concatenation of \f$a_{j,(i-1)k+1}, ... , a_{j,ik}\f$.
   *
   * Step 3. for \f$h = 1,2, ... , c\f$ do
   *   calculate \f$C_{jh} = C_{jh} + T_{xh}\f$.
   */
  assert(A->offset == 0);
  assert(B->offset == 0);
  assert(C->offset == 0);
  rci_t x1, x2, x3, x4;
  word *t1, *t2, *t3, *t4;

#ifdef M4RM_GRAY8
  rci_t x5, x6, x7, x8;
  word *t5, *t6, *t7, *t8;
#endif

  word *c;

  rci_t const a_nr = A->nrows;
  rci_t const a_nc = A->ncols;
  rci_t const b_nc = B->ncols;

  if (b_nc < RADIX-10 || a_nr < 16) {
    if(clear)
      return mzd_mul_naive(C, A, B);
    else
      return mzd_addmul_naive(C, A, B);
  }

  wi_t wide = C->width;

  /* clear first */
  if (clear) {
    mzd_set_ui(C, 0);
  }

  int const blocksize = MZD_MUL_BLOCKSIZE;

  if (k == 0) {
    k = m4ri_opt_k(blocksize, a_nc, b_nc);
#ifdef M4RM_GRAY8
    if (k > 3)
      k -= 2;
    /* reduce k further if that has a chance of hitting L1 */
    size_t const tsize = 0.8 * TWOPOW(k) * b_nc;
    if(CPU_L1_CACHE < tsize && tsize <= 2 * CPU_L1_CACHE)
      k -= 1;
#else
    if (k > 2)
      k -= 1;
#endif
  }

#ifndef M4RM_GRAY8
  rci_t *buffer = (rci_t*)m4ri_mm_malloc(4 * TWOPOW(k) * sizeof(rci_t));
#else
  rci_t *buffer = (rci_t*)m4ri_mm_malloc(8 * TWOPOW(k) * sizeof(rci_t));
#endif

  mzd_t *T1 = mzd_init(TWOPOW(k), b_nc);
  rci_t *L1 = buffer;
  mzd_t *T2 = mzd_init(TWOPOW(k), b_nc);
  rci_t *L2 = buffer + 1*TWOPOW(k);
  mzd_t *T3 = mzd_init(TWOPOW(k), b_nc);
  rci_t *L3 = buffer + 2*TWOPOW(k);
  mzd_t *T4 = mzd_init(TWOPOW(k), b_nc);
  rci_t *L4 = buffer + 3*TWOPOW(k);

#ifdef M4RM_GRAY8
  mzd_t *T5 = mzd_init(TWOPOW(k), b_nc);
  rci_t *L5 = buffer + 4*TWOPOW(k);
  mzd_t *T6 = mzd_init(TWOPOW(k), b_nc);
  rci_t *L6 = buffer + 5*TWOPOW(k);
  mzd_t *T7 = mzd_init(TWOPOW(k), b_nc);
  rci_t *L7 = buffer + 6*TWOPOW(k);
  mzd_t *T8 = mzd_init(TWOPOW(k), b_nc);
  rci_t *L8 = buffer + 7*TWOPOW(k);
#endif

  /* process stuff that fits into multiple of k first, but blockwise (babystep-giantstep)*/
#ifdef M4RM_GRAY8
  int const kk = 8 * k;
#else
  int const kk = 4 * k;
#endif
  rci_t const end = a_nc / kk;

  rci_t giantstep = 0;
  for (; giantstep + blocksize <= a_nr; giantstep += blocksize) {
    for(rci_t i = 0; i < end; ++i) {
      mzd_make_table( B, kk*i, 0, k, T1, L1);
      mzd_make_table( B, kk*i+k, 0, k, T2, L2);
      mzd_make_table( B, kk*i+k+k, 0, k, T3, L3);
      mzd_make_table( B, kk*i+k+k+k, 0, k, T4, L4);
#ifdef M4RM_GRAY8
      mzd_make_table( B, kk*i+k+k+k+k, 0, k, T5, L5);
      mzd_make_table( B, kk*i+k+k+k+k+k, 0, k, T6, L6);
      mzd_make_table( B, kk*i+k+k+k+k+k+k, 0, k, T7, L7);
      mzd_make_table( B, kk*i+k+k+k+k+k+k+k, 0, k, T8, L8);
#endif   

      for(int babystep = 0; babystep < blocksize; ++babystep) {
        rci_t j = giantstep + babystep;
        x1 = L1[ mzd_read_bits_int(A, j, kk*i, k) ];
        x2 = L2[ mzd_read_bits_int(A, j, kk*i+k, k) ];
        x3 = L3[ mzd_read_bits_int(A, j, kk*i+k+k, k) ];
        x4 = L4[ mzd_read_bits_int(A, j, kk*i+k+k+k, k) ];
#ifdef M4RM_GRAY8
        x5 = L5[ mzd_read_bits_int(A, j, kk*i+k+k+k+k, k) ];
        x6 = L6[ mzd_read_bits_int(A, j, kk*i+k+k+k+k+k, k) ];
        x7 = L7[ mzd_read_bits_int(A, j, kk*i+k+k+k+k+k+k, k) ];
        x8 = L8[ mzd_read_bits_int(A, j, kk*i+k+k+k+k+k+k+k, k) ];
#endif
        c = C->rows[j];
        t1 = T1->rows[x1];
        t2 = T2->rows[x2];
        t3 = T3->rows[x3];
        t4 = T4->rows[x4];
#ifdef M4RM_GRAY8
        t5 = T5->rows[x5];
        t6 = T6->rows[x6];
        t7 = T7->rows[x7];
        t8 = T8->rows[x8];
#endif
        _MZD_COMBINE;
      }
    }
  }
  
  for(rci_t i = 0; i < end; ++i) {
    mzd_make_table( B, kk*i, 0, k, T1, L1);
    mzd_make_table( B, kk*i+k, 0, k, T2, L2);
    mzd_make_table( B, kk*i+k+k, 0, k, T3, L3);
    mzd_make_table( B, kk*i+k+k+k, 0, k, T4, L4);
#ifdef M4RM_GRAY8
    mzd_make_table( B, kk*i+k+k+k+k, 0, k, T5, L5);
    mzd_make_table( B, kk*i+k+k+k+k+k, 0, k, T6, L6);
    mzd_make_table( B, kk*i+k+k+k+k+k+k, 0, k, T7, L7);
    mzd_make_table( B, kk*i+k+k+k+k+k+k+k, 0, k, T8, L8);
#endif
    for(int babystep = 0; babystep < a_nr - giantstep; ++babystep) {
      rci_t j = giantstep + babystep;
      x1 = L1[ mzd_read_bits_int(A, j, kk*i, k) ];
      x2 = L2[ mzd_read_bits_int(A, j, kk*i+k, k) ];
      x3 = L3[ mzd_read_bits_int(A, j, kk*i+k+k, k) ];
      x4 = L4[ mzd_read_bits_int(A, j, kk*i+k+k+k, k) ];
#ifdef M4RM_GRAY8
      x5 = L5[ mzd_read_bits_int(A, j, kk*i+k+k+k+k, k) ];
      x6 = L6[ mzd_read_bits_int(A, j, kk*i+k+k+k+k+k, k) ];
      x7 = L7[ mzd_read_bits_int(A, j, kk*i+k+k+k+k+k+k, k) ];
      x8 = L8[ mzd_read_bits_int(A, j, kk*i+k+k+k+k+k+k+k, k) ];
#endif
      c = C->rows[j];
      t1 = T1->rows[x1];
      t2 = T2->rows[x2];
      t3 = T3->rows[x3];
      t4 = T4->rows[x4];
#ifdef M4RM_GRAY8
      t5 = T5->rows[x5];
      t6 = T6->rows[x6];
      t7 = T7->rows[x7];
      t8 = T8->rows[x8];
#endif
      _MZD_COMBINE;
    }
  }

  /* handle stuff that doesn't fit into multiple of kk */
  if (a_nc%kk) {
    rci_t i;
    for (i = kk / k * end; i < a_nc / k; ++i) {
      mzd_make_table( B, k*i, 0, k, T1, L1);
      for(rci_t j = 0; j < a_nr; ++j) {
        x1 = L1[ mzd_read_bits_int(A, j, k*i, k) ];
        c = C->rows[j];
        t1 = T1->rows[x1];
        for(wi_t ii = 0; ii < wide; ++ii) {
          c[ii] ^= t1[ii];
        }
      }
    }
    /* handle stuff that doesn't fit into multiple of k */
    if (a_nc%k) {
      mzd_make_table( B, k*(a_nc/k), 0, a_nc%k, T1, L1);
      for(rci_t j = 0; j < a_nr; ++j) {
        x1 = L1[ mzd_read_bits_int(A, j, k*i, a_nc%k) ];
        c = C->rows[j];
        t1 = T1->rows[x1];
        for(wi_t ii = 0; ii < wide; ++ii) {
          c[ii] ^= t1[ii];
        }
      }
    }
  }

  mzd_free(T1);
  mzd_free(T2);
  mzd_free(T3);
  mzd_free(T4);
#ifdef M4RM_GRAY8
  mzd_free(T5);
  mzd_free(T6);
  mzd_free(T7);
  mzd_free(T8);
#endif
  m4ri_mm_free(buffer);
  return C;
}

/* TRSM */

void _mzd_trsm_upper_left_even_submatrix(mzd_t *U, mzd_t *B, rci_t const start_row, int const k, word const mask_begin, word const mask_end) {
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < i; ++j) {
      if (mzd_read_bit(U, start_row+(k-i-1), start_row+(k-i)+j)) {
        word *a = B->rows[start_row+(k-i-1)];
        word *b = B->rows[start_row+(k-i)+j];

        *a++ ^= *b++ & mask_begin; 
	wi_t ii;
        for(ii = 1; ii + 8 <= B->width - 1; ii += 8) {
          *a++ ^= *b++;
          *a++ ^= *b++;
          *a++ ^= *b++;
          *a++ ^= *b++;
          *a++ ^= *b++;
          *a++ ^= *b++;
          *a++ ^= *b++;
          *a++ ^= *b++;
        }
        switch(B->width - ii) {
        case 8:  *a++ ^= *b++;
        case 7:  *a++ ^= *b++;
        case 6:  *a++ ^= *b++;
        case 5:  *a++ ^= *b++;
        case 4:  *a++ ^= *b++;
        case 3:  *a++ ^= *b++;
        case 2:  *a++ ^= *b++;
        case 1:  *a++ ^= (*b++ & mask_end);
        }
      }
    }
  }
}

//#undef M4RM_GRAY8

void _mzd_trsm_upper_left_even_m4r(mzd_t *U, mzd_t *B, int k) {
  wi_t const wide = B->width;
  int const blocksize = MZD_MUL_BLOCKSIZE;

  word mask_begin = RIGHT_BITMASK(RADIX - B->offset);
  word mask_end = LEFT_BITMASK((B->ncols + B->offset) % RADIX);
  
  if (B->width == 1)
    mask_begin = mask_begin & mask_end;

  if (k == 0) {
    k = m4ri_opt_k(blocksize, B->nrows, B->ncols);
#ifdef M4RM_GRAY8
    if (k > 3)
      k -= 2;
    /* reduce k further if that has a chance of hitting L1 */
    size_t const tsize = (int)(0.8 * (TWOPOW(k) * B->nrows));
    if(CPU_L1_CACHE < tsize && tsize <= 2 * CPU_L1_CACHE)
      k -= 1;
#else
    if (k > 2)
      k -= 1;
#endif
  }

  mzd_t *T0 = mzd_init(TWOPOW(k), B->ncols + B->offset);
  mzd_t *T1 = mzd_init(TWOPOW(k), B->ncols + B->offset);
  mzd_t *T2 = mzd_init(TWOPOW(k), B->ncols + B->offset);
  mzd_t *T3 = mzd_init(TWOPOW(k), B->ncols + B->offset);
#ifdef M4RM_GRAY8
  mzd_t *T4 = mzd_init(TWOPOW(k), B->ncols + B->offset);
  mzd_t *T5 = mzd_init(TWOPOW(k), B->ncols + B->offset);
  mzd_t *T6 = mzd_init(TWOPOW(k), B->ncols + B->offset);
  mzd_t *T7 = mzd_init(TWOPOW(k), B->ncols + B->offset);
#endif

  rci_t *L0 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L1 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L2 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L3 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
#ifdef M4RM_GRAY8
  rci_t *L4 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L5 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L6 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *L7 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
#endif

#ifdef M4RM_GRAY8
  int kk = 8 * k;
#else
  int kk = 4 * k;
#endif

  rci_t i = 0;
  for (; i < B->nrows - kk; i += kk) {

    _mzd_trsm_upper_left_even_submatrix(U, B, B->nrows-i-kk, kk, mask_begin, mask_end);

#ifdef M4RM_GRAY8
    mzd_make_table(B, B->nrows - i - 8*k, 0, k, T7, L7);
    mzd_make_table(B, B->nrows - i - 7*k, 0, k, T6, L6);
    mzd_make_table(B, B->nrows - i - 6*k, 0, k, T5, L5);
    mzd_make_table(B, B->nrows - i - 5*k, 0, k, T4, L4);
#endif
    mzd_make_table(B, B->nrows - i - 4*k, 0, k, T3, L3);
    mzd_make_table(B, B->nrows - i - 3*k, 0, k, T2, L2);
    mzd_make_table(B, B->nrows - i - 2*k, 0, k, T1, L1);
    mzd_make_table(B, B->nrows - i - 1*k, 0, k, T0, L0);

    for(rci_t j = 0; j < B->nrows - i - kk; ++j) {
#ifdef M4RM_GRAY8
      rci_t const x7 = L7[ mzd_read_bits_int(U, j, B->nrows - i - 8*k, k) ];
      rci_t const x6 = L6[ mzd_read_bits_int(U, j, B->nrows - i - 7*k, k) ];
      rci_t const x5 = L5[ mzd_read_bits_int(U, j, B->nrows - i - 6*k, k) ];
      rci_t const x4 = L4[ mzd_read_bits_int(U, j, B->nrows - i - 5*k, k) ];
#endif
      rci_t const x3 = L3[ mzd_read_bits_int(U, j, B->nrows - i - 4*k, k) ];
      rci_t const x2 = L2[ mzd_read_bits_int(U, j, B->nrows - i - 3*k, k) ];
      rci_t const x1 = L1[ mzd_read_bits_int(U, j, B->nrows - i - 2*k, k) ];
      rci_t const x0 = L0[ mzd_read_bits_int(U, j, B->nrows - i - 1*k, k) ];


      word *b = B->rows[j];
#ifdef M4RM_GRAY8
      word *t7 = T7->rows[x7];
      word *t6 = T6->rows[x6];
      word *t5 = T5->rows[x5];
      word *t4 = T4->rows[x4];
#endif
      word *t3 = T3->rows[x3];
      word *t2 = T2->rows[x2];
      word *t1 = T1->rows[x1];
      word *t0 = T0->rows[x0];

#ifdef M4RM_GRAY8
      _mzd_combine8(b, t0, t1, t2, t3, t4, t5, t6, t7, wide);
      //b[wide-1] ^= (t0[wide-1] ^ t1[wide-1] ^ t2[wide-1] ^ t3[wide-1] ^ t4[wide-1] ^ t5[wide-1] ^ t6[wide-1] ^ t7[wide-1]) & mask_end;
#else
      _mzd_combine4(b, t0, t1, t2, t3, wide);
      //b[wide-1] ^= (t0[wide-1] ^ t1[wide-1] ^ t2[wide-1] ^ t3[wide-1]) & mask_end;
#endif
      
    }
  }

  /* handle stuff that doesn't fit in multiples of kk */
  for ( ;i < B->nrows; i += k) {
    if (i > B->nrows - k)
      k = B->nrows - i;

    _mzd_trsm_upper_left_even_submatrix(U, B, B->nrows-i-k, k, mask_begin, mask_end);

    mzd_make_table(B, B->nrows - i - 1*k, 0, k, T0, L0);

    for(rci_t j = 0; j < B->nrows - i - k; ++j) {
      rci_t const x0 = L0[ mzd_read_bits_int(U, j, B->nrows - i - 1*k, k) ];

      word *b = B->rows[j];
      word *t0 = T0->rows[x0];

      for (wi_t ii = 0; ii < wide; ++ii)
        b[ii] ^= t0[ii];
    }
  }
  
  mzd_free(T0);
  mzd_free(T1);
  mzd_free(T2);
  mzd_free(T3);
#ifdef M4RM_GRAY8
  mzd_free(T4);
  mzd_free(T5);
  mzd_free(T6);
  mzd_free(T7);
#endif

  m4ri_mm_free(L0);
  m4ri_mm_free(L1);
  m4ri_mm_free(L2);
  m4ri_mm_free(L3);
#ifdef M4RM_GRAY8
  m4ri_mm_free(L4);
  m4ri_mm_free(L5);
  m4ri_mm_free(L6);
  m4ri_mm_free(L7);
#endif
}
