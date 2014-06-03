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

#include "brilliantrussian.h"
#include "xor.h"
#include "graycode.h"
#include "echelonform.h"
#include "ple_russian.h"

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
  assert(k <= m4ri_radix);
  rci_t start_row = r;
  rci_t j;
  for (j = c; j < c + k; ++j) {
    int found = 0;
    for (rci_t i = start_row; i < end_row; ++i) {
      /* first we need to clear the first columns */
      word const tmp = mzd_read_bits(A, i, c, j - c + 1);
      if(tmp) {
        for (int l = 0; l < j - c; ++l)
          if (__M4RI_GET_BIT(tmp, l))
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
      break;
    }
  }
  __M4RI_DD_MZD(A);
  __M4RI_DD_INT(j - c);
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
      break;
    }
  }
  __M4RI_DD_MZD(A);
  __M4RI_DD_INT(j - c);
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
  __M4RI_DD_MZD(A);
  __M4RI_DD_INT(k);
  return k;
}

static inline void _mzd_copy_back_rows(mzd_t *A, mzd_t const *U, rci_t r, rci_t c, int k) {
  wi_t const startblock = c / m4ri_radix;
  wi_t const width = A->width - startblock;
  for (int i = 0; i < k; ++i) {
    word const *const src = U->rows[i] + startblock;
    word *const dst = A->rows[r+i] + startblock;
    for (wi_t j = 0; j < width; ++j) {
      dst[j] = src[j];
    }
  }
  __M4RI_DD_MZD(A);
}

void mzd_make_table(mzd_t const *M, rci_t r, rci_t c, int k, mzd_t *T, rci_t *L)
{
  wi_t const homeblock = c / m4ri_radix;
  word const mask_end = __M4RI_LEFT_BITMASK(M->ncols % m4ri_radix);
  word const pure_mask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix - (c % m4ri_radix));
  word const mask_begin = (M->width - homeblock != 1) ? pure_mask_begin : pure_mask_begin & mask_end;
  wi_t const wide = M->width - homeblock;

  int const twokay = __M4RI_TWOPOW(k);
  L[0] = 0;
  for (rci_t i = 1; i < twokay; ++i) {
    word *ti = T->rows[i] + homeblock;
    word *ti1 = T->rows[i-1] + homeblock;

    rci_t const rowneeded = r + m4ri_codebook[k]->inc[i - 1];
    int const id = m4ri_codebook[k]->ord[i];
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
  __M4RI_DD_MZD(T);
  __M4RI_DD_RCI_ARRAY(L, twokay);
}

void mzd_process_rows(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k, mzd_t const *T, rci_t const *L) {
  wi_t const block = startcol / m4ri_radix;
  wi_t const wide = M->width - block;
  wi_t const count = (wide + 7) / 8;	/* Unrolled loop count */
  int const entry_point = wide % 8;	/* Unrolled loop entry point */

  if(k == 1) {
    word const bm = m4ri_one << (startcol % m4ri_radix);

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
    __M4RI_DD_MZD(M);
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

  __M4RI_DD_MZD(M);
}

void mzd_process_rows2(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k,
                       mzd_t const *T0, rci_t const *L0, mzd_t const *T1, rci_t const *L1) {
  assert(k <= m4ri_radix);
  wi_t const blocknum = startcol / m4ri_radix;
  wi_t const wide = M->width - blocknum;

  int const ka = k / 2;
  int const kb = k - k / 2;

  rci_t r;

  word const ka_bm = __M4RI_LEFT_BITMASK(ka);
  word const kb_bm = __M4RI_LEFT_BITMASK(kb);

#if __M4RI_HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(static,512) // MAX((__M4RI_CPU_L1_CACHE >> 3) / wide,
#endif
  for(r = startrow; r < stoprow; ++r) {
    word bits = mzd_read_bits(M, r, startcol, k);
    rci_t const x0 = L0[ bits & ka_bm ]; bits>>=ka;
    rci_t const x1 = L1[ bits & kb_bm ];
    if((x0 | x1) == 0)	// x0 == 0 && x1 == 0
      continue;
    word *m0 = M->rows[r] + blocknum;
    word const *t[2];
    t[0] = T0->rows[x0] + blocknum;
    t[1] = T1->rows[x1] + blocknum;

    _mzd_combine_2( m0, t, wide);
  }

  __M4RI_DD_MZD(M);
}

void mzd_process_rows3(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k,
                       mzd_t const *T0, rci_t const *L0, mzd_t const *T1, rci_t const *L1, mzd_t const *T2, rci_t const *L2) {
  assert(k <= m4ri_radix);
  wi_t const blocknum = startcol / m4ri_radix;
  wi_t const wide = M->width - blocknum;

  int rem = k % 3;

  int const ka = k / 3 + ((rem >= 2) ? 1 : 0);
  int const kb = k / 3 + ((rem >= 1) ? 1 : 0);
  int const kc = k / 3;

  rci_t r;

  word const ka_bm = __M4RI_LEFT_BITMASK(ka);
  word const kb_bm = __M4RI_LEFT_BITMASK(kb);
  word const kc_bm = __M4RI_LEFT_BITMASK(kc);

#if __M4RI_HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(static,512) //if(stoprow-startrow > 128)
#endif
  for(r= startrow; r < stoprow; ++r) {
    word bits = mzd_read_bits(M, r, startcol, k);
    rci_t const x0 = L0[ bits & ka_bm ]; bits>>=ka;
    rci_t const x1 = L1[ bits & kb_bm ]; bits>>=kb;
    rci_t const x2 = L2[ bits & kc_bm ];
    if((x0 | x1 | x2) == 0) // x0 == 0 && x1 == 0 && x2 == 0
      continue;

    word *m0 = M->rows[r] + blocknum;
    word const *t[3];
    t[0] = T0->rows[x0] + blocknum;
    t[1] = T1->rows[x1] + blocknum;
    t[2] = T2->rows[x2] + blocknum;

    _mzd_combine_3( m0, t, wide);
  }

  __M4RI_DD_MZD(M);
}

void mzd_process_rows4(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k,
                       mzd_t const *T0, rci_t const *L0, mzd_t const *T1, rci_t const *L1, mzd_t const *T2, rci_t const *L2, 
                       mzd_t const *T3, rci_t const *L3) {
  assert(k <= m4ri_radix);
  wi_t const blocknum = startcol / m4ri_radix;
  wi_t const wide = M->width - blocknum;

  int const rem = k % 4;

  int const ka = k / 4 + ((rem >= 3) ? 1 : 0);
  int const kb = k / 4 + ((rem >= 2) ? 1 : 0);
  int const kc = k / 4 + ((rem >= 1) ? 1 : 0);
  int const kd = k / 4;

  rci_t r;

  word const ka_bm = __M4RI_LEFT_BITMASK(ka);
  word const kb_bm = __M4RI_LEFT_BITMASK(kb);
  word const kc_bm = __M4RI_LEFT_BITMASK(kc);
  word const kd_bm = __M4RI_LEFT_BITMASK(kd);

#if __M4RI_HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(static,512) //if(stoprow-startrow > 128)
#endif
  for(r = startrow; r < stoprow; ++r) {
    word bits = mzd_read_bits(M, r, startcol, k);
    rci_t const x0 = L0[ bits & ka_bm ]; bits>>=ka;
    rci_t const x1 = L1[ bits & kb_bm ]; bits>>=kb;
    rci_t const x2 = L2[ bits & kc_bm ]; bits>>=kc;
    rci_t const x3 = L3[ bits & kd_bm ];
    if(((x0 | x1) | (x2 | x3)) == 0) // x0 == 0 && x1 == 0 && x2 == 0 && x3 == 0
      continue;

    word *m0 = M->rows[r] + blocknum;
    word const *t[4];
    t[0] = T0->rows[x0] + blocknum;
    t[1] = T1->rows[x1] + blocknum;
    t[2] = T2->rows[x2] + blocknum;
    t[3] = T3->rows[x3] + blocknum;

    _mzd_combine_4( m0, t, wide);
  }

  __M4RI_DD_MZD(M);
}

void mzd_process_rows5(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k,
                       mzd_t const *T0, rci_t const *L0, mzd_t const *T1, rci_t const *L1, mzd_t const *T2, rci_t const *L2,
		       mzd_t const *T3, rci_t const *L3, mzd_t const *T4, rci_t const *L4) {
  assert(k <= m4ri_radix);
  wi_t const blocknum = startcol / m4ri_radix;
  wi_t const wide = M->width - blocknum;
  int rem = k % 5;

  int const ka = k / 5 + ((rem >= 4) ? 1 : 0);
  int const kb = k / 5 + ((rem >= 3) ? 1 : 0);
  int const kc = k / 5 + ((rem >= 2) ? 1 : 0);
  int const kd = k / 5 + ((rem >= 1) ? 1 : 0);
  int const ke = k / 5;

  rci_t r;

  word const ka_bm = __M4RI_LEFT_BITMASK(ka);
  word const kb_bm = __M4RI_LEFT_BITMASK(kb);
  word const kc_bm = __M4RI_LEFT_BITMASK(kc);
  word const kd_bm = __M4RI_LEFT_BITMASK(kd);
  word const ke_bm = __M4RI_LEFT_BITMASK(ke);

#if __M4RI_HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(static,512) //if(stoprow-startrow > 128)
#endif
  for(r = startrow; r < stoprow; ++r) {
    word bits = mzd_read_bits(M, r, startcol, k);
    rci_t const x0 = L0[ bits & ka_bm ]; bits>>=ka;
    rci_t const x1 = L1[ bits & kb_bm ]; bits>>=kb;
    rci_t const x2 = L2[ bits & kc_bm ]; bits>>=kc;
    rci_t const x3 = L3[ bits & kd_bm ]; bits>>=kd;
    rci_t const x4 = L4[ bits & ke_bm ];

    if(((x0 | x1 | x2) | (x3 | x4)) == 0) // x0 == 0 && x1 == 0 && x2 == 0 && x3 == 0 && x4 == 0
      continue;

    word *m0 = M->rows[r] + blocknum;
    word const *t[5];
    t[0] = T0->rows[x0] + blocknum;
    t[1] = T1->rows[x1] + blocknum;
    t[2] = T2->rows[x2] + blocknum;
    t[3] = T3->rows[x3] + blocknum;
    t[4] = T4->rows[x4] + blocknum;

    _mzd_combine_5( m0, t, wide);
  }

  __M4RI_DD_MZD(M);
}

void mzd_process_rows6(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k,
                       mzd_t const *T0, rci_t const *L0, mzd_t const *T1, rci_t const *L1, mzd_t const *T2,
		       rci_t const *L2, mzd_t const *T3, rci_t const *L3, mzd_t const *T4, rci_t const *L4,
		       mzd_t const *T5, rci_t const *L5) {
  assert(k <= m4ri_radix);
  wi_t const blocknum = startcol / m4ri_radix;
  wi_t const wide = M->width - blocknum;

  int const rem = k % 6;

  int const ka = k / 6 + ((rem >= 5) ? 1 : 0);
  int const kb = k / 6 + ((rem >= 4) ? 1 : 0);
  int const kc = k / 6 + ((rem >= 3) ? 1 : 0);
  int const kd = k / 6 + ((rem >= 2) ? 1 : 0);
  int const ke = k / 6 + ((rem >= 1) ? 1 : 0);;
  int const kf = k / 6;

  rci_t r;

  word const ka_bm = __M4RI_LEFT_BITMASK(ka);
  word const kb_bm = __M4RI_LEFT_BITMASK(kb);
  word const kc_bm = __M4RI_LEFT_BITMASK(kc);
  word const kd_bm = __M4RI_LEFT_BITMASK(kd);
  word const ke_bm = __M4RI_LEFT_BITMASK(ke);
  word const kf_bm = __M4RI_LEFT_BITMASK(kf);

#if __M4RI_HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(static,512) //if(stoprow-startrow > 128)
#endif
  for(r = startrow; r < stoprow; ++r) {
    word bits = mzd_read_bits(M, r, startcol, k);
    rci_t const x0 = L0[ bits & ka_bm ]; bits>>=ka;
    rci_t const x1 = L1[ bits & kb_bm ]; bits>>=kb;
    rci_t const x2 = L2[ bits & kc_bm ]; bits>>=kc;
    rci_t const x3 = L3[ bits & kd_bm ]; bits>>=kd;
    rci_t const x4 = L4[ bits & ke_bm ]; bits>>=ke;
    rci_t const x5 = L5[ bits & kf_bm ];

    /* Waste three clocks on OR-ing (modern CPU can do three in
     * parallel) to avoid possible multiple conditional jumps. */
    if(((x0 | x1) | (x2 | x3) | (x4 | x5)) == 0) // x0 == 0 && x1 == 0 && x2 == 0 && x3 == 0 && x4 == 0 && x5 == 0
      continue;

    word *m0 = M->rows[r] + blocknum;
    word const *t[6];
    t[0] = T0->rows[x0] + blocknum;
    t[1] = T1->rows[x1] + blocknum;
    t[2] = T2->rows[x2] + blocknum;
    t[3] = T3->rows[x3] + blocknum;
    t[4] = T4->rows[x4] + blocknum;
    t[5] = T5->rows[x5] + blocknum;

    _mzd_combine_6( m0, t, wide);
  }

  __M4RI_DD_MZD(M);
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
    if (0.75 * __M4RI_TWOPOW(k) * ncols > __M4RI_CPU_L3_CACHE / 2.0)
      k -= 1;
  }
  int kk = 6 * k;

  mzd_t *U  = mzd_init(kk, ncols);
  mzd_t *T0 = mzd_init(__M4RI_TWOPOW(k), ncols);
  mzd_t *T1 = mzd_init(__M4RI_TWOPOW(k), ncols);
  mzd_t *T2 = mzd_init(__M4RI_TWOPOW(k), ncols);
  mzd_t *T3 = mzd_init(__M4RI_TWOPOW(k), ncols);
  mzd_t *T4 = mzd_init(__M4RI_TWOPOW(k), ncols);
  mzd_t *T5 = mzd_init(__M4RI_TWOPOW(k), ncols);
  rci_t *L0 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L1 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L2 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L3 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L4 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L5 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));

  rci_t last_check = 0;
  rci_t r = 0;
  rci_t c = 0;

  if (heuristic) {
    if (c < ncols && r < A->nrows && _mzd_density(A, 32, 0, 0) >= threshold) {
      wi_t const tmp = c / m4ri_radix;
      rci_t const tmp2 = tmp * m4ri_radix;
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
        mzd_t *Abar = mzd_init_window(A, r, (c / m4ri_radix) * m4ri_radix, A->nrows, ncols);
        if (!full) {
          r += mzd_echelonize_pluq(Abar, full);
        } else {
          rci_t r2 = mzd_echelonize_pluq(Abar, full);
          if (r > 0)
            _mzd_top_echelonize_m4ri(A, 0, r, c, r);
          r += r2;
        }
        mzd_free(Abar);
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

  __M4RI_DD_MZD(A);
  __M4RI_DD_RCI(r);
  return r;
}

rci_t _mzd_top_echelonize_m4ri(mzd_t *A, int k, rci_t r, rci_t c, rci_t max_r) {
  rci_t const ncols = A->ncols;
  int kbar = 0;

  if (k == 0) {
    k = m4ri_opt_k(max_r, A->ncols, 0);
    if (k >= 7)
      k = 7;
    if (0.75 * __M4RI_TWOPOW(k) *A->ncols > __M4RI_CPU_L3_CACHE / 2.0)
      k -= 1;
  }
  int kk = 6 * k;

  mzd_t *U  = mzd_init(kk, A->ncols);
  mzd_t *T0 = mzd_init(__M4RI_TWOPOW(k), A->ncols);
  mzd_t *T1 = mzd_init(__M4RI_TWOPOW(k), A->ncols);
  mzd_t *T2 = mzd_init(__M4RI_TWOPOW(k), A->ncols);
  mzd_t *T3 = mzd_init(__M4RI_TWOPOW(k), A->ncols);
  mzd_t *T4 = mzd_init(__M4RI_TWOPOW(k), A->ncols);
  mzd_t *T5 = mzd_init(__M4RI_TWOPOW(k), A->ncols);
  rci_t *L0 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L1 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L2 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L3 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L4 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L5 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));

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

  __M4RI_DD_MZD(A);
  __M4RI_DD_RCI(r);
  return r;
}

void mzd_top_echelonize_m4ri(mzd_t *M, int k) {
  _mzd_top_echelonize_m4ri(M,k,0,0,M->nrows);
}

mzd_t *mzd_inv_m4ri(mzd_t *B, mzd_t const* A, int k) {
  assert(A->nrows == A->ncols);
  if(B == NULL) {
    B = mzd_init(A->nrows, A->ncols);
  } else {
    assert(B->ncols == A->ncols && B->nrows && A->ncols);
  }

  const rci_t n  = A->nrows;
  const rci_t nr = m4ri_radix * A->width;
  mzd_t *C = mzd_init(n, 2*nr);

  mzd_t *AW = mzd_init_window(C, 0, 0,  n, n);
  mzd_t *BW = mzd_init_window(C, 0, nr, n, nr+n);

  mzd_copy(AW, A);
  mzd_set_ui(BW, 1);

  mzd_echelonize_m4ri(C, TRUE, 0);

  mzd_copy(B, BW);
  mzd_free_window(AW);
  mzd_free_window(BW);
  mzd_free(C);
  __M4RI_DD_MZD(B);
  return B;
}


mzd_t *mzd_mul_m4rm(mzd_t *C, mzd_t const *A, mzd_t const *B, int k) {
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

mzd_t *mzd_addmul_m4rm(mzd_t *C, mzd_t const *A, mzd_t const *B, int k) {
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

#define __M4RI_M4RM_NTABLES 8

mzd_t *_mzd_mul_m4rm(mzd_t *C, mzd_t const *A, mzd_t const *B, int k, int clear) {
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

  rci_t        x[__M4RI_M4RM_NTABLES];
  rci_t       *L[__M4RI_M4RM_NTABLES];
  word  const *t[__M4RI_M4RM_NTABLES];
  mzd_t       *T[__M4RI_M4RM_NTABLES];
#ifdef __M4RI_HAVE_SSE2
  mzd_t  *Talign[__M4RI_M4RM_NTABLES];
  int c_align = (__M4RI_ALIGNMENT(C->rows[0], 16) == 8);
#endif

  word *c;

  rci_t const a_nr = A->nrows;
  rci_t const a_nc = A->ncols;
  rci_t const b_nc = B->ncols;

  if (b_nc < m4ri_radix-10 || a_nr < 16) {
    if(clear)
      return mzd_mul_naive(C, A, B);
    else
      return mzd_addmul_naive(C, A, B);
  }

  /* clear first */
  if (clear) {
    mzd_set_ui(C, 0);
  }

  const int blocksize = __M4RI_MUL_BLOCKSIZE;

  if(k==0) {
    /* __M4RI_CPU_L2_CACHE == 2^k * B->width * 8 * 8 */
    k = (int)log2((__M4RI_CPU_L2_CACHE/64)/(double)B->width);
    if ((__M4RI_CPU_L2_CACHE - 64*__M4RI_TWOPOW(k)*B->width) > (64*__M4RI_TWOPOW(k+1)*B->width - __M4RI_CPU_L2_CACHE))
      k++;

    rci_t const klog = round(0.75 * log2_floor(MIN(MIN(a_nr,a_nc),b_nc)));

    if(klog < k)
      k = klog;

    if (k<2)
      k=2;
    else if(k>6)
      k=6;
  }

  const wi_t wide = C->width;
  const word bm = __M4RI_TWOPOW(k)-1;

  rci_t *buffer = (rci_t*)m4ri_mm_malloc(__M4RI_M4RM_NTABLES * __M4RI_TWOPOW(k) * sizeof(rci_t));
  for(int z=0; z<__M4RI_M4RM_NTABLES; z++) {
    L[z] = buffer + z*__M4RI_TWOPOW(k);
#ifdef __M4RI_HAVE_SSE2
    /* we make sure that T are aligned as C */
    Talign[z] = mzd_init(__M4RI_TWOPOW(k), b_nc+m4ri_radix);
    T[z] = mzd_init_window(Talign[z], 0, c_align*m4ri_radix, Talign[z]->nrows, b_nc + c_align*m4ri_radix);
#else
    T[z] = mzd_init(__M4RI_TWOPOW(k), b_nc);
#endif
  }

  /* process stuff that fits into multiple of k first, but blockwise (babystep-giantstep)*/
  int const kk = __M4RI_M4RM_NTABLES * k;
  assert(kk <= m4ri_radix);
  rci_t const end = a_nc / kk;

  for (rci_t giantstep = 0; giantstep < a_nr; giantstep += blocksize) {
    for(rci_t i = 0; i < end; ++i) {
#if __M4RI_HAVE_OPENMP
#pragma omp parallel for schedule(static,1)
#endif
      for(int z=0; z<__M4RI_M4RM_NTABLES; z++) {
        mzd_make_table( B, kk*i + k*z, 0, k, T[z], L[z]);
      }

      const rci_t blockend = MIN(giantstep+blocksize, a_nr);
#if __M4RI_HAVE_OPENMP
#pragma omp parallel for schedule(static,512) private(x,t)
#endif
      for(rci_t j = giantstep; j < blockend; j++) {
        const word a = mzd_read_bits(A, j, kk*i, kk);

        switch(__M4RI_M4RM_NTABLES) {
        case 8: t[7] = T[ 7]->rows[ L[7][ (a >> 7*k) & bm ] ];
        case 7: t[6] = T[ 6]->rows[ L[6][ (a >> 6*k) & bm ] ];
        case 6: t[5] = T[ 5]->rows[ L[5][ (a >> 5*k) & bm ] ];
        case 5: t[4] = T[ 4]->rows[ L[4][ (a >> 4*k) & bm ] ];
        case 4: t[3] = T[ 3]->rows[ L[3][ (a >> 3*k) & bm ] ];
        case 3: t[2] = T[ 2]->rows[ L[2][ (a >> 2*k) & bm ] ];
        case 2: t[1] = T[ 1]->rows[ L[1][ (a >> 1*k) & bm ] ];
        case 1: t[0] = T[ 0]->rows[ L[0][ (a >> 0*k) & bm ] ];
          break;
        default:
          m4ri_die("__M4RI_M4RM_NTABLES must be <= 8 but got %d", __M4RI_M4RM_NTABLES);
        }

        c = C->rows[j];

        switch(__M4RI_M4RM_NTABLES) {
        case 8: _mzd_combine_8(c, t, wide); break;
        case 7: _mzd_combine_7(c, t, wide); break;
        case 6: _mzd_combine_6(c, t, wide); break;
        case 5: _mzd_combine_5(c, t, wide); break;
        case 4: _mzd_combine_4(c, t, wide); break;
        case 3: _mzd_combine_3(c, t, wide); break;
        case 2: _mzd_combine_2(c, t, wide); break;
        case 1: _mzd_combine(c, t[0], wide);
          break;
        default:
          m4ri_die("__M4RI_M4RM_NTABLES must be <= 8 but got %d", __M4RI_M4RM_NTABLES);
        }
      }
    }
  }

  /* handle stuff that doesn't fit into multiple of kk */
  if (a_nc%kk) {
    rci_t i;
    for (i = kk / k * end; i < a_nc / k; ++i) {
      mzd_make_table( B, k*i, 0, k, T[0], L[0]);
      for(rci_t j = 0; j < a_nr; ++j) {
        x[0] = L[0][ mzd_read_bits_int(A, j, k*i, k) ];
        c = C->rows[j];
        t[0] = T[0]->rows[x[0]];
        for(wi_t ii = 0; ii < wide; ++ii) {
          c[ii] ^= t[0][ii];
        }
      }
    }
    /* handle stuff that doesn't fit into multiple of k */
    if (a_nc%k) {
      mzd_make_table( B, k*(a_nc/k), 0, a_nc%k, T[0], L[0]);
      for(rci_t j = 0; j < a_nr; ++j) {
        x[0] = L[0][ mzd_read_bits_int(A, j, k*i, a_nc%k) ];
        c = C->rows[j];
        t[0] = T[0]->rows[x[0]];
        for(wi_t ii = 0; ii < wide; ++ii) {
          c[ii] ^= t[0][ii];
        }
      }
    }
  }

  for(int j=0; j<__M4RI_M4RM_NTABLES; j++) {
    mzd_free(T[j]);
#ifdef __M4RI_HAVE_SSE2
    mzd_free(Talign[j]);
#endif
  }
  m4ri_mm_free(buffer);

  __M4RI_DD_MZD(C);
  return C;
}

