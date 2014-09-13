
/*******************************************************************
*
*                 M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2014 Martin Albrecht <martinralbrecht@googlemail.com>
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

#include "m4ri_config.h"
#include "misc.h"
#include "mp.h"
#include "brilliantrussian.h"
#include "strassen.h"

#if __M4RI_HAVE_OPENMP

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#include <omp.h>

// Returns true if a is closer to cutoff than a/2.
static inline int closer(rci_t a, int cutoff) {
  return 3 * a < 4 * cutoff;
}


mzd_t *_mzd_addmul_mp4(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff) {
  /**
   * \todo make sure not to overwrite crap after ncols and before width * m4ri_radix
   */
  rci_t a = A->nrows;
  rci_t b = A->ncols;
  rci_t c = B->ncols;
  /* handle case first, where the input matrices are too small already */
  if (closer(A->nrows, cutoff) || closer(A->ncols, cutoff) || closer(B->ncols, cutoff)) {
    /* we copy the matrix first since it is only constant memory
       overhead and improves data locality, if you remove it make sure
       there are no speed regressions */
    /* C = _mzd_mul_m4rm(C, A, B, 0, TRUE); */
    mzd_t *Cbar = mzd_init(C->nrows, C->ncols);
    Cbar = _mzd_mul_m4rm(Cbar, A, B, 0, FALSE);
    mzd_add(C, C, Cbar);
    mzd_free(Cbar);
    return C;
  }

  /* adjust cutting numbers to work on words */
  {
    rci_t mult = 2 * m4ri_radix;
    a -= a % mult;
    b -= b % mult;
    c -= c % mult;
  }

  rci_t anr = ((a / m4ri_radix) >> 1) * m4ri_radix;
  rci_t anc = ((b / m4ri_radix) >> 1) * m4ri_radix;
  rci_t bnr = anc;
  rci_t bnc = ((c / m4ri_radix) >> 1) * m4ri_radix;

  mzd_t const *A00 = mzd_init_window_const(A,   0,   0,   anr,   anc);
  mzd_t const *A01 = mzd_init_window_const(A,   0, anc,   anr, 2*anc);
  mzd_t const *A10 = mzd_init_window_const(A, anr,   0, 2*anr,   anc);
  mzd_t const *A11 = mzd_init_window_const(A, anr, anc, 2*anr, 2*anc);

  mzd_t const *B00 = mzd_init_window_const(B,   0,   0,   bnr,   bnc);
  mzd_t const *B01 = mzd_init_window_const(B,   0, bnc,   bnr, 2*bnc);
  mzd_t const *B10 = mzd_init_window_const(B, bnr,   0, 2*bnr,   bnc);
  mzd_t const *B11 = mzd_init_window_const(B, bnr, bnc, 2*bnr, 2*bnc);

  mzd_t *C00 = mzd_init_window(C,   0,   0,   anr,   bnc);
  mzd_t *C01 = mzd_init_window(C,   0, bnc,   anr, 2*bnc);
  mzd_t *C10 = mzd_init_window(C, anr,   0, 2*anr,   bnc);
  mzd_t *C11 = mzd_init_window(C, anr, bnc, 2*anr, 2*bnc);

#pragma omp parallel sections
  {
#pragma omp section
    {
      _mzd_addmul_even(C00, A00, B00, cutoff);
      _mzd_addmul_even(C00, A01, B10, cutoff);
    }
#pragma omp section
    {
      _mzd_addmul_even(C01, A00, B01, cutoff);
      _mzd_addmul_even(C01, A01, B11, cutoff);
    }
#pragma omp section
    {
      _mzd_addmul_even(C10, A10, B00, cutoff);
      _mzd_addmul_even(C10, A11, B10, cutoff);
    }
#pragma omp section
    {
      _mzd_addmul_even(C11, A10, B01, cutoff);
      _mzd_addmul_even(C11, A11, B11, cutoff);
    }
  }

  /* deal with rest */
  if (B->ncols > 2 * bnc) {
    mzd_t const *B_last_col = mzd_init_window_const(B, 0, 2*bnc, A->ncols, B->ncols);
    mzd_t *C_last_col = mzd_init_window(C, 0, 2*bnc, A->nrows, C->ncols);
    mzd_addmul_m4rm(C_last_col, A, B_last_col, 0);
    mzd_free_window((mzd_t*)B_last_col);
    mzd_free_window(C_last_col);
  }
  if (A->nrows > 2 * anr) {
    mzd_t const *A_last_row = mzd_init_window_const(A, 2*anr, 0, A->nrows, A->ncols);
    mzd_t const *B_bulk = mzd_init_window_const(B, 0, 0, B->nrows, 2*bnc);
    mzd_t *C_last_row = mzd_init_window(C, 2*anr, 0, C->nrows, 2*bnc);
    mzd_addmul_m4rm(C_last_row, A_last_row, B_bulk, 0);
    mzd_free_window((mzd_t*)A_last_row);
    mzd_free_window((mzd_t*)B_bulk);
    mzd_free_window(C_last_row);
  }
  if (A->ncols > 2 * anc) {
    mzd_t const *A_last_col = mzd_init_window_const(A,     0, 2*anc, 2*anr, A->ncols);
    mzd_t const *B_last_row = mzd_init_window_const(B, 2*bnr,     0, B->nrows, 2*bnc);
    mzd_t *C_bulk = mzd_init_window(C, 0, 0, 2*anr, 2*bnc);
    mzd_addmul_m4rm(C_bulk, A_last_col, B_last_row, 0);
    mzd_free_window((mzd_t*)A_last_col);
    mzd_free_window((mzd_t*)B_last_row);
    mzd_free_window(C_bulk);
  }

  /* clean up */
  mzd_free_window((mzd_t*)A00); mzd_free_window((mzd_t*)A01);
  mzd_free_window((mzd_t*)A10); mzd_free_window((mzd_t*)A11);

  mzd_free_window((mzd_t*)B00); mzd_free_window((mzd_t*)B01);
  mzd_free_window((mzd_t*)B10); mzd_free_window((mzd_t*)B11);

  mzd_free_window(C00); mzd_free_window(C01);
  mzd_free_window(C10); mzd_free_window(C11);

  __M4RI_DD_MZD(C);
  return C;
}

mzd_t *_mzd_mul_mp4(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff) {
  /**
   * \todo make sure not to overwrite crap after ncols and before width * m4ri_radix
   */
  rci_t a = A->nrows;
  rci_t b = A->ncols;
  rci_t c = B->ncols;
  /* handle case first, where the input matrices are too small already */
  if (closer(A->nrows, cutoff) || closer(A->ncols, cutoff) || closer(B->ncols, cutoff)) {
    /* we copy the matrix first since it is only constant memory
       overhead and improves data locality, if you remove it make sure
       there are no speed regressions */
    /* C = _mzd_mul_m4rm(C, A, B, 0, TRUE); */
    mzd_t *Cbar = mzd_init(C->nrows, C->ncols);
    Cbar = _mzd_mul_m4rm(Cbar, A, B, 0, FALSE);
    mzd_copy(C, Cbar);
    mzd_free(Cbar);
    return C;
  }

  /* adjust cutting numbers to work on words */
  {
    rci_t mult = 2 * m4ri_radix;
    a -= a % mult;
    b -= b % mult;
    c -= c % mult;
  }

  rci_t anr = ((a / m4ri_radix) >> 1) * m4ri_radix;
  rci_t anc = ((b / m4ri_radix) >> 1) * m4ri_radix;
  rci_t bnr = anc;
  rci_t bnc = ((c / m4ri_radix) >> 1) * m4ri_radix;

  mzd_t const *A00 = mzd_init_window_const(A,   0,   0,   anr,   anc);
  mzd_t const *A01 = mzd_init_window_const(A,   0, anc,   anr, 2*anc);
  mzd_t const *A10 = mzd_init_window_const(A, anr,   0, 2*anr,   anc);
  mzd_t const *A11 = mzd_init_window_const(A, anr, anc, 2*anr, 2*anc);

  mzd_t const *B00 = mzd_init_window_const(B,   0,   0,   bnr,   bnc);
  mzd_t const *B01 = mzd_init_window_const(B,   0, bnc,   bnr, 2*bnc);
  mzd_t const *B10 = mzd_init_window_const(B, bnr,   0, 2*bnr,   bnc);
  mzd_t const *B11 = mzd_init_window_const(B, bnr, bnc, 2*bnr, 2*bnc);

  mzd_t *C00 = mzd_init_window(C,   0,   0,   anr,   bnc);
  mzd_t *C01 = mzd_init_window(C,   0, bnc,   anr, 2*bnc);
  mzd_t *C10 = mzd_init_window(C, anr,   0, 2*anr,   bnc);
  mzd_t *C11 = mzd_init_window(C, anr, bnc, 2*anr, 2*bnc);

#pragma omp parallel sections
  {
#pragma omp section
    {
      _mzd_mul_even(C00, A00, B00, cutoff);
      _mzd_addmul_even(C00, A01, B10, cutoff);
    }
#pragma omp section
    {
      _mzd_mul_even(C01, A00, B01, cutoff);
      _mzd_addmul_even(C01, A01, B11, cutoff);
    }
#pragma omp section
    {
      _mzd_mul_even(C10, A10, B00, cutoff);
      _mzd_addmul_even(C10, A11, B10, cutoff);
    }
#pragma omp section
    {
      _mzd_mul_even(C11, A10, B01, cutoff);
      _mzd_addmul_even(C11, A11, B11, cutoff);
    }
  }

  /* deal with rest */
  if (B->ncols > 2 * bnc) {
    mzd_t const *B_last_col = mzd_init_window_const(B, 0, 2*bnc, A->ncols, B->ncols);
    mzd_t *C_last_col = mzd_init_window(C, 0, 2*bnc, A->nrows, C->ncols);
    mzd_addmul_m4rm(C_last_col, A, B_last_col, 0);
    mzd_free_window((mzd_t*)B_last_col);
    mzd_free_window(C_last_col);
  }
  if (A->nrows > 2 * anr) {
    mzd_t const *A_last_row = mzd_init_window_const(A, 2*anr, 0, A->nrows, A->ncols);
    mzd_t const *B_bulk = mzd_init_window_const(B, 0, 0, B->nrows, 2*bnc);
    mzd_t *C_last_row = mzd_init_window(C, 2*anr, 0, C->nrows, 2*bnc);
    mzd_addmul_m4rm(C_last_row, A_last_row, B_bulk, 0);
    mzd_free_window((mzd_t*)A_last_row);
    mzd_free_window((mzd_t*)B_bulk);
    mzd_free_window(C_last_row);
  }
  if (A->ncols > 2 * anc) {
    mzd_t const *A_last_col = mzd_init_window_const(A,     0, 2*anc, 2*anr, A->ncols);
    mzd_t const *B_last_row = mzd_init_window_const(B, 2*bnr,     0, B->nrows, 2*bnc);
    mzd_t *C_bulk = mzd_init_window(C, 0, 0, 2*anr, 2*bnc);
    mzd_addmul_m4rm(C_bulk, A_last_col, B_last_row, 0);
    mzd_free_window((mzd_t*)A_last_col);
    mzd_free_window((mzd_t*)B_last_row);
    mzd_free_window(C_bulk);
  }

  /* clean up */
  mzd_free_window((mzd_t*)A00); mzd_free_window((mzd_t*)A01);
  mzd_free_window((mzd_t*)A10); mzd_free_window((mzd_t*)A11);

  mzd_free_window((mzd_t*)B00); mzd_free_window((mzd_t*)B01);
  mzd_free_window((mzd_t*)B10); mzd_free_window((mzd_t*)B11);

  mzd_free_window(C00); mzd_free_window(C01);
  mzd_free_window(C10); mzd_free_window(C11);

  __M4RI_DD_MZD(C);
  return C;
}

mzd_t *mzd_mul_mp(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff) {
  if(A->ncols != B->nrows)
    m4ri_die("mzd_mul_mp: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);

  if (cutoff < 0)
    m4ri_die("mzd_mul_mp: cutoff must be >= 0.\n");

  if(cutoff == 0) {
    cutoff = __M4RI_STRASSEN_MUL_CUTOFF;
  }

  cutoff = cutoff / m4ri_radix * m4ri_radix;
  if (cutoff < m4ri_radix) {
    cutoff = m4ri_radix;
  };

  if (C == NULL) {
    C = mzd_init(A->nrows, B->ncols);
  } else if (C->nrows != A->nrows || C->ncols != B->ncols){
    m4ri_die("mzd_mul_mp: C (%d x %d) has wrong dimensions, expected (%d x %d)\n",
	     C->nrows, C->ncols, A->nrows, B->ncols);
  }

  _mzd_mul_mp4(C, A, B, cutoff);
  return C;
}


mzd_t *mzd_addmul_mp(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff) {
  if(A->ncols != B->nrows)
    m4ri_die("mzd_addmul_mp: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);

  if (cutoff < 0)
    m4ri_die("mzd_addmul_mp: cutoff must be >= 0.\n");

  if(cutoff == 0) {
    cutoff = __M4RI_STRASSEN_MUL_CUTOFF;
  }

  cutoff = cutoff / m4ri_radix * m4ri_radix;
  if (cutoff < m4ri_radix) {
    cutoff = m4ri_radix;
  };

  if (C == NULL) {
    C = mzd_init(A->nrows, B->ncols);
  } else if (C->nrows != A->nrows || C->ncols != B->ncols){
    m4ri_die("mzd_addmul_mp: C (%d x %d) has wrong dimensions, expected (%d x %d)\n",
	     C->nrows, C->ncols, A->nrows, B->ncols);
  }
  if(A->nrows == 0 || A->ncols == 0 || B->ncols == 0) {
    __M4RI_DD_MZD(C);
    return C;
  }

  C = _mzd_addmul_mp4(C, A, B, cutoff);
  __M4RI_DD_MZD(C);
  return C;
}


#endif //__M4RI_HAVE_OPENMP
