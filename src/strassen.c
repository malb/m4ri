/*******************************************************************
*
*                 M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2008 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
*    Copyright (C) 2008 Clement Pernet <pernet@math.washington.edu>
*    Copyright (C) 2008 Marco Bodrato <bodrato@mail.dm.unipi.it>
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

#include "graycode.h"
#include "strassen.h"
#include "parity.h"
#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#if __M4RI_HAVE_OPENMP
#include <omp.h>
#endif

// Returns true if a is closer to cutoff than a/2.
static inline int closer(rci_t a, int cutoff) {
  return 3 * a < 4 * cutoff;
}


mzd_t *_mzd_mul_even(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff) {
  rci_t mmm, kkk, nnn;

  if(C->nrows == 0 || C->ncols == 0)
    return C;

  rci_t m = A->nrows;
  rci_t k = A->ncols;
  rci_t n = B->ncols;

  /* handle case first, where the input matrices are too small already */
  if (closer(m, cutoff) || closer(k, cutoff) || closer(n, cutoff)) {
    /* we copy the matrices first since it is only constant memory overhead and improves data
       locality */
    if(mzd_is_windowed(A)|mzd_is_windowed(B)|mzd_is_windowed(C)) {
      mzd_t *Abar = mzd_copy(NULL, A);
      mzd_t *Bbar = mzd_copy(NULL, B);
      mzd_t *Cbar = mzd_init(m, n);
      _mzd_mul_m4rm(Cbar, Abar, Bbar, 0, FALSE);
      mzd_copy(C, Cbar);
      mzd_free(Cbar);
      mzd_free(Bbar);
      mzd_free(Abar);
    } else {
      _mzd_mul_m4rm(C, A, B, 0, TRUE);
    }
    return C;
  }

  /* adjust cutting numbers to work on words */

  rci_t mult = m4ri_radix;
  rci_t width = MIN(MIN(m, n), k) / 2;
  while (width > cutoff) {
    width /= 2;
    mult *= 2;
  }

  mmm = (((m - m % mult) / m4ri_radix) >> 1) * m4ri_radix;
  kkk = (((k - k % mult) / m4ri_radix) >> 1) * m4ri_radix;
  nnn = (((n - n % mult) / m4ri_radix) >> 1) * m4ri_radix;

  /*         |A |   |B |   |C |
   * Compute |  | x |  | = |  | */
  {
    mzd_t const *A11 = mzd_init_window_const(A,   0,   0,   mmm,   kkk);
    mzd_t const *A12 = mzd_init_window_const(A,   0, kkk,   mmm, 2*kkk);
    mzd_t const *A21 = mzd_init_window_const(A, mmm,   0, 2*mmm,   kkk);
    mzd_t const *A22 = mzd_init_window_const(A, mmm, kkk, 2*mmm, 2*kkk);

    mzd_t const *B11 = mzd_init_window_const(B,   0,   0,   kkk,   nnn);
    mzd_t const *B12 = mzd_init_window_const(B,   0, nnn,   kkk, 2*nnn);
    mzd_t const *B21 = mzd_init_window_const(B, kkk,   0, 2*kkk,   nnn);
    mzd_t const *B22 = mzd_init_window_const(B, kkk, nnn, 2*kkk, 2*nnn);

    mzd_t *C11 = mzd_init_window(C,   0,   0,   mmm,   nnn);
    mzd_t *C12 = mzd_init_window(C,   0, nnn,   mmm, 2*nnn);
    mzd_t *C21 = mzd_init_window(C, mmm,   0, 2*mmm,   nnn);
    mzd_t *C22 = mzd_init_window(C, mmm, nnn, 2*mmm, 2*nnn);

    /**
     * \note See Marco Bodrato; "A Strassen-like Matrix Multiplication
     * Suited for Squaring and Highest Power Computation";
     * http://bodrato.it/papres/#CIVV2008 for reference on the used
     * sequence of operations.
     */

    /* change this to mzd_init(mmm, MAX(nnn,kkk)) to fix the todo below */
    mzd_t *Wmk = mzd_init(mmm, kkk);
    mzd_t *Wkn = mzd_init(kkk, nnn);

    _mzd_add(Wkn, B22, B12);		 /* Wkn = B22 + B12 */
    _mzd_add(Wmk, A22, A12);		 /* Wmk = A22 + A12 */
    _mzd_mul_even(C21, Wmk, Wkn, cutoff);/* C21 = Wmk * Wkn */

    _mzd_add(Wmk, A22, A21);		 /* Wmk = A22 - A21 */
    _mzd_add(Wkn, B22, B21);		 /* Wkn = B22 - B21 */
    _mzd_mul_even(C22, Wmk, Wkn, cutoff);/* C22 = Wmk * Wkn */

    _mzd_add(Wkn, Wkn, B12);		 /* Wkn = Wkn + B12 */
    _mzd_add(Wmk, Wmk, A12);		 /* Wmk = Wmk + A12 */
    _mzd_mul_even(C11, Wmk, Wkn, cutoff);/* C11 = Wmk * Wkn */

    _mzd_add(Wmk, Wmk, A11);		 /* Wmk = Wmk - A11 */
    _mzd_mul_even(C12, Wmk, B12, cutoff);/* C12 = Wmk * B12 */
    _mzd_add(C12, C12, C22);		 /* C12 = C12 + C22 */

    /**
     * \todo ideally we would use the same Wmk throughout the function
     * but some called function doesn't like that and we end up with a
     * wrong result if we use virtual Wmk matrices. Ideally, this should
     * be fixed not worked around. The check whether the bug has been
     * fixed, use only one Wmk and check if mzd_mul(4096, 3528,
     * 4096, 2124) still returns the correct answer.
     */

    mzd_free(Wmk);
    Wmk = mzd_mul(NULL, A12, B21, cutoff);/*Wmk = A12 * B21 */

    _mzd_add(C11, C11, Wmk);		  /* C11 = C11 + Wmk */
    _mzd_add(C12, C11, C12);		  /* C12 = C11 - C12 */
    _mzd_add(C11, C21, C11);		  /* C11 = C21 - C11 */
    _mzd_add(Wkn, Wkn, B11);		  /* Wkn = Wkn - B11 */
    _mzd_mul_even(C21, A21, Wkn, cutoff); /* C21 = A21 * Wkn */
    mzd_free(Wkn);

    _mzd_add(C21, C11, C21);		  /* C21 = C11 - C21 */
    _mzd_add(C22, C22, C11);		  /* C22 = C22 + C11 */
    _mzd_mul_even(C11, A11, B11, cutoff); /* C11 = A11 * B11 */

    _mzd_add(C11, C11, Wmk);		  /* C11 = C11 + Wmk */

    /* clean up */
    mzd_free_window((mzd_t*)A11); mzd_free_window((mzd_t*)A12);
    mzd_free_window((mzd_t*)A21); mzd_free_window((mzd_t*)A22);

    mzd_free_window((mzd_t*)B11); mzd_free_window((mzd_t*)B12);
    mzd_free_window((mzd_t*)B21); mzd_free_window((mzd_t*)B22);

    mzd_free_window(C11); mzd_free_window(C12);
    mzd_free_window(C21); mzd_free_window(C22);

    mzd_free(Wmk);
  }
  /* deal with rest */
  nnn *= 2;
  if (n > nnn) {
    /*         |AA|   | B|   | C|
     * Compute |AA| x | B| = | C| */
    mzd_t const *B_last_col = mzd_init_window_const(B, 0, nnn, k, n);
    mzd_t *C_last_col = mzd_init_window(C, 0, nnn, m, n);
    _mzd_mul_m4rm(C_last_col, A, B_last_col, 0, TRUE);
    mzd_free_window((mzd_t*)B_last_col);
    mzd_free_window(C_last_col);
  }
  mmm *= 2;
  if (m > mmm) {
    /*         |  |   |B |   |  |
     * Compute |AA| x |B | = |C | */
    mzd_t const *A_last_row = mzd_init_window_const(A, mmm, 0, m, k);
    mzd_t const *B_first_col= mzd_init_window_const(B,   0, 0, k, nnn);
    mzd_t *C_last_row = mzd_init_window(C, mmm, 0, m, nnn);
    _mzd_mul_m4rm(C_last_row, A_last_row, B_first_col, 0, TRUE);
    mzd_free_window((mzd_t*)A_last_row);
    mzd_free_window((mzd_t*)B_first_col);
    mzd_free_window(C_last_row);
  }
  kkk *= 2;
  if (k > kkk) {
    /* Add to  |  |   | B|   |C |
     * result  |A | x |  | = |  | */
    mzd_t const *A_last_col = mzd_init_window_const(A,   0, kkk, mmm, k);
    mzd_t const *B_last_row = mzd_init_window_const(B, kkk,   0,   k, nnn);
    mzd_t *C_bulk = mzd_init_window(C, 0, 0, mmm, nnn);
    mzd_addmul_m4rm(C_bulk, A_last_col, B_last_row, 0);
    mzd_free_window((mzd_t*)A_last_col);
    mzd_free_window((mzd_t*)B_last_row);
    mzd_free_window(C_bulk);
  }

  __M4RI_DD_MZD(C);
  return C;
}

mzd_t *_mzd_sqr_even(mzd_t *C, mzd_t const *A, int cutoff) {
  rci_t m;

  m = A->nrows;
  /* handle case first, where the input matrices are too small already */
  if (closer(m, cutoff)) {
    /* we copy the matrices first since it is only constant memory overhead and improves data
       locality */
    if(mzd_is_windowed(A)|mzd_is_windowed(C)) {
      mzd_t *Abar = mzd_copy(NULL, A);
      mzd_t *Cbar = mzd_init(m, m);
      _mzd_mul_m4rm(Cbar, Abar, Abar, 0, FALSE);
      mzd_copy(C, Cbar);
      mzd_free(Cbar);
      mzd_free(Abar);
    } else {
      _mzd_mul_m4rm(C, A, A, 0, TRUE);
    }
    return C;
  }

  /* adjust cutting numbers to work on words */
  rci_t mmm;
  {
    rci_t mult = m4ri_radix;
    rci_t width = m / 2;
    while (width > cutoff) {
      width /= 2;
      mult *= 2;
    }
    mmm = (((m - m % mult) / m4ri_radix) >> 1) * m4ri_radix;
  }
  /*         |A |   |A |   |C |
   * Compute |  | x |  | = |  | */
  {
    mzd_t const *A11 = mzd_init_window_const(A,   0,   0,   mmm,   mmm);
    mzd_t const *A12 = mzd_init_window_const(A,   0, mmm,   mmm, 2*mmm);
    mzd_t const *A21 = mzd_init_window_const(A, mmm,   0, 2*mmm,   mmm);
    mzd_t const *A22 = mzd_init_window_const(A, mmm, mmm, 2*mmm, 2*mmm);

    mzd_t *C11 = mzd_init_window(C,   0,   0,   mmm,   mmm);
    mzd_t *C12 = mzd_init_window(C,   0, mmm,   mmm, 2*mmm);
    mzd_t *C21 = mzd_init_window(C, mmm,   0, 2*mmm,   mmm);
    mzd_t *C22 = mzd_init_window(C, mmm, mmm, 2*mmm, 2*mmm);

    /**
     * \note See Marco Bodrato; "A Strassen-like Matrix Multiplication
     * Suited for Squaring and Highest Power Computation";
     * http://bodrato.it/papres/#CIVV2008 for reference on the used
     * sequence of operations.
     */

    mzd_t *Wmk;
    mzd_t *Wkn = mzd_init(mmm, mmm);

    _mzd_add(Wkn, A22, A12);                 /* Wkn = A22 + A12 */
    _mzd_sqr_even(C21, Wkn, cutoff);     /* C21 = Wkn^2 */

    _mzd_add(Wkn, A22, A21);                 /* Wkn = A22 - A21 */
    _mzd_sqr_even(C22, Wkn, cutoff);     /* C22 = Wkn^2 */

    _mzd_add(Wkn, Wkn, A12);                 /* Wkn = Wkn + A12 */
    _mzd_sqr_even(C11, Wkn, cutoff);     /* C11 = Wkn^2 */

    _mzd_add(Wkn, Wkn, A11);                 /* Wkn = Wkn - A11 */
    _mzd_mul_even(C12, Wkn, A12, cutoff);/* C12 = Wkn * A12 */
    _mzd_add(C12, C12, C22);		  /* C12 = C12 + C22 */

    Wmk = mzd_mul(NULL, A12, A21, cutoff);/*Wmk = A12 * A21 */

    _mzd_add(C11, C11, Wmk);		  /* C11 = C11 + Wmk */
    _mzd_add(C12, C11, C12);		  /* C12 = C11 - C12 */
    _mzd_add(C11, C21, C11);		  /* C11 = C21 - C11 */
    _mzd_mul_even(C21, A21, Wkn, cutoff);/* C21 = A21 * Wkn */
    mzd_free(Wkn);

    _mzd_add(C21, C11, C21);		  /* C21 = C11 - C21 */
    _mzd_add(C22, C22, C11);		  /* C22 = C22 + C11 */
    _mzd_sqr_even(C11, A11, cutoff);     /* C11 = A11^2 */

    _mzd_add(C11, C11, Wmk);		  /* C11 = C11 + Wmk */

    /* clean up */
    mzd_free_window((mzd_t*)A11); mzd_free_window((mzd_t*)A12);
    mzd_free_window((mzd_t*)A21); mzd_free_window((mzd_t*)A22);

    mzd_free_window(C11); mzd_free_window(C12);
    mzd_free_window(C21); mzd_free_window(C22);

    mzd_free(Wmk);
  }
  /* deal with rest */
  mmm *= 2;
  if (m > mmm) {
    /*         |AA|   | A|   | C|
     * Compute |AA| x | A| = | C| */
    {
      mzd_t const *A_last_col = mzd_init_window_const(A, 0, mmm, m, m);
      mzd_t *C_last_col = mzd_init_window(C, 0, mmm, m, m);
      _mzd_mul_m4rm(C_last_col, A, A_last_col, 0, TRUE);
      mzd_free_window((mzd_t*)A_last_col);
      mzd_free_window(C_last_col);
    }
    /*         |  |   |A |   |  |
     * Compute |AA| x |A | = |C | */
    {
      mzd_t const *A_last_row = mzd_init_window_const(A, mmm, 0, m, m);
      mzd_t const *A_first_col= mzd_init_window_const(A,   0, 0, m, mmm);
      mzd_t *C_last_row = mzd_init_window(C, mmm, 0, m, mmm);
      _mzd_mul_m4rm(C_last_row, A_last_row, A_first_col, 0, TRUE);
      mzd_free_window((mzd_t*)A_last_row);
      mzd_free_window((mzd_t*)A_first_col);
      mzd_free_window(C_last_row);
    }
    /* Add to  |  |   | A|   |C |
     * result  |A | x |  | = |  | */
    {
      mzd_t const *A_last_col = mzd_init_window_const(A,   0, mmm, mmm, m);
      mzd_t const *A_last_row = mzd_init_window_const(A, mmm,   0,   m, mmm);
      mzd_t *C_bulk = mzd_init_window(C, 0, 0, mmm, mmm);
      mzd_addmul_m4rm(C_bulk, A_last_col, A_last_row, 0);
      mzd_free_window((mzd_t*)A_last_col);
      mzd_free_window((mzd_t*)A_last_row);
      mzd_free_window(C_bulk);
    }
  }

  __M4RI_DD_MZD(C);
  return C;
}



mzd_t *mzd_mul(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff) {
  if(A->ncols != B->nrows)
    m4ri_die("mzd_mul: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);

  if (cutoff < 0)
    m4ri_die("mzd_mul: cutoff must be >= 0.\n");

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
    m4ri_die("mzd_mul: C (%d x %d) has wrong dimensions, expected (%d x %d)\n",
	     C->nrows, C->ncols, A->nrows, B->ncols);
  }

  C = (A == B) ? _mzd_sqr_even(C, A, cutoff) : _mzd_mul_even(C, A, B, cutoff);
  return C;
}

mzd_t *_mzd_addmul_even(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff) {
  /**
   * \todo make sure not to overwrite crap after ncols and before width * m4ri_radix
   */
  if(C->nrows == 0 || C->ncols == 0)
    return C;

  rci_t m = A->nrows;
  rci_t k = A->ncols;
  rci_t n = B->ncols;

  /* handle case first, where the input matrices are too small already */
  if (closer(m, cutoff) || closer(k, cutoff) || closer(n, cutoff)) {
    /* we copy the matrices first since it is only constant memory overhead and improves data
       locality */
    if(mzd_is_windowed(A)|mzd_is_windowed(B)|mzd_is_windowed(C)) {
      mzd_t *Abar = mzd_copy(NULL, A);
      mzd_t *Bbar = mzd_copy(NULL, B);
      mzd_t *Cbar = mzd_copy(NULL, C);
      mzd_addmul_m4rm(Cbar, Abar, Bbar, 0);
      mzd_copy(C, Cbar);
      mzd_free(Cbar);
      mzd_free(Bbar);
      mzd_free(Abar);
    } else {
      mzd_addmul_m4rm(C, A, B, 0);
    }
    return C;
  }

  /* adjust cutting numbers to work on words */
  rci_t mmm, kkk, nnn;
  {
    rci_t mult = m4ri_radix;
    rci_t width = MIN(MIN(m, n), k) / 2;
    while (width > cutoff) {
      width /= 2;
      mult *= 2;
    }

    mmm = (((m - m % mult) / m4ri_radix) >> 1) * m4ri_radix;
    kkk = (((k - k % mult) / m4ri_radix) >> 1) * m4ri_radix;
    nnn = (((n - n % mult) / m4ri_radix) >> 1) * m4ri_radix;
  }

  /*         |C |    |A |   |B |
   * Compute |  | += |  | x |  |  */
  {
    mzd_t const *A11 = mzd_init_window_const(A,   0,   0,   mmm,   kkk);
    mzd_t const *A12 = mzd_init_window_const(A,   0, kkk,   mmm, 2*kkk);
    mzd_t const *A21 = mzd_init_window_const(A, mmm,   0, 2*mmm,   kkk);
    mzd_t const *A22 = mzd_init_window_const(A, mmm, kkk, 2*mmm, 2*kkk);

    mzd_t const *B11 = mzd_init_window_const(B,   0,   0,   kkk,   nnn);
    mzd_t const *B12 = mzd_init_window_const(B,   0, nnn,   kkk, 2*nnn);
    mzd_t const *B21 = mzd_init_window_const(B, kkk,   0, 2*kkk,   nnn);
    mzd_t const *B22 = mzd_init_window_const(B, kkk, nnn, 2*kkk, 2*nnn);

    mzd_t *C11 = mzd_init_window(C,   0,   0,   mmm,   nnn);
    mzd_t *C12 = mzd_init_window(C,   0, nnn,   mmm, 2*nnn);
    mzd_t *C21 = mzd_init_window(C, mmm,   0, 2*mmm,   nnn);
    mzd_t *C22 = mzd_init_window(C, mmm, nnn, 2*mmm, 2*nnn);

    /**
     * \note See Marco Bodrato; "A Strassen-like Matrix Multiplication
     * Suited for Squaring and Highest Power Computation";
     * http://bodrato.it/papres/#CIVV2008 for reference on the used
     * sequence of operations.
     */

    mzd_t *S = mzd_init(mmm, kkk);
    mzd_t *T = mzd_init(kkk, nnn);
    mzd_t *U = mzd_init(mmm, nnn);

    _mzd_add(S, A22, A21);                   /* 1  S = A22 - A21       */
    _mzd_add(T, B22, B21);                   /* 2  T = B22 - B21       */
    _mzd_mul_even(U, S, T, cutoff);          /* 3  U = S*T             */
    _mzd_add(C22, U, C22);                   /* 4  C22 = U + C22       */
    _mzd_add(C12, U, C12);                   /* 5  C12 = U + C12       */

    _mzd_mul_even(U, A12, B21, cutoff);      /* 8  U = A12*B21         */
    _mzd_add(C11, U, C11);                   /* 9  C11 = U + C11       */

    _mzd_addmul_even(C11, A11, B11, cutoff); /* 11 C11 = A11*B11 + C11 */

    _mzd_add(S, S, A12);                     /* 6  S = S - A12         */
    _mzd_add(T, T, B12);                     /* 7  T = T - B12         */
    _mzd_addmul_even(U, S, T, cutoff);       /* 10 U = S*T + U         */
    _mzd_add(C12, C12, U);                   /* 15 C12 = U + C12       */

    _mzd_add(S, A11, S);                     /* 12 S = A11 - S         */
    _mzd_addmul_even(C12, S, B12, cutoff);   /* 14 C12 = S*B12 + C12   */

    _mzd_add(T, B11, T);                     /* 13 T = B11 - T         */
    _mzd_addmul_even(C21, A21, T, cutoff);   /* 16 C21 = A21*T + C21   */

    _mzd_add(S, A22, A12);                   /* 17 S = A22 + A21       */
    _mzd_add(T, B22, B12);                   /* 18 T = B22 + B21       */
    _mzd_addmul_even(U, S, T, cutoff);       /* 19 U = U - S*T         */
    _mzd_add(C21, C21, U);                   /* 20 C21 = C21 - U       */
    _mzd_add(C22, C22, U);                   /* 21 C22 = C22 - U       */

    /* clean up */
    mzd_free_window((mzd_t*)A11); mzd_free_window((mzd_t*)A12);
    mzd_free_window((mzd_t*)A21); mzd_free_window((mzd_t*)A22);

    mzd_free_window((mzd_t*)B11); mzd_free_window((mzd_t*)B12);
    mzd_free_window((mzd_t*)B21); mzd_free_window((mzd_t*)B22);

    mzd_free_window(C11); mzd_free_window(C12);
    mzd_free_window(C21); mzd_free_window(C22);

    mzd_free(S);
    mzd_free(T);
    mzd_free(U);
  }
  /* deal with rest */
  nnn *= 2;
  if (n > nnn) {
    /*         | C|    |AA|   | B|
     * Compute | C| += |AA| x | B| */
    mzd_t const *B_last_col = mzd_init_window_const(B, 0, nnn, k, n);
    mzd_t *C_last_col = mzd_init_window(C, 0, nnn, m, n);
    mzd_addmul_m4rm(C_last_col, A, B_last_col, 0);
    mzd_free_window((mzd_t*)B_last_col);
    mzd_free_window(C_last_col);
  }
  mmm *= 2;
  if (m > mmm) {
    /*         |  |    |  |   |B |
     * Compute |C | += |AA| x |B | */
    mzd_t const *A_last_row = mzd_init_window_const(A, mmm, 0, m, k);
    mzd_t const *B_first_col= mzd_init_window_const(B,   0, 0, k, nnn);
    mzd_t *C_last_row = mzd_init_window(C, mmm, 0, m, nnn);
    mzd_addmul_m4rm(C_last_row, A_last_row, B_first_col, 0);
    mzd_free_window((mzd_t*)A_last_row);
    mzd_free_window((mzd_t*)B_first_col);
    mzd_free_window(C_last_row);
  }
  kkk *= 2;
  if (k > kkk) {
    /* Add to  |  |   | B|   |C |
     * result  |A | x |  | = |  | */
    mzd_t const *A_last_col = mzd_init_window_const(A,   0, kkk, mmm, k);
    mzd_t const *B_last_row = mzd_init_window_const(B, kkk,   0,   k, nnn);
    mzd_t *C_bulk = mzd_init_window(C, 0, 0, mmm, nnn);
    mzd_addmul_m4rm(C_bulk, A_last_col, B_last_row, 0);
    mzd_free_window((mzd_t*)A_last_col);
    mzd_free_window((mzd_t*)B_last_row);
    mzd_free_window(C_bulk);
  }

  __M4RI_DD_MZD(C);
  return C;
}

mzd_t *_mzd_addsqr_even(mzd_t *C, mzd_t const *A, int cutoff) {
  /**
   * \todo make sure not to overwrite crap after ncols and before width * m4ri_radix
   */
  if(C->nrows == 0)
    return C;

  rci_t m = A->nrows;

  /* handle case first, where the input matrices are too small already */
  if (closer(m, cutoff)) {
    /* we copy the matrices first since it is only constant memory overhead and improves data
       locality */
    if(mzd_is_windowed(A)|mzd_is_windowed(C)) {
      mzd_t *Cbar = mzd_copy(NULL, C);
      mzd_t *Abar = mzd_copy(NULL, A);
      mzd_addmul_m4rm(Cbar, Abar, Abar, 0);
      mzd_copy(C, Cbar);
      mzd_free(Cbar);
      mzd_free(Abar);
    } else {
      mzd_addmul_m4rm(C, A, A, 0);
    }
    return C;
  }

  /* adjust cutting numbers to work on words */
  rci_t mmm;
  {
    rci_t mult = m4ri_radix;
    rci_t width = m / 2;
    while (width > cutoff) {
      width /= 2;
      mult *= 2;
    }

    mmm = (((m - m % mult) / m4ri_radix) >> 1) * m4ri_radix;
  }

  /*         |C |    |A |   |B |
   * Compute |  | += |  | x |  |  */
  {
    mzd_t const *A11 = mzd_init_window_const(A,   0,   0,   mmm,   mmm);
    mzd_t const *A12 = mzd_init_window_const(A,   0, mmm,   mmm, 2*mmm);
    mzd_t const *A21 = mzd_init_window_const(A, mmm,   0, 2*mmm,   mmm);
    mzd_t const *A22 = mzd_init_window_const(A, mmm, mmm, 2*mmm, 2*mmm);

    mzd_t *C11 = mzd_init_window(C,   0,   0,   mmm,   mmm);
    mzd_t *C12 = mzd_init_window(C,   0, mmm,   mmm, 2*mmm);
    mzd_t *C21 = mzd_init_window(C, mmm,   0, 2*mmm,   mmm);
    mzd_t *C22 = mzd_init_window(C, mmm, mmm, 2*mmm, 2*mmm);

    /**
     * \note See Marco Bodrato; "A Strassen-like Matrix Multiplication
     * Suited for Squaring and Highest Power Computation"; on-line v.
     * http://bodrato.it/papres/#CIVV2008 for reference on the used
     * sequence of operations.
     */

    mzd_t *S = mzd_init(mmm, mmm);
    mzd_t *U = mzd_init(mmm, mmm);

    _mzd_add(S, A22, A21);                   /* 1  S = A22 - A21       */
    _mzd_sqr_even(U, S, cutoff);             /* 3  U = S^2             */
    _mzd_add(C22, U, C22);                   /* 4  C22 = U + C22       */
    _mzd_add(C12, U, C12);                   /* 5  C12 = U + C12       */

    _mzd_mul_even(U, A12, A21, cutoff);      /* 8  U = A12*A21         */
    _mzd_add(C11, U, C11);                   /* 9  C11 = U + C11       */

    _mzd_addsqr_even(C11, A11, cutoff);      /* 11 C11 = A11^2 + C11   */

    _mzd_add(S, S, A12);                     /* 6  S = S + A12         */
    _mzd_addsqr_even(U, S, cutoff);          /* 10 U = S^2 + U         */
    _mzd_add(C12, C12, U);                   /* 15 C12 = U + C12       */

    _mzd_add(S, A11, S);                     /* 12 S = A11 - S         */
    _mzd_addmul_even(C12, S, A12, cutoff);   /* 14 C12 = S*B12 + C12   */

    _mzd_addmul_even(C21, A21, S, cutoff);   /* 16 C21 = A21*T + C21   */

    _mzd_add(S, A22, A12);                   /* 17 S = A22 + A21       */
    _mzd_addsqr_even(U, S, cutoff);          /* 19 U = U - S^2         */
    _mzd_add(C21, C21, U);                   /* 20 C21 = C21 - U3      */
    _mzd_add(C22, C22, U);                   /* 21 C22 = C22 - U3      */

    /* clean up */
    mzd_free_window((mzd_t*)A11); mzd_free_window((mzd_t*)A12);
    mzd_free_window((mzd_t*)A21); mzd_free_window((mzd_t*)A22);

    mzd_free_window(C11); mzd_free_window(C12);
    mzd_free_window(C21); mzd_free_window(C22);

    mzd_free(S);
    mzd_free(U);
  }
  /* deal with rest */
  mmm *= 2;
  if (m > mmm) {
    /*         | C|    |AA|   | B|
     * Compute | C| += |AA| x | B| */
    {
      mzd_t const *A_last_col = mzd_init_window_const(A, 0, mmm, m, m);
      mzd_t *C_last_col = mzd_init_window(C, 0, mmm, m, m);
      mzd_addmul_m4rm(C_last_col, A, A_last_col, 0);
      mzd_free_window((mzd_t*)A_last_col);
      mzd_free_window(C_last_col);
    }
    /*         |  |    |  |   |B |
     * Compute |C | += |AA| x |B | */
    {
      mzd_t const *A_last_row = mzd_init_window_const(A, mmm, 0, m, m);
      mzd_t const *A_first_col= mzd_init_window_const(A,   0, 0, m, mmm);
      mzd_t *C_last_row = mzd_init_window(C, mmm, 0, m, mmm);
      mzd_addmul_m4rm(C_last_row, A_last_row, A_first_col, 0);
      mzd_free_window((mzd_t*)A_last_row);
      mzd_free_window((mzd_t*)A_first_col);
      mzd_free_window(C_last_row);
    }
    /* Add to  |  |   | B|   |C |
     * result  |A | x |  | = |  | */
    {
      mzd_t const *A_last_col = mzd_init_window_const(A,   0, mmm, mmm, m);
      mzd_t const *A_last_row = mzd_init_window_const(A, mmm,   0,   m, mmm);
      mzd_t *C_bulk = mzd_init_window(C, 0, 0, mmm, mmm);
      mzd_addmul_m4rm(C_bulk, A_last_col, A_last_row, 0);
      mzd_free_window((mzd_t*)A_last_col);
      mzd_free_window((mzd_t*)A_last_row);
      mzd_free_window(C_bulk);
    }
  }

  __M4RI_DD_MZD(C);
  return C;
}

mzd_t *_mzd_addmul(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff) {
  /**
   * Assumes that B and C are aligned in the same manner (as in a Schur complement)
   */

  return (A == B) ? _mzd_addsqr_even(C, A, cutoff) : _mzd_addmul_even(C, A, B, cutoff);
}

mzd_t *mzd_addmul(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff) {
  if(A->ncols != B->nrows)
    m4ri_die("mzd_addmul: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);

  if (cutoff < 0)
    m4ri_die("mzd_addmul: cutoff must be >= 0.\n");

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
    m4ri_die("mzd_addmul: C (%d x %d) has wrong dimensions, expected (%d x %d)\n",
	     C->nrows, C->ncols, A->nrows, B->ncols);
  }
  if(A->nrows == 0 || A->ncols == 0 || B->ncols == 0) {
    __M4RI_DD_MZD(C);
    return C;
  }

  C = _mzd_addmul(C, A, B, cutoff);
  __M4RI_DD_MZD(C);
  return C;
}
