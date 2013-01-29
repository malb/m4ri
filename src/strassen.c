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
static inline int CLOSER(rci_t a, int cutoff)
{
  return 3 * a < 4 * cutoff;
}

/**
 * Simple blockwise product
 */
mzd_t *_mzd_addmul_mp_even(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff);

mzd_t *_mzd_mul_even(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff) {
  rci_t mmm, kkk, nnn;

  if(C->nrows == 0 || C->ncols == 0)
    return C;

  rci_t m = A->nrows;
  rci_t k = A->ncols;
  rci_t n = B->ncols;

  /* handle case first, where the input matrices are too small already */
  if (CLOSER(m, cutoff) || CLOSER(k, cutoff) || CLOSER(n, cutoff)) {
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
  if (CLOSER(m, cutoff)) {
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


#if __M4RI_HAVE_OPENMP
mzd_t *_mzd_addmul_mp_even(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff) {
  /**
   * \todo make sure not to overwrite crap after ncols and before width * m4ri_radix
   */
  rci_t a = A->nrows;
  rci_t b = A->ncols;
  rci_t c = B->ncols;
  /* handle case first, where the input matrices are too small already */
  if (CLOSER(A->nrows, cutoff) || CLOSER(A->ncols, cutoff) || CLOSER(B->ncols, cutoff)) {
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
/*    rci_t width = a; */
/*    while (width > 2 * cutoff) { */
/*      width /= 2; */
/*      mult *= 2; */
/*    } */
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

#pragma omp parallel sections num_threads(4)
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
#endif

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

  if(A->offset || B->offset || C->offset) {
    mzd_set_ui(C, 0);
    mzd_addmul(C, A, B, cutoff);
    return C;
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
  if (CLOSER(m, cutoff) || CLOSER(k, cutoff) || CLOSER(n, cutoff)) {
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
  if (CLOSER(m, cutoff)) {
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

  if (!A->offset){
    if (!B->offset) /* A even, B even */
      return (A == B) ? _mzd_addsqr_even(C, A, cutoff) : _mzd_addmul_even(C, A, B, cutoff);
    else {          /* A even, B weird */
      int const bnc = m4ri_radix - B->offset;
      if (B->ncols <= bnc){
	_mzd_addmul_even_weird  (C,  A, B, cutoff);
      } else {
	mzd_t const *B0 = mzd_init_window_const (B, 0, 0, B->nrows, bnc);
	mzd_t *C0 = mzd_init_window (C, 0, 0, C->nrows, bnc);
	mzd_t const *B1 = mzd_init_window_const (B, 0, bnc, B->nrows, B->ncols);
	mzd_t *C1 = mzd_init_window (C, 0, bnc, C->nrows, C->ncols);
	_mzd_addmul_even_weird  (C0,  A, B0, cutoff);
	_mzd_addmul_even(C1, A, B1, cutoff);
	mzd_free_window ((mzd_t*)B0); mzd_free_window ((mzd_t*)B1);
	mzd_free_window (C0); mzd_free_window (C1);
      }
    }
  } else if (B->offset) { /* A weird, B weird */
    int const anc = m4ri_radix - A->offset;
    int const bnc = m4ri_radix - B->offset;
    if (B->ncols <= bnc){
      if (A->ncols <= anc)
	_mzd_addmul_weird_weird (C, A, B);
      else {
	mzd_t const *A0  = mzd_init_window_const (A, 0, 0, A->nrows, anc);
	mzd_t const *A1  = mzd_init_window_const (A, 0, anc, A->nrows, A->ncols);
	mzd_t const *B0  = mzd_init_window_const (B, 0, 0, anc, B->ncols);
	mzd_t const *B1  = mzd_init_window_const (B, anc, 0, B->nrows, B->ncols);
	_mzd_addmul_weird_weird (C, A0, B0);
	_mzd_addmul_even_weird  (C, A1, B1, cutoff);
	mzd_free_window ((mzd_t*)A0);  mzd_free_window ((mzd_t*)A1);
	mzd_free_window ((mzd_t*)B0);  mzd_free_window ((mzd_t*)B1);
      }
    } else if (A->ncols <= anc) {
      mzd_t const *B0 = mzd_init_window_const (B, 0, 0, B->nrows, bnc);
      mzd_t const *B1 = mzd_init_window_const (B, 0, bnc, B->nrows, B->ncols);
      mzd_t *C0 = mzd_init_window (C, 0, 0, C->nrows, bnc);
      mzd_t *C1 = mzd_init_window (C, 0, bnc, C->nrows, C->ncols);
      _mzd_addmul_weird_weird (C0, A, B0);
      _mzd_addmul_weird_even  (C1, A, B1, cutoff);
      mzd_free_window ((mzd_t*)B0); mzd_free_window ((mzd_t*)B1);
      mzd_free_window (C0); mzd_free_window (C1);
    } else {
      mzd_t const *A0  = mzd_init_window_const (A, 0, 0, A->nrows, anc);
      mzd_t const *A1  = mzd_init_window_const (A, 0, anc, A->nrows, A->ncols);
      mzd_t const *B00 = mzd_init_window_const (B, 0, 0, anc, bnc);
      mzd_t const *B01 = mzd_init_window_const (B, 0, bnc, anc, B->ncols);
      mzd_t const *B10 = mzd_init_window_const (B, anc, 0, B->nrows, bnc);
      mzd_t const *B11 = mzd_init_window_const (B, anc, bnc, B->nrows, B->ncols);
      mzd_t *C0 = mzd_init_window (C, 0, 0, C->nrows, bnc);
      mzd_t *C1 = mzd_init_window (C, 0, bnc, C->nrows, C->ncols);

      _mzd_addmul_weird_weird (C0, A0, B00);
      _mzd_addmul_even_weird  (C0,  A1, B10, cutoff);
      _mzd_addmul_weird_even  (C1,  A0, B01, cutoff);
      _mzd_addmul_even  (C1,  A1, B11, cutoff);

      mzd_free_window ((mzd_t*)A0);  mzd_free_window ((mzd_t*)A1);
      mzd_free_window (C0);  mzd_free_window (C1);
      mzd_free_window ((mzd_t*)B00); mzd_free_window ((mzd_t*)B01);
      mzd_free_window ((mzd_t*)B10); mzd_free_window ((mzd_t*)B11);
    }
  } else { /* A weird, B even */
    int const anc = m4ri_radix - A->offset;
    if (A->ncols <= anc){
      _mzd_addmul_weird_even  (C,  A, B, cutoff);
    } else {
      mzd_t const *A0  = mzd_init_window_const (A, 0, 0, A->nrows, anc);
      mzd_t const *A1  = mzd_init_window_const (A, 0, anc, A->nrows, A->ncols);
      mzd_t const *B0  = mzd_init_window_const (B, 0, 0, anc, B->ncols);
      mzd_t const *B1  = mzd_init_window_const (B, anc, 0, B->nrows, B->ncols);
      _mzd_addmul_weird_even (C, A0, B0, cutoff);
      _mzd_addmul_even  (C, A1, B1, cutoff);
      mzd_free_window ((mzd_t*)A0); mzd_free_window ((mzd_t*)A1);
      mzd_free_window ((mzd_t*)B0); mzd_free_window ((mzd_t*)B1);
    }
  }

  __M4RI_DD_MZD(C);
  return C;
}

mzd_t *_mzd_addmul_weird_even (mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff) {
  mzd_t *tmp = mzd_init(A->nrows, MIN((rci_t)m4ri_radix - A->offset, A->ncols));
  for (rci_t i = 0; i < A->nrows; ++i){
    tmp->rows[i][0] = (A->rows[i][0] >> A->offset);
  }
  _mzd_addmul_even (C, tmp, B, cutoff);
  mzd_free(tmp);
  return C;
}

mzd_t *_mzd_addmul_even_weird (mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff) {
  assert(B->width == 1 && C->width == 1);
  assert(mzd_is_windowed(B));	// Otherwise the whole copying makes no sense.
  // D will contain the first 64 columns of B.
  mzd_t *D = mzd_init(B->nrows, (rci_t)m4ri_radix);
  // Make a backup of the shape of C.
  int const offset = C->offset;
  rci_t const cncols = C->ncols;
  int const flags = C->flags;
  word const bitmask = C->low_bitmask;
  // Extend it's columns to the right to include the whole first word, and only the first word.
  C->offset = 0;
  C->ncols = m4ri_radix;
  C->flags = (flags & (mzd_flag_multiple_blocks | mzd_flag_windowed_ownsblocks)) | (mzd_flag_windowed_zerooffset | mzd_flag_windowed_zeroexcess);
  C->low_bitmask = C->high_bitmask = m4ri_ffff;
  // Run Bw and Dw over all (now single word) rows of B and D respectively.
  word* RESTRICT Bw = mzd_first_row(B);
  word* RESTRICT Dw = mzd_first_row(D);
  int Bblock = 0;
  int Dblock = 0;
  int Bcount = mzd_rows_in_block(B, 0);
  int Dcount = mzd_rows_in_block(D, 0);
  // This is true because D->row_offset == 0, and D contains much less columns
  // (well, at least, it's rowstride is less than or equal B's rowstride).
  // It can at most contain more rows in the first block than B.
  assert(Bcount <= Dcount);
  int count = Bcount;
  // Already substract 'count' from Bcount and Dcount.
  Dcount -= Bcount;
  Bcount = 0;
  wi_t const Browstride = B->rowstride;
  word const mask = B->low_bitmask;
  while(1) {
    // Make count even.
    if ((count & 1)) {
      *Dw = *Bw & mask;
      Bw += Browstride;
      Dw += 1;
    }
    // Unroll the loop a factor of two.
    count >>= 1;
    while (count--) {
      // Inner loop. Copy the first word of B to D, setting the extra columns to zero.
      Dw[0] = *Bw & mask;
      Bw += Browstride;
      Dw[1] = *Bw & mask;
      Bw += Browstride;
      Dw += 2;			// D->rowstride == 1
    }
    // Unless we have more than one block, we're done.
    // This is always the case the first time we get here, and almost always true subsequent times.
    if (__M4RI_LIKELY(Bcount == 0)) {
      // This is true if we just processed the last block; optimize for the case
      // of a single block matrix and mark it as likely.
      if (__M4RI_LIKELY((Bcount = mzd_rows_in_block(B, ++Bblock)) <= 0))
	break;
      // Put Bw at the start of the next block.
      Bw = mzd_first_row_next_block(B, Bblock);
    } else {
      // then Dcount == 0. Do the same as above but for D.
      if ((Dcount = mzd_rows_in_block(D, ++Dblock)) <= 0)
	break;
      Dw = mzd_first_row_next_block(D, Dblock);
    }
    count = MIN(Bcount, Dcount);
    Bcount -= count;
    Dcount -= count;
  }
  _mzd_addmul_even (C, A, D, cutoff);
  C->offset = offset;
  C->ncols = cncols;
  C->flags = flags;
  C->low_bitmask = C->high_bitmask = bitmask;
  mzd_free(D);
  return C;
}

mzd_t *_mzd_addmul_weird_weird (mzd_t *C, mzd_t const *A, mzd_t const *B) {
  mzd_t *BT = mzd_init( B->ncols, B->nrows );

  for (rci_t i = 0; i < B->ncols; ++i) {
    word *dstp = BT->rows[i];
    wi_t const ii = (i + B->offset) / m4ri_radix;
    int const is = (i + B->offset) % m4ri_radix;
    word const mask = m4ri_one << is;
    int const ke = (is - A->offset < 0) ? -1 : is - A->offset;
    rci_t k = B->nrows - 1;
    while (k > ke)
    {
      *dstp |= (B->rows[k][ii] & mask) << (k - ke);
      --k;
    }
    while (k >= 0)
    {
      *dstp |= (B->rows[k][ii] & mask) >> (ke - k);
      --k;
    }
  }

  assert(C->offset + C->ncols - 1 < 64);
  word parity[64];
  memset(parity, 0, sizeof(parity));
#ifdef M4RI_WRAPWORD
  word::init_array(parity, 64);
#endif
  for (rci_t i = 0; i < A->nrows; ++i) {
    word *a = A->rows[i];
    word *c = C->rows[i];
    for (rci_t k = 0; k < C->ncols; ++k) {
      word *b = BT->rows[k];
      parity[(k + C->offset)] = (*a) & (*b);
    }
    word par = m4ri_parity64(parity);
    *c ^= par;//m4ri_parity64(parity);
  }
  mzd_free (BT);

  __M4RI_DD_MZD(C);
  return C;
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
