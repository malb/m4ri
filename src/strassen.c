 /*******************************************************************
 *
 *            M4RI: Method of the Four Russians Inversion
 *
 *       Copyright (C) 2008 Martin Albrecht <M.R.Albrecht@rhu.ac.uk>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
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

#include "grayflex.h"
#include "strassen.h"
#include "misc.h"
#include "parity.h"
#define CLOSER(a,b,target) (abs((long)a-(long)target)<abs((long)b-(long)target))
#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#ifdef HAVE_OPENMP
#include <omp.h>
#endif



packedmatrix *_mzd_mul_even(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  size_t a,b,c;
  size_t anr, anc, bnr, bnc;
  
  a = A->nrows;
  b = A->ncols;
  c = B->ncols;

  /* handle case first, where the input matrices are too small already */
  if (CLOSER(A->nrows, A->nrows/2, cutoff) || CLOSER(A->ncols, A->ncols/2, cutoff) || CLOSER(B->ncols, B->ncols/2, cutoff)) {
    /* we copy the matrix first since it is only constant memory
       overhead and improves data locality, if you remove it make sure
       there are no speed regressions */
    /* C = _mzd_mul_m4rm(C, A, B, 0, TRUE); */
    /* return C; */
    packedmatrix *Cbar = mzd_init(C->nrows, C->ncols);
    Cbar = _mzd_mul_m4rm(Cbar, A, B, 0, FALSE);
    mzd_copy(C, Cbar);
    mzd_free(Cbar);
    return C;
  }

  /* adjust cutting numbers to work on words */
  unsigned long mult = 1;
  long width = MIN(MIN(a,b),c);
  while (width > 2*cutoff) {
    width/=2;
    mult*=2;
  }
  a -= a%(RADIX*mult);
  b -= b%(RADIX*mult);
  c -= c%(RADIX*mult);

  anr = ((a/RADIX) >> 1) * RADIX;
  anc = ((b/RADIX) >> 1) * RADIX;
  bnr = anc;
  bnc = ((c/RADIX) >> 1) * RADIX;

  packedmatrix *A00 = mzd_init_window(A,   0,   0,   anr,   anc);
  packedmatrix *A01 = mzd_init_window(A,   0, anc,   anr, 2*anc);
  packedmatrix *A10 = mzd_init_window(A, anr,   0, 2*anr,   anc);
  packedmatrix *A11 = mzd_init_window(A, anr, anc, 2*anr, 2*anc);

  packedmatrix *B00 = mzd_init_window(B,   0,   0,   bnr,   bnc);
  packedmatrix *B01 = mzd_init_window(B,   0, bnc,   bnr, 2*bnc);
  packedmatrix *B10 = mzd_init_window(B, bnr,   0, 2*bnr,   bnc);
  packedmatrix *B11 = mzd_init_window(B, bnr, bnc, 2*bnr, 2*bnc);

  packedmatrix *C00 = mzd_init_window(C,   0,   0,   anr,   bnc);
  packedmatrix *C01 = mzd_init_window(C,   0, bnc,   anr, 2*bnc);
  packedmatrix *C10 = mzd_init_window(C, anr,   0, 2*anr,   bnc);
  packedmatrix *C11 = mzd_init_window(C, anr, bnc, 2*anr, 2*bnc);
  
  /**
   * \note See Jean-Guillaume Dumas, Clement Pernet, Wei Zhou; "Memory
   * efficient scheduling of Strassen-Winograd's matrix multiplication
   * algorithm"; http://arxiv.org/pdf/0707.2347v3 for reference on the
   * used operation scheduling.
   */

  /* change this to mzd_init(anr, MAX(bnc,anc)) to fix the todo below */
  packedmatrix *X0 = mzd_init(anr, anc);
  packedmatrix *X1 = mzd_init(bnr, bnc);
  
  _mzd_add(X0, A00, A10);              /*1    X0 = A00 + A10 */
  _mzd_add(X1, B11, B01);              /*2    X1 = B11 + B01 */
  _mzd_mul_even(C10, X0, X1, cutoff);  /*3   C10 = X0*X1 */

  _mzd_add(X0, A10, A11);              /*4    X0 = A10 + A11 */
  _mzd_add(X1, B01, B00);              /*5    X1 = B01 + B00*/
  _mzd_mul_even(C11, X0, X1, cutoff);  /*6   C11 = X0*X1 */

  _mzd_add(X0, X0, A00);               /*7    X0 = X0 + A00 */
  _mzd_add(X1, X1, B11);               /*8    X1 = B11 + X1 */
  _mzd_mul_even(C01, X0, X1, cutoff);  /*9   C01 = X0*X1 */

  _mzd_add(X0, X0, A01);               /*10   X0 = A01 + X0 */
  _mzd_mul_even(C00, X0, B11, cutoff); /*11  C00 = X0*B11 */

  /**
   * \todo ideally we would use the same X0 throughout the function
   * but some called function doesn't like that and we end up with a
   * wrong result if we use virtual X0 matrices. Ideally, this should
   * be fixed not worked around. The check whether the bug has been
   * fixed, use only one X0 and check if mzd_mul(4096, 3528, 4096,
   * 1024) still returns the correct answer.
   */

  mzd_free(X0);
  X0 = mzd_mul(NULL, A00, B00, cutoff);/*12  X0 = A00*B00*/

  _mzd_add(C01, X0, C01);              /*13  C01 =  X0 + C01 */
  _mzd_add(C10, C01, C10);             /*14  C10 = C01 + C10 */
  _mzd_add(C01, C01, C11);             /*15  C01 = C01 + C11 */
  _mzd_add(C11, C10, C11);             /*16  C11 = C10 + C11 */
  _mzd_add(C01, C01, C00);             /*17  C01 = C01 + C00 */
  _mzd_add(X1, X1, B10);               /*18   X1 = X1 + B10 */
  _mzd_mul_even(C00, A11, X1, cutoff); /*19  C00 = A11*X1 */

  _mzd_add(C10, C10, C00);             /*20  C10 = C10 + C00 */
  _mzd_mul_even(C00, A01, B10, cutoff);/*21  C00 = A01*B10 */

  _mzd_add(C00, C00, X0);              /*22  C00 = X0 + C00 */

  /* deal with rest */
  if (B->ncols > (int)(2*bnc)) {
    packedmatrix *B_last_col = mzd_init_window(B, 0, 2*bnc, A->ncols, B->ncols); 
    packedmatrix *C_last_col = mzd_init_window(C, 0, 2*bnc, A->nrows, C->ncols);
    _mzd_mul_m4rm(C_last_col, A, B_last_col, 0, TRUE);
    mzd_free_window(B_last_col);
    mzd_free_window(C_last_col);
  }
  if (A->nrows > (int)(2*anr)) {
    packedmatrix *A_last_row = mzd_init_window(A, 2*anr, 0, A->nrows, A->ncols);
    packedmatrix *C_last_row = mzd_init_window(C, 2*anr, 0, C->nrows, C->ncols);
    _mzd_mul_m4rm(C_last_row, A_last_row, B, 0, TRUE);
    mzd_free_window(A_last_row);
    mzd_free_window(C_last_row);
  }
  if (A->ncols > (int)(2*anc)) {
    packedmatrix *A_last_col = mzd_init_window(A,     0, 2*anc, 2*anr, A->ncols);
    packedmatrix *B_last_row = mzd_init_window(B, 2*bnr,     0, B->nrows, 2*bnc);
    packedmatrix *C_bulk = mzd_init_window(C, 0, 0, 2*anr, bnc*2);
    mzd_addmul_m4rm(C_bulk, A_last_col, B_last_row, 0);
    mzd_free_window(A_last_col);
    mzd_free_window(B_last_row);
    mzd_free_window(C_bulk);
  }

  /* clean up */
  mzd_free_window(A00); mzd_free_window(A01);
  mzd_free_window(A10); mzd_free_window(A11);

  mzd_free_window(B00); mzd_free_window(B01);
  mzd_free_window(B10); mzd_free_window(B11);

  mzd_free_window(C00); mzd_free_window(C01);
  mzd_free_window(C10); mzd_free_window(C11);
  
  mzd_free(X0);
  mzd_free(X1);

  return C;
}

#ifdef HAVE_OPENMP
packedmatrix *_mzd_mul_mp_even(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  /**
   * \todo: make sure not to overwrite crap after ncols and before width*RADIX
   */
  size_t a,b,c;
  size_t anr, anc, bnr, bnc;
  
  a = A->nrows;
  b = A->ncols;
  c = B->ncols;
  /* handle case first, where the input matrices are too small already */
  if (CLOSER(A->nrows, A->nrows/2, cutoff) || CLOSER(A->ncols, A->ncols/2, cutoff) || CLOSER(B->ncols, B->ncols/2, cutoff)) {
    /* we copy the matrix first since it is only constant memory
       overhead and improves data locality, if you remove it make sure
       there are no speed regressions */
    /* C = _mzd_mul_m4rm(C, A, B, 0, TRUE); */
    packedmatrix *Cbar = mzd_init(C->nrows, C->ncols);
    Cbar = _mzd_mul_m4rm(Cbar, A, B, 0, FALSE);
    mzd_copy(C, Cbar);
    mzd_free(Cbar);
    return C;
  }

  /* adjust cutting numbers to work on words */
  unsigned long mult = 2;
/*   long width = a; */
/*   while (width > 2*cutoff) { */
/*     width/=2; */
/*     mult*=2; */
/*   } */
  a -= a%(RADIX*mult);
  b -= b%(RADIX*mult);
  c -= c%(RADIX*mult);

  anr = ((a/RADIX) >> 1) * RADIX;
  anc = ((b/RADIX) >> 1) * RADIX;
  bnr = anc;
  bnc = ((c/RADIX) >> 1) * RADIX;

  packedmatrix *A00 = mzd_init_window(A,   0,   0,   anr,   anc);
  packedmatrix *A01 = mzd_init_window(A,   0, anc,   anr, 2*anc);
  packedmatrix *A10 = mzd_init_window(A, anr,   0, 2*anr,   anc);
  packedmatrix *A11 = mzd_init_window(A, anr, anc, 2*anr, 2*anc);

  packedmatrix *B00 = mzd_init_window(B,   0,   0,   bnr,   bnc);
  packedmatrix *B01 = mzd_init_window(B,   0, bnc,   bnr, 2*bnc);
  packedmatrix *B10 = mzd_init_window(B, bnr,   0, 2*bnr,   bnc);
  packedmatrix *B11 = mzd_init_window(B, bnr, bnc, 2*bnr, 2*bnc);

  packedmatrix *C00 = mzd_init_window(C,   0,   0,   anr,   bnc);
  packedmatrix *C01 = mzd_init_window(C,   0, bnc,   anr, 2*bnc);
  packedmatrix *C10 = mzd_init_window(C, anr,   0, 2*anr,   bnc);
  packedmatrix *C11 = mzd_init_window(C, anr, bnc, 2*anr, 2*bnc);
  
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
  if (B->ncols > (int)(2*bnc)) {
    packedmatrix *B_last_col = mzd_init_window(B, 0, 2*bnc, A->ncols, B->ncols); 
    packedmatrix *C_last_col = mzd_init_window(C, 0, 2*bnc, A->nrows, C->ncols);
    _mzd_mul_m4rm(C_last_col, A, B_last_col, 0, TRUE);
    mzd_free_window(B_last_col);
    mzd_free_window(C_last_col);
  }
  if (A->nrows > (int)(2*anr)) {
    packedmatrix *A_last_row = mzd_init_window(A, 2*anr, 0, A->nrows, A->ncols);
    packedmatrix *C_last_row = mzd_init_window(C, 2*anr, 0, C->nrows, C->ncols);
    _mzd_mul_m4rm(C_last_row, A_last_row, B, 0, TRUE);
    mzd_free_window(A_last_row);
    mzd_free_window(C_last_row);
  }
  if (A->ncols > (int)(2*anc)) {
    packedmatrix *A_last_col = mzd_init_window(A,     0, 2*anc, 2*anr, A->ncols);
    packedmatrix *B_last_row = mzd_init_window(B, 2*bnr,     0, B->nrows, 2*bnc);
    packedmatrix *C_bulk = mzd_init_window(C, 0, 0, 2*anr, bnc*2);
    mzd_addmul_m4rm(C_bulk, A_last_col, B_last_row, 0);
    mzd_free_window(A_last_col);
    mzd_free_window(B_last_row);
    mzd_free_window(C_bulk);
  }

  /* clean up */
  mzd_free_window(A00); mzd_free_window(A01);
  mzd_free_window(A10); mzd_free_window(A11);

  mzd_free_window(B00); mzd_free_window(B01);
  mzd_free_window(B10); mzd_free_window(B11);

  mzd_free_window(C00); mzd_free_window(C01);
  mzd_free_window(C10); mzd_free_window(C11);
  
  return C;
}
#endif

packedmatrix *mzd_mul(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  if(A->ncols != B->nrows)
    m4ri_die("mzd_mul: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);
  
  if (cutoff < 0)
    m4ri_die("mzd_mul: cutoff must be >= 0.\n");

  if(cutoff == 0) {
    cutoff = STRASSEN_MUL_CUTOFF;
  }

  cutoff = cutoff/RADIX * RADIX;
  if (cutoff < RADIX) {
    cutoff = RADIX;
  };

  if (C == NULL) {
    C = mzd_init(A->nrows, B->ncols);
  } else if (C->nrows != A->nrows || C->ncols != B->ncols){
    m4ri_die("mzd_mul: C (%d x %d) has wrong dimensions, expected (%d x %d)\n",
	     C->nrows, C->ncols, A->nrows, B->ncols);
  }
#ifdef HAVE_OPENMP
  /* this one isn't optimal */
  if (omp_get_max_threads() > 1) {
    C = _mzd_mul_mp_even(C, A, B, cutoff);
  } else {
    C = _mzd_mul_even(C, A, B, cutoff);
  }
#else
  C = _mzd_mul_even(C, A, B, cutoff);
#endif  
  return C;
}

packedmatrix *_mzd_addmul_even(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  /**
   * \todo: make sure not to overwrite crap after ncols and before width*RADIX
   */

  size_t a,b,c;
  size_t anr, anc, bnr, bnc;
  
  a = A->nrows;
  b = A->ncols;
  c = B->ncols;
  /* handle case first, where the input matrices are too small already */
  if (CLOSER(A->nrows, A->nrows/2, cutoff) || CLOSER(A->ncols, A->ncols/2, cutoff) || CLOSER(B->ncols, B->ncols/2, cutoff)) {
    /* we copy the matrix first since it is only constant memory
       overhead and improves data locality, if you remove it make sure
       there are no speed regressions */
    packedmatrix *Cbar = mzd_copy (NULL, C);
    mzd_addmul_m4rm (Cbar, A, B, 0);
    mzd_copy(C, Cbar);
    mzd_free(Cbar);
    return C;
  }

  /* adjust cutting numbers to work on words */
  unsigned long mult = 1;
  long width = MIN(MIN(a,b),c);
  while (width > 2*cutoff) {
    width/=2;
    mult*=2;
  }
  a -= a%(RADIX*mult);
  b -= b%(RADIX*mult);
  c -= c%(RADIX*mult);

  anr = ((a/RADIX) >> 1) * RADIX;
  anc = ((b/RADIX) >> 1) * RADIX;
  bnr = anc;
  bnc = ((c/RADIX) >> 1) * RADIX;

  packedmatrix *A00 = mzd_init_window(A,   0,   0,   anr,   anc);
  packedmatrix *A01 = mzd_init_window(A,   0, anc,   anr, 2*anc);
  packedmatrix *A10 = mzd_init_window(A, anr,   0, 2*anr,   anc);
  packedmatrix *A11 = mzd_init_window(A, anr, anc, 2*anr, 2*anc);

  packedmatrix *B00 = mzd_init_window(B,   0,   0,   bnr,   bnc);
  packedmatrix *B01 = mzd_init_window(B,   0, bnc,   bnr, 2*bnc);
  packedmatrix *B10 = mzd_init_window(B, bnr,   0, 2*bnr,   bnc);
  packedmatrix *B11 = mzd_init_window(B, bnr, bnc, 2*bnr, 2*bnc);

  packedmatrix *C00 = mzd_init_window(C,   0,   0,   anr,   bnc);
  packedmatrix *C01 = mzd_init_window(C,   0, bnc,   anr, 2*bnc);
  packedmatrix *C10 = mzd_init_window(C, anr,   0, 2*anr,   bnc);
  packedmatrix *C11 = mzd_init_window(C, anr, bnc, 2*anr, 2*bnc);

  /**
   * \note See Jean-Guillaume Dumas, Clement Pernet, Wei Zhou; "Memory
   * efficient scheduling of Strassen-Winograd's matrix multiplication
   * algorithm"; http://arxiv.org/pdf/0707.2347v3 for reference on the
   * used operation scheduling.
   */

  packedmatrix *X0 = mzd_init(anr, anc);
  packedmatrix *X1 = mzd_init(bnr, bnc);
  packedmatrix *X2 = mzd_init(anr, bnc);
  
  _mzd_add(X0, A10, A11);                  /* 1  S1 = A21 + A22        X1 */
  _mzd_add(X1, B01, B00);                  /* 2  T1 = B12 - B11        X2 */
  _mzd_mul_even(X2, X0, X1, cutoff);       /* 3  P5 = S1 T1            X3 */
  
  _mzd_add(C11, X2, C11);                  /* 4  C22 = P5 + C22       C22 */
  _mzd_add(C01, X2, C01);                  /* 5  C12 = P5 + C12       C12 */
  _mzd_add(X0, X0, A00);                   /* 6  S2 = S1 - A11         X1 */
  _mzd_add(X1, B11, X1);                   /* 7  T2 = B22 - T1         X2 */
  _mzd_mul_even(X2, A00, B00, cutoff);     /* 8  P1 = A11 B11          X3 */
  
  _mzd_add(C00, X2, C00);                  /* 9  C11 = P1 + C11       C11 */
  _mzd_addmul_even(X2, X0, X1, cutoff);    /* 10 U2 = S2 T2 + P1       X3 */

  _mzd_addmul_even(C00, A01, B10, cutoff); /* 11 U1 = A12 B21 + C11   C11 */
  
  _mzd_add(X0, A01, X0);                   /* 12 S4 = A12 - S2         X1 */
  _mzd_add(X1, X1, B10);                   /* 13 T4 = T2 - B21         X2 */
  _mzd_addmul_even(C01, X0, B11, cutoff);  /* 14 C12 = S4 B22 + C12   C12 */
  
  _mzd_add(C01, X2, C01);                  /* 15 U5 = U2 + C12        C12 */
  _mzd_addmul_even(C10, A11, X1, cutoff);  /* 16 P4 = A22 T4 - C21    C21 */
  
  _mzd_add(X0, A00, A10);                  /* 17 S3 = A11 - A21        X1 */
  _mzd_add(X1, B11, B01);                  /* 18 T3 = B22 - B12        X2 */
  _mzd_addmul_even(X2, X0, X1, cutoff);    /* 19 U3 = S3 T3 + U2       X3 */
  
  _mzd_add(C11, X2, C11);                  /* 20 U7 = U3 + C22        C22 */
  _mzd_add(C10, X2, C10);                  /* 21 U6 = U3 - C21        C21 */

  /* deal with rest */
  if (B->ncols > 2*bnc) {
    packedmatrix *B_last_col = mzd_init_window(B, 0, 2*bnc, A->ncols, B->ncols); 
    packedmatrix *C_last_col = mzd_init_window(C, 0, 2*bnc, A->nrows, C->ncols);
    mzd_addmul_m4rm(C_last_col, A, B_last_col, 0);
    mzd_free_window(B_last_col);
    mzd_free_window(C_last_col);
  }
  if (A->nrows > 2*anr) {
    packedmatrix *A_last_row = mzd_init_window(A, 2*anr, 0, A->nrows, A->ncols);
    packedmatrix *B_bulk = mzd_init_window(B, 0, 0, B->nrows, 2*bnc);
    packedmatrix *C_last_row = mzd_init_window(C, 2*anr, 0, C->nrows, 2*bnc);
    mzd_addmul_m4rm(C_last_row, A_last_row, B_bulk, 0);
    mzd_free_window(A_last_row);
    mzd_free_window(B_bulk);
    mzd_free_window(C_last_row);
  }
  if (A->ncols > 2*anc) {
    packedmatrix *A_last_col = mzd_init_window(A,     0, 2*anc, 2*anr, A->ncols);
    packedmatrix *B_last_row = mzd_init_window(B, 2*bnr,     0, B->nrows, 2*bnc);
    packedmatrix *C_bulk = mzd_init_window(C, 0, 0, 2*anr, bnc*2);
    mzd_addmul_m4rm(C_bulk, A_last_col, B_last_row, 0);
    mzd_free_window(A_last_col);
    mzd_free_window(B_last_row);
    mzd_free_window(C_bulk);
  }

  /* clean up */
  mzd_free_window(A00); mzd_free_window(A01);
  mzd_free_window(A10); mzd_free_window(A11);

  mzd_free_window(B00); mzd_free_window(B01);
  mzd_free_window(B10); mzd_free_window(B11);

  mzd_free_window(C00); mzd_free_window(C01);
  mzd_free_window(C10); mzd_free_window(C11);
  
  mzd_free(X0);
  mzd_free(X1);
  mzd_free(X2);

  return C;
}

packedmatrix *_mzd_addmul(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff){
  /**
   * Assumes that B and C are aligned in the same manner (as in a Schur complement)
   */
  
  if (!A->offset){
    if (!B->offset) /* A even, B even */
      return _mzd_addmul_even (C, A, B, cutoff);
    else {  /* A even, B weird */
      size_t bnc = RADIX - B->offset;
      if (B->ncols <= bnc){
	_mzd_addmul_even_weird  (C,  A, B, cutoff);
      } else {
	packedmatrix * B0 = mzd_init_window (B, 0, 0, B->nrows, bnc);
	packedmatrix * C0 = mzd_init_window (C, 0, 0, C->nrows, bnc);
	packedmatrix * B1 = mzd_init_window (B, 0, bnc, B->nrows, B->ncols);
	packedmatrix * C1 = mzd_init_window (C, 0, bnc, C->nrows, C->ncols);
	_mzd_addmul_even_weird  (C0,  A, B0, cutoff);
	_mzd_addmul_even (C1, A, B1, cutoff);
	mzd_free_window (B0); mzd_free_window (B1);
	mzd_free_window (C0); mzd_free_window (C1);
      }
    }
  } else if (B->offset) { /* A weird, B weird */
    size_t anc = RADIX - A->offset;
    size_t bnc = RADIX - B->offset;
    if (B->ncols <= bnc){
      if (A->ncols <= anc)
	_mzd_addmul_weird_weird (C, A, B, cutoff);
      else {
	packedmatrix * A0  = mzd_init_window (A, 0, 0, A->nrows, anc);
	packedmatrix * A1  = mzd_init_window (A, 0, anc, A->nrows, A->ncols);
	packedmatrix * B0  = mzd_init_window (B, 0, 0, anc, B->ncols);
	packedmatrix * B1  = mzd_init_window (B, anc, 0, B->nrows, B->ncols);
	_mzd_addmul_weird_weird (C, A0, B0, cutoff);
	_mzd_addmul_even_weird  (C, A1, B1, cutoff);
	mzd_free_window (A0);  mzd_free_window (A1);
	mzd_free_window (B0);  mzd_free_window (B1);
      }
    } else if (A->ncols <= anc) {
      packedmatrix * B0 = mzd_init_window (B, 0, 0, B->nrows, bnc);
      packedmatrix * B1 = mzd_init_window (B, 0, bnc, B->nrows, B->ncols);
      packedmatrix * C0 = mzd_init_window (C, 0, 0, C->nrows, bnc);
      packedmatrix * C1 = mzd_init_window (C, 0, bnc, C->nrows, C->ncols);
      _mzd_addmul_weird_weird (C0, A, B0, cutoff);
      _mzd_addmul_weird_even  (C1, A, B1, cutoff);
      mzd_free_window (B0); mzd_free_window (B1);
      mzd_free_window (C0); mzd_free_window (C1);
    } else {
      packedmatrix * A0  = mzd_init_window (A, 0, 0, A->nrows, anc);
      packedmatrix * A1  = mzd_init_window (A, 0, anc, A->nrows, A->ncols);
      packedmatrix * B00 = mzd_init_window (B, 0, 0, anc, bnc);
      packedmatrix * B01 = mzd_init_window (B, 0, bnc, anc, B->ncols);
      packedmatrix * B10 = mzd_init_window (B, anc, 0, B->nrows, bnc);
      packedmatrix * B11 = mzd_init_window (B, anc, bnc, B->nrows, B->ncols);
      packedmatrix * C0 = mzd_init_window (C, 0, 0, C->nrows, bnc);
      packedmatrix * C1 = mzd_init_window (C, 0, bnc, C->nrows, C->ncols);
      
      _mzd_addmul_weird_weird (C0, A0, B00, cutoff);
      _mzd_addmul_even_weird  (C0,  A1, B10, cutoff);
      _mzd_addmul_weird_even  (C1,  A0, B01, cutoff);
      _mzd_addmul_even  (C1,  A1, B11, cutoff);

      mzd_free_window (A0);  mzd_free_window (A1);
      mzd_free_window (C0);  mzd_free_window (C1);
      mzd_free_window (B00); mzd_free_window (B01);
      mzd_free_window (B10); mzd_free_window (B11);
    }
  } else { /* A weird, B even */
    int anc = RADIX - A->offset;
    if (A->ncols <= anc){
      _mzd_addmul_weird_even  (C,  A, B, cutoff);
    } else {
      packedmatrix * A0  = mzd_init_window (A, 0, 0, A->nrows, anc);
      packedmatrix * A1  = mzd_init_window (A, 0, anc, A->nrows, A->ncols);
      packedmatrix * B0  = mzd_init_window (B, 0, 0, anc, B->ncols);
      packedmatrix * B1  = mzd_init_window (B, anc, 0, B->nrows, B->ncols);
      _mzd_addmul_weird_even (C, A0, B0, cutoff);
      _mzd_addmul_even  (C, A1, B1, cutoff);
      mzd_free_window (A0); mzd_free_window (A1);
      mzd_free_window (B0); mzd_free_window (B1);
    }
  }
  return C;
}

packedmatrix *_mzd_addmul_weird_even (packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff){
  packedmatrix * tmp = mzd_init (A->nrows, MIN(RADIX- A->offset, A->ncols));
  for (size_t i=0; i < A->nrows; ++i){
    tmp->values [tmp->rowswap[i]] = (A->values [A->rowswap [i]] << A->offset);
  }
  _mzd_addmul_even (C, tmp, B, cutoff);
  mzd_free(tmp);
  return C;
}

 packedmatrix *_mzd_addmul_even_weird (packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff){
   packedmatrix * tmp = mzd_init (B->nrows, RADIX);
   size_t offset = C->offset;
   size_t cncols = C->ncols;
   C->offset=0;
   C->ncols = RADIX;
   word mask = ((ONE << B->ncols) - 1) << (RADIX-B->offset - B->ncols);
   for (size_t i=0; i < B->nrows; ++i)
     tmp->values [tmp->rowswap[i]] = B->values [B->rowswap [i]] & mask;
   _mzd_addmul_even (C, A, tmp, cutoff);
   C->offset=offset;
   C->ncols = cncols;
   mzd_free (tmp);
   return C;
}

 packedmatrix* _mzd_addmul_weird_weird (packedmatrix* C, packedmatrix* A, packedmatrix *B, int cutoff){
   packedmatrix *BT;
   word* temp;
   BT = mzd_init( B->ncols, B->nrows );
   
   for (int i = 0; i < B->ncols; ++i) {
     temp = BT->values + BT->rowswap[i];
     for (size_t k = 0; k < B->nrows; k++) {
      *temp |= ((word)mzd_read_bit (B, k, i)) << (RADIX-1-k-A->offset);

     }
   }
   
   word parity[64];
   for (size_t i = 0; i < 64; i++) {
     parity[i] = 0;
   }
   for (size_t i = 0; i < A->nrows; ++i) {
     word * a = A->values + A->rowswap[i];
     word * c = C->values + C->rowswap[i];
     for (size_t k=0; k< C->ncols; k++) {
       word *b = BT->values + BT->rowswap[k];
       parity[k+C->offset] = (*a) & (*b);
     }
     word par = parity64(parity);
     *c ^= par;//parity64(parity);
   }
   mzd_free (BT);
   return C;
 }

 packedmatrix *mzd_addmul(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  if(A->ncols != B->nrows)
    m4ri_die("mzd_addmul: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);
  
  if (cutoff < 0)
    m4ri_die("mzd_addmul: cutoff must be >= 0.\n");

  if(cutoff == 0) {
    cutoff = STRASSEN_MUL_CUTOFF;
  }
  
  cutoff = cutoff/RADIX * RADIX;
  if (cutoff < RADIX) {
    cutoff = RADIX;
  };

  if (C == NULL) {
    C = mzd_init(A->nrows, B->ncols);
  } else if (C->nrows != A->nrows || C->ncols != B->ncols){
    m4ri_die("mzd_addmul: C (%d x %d) has wrong dimensions, expected (%d x %d)\n",
	     C->nrows, C->ncols, A->nrows, B->ncols);
  }
  C = _mzd_addmul(C, A, B, cutoff);
  return C;
}


