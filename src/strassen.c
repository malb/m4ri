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


#define CLOSER(a,b,target) (abs((long)a-(long)target)<abs((long)b-(long)target))

packedmatrix *_mzd_mul_strassen_impl_even(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  /**
   * \todo: make sure not to overwrite crap after ncols and before width*RADIX
   */
  unsigned int a,b,c;
  unsigned int anr, anc, bnr, bnc;
  
  a = A->nrows;
  b = A->ncols;
  c = B->ncols;
  /* handle case first, where the input matrices are too small already */
  if (CLOSER(A->nrows, A->nrows/2, cutoff) || CLOSER(A->ncols, A->ncols/2, cutoff) || CLOSER(B->ncols, B->ncols/2, cutoff)) {
    /* we copy the matrix first since it is only constant memory
       overhead and improves data locality, if you remove it make sure
       there are no speed regressions */
    /* C = _mzd_mul_m4rm_impl(C, A, B, 0, TRUE); */
    packedmatrix *Cbar = mzd_init(C->nrows, C->ncols);
    Cbar = _mzd_mul_m4rm_impl(Cbar, A, B, 0, FALSE);
    mzd_copy(C, Cbar);
    mzd_free(Cbar);
    return C;
  }

  /* adjust cutting numbers to work on words */
  unsigned long mult = 1;
  long width = a;
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
  
  _mzd_add_impl(X0, A00, A10);                  /*1    X0 = A00 + A10 */
  _mzd_add_impl(X1, B11, B01);                  /*2    X1 = B11 + B01 */
  _mzd_mul_strassen_impl_even(C10, X0, X1, cutoff);  /*3   C10 = X0*X1 */

  _mzd_add_impl(X0, A10, A11);                  /*4    X0 = A10 + A11 */
  _mzd_add_impl(X1, B01, B00);                  /*5    X1 = B01 + B00*/
  _mzd_mul_strassen_impl_even(C11, X0, X1, cutoff);  /*6   C11 = X0*X1 */

  _mzd_add_impl(X0, X0, A00);                   /*7    X0 = X0 + A00 */
  _mzd_add_impl(X1, X1, B11);                   /*8    X1 = B11 + X1 */
  _mzd_mul_strassen_impl_even(C01, X0, X1, cutoff);  /*9   C01 = X0*X1 */

  _mzd_add_impl(X0, X0, A01);                   /*10   X0 = A01 + X0 */
  _mzd_mul_strassen_impl_even(C00, X0, B11, cutoff); /*11  C00 = X0*B11 */

  /**
   * \todo ideally we would use the same X0 throughout the function
   * but some called function doesn't like that and we end up with a
   * wrong result if we use virtual X0 matrices. Ideally, this should
   * be fixed not worked around. The check whether the bug has been
   * fixed, use only one X0 and check if mzd_mul_strassen(4096, 3528,
   * 4096, 1024) still returns the correct answer.
   */

  mzd_free(X0);
  X0 = mzd_mul(NULL, A00, B00, cutoff);/*12  X0 = A00*B00*/

  _mzd_add_impl(C01, X0, C01);                  /*13  C01 =  X0 + C01 */
  _mzd_add_impl(C10, C01, C10);                 /*14  C10 = C01 + C10 */
  _mzd_add_impl(C01, C01, C11);                 /*15  C01 = C01 + C11 */
  _mzd_add_impl(C11, C10, C11);                 /*16  C11 = C10 + C11 */
  _mzd_add_impl(C01, C01, C00);                 /*17  C01 = C01 + C00 */
  _mzd_add_impl(X1, X1, B10);                   /*18   X1 = X1 + B10 */
  _mzd_mul_strassen_impl_even(C00, A11, X1, cutoff); /*19  C00 = A11*X1 */

  _mzd_add_impl(C10, C10, C00);                 /*20  C10 = C10 + C00 */
  _mzd_mul_strassen_impl_even(C00, A01, B10, cutoff);/*21  C00 = A01*B10 */

  _mzd_add_impl(C00, C00, X0);                  /*22  C00 = X0 + C00 */

  /* deal with rest */
  if (B->ncols > (int)(2*bnc)) {
    packedmatrix *B_last_col = mzd_init_window(B, 0, 2*bnc, A->ncols, B->ncols); 
    packedmatrix *C_last_col = mzd_init_window(C, 0, 2*bnc, A->nrows, C->ncols);
    _mzd_mul_m4rm_impl(C_last_col, A, B_last_col, 0, TRUE);
    mzd_free_window(B_last_col);
    mzd_free_window(C_last_col);
  }
  if (A->nrows > (int)(2*anr)) {
    packedmatrix *A_last_row = mzd_init_window(A, 2*anr, 0, A->nrows, A->ncols);
    packedmatrix *C_last_row = mzd_init_window(C, 2*anr, 0, C->nrows, C->ncols);
    _mzd_mul_m4rm_impl(C_last_row, A_last_row, B, 0, TRUE);
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
packedmatrix *_mzd_mul_strassen_mp_impl_even(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  int a,b,c;
  int anr, anc, bnr, bnc;
  
  a = A->nrows;
  b = A->ncols;
  c = B->ncols;

  /* handle case first, where the input matrices are too small already */
  if (A->nrows <= cutoff || A->ncols <= cutoff || B->ncols <= cutoff) {
    C = _mzd_mul_m4rm_impl(C, A, B, 0, TRUE);
    return C;
  }

  /* adjust cutting numbers to work on words */
  a -= (a%RADIX);
  b -= (b%RADIX);
  c -= (c%RADIX);

  anr = (a/RADIX >> 1) * RADIX;
  anc = (b/RADIX >> 1) * RADIX;
  bnr = anc;
  bnc = (c/RADIX >> 1) * RADIX;

  packedmatrix *A00 = mzd_init_window(A,   0,   0,   anr,   anc);
  packedmatrix *A01 = mzd_init_window(A,   0, anc,   anr, 2*anc);
  packedmatrix *A10 = mzd_init_window(A, anr,   0, 2*anr,   anc);
  packedmatrix *A11 = mzd_init_window(A, anr, anc, 2*anr, 2*anc);

  packedmatrix *B00 = mzd_init_window(B,   0,   0,   bnr,   bnc);
  packedmatrix *B01 = mzd_init_window(B,   0, bnc,   bnr, 2*bnc);
  packedmatrix *B10 = mzd_init_window(B, bnr,   0, 2*bnr,   bnc);
  packedmatrix *B11 = mzd_init_window(B, bnr, bnc, 2*bnr, 2*bnc);

  packedmatrix *U0 = mzd_init_window(C,   0,   0,   anr,   bnc);
  packedmatrix *U6 = mzd_init_window(C,   0, bnc,   anr, 2*bnc);
  packedmatrix *U3 = mzd_init_window(C, anr,   0, 2*anr,   bnc);
  packedmatrix *U4 = mzd_init_window(C, anr, bnc, 2*anr, 2*bnc);

  /**
   * \todo this is way too much memory
   */
  packedmatrix *S0 = mzd_init(anr, anc);
  packedmatrix *S1 = mzd_init(anr, anc);
  packedmatrix *S2 = mzd_init(anr, anc);
  packedmatrix *S3 = mzd_init(anr, anc);

  packedmatrix *T0 = mzd_init(bnr, bnc);
  packedmatrix *T1 = mzd_init(bnr, bnc);
  packedmatrix *T2 = mzd_init(bnr, bnc);
  packedmatrix *T3 = mzd_init(bnr, bnc);
  
  packedmatrix *Q0 = mzd_init(anr, bnc);
  packedmatrix *Q1 = mzd_init(anr, bnc);
  packedmatrix *Q2 = mzd_init(anr, bnc);
  packedmatrix *Q3 = mzd_init(anr, bnc);

#pragma omp parallel sections
  {
#pragma omp section
    {
      _mzd_add_impl(S0, A10, A11);
      _mzd_add_impl(S1,  S0, A00);
      _mzd_add_impl(S2, A00, A10);
      _mzd_add_impl(S3, A01,  S1);
    }
#pragma omp section 
    {
      _mzd_add_impl(T0, B01, B00);
      _mzd_add_impl(T1, B11,  T0);
      _mzd_add_impl(T2, B11, B01);
      _mzd_add_impl(T3, B10,  T1);
    }
  }
  
#pragma omp parallel sections
  {
#pragma omp section
    {
      _mzd_mul_strassen_mp_impl_even(Q0, A00, B00, cutoff); /* now Q0 holds P0 */
    }
#pragma omp section 
    {
      _mzd_mul_strassen_mp_impl_even(Q1, A01, B10, cutoff); /* now Q1 holds P1 */
    }  
  }
  _mzd_add_impl(U0, Q0, Q1); /* now U0 is correct */
  
#pragma omp parallel sections
  {
#pragma omp section
    {
      packedmatrix *S1T1 = mzd_mul_strassen(NULL, S1, T1, cutoff);
      _mzd_add_impl(Q0, Q0, S1T1); /* now Q0 holds U1 */
      mzd_free(S1T1);
    }    
#pragma omp section
    {
      _mzd_mul_strassen_mp_impl_even(Q1, S2, T2, cutoff); /* now Q1 holds P4 */
    }
  }  
  _mzd_add_impl(Q1, Q1, Q0); /* now Q1 holds U2 */

#pragma omp parallel sections
  {
#pragma omp section
    {
      _mzd_mul_strassen_mp_impl_even(Q2, A11, T3, cutoff); /* now Q2 holds P6 */
      _mzd_add_impl(U3, Q1, Q2); /* now U3 is correct */
    }  
#pragma omp section
    {
      _mzd_mul_strassen_mp_impl_even(Q3, S0, T0, cutoff); /* now Q3 holds P2 */
      _mzd_add_impl(U4, Q1, Q3); /* now U4 is correct */
    }
  }  
  _mzd_add_impl(Q0, Q0, Q3); /* now Q0 holds U5 */
  _mzd_mul_strassen_mp_impl_even(Q2, S3, B11, cutoff); /* now Q2 holds P5 */
  _mzd_add_impl(U6, Q0, Q2); /* now U6 is correct */

  /* deal with rest */
  if (B->ncols > 2*bnc) {
    packedmatrix *B_last_col = mzd_init_window(B, 0, 2*bnc, A->ncols, B->ncols); 
    packedmatrix *C_last_col = mzd_init_window(C, 0, 2*bnc, A->nrows, C->ncols);
    _mzd_mul_m4rm_impl(C_last_col, A, B_last_col, 0, TRUE);
    mzd_free_window(B_last_col);
    mzd_free_window(C_last_col);

  }
  if (A->nrows > 2*anr) {
    packedmatrix *A_last_row = mzd_init_window(A, 2*anr, 0, A->nrows, A->ncols);
    packedmatrix *C_last_row = mzd_init_window(C, 2*anr, 0, C->nrows, C->ncols);
    _mzd_mul_m4rm_impl(C_last_row, A_last_row, B, 0, TRUE);
    mzd_free_window(A_last_row);
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

  mzd_free_window(U0); mzd_free_window(U6);
  mzd_free_window(U3); mzd_free_window(U4);
  
  mzd_free(S0); mzd_free(S1);
  mzd_free(S2); mzd_free(S3);

  mzd_free(T0); mzd_free(T1);
  mzd_free(T2); mzd_free(T3);

  mzd_free(Q0); mzd_free(Q1);
  mzd_free(Q2); mzd_free(Q3);

  return C;
}
#endif

packedmatrix *mzd_mul(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  if(A->ncols != B->nrows)
    m4ri_die("mzd_mul_strassen: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);
  
  if (cutoff < 0)
    m4ri_die("mzd_mul_strassen: cutoff must be > 0.\n");

  if(cutoff == 0) {
    cutoff = STRASSEN_MUL_CUTOFF;
  }

  cutoff = cutoff/RADIX * RADIX;
  if (cutoff == 0) {
    cutoff = RADIX;
  };

  if (C == NULL) {
    C = mzd_init(A->nrows, B->ncols);
  } else if (C->nrows != A->nrows || C->ncols != B->ncols){
    m4ri_die("mzd_mul_strassen: C (%d x %d) has wrong dimensions, expected (%d x %d)\n",
	     C->nrows, C->ncols, A->nrows, B->ncols);
  }
#ifdef HAVE_OPENMP
  /* this one isn't optimal */
  return _mzd_mul_strassen_mp_impl_even(C, A, B, cutoff);
#else
  return _mzd_mul_strassen_impl_even(C, A, B, cutoff);
#endif  
}

packedmatrix *_mzd_addmul_strassen_impl_even(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  /**
   * \todo: make sure not to overwrite crap after ncols and before width*RADIX
   */

  unsigned int a,b,c;
  unsigned int anr, anc, bnr, bnc;
  
  a = A->nrows;
  b = A->ncols;
  c = B->ncols;
  /* handle case first, where the input matrices are too small already */
  if (CLOSER(A->nrows, A->nrows/2, cutoff) || CLOSER(A->ncols, A->ncols/2, cutoff) || CLOSER(B->ncols, B->ncols/2, cutoff)) {
    /* we copy the matrix first since it is only constant memory
       overhead and improves data locality, if you remove it make sure
       there are no speed regressions */
    packedmatrix *Cbar = mzd_copy(NULL, C);
    Cbar = mzd_addmul_m4rm(Cbar, A, B, 0);
    mzd_copy(C, Cbar);
    mzd_free(Cbar);
    return C;
  }

  /* adjust cutting numbers to work on words */
  unsigned long mult = 1;
  long width = a;
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
  
  _mzd_add_impl(X0, A10, A11);                      /* 1  S1 = A21 + A22        X1 */
  _mzd_add_impl(X1, B01, B00);                      /* 2  T1 = B12 - B11        X2 */
  _mzd_mul_strassen_impl_even(X2, X0, X1, cutoff);       /* 3  P5 = S1 T1            X3 */
  
  _mzd_add_impl(C11, X2, C11);                      /* 4  C22 = P5 + C22       C22 */
  _mzd_add_impl(C01, X2, C01);                      /* 5  C12 = P5 + C12       C12 */
  _mzd_add_impl(X0, X0, A00);                       /* 6  S2 = S1 - A11         X1 */
  _mzd_add_impl(X1, B11, X1);                       /* 7  T2 = B22 - T1         X2 */
  _mzd_mul_strassen_impl_even(X2, A00, B00, cutoff);     /* 8  P1 = A11 B11          X3 */
  
  _mzd_add_impl(C00, X2, C00);                      /* 9  C11 = P1 + C11       C11 */
  _mzd_addmul_strassen_impl_even(X2, X0, X1, cutoff);    /* 10 U2 = S2 T2 + P1       X3 */

  _mzd_addmul_strassen_impl_even(C00, A01, B10, cutoff); /* 11 U1 = A12 B21 + C11   C11 */
  
  _mzd_add_impl(X0, A01, X0);                       /* 12 S4 = A12 - S2         X1 */
  _mzd_add_impl(X1, X1, B10);                       /* 13 T4 = T2 - B21         X2 */
  _mzd_addmul_strassen_impl_even(C01, X0, B11, cutoff);  /* 14 C12 = S4 B22 + C12   C12 */
  
  _mzd_add_impl(C01, X2, C01);                      /* 15 U5 = U2 + C12        C12 */
  _mzd_addmul_strassen_impl_even(C10, A11, X1, cutoff);  /* 16 P4 = A22 T4 - C21    C21 */
  
  _mzd_add_impl(X0, A00, A10);                      /* 17 S3 = A11 - A21        X1 */
  _mzd_add_impl(X1, B11, B01);                      /* 18 T3 = B22 - B12        X2 */
  _mzd_addmul_strassen_impl_even(X2, X0, X1, cutoff);    /* 19 U3 = S3 T3 + U2       X3 */
  
  _mzd_add_impl(C11, X2, C11);                      /* 20 U7 = U3 + C22        C22 */
  _mzd_add_impl(C10, X2, C10);                      /* 21 U6 = U3 - C21        C21 */

  /* deal with rest */
  if (B->ncols > (int)(2*bnc)) {
    packedmatrix *B_last_col = mzd_init_window(B, 0, 2*bnc, A->ncols, B->ncols); 
    packedmatrix *C_last_col = mzd_init_window(C, 0, 2*bnc, A->nrows, C->ncols);
    mzd_addmul_m4rm(C_last_col, A, B_last_col, 0);
    mzd_free_window(B_last_col);
    mzd_free_window(C_last_col);
  }
  if (A->nrows > (int)(2*anr)) {
    packedmatrix *A_last_row = mzd_init_window(A, 2*anr, 0, A->nrows, A->ncols);
    packedmatrix *B_bulk = mzd_init_window(B, 0, 0, B->nrows, 2*bnc);
    packedmatrix *C_last_row = mzd_init_window(C, 2*anr, 0, C->nrows, 2*bnc);//C->ncols);
    mzd_addmul_m4rm(C_last_row, A_last_row, B_bulk, 0);
    mzd_free_window(A_last_row);
    mzd_free_window(B_bulk);
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

#define _mzd_addmul_strassen_impl _mzd_addmul_strassen_impl_even

packedmatrix *mzd_addmul(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  if(A->ncols != B->nrows)
    m4ri_die("mzd_addmul_strassen: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);


  
  if (cutoff < 0)
    m4ri_die("mzd_addmul_strassen: cutoff must be >= 0.\n");

  if(cutoff == 0) {
    cutoff = STRASSEN_MUL_CUTOFF;
  }
  
  cutoff = cutoff/RADIX * RADIX;
  if (cutoff == 0) {
    cutoff = RADIX;
  };

  if (C == NULL) {
    C = mzd_init(A->nrows, B->ncols);
  } else if (C->nrows != A->nrows || C->ncols != B->ncols){
    m4ri_die("mzd_addmul_strassen: C (%d x %d) has wrong dimensions, expected (%d x %d)\n",
	     C->nrows, C->ncols, A->nrows, B->ncols);
  }
  return _mzd_addmul_strassen_impl(C, A, B, cutoff);
}


