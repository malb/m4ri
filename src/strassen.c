 /*******************************************************************
 *
 *            M4RI: Method of the Four Russians Inversion
 *
 *       Copyright (C) 2007, 2008 Gregory Bard <bard@fordham.edu>
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


packedmatrix *_mzd_mul_strassen_impl(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  unsigned int a,b,c;
  unsigned int anr, anc, bnr, bnc;
  
  a = A->nrows;
  b = A->ncols;
  c = B->ncols;
  /* handle case first, where the input matrices are too small already */
  if (A->nrows <= cutoff || A->ncols <= cutoff || B->ncols <= cutoff) {
    C = _mzd_mul_m4rm_impl(C, A, B, 0, NULL, NULL, TRUE);
    return C;
  }

  /* adjust cutting numbers to work on words */
  a -= a%RADIX;
  b -= b%RADIX;
  c -= c%RADIX;

  anr = ((a/RADIX >> 1) * RADIX);
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
  packedmatrix *X0 = mzd_init(anr, MAX(bnc,anc));
  packedmatrix *X1 = mzd_init(anc, bnc);
  
  _mzd_add_impl(X0, A00, A10);                  /*  X0 = A00 + A10 */
  _mzd_add_impl(X1, B11, B01);                  /*  X1 = B11 + B01 */
  _mzd_mul_strassen_impl(C10, X0, X1, cutoff);  /* C10 = X0*X1 */

  _mzd_add_impl(X0, A10, A11);                  /*  X0 = A10 + A11 */
  _mzd_add_impl(X1, B01, B00);                  /*  X1 = B01 + B00*/
  _mzd_mul_strassen_impl(C11, X0, X1, cutoff);  /* C11 = X0*X1 */

  _mzd_add_impl(X0, X0, A00);                   /*  X0 = X0 + A00 */
  _mzd_add_impl(X1, X1, B11);                   /*  X1 = B11 + X1 */
  _mzd_mul_strassen_impl(C01, X0, X1, cutoff);  /* C01 = X0*X1 */

  _mzd_add_impl(X0, X0, A01);                   /*  X0 = A01 + X0 */
  _mzd_mul_strassen_impl(C00, X0, B11, cutoff); /* C00 = X0*B11 */

  _mzd_mul_strassen_impl(X0, A00, B00, cutoff); /* X0 = A00*B00*/

  _mzd_add_impl(C01, X0, C01);                  /* C01 =  X0 + C01 */
  _mzd_add_impl(C10, C01, C10);                 /* C10 = C01 + C10 */
  _mzd_add_impl(C01, C01, C11);                 /* C01 = C01 + C11 */
  _mzd_add_impl(C11, C10, C11);                 /* C11 = C10 + C11 */
  _mzd_add_impl(C01, C01, C00);                 /* C01 = C01 + C00 */
  _mzd_add_impl(X1, X1, B10);                   /*  X1 = X1 + B10 */
  _mzd_mul_strassen_impl(C00, A11, X1, cutoff); /* C00 = A11*X1 */

  _mzd_add_impl(C10, C10, C00);                 /* C10 = C10 + C00 */
  _mzd_mul_strassen_impl(C00, A01, B10, cutoff);/* C00 = A01*B10 */

  _mzd_add_impl(C00, C00, X0);                  /* C00 = X0 + C00 */

  /* deal with rest */
  if (B->ncols > 2*bnc) {
    packedmatrix *B_last_col = mzd_init_window(B, 0, 2*bnc, A->ncols, B->ncols); 
    packedmatrix *C_last_col = mzd_init_window(C, 0, 2*bnc, A->nrows, C->ncols);
    _mzd_mul_m4rm_impl(C_last_col, A, B_last_col, 0, NULL, NULL, TRUE);
    mzd_free_window(B_last_col);
    mzd_free_window(C_last_col);

  }
  if (A->nrows > 2*anr) {
    packedmatrix *A_last_row = mzd_init_window(A, 2*anr, 0, A->nrows, A->ncols);
    packedmatrix *C_last_row = mzd_init_window(C, 2*anr, 0, C->nrows, C->ncols);
    _mzd_mul_m4rm_impl(C_last_row, A_last_row, B, 0, NULL, NULL, TRUE);
    mzd_free_window(A_last_row);
    mzd_free_window(C_last_row);

  }
  if (A->ncols > 2*anc) {
    packedmatrix *A_last_col = mzd_init_window(A,     0, 2*anc, 2*anr, A->ncols);
    packedmatrix *B_last_row = mzd_init_window(B, 2*bnr,     0, B->nrows, 2*bnc);
    packedmatrix *C_bulk = mzd_init_window(C, 0, 0, 2*anr, bnc*2);
    mzd_addmul_m4rm(C_bulk, A_last_col, B_last_row, 0, NULL, NULL);
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

packedmatrix *_mzd_mul_strassen_impl_old_sched(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  int a,b,c;
  int k;
  int anr, anc, bnr, bnc;
  
  a = A->nrows;
  b = A->ncols;
  c = B->ncols;

  /* handle case first, where the input matrices are too small already */
  if (A->nrows <= cutoff || A->ncols <= cutoff || B->ncols <= cutoff) {
    C = _mzd_mul_m4rm_impl(C, A, B, 0, NULL, NULL, TRUE);
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

  packedmatrix *tmp = mzd_init(7*anr + 4*anc, MAX(anc,bnc));

  int start_row = 0;
  packedmatrix *S0 = mzd_init_window(tmp, start_row, 0, start_row + anr, anc);

  start_row += anr;
  packedmatrix *S1 = mzd_init_window(tmp, start_row, 0, start_row + anr, anc);

  start_row += anr;
  packedmatrix *S2 = mzd_init_window(tmp, start_row, 0, start_row + anr, anc);
  
  start_row += anr;
  packedmatrix *S3 = mzd_init_window(tmp, start_row, 0, start_row + anr, anc);

  start_row += anr;
  packedmatrix *T0 = mzd_init_window(tmp, start_row, 0, start_row + anc, bnc);

  start_row += anc;
  packedmatrix *T1 = mzd_init_window(tmp, start_row, 0, start_row + anc, bnc);
  
  start_row += anc;
  packedmatrix *T2 = mzd_init_window(tmp, start_row, 0, start_row + anc, bnc);

  start_row += anc;
  packedmatrix *T3 = mzd_init_window(tmp, start_row, 0, start_row + anc, bnc);
  
  start_row += anc;
  packedmatrix *Q0 = mzd_init_window(tmp, start_row, 0, start_row + anr, bnc);

  start_row += anr;
  packedmatrix *Q1 = mzd_init_window(tmp, start_row, 0, start_row + anr, bnc);
  
  start_row += anr;
  packedmatrix *Q2 = mzd_init_window(tmp, start_row, 0, start_row + anr, bnc);

  _mzd_add_impl(S0, A10, A11);
  _mzd_add_impl(S1,  S0, A00);
  _mzd_add_impl(S2, A00, A10);
  _mzd_add_impl(S3, A01,  S1);

  _mzd_add_impl(T0, B01, B00);
  _mzd_add_impl(T1, B11,  T0);
  _mzd_add_impl(T2, B11, B01);
  _mzd_add_impl(T3, B10,  T1);

  _mzd_mul_strassen_impl(Q0, A00, B00, cutoff); /* now Q0 holds P0 */
  _mzd_mul_strassen_impl(Q1, A01, B10, cutoff); /* now Q1 holds P1 */
  
  _mzd_add_impl(U0, Q0, Q1); /* now U0 is correct */
  
  packedmatrix *S1T1 = mzd_mul_strassen(NULL, S1, T1, cutoff);
  _mzd_add_impl(Q0, Q0, S1T1); /* now Q0 holds U1 */
  mzd_free(S1T1);
    
  _mzd_mul_strassen_impl(Q1, S2, T2, cutoff); /* now Q1 holds P4 */
  
  _mzd_add_impl(Q1, Q1, Q0); /* now Q1 holds U2 */
  _mzd_mul_strassen_impl(Q2, A11, T3, cutoff); /* now Q2 holds P6 */
  _mzd_add_impl(U3, Q1, Q2); /* now U3 is correct */
  
  _mzd_mul_strassen_impl(Q2, S0, T0, cutoff); /* now Q2 holds P2 */
  _mzd_add_impl(U4, Q2, Q1); /* now U4 is correct */
  
  _mzd_add_impl(Q0, Q0, Q2); /* now Q0 holds U5 */
  _mzd_mul_strassen_impl(Q2, S3, B11, cutoff); /* now Q2 holds P5 */
  _mzd_add_impl(U6, Q0, Q2); /* now U6 is correct */

  /* deal with rest */
  if (B->ncols > 2*bnc) {
    packedmatrix *B_last_col = mzd_init_window(B, 0, 2*bnc, A->ncols, B->ncols); 
    packedmatrix *C_last_col = mzd_init_window(C, 0, 2*bnc, A->nrows, C->ncols);
    _mzd_mul_m4rm_impl(C_last_col, A, B_last_col, 0, NULL, NULL, TRUE);
    mzd_free_window(B_last_col);
    mzd_free_window(C_last_col);

  }
  if (A->nrows > 2*anr) {
    packedmatrix *A_last_row = mzd_init_window(A, 2*anr, 0, A->nrows, A->ncols);
    packedmatrix *C_last_row = mzd_init_window(C, 2*anr, 0, C->nrows, C->ncols);
    _mzd_mul_m4rm_impl(C_last_row, A_last_row, B, 0, NULL, NULL, TRUE);
    mzd_free_window(A_last_row);
    mzd_free_window(C_last_row);

  }
  if (A->ncols > 2*anc) {
    packedmatrix *A_last_col = mzd_init_window(A,     0, 2*anc, 2*anr, A->ncols);
    packedmatrix *B_last_row = mzd_init_window(B, 2*bnr,     0, B->nrows, 2*bnc);
    packedmatrix *C_bulk = mzd_init_window(C, 0, 0, 2*anr, bnc*2);
    mzd_addmul_m4rm(C_bulk, A_last_col, B_last_row, 0, NULL, NULL);
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
  
  mzd_free_window(S0); mzd_free_window(S1);
  mzd_free_window(S2); mzd_free_window(S3);

  mzd_free_window(T0); mzd_free_window(T1);
  mzd_free_window(T2); mzd_free_window(T3);

  mzd_free_window(Q0); mzd_free_window(Q1);
  mzd_free_window(Q2);

  mzd_free(tmp);

  return C;
}

packedmatrix *mzd_mul_strassen(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  if(A->ncols != B->nrows)
    m4ri_die("mzd_mul_strassen: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);
  
  if (cutoff <= 0)
    m4ri_die("mzd_mul_strassen: cutoff must be > 0.\n");
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

  return _mzd_mul_strassen_impl(C, A, B, cutoff);
}
