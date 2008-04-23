/**
 * \file packedmatrix.h
 * \brief Dense matrices over GF(2) represented as a bit field.
 *
 * \author Gregory Bard <bard@fordham.edu>
 * \author Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
 */

#ifndef PACKEDMATRIX_H
#define PACKEDMATRIX_H
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

#include "misc.h"
#include <stdio.h>

#define SAFECHAR (int)(RADIX*1.3)

/**
 * \brief Dense matrices over GF(2). 
 * 
 * The most fundamental data type in this library.
 */

typedef struct {
  /**
   * Contains the actual values packed into words of size RADIX.
   */

  word *values;

  /**
   * Number of rows.
   */

  int nrows;

  /**
   * Number of columns.
   */

  int ncols;

  /**
   * width = ceil(nrows/RADIX)
   */
  int width; 

  /**
   * Offsets to each row, so e.g. the first word of the \f$i\f$-th row
   * is m->values[m->rowswap[i]]
   */

  int *rowswap;

} packedmatrix;

/**
 * Create a new matrix of dimension r x c.
 *
 * Use mzd_free to kill it.
 *
 * \param r number of rows
 * \param c number of columns
 *
 */

packedmatrix *mzd_init(int r, int c);

/**
 * Free a matrix created with mzd_init.
 * 
 * \param condemned matrix
 */

void mzd_free( packedmatrix *condemned );

/**
 * Create a window/view into the matrix m.
 *
 * A matrix window for m is a meta structure on the matrix m. It is
 * setup to point into the matrix so m MUST NOT be freed while the
 * matrix window is used.
 *
 * This function puts restrictions on the provided parameters which
 * are not enforced currently.
 *
 *  - lowc must be divisible by RADIX
 *  - highc must be divisible by RADIX
 *  - all parameters must be within range for m
 *
 * Use mzd_free_free to free the window.
 *
 * \param m a matrix
 * \param lowr starting row (inclusive)
 * \param lowc starting column (inclusive)
 * \param highr end row (exclusive)
 * \param highc end column (exclusive)
 *
 */

packedmatrix *mzd_init_window(packedmatrix *m, int lowr, int lowc, int highr, int highc);

/**
 * Frees a matrix window created with mzd_init_window.
 * 
 * \param condemned matrix
 */

void mzd_free_window( packedmatrix *condemned);

/**
 * Swap the two rows rowa and rowb.
 * 
 * \param m Matrix
 * \param rowa Row index
 * \param rowb Row index
 */
 
static inline void mzd_row_swap( packedmatrix *m, int rowa, int rowb ) {
  int temp=m->rowswap[rowa];
  m->rowswap[rowa]=m->rowswap[rowb];
  m->rowswap[rowb]=temp;
}

/**
 * Read the bit at position m[row,col]
 * 
 * \param m Matrix
 * \param row Row index
 * \param col Column index
 */

static inline BIT mzd_read_bit( packedmatrix *m, int row, int col ) {
  return GET_BIT(m->values[ m->rowswap[row] + col/RADIX ], col%RADIX);
}

/**
 * Write the bit value to position m[row,col]
 * 
 * \param m Matrix
 * \param row Row index
 * \param col Column index
 * \param value Either 0 or 1 
 */

static inline void mzd_write_bit( packedmatrix *m, int row, int col, BIT value) {
  if (value==1)
    SET_BIT(m->values[ m->rowswap[row] + col/RADIX ], col % RADIX);
  else
    CLR_BIT(m->values[ m->rowswap[row] + col/RADIX ], col % RADIX);
}

/**
 * \note Keep in mind that the row, col refer to a row and column (of
 *  bits), and you can address the block by any of the RADIX (usually
 *  64) & \f$A_{i,j}\f$ there.
 */

static inline void mzd_xor_block( packedmatrix *m, int row, int col, word value) {
  int block=col/RADIX;
  int truerow=m->rowswap[row];

  word *entry=m->values + block + truerow;
  *entry ^= value;
}

/**
 * \note Keep in mind that the row, col refer to a row and column (of
 * bits), and you can address the block by any of the RADIX (usually
 * 64) \f$A_{i,j}\f$ there.
 */

static inline void mzd_write_block( packedmatrix *m, int row, int col, word value) {
  m->values[ m->rowswap[row] + col/RADIX ] = value;
}

/**
 * \note Keep in mind that the row, col refer to a row and column (of
 * bits), and you can address the block by any of the RADIX (usually
 * 64) \f$A_{i,j}\f$ there.
 */

static inline word mzd_read_block( packedmatrix *m, int row, int col ) {
  return m->values[ m->rowswap[row] + col/RADIX ];
}

/**
 * Print the matrix to stdout. The output will contain colons between
 * every 4-th column.
 *
 * \param m Matrix
 */

void mzd_print_matrix( packedmatrix *m );

/**
 * Print the matrix to stdout.
 *
 * \param m Matrix
 */

void mzd_print_matrix_tight( packedmatrix *m );

/**
 * Add the rows sourcerow and destrow and stores the total in the row
 * destrow, but only begins at the column coloffset.
 *
 * \param m Matrix
 * \param sourcerow Index of source row
 * \param destrow Index of target row
 * \param coloffset Column offset
 *
 * \note this can be done much faster with mzd_combine.
 */

void mzd_row_add_offset( packedmatrix *m, int sourcerow, int destrow, 
			 int coloffset );

/**
 * Clears the given row, but only begins at the column coloffset.
 *
 * \param m Matrix
 * \param row Index of row
 * \param coloffset Column offset
 *
 */

void mzd_row_clear_offset(packedmatrix *m, int row, int coloffset);

/**
 * Add the rows sourcerow and destrow and stores the total in the row
 * destrow.
 *
 * \param m Matrix
 * \param sourcerow Index of source row
 * \param destrow Index of target row
 *
 * \note this can be done much faster with mzd_combine.
 */

void mzd_row_add( packedmatrix *m, int sourcerow, int destrow);

/**
 * Transpose a matrix.
 *
 * This is not efficient, but it is quadratic time, so who cares?
 * Efficient, would be to use the fact that:
 *
\verbatim
   [ A B ]T    [AT CT]
   [ C D ]  =  [BT DT] 
 \endverbatim 
 * and thus rearrange the blocks recursively. 
 *
 * \param newmatrix return matrix (may be NULL)
 * \param data matrix
 */

packedmatrix *mzd_transpose(packedmatrix *newmatrix, packedmatrix *data );

/**
 * Naive cubic matrix multiplication, i.e. compute \f$C\f$ such that
 * \f$C == AB\f$. Where B has to be provided transposed as bT.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A
 * \param bT Input matrix \f$B^T\f$
 *
 * \note Normally, if you will multiply several times by b, it is
 * smarter to calculate bT yourself, and keep it, and then use the
 * function called matrixTimesMatrixTranspose
 */

packedmatrix *mzd_mul_naiv_t(packedmatrix *C, packedmatrix *A, packedmatrix *bT);
  
/**
 * Naive cubic matrix multiplication, i.e. compute \f$C\f$ such that
 * \f$C == AB\f$.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A
 * \param B Input matrix B
 *
 * \note Normally, if you will multiply several times by b, it is
 * smarter to calculate bT yourself, and keep it, and then use the
 * function called matrixTimesMatrixTranspose
 */
packedmatrix *mzd_mul_naiv(packedmatrix *C, packedmatrix *A, packedmatrix *B);

/**
 * Fill matrix a with uniformly distributed bits.
 *
 * \param a Matrix
 */

void mzd_randomize( packedmatrix *a );

/**
 * Set the matrix a to the value equivalent to the integer value
 * provided. Specifically, this function does nothing if \f$value\%2 ==
 * 0\f$ and returns the identity matrix if \f$value\%2 == 1\f$.
 *
 * If the matrix is not square then the largest possible square
 * submatrix is set to the identity matrix.
 *
 * \param a Matrix
 * \param value Either 0 or 1
 */

void mzd_set_ui(packedmatrix *a, unsigned value);

/**
 * This will do Gaussian elimination on the matrix m but will start
 * not at column 0 necc but at column "startcol". If full=NO, then it
 * will do triangular style elimination, and if full=TRUE, it will do
 * Gauss-Jordan style, or full elimination.
 */

int mzd_gauss_delayed(packedmatrix *m, int startcol, int full);

/**
 * This will do Gaussian elimination on the matrix m.  If  full =
 *  FALSE, then it will do triangular style elimination, and if 
 *  full = TRUE, it will do Gauss-Jordan style, or full elimination.
 */
int mzd_reduce_naiv(packedmatrix *m, int full);

/**
 * Return TRUE if  a ==  b.
 *
 * \param a Matrix.
 * \param b Matrix.
 */

BIT mzd_equal( packedmatrix *a, packedmatrix *b );

/**
 * Return -1,0,1 if if  a < ,  a ==  b or  a >
 *  b respectively.
 *
 * \param a Matrix.
 * \param b Matrix.
 *
 * \note This comparison is not well defined mathematically and
 * relatively arbitrary since elements of GF(2) don't have an
 * ordering.
 */

int mzd_cmp(packedmatrix *a, packedmatrix *b);

/**
 * Copy matrix  p to  dst.
 *
 * \param dst May be NULL for automatic creation.
 * \param p Source matrix.
 */

packedmatrix *mzd_copy(packedmatrix *dst, packedmatrix *p);

/**
 * \brief Concatenate  b to  a.
 * 
 * That is,
 *
\verbatim
[ A ], [ B ] -> [ A  B ]
\endverbatim
 *
 *  The inputs are not modified but a new matrix is created.
 *
 * \param a Matrix
 * \param b Matrix
 *
 * \note This is sometimes called augment.
 */
packedmatrix *mzd_concat( packedmatrix *a, packedmatrix *b);

/**
 * \brief Stack  a on to of a.
 *
 * That is, 
 *
\verbatim
[ A ], [ B ] -> [ A ]
                [ B ]
\endverbatim
 *
 *  The inputs are not modified but a new matrix is created.
 *
 * \param a Matrix
 * \param b Matrix
 */

packedmatrix *mzd_stack( packedmatrix *a, packedmatrix *b);

/**
 * Copy a submatrix.
 * 
 * Note that the upper bounds are not included.
 *
 * \param a matrix
 * \param lowr start rows
 * \param lowc start column
 * \param highr stop row (this row is NOT included)
 * \param highc stop column (this column is NOT included)
 */
packedmatrix *mzd_submatrix(packedmatrix *a, int lowr, int lowc,
			    int highr, int highc);

/**
 * Invert the matrix target using Gaussian elimination. To avoid
 * recomputing the identity matrix over and over again, I may be
 * passed in as identity parameter.
 *
 * \param target Matrix to be reduced.
 * \param identity Identity matrix.
 *
 */

packedmatrix *mzd_invert_naiv(packedmatrix *target, packedmatrix *identity);

/**
 * Add left to right and write the result to ret. ret is also
 * returned.
 *
 * If ret is NULL then a new matrix is created which must be freed by
 * mzd_free.
 *
 * \param ret Matrix may be NULL for automatic creation.
 * \param left Matrix
 * \param right Matrix
 */

packedmatrix *mzd_add(packedmatrix *ret, packedmatrix *left, packedmatrix *right);

/**
 * Same as mzd_add but without any checks on the input.
 *
 * \param ret Matrix
 * \param left Matrix
 * \param right Matrix
 */

packedmatrix *_mzd_add_impl(packedmatrix *ret, packedmatrix *left, packedmatrix *right);

#endif //PACKEDMATRIX_H
