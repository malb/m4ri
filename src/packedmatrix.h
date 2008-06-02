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
   * Offsets to each row, so e.g. the first word of the i-th row
   * is m->values[m->rowswap[i]]
   */

  int *rowswap;

} packedmatrix;

/**
 * \brief Create a new matrix of dimension r x c.
 *
 * Use mzd_free to kill it.
 *
 * \param r Number of rows
 * \param c Number of columns
 *
 */

packedmatrix *mzd_init(const int r, const int c);

/**
 * \brief Free a matrix created with mzd_init.
 * 
 * \param A Matrix
 */

void mzd_free(packedmatrix *A);

/**
 * \brief Create a window/view into the matrix M.
 *
 * A matrix window for m is a meta structure on the matrix M. It is
 * setup to point into the matrix so M \em must \em not be freed while the
 * matrix window is used.
 *
 * This function puts restrictions on the provided parameters which
 * are not enforced currently.
 *
 *  - lowc must be divisible by RADIX
 *  - highc must be divisible by RADIX
 *  - all parameters must be within range for M
 *
 * Use mzd_free_free to free the window.
 *
 * \param M Matrix
 * \param lowr Starting row (inclusive)
 * \param lowc Starting column (inclusive)
 * \param highr End row (exclusive)
 * \param highc End column (exclusive)
 *
 */

packedmatrix *mzd_init_window(const packedmatrix *M, const int lowr, const int lowc, const int highr, const int highc);

/**
 * \brief Free a matrix window created with mzd_init_window.
 * 
 * \param A Matrix
 */

void mzd_free_window(packedmatrix *A);

/**
 * \brief Swap the two rows rowa and rowb.
 * 
 * \param M Matrix
 * \param rowa Row index
 * \param rowb Row index
 */
 
static inline void mzd_row_swap(packedmatrix *M, const int rowa, const int rowb) {
  /**
   * \todo it might be better to actually copy stuff around if the
   * rows are far apart to improve data locality.
   */
  int temp=M->rowswap[rowa];
  M->rowswap[rowa]=M->rowswap[rowb];
  M->rowswap[rowb]=temp;
}

/**
 * \brief Read the bit at position M[row,col].
 *
 * \param M Matrix
 * \param row Row index
 * \param col Column index
 *
 * \note No bounds checks whatsoever are performed.
 */

static inline BIT mzd_read_bit(const packedmatrix *M, const int row, const int col ) {
  return GET_BIT(M->values[ M->rowswap[row] + col/RADIX ], col%RADIX);
}

/**
 * \brief Write the bit value to position M[row,col]
 * 
 * \param M Matrix
 * \param row Row index
 * \param col Column index
 * \param value Either 0 or 1 
 *
 * \note No bounds checks whatsoever are performed.
 */

static inline void mzd_write_bit(packedmatrix *M, const int row, const int col, const BIT value) {
  if (value==1)
    SET_BIT(M->values[ M->rowswap[row] + col/RADIX ], col % RADIX);
  else
    CLR_BIT(M->values[ M->rowswap[row] + col/RADIX ], col % RADIX);
}

/**
 * \brief Add value to the word at position M[row,col].
 *
 * \param M Matrix
 * \param row Row index
 * \param col Column index
 * \param value Word of BITs.
 *
 * \note Keep in mind that the row, col refer to a row and column (of
 *  bits), and you can address the block by any of the RADIX (usually
 *  64) & A[i,j] there.
 *
 * \note No bounds checks whatsoever are performed.
 */

static inline void mzd_xor_block(packedmatrix *M, const int row, const int col, const word value) {
  int block=col/RADIX;
  int truerow=M->rowswap[row];

  word *entry=M->values + block + truerow;
  *entry ^= value;
}

/**
 * \brief Write value to the word at position M[row,col].
 *
 * \param M Matrix
 * \param row Row index
 * \param col Column index
 * \param value Word of BITs.
 *
 * \note Keep in mind that the row, col refer to a row and column (of
 * bits), and you can address the block by any of the RADIX (usually
 * 64) A[i,j] there.
 *
 * \note No bounds checks whatsoever are performed.
 */

static inline void mzd_write_block(packedmatrix *M, const int row, const int col, const word value) {
  M->values[ M->rowswap[row] + col/RADIX ] = value;
}

/**
 * \brief Read the word  at position M[row,col].
 *
 * \param M Matrix
 * \param row Row index
 * \param col Column index
 *
 * \note Keep in mind that the row, col refer to a row and column (of
 * bits), and you can address the block by any of the RADIX (usually
 * 64) A[i,j] there.
 *
 * \note No bounds checks whatsoever are performed.
 */

static inline word mzd_read_block(const packedmatrix *M, const int row, const int col ) {
  return M->values[ M->rowswap[row] + col/RADIX ];
}

/**
 * \brief Print a matrix to stdout. 
 *
 * The output will contain colons between every 4-th column.
 *
 * \param M Matrix
 *
 */

void mzd_print_matrix(const packedmatrix *M );

/**
 * \brief Print the matrix to stdout.
 *
 * \param M Matrix
 */

void mzd_print_matrix_tight(const packedmatrix *M );

/**
 * \brief Add the rows sourcerow and destrow and stores the total in the row
 * destrow, but only begins at the column coloffset.
 *
 * \param M Matrix
 * \param sourcerow Index of source row
 * \param destrow Index of target row
 * \param coloffset Column offset
 */

void mzd_row_add_offset(packedmatrix *M, const int sourcerow, const int destrow, const int coloffset );

/**
 * \brief Clear the given row, but only begins at the column coloffset.
 *
 * \param M Matrix
 * \param row Index of row
 * \param coloffset Column offset
 *
 */

void mzd_row_clear_offset(packedmatrix *M, const int row, const int coloffset);

/**
 * \brief Add the rows sourcerow and destrow and stores the total in
 * the row destrow.
 *
 * \param M Matrix
 * \param sourcerow Index of source row
 * \param destrow Index of target row
 *
 * \note this can be done much faster with mzd_combine.
 */

void mzd_row_add(packedmatrix *M, const int sourcerow, const int destrow);

/**
 * \brief Transpose a matrix.
 *
 * This function uses the fact that:
\verbatim
   [ A B ]T    [AT CT]
   [ C D ]  =  [BT DT] 
 \endverbatim 
 * and thus rearranges the blocks recursively. 
 *
 * \param DST Preallocated return matrix, may be NULL for automatic creation.
 * \param A Matrix
 */

packedmatrix *mzd_transpose(packedmatrix *DST, const packedmatrix *A );

/**
 * \brief Naive cubic matrix multiplication.
 *
 * That is, compute C such that C == AB.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A
 * \param B Input matrix B
 *
 * \note Normally, if you will multiply several times by b, it is
 * smarter to calculate bT yourself, and keep it, and then use the
 * function called matrixTimesMatrixTranspose
 */
packedmatrix *mzd_mul_naiv(packedmatrix *C, const packedmatrix *A, const packedmatrix *B);

/**
 * \brief Fill matrix M with uniformly distributed bits.
 *
 * \param M Matrix
 *
 * \todo Allow the user to provide a RNG callback.
 */

void mzd_randomize(packedmatrix *M );

/**
 * \brief Set the matrix M to the value equivalent to the integer
 * value provided.
 *
 * Specifically, this function does nothing if value%2 == 0 and
 * returns the identity matrix if value%2 == 1.
 *
 * If the matrix is not square then the largest possible square
 * submatrix is set to the identity matrix.
 *
 * \param M Matrix
 * \param value Either 0 or 1
 */

void mzd_set_ui(packedmatrix *M, const unsigned value);

/**
 * \brief Gaussian elimination.
 * 
 * This will do Gaussian elimination on the matrix m but will start
 * not at column 0 necc but at column startcol. If full=FALSE, then it
 * will do triangular style elimination, and if full=TRUE, it will do
 * Gauss-Jordan style, or full elimination.
 * 
 * \param M Matrix
 * \param startcol First column to consider for reduction.
 * \param full Gauss-Jordan style or upper triangular form only.
 */

int mzd_gauss_delayed(packedmatrix *M, const int startcol, const int full);

/**
 * \brief Gaussian elimination.
 * 
 * This will do Gaussian elimination on the matrix m.  If  full =
 *  FALSE, then it will do triangular style elimination, and if 
 *  full = TRUE, it will do Gauss-Jordan style, or full elimination.
 *
 * \param M Matrix
 * \param full Gauss-Jordan style or upper triangular form only.
 */
int mzd_reduce_naiv(packedmatrix *M, const int full);

/**
 * \brief Return TRUE if A == B.
 *
 * \param A Matrix
 * \param B Matrix
 */

BIT mzd_equal(const packedmatrix *A, const packedmatrix *B );

/**
 * \brief Return -1,0,1 if if A < B, A == B or A > B respectively.
 *
 * \param A Matrix.
 * \param B Matrix.
 *
 * \note This comparison is not well defined mathematically and
 * relatively arbitrary since elements of GF(2) don't have an
 * ordering.
 */

int mzd_cmp(const packedmatrix *A, const packedmatrix *B);

/**
 * \brief Copy matrix  A to DST.
 *
 * \param DST May be NULL for automatic creation.
 * \param A Source matrix.
 */

packedmatrix *mzd_copy(packedmatrix *DST, const packedmatrix *A);

/**
 * \brief Concatenate B to A and write the result to C.
 * 
 * That is,
 *
\verbatim
[ A ], [ B ] -> [ A  B ] = C
\endverbatim
 *
 * The inputs are not modified but a new matrix is created.
 *
 * \param C Matrix, may be NULL for automatic creation
 * \param A Matrix
 * \param B Matrix
 *
 * \note This is sometimes called augment.
 */
packedmatrix *mzd_concat(packedmatrix *C, const packedmatrix *A, const packedmatrix *B);

/**
 * \brief Stack A on top of B and write the result to C.
 *
 * That is, 
 *
\verbatim
[ A ], [ B ] -> [ A ] = C
                [ B ]
\endverbatim
 *
 * The inputs are not modified but a new matrix is created.
 *
 * \param C Matrix, may be NULL for automatic creation
 * \param A Matrix
 * \param B Matrix
 */

packedmatrix *mzd_stack(packedmatrix *C, const packedmatrix *A, const packedmatrix *B);

/**
 * \brief Copy a submatrix.
 * 
 * Note that the upper bounds are not included.
 *
 * \param S Preallocated space for submatrix, may be NULL for automatic creation.
 * \param M Matrix
 * \param lowr start rows
 * \param lowc start column
 * \param highr stop row (this row is \em not included)
 * \param highc stop column (this column is \em not included)
 */
packedmatrix *mzd_submatrix(packedmatrix *S, const packedmatrix *M, const int lowr, const int lowc, const int highr, const int highc);

/**
 * \brief Invert the matrix target using Gaussian elimination. 
 *
 * To avoid recomputing the identity matrix over and over again, I may
 * be passed in as identity parameter.
 *
 * \param INV Preallocated space for inversion matrix, may be NULL for automatic creation.
 * \param A Matrix to be reduced.
 * \param I Identity matrix.
 *
 */

packedmatrix *mzd_invert_naiv(packedmatrix *INV, packedmatrix *A, const packedmatrix *I);

/**
 * \brief Set C = A+B.
 *
 * C is also returned. If C is NULL then a new matrix is created which
 * must be freed by mzd_free.
 *
 * \param C Preallocated sum matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 */

packedmatrix *mzd_add(packedmatrix *C, const packedmatrix *A, const packedmatrix *B);

/**
 * \brief Same as mzd_add but without any checks on the input.
 *
 * \param C Preallocated sum matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 */

packedmatrix *_mzd_add_impl(packedmatrix *C, const packedmatrix *A, const packedmatrix *B);

/**
 * \brief Same as mzd_add.
 *
 * \param C Preallocated difference matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 */

#define mzd_sub mzd_add

/**
 * \brief Same as mzd_sub but without any checks on the input.
 *
 * \param C Preallocated difference matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 */

#define _mzd_sub_impl _mzd_add_impl

/**
 * \brief row3[col3:] = row1[col1:] + row2[col2:]
 * 
 * Adds row1 of SC1, starting with startblock1 to the end, to
 * row2 of SC2, starting with startblock2 to the end. This gets stored
 * in DST, in row3, starting with startblock3.
 *
 * \param DST destination matrix
 * \param row3 destination row for matrix dst
 * \param startblock3 starting block to work on in matrix dst
 * \param SC1 source matrix
 * \param row1 source row for matrix sc1
 * \param startblock1 starting block to work on in matrix sc1
 * \param SC2 source matrix
 * \param startblock2 starting block to work on in matrix sc2
 * \param row2 source row for matrix sc2
 */

void mzd_combine(packedmatrix * DST, const int row3, const int startblock3,
		 const packedmatrix * SC1, const int row1, const int startblock1, 
		 const packedmatrix * SC2, const int row2, const int startblock2);

#ifdef HAVE_SSE2
/**
 * Cutoff in words after which row length SSE2 instructions should be
 * used.
 */

#define SSE2_CUTOFF 20
#endif

/**
 * Defines the number of rows of the matrix A that are processed as
 * one block during the execution of a multiplication algorithm.
 */

#define MZD_MUL_BLOCKSIZE 768


#endif //PACKEDMATRIX_H
