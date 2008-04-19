#ifndef BRILLIANTRUSSIAN_H
#define BRILLIANTRUSSIAN_H

/**
 *
 *            M4RI: Method of the Four Russians Inversion
 *
 *       Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org> 
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
 *  
 * See http://eprint.iacr.org/2006/251.pdf for details of the used
 * algorithm.
 */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "misc.h"
#include "packedmatrix.h"

/**
 * Finds a pivot row between xstart and xstop. The column where this
 * pivot is search is y. Returns YES if such a pivot row was
 * found. Also, the appropriate row is swapped to the top (== xstart).
 *
 * 
 * @param m matrix to operate on
 * @param xstart start row
 * @param xstop stop row (including)
 * @param y column to read
 *
 * @return True if a pivot row was found
 */

int forceNonZero(packedmatrix *m, int xstart, int xstop, int y);


/***************************************************/

/**
 * Performs Gaussian elimination on a submatrix of 3k x k starting at
 * point (homepoint, homepoint) of m.
 *
 * @param m matrix to operate on
 * @param homepoint row,col where to start
 * @param k
 *
 * @return rank of 3k x k submatrix.
 */

int prep(packedmatrix *m, int ai, int k);

/**
 * Adds row1 of s1, starting with startblock1 to the end, to row2 of
 * s2, starting with startblock2 to the end. This gets stored in dest,
 * in row3, starting with startblock3.
 *
 * Almost all computation time is spent in this function and thus
 * implementation improvements focus here.
 *
 * @param sc1 source matrix
 * @param row1 sourc row for matrix sc1
 * @param startblock1 starting block to work on in matrix sc1
 * @param sc2 source matrix
 * @param startblock2 starting block to work on in matrix sc2
 * @param row2 sourc row for matrix sc2
 * @param dst destination matrix
 * @param row3 sourc row for matrix dst
 * @param startblock3 starting block to work on in matrix dst
 *
 * 
@verbatim
  row3[col3:] = row1[col1:] + row2[col2:]
@endverbatim
 * 
 */

void m2t_combine( packedmatrix * sc1, int row1, int startblock1, 
		  packedmatrix * sc2, int row2, int startblock2,
		  packedmatrix * dst, int row3, int startblock3 );

/**
 * Constructs all possible $2^k$ row combinations using the gray code
 * table.
 * 
 * @param m matrix to operate on
 * @param ai the starting position
 * @param k
 * @param tablepacked prealloced matrix of dimension $2^k$ x m->ncols
 * @param lookuppacked prealloced table of length $2^k$
 * @param full touch columns before ai?
 *
 */

void makeTable( packedmatrix *m, int ai, int k, packedmatrix *tablepacked, int *lookuppacked, int full);


/**
 * Adds the correct row from tablepacked to the row 'row' in 'm'
 * starting at homecol.
 * 
 * @param m matrix to operate on
 * @param row the row which is operated on
 * @param homecol starting column for addition
 * @param k
 * @param tablepacked contains the correct row to be added
 * @param lookuptable contains row number to be addede
 */

void processRow(packedmatrix *m, int row, int homecol, int k, packedmatrix *tablepacked, int *lookuppacked);

/**
 *  Iterates proccessRow from startrow to stoprow.
 */

void process(packedmatrix *m, int startrow, int stoprow, int startcol, int k, packedmatrix *tablepacked, int *lookuppacked);

/**
 * This is the actual heart of the M4RI algorithm. 
 *
 * @param m the matrix
 * @param full perform full reduction
 * @param k 
 * @param ai start column
 * @param tablepacked buffer
 * @param lookuppacked buffer
 */

int stepM4RI(packedmatrix *m, int full, int k, int ai, packedmatrix *tablepacked, int *lookuppacked);

/**
 * Perform matrix reduction using the 'Method of the Four Russians' or
 * Kronrod-Method.
 * 
 * @param m the matrix to be reduced
 * @param full return the reduced row echelon form, not only the row echelon form
 * @param k
 * @param tablepacked preallocated table, may be NULL for automatic creation
 * @param lookuppacked preallocated lookup table, may be NULL for automatic creation
 */

int reduceM4RI(packedmatrix *m, int full, int k, packedmatrix *tablepacked, int *lookuppacked);

/**
 * Given a matrix in row echelon form compute the reduced row echelon
 * form of that matrix.
 * 
 * @param m the matrix to be reduced
 * @param full return the reduced row echelon form, not only the row echelon form
 * @param k
 * @param tablepacked preallocated table, may be NULL for automatic creation
 * @param lookuppacked preallocated lookup table, may be NULL for automatic creation
 */

void topReduceM4RI(packedmatrix *m, int k, packedmatrix *tablepacked, int *lookuppacked);

/**
 * Inverts the matrix m using Konrod's method. To avoid recomputing
 * the identity matrix over and over again, I may be passed in as
 * identity parameter.
 */

packedmatrix *invertM4RI(packedmatrix *m, packedmatrix *identity, int k);


/**
 * Matrix multiplication using Gray tables.
 *
 */

packedmatrix *m2t_mul_m4rm(packedmatrix *ret, packedmatrix *A, packedmatrix *B, int k, packedmatrix *tablepacked, int *lookuppacked);

packedmatrix *m2t_mul_m4rm_t(packedmatrix *ret, packedmatrix *A, packedmatrix *B, int k);

/**
 * Matrix multiplication via Strasen + Gray codes.
 *
 */

packedmatrix *m2t_mul_strassen(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff);

#endif //BRILLIANTRUSSIAN_H
