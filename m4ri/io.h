/**
 * \file io.h
 * \brief Input/output routines for matrices
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RI_IO_H
#define M4RI_IO_H

/*******************************************************************
 *
 *                M4RI: Linear Algebra over GF(2)
 *
 *    Copyright (C) 2011 Martin Albrecht <martinralbrecht@googlemail.com>
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

#include <m4ri/m4ri_config.h>
#include <m4ri/mzd.h>

/**
 * \brief Print row i of M to an output stream.
 *
 * The output will contain colons between every 4-th column.
 *
 * \param stream Output stream
 * \param M Matrix
 * \param i Row to print
 */

void mzd_fprint_row(FILE *stream, mzd_t const *M, const rci_t i);

/**
 * \brief Print row i of M to stdout.
 *
 * The output will contain colons between every 4-th column.
 *
 * \param M Matrix
 * \param i Row to print
 */

static inline void mzd_print_row(mzd_t const *M, const rci_t i) { mzd_fprint_row(stdout, M, i); }

/**
 * \brief Print a matrix to an output stream.
 *
 * The output will contain colons between every 4-th column.
 *
 * \param stream Output stream
 * \param M Matrix
 */

static inline void mzd_fprint(FILE *stream, mzd_t const *M) {
  for (rci_t i = 0; i < M->nrows; ++i) { mzd_fprint_row(stream, M, i); }
}

/**
 * \brief Print a matrix to stdout.
 *
 * The output will contain colons between every 4-th column.
 *
 * \param M Matrix
 */

static inline void mzd_print(mzd_t const *M) { mzd_fprint(stdout, M); }

/**
 * \brief Print compact information about the matrix to stdout.
 *
 * Prints number of rows, number of columns, density (and rank).
 *
 * \param A Matrix
 * \param do_rank Also display the rank (expensive)
 */

void mzd_info(const mzd_t *A, int do_rank);

#if __M4RI_HAVE_LIBPNG

/**
 * \brief Read matrix from 1-bit PNG image.
 *
 * This function returns a matrix on success and NULL otherwise. 1-bit
 * Grayscale and 1-bit Palette images are supported.
 *
 * \param fn Filename
 * \param verbose Print error message to stdout if != 0
 */

mzd_t *mzd_from_png(const char *fn, int verbose);

/**
 * \brief Write matrix to 1-bit PNG image.
 *
 * This function returns zero on success and some value != 0
 * otherwise. The parameter compression_level takes a zlib compression
 * level, i.e., an integer betweeen -1 and 9 (inclusive) such that
 *
\verbatim
#define Z_NO_COMPRESSION         0
#define Z_BEST_SPEED             1
#define Z_BEST_COMPRESSION       9
#define Z_DEFAULT_COMPRESSION  (-1)
\endverbatim
 *
 * The optional comment string is written as a PNG comment.
 *
 *
 * \param A Matrix
 * \param fn Filename (must have write permission)
 * \param compression_level Zlib compression level (see above)
 * \param comment Optional comment (may be NULL)
 * \param verbose Print error message to stdout if != 0
 */

int mzd_to_png(const mzd_t *A, const char *fn, int compression_level, const char *comment,
               int verbose);

#endif  //__M4RI_HAVE_LIBPNG

/**
 * \brief Read matrix from ASCII file in JCF format.
 *
 * The format is as follows:
\verbatim
nrows ncols modulus
nonzero_entries_upper_bound
column_index
\endverbatim
 *
 * where a negative column_index indicates a row_index increase by one and a non-zero entry at index
 * -column_index.
 *
 * \note the JCF format is one-based in contrast to everything else in this library which is
 * zero-based.
 *
 * For example, a valid input is:
\verbatim
3 2 2
3

-2
-1
-2
\endverbatim
 *
 * which produces the matrix
\verbatim
[0 1]
[1 0]
[0 1]
\endverbatim
 *
 *
 * \param fn Filename
 * \param verbose Print error message to stdout if != 0
 */

mzd_t *mzd_from_jcf(const char *fn, int verbose);

/**
 * \brief Create matrix from dense ASCII string
 *
 * The provided string is parsed in row major ordering, i.e. the first entry is
 * writen to A[0,0], the second entry to A[0,1] etc.
 *
 * For example, calling
\verbatim
mzd_t *A = mzd_from_str(4, 4, "1000010000100001");
\endverbatim
 *
 * would create a 4 x 4 identity matrix.
 *
 * \param m Number of rows
 * \param n Nimber of columns
 * \param str String containing ASCII zeros and ones of length m*n
 */

mzd_t *mzd_from_str(rci_t m, rci_t n, const char *str);

#endif  // M4RI_IO_H
