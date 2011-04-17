/******************************************************************************
*
*            M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2011 Carlo Wood <carlo@alinoe.com>
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
******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "packedmatrix.h"
#include "permutation.h"

#if __M4RI_DEBUG_DUMP

static unsigned long dd_sequence_number = 0;

static void entry(char const* function, char const* file, int line)
{
  printf("Sequence#: %ld; %s @ %s:%d; ", dd_sequence_number, function, file, line);
  ++dd_sequence_number;
}

static word calculate_hash(word const* rowptr, wi_t wide)
{
  unsigned long long hash = 0;
  for (word const* ptr = rowptr; ptr < rowptr + wide; ++ptr)
    hash ^= *ptr;
  return hash;
}

static inline word rotate_word(word w, int shift)
{
  return (w << shift) | (w >> (m4ri_radix - w));
}

void m4ri_dd_int(char const* function, char const* file, int line, int i)
{
  entry(function, file, line);
  printf("int: %d\n", i);
}

void m4ri_dd_rci(char const* function, char const* file, int line, rci_t rci)
{
  entry(function, file, line);
  printf("rci: %d\n", rci);
}

void m4ri_dd_rci_array(char const* function, char const* file, int line, rci_t *rciptr, int len)
{
  entry(function, file, line);
  unsigned long long hash = 0;
  for (int i = 0; i < len; ++i)
    hash ^= rotate_word(rciptr[i], i % m4ri_radix);
  printf("rci array (size %d) hash: %llx\n", len, hash);
}

void m4ri_dd_rawrow(char const* function, char const* file, int line, word const* rowptr, wi_t wide)
{
  entry(function, file, line);
  unsigned long long hash = calculate_hash(rowptr, wide);
  printf("raw row (%d words) hash: %llx\n", wide, hash);
}

void m4ri_dd_row(char const* function, char const* file, int line, mzd_t const* M, rci_t row)
{
  entry(function, file, line);
  unsigned long long hash = calculate_hash(M->rows[row], M->width);
  printf("row %d hash: %llx\n", row, hash);
}

void m4ri_dd_mzd(char const* function, char const* file, int line, mzd_t const* M)
{
  entry(function, file, line);
  unsigned long long hash = 0;
  for (rci_t r = 0; r < M->nrows; ++r)
    hash ^= rotate_word(calculate_hash(M->rows[r], M->width), r % m4ri_radix);
  printf("mzd hash: %llx\n", hash);
}

void m4ri_dd_mzp(char const* function, char const* file, int line, mzp_t const* P)
{
  entry(function, file, line);
  unsigned long long hash = 0;
  for (rci_t i = 0; i < P->length; ++i)
    hash ^= rotate_word(P->values[i], i % m4ri_radix);
  printf("mzp hash: %llx\n", hash);
}

#endif
