/*
 * To enable dumping of output per function, configure the library with --enable-debug-dump.
 */

/******************************************************************************
*
*                 M4RI: Linear Algebra over GF(2)
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

#ifndef M4RI_DEBUG_DUMP
#define M4RI_DEBUG_DUMP

static inline word calculate_hash(word const* rowptr, wi_t wide) {
  word hash = 0;
  for (word const* ptr = rowptr; ptr < rowptr + wide; ++ptr) {
    hash ^= *ptr;
}
  return hash;
}

static inline word rotate_word(word w, int shift) {
  return (w << shift) | (w >> (m4ri_radix - w));
}

#if __M4RI_DEBUG_DUMP

struct mzd_t;
struct mzp_t;

extern void m4ri_dd_int(char const* function, char const* file, int line, int i);
extern void m4ri_dd_rci(char const* function, char const* file, int line, rci_t rci);
extern void m4ri_dd_rci_array(char const* function, char const* file, int line, rci_t *rciptr, int len);
extern void m4ri_dd_rawrow(char const* function, char const* file, int line, word const* rowptr, wi_t wide);
extern void m4ri_dd_row(char const* function, char const* file, int line, struct mzd_t const* M, rci_t row);
extern void m4ri_dd_mzd(char const* function, char const* file, int line, struct mzd_t const* M);
extern void m4ri_dd_mzp(char const* function, char const* file, int line, struct mzp_t const* P);

#define __M4RI_DD_INT(i) m4ri_dd_int(__FUNCTION__, __FILE__, __LINE__, i)
#define __M4RI_DD_RCI(rci) m4ri_dd_rci(__FUNCTION__, __FILE__, __LINE__, rci)
#define __M4RI_DD_RCI_ARRAY(rciptr, len) m4ri_dd_rci_array(__FUNCTION__, __FILE__, __LINE__, rciptr, len)
#define __M4RI_DD_RAWROW(rowptr, wide) m4ri_dd_rawrow(__FUNCTION__, __FILE__, __LINE__, rowptr, wide)
#define __M4RI_DD_ROW(M, row) m4ri_dd_row(__FUNCTION__, __FILE__, __LINE__, M, row)
#define __M4RI_DD_MZD(M) m4ri_dd_mzd(__FUNCTION__, __FILE__, __LINE__, M)
#define __M4RI_DD_MZP(P) m4ri_dd_mzp(__FUNCTION__, __FILE__, __LINE__, P)

#else // __M4RI_DEBUG_DUMP

#define __M4RI_DD_INT(i)
#define __M4RI_DD_RCI(rci)
#define __M4RI_DD_RCI_ARRAY(rciptr, len)
#define __M4RI_DD_RAWROW(rowptr, wide)
#define __M4RI_DD_ROW(M, row)
#define __M4RI_DD_MZD(M)
#define __M4RI_DD_MZP(P)

#endif // __M4RI_DEBUG_DUMP

#endif // M4RI_DEBUG_DUMP
