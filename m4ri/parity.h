#ifndef M4RI_PARITY_H
#define M4RI_PARITY_H

/*******************************************************************
 *
 *                 M4RI: Linear Algebra over GF(2)
 *
 *    Copyright (C) 2008 David Harvey <dmharvey@cims.nyu.edu>
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

/**
 * \file parity.h
 *
 * \brief Compute the parity of 64 words in parallel.
 *
 * \author David Harvey
 */

#include <m4ri/misc.h>

/**
 * \brief Step for mixing two 64-bit words to compute their parity.
 */

#define __M4RI_MIX32(a, b) (((((a) >> 32) ^ (a)) << 32) | ((((b) << 32) ^ (b)) >> 32))

/**
 * \brief Step for mixing two 64-bit words to compute their parity.
 */

#define __M4RI_MIX16(a, b)                                                                         \
  (((((a) << 16) ^ (a)) & __M4RI_CONVERT_TO_WORD(0xFFFF0000FFFF0000ull)) |                         \
   ((((b) >> 16) ^ (b)) & __M4RI_CONVERT_TO_WORD(0x0000FFFF0000FFFFull)));
/**
 * \brief Step for mixing two 64-bit words to compute their parity.
 */

#define __M4RI_MIX8(a, b)                                                                          \
  (((((a) << 8) ^ (a)) & __M4RI_CONVERT_TO_WORD(0xFF00FF00FF00FF00ull)) |                          \
   ((((b) >> 8) ^ (b)) & __M4RI_CONVERT_TO_WORD(0x00FF00FF00FF00FFull)));
/**
 * \brief Step for mixing two 64-bit words to compute their parity.
 */

#define __M4RI_MIX4(a, b)                                                                          \
  (((((a) << 4) ^ (a)) & __M4RI_CONVERT_TO_WORD(0xF0F0F0F0F0F0F0F0ull)) |                          \
   ((((b) >> 4) ^ (b)) & __M4RI_CONVERT_TO_WORD(0x0F0F0F0F0F0F0F0Full)));
/**
 * \brief Step for mixing two 64-bit words to compute their parity.
 */

#define __M4RI_MIX2(a, b)                                                                          \
  (((((a) << 2) ^ (a)) & __M4RI_CONVERT_TO_WORD(0xCCCCCCCCCCCCCCCCull)) |                          \
   ((((b) >> 2) ^ (b)) & __M4RI_CONVERT_TO_WORD(0x3333333333333333ull)));
/**
 * \brief Step for mixing two 64-bit words to compute their parity.
 */

#define __M4RI_MIX1(a, b)                                                                          \
  (((((a) << 1) ^ (a)) & __M4RI_CONVERT_TO_WORD(0xAAAAAAAAAAAAAAAAull)) |                          \
   ((((b) >> 1) ^ (b)) & __M4RI_CONVERT_TO_WORD(0x5555555555555555ull)));

/**
 * \brief See parity64.
 */

static inline word m4ri_parity64_helper(word *buf) {
  word a0, a1, b0, b1, c0, c1;

  a0 = __M4RI_MIX32(buf[0x20], buf[0x00]);
  a1 = __M4RI_MIX32(buf[0x30], buf[0x10]);
  b0 = __M4RI_MIX16(a1, a0);

  a0 = __M4RI_MIX32(buf[0x28], buf[0x08]);
  a1 = __M4RI_MIX32(buf[0x38], buf[0x18]);
  b1 = __M4RI_MIX16(a1, a0);

  c0 = __M4RI_MIX8(b1, b0);

  a0 = __M4RI_MIX32(buf[0x24], buf[0x04]);
  a1 = __M4RI_MIX32(buf[0x34], buf[0x14]);
  b0 = __M4RI_MIX16(a1, a0);

  a0 = __M4RI_MIX32(buf[0x2C], buf[0x0C]);
  a1 = __M4RI_MIX32(buf[0x3C], buf[0x1C]);
  b1 = __M4RI_MIX16(a1, a0);

  c1 = __M4RI_MIX8(b1, b0);

  return __M4RI_MIX4(c1, c0);
}

/**
 * \brief Computes parity of each of buf[0], buf[1], ..., buf[63].
 * Returns single word whose bits are the parities of buf[0], ...,
 * buf[63].
 *
 * \param buf buffer of words of length 64
 */
static inline word m4ri_parity64(word *buf) {
  word d0, d1, e0, e1;

  d0 = m4ri_parity64_helper(buf);
  d1 = m4ri_parity64_helper(buf + 2);
  e0 = __M4RI_MIX2(d1, d0);

  d0 = m4ri_parity64_helper(buf + 1);
  d1 = m4ri_parity64_helper(buf + 3);
  e1 = __M4RI_MIX2(d1, d0);

  return __M4RI_MIX1(e1, e0);
}

#endif  // M4RI_PARITY_H
