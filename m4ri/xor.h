/*
 * Functions for adding vectors.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 *
 */

#ifndef M4RI_XOR_H
#define M4RI_XOR_H

/*******************************************************************
 *
 *                 M4RI:  Linear Algebra over GF(2)
 *
 *    Copyright (C) 2008-2013  Martin Albrecht <martinralbrecht@googlemail.com>
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

#if __M4RI_HAVE_SSE2
#include <emmintrin.h>
#endif

#include <m4ri/misc.h>

/**
 * Compute c[i] += t1[i] for 0 <= i < wide
 *
 */

static inline void _mzd_combine(word *c, word const *t1, wi_t wide_in) {
  wi_t wide = wide_in;
#if __M4RI_HAVE_SSE2
  /* assuming c, t1 are alligned the same way */

  if (__M4RI_ALIGNMENT(c, 16) == 8 && wide) {
    *c++ ^= *t1++;
    wide--;
  }

  __m128i *__c       = (__m128i *)c;
  __m128i *__t1      = (__m128i *)t1;
  const __m128i *eof = (__m128i *)((unsigned long)(c + wide) & ~0xFUL);
  __m128i xmm1;

  while (__c < eof - 1) {
    xmm1   = _mm_xor_si128(*__c, *__t1++);
    *__c++ = xmm1;
    xmm1   = _mm_xor_si128(*__c, *__t1++);
    *__c++ = xmm1;
  }

  if (__c < eof) {
    xmm1   = _mm_xor_si128(*__c, *__t1++);
    *__c++ = xmm1;
  }

  c    = (word *)__c;
  t1   = (word *)__t1;
  wide = ((sizeof(word) * wide) % 16) / sizeof(word);

  if (!wide) {
    __M4RI_DD_RAWROW(c, wide_in);
    return;
  }
#endif  // __M4RI_HAVE_SSE2

  wi_t n = (wide + 7) / 8;
  switch (wide % 8) {
  case 0: do { *c++ ^= *t1++;
    case 7: *c++ ^= *t1++;
    case 6: *c++ ^= *t1++;
    case 5: *c++ ^= *t1++;
    case 4: *c++ ^= *t1++;
    case 3: *c++ ^= *t1++;
    case 2: *c++ ^= *t1++;
    case 1: *c++ ^= *t1++;
    } while (--n > 0);
  }
  __M4RI_DD_RAWROW(c, wide_in);
}

#define N 2
#include "xor_template.h"
#undef N

#define N 3
#include "xor_template.h"
#undef N

#define N 4
#include "xor_template.h"
#undef N

#define N 5
#include "xor_template.h"
#undef N

#define N 6
#include "xor_template.h"
#undef N

#define N 7
#include "xor_template.h"
#undef N

#define N 8
#include "xor_template.h"
#undef N

#endif  // M4RI_XOR_H
