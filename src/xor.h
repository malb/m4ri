/**
 * \file xor.h
 * \brief Functions for adding vectors.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 *
 * \todo start counting at 0!
 */

#ifndef M4RI_XOR_H
#define M4RI_XOR_H

 /*******************************************************************
 *
 *                 M4RI:  Linear Algebra over GF(2)
 *
 *    Copyright (C) 2008-2010  Martin Albrecht <martinralbrecht@googlemail.com>
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
 * Compute c[i] += t1[i] + t2[i] + t3[i] + t4[i] + t5[i] + t6[i] + t7[i] + t8[i] for 0 <= i < wide
 *
 * \todo the non SSE2 version of this code is slow, replace by code
 * from mzd_process_rows8
 */

static inline void _mzd_combine8(word *c, word const *t1, word const *t2, word const *t3, word const *t4,
                                 word const *t5, word const *t6, word const *t7, word const *t8, wi_t wide_in) {
  wi_t wide = wide_in;
#if __M4RI_HAVE_SSE2
  /* assuming t1 ... t8 are aligned, but c might not be */
  if (__M4RI_ALIGNMENT(c,16)==0) {
    __m128i *__c = (__m128i*)c;
    __m128i *__t1 = (__m128i*)t1;
    __m128i *__t2 = (__m128i*)t2;
    __m128i *__t3 = (__m128i*)t3;
    __m128i *__t4 = (__m128i*)t4;
    __m128i *__t5 = (__m128i*)t5;
    __m128i *__t6 = (__m128i*)t6;
    __m128i *__t7 = (__m128i*)t7;
    __m128i *__t8 = (__m128i*)t8;
    const __m128i *eof = (__m128i*)((unsigned long)(c + wide) & ~0xFUL);
    __m128i xmm1;
    
    while(__c < eof) {
      xmm1 = _mm_xor_si128(*__c, *__t1++);
      xmm1 = _mm_xor_si128(xmm1, *__t2++);
      xmm1 = _mm_xor_si128(xmm1, *__t3++);
      xmm1 = _mm_xor_si128(xmm1, *__t4++);
      xmm1 = _mm_xor_si128(xmm1, *__t5++);
      xmm1 = _mm_xor_si128(xmm1, *__t6++);
      xmm1 = _mm_xor_si128(xmm1, *__t7++);
      xmm1 = _mm_xor_si128(xmm1, *__t8++);
      *__c++ = xmm1;
    }
    c  = (word*)__c;
    t1 = (word*)__t1;
    t2 = (word*)__t2;
    t3 = (word*)__t3;
    t4 = (word*)__t4;
    t5 = (word*)__t5;
    t6 = (word*)__t6;
    t7 = (word*)__t7;
    t8 = (word*)__t8;
    wide = ((sizeof(word) * wide) % 16) / sizeof(word);
  }
#endif
  for(wi_t i = 0; i < wide; ++i) {
    c[i] ^= t1[i] ^ t2[i] ^ t3[i] ^ t4[i] ^ t5[i] ^ t6[i] ^ t7[i] ^ t8[i];
  }

  __M4RI_DD_RAWROW(c, wide_in);
}
/**
 * Compute c[i] += t1[i] + t2[i] + t3[i] + t4[i] + t5[i] + t6[i] for 0 <= i < wide
 *
 */

static inline void _mzd_combine6(word *c, word const *t1, word const *t2, word const *t3, word const *t4, word const *t5, word const *t6, wi_t wide_in) {
  wi_t wide = wide_in;
#if __M4RI_HAVE_SSE2
  /* assuming t1 ... t4 are aligned, but c might not be */
  assert(__M4RI_ALIGNMENT(c,16) == 8 || __M4RI_ALIGNMENT(c,16) == 0);
  assert(__M4RI_ALIGNMENT(c,16) == __M4RI_ALIGNMENT(t1,16));

  if (__M4RI_ALIGNMENT(c,16) == 8) {
    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++ ^ *t6++;
    wide--;
  }

  __m128i *__c = (__m128i*)c;
  __m128i *__t1 = (__m128i*)t1;
  __m128i *__t2 = (__m128i*)t2;
  __m128i *__t3 = (__m128i*)t3;
  __m128i *__t4 = (__m128i*)t4;
  __m128i *__t5 = (__m128i*)t5;
  __m128i *__t6 = (__m128i*)t6;
  const __m128i *eof = (__m128i*)((unsigned long)(c + wide) & ~0xFUL);
  __m128i xmm1;

  while(__c < eof) {
    xmm1 = _mm_xor_si128(*__c, *__t1++);
    xmm1 = _mm_xor_si128(xmm1, *__t2++);
    xmm1 = _mm_xor_si128(xmm1, *__t3++);
    xmm1 = _mm_xor_si128(xmm1, *__t4++);
    xmm1 = _mm_xor_si128(xmm1, *__t5++);
    xmm1 = _mm_xor_si128(xmm1, *__t6++);
    *__c++ = xmm1;
  }
  c  = (word*)__c;
  t1 = (word*)__t1;
  t2 = (word*)__t2;
  t3 = (word*)__t3;
  t4 = (word*)__t4;
  t5 = (word*)__t5;
  t6 = (word*)__t6;
  wide = ((sizeof(word) * wide) % 16) / sizeof(word);

  if(wide) {
    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++ ^ *t6++;
  }
  __M4RI_DD_RAWROW(c, wide_in);
  return;
#else
  wi_t n = (wide + 7) / 8;
  switch (wide % 8) {
  case 0: do { *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++ ^ *t6++;
    case 7:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++ ^ *t6++;
    case 6:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++ ^ *t6++;
    case 5:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++ ^ *t6++;
    case 4:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++ ^ *t6++;
    case 3:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++ ^ *t6++;
    case 2:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++ ^ *t6++;
    case 1:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++ ^ *t6++;
    } while (--n > 0);
  }
  __M4RI_DD_RAWROW(c, wide_in);
  return;
#endif // __M4RI_HAVE_SSE2
}


/**
 * Compute c[i] += t1[i] + t2[i] + t3[i] + t4[i] + t5[i] for 0 <= i < wide
 *
 */

static inline void _mzd_combine5(word *c, word const *t1, word const *t2, word const *t3, word const *t4, word const *t5, wi_t wide_in) {
  wi_t wide = wide_in;
#if __M4RI_HAVE_SSE2
  /* assuming t1 ... t4 are aligned, but c might not be */
  assert(__M4RI_ALIGNMENT(c,16) == 8 || __M4RI_ALIGNMENT(c,16) == 0);
  assert(__M4RI_ALIGNMENT(c,16) == __M4RI_ALIGNMENT(t1,16));

  if (__M4RI_ALIGNMENT(c,16) == 8) {
    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
    wide--;
  }

  __m128i *__c = (__m128i*)c;
  __m128i *__t1 = (__m128i*)t1;
  __m128i *__t2 = (__m128i*)t2;
  __m128i *__t3 = (__m128i*)t3;
  __m128i *__t4 = (__m128i*)t4;
  __m128i *__t5 = (__m128i*)t5;
  const __m128i *eof = (__m128i*)((unsigned long)(c + wide) & ~0xFUL);
  __m128i xmm1;

  while(__c < eof) {
    xmm1 = _mm_xor_si128(*__c, *__t1++);
    xmm1 = _mm_xor_si128(xmm1, *__t2++);
    xmm1 = _mm_xor_si128(xmm1, *__t3++);
    xmm1 = _mm_xor_si128(xmm1, *__t4++);
    xmm1 = _mm_xor_si128(xmm1, *__t5++);
    *__c++ = xmm1;
  }
  c  = (word*)__c;
  t1 = (word*)__t1;
  t2 = (word*)__t2;
  t3 = (word*)__t3;
  t4 = (word*)__t4;
  t5 = (word*)__t5;
  wide = ((sizeof(word) * wide) % 16) / sizeof(word);

  if(wide) {
    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
  }
  __M4RI_DD_RAWROW(c, wide_in);
  return;
#else
  wi_t n = (wide + 7) / 8;
  switch (wide % 8) {
  case 0: do { *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
    case 7:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
    case 6:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
    case 5:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
    case 4:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
    case 3:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
    case 2:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
    case 1:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
    } while (--n > 0);
  }
  __M4RI_DD_RAWROW(c, wide_in);
  return;
#endif // __M4RI_HAVE_SSE2
}


/**
 * Compute c[i] += t1[i] + t2[i] + t3[i] + t4[i] for 0 <= i < wide
 *
 */

static inline void _mzd_combine4(word *c, word const *t1, word const *t2, word const *t3, word const *t4, wi_t wide_in) {
  wi_t wide = wide_in;
#if __M4RI_HAVE_SSE2
  /* assuming t1 ... t4 are aligned, but c might not be */
  assert(__M4RI_ALIGNMENT(c,16) == 8 || __M4RI_ALIGNMENT(c,16) == 0);
  assert(__M4RI_ALIGNMENT(c,16) == __M4RI_ALIGNMENT(t1,16));

  if (__M4RI_ALIGNMENT(c,16) == 8) {
    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
    wide--;
  }

  __m128i *__c = (__m128i*)c;
  __m128i *__t1 = (__m128i*)t1;
  __m128i *__t2 = (__m128i*)t2;
  __m128i *__t3 = (__m128i*)t3;
  __m128i *__t4 = (__m128i*)t4;
  const __m128i *eof = (__m128i*)((unsigned long)(c + wide) & ~0xFUL);
  __m128i xmm1;

  while(__c < eof) {
    xmm1 = _mm_xor_si128(*__c, *__t1++);
    xmm1 = _mm_xor_si128(xmm1, *__t2++);
    xmm1 = _mm_xor_si128(xmm1, *__t3++);
    xmm1 = _mm_xor_si128(xmm1, *__t4++);
    *__c++ = xmm1;
  }
  c  = (word*)__c;
  t1 = (word*)__t1;
  t2 = (word*)__t2;
  t3 = (word*)__t3;
  t4 = (word*)__t4;
  wide = ((sizeof(word) * wide) % 16) / sizeof(word);

  if(wide) {
    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
  }
  __M4RI_DD_RAWROW(c, wide_in);
  return;
#else
  wi_t n = (wide + 7) / 8;
  switch (wide % 8) {
  case 0: do { *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
    case 7:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
    case 6:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
    case 5:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
    case 4:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
    case 3:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
    case 2:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
    case 1:    *c++ ^= *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
    } while (--n > 0);
  }
  __M4RI_DD_RAWROW(c, wide_in);
  return;
#endif // __M4RI_HAVE_SSE2
}

/**
 * Compute c[i] += t1[i] + t2[i] + t3[i] for 0 <= i < wide
 *
 */

static inline void _mzd_combine3(word *c, word const *t1, word const *t2, word const *t3, wi_t wide_in) {
  wi_t wide = wide_in;
#if __M4RI_HAVE_SSE2
  /* assuming t1 ... t3 are aligned, but c might not be */
  if (__M4RI_ALIGNMENT(c,16)==0) {
    __m128i *__c = (__m128i*)c;
    __m128i *__t1 = (__m128i*)t1;
    __m128i *__t2 = (__m128i*)t2;
    __m128i *__t3 = (__m128i*)t3;
    const __m128i *eof = (__m128i*)((unsigned long)(c + wide) & ~0xFUL);
    __m128i xmm1;
    
    while(__c < eof) {
      xmm1 = _mm_xor_si128(*__c, *__t1++);
      xmm1 = _mm_xor_si128(xmm1, *__t2++);
      xmm1 = _mm_xor_si128(xmm1, *__t3++);
      *__c++ = xmm1;
    }
    c  = (word*)__c;
    t1 = (word*)__t1;
    t2 = (word*)__t2;
    t3 = (word*)__t3;
    wide = ((sizeof(word) * wide) % 16) / sizeof(word);
  }
  if(!wide) {
    __M4RI_DD_RAWROW(c, wide_in);
    return;
  }
#endif // __M4RI_HAVE_SSE2
  wi_t n = (wide + 7) / 8;
  switch (wide % 8) {
  case 0: do { *c++ ^= *t1++ ^ *t2++ ^ *t3++;
    case 7:    *c++ ^= *t1++ ^ *t2++ ^ *t3++;
    case 6:    *c++ ^= *t1++ ^ *t2++ ^ *t3++;
    case 5:    *c++ ^= *t1++ ^ *t2++ ^ *t3++;
    case 4:    *c++ ^= *t1++ ^ *t2++ ^ *t3++;
    case 3:    *c++ ^= *t1++ ^ *t2++ ^ *t3++;
    case 2:    *c++ ^= *t1++ ^ *t2++ ^ *t3++;
    case 1:    *c++ ^= *t1++ ^ *t2++ ^ *t3++;
    } while (--n > 0);
  }
  __M4RI_DD_RAWROW(c, wide_in);
}


/**
 * Compute c[i] += t1[i] + t2[i] for 0 <= i < wide
 *
 */

static inline void _mzd_combine2(word *c, word const *t1, word const *t2, wi_t wide_in) {
  wi_t wide = wide_in;
#if __M4RI_HAVE_SSE2
  /* assuming t1 ... t2 are aligned, but c might not be */
  if (__M4RI_ALIGNMENT(c,16)==0) {
    __m128i *__c = (__m128i*)c;
    __m128i *__t1 = (__m128i*)t1;
    __m128i *__t2 = (__m128i*)t2;
    const __m128i *eof = (__m128i*)((unsigned long)(c + wide) & ~0xFUL);
    __m128i xmm1;
    
    while(__c < eof) {
      xmm1 = _mm_xor_si128(*__c, *__t1++);
      xmm1 = _mm_xor_si128(xmm1, *__t2++);
      *__c++ = xmm1;
    }
    c  = (word*)__c;
    t1 = (word*)__t1;
    t2 = (word*)__t2;
    wide = ((sizeof(word) * wide) % 16) / sizeof(word);
  }
  if(!wide) {
    __M4RI_DD_RAWROW(c, wide_in);
    return;
  }
#endif // __M4RI_HAVE_SSE2
  wi_t n = (wide + 7) / 8;
  switch (wide % 8) {
  case 0: do { *c++ ^= *t1++ ^ *t2++;
    case 7:    *c++ ^= *t1++ ^ *t2++;
    case 6:    *c++ ^= *t1++ ^ *t2++;
    case 5:    *c++ ^= *t1++ ^ *t2++;
    case 4:    *c++ ^= *t1++ ^ *t2++;
    case 3:    *c++ ^= *t1++ ^ *t2++;
    case 2:    *c++ ^= *t1++ ^ *t2++;
    case 1:    *c++ ^= *t1++ ^ *t2++;
    } while (--n > 0);
  }
  __M4RI_DD_RAWROW(c, wide_in);
}

/**
 * Compute c[i] += t1[i] for 0 <= i < wide
 *
 */

static inline void _mzd_combine(word *c, word const *t1, wi_t wide_in) {
  wi_t wide = wide_in;
#if __M4RI_HAVE_SSE2
  /* assuming c, t1 are alligned the same way */

  if (__M4RI_ALIGNMENT(c,16)==8 && wide) {
    *c++ ^= *t1++;
    wide--;
  }

  __m128i *__c = (__m128i*)c;
  __m128i *__t1 = (__m128i*)t1;
  const __m128i *eof = (__m128i*)((unsigned long)(c + wide) & ~0xFUL);
  __m128i xmm1;
  
  
  while(__c < eof-1) {
    xmm1 = _mm_xor_si128(*__c, *__t1++);
    *__c++ = xmm1;
    xmm1 = _mm_xor_si128(*__c, *__t1++);
    *__c++ = xmm1;
  }

  if(__c < eof) {
    xmm1 = _mm_xor_si128(*__c, *__t1++); 
    *__c++ = xmm1;      
  }
  
  c  = (word*)__c;
  t1 = (word*)__t1;
  wide = ((sizeof(word) * wide) % 16) / sizeof(word);

  if(!wide) {
    __M4RI_DD_RAWROW(c, wide_in);
    return;
  }
#endif // __M4RI_HAVE_SSE2

  wi_t n = (wide + 7) / 8;
  switch (wide % 8) {
  case 0: do { *c++ ^= *t1++;
    case 7:    *c++ ^= *t1++;
    case 6:    *c++ ^= *t1++;
    case 5:    *c++ ^= *t1++;
    case 4:    *c++ ^= *t1++;
    case 3:    *c++ ^= *t1++;
    case 2:    *c++ ^= *t1++;
    case 1:    *c++ ^= *t1++;
    } while (--n > 0);
  }
  __M4RI_DD_RAWROW(c, wide_in);
}

#endif // M4RI_XOR_H
