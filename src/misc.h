
/**
 * \file misc.h
 * \brief Helper functions.
 *
 * \author Gregory Bard <bard@fordham.edu>
 * \author Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
 * \author Carlo Wood <carlo@alinoe.com>
 */

#ifndef M4RI_MISC_H
#define M4RI_MISC_H

/*******************************************************************
*
*                 M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2007, 2008 Gregory Bard <bard@fordham.edu>
*    Copyright (C) 2008 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
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
*
********************************************************************/

#include <m4ri/m4ri_config.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if __M4RI_USE_MM_MALLOC
#include <mm_malloc.h>
#endif

#include <stdlib.h>
#include <assert.h>
#include <string.h>
/// @cond INTERNAL
#define __STDC_LIMIT_MACROS
/// @endcond
#include <stdint.h>

/*
 * These define entirely the word width used in the library.
 */

/**
 * \brief Pretty for a boolean int.
 *
 * The value of a BIT is either 0 or 1.
 */

typedef int BIT;

/**
 * \brief Type of row and column indexes.
 *
 * This type is used for integer values that hold row/colum sized values.
 */

typedef int rci_t;

/**
 * \brief Type of word indexes.
 *
 * This type is used for the array of words that make up a row.
 */

typedef int wi_t;


/**
 * \brief A word is the typical packed data structure to represent packed bits.
 */

typedef uint64_t word;

/**
 * \brief Explicit conversion macro.
 *
 * Explicit conversion of a word, representing 64 columns, to an integer
 * to be used as index into an array. This is used for Gray codes.
 * No error checking is done that the most significant bits in w are zero.
 *
 * \note This is a no-op. It's purpose it to track intention.
 */

#define __M4RI_CONVERT_TO_INT(w) ((int)(w))

/**
 * \brief Explicit conversion macro.
 *
 * Explicit conversion of a word, representing 64 columns, to a BIT
 * to be used as boolean: this is an int with value 0 (false) or 1 (true).
 * No error checking is done that only the least significant bit is set (if any).
 *
 * \note This is a no-op. It's purpose it to track intention.
 */

#define __M4RI_CONVERT_TO_BIT(w) ((BIT)(w))

/**
 * \brief Explicit conversion macro.
 *
 * Explicit conversion of a word, representing 64 columns, to an uint64_t.
 *
 * The returned value is the underlaying integer representation of these 64 columns,
 * meaning in particular that if val is an uint64_t then
 * __M4RI_CONVERT_TO_UINT64_T(__M4RI_CONVERT_TO_WORD(val)) == val.
 *
 * \note This is a no-op. It's purpose it to track intention.
 */

#define __M4RI_CONVERT_TO_UINT64_T(w) (w)

/**
 * \brief Explicit conversion macro.
 *
 * Explicit conversion of an integer to a word.
 *
 * \note This is a no-op. It's purpose it to track intention.
 */

#define __M4RI_CONVERT_TO_WORD(i) ((word)(i))

/**
 * \brief The number of bits in a word.
 */

static int const m4ri_radix = 64;

/**
 * \brief The number one as a word.
 */

static word const m4ri_one = __M4RI_CONVERT_TO_WORD(1);

/**
 * \brief A word with all bits set.
 */

static word const m4ri_ffff = __M4RI_CONVERT_TO_WORD(-1);

/**
 * \brief Return the maximal element of x and y
 *
 * \param x Word
 * \param y Word
 */

#ifndef MAX
#define MAX(x,y) (((x) > (y))?(x):(y))
#endif

/**
 * \brief Return the minimal element of x and y
 *
 * \param x Word
 * \param y Word
 */

#ifndef MIN
#define MIN(x,y) (((x) < (y))?(x):(y))
#endif

/**
 *\brief Pretty for 1.
 */ 

#ifndef TRUE
#define TRUE 1
#endif

/**
 *\brief Pretty for 0.
 */ 

#ifndef FALSE
#define FALSE 0
#endif

/**
 * \brief $2^i$
 *
 * \param i Integer.
 */ 

#define __M4RI_TWOPOW(i) ((uint64_t)1 << (i))

/**
 * \brief Clear the bit spot (counting from the left) in the word w
 * 
 * \param w Word
 * \param spot Integer with 0 <= spot < m4ri_radix
 */

#define __M4RI_CLR_BIT(w, spot) ((w) &= ~(m4ri_one << (spot))

/**
 * \brief Set the bit spot (counting from the left) in the word w
 * 
 * \param w Word
 * \param spot Integer with 0 <= spot < m4ri_radix
 */

#define __M4RI_SET_BIT(w, spot) ((w) |= (m4ri_one << (spot)))

/**
 * \brief Get the bit spot (counting from the left) in the word w
 * 
 * \param w Word
 * \param spot Integer with 0 <= spot < m4ri_radix
 */

#define __M4RI_GET_BIT(w, spot) __M4RI_CONVERT_TO_BIT(((w) >> (spot)) & m4ri_one)

/**
 * \brief Write the value to the bit spot in the word w
 * 
 * \param w Word.
 * \param spot Integer with 0 <= spot < m4ri_radix.
 * \param value Either 0 or 1.
 */

#define __M4RI_WRITE_BIT(w, spot, value) ((w) = (((w) & ~(m4ri_one << (spot))) | (-__M4RI_CONVERT_TO_WORD(value) & (m4ri_one << (spot)))))

/**
 * \brief Flip the spot in the word w
 *
 * \param w Word.
 * \param spot Integer with 0 <= spot < m4ri_radix.
 */

#define __M4RI_FLIP_BIT(w, spot) ((w) ^= (m4ri_one << (spot)))

/**
* \brief create a bit mask to zero out all but the (n - 1) % m4ri_radix + 1 leftmost bits.
*
* This function returns 1..64 bits, never zero bits.
* This mask is mainly used to mask the valid bits in the most significant word,
* by using __M4RI_LEFT_BITMASK((M->ncols + M->offset) % m4ri_radix).
* In other words, the set bits represent the columns with the lowest index in the word.
*
*  Thus,
*
*  n	Output
*  0=64 1111111111111111111111111111111111111111111111111111111111111111
*  1	0000000000000000000000000000000000000000000000000000000000000001
*  2    0000000000000000000000000000000000000000000000000000000000000011
*  .                                   ...
* 62    0011111111111111111111111111111111111111111111111111111111111111
* 63	0111111111111111111111111111111111111111111111111111111111111111
*
* Note that n == 64 is only passed from __M4RI_MIDDLE_BITMASK, and still works
* (behaves the same as n == 0): the input is modulo 64.
*
* \param n Integer with 0 <= n <= m4ri_radix
*/

#define __M4RI_LEFT_BITMASK(n) (m4ri_ffff >> (m4ri_radix - (n)) % m4ri_radix)

/**
* \brief create a bit mask to zero out all but the n rightmost bits.
*
* This function returns 1..64 bits, never zero bits.
* This mask is mainly used to mask the n valid bits in the least significant word
* with valid bits by using __M4RI_RIGHT_BITMASK(m4ri_radix - M->offset).
* In other words, the set bits represent the columns with the highest index in the word.
*
*  Thus,
*
*  n	Output
*  1	1000000000000000000000000000000000000000000000000000000000000000
*  2    1100000000000000000000000000000000000000000000000000000000000000
*  3    1110000000000000000000000000000000000000000000000000000000000000
*  .                                   ...
* 63	1111111111111111111111111111111111111111111111111111111111111110
* 64	1111111111111111111111111111111111111111111111111111111111111111
*
* Note that n == 0 is never passed and would fail.
*
* \param n Integer with 0 < n <= m4ri_radix
*/

#define __M4RI_RIGHT_BITMASK(n) (m4ri_ffff << (m4ri_radix - (n)))

/**
* \brief create a bit mask that is the combination of __M4RI_LEFT_BITMASK and __M4RI_RIGHT_BITMASK.
*
* This function returns 1..64 bits, never zero bits.
* This mask is mainly used to mask the n valid bits in the only word with valid bits,
* when M->ncols + M->offset <= m4ri_radix), by using __M4RI_MIDDLE_BITMASK(M->ncols, M->offset).
* It is equivalent to __M4RI_LEFT_BITMASK(n + offset) & __M4RI_RIGHT_BITMASK(m4ri_radix - offset).
* In other words, the set bits represent the valid columns in the word.
*
* Note that when n == m4ri_radix (and thus offset == 0) then __M4RI_LEFT_BITMASK is called with n == 64.
*
* \param n Integer with 0 < n <= m4ri_radix - offset
* \param offset Column offset, with 0 <= offset < m4ri_radix
*/

#define __M4RI_MIDDLE_BITMASK(n, offset) (__M4RI_LEFT_BITMASK(n) << (offset))

/**
 * \brief swap bits in the word v
 *
 * \param v The word whose bits need to be reversed.
 */

static inline word m4ri_swap_bits(word v) {
  v = ((v >>  1) & 0x5555555555555555ULL) | ((v & 0x5555555555555555ULL) << 1);
  v = ((v >>  2) & 0x3333333333333333ULL) | ((v & 0x3333333333333333ULL) << 2);
  v = ((v >>  4) & 0x0F0F0F0F0F0F0F0FULL) | ((v & 0x0F0F0F0F0F0F0F0FULL) << 4);
  v = ((v >>  8) & 0x00FF00FF00FF00FFULL) | ((v & 0x00FF00FF00FF00FFULL) << 8);
  v = ((v >> 16) & 0x0000FFFF0000FFFFULL) | ((v & 0x0000FFFF0000FFFFULL) << 16);
  v =  (v >> 32)                          |  (v                          << 32);
  return v;
}

/**
 * \brief pack bits (inverse of m4ri_spread_bits)
 *
 * \param from bitstring
 * \param Q array with bit positions
 * \param length bitsize of the output
 * \param base subtracted from every value in Q 
 *
 * \returns inverse of m4ri_spread_bits)
 *
 * \see m4ri_spread_bits
 */

static inline word m4ri_shrink_bits(word const from, rci_t* const Q, int const length, int const base) {
  word to = 0;
  switch(length-1) {
  case 15: to |= (from & (m4ri_one << (Q[15] - base))) >> (Q[15] - 15 - base);
  case 14: to |= (from & (m4ri_one << (Q[14] - base))) >> (Q[14] - 14 - base);
  case 13: to |= (from & (m4ri_one << (Q[13] - base))) >> (Q[13] - 13 - base);
  case 12: to |= (from & (m4ri_one << (Q[12] - base))) >> (Q[12] - 12 - base);
  case 11: to |= (from & (m4ri_one << (Q[11] - base))) >> (Q[11] - 11 - base);
  case 10: to |= (from & (m4ri_one << (Q[10] - base))) >> (Q[10] - 10 - base);
  case  9: to |= (from & (m4ri_one << (Q[ 9] - base))) >> (Q[ 9] -  9 - base);
  case  8: to |= (from & (m4ri_one << (Q[ 8] - base))) >> (Q[ 8] -  8 - base);
  case  7: to |= (from & (m4ri_one << (Q[ 7] - base))) >> (Q[ 7] -  7 - base);
  case  6: to |= (from & (m4ri_one << (Q[ 6] - base))) >> (Q[ 6] -  6 - base);
  case  5: to |= (from & (m4ri_one << (Q[ 5] - base))) >> (Q[ 5] -  5 - base);
  case  4: to |= (from & (m4ri_one << (Q[ 4] - base))) >> (Q[ 4] -  4 - base);
  case  3: to |= (from & (m4ri_one << (Q[ 3] - base))) >> (Q[ 3] -  3 - base);
  case  2: to |= (from & (m4ri_one << (Q[ 2] - base))) >> (Q[ 2] -  2 - base); 
  case  1: to |= (from & (m4ri_one << (Q[ 1] - base))) >> (Q[ 1] -  1 - base);
  case  0: to |= (from & (m4ri_one << (Q[ 0] - base))) >> (Q[ 0] -  0 - base);
    break;
  default:
    abort();
  }
  return to;
}

/**
 * \brief spread bits
 *
 * Given a bitstring 'from' and a spreading table Q, return a
 * bitstring where the bits of 'from' are in the positions indicated
 * by Q.
 * 
 * \param from bitstring of length 'length' stored in a word
 * \param Q table with new bit positions
 * \param length bitsize of input
 * \param base subtracted from every value in Q 
 * 
 * \returns bitstring having the same bits as from but spread using Q
 * 
 * \see m4ri_shrink_bits
 */

static inline word m4ri_spread_bits(word const from, rci_t* const Q, int const length, int const base) {
  word to = 0;
  switch(length-1) {
  case 15: to |= (from & (m4ri_one << (15))) << (Q[15]-15-base);
  case 14: to |= (from & (m4ri_one << (14))) << (Q[14]-14-base);
  case 13: to |= (from & (m4ri_one << (13))) << (Q[13]-13-base);
  case 12: to |= (from & (m4ri_one << (12))) << (Q[12]-12-base);
  case 11: to |= (from & (m4ri_one << (11))) << (Q[11]-11-base);
  case 10: to |= (from & (m4ri_one << (10))) << (Q[10]-10-base);
  case  9: to |= (from & (m4ri_one << ( 9))) << (Q[ 9]- 9-base);
  case  8: to |= (from & (m4ri_one << ( 8))) << (Q[ 8]- 8-base);
  case  7: to |= (from & (m4ri_one << ( 7))) << (Q[ 7]- 7-base);
  case  6: to |= (from & (m4ri_one << ( 6))) << (Q[ 6]- 6-base);
  case  5: to |= (from & (m4ri_one << ( 5))) << (Q[ 5]- 5-base);
  case  4: to |= (from & (m4ri_one << ( 4))) << (Q[ 4]- 4-base);
  case  3: to |= (from & (m4ri_one << ( 3))) << (Q[ 3]- 3-base);
  case  2: to |= (from & (m4ri_one << ( 2))) << (Q[ 2]- 2-base);
  case  1: to |= (from & (m4ri_one << ( 1))) << (Q[ 1]- 1-base);
  case  0: to |= (from & (m4ri_one << ( 0))) << (Q[ 0]- 0-base);
    break;
  default:
    abort();
  }
  return to;
}

/**
 * \brief Return alignment of addr w.r.t. n. For example the address
 * 17 would be 1 aligned w.r.t. 16.
 *
 * \param addr
 * \param n
 */

#define __M4RI_ALIGNMENT(addr, n) (((unsigned long)(addr))%(n))

/**
 * \brief Test for gcc >= maj.min, as per __GNUC_PREREQ in glibc
 *
 * \param maj The major version.
 * \param min The minor version.
 * \return TRUE iff we are using a GNU compile of at least version maj.min.
 */
#if defined(__GNUC__) && defined(__GNUC_MINOR__)
#define __M4RI_GNUC_PREREQ(maj, min) ((__GNUC__ << 16) + __GNUC_MINOR__ >= ((maj) << 16) + (min))
#else
#define __M4RI_GNUC_PREREQ(maj, min) FALSE
#endif

/* __builtin_expect is in gcc 3.0, and not in 2.95. */
#if __M4RI_GNUC_PREREQ(3,0) || defined(M4RI_DOXYGEN)

/**
 * \brief Macro to help with branch prediction.
 */

#define __M4RI_LIKELY(cond)    __builtin_expect ((cond) != 0, 1)

/**
 * \brief Macro to help with branch prediction.
 */

#define __M4RI_UNLIKELY(cond)  __builtin_expect ((cond) != 0, 0)

#else
#define __M4RI_LIKELY(cond)    (cond)
#define __M4RI_UNLIKELY(cond)  (cond)
#endif

/**
 * Return true if a's least significant bit is smaller than b's least significant bit.
 *
 * return true if LSBI(a) < LSBI(b),
 * where LSBI(w) is the index of the least significant bit that is set in w, or 64 if w is zero.
 *
 * \param a Word
 * \param b Word
 */

static inline int m4ri_lesser_LSB(word a, word b)
{
  uint64_t const ia = __M4RI_CONVERT_TO_UINT64_T(a);
  uint64_t const ib = __M4RI_CONVERT_TO_UINT64_T(b);
  /*
   * If a is zero then we should always return false, otherwise
   * if b is zero we should return true iff a has at least one bit set.
   */
  return !(ib ? ((ia - 1) ^ ia) & ib : !ia);
}


/**** Error Handling *****/

/**
 * \brief Print error message and abort(). 
 * 
 * The function accepts additional
 * parameters like printf, so e.g. m4ri_die("foo %d bar %f\n",1 ,2.0)
 * is valid and will print the string "foo 1 bar 2.0" before dying.
 *
 * \param errormessage a string to be printed.
 *
 * \todo Allow user to register callback which is called on
 * m4ri_die().
 *
 * \warning The provided string is not free'd.
 */

void m4ri_die(const char *errormessage, ...);

/**** IO *****/

/**
 * \brief Write a sting representing the word data to destination. 
 * 
 * \param destination Address of buffer of length at least m4ri_radix*1.3
 * \param data Source word
 * \param colon Insert a Colon after every 4-th bit. 
 * \warning Assumes destination has m4ri_radix*1.3 bytes available 
 */
void m4ri_word_to_str( char *destination, word data, int colon);

/**
 * \brief Return 1 or 0 uniformly randomly distributed.
 *
 * \todo Allow user to provide her own random() function.
 */

static inline BIT m4ri_coin_flip() {
  if (rand() < RAND_MAX/2) {
    return 0;
  }  else {
    return 1;
  }
}

/**
 * \brief Return uniformly randomly distributed random word.
 *
 * \todo Allow user to provide her own random() function.
 */

word m4ri_random_word();

/***** Initialization *****/

/**
 * \brief Initialize global data structures for the M4RI library.
 *
 * On Linux/Solaris this is called automatically when the shared
 * library is loaded, but it doesn't harm if it is called twice.
 */

#if defined(__GNUC__)
void __attribute__ ((constructor)) m4ri_init(void);
#else
void m4ri_init(void);
#endif

#ifdef __SUNPRO_C
#pragma init(m4ri_init)
#endif

/**
 * \brief De-initialize global data structures from the M4RI library. 
 *
 * On Linux/Solaris this is called automatically when the shared
 * library is unloaded, but it doesn't harm if it is called twice.
 */

#if defined(__GNUC__)
void __attribute__ ((destructor)) m4ri_fini(void);
#else
void m4ri_fini(void);
#endif

#ifdef __SUNPRO_C
#pragma fini(m4ri_fini)
#endif

/***** Memory Management *****/

/// @cond INTERNAL

#if __M4RI_CPU_L3_CACHE == 0
/*
 * Fix some standard value for L3 cache size if it couldn't be
 * determined by configure.
 */

#undef __M4RI_CPU_L3_CACHE
#if __M4RI_CPU_L2_CACHE
#define __M4RI_CPU_L3_CACHE __M4RI_CPU_L2_CACHE
#else
#define __M4RI_CPU_L3_CACHE 4194304
#endif // __M4RI_CPU_L2_CACHE
#endif // __M4RI_CPU_L3_CACHE

#if __M4RI_CPU_L2_CACHE == 0
/*
 * Fix some standard value for L2 cache size if it couldn't be
 * determined by configure.
 */
#undef __M4RI_CPU_L2_CACHE
#define __M4RI_CPU_L2_CACHE 262144
#endif // __M4RI_CPU_L2_CACHE


#if __M4RI_CPU_L1_CACHE == 0
/*
 * Fix some standard value for L1 cache size if it couldn't be
 * determined by configure.
 */
#undef __M4RI_CPU_L1_CACHE
#define __M4RI_CPU_L1_CACHE 16384
#endif // __M4RI_CPU_L1_CACHE

/// @endcond

/**
 * \brief Calloc wrapper.
 *
 * \param count Number of elements.
 * \param size Size of each element.
 *
 * \return pointer to allocated memory block.
 *
 * \todo Allow user to register calloc function.
 */

static inline void *m4ri_mm_calloc(size_t count, size_t size) {
  void *newthing;
#if __M4RI_USE_MM_MALLOC
  newthing = _mm_malloc(count * size, 64);
#elif __M4RI_USE_POSIX_MEMALIGN
  int error = posix_memalign(&newthing, 64, count * size);
  if (error) newthing = NULL;
#else
  newthing = calloc(count, size);
#endif

  if (newthing == NULL) {
    m4ri_die("m4ri_mm_calloc: calloc returned NULL\n");
    return NULL; /* unreachable. */
  }
#if __M4RI_USE_MM_MALLOC || __M4RI_USE_POSIX_MEMALIGN
  char *b = (char*)newthing;
  memset(b, 0, count * size);
#endif
  return newthing;
}

/**
 * \brief Aligned malloc wrapper.
 *
 * This function will attempt to align memory, but does not guarantee
 * success in case neither _mm_malloc nor posix_memalign are available.
 *
 * \param size Size in bytes.
 * \param alignment Alignment (16,64,...).
 *
 * \return pointer to allocated memory block.
 *
 * \todo Allow user to register malloc function.
 */

static inline void *m4ri_mm_malloc_aligned(size_t size, size_t alignment) {
  void *newthing;

#if __M4RI_USE_MM_MALLOC
  newthing = _mm_malloc(size, alignment);
#elif __M4RI_USE_POSIX_MEMALIGN
  int error = posix_memalign(&newthing, alignment, size);
  if (error)
    newthing = NULL;
#else
  newthing = malloc(size);
#endif

 if (newthing==NULL && (size>0)) {
   m4ri_die("m4ri_mm_malloc: malloc returned NULL\n");
   return NULL; /* unreachable */
 }
 else return newthing;
}

/**
 * \brief Malloc wrapper.
 *
 * \param size Size in bytes.
 *
 * \return pointer to allocated memory block.
 *
 * \todo Allow user to register malloc function.
 */

static inline void *m4ri_mm_malloc(size_t size) {
  void *newthing;
#if __M4RI_USE_MM_MALLOC
    newthing = _mm_malloc(size, 64);
#elif __M4RI_USE_POSIX_MEMALIGN
    int error = posix_memalign(&newthing, 64, size);
    if (error) newthing = NULL;
#else
    newthing = malloc(size);
#endif //__M4RI_USE_MM_MALLOC
  if (newthing==NULL && (size>0)) {
    m4ri_die("m4ri_mm_malloc: malloc returned NULL\n");
    return NULL; /* unreachable */
  }
  else return newthing;
}


/**
 * \brief Free wrapper.
 *
 * \param condemned Pointer.
 *
 * \todo Allow user to register free function.
 */

/* void m4ri_mm_free(void *condemned, ...); */
static inline void m4ri_mm_free(void *condemned, ...) { 
#if __M4RI_USE_MM_MALLOC
    _mm_free(condemned);
#else
    free(condemned);
#endif
}

/// @cond INTERNAL

/*
 * MSVC does not understand the restrict keyword 
 */

#if defined (__GNUC__)
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif



/*
 * Macros for template expansion.
 */

#define __M4RI_TEMPLATE_EXPAND0(x,y) x ## _ ## y
#define __M4RI_TEMPLATE_EXPAND1(x,y) __M4RI_TEMPLATE_EXPAND0(x,y)
#define __M4RI_TEMPLATE_NAME(fun) __M4RI_TEMPLATE_EXPAND1(fun, N)

//// @endcond

#endif // M4RI_MISC_H
