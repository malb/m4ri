/**
 * \file misc.h
 * \brief Helper functions.
 *
 * \author Gregory Bard <bard@fordham.edu>
 * \author Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
 */

#ifndef MISC_H
#define MISC_H
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

#include <stdlib.h>

/*
 * These define entirely the word width used in the library.
 */

/**
 * A word is the typical packed data structure to represent packed
 * bits.
 */

typedef unsigned long long word;

/**
 * The number of bits in a word.
 */

#define RADIX (sizeof(word)<<3)

/**
 * The number one as a word.
 */

#define ONE ((word)1)


/**
 * \brief Return the maximal element of x and y
 *
 * \param x Word
 * \param y Word
 */

#ifndef MAX
#define MAX(x,y) ((x > y)?x:y)
#endif

/**
 * \brief Return the minimal element of x and y
 *
 * \param x Word
 * \param y Word
 */

#ifndef MIN
#define MIN(x,y) ((x < y)?x:y)
#endif

/**
 * \brief Return r such that x elements fit into r blocks of length y.
 *
 * \param x Number of elements
 * \param y Block size
 */

#define DIV_CEIL(x,y) ((x%y)?x/y+1:x/y)

/**
 *\brief Pretty for 1.
 */ 

#define TRUE 1

/**
 *\brief Pretty for 0.
 */ 

#define FALSE 0

/**
 * \brief 2^i
 *
 * \param i Integer.
 */ 

#define TWOPOW(i) (1<<(i))

/**
 * Pretty for unsigned char.
 */

typedef unsigned char BIT;

/**
 * \brief Clear the bit spot in the word w
 * 
 * \param w Word
 * \param spot Integer with 0 <= spot < RADIX
 */

#define CLR_BIT(w, spot) (w &= ~(ONE<<(RADIX - spot - 1)))

/**
 * \brief Set the bit spot in the word w
 * 
 * \param w Word
 * \param spot Integer with 0 <= spot < RADIX
 */

#define SET_BIT(w, spot) (w |= (ONE<<(RADIX - spot - 1)))

/**
 * \brief Get the bit spot in the word w
 * 
 * \param w Word
 * \param spot Integer with 0 <= spot < RADIX
 */

#define GET_BIT(w, spot) ((w & (ONE<<(RADIX - spot - 1))) >> (RADIX - spot - 1))

/**
* \brief Return the n leftmost bits of the word w.
*
* \param w Word
* \param n Integer with 0 <= spot < RADIX
*/

#define LEFTMOST_BITS(w, n)  (w & ~((ONE<<(RADIX-n))-1))>>(RADIX-n)

/**
* \brief Return the n rightmost bits of the word w.
*
* \param w Word
* \param n Integer with 0 <= spot < RADIX
*/

#define RIGHTMOST_BITS(w, n) (w &  ((ONE<<n)-1))

/**
 * \brief return alignment of addr w.r.t. n. For example the address
 * 17 would be 1 aligned w.r.t. 16.
 *
 * \param addr
 * \param n
 */

#define ALIGNMENT(addr, n) (((unsigned long)(addr))%n)

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

void m4ri_die(char *errormessage, ...);

/**** IO *****/

/**
 * \brief Write a sting representing the word data to destination. 
 * 
 * \param destination Address of buffer of length at least RADIX*1.3
 * \param data Source word
 * \param colon Insert a Colon after every 4-th bit. 
 * \warning Assumes destination has RADIX*1.3 bytes available 
 */
void m4ri_word_to_str( char *destination, word data, int colon);

/***** Memory Management *****/

/**
 * \brief Calloc wrapper.
 *
 * \param count Number of elements.
 * \param size Size of each element.
 *
 * \todo Allow user to register calloc function.
 */

void *m4ri_mm_calloc( int count, int size );

/**
 * \brief Malloc wrapper.
 *
 * \param size Size in bytes.
 *
 * \todo Allow user to register malloc function.
 */

void *m4ri_mm_malloc( int size );

/**
 * \brief Free wrapper.
 *
 * \param condemned Pointer.
 *
 * \todo Allow user to register free function.
 */

void m4ri_mm_free(void *condemned);

/**
 * \brief Return 1 or 0 uniformly randomly distributed.
 *
 * \todo Allow user to provide her own random() function.
 */

BIT m4ri_coin_flip();

#endif //MISC_H
