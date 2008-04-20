#ifndef MATRIX_H
#define MATRIX_H

/******************************************************************************
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
******************************************************************************/

//#include "grayflex.h"
#include <stdlib.h>

typedef unsigned long long word;

#define RADIX 64
#define MAX(x,y) ((x > y)?x:y)
#define MIN(x,y) ((x < y)?x:y)
#define DIV_CEIL(x,y) ((x%y)?x/y+1:x/y)

#define TRUE 1
#define FALSE 0

#define TWOPOW(i) (1<<(i))

typedef unsigned char BIT;

/**** bitmasks *****/

extern word packingmask[RADIX]; //deprecated
extern word bytemask[RADIX/8]; //deprecated
extern word sixteenmask[RADIX/16]; //deprecated

/**
 * Clear the bit spot in the word w
 * 
 * @param w Word
 * @param spot Integer with 0 <= spot < RADIX
 */

#define CLR_BIT(w, spot) (w &= ~(1ULL<<(RADIX - spot - 1)))
#define SET_BIT(w, spot) (w |= (1ULL<<(RADIX - spot - 1)))
#define GET_BIT(w, spot) ((w & (1ULL<<(RADIX - spot - 1))) >> (RADIX - spot - 1))
#define LEFTMOST_BITS(w, spot)  (w & ~((1ULL<<(RADIX-spot))-1))>>(RADIX-spot)
#define RIGHTMOST_BITS(w, spot) (w &  ((1ULL<<spot)-1))

void m4ri_setup_packing_masks();

/**** Error Handling *****/

/**
 * Print error message and abort().
 *
 * \param errormessage a string to be printed
 *
 * \todo Allow user to register callback.
 *
 * \warning The provided string is not free'd.
 */

void m4ri_die(char *errormessage);

void m4ri_print_bit_string(int number, int length);

/***** Memory Management *****/

/**
 * \brief Calloc wrapper.
 *
 * \param count Number of elements.
 * \param size Size of each element.
 */

void *m4ri_mm_calloc( int count, int size );

/**
 * \brief Malloc wrapper.
 *
 * \param size Size in bytes.
 */

void *m4ri_mm_malloc( int size );

/**
 * \brief Free wrapper.
 *
 * \param condemned Pointer.
 */

static inline void m4ri_mm_free(void *condemned) { 
  free(condemned); 
};

BIT m4ri_coin_flip();

#endif //MATRIX_H
