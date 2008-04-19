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

#include "grayflex.h"

#define HAVE_SSE2

#define max(x,y) ((x > y)?x:y)
#define min(x,y) ((x < y)?x:y)

#define YES 1
#define NO 0

#define LEFTMOST_BITS(w, spot)  (w & ~((1ULL<<(RADIX-spot))-1))>>(RADIX-spot)
#define RIGHTMOST_BITS(w, spot) (w &  ((1ULL<<spot)-1))

typedef unsigned char /* renamed as */ BIT;

struct matrixstruct {
  BIT *cells;
  int *rowswap;
  int *colswap;
  int nrows;
  int ncols;
};

typedef struct matrixstruct /* renamed as */ matrix;

extern BIT *table;
extern int *lookup;
extern int numcols;

void die(char *errormessage);

/* MEMLEAK, use free */
void *safeCalloc( int count, int size );

/* MEMLEAK, use free */
void *safeMalloc( int count, int size );

BIT coinFlip();

#endif //MATRIX_H
