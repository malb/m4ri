#ifndef PACKEDMATRIX_H
#define PACKEDMATRIX_H

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

#include "misc.h"
#include <stdio.h>

#define RADIX 64
#define SAFECHAR 85 /* Radix*1.25 plus a few. */
#define ONE 1ULL


typedef unsigned long long word;

struct packedmatrixstruct {
  word *values;

  int nrows;
  int ncols;
  int width; /* rounded as floor(rows/RADIX)+1. */

  int *rowswap;

};

typedef struct packedmatrixstruct packedmatrix;

extern word packingmask[RADIX];
extern word bytemask[RADIX/8];
extern word sixteenmask[RADIX/16];

packedmatrix *createMatrix(int r, int c);

void destroyMatrix( packedmatrix *condemned );

static inline void rowSwap( packedmatrix *m, int rowa, int rowb ) {
  int temp=m->rowswap[rowa];
  m->rowswap[rowa]=m->rowswap[rowb];
  m->rowswap[rowb]=temp;
}

/* Internal: do not call */
void setupPackingMasks();

static inline BIT readCell( packedmatrix *m, int row, int col ) {
  int block=col/RADIX;
  int spot=col % RADIX;
  int truerow=m->rowswap[row];
  
  word entry=m->values[ block + truerow ];
 
  word resolved=entry & ((ONE)<<(RADIX - spot - 1));
  
  return (resolved >> (RADIX - spot -1));
}

static inline void writeCell( packedmatrix *m, int row, int col, BIT value) {
  int block=col/RADIX;
  int spot=col % RADIX;
  int truerow=m->rowswap[row];

  if (value==0) {
    m->values[ block + truerow ] &= ~((ONE) <<(RADIX - spot - 1));
  } else {
    m->values[ block + truerow ] |= ((ONE)<<(RADIX - spot - 1));
  }
}


/* Keep in mind that the row, col refer to a row and column (of bits), and
   you can address the block by any of the RADIX (usually 64) A_ijs there. */
static inline void xorBlock( packedmatrix *m, int row, int col, word value) {
  int block=col/RADIX;
  int truerow=m->rowswap[row];

  word *entry=m->values + block + truerow;
  *entry ^= value;
}

/* Keep in mind that the row, col refer to a row and column (of bits), and
   you can address the block by any of the RADIX (usually 64) A_ijs there. */
static inline void writeBlock( packedmatrix *m, int row, int col, word value) {
  int block=col/RADIX;
  int truerow=m->rowswap[row];

  m->values[ block + truerow] = value;
}



word fetchByte( word data, int which );

/* The RADIX-bit word can be divided into RADIX/8 octets or 8 bit bytes */
/* The most significant byte is numbered 0. */
word fetch16( word data, int which );

/* Warning: I assume *destination has RADIX+1 bytes available */

void wordToString( char *destination, word data);

/* Warning: I assume *destination has 9 bytes available */

void byteToString( char *destination, word data);

/* Warning: I assume *destination has RADIX*1.25 bytes available */

void wordToStringComma( char *destination, word data);
  
/* Important note: You can specify any of the RADIX bits (64 bits usually),
   inside of the block, and it will still return the correct entire block */
/* Important note: You can specify any of the RADIX bits (64 bits usually),
   inside of the block, and it will still return the correct entire block */
static inline word readBlock( packedmatrix *m, int row, int col ) {
  int block=col/RADIX;
  int truerow=m->rowswap[row];

  word entry=m->values[ block + truerow ];
  
  return entry;
}

void printMatrix( packedmatrix *m );

void printMatrixTight( packedmatrix *m );

/**********************************************************************/
/* this adds rows sourcerow and destrow and stores the total in row
   destrow, but only begins at the column coloffset */

void rowAddOffset( packedmatrix *m, int sourcerow, int destrow, 
			 int coloffset );

void rowClearOffset(packedmatrix *m, int row, int coloffset);

void rowAdd( packedmatrix *m, int sourcerow, int destrow);

int reduceGaussianDelayed(packedmatrix *m, int startcol, int full);

int reduceGaussian(packedmatrix *m, int full);

static inline BIT dotProduct( word a, word b ) {
  word temp=a & b;
  //int i, 
  int total=0;

/*   for (i=0; i<RADIX; i++) { */
/*     if ((temp & packingmask[i])!=0) total++; */
/*   } */
/*   return (total % 2); */
  while (temp)  {
    total = !total;
    temp = temp & (temp - 1);
  }
  return total;
}

packedmatrix *transpose( packedmatrix *data );

BIT bigDotProduct( packedmatrix *a, packedmatrix *bT, int rowofa,
			 int rowofb );

/* MEMLEAK use destroyMatrix */

packedmatrix *matrixTimesMatrixTranspose( packedmatrix *a, 
						packedmatrix *bT );
  
/* MEMLEAK: use destroyMatrix */
/* Normally, if you will multiply several times by b, it is smarter to
  calculate bT yourself, and keep it, and then use the function called
  matrixTimesMatrixTranspose */
packedmatrix *matrixTimesMatrix( packedmatrix *a, 
				       packedmatrix *b);
word randomWord();

void fillRandomly( packedmatrix *a );

void makeIdentity( packedmatrix *a );

matrix *unpackMatrix( packedmatrix *pm );

packedmatrix *packMatrix( matrix *m );

static inline BIT isEqualWord( word a, word b) {
  if (a==b) return YES;
  else return NO;
}

BIT equalMatrix( packedmatrix *a, packedmatrix *b );

int compareMatrix(packedmatrix *a, packedmatrix *b);

/* MEMLEAK: use destroyMatrix */
packedmatrix *clone( packedmatrix *p);

/* MEMLEAK: use destroyMatrix */
/* This is sometimes called augment */
packedmatrix *concat( packedmatrix *a, packedmatrix *b);


/* MEMLEAK: use destroyMatrix */
packedmatrix *copySubMatrix( packedmatrix *a, int lowr, int lowc,
				   int highr, int highc);

/* MEMLEAK: use destroyMatrix */
packedmatrix *invertGaussian(packedmatrix *target, 
				   packedmatrix *identity);

/* MEMLEAK: use destroyMatrix */
packedmatrix *add(packedmatrix *left, packedmatrix *right);

void lazyPrint(packedmatrix *a);


#endif //PACKEDMATRIX_H
