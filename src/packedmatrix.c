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

#include <stdlib.h>
#include <string.h>
#include "packedmatrix.h"

word packingmask[RADIX];
word bytemask[RADIX/8];
word sixteenmask[RADIX/16];

/****************************************/

packedmatrix *m2t_init(int r, int c) {
  packedmatrix *newmatrix;
  int i;

  newmatrix=(packedmatrix *)calloc(1, sizeof(packedmatrix));

  if ((c % RADIX)==0) 
    newmatrix->width=(c/RADIX);
  else 
    newmatrix->width=(c/RADIX) + 1;

  newmatrix->ncols=c;
  newmatrix->nrows=r;

  newmatrix->values=(word *)safeCalloc( (newmatrix->width)*r, sizeof(word));

  newmatrix->rowswap=(int *)safeMalloc( r, sizeof(int));

  // Rowswap does not contain the rowswap index i but the correct
  // offset in the values table. Rowswap is exclusively used to access
  // elements in that table and this speeds up computation a little. (malb)

  for (i=0; i<r; i++) { newmatrix->rowswap[i]=i*(newmatrix->width); }

  return newmatrix;
}

/** We don't perform any sanity checks! **/
packedmatrix *m2t_init_window(packedmatrix *m, int lowr, int lowc, int highr, int highc) {
  int nrows, ncols, i, offset; 
  packedmatrix *window = (packedmatrix *)calloc(1, sizeof(packedmatrix));
  nrows = highr - lowr;
  ncols = highc - lowc;
  
  window->ncols = ncols;
  window->nrows = nrows;
  window->width = ncols/RADIX;
  window->values = m->values;
  window->rowswap = (int *)safeCalloc( nrows, sizeof(int));

  offset = lowc / RADIX;

  for(i=0; i<nrows; i++) {
    window->rowswap[i] = m->rowswap[lowr + i] + offset;
  }
  
  return window;
}

void m2t_free( packedmatrix *condemned) {
  free(condemned->values);
  free(condemned->rowswap);
  free(condemned);
}

void m2t_free_window( packedmatrix *condemned) {
  free(condemned->rowswap);
  free(condemned);
}

/************************************************************/

/* Internal: do not call */
void setupPackingMasks() {
  int i, j;
  word x=1;

  for (i=RADIX-1; i>=0; i--) {
    packingmask[i]=x;
    x<<=1;
  }

  for (i=0; i<RADIX/8; i++) {
    x=0;
    for (j=0; j<8; j++) {
      x|=packingmask[j+i*8];
    }
    bytemask[i]=x; 
  }

  for (i=0; i<RADIX/16; i++) {
    x=0;
    for (j=0; j<16; j++) {
      x|=packingmask[j+i*16];
    }
    sixteenmask[i]=x; 
  }
} 

/**********************************************************************/
/* The RADIX-bit word can be divided into RADIX/8 octets or 8 bit bytes */
/* The most significant byte is numbered 0. */
word fetchByte( word data, int which ) {
  word masked=data & bytemask[which];
  int chunks=RADIX/8;
  int moved=(chunks-which-1)*8;

  masked>>=moved;

  return masked;
}

/**********************************************************************/

/* Warning: I assume *destination has RADIX+1 bytes available */

void wordToString( char *destination, word data) {
  int i;

  for (i=0; i<RADIX; i++) {
    if ((data & packingmask[i]) == 0) {
      destination[i]='0';
    } else destination[i]='1';
  }
  
  destination[RADIX]='\0';
}

/* Warning: I assume *destination has 9 bytes available */

void byteToString( char *destination, word data) {
  int i;

  for (i=0; i<8; i++) {
    if ((data & packingmask[i+RADIX-8]) == 0) {
      destination[i]='0';
    } else destination[i]='1';
  }
  
  destination[8]='\0';
}

/**********************************************************************/

/* Warning: I assume *destination has RADIX*1.25 bytes available */
void wordToStringComma( char *destination, word data) {
  int i;
  int j=0;
  
  for (i=0; i<RADIX; i++) {
    if ((data & packingmask[i]) == 0) {
      destination[j]='0';
    } else destination[j]='1';
    j++;
      
    if (((i % 4)==3) && (i!=RADIX-1)) {
      destination[j]=':';
      j++;
    }
  }
  
  destination[(int)(RADIX*1.25)-1]='\0';
}
  
/**********************************************************************/

/* Works */
void printMatrix( packedmatrix *m ) {
  int i, j;
  char temp[SAFECHAR];
  word block;

  for (i=0; i< m->nrows; i++ ) {
    printf("[ ");

    for (j=0; j< m->ncols; j+=RADIX) {
      block=readBlock(m, i, j);
      wordToStringComma(temp, block);
      printf("%s ", temp);
    }
    printf("]\n");
  }
}

/**********************************************************************/

/* Works */
void printMatrixTight( packedmatrix *m ) {
  int i, j;
  char temp[SAFECHAR];
  word block;

  for (i=0; i< m->nrows; i++ ) {
    printf("[");

    for (j=0; j< m->ncols; j+=RADIX) {
      block=readBlock(m, i, j);
      wordToString(temp, block);
      printf("%s", temp);
    }
    printf("]\n");
  }

  printf("\n\n\n");
}

/**********************************************************************/
/* this clears the row, but only begins at the column coloffset */
void rowClearOffset(packedmatrix *m, int row, int coloffset) {
  int startblock= coloffset/RADIX;
  int i;
  word temp;
  
  // make sure to start clearing at coloffset
  if (coloffset%RADIX) {
    temp=readBlock(m, row, coloffset);
    temp &=  ~(((ONE<<(RADIX-coloffset%RADIX))) - ONE);
  } else {
    temp = 0;
  }
  writeBlock(m, row, coloffset, temp);

  temp=0;

  for ( i=startblock+1; i < (m->width); i++ ) {
    writeBlock(m, row, i*RADIX, temp);
  }
}


/**********************************************************************/
/* this adds rows sourcerow and destrow and stores the total in row
   destrow, but only begins at the column coloffset */

void rowAddOffset( packedmatrix *m, int sourcerow, int destrow, 
		   int coloffset ) {

  int startblock= coloffset/RADIX;
  int i;
  word temp;
  
  // make sure to start adding at coloffset
  temp=readBlock(m, sourcerow, startblock*RADIX);
  if (coloffset%RADIX)
    temp &= (ONE<<(RADIX - (coloffset%RADIX))) - ONE;
  xorBlock(m, destrow, startblock*RADIX, temp);

  for ( i=startblock+1; i < (m->width); i++ ) {
    temp=readBlock(m, sourcerow, i*RADIX);
    xorBlock(m, destrow, i*RADIX, temp);
  }
}


void rowAdd( packedmatrix *m, int sourcerow, int destrow) {
  rowAddOffset(m, sourcerow, destrow, 0);
}

/**********************************************************************/
/* This will do Gaussian Elimination on the matrix m but will start not
 at column 0 necc but at column "startcol". If full=NO, then it will do
 triangular style elimination, and if full=YES, it will do Gauss-Jordan style,
 or full elimination.*/

int reduceGaussianDelayed(packedmatrix *m, int startcol, int full) {
  int i,j;
  int start; 

  int startrow = startcol;
  int ii;
  int pivots = 0;
  for (i=startcol ; i<m->ncols ; i++) {

    for(j=startrow ; j < m->nrows; j++) {
      if (readCell(m,j,i)) {
	rowSwap(m,startrow,j);
	pivots++;

	if (full==YES) start=0; else start=i+1;

	for(ii=start ;  ii < m->nrows ; ii++) {
	  if (ii != startrow) {
	    if (readCell(m, ii, i)) {
	      rowAddOffset(m, startrow, ii, i);
	    }
	  }
	}
	startrow = startrow + 1;
	break;
      }
    }
  }

  return pivots;
}


/* This will do Gaussian Elimination on the matrix m. 
   If full=NO, then it will do
 triangular style elimination, and if full=YES, it will do Gauss-Jordan style,
 or full elimination.*/
int reduceGaussian(packedmatrix *m, int full) { 
  return reduceGaussianDelayed(m,0, full); 
}


packedmatrix *m2t_transpose(packedmatrix *newmatrix, packedmatrix *data) {
  int i,j,k;
  word temp;

  if (newmatrix == NULL) {
    newmatrix = m2t_init( data->ncols, data->nrows );
  } else {
    if (newmatrix->nrows != data->ncols || newmatrix->ncols != data->nrows) {
      die("Wrong size for return matrix.\n");
    }
  }

  for (i=0; i<newmatrix->nrows; i++) {
    for (j=0; j<newmatrix->width; j++) {
      temp=(word)0;
      for (k=0; k<RADIX; k++) {
	if (  (j*RADIX+k) < data->nrows ) { 
	  if (readCell(data, j*RADIX+k, i)==1)
	    temp=temp | packingmask[k]; 
	}
      }
      writeBlock(newmatrix, i, j*RADIX, temp);
    }
  }
	
  return newmatrix;
}

/********************************************************/

/* Works */
/* Internal to naive matrix mult */
BIT bigDotProduct( packedmatrix *a, packedmatrix *bT, int rowofa,
			 int rowofb ) {
  /* ``a slot'' is a row of A, and a column of B when calcing AB */
  /* but since we use B^T so that we are working only with rows, */
  /* ``a slot'' of A is a row, ``a slot'' of B is a row of B^T */
  int total, i;

  total=0;
  for (i=0; i< a->width; i++) {
    if (  (i*RADIX) < a->nrows )  
    total+=dotProduct( readBlock(a, rowofa, i*RADIX), 
		       readBlock(bT, rowofb, i*RADIX) );
  }

  return (BIT)(total % 2);
}

/********************************************************/

packedmatrix *m2t_mul_naiv_t(packedmatrix *product, packedmatrix *a, 
			     packedmatrix *bT ) {
  int i, j;
  int newrows=a->nrows;
  int newcols=bT->nrows;

  if (product==NULL) {
    product=m2t_init(newrows, newcols);
  } else {
    if (product->nrows != newrows || product->ncols != newcols) {
      die("Provided return matrix has wrong dimensions.\n");
    }
  }

  for (i=0; i<newrows; i++) {
    for (j=0; j<newcols; j++) {
      writeCell(product, i, j, bigDotProduct( a, bT, i, j ) );
    }
  }

  return product;
}
  
packedmatrix *m2t_mul_naiv(packedmatrix *product, packedmatrix *a, 
			   packedmatrix *b) {
  packedmatrix *bT = m2t_transpose(NULL, b);
  product = m2t_mul_naiv_t(product, a, bT );
  m2t_free(bT);

  return product;
}

/********************************************************/

/* Works */
word randomWord() {
  int i;
  word temp=0;

  for (i=0; i<RADIX; i++) {
    if (coinFlip()==1) {
      temp|=packingmask[i];
    }
  }

  return temp;
}

/********************************************************/

/* Works */
void fillRandomly( packedmatrix *a ) {
  int i, j;
  
  for (i=0; i < (a->nrows); i++) {
    for (j=0; j < (a->ncols); j++) {
      writeCell(a, i, j, coinFlip() );
    }
  }
}

/********************************************************/

/* Works */
void makeIdentity( packedmatrix *a ) {
  int i,j;
  int stop=min(a->nrows, a->ncols);

  for (i=0; i< (a->nrows); i++) {
    for (j=0; j< (a->width); j++) {
      
      writeBlock(a, i, j*RADIX, 0);
    }
  }

  for (i=0; i<stop; i++) {
    writeCell(a, i, i, 1);
  }
}

/********************************************************/

/********************************************************/

/* Works */
BIT equalMatrix( packedmatrix *a, packedmatrix *b ) {
  int i, j;
  word block1, block2;

  if (a->nrows!=b->nrows) return NO;
  if (a->ncols!=b->ncols) return NO;

  for (i=0; i< a->nrows; i++) {
    for (j=0; j< a->width; j++) {
      block1=readBlock(a, i, j*RADIX);
      block2=readBlock(b, i, j*RADIX);
      if (isEqualWord(block1, block2)==NO) return NO;
    }
  }
  return YES;
}

int compareMatrix(packedmatrix *a, packedmatrix *b) {

  int i,j;

  if(a->nrows < b->nrows) return -1;
  if(b->nrows < a->nrows) return 1;
  if(a->ncols < b->ncols) return -1;
  if(b->ncols < a->ncols) return 1;

  for(i=0; i < a->nrows ; i++) {
    for(j=0 ; j< a->width ; j++) {
      if ( a->values[a->rowswap[i] + j] < b->values[b->rowswap[i] + j])
	return -1;
      else if( a->values[a->rowswap[i] + j] > b->values[b->rowswap[i] + j])
	return 1;
    }
  }
  return 0;
}

/********************************************************/

/* MEMLEAK: use m2t_free */
packedmatrix *cloneMatrix( packedmatrix *p) {
  packedmatrix *n = m2t_init(p->nrows, p->ncols);
  int i, j, p_truerow, n_truerow;
  
  for (i=0; i<p->nrows; i++) {
    p_truerow = p->rowswap[i];
    n_truerow = n->rowswap[i];
    for (j=0; j<p->width; j++) {
      n->values[n_truerow + j] = p->values[p_truerow + j];
    }
  }

  return n;
}


/********************************************************/

/* MEMLEAK: use m2t_free */
/* This is sometimes called augment */
packedmatrix *concat( packedmatrix *a, packedmatrix *b) {
  packedmatrix *newmatrix;
  int i, j, truerow, width;
  //word entry;
  
  if (a->nrows!=b->nrows) {
    die("Bad arguments to concat!\n");
  }

  newmatrix=m2t_init(a->nrows, a->ncols + b->ncols);
  width = newmatrix->width;

  for (i=0; i<a->nrows; i++) {
    truerow = a->rowswap[i];
    for (j=0; j<a->width; j++) {
      newmatrix->values[i*width + j] = a->values[truerow + j];
    }
  }

  for (i=0; i<b->nrows; i++) {
    for (j=0; j<b->ncols; j++) {
      writeCell(newmatrix, i, j+(a->ncols), 
		readCell(b, i, j) );
    }
  }

  return newmatrix;
}

packedmatrix *stack(packedmatrix *a, packedmatrix *b) {
  packedmatrix *newmatrix;
  int i, j, offset, truerow, width;

  if (a->ncols != b->ncols) {
    die("Bad arguments to stack!\n");
  }
  newmatrix = m2t_init(a->nrows + b->nrows, a->ncols);
  width = newmatrix->width;
  
  for(i=0; i<a->nrows; i++) {
    truerow = a->rowswap[i];
    for (j=0; j<a->width; j++) {
      newmatrix->values[i*width + j] = a->values[truerow + j]; 
    }
  }
  offset = a->nrows * a->width;

  for(i=0; i<b->nrows; i++) {
    truerow = b->rowswap[i];
    for (j=0; j<b->width; j++) {
      newmatrix->values[offset + i*width + j] = b->values[truerow + j]; 
    }
  }
  return newmatrix;
}

/*********************************************************/
/* MEMLEAK: use m2t_free */
packedmatrix *invertGaussian(packedmatrix *target, 
				   packedmatrix *identity) {
  packedmatrix *huge, *inverse;
  int x;

  huge=concat(target, identity);

  x=reduceGaussian(huge, YES);

  if (x==NO) { m2t_free(huge); return NULL; }
  
  inverse=copySubMatrix(huge, 0, target->ncols, target->nrows, 
			target->ncols*2);

  m2t_free(huge);
  return inverse;
}

packedmatrix *m2t_add(packedmatrix *ret, packedmatrix *left, packedmatrix *right) {
  if (left->nrows != right->nrows || left->ncols != right->ncols) {
    die("rows and columns must match");
  }
  if (ret == NULL) {
    ret = cloneMatrix(left);
  } else if (ret != left) {
    if (ret->nrows != left->nrows || ret->ncols != left->ncols) {
      die("rows and columns of returned matrix must match");
    }
  }
  return _m2t_add_impl(ret, left, right);
}

packedmatrix *_m2t_add_impl(packedmatrix *ret, packedmatrix *left, packedmatrix *right) {
  int i,j,left_truerow, ret_truerow;

  if (ret != left) {
    for (i=0; i<left->nrows; i++) {
      left_truerow = left->rowswap[i];
      ret_truerow = ret->rowswap[i];
      for (j=0; j < left->width; j++) {
	ret->values[ret_truerow + j] = left->values[left_truerow + j];
      }
    }
  }
  
  for(i=0; i < right->nrows; i++) {
    for(j=0; j < left->width; j++) {
      ret->values[  ret->rowswap[i] + j ]  ^= right->values[ right->rowswap[i] + j];
    }
  }
  return ret;
}

packedmatrix *copySubMatrix(packedmatrix *m, int startrow, int startcol, int endrow, int endcol) {
  int nrows, ncols, truerow, i, colword, x, y, block, spot, startword;
  word temp  = 0;
  
  nrows = endrow - startrow;
  ncols = endcol - startcol;

  packedmatrix *newmatrix = m2t_init(nrows, ncols);

  startword = startcol / RADIX;

  /** we start at the beginning of a word **/
  if (startcol%RADIX == 0) {
    for(x = startrow, i=0; i<nrows; i++, x+=1) {
      truerow = m->rowswap[x];

      /** process full words first **/
      for(y = startcol, colword=0; colword<ncols/RADIX; colword++, y+=RADIX) {
	block = truerow + colword + startword;
	temp = m->values[block];
	newmatrix->values[newmatrix->rowswap[i] + colword] = temp;
      }

      /** process remaining bits **/
      if (ncols%RADIX) {
	colword = ncols/RADIX;
	block = truerow + colword;
	temp = m->values[block] & ~((1ULL<<(RADIX-ncols%RADIX))-1);
	newmatrix->values[newmatrix->rowswap[i] + colword] = temp;
      } 
    }

    /** startcol is not the beginning of a word **/
  } else { 
    spot = startcol % RADIX;
    for(x = startrow, i=0; i<nrows; i++, x+=1) {
      truerow = m->rowswap[x];

      /** process full words first **/
      for(y = startcol, colword=0; colword<ncols/RADIX; colword++, y+=RADIX) {
	block = truerow + colword + startword;
	temp = (m->values[block] << (spot)) | (m->values[block + 1] >> (RADIX-spot) ); 
	newmatrix->values[newmatrix->rowswap[i] + colword] = temp;
      }
      /** process remaining bits (lazy)**/
      colword = ncols/RADIX;
      for (y=0; y < ncols%RADIX; y++) {
	temp = readCell(m, x, startcol + colword*RADIX + y);
	writeCell(newmatrix, i, colword*RADIX + y, temp);
      }
    }
  }
  return newmatrix;
}
