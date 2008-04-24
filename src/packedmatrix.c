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

static BIT mzd_big_dot_product( packedmatrix *a, packedmatrix *bT, int rowofa, int rowofb );
static inline BIT mzd_dot_product( word a, word b );

packedmatrix *mzd_init(int r, int c) {
  packedmatrix *newmatrix;
  int i;

  newmatrix=(packedmatrix *)m4ri_mm_calloc(1, sizeof(packedmatrix));

  if ((c % RADIX)==0) 
    newmatrix->width=(c/RADIX);
  else 
    newmatrix->width=(c/RADIX) + 1;

  newmatrix->ncols=c;
  newmatrix->nrows=r;

  newmatrix->values=(word *)m4ri_mm_calloc( (newmatrix->width)*r, sizeof(word));

  newmatrix->rowswap=(int *)m4ri_mm_malloc( r * sizeof(int));

  /* Rowswap does not contain the rowswap index i but the correct
   * offset in the values table. Rowswap is exclusively used to access
   * elements in that table and this speeds up computation a
   * little.
   */

  for (i=0; i<r; i++) { 
    newmatrix->rowswap[i]=i*(newmatrix->width); 
  }

  return newmatrix;
}

/* We don't perform any sanity checks! */
packedmatrix *mzd_init_window(packedmatrix *m, int lowr, int lowc, int highr, int highc) {
  int nrows, ncols, i, offset; 
  packedmatrix *window = (packedmatrix *)m4ri_mm_calloc(1, sizeof(packedmatrix));
  nrows = MIN(highr - lowr, m->nrows - lowr);
  ncols = highc - lowc;
  
  window->ncols = ncols;
  window->nrows = nrows;
  window->width = ncols/RADIX;
  if (ncols%RADIX)
    window->width++;
  window->values = m->values;
  window->rowswap = (int *)m4ri_mm_calloc( nrows, sizeof(int));

  offset = lowc / RADIX;

  for(i=0; i<nrows; i++) {
    window->rowswap[i] = m->rowswap[lowr + i] + offset;
  }
  
  return window;
}

void mzd_free( packedmatrix *condemned) {
  m4ri_mm_free(condemned->values);
  m4ri_mm_free(condemned->rowswap);
  m4ri_mm_free(condemned);
}

void mzd_free_window( packedmatrix *condemned) {
  m4ri_mm_free(condemned->rowswap);
  m4ri_mm_free(condemned);
}

/* Warning: I assume *destination has RADIX*1.25 bytes available */
void wordToStringComma( char *destination, word data) {
  int i;
  int j=0;
  
  for (i=0; i<RADIX; i++) {
    if (GET_BIT(data,i) == 0) {
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
  
void mzd_print_matrix( packedmatrix *m ) {
  int i, j;
  char temp[SAFECHAR];
  word block;

  for (i=0; i< m->nrows; i++ ) {
    printf("[ ");

    for (j=0; j< m->ncols; j+=RADIX) {
      block=mzd_read_block(m, i, j);
      m4ri_word_to_str(temp, block, 1);
      printf("%s ", temp);
    }
    printf("]\n");
  }
}

void mzd_print_matrix_tight( packedmatrix *m ) {
  int i, j;
  char temp[SAFECHAR];
  word block;

  for (i=0; i< m->nrows; i++ ) {
    printf("[");

    for (j=0; j< m->ncols; j+=RADIX) {
      block=mzd_read_block(m, i, j);
      m4ri_word_to_str(temp, block, 0);
      printf("%s", temp);
    }
    printf("]\n");
  }

  printf("\n\n\n");
}

void mzd_row_clear_offset(packedmatrix *m, int row, int coloffset) {
  int startblock= coloffset/RADIX;
  int i;
  word temp;
  
  // make sure to start clearing at coloffset
  if (coloffset%RADIX) {
    temp=mzd_read_block(m, row, coloffset);
    temp &=  ~(((ONE<<(RADIX-coloffset%RADIX))) - ONE);
  } else {
    temp = 0;
  }
  mzd_write_block(m, row, coloffset, temp);

  temp=0;

  for ( i=startblock+1; i < (m->width); i++ ) {
    mzd_write_block(m, row, i*RADIX, temp);
  }
}


void mzd_row_add_offset( packedmatrix *m, int sourcerow, int destrow, 
		   int coloffset ) {

  int startblock= coloffset/RADIX;
  int i;
  word temp;
  
  // make sure to start adding at coloffset
  temp=mzd_read_block(m, sourcerow, startblock*RADIX);
  if (coloffset%RADIX)
    temp &= (ONE<<(RADIX - (coloffset%RADIX))) - ONE;
  mzd_xor_block(m, destrow, startblock*RADIX, temp);

  for ( i=startblock+1; i < (m->width); i++ ) {
    temp=mzd_read_block(m, sourcerow, i*RADIX);
    mzd_xor_block(m, destrow, i*RADIX, temp);
  }
}


void mzd_row_add( packedmatrix *m, int sourcerow, int destrow) {
  mzd_row_add_offset(m, sourcerow, destrow, 0);
}

int mzd_gauss_delayed(packedmatrix *m, int startcol, int full) {
  int i,j;
  int start; 

  int startrow = startcol;
  int ii;
  int pivots = 0;
  for (i=startcol ; i<m->ncols ; i++) {

    for(j=startrow ; j < m->nrows; j++) {
      if (mzd_read_bit(m,j,i)) {
	mzd_row_swap(m,startrow,j);
	pivots++;

	if (full==TRUE) start=0; else start=i+1;

	for(ii=start ;  ii < m->nrows ; ii++) {
	  if (ii != startrow) {
	    if (mzd_read_bit(m, ii, i)) {
	      mzd_row_add_offset(m, startrow, ii, i);
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

int mzd_reduce_naiv(packedmatrix *m, int full) { 
  return mzd_gauss_delayed(m,0, full); 
}

packedmatrix *mzd_transpose(packedmatrix *newmatrix, packedmatrix *data) {
  int i,j,k;
  word temp;

  if (newmatrix == NULL) {
    newmatrix = mzd_init( data->ncols, data->nrows );
  } else {
    if (newmatrix->nrows != data->ncols || newmatrix->ncols != data->nrows) {
      m4ri_die("Wrong size for return matrix.\n");
    }
  }

  for (i=0; i<newmatrix->nrows; i++) {
    for (j=0; j<newmatrix->width; j++) {
      temp=(word)0;
      for (k=0; k<RADIX; k++) {
	if (  (j*RADIX+k) < data->nrows ) { 
	  if (mzd_read_bit(data, j*RADIX+k, i)==1)
	    SET_BIT(temp,k);
	}
      }
      mzd_write_block(newmatrix, i, j*RADIX, temp);
    }
  }
	
  return newmatrix;
}

static inline BIT _mzd_dot_product( word a, word b ) {
  word temp=a & b;
  //int i, 
  int total=0;

  while (temp)  {
    total = !total;
    temp = temp & (temp - 1);
  }
  return total;
}

/* Internal to naive matrix mult */
static BIT _mzd_big_dot_product( packedmatrix *a, packedmatrix *bT, int rowofa,
				 int rowofb ) {
  /* ``a slot'' is a row of A, and a column of B when calcing AB */
  /* but since we use B^T so that we are working only with rows, */
  /* ``a slot'' of A is a row, ``a slot'' of B is a row of B^T */
  int total, i;

  total=0;
  for (i=0; i< a->width; i++) {
    //if (  (i*RADIX) < a->nrows )  
    total+=_mzd_dot_product( mzd_read_block(a, rowofa, i*RADIX), 
			     mzd_read_block(bT, rowofb, i*RADIX) );
  }

  return (BIT)(total % 2);
}

packedmatrix *mzd_mul_naiv_t(packedmatrix *C, packedmatrix *A, 
			     packedmatrix *bT ) {
  int i, j;
  int newrows=A->nrows;
  int newcols=bT->nrows;

  if (C==NULL) {
    C=mzd_init(newrows, newcols);
  } else {
    if (C->nrows != newrows || C->ncols != newcols) {
      m4ri_die("Provided return matrix has wrong dimensions.\n");
    }
  }

  for (i=0; i<newrows; i++) {
    for (j=0; j<newcols; j++) {
      mzd_write_bit(C, i, j, _mzd_big_dot_product( A, bT, i, j ) );
    }
  }

  return C;
}
  
packedmatrix *mzd_mul_naiv(packedmatrix *C, packedmatrix *A,  packedmatrix *B) {
  packedmatrix *bT = mzd_transpose(NULL, B);
  C = mzd_mul_naiv_t(C, A, bT );
  mzd_free(bT);

  return C;
}

void mzd_randomize( packedmatrix *a ) {
  int i, j;
  for (i=0; i < (a->nrows); i++) {
    for (j=0; j < (a->ncols); j++) {
      mzd_write_bit(a, i, j, m4ri_coin_flip() );
    }
  }
}

void mzd_set_ui( packedmatrix *a, unsigned int value) {
  if (value%2 == 0) {
    return;
  }

  int i,j;
  int stop = MIN(a->nrows, a->ncols);

  for (i=0; i< (a->nrows); i++) {
    for (j=0; j< (a->width); j++) {
      
      mzd_write_block(a, i, j*RADIX, 0);
    }
  }

  for (i=0; i<stop; i++) {
    mzd_write_bit(a, i, i, 1);
  }
}

BIT mzd_equal( packedmatrix *a, packedmatrix *b ) {
  int i, j;
  word block1, block2;

  if (a->nrows!=b->nrows) return FALSE;
  if (a->ncols!=b->ncols) return FALSE;

  for (i=0; i< a->nrows; i++) {
    for (j=0; j< a->width; j++) {
      block1=mzd_read_block(a, i, j*RADIX);
      block2=mzd_read_block(b, i, j*RADIX);
      if (block1 != block2)
	return FALSE;
    }
  }
  return TRUE;
}

int mzd_cmp(packedmatrix *a, packedmatrix *b) {

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

packedmatrix *mzd_copy(packedmatrix *n, packedmatrix *p) {
  if (n == NULL) {
    n = mzd_init(p->nrows, p->ncols);
  } else {
    if (n == p) {
      return p;
    } else if (n->nrows < p->nrows || n->ncols < p->ncols) {
      m4ri_die("Target matrix is too small.");
    }
  }
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

/* This is sometimes called augment */
packedmatrix *mzd_concat( packedmatrix *a, packedmatrix *b) {
  packedmatrix *newmatrix;
  int i, j, truerow, width;
  //word entry;
  
  if (a->nrows!=b->nrows) {
    m4ri_die("Bad arguments to concat!\n");
  }

  newmatrix=mzd_init(a->nrows, a->ncols + b->ncols);
  width = newmatrix->width;

  for (i=0; i<a->nrows; i++) {
    truerow = a->rowswap[i];
    for (j=0; j<a->width; j++) {
      newmatrix->values[i*width + j] = a->values[truerow + j];
    }
  }

  for (i=0; i<b->nrows; i++) {
    for (j=0; j<b->ncols; j++) {
      mzd_write_bit(newmatrix, i, j+(a->ncols), 
		mzd_read_bit(b, i, j) );
    }
  }

  return newmatrix;
}

packedmatrix *mzd_stack(packedmatrix *a, packedmatrix *b) {
  packedmatrix *newmatrix;
  int i, j, offset, truerow, width;

  if (a->ncols != b->ncols) {
    m4ri_die("Bad arguments to stack!\n");
  }
  newmatrix = mzd_init(a->nrows + b->nrows, a->ncols);
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

packedmatrix *mzd_invert_naiv(packedmatrix *target, packedmatrix *identity) {
  packedmatrix *huge, *inverse;
  int x;

  huge=mzd_concat(target, identity);

  x = mzd_reduce_naiv(huge, TRUE);

  if (x == FALSE) { mzd_free(huge); return NULL; }
  
  inverse=mzd_submatrix(huge, 0, target->ncols, target->nrows, 
			target->ncols*2);

  mzd_free(huge);
  return inverse;
}

packedmatrix *mzd_add(packedmatrix *ret, packedmatrix *left, packedmatrix *right) {
  if (left->nrows != right->nrows || left->ncols != right->ncols) {
    m4ri_die("rows and columns must match");
  }
  if (ret == NULL) {
    ret = mzd_copy(ret, left);
  } else if (ret != left) {
    if (ret->nrows != left->nrows || ret->ncols != left->ncols) {
      m4ri_die("rows and columns of returned matrix must match");
    }
  }
  return _mzd_add_impl(ret, left, right);
}

packedmatrix *_mzd_add_impl(packedmatrix *ret, packedmatrix *left, packedmatrix *right) {
  int i,j,left_truerow, ret_truerow;

  if (ret != left) {
    for (i=0; i < MIN(ret->nrows, left->nrows); i++) {
      left_truerow = left->rowswap[i];
      ret_truerow = ret->rowswap[i];
      for (j=0; j < left->width; j++) {
	ret->values[ret_truerow + j] = left->values[left_truerow + j];
      }
    }
  }
  
  for(i=0; i < MIN(ret->nrows, right->nrows); i++) {
    for(j=0; j < left->width; j++) {
      ret->values[  ret->rowswap[i] + j ]  ^= right->values[ right->rowswap[i] + j];
    }
  }
  return ret;
}

packedmatrix *mzd_submatrix(packedmatrix *m, int startrow, int startcol, int endrow, int endcol) {
  int nrows, ncols, truerow, i, colword, x, y, block, spot, startword;
  word temp  = 0;
  
  nrows = endrow - startrow;
  ncols = endcol - startcol;

  packedmatrix *newmatrix = mzd_init(nrows, ncols);

  startword = startcol / RADIX;

  /* we start at the beginning of a word */
  if (startcol%RADIX == 0) {
    for(x = startrow, i=0; i<nrows; i++, x+=1) {
      truerow = m->rowswap[x];

      /* process full words first */
      for(y = startcol, colword=0; colword<ncols/RADIX; colword++, y+=RADIX) {
	block = truerow + colword + startword;
	temp = m->values[block];
	newmatrix->values[newmatrix->rowswap[i] + colword] = temp;
      }

      /* process remaining bits */
      if (ncols%RADIX) {
	colword = ncols/RADIX;
	block = truerow + colword;
	temp = m->values[block] & ~((ONE<<(RADIX-ncols%RADIX))-1);
	newmatrix->values[newmatrix->rowswap[i] + colword] = temp;
      } 
    }

    /* startcol is not the beginning of a word */
  } else { 
    spot = startcol % RADIX;
    for(x = startrow, i=0; i<nrows; i++, x+=1) {
      truerow = m->rowswap[x];

      /* process full words first */
      for(y = startcol, colword=0; colword<ncols/RADIX; colword++, y+=RADIX) {
	block = truerow + colword + startword;
	temp = (m->values[block] << (spot)) | (m->values[block + 1] >> (RADIX-spot) ); 
	newmatrix->values[newmatrix->rowswap[i] + colword] = temp;
      }
      /* process remaining bits (lazy)*/
      colword = ncols/RADIX;
      for (y=0; y < ncols%RADIX; y++) {
	temp = mzd_read_bit(m, x, startcol + colword*RADIX + y);
	mzd_write_bit(newmatrix, i, colword*RADIX + y, temp);
      }
    }
  }
  return newmatrix;
}

