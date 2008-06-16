/******************************************************************************
*
*            M4BI: Method of the Four Russians Inversion
*
*       Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org> 
*
*  Distributed under the terms of the GNU General Public License (GEL)
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
#include "parity.h"

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

#define SAFECHAR (int)(RADIX+RADIX/3)

packedmatrix *mzd_init(int r, int c) {
  packedmatrix *newmatrix;
  int i;

  newmatrix=(packedmatrix *)m4ri_mm_malloc(sizeof(packedmatrix));
  newmatrix->width=DIV_CEIL(c,RADIX);

#ifdef HAVE_SSE2
  int incw = 0;
  /* make sure each row is 16-byte aligned */
  if (newmatrix->width & 1) {
    newmatrix->width++;
    incw = 1;
  }
#endif

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

#ifdef HAVE_SSE2
  if (incw) {
    newmatrix->width--;
  }
#endif

  return newmatrix;
}

/* We don't perform any sanity checks! */
packedmatrix *mzd_init_window(const packedmatrix *m, int lowr, int lowc, int highr, int highc) {
  int nrows, ncols, i, offset; 
  packedmatrix *window = (packedmatrix *)m4ri_mm_malloc(sizeof(packedmatrix));
  nrows = MIN(highr - lowr, m->nrows - lowr);
  ncols = highc - lowc;
  
  window->ncols = ncols;
  window->nrows = nrows;
  window->width = ncols/RADIX;
  if (ncols%RADIX)
    window->width++;
  window->values = m->values;
  window->rowswap = (int *)m4ri_mm_malloc( nrows * sizeof(int));

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

void mzd_print_matrix( const packedmatrix *M ) {
  int i, j;
  char temp[SAFECHAR];
  word *row;

  for (i=0; i< M->nrows; i++ ) {
    printf("[ ");
    row = M->values + M->rowswap[i];
    for (j=0; j< M->ncols/RADIX; j++) {
      m4ri_word_to_str(temp, row[j], 1);
      printf("%s ", temp);
    }
    row = row + M->width - 1;
    for (j=0; j< (M->ncols%RADIX); j++) {
      printf("%d", (int)GET_BIT(*row, j));
      if (((j % 4)==3) && (j!=RADIX-1))
        printf(":");
    }
    if (M->ncols%RADIX)
      printf(" ]\n");
    else
      printf("]\n");
  }
}

void mzd_print_matrix_tight( const packedmatrix *M ) {
  int i, j;
  char temp[SAFECHAR];
  word *row;

  for (i=0; i< M->nrows; i++ ) {
    printf("[");
    row = M->values + M->rowswap[i];
    for (j=0; j< M->ncols/RADIX; j++) {
      m4ri_word_to_str(temp, row[j], 0);
      printf("%s", temp);
    }
    row = row + M->width - 1;
    for (j=0; j< (M->ncols%RADIX); j++) {
      printf("%d", (int)GET_BIT(*row, j));
    }
    printf("]\n");
  }
}

void mzd_row_clear_offset(packedmatrix *m, int row, int coloffset) {
  int startblock= coloffset/RADIX;
  int i;
  word temp;
  
  /* make sure to start clearing at coloffset */
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


void mzd_row_add_offset( packedmatrix *M, int srcrow, int dstrow, int coloffset) {
  int startblock= coloffset/RADIX;
  int i;
  word temp;
  
  /* make sure to start adding at coloffset */
  temp=mzd_read_block(M, srcrow, startblock*RADIX);
  if (coloffset%RADIX)
    temp &= (ONE<<(RADIX - (coloffset%RADIX))) - ONE;
  mzd_xor_block(M, dstrow, startblock*RADIX, temp);

  word *src = M->values + M->rowswap[srcrow];
  word *dst = M->values + M->rowswap[dstrow];
  for ( i=startblock+1; i < M->width; i++ ) {
    dst[i] ^= src[i];
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

	if (full==TRUE) 
          start=0; 
        else 
          start=startrow+1;

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

static inline packedmatrix *_mzd_transpose_direct(packedmatrix *DST, const packedmatrix *A) {
  int i,j,k, eol;
  word *temp;


  if(DST->ncols%RADIX) {
    eol = RADIX*(DST->width-1);
  } else {
    eol = RADIX*(DST->width);
  }

  for (i=0; i<DST->nrows; i++) {
    temp = DST->values + DST->rowswap[i];
    for (j=0; j < eol; j+=RADIX) {
      for (k=0; k<RADIX; k++) {
        *temp |= ((word)mzd_read_bit(A, j+k, i))<<(RADIX-1-k);
      }
      temp++;
    }
    j = A->nrows - (A->nrows%RADIX);
    for (k=0; k<(int)(A->nrows%RADIX); k++) {
      *temp |= ((word)mzd_read_bit(A, j+k, i))<<(RADIX-1-k);
    }
  }
  return DST;
}

static inline packedmatrix *_mzd_transpose_impl(packedmatrix *DST, const packedmatrix *X) {
  const int nr = X->nrows;
  const int nc = X->ncols;
  const int cutoff = 256; /* 256 seems optimal */

  if(nr <= cutoff || nc <= cutoff) {
    packedmatrix *x = mzd_copy(NULL, X);
    _mzd_transpose_direct(DST, x);
    mzd_free(x);
    return DST;
  }

  const int nr2 = RADIX*(X->nrows/(2*RADIX));
  const int nc2 = RADIX*(X->ncols/(2*RADIX));

  packedmatrix *A = mzd_init_window(X,    0,   0, nr2, nc2);
  packedmatrix *B = mzd_init_window(X,    0, nc2, nr2,  nc);
  packedmatrix *C = mzd_init_window(X,  nr2,   0,  nr, nc2);
  packedmatrix *D = mzd_init_window(X,  nr2, nc2,  nr,  nc);

  packedmatrix *AT = mzd_init_window(DST,   0,   0, nc2, nr2);
  packedmatrix *CT = mzd_init_window(DST,   0, nr2, nc2,  nr);
  packedmatrix *BT = mzd_init_window(DST, nc2,   0,  nc, nr2);
  packedmatrix *DT = mzd_init_window(DST, nc2, nr2,  nc,  nr);

  _mzd_transpose_impl(AT, A);
  _mzd_transpose_impl(BT, B);
  _mzd_transpose_impl(CT, C);
  _mzd_transpose_impl(DT, D);

  mzd_free_window(A); mzd_free_window(B);
  mzd_free_window(C); mzd_free_window(D);

  mzd_free_window(AT); mzd_free_window(CT);
  mzd_free_window(BT); mzd_free_window(DT);
  
  return DST;
}

packedmatrix *mzd_transpose(packedmatrix *DST, const packedmatrix *A) {
  if (DST == NULL) {
    DST = mzd_init( A->ncols, A->nrows );
  } else {
    if (DST->nrows != A->ncols || DST->ncols != A->nrows) {
      m4ri_die("mzd_transpose: Wrong size for return matrix.\n");
    }
  }
  return _mzd_transpose_impl(DST, A);
}

packedmatrix *mzd_mul_naiv(packedmatrix *C, const packedmatrix *A, const packedmatrix *B) {
  int i, j, k, ii, eol;
  packedmatrix *BT = mzd_transpose(NULL, B);
  word *a, *b, *c;

  if (C==NULL) {
    C=mzd_init(A->nrows, B->ncols);
  } else {
    if (C->nrows != A->nrows || C->ncols != B->ncols) {
      m4ri_die("mzd_mul_naiv_t: Provided return matrix has wrong dimensions.\n");
    }
  }

  if(C->ncols%RADIX) {
    eol = (C->width-1);
  } else {
    eol = (C->width);
  }

  word parity[64];
  for (i=0; i<64; i++) {
    parity[i] = 0;
  }
  const int wide = A->width;
  const int blocksize = MZD_MUL_BLOCKSIZE;
  int start;
  for (start = 0; start + blocksize <= C->nrows; start += blocksize) {
    for (i=start; i<start+blocksize; i++) {
      a = A->values + A->rowswap[i];
      c = C->values + C->rowswap[i];
      for (j=RADIX*(eol-1); j>=0; j-=RADIX) {
        for (k=RADIX-1; k>=0; k--) {
          b = BT->values + BT->rowswap[j+k];
          parity[k] = a[0] & b[0];
          for (ii=wide-1; ii>=1; ii--)
          parity[k] ^= a[ii] & b[ii];
        }
        c[j/RADIX] ^= parity64(parity);
      }
      
      if (eol != C->width) {
        for (k=0; k<(int)(C->ncols%RADIX); k++) {
          b = BT->values + BT->rowswap[RADIX*eol+k];
          parity[k] = a[0] & b[0];
          for (ii=1; ii<A->width; ii++)
            parity[k] ^= a[ii] & b[ii];
        }
        c[eol] ^= parity64(parity) & ~((ONE<<(RADIX-(C->ncols%RADIX)))-1);
      }
    }
  }

  for (i=C->nrows - (C->nrows%blocksize); i<C->nrows; i++) {
    a = A->values + A->rowswap[i];
    c = C->values + C->rowswap[i];
    for (j=RADIX*(eol-1); j>=0; j-=RADIX) {
      for (k=RADIX-1; k>=0; k--) {
        b = BT->values + BT->rowswap[j+k];
        parity[k] = a[0] & b[0];
        for (ii=wide-1; ii>=1; ii--)
          parity[k] ^= a[ii] & b[ii];
      }
      c[j/RADIX] ^= parity64(parity);
      }
    
    if (eol != C->width) {
      for (k=0; k<(int)(C->ncols%RADIX); k++) {
        b = BT->values + BT->rowswap[RADIX*eol+k];
        parity[k] = a[0] & b[0];
        for (ii=1; ii<A->width; ii++)
          parity[k] ^= a[ii] & b[ii];
      }
      c[eol] ^= parity64(parity) & ~((ONE<<(RADIX-(C->ncols%RADIX)))-1);
    }
  }

  mzd_free(BT);
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

  int i,j;
  int stop = MIN(a->nrows, a->ncols);

  for (i=0; i< (a->nrows); i++) {
    for (j=0; j< (a->width); j++) {
      
      mzd_write_block(a, i, j*RADIX, 0);
    }
  }

  if(value%2 == 0)
    return;

  for (i=0; i<stop; i++) {
    mzd_write_bit(a, i, i, 1);
  }
}

BIT mzd_equal(const packedmatrix *a, const packedmatrix *b ) {
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

int mzd_cmp(const packedmatrix *a, const packedmatrix *b) {

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

packedmatrix *mzd_copy(packedmatrix *n, const packedmatrix *p) {
  if (n == NULL) {
    n = mzd_init(p->nrows, p->ncols);
  } else {
    if (n == p) {
      return n;
    } else if (n->nrows < p->nrows || n->ncols < p->ncols) {
      m4ri_die("mzd_copy: Target matrix is too small.");
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
packedmatrix *mzd_concat(packedmatrix *C, const packedmatrix *A, const packedmatrix *B) {
  int i, j, src_truerow, dst_truerow;
  
  if (A->nrows != B->nrows) {
    m4ri_die("mzd_concat: Bad arguments to concat!\n");
  }

  if (C == NULL) {
    C = mzd_init(A->nrows, A->ncols + B->ncols);
  } else if (C->nrows != A->nrows || C->ncols != (A->ncols + B->ncols)) {
    m4ri_die("mzd_concat: C has wrong dimension!\n");
  }

  for (i=0; i<A->nrows; i++) {
    dst_truerow = C->rowswap[i];
    src_truerow = A->rowswap[i];
    for (j=0; j <A->width; j++) {
      C->values[dst_truerow + j] = A->values[src_truerow + j];
    }
  }

  for (i=0; i<B->nrows; i++) {
    for (j=0; j<B->ncols; j++) {
      mzd_write_bit(C, i, j+(A->ncols), mzd_read_bit(B, i, j) );
    }
  }

  return C;
}

packedmatrix *mzd_stack(packedmatrix *C, const packedmatrix *A, const packedmatrix *B) {
  int i, j, src_truerow, dst_truerow;

  if (A->ncols != B->ncols) {
    m4ri_die("mzd_stack: A->ncols (%d) != B->ncols (%d)!\n",A->ncols, B->ncols);
  }

  if (C == NULL) {
    C = mzd_init(A->nrows + B->nrows, A->ncols);
  } else if (C->nrows != (A->nrows + B->nrows) || C->ncols != A->ncols) {
    m4ri_die("mzd_stack: C has wrong dimension!\n");
  }
  
  for(i=0; i<A->nrows; i++) {
    src_truerow = A->rowswap[i];
    dst_truerow = C->rowswap[i];
    for (j=0; j<A->width; j++) {
      C->values[dst_truerow + j] = A->values[src_truerow + j]; 
    }
  }

  for(i=0; i<B->nrows; i++) {
    dst_truerow = C->rowswap[A->nrows + i];
    src_truerow = B->rowswap[i];
    for (j=0; j<B->width; j++) {
      C->values[dst_truerow + j] = B->values[src_truerow + j]; 
    }
  }
  return C;
}

packedmatrix *mzd_invert_naiv(packedmatrix *INV, packedmatrix *A, const packedmatrix *I) {
  packedmatrix *H;
  int x;

  H = mzd_concat(NULL, A, I);

  x = mzd_reduce_naiv(H, TRUE);

  if (x == FALSE) { 
    mzd_free(H); 
    return NULL; 
  }
  
  INV = mzd_submatrix(INV, H, 0, A->ncols, A->nrows, A->ncols*2);

  mzd_free(H);
  return INV;
}

packedmatrix *mzd_add(packedmatrix *ret, const packedmatrix *left, const packedmatrix *right) {
  if (left->nrows != right->nrows || left->ncols != right->ncols) {
    m4ri_die("mzd_add: rows and columns must match.\n");
  }
  if (ret == NULL) {
    ret = mzd_copy(ret, left);
  } else if (ret != left) {
    if (ret->nrows != left->nrows || ret->ncols != left->ncols) {
      m4ri_die("mzd_add: rows and columns of returned matrix must match.\n");
    }
  }
  return _mzd_add_impl(ret, left, right);
}

packedmatrix *_mzd_add_impl(packedmatrix *C, const packedmatrix *A, const packedmatrix *B) {
  int i;
  int nrows = MIN(MIN(A->nrows, B->nrows), C->nrows);
  const packedmatrix *tmp;

  if (C == B) { //swap
    tmp = A;
    A = B;
    B = tmp;
  }
  
  for(i=nrows-1; i>=0; i--) {
    mzd_combine(C,i,0, A,i,0, B,i,0);
  }
  return C;
}

packedmatrix *mzd_submatrix(packedmatrix *S, const packedmatrix *m, const int startrow, const int startcol, const int endrow, const int endcol) {
  int nrows, ncols, i, colword, x, y, block, spot, startword;
  unsigned int truerow;
  word temp  = 0;
  
  nrows = endrow - startrow;
  ncols = endcol - startcol;

  if (S == NULL) {
    S = mzd_init(nrows, ncols);
  } else if(S->nrows < nrows || S->ncols < ncols) {
    m4ri_die("mzd_submatrix: got S with dimension %d x %d but expected %d x %d\n",S->nrows,S->ncols,nrows,ncols);
  }

  startword = startcol / RADIX;

  /* we start at the beginning of a word */
  if (startcol%RADIX == 0) {
    if(ncols/RADIX) {
      for(x = startrow, i=0; i<nrows; i++, x++) {
        memcpy(S->values + S->rowswap[i], m->values + m->rowswap[x] + startword, 8*(ncols/RADIX));
      }
    }
    if (ncols%RADIX) {
      for(x = startrow, i=0; i<nrows; i++, x++) {
        /* process remaining bits */
	temp = m->values[m->rowswap[x] + startword + ncols/RADIX] & ~((ONE<<(RADIX-ncols%RADIX))-1);
	S->values[S->rowswap[i] + ncols/RADIX] = temp;
      } 
    }
    /* startcol is not the beginning of a word */
  } else { 
    spot = startcol % RADIX;
    for(x = startrow, i=0; i<nrows; i++, x+=1) {
      truerow = m->rowswap[x];

      /* process full words first */
      for(y = startcol, colword=0; colword<(int)(ncols/RADIX); colword++, y+=RADIX) {
	block = truerow + colword + startword;
	temp = (m->values[block] << (spot)) | (m->values[block + 1] >> (RADIX-spot) ); 
	S->values[S->rowswap[i] + colword] = temp;
      }
      /* process remaining bits (lazy)*/
      colword = ncols/RADIX;
      for (y=0; y < (int)(ncols%RADIX); y++) {
	temp = mzd_read_bit(m, x, startcol + colword*RADIX + y);
	mzd_write_bit(S, i, colword*RADIX + y, (BIT)temp);
      }
    }
  }
  return S;
}

void mzd_combine( packedmatrix * C, const int c_row, const int c_startblock,
		  const packedmatrix * A, const int a_row, const int a_startblock, 
		  const packedmatrix * B, const int b_row, const int b_startblock) {
  int i;
  int wide = A->width - a_startblock;

  word *a = A->values + a_startblock + A->rowswap[a_row];
  word *b = B->values + b_startblock + B->rowswap[b_row];
  
  if( C == A && a_row == c_row && a_startblock == c_startblock) {
#ifdef HAVE_SSE2
    if(wide > SSE2_CUTOFF) {
      /** check alignments **/
      if (ALIGNMENT(a,16)) {
        *a++ ^= *b++;
        wide--;
      }

      if (ALIGNMENT(a,16)==0 && ALIGNMENT(b,16)==0) {
	__m128i *a128 = (__m128i*)a;
	__m128i *b128 = (__m128i*)b;
	const __m128i *eof = (__m128i*)((unsigned long)(a + wide) & ~0xF);

	do {
	  *a128 = _mm_xor_si128(*a128, *b128);
	  ++b128;
	  ++a128;
	} while(a128 < eof);
	
	a = (word*)a128;
	b = (word*)b128;
	wide = ((sizeof(word)*wide)%16)/sizeof(word);
      }
    }
#endif //HAVE_SSE2
    for(i=wide-1; i >= 0; i--)
      a[i] ^= b[i];
    return;
    
  } else { /* C != A */
    word *c = C->values + c_startblock + C->rowswap[c_row];

    /* this is a corner case triggered by Strassen multiplication
       which assumes certain (virtual) matrix sizes */
    if (a_row >= A->nrows) {
      for(i = wide - 1 ; i >= 0 ; i--) {
        c[i] = b[i];
      }
    } else {
#ifdef HAVE_SSE2
    if(wide > SSE2_CUTOFF) {
      /** check alignments **/
      if (ALIGNMENT(a,16)) {
        *c++ = *b++ ^ *a++;
        wide--;
      }

      if ((ALIGNMENT(b,16)==0) && (ALIGNMENT(c,16)==0)) {
	__m128i *a128 = (__m128i*)a;
	__m128i *b128 = (__m128i*)b;
	__m128i *c128 = (__m128i*)c;
	const __m128i *eof = (__m128i*)((unsigned long)(a + wide) & ~0xF);
	
	do {
          *c128 = _mm_xor_si128(*a128, *b128);
	  ++c128;
	  ++b128;
	  ++a128;
	} while(a128 < eof);
	
	a = (word*)a128;
	b = (word*)b128;
	c = (word*)c128;
	wide = ((sizeof(word)*wide)%16)/sizeof(word);
      }
    }
#endif //HAVE_SSE2
    for(i = wide - 1 ; i >= 0 ; i--) {
      c[i] = a[i] ^ b[i];
    }
    return;
    }
  }
}
