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

packedmatrix *mzd_init(size_t r, size_t c) {
  packedmatrix *newmatrix;
  size_t i;

#ifdef HAVE_OPENMP
#pragma omp critical
{
#endif
  newmatrix=(packedmatrix *)m4ri_mmc_malloc(sizeof(packedmatrix));
#ifdef HAVE_OPENMP
 }
#endif

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
  newmatrix->offset = 0;
#ifdef HAVE_OPENMP
#pragma omp critical
{
#endif
  newmatrix->values=(word *)m4ri_mmc_calloc( (newmatrix->width)*r, sizeof(word) );
#ifdef HAVE_OPENMP
 }
#endif

#ifdef HAVE_OPENMP
#pragma omp critical
{
#endif
  newmatrix->rowswap=(size_t *)m4ri_mmc_malloc( r * sizeof(size_t) );
#ifdef HAVE_OPENMP
 }
#endif

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

packedmatrix *mzd_init_window (const packedmatrix *m, size_t lowr, size_t lowc, size_t highr, size_t highc) {
  size_t nrows, ncols, i, offset; 
  packedmatrix *window;
#ifdef HAVE_OPENMP
#pragma omp critical
{
#endif
  window = (packedmatrix *)m4ri_mmc_malloc(sizeof(packedmatrix));
#ifdef HAVE_OPENMP
}
#endif
  nrows = MIN(highr - lowr, m->nrows - lowr);
  ncols = highc - lowc;
  
  window->ncols = ncols;
  window->nrows = nrows;

  window->offset = (m->offset + lowc) % RADIX;
  offset = (m->offset + lowc) / RADIX;
  
  window->width = (window->offset + ncols) / RADIX;
  if ((window->offset + ncols) % RADIX)
    window->width++;
  window->values = m->values;

#ifdef HAVE_OPENMP
#pragma omp critical
{
#endif
  window->rowswap = (size_t *)m4ri_mmc_malloc( nrows * sizeof(size_t));
#ifdef HAVE_OPENMP
}
#endif
  for(i=0; i<nrows; i++) {
    window->rowswap[i] = m->rowswap[lowr + i] + offset;
  }
  
  return window;
}

permutation *mzd_init_permutation_window (permutation* P, size_t begin, size_t end){
  permutation *window = (permutation *)m4ri_mm_malloc(sizeof(permutation));
  window->values = P->values + begin;
  window->length = begin-end;
  return window;
}

void mzd_free( packedmatrix *condemned) {
#ifdef HAVE_OPENMP
#pragma omp critical
{
#endif
  m4ri_mmc_free(condemned->values, condemned->width*condemned->nrows*sizeof(word));
  m4ri_mmc_free(condemned->rowswap, condemned->nrows * sizeof(size_t));
  m4ri_mmc_free(condemned, sizeof(packedmatrix));
#ifdef HAVE_OPENMP
}
#endif
}

void mzd_free_window( packedmatrix *condemned) {
#ifdef HAVE_OPENMP
#pragma omp critical
{
#endif
  m4ri_mmc_free(condemned->rowswap, condemned->nrows * sizeof(size_t));
  m4ri_mmc_free(condemned, sizeof(packedmatrix));
#ifdef HAVE_OPENMP
}
#endif
}

void mzd_free_permutation_window (permutation* condemned){
  m4ri_mm_free(condemned);
}

void mzd_print_matrix( const packedmatrix *M ) {
  size_t i, j;
  char temp[SAFECHAR];
  word *row;

  for (i=0; i< M->nrows; i++ ) {
    printf("[ ");
    row = M->values + M->rowswap[i];
    /* TODO: This is not correct */
    for (j=0; j< (M->ncols+M->offset)/RADIX; j++) {
      m4ri_word_to_str(temp, row[j], 1);
      printf("%s ", temp);
    }
    row = row + M->width - 1;
    for (j=0; j< (size_t)((M->ncols+M->offset)%RADIX); j++) {
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
  assert(M->offset == 0);

  size_t i, j;
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
    for (j=0; j< (int)(M->ncols%RADIX); j++) {
      printf("%d", (int)GET_BIT(*row, j));
    }
    printf("]\n");
  }
}

void mzd_row_clear_offset(packedmatrix *M, size_t row, size_t coloffset) {
  assert(M->offset == 0);

  size_t startblock= coloffset/RADIX;
  size_t i;
  word temp;
  
  /* make sure to start clearing at coloffset */
  if (coloffset%RADIX) {
    temp=mzd_read_block(M, row, coloffset);
    temp &= RIGHT_BITMASK(RADIX-coloffset);
  } else {
    temp = 0;
  }
  mzd_write_block(M, row, coloffset, temp);

  temp=0;

  for ( i=startblock+1; i < (M->width); i++ ) {
    mzd_write_block(M, row, i*RADIX, temp);
  }
}


void mzd_row_add_offset( packedmatrix *M, size_t dstrow, size_t srcrow, size_t coloffset) {
  assert(M->offset == 0);

  size_t startblock= coloffset/RADIX;
  size_t i;
  
  /* make sure to start adding at coloffset */
  word *src = M->values + M->rowswap[srcrow];
  word *dst = M->values + M->rowswap[dstrow];
  word temp = src[startblock];

  if (coloffset%RADIX)
    temp = RIGHTMOST_BITS(temp, (RADIX-(coloffset%RADIX)-1));

  dst[startblock] ^= temp;

  for ( i=startblock+1; i < M->width; i++ ) {
    dst[i] ^= src[i];
  }
}


void mzd_row_add( packedmatrix *m, size_t sourcerow, size_t destrow) {
  mzd_row_add_offset(m, destrow, sourcerow, 0);
}

int mzd_gauss_delayed(packedmatrix *M, size_t startcol, int full) {
  assert(M->offset == 0);
  size_t i,j;
  size_t start; 

  size_t startrow = startcol;
  size_t ii;
  size_t pivots = 0;
  for (i=startcol ; i<M->ncols ; i++) {

    for(j=startrow ; j < M->nrows; j++) {
      if (mzd_read_bit(M,j,i)) {
	mzd_row_swap(M,startrow,j);
	pivots++;

	if (full==TRUE) 
          start=0; 
        else 
          start=startrow+1;

	for(ii=start ;  ii < M->nrows ; ii++) {
	  if (ii != startrow) {
	    if (mzd_read_bit(M, ii, i)) {
	      mzd_row_add_offset(M, ii, startrow, i);
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
  return mzd_gauss_delayed(m, 0, full); 
}

static inline packedmatrix *_mzd_transpose_direct(packedmatrix *DST, const packedmatrix *A) {
  assert(A->offset == 0);

  size_t i,j,k, eol;
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
        *temp |= ((word)mzd_read_bit(A, j+k, i+A->offset))<<(RADIX-1-k);
      }
      temp++;
    }
    j = A->nrows - (A->nrows%RADIX);
    for (k=0; k<(size_t)(A->nrows%RADIX); k++) {
      *temp |= ((word)mzd_read_bit(A, j+k, i+A->offset))<<(RADIX-1-k);
    }
  }
  DST->offset = 0;
  return DST;
}

static inline packedmatrix *_mzd_transpose(packedmatrix *DST, const packedmatrix *X) {
  assert(X->offset == 0);

  const size_t nr = X->nrows;
  const size_t nc = X->ncols;
  const size_t cutoff = 256; /* 256 seems optimal */

  if(nr <= cutoff || nc <= cutoff) {
    packedmatrix *x = mzd_copy(NULL, X);
    _mzd_transpose_direct(DST, x);
    mzd_free(x);
    return DST;
  }

  const size_t nr2 = RADIX*(X->nrows/(2*RADIX));
  const size_t nc2 = RADIX*(X->ncols/(2*RADIX));

  packedmatrix *A = mzd_init_window(X,    0,   0, nr2, nc2);
  packedmatrix *B = mzd_init_window(X,    0, nc2, nr2,  nc);
  packedmatrix *C = mzd_init_window(X,  nr2,   0,  nr, nc2);
  packedmatrix *D = mzd_init_window(X,  nr2, nc2,  nr,  nc);

  packedmatrix *AT = mzd_init_window(DST,   0,   0, nc2, nr2);
  packedmatrix *CT = mzd_init_window(DST,   0, nr2, nc2,  nr);
  packedmatrix *BT = mzd_init_window(DST, nc2,   0,  nc, nr2);
  packedmatrix *DT = mzd_init_window(DST, nc2, nr2,  nc,  nr);

  _mzd_transpose(AT, A);
  _mzd_transpose(BT, B);
  _mzd_transpose(CT, C);
  _mzd_transpose(DT, D);

  mzd_free_window(A); mzd_free_window(B);
  mzd_free_window(C); mzd_free_window(D);

  mzd_free_window(AT); mzd_free_window(CT);
  mzd_free_window(BT); mzd_free_window(DT);
  
  return DST;
}

packedmatrix *mzd_transpose(packedmatrix *DST, const packedmatrix *A) {
  assert(A->offset == 0);

  if (DST == NULL) {
    DST = mzd_init( A->ncols, A->nrows );
  } else {
    if (DST->nrows != A->ncols || DST->ncols != A->nrows) {
      m4ri_die("mzd_transpose: Wrong size for return matrix.\n");
    }
  }
  return _mzd_transpose(DST, A);
}

packedmatrix *mzd_mul_naiv(packedmatrix *C, const packedmatrix *A, const packedmatrix *B) {
  packedmatrix *BT = mzd_transpose(NULL, B);

  if (C==NULL) {
    C=mzd_init(A->nrows, B->ncols);
  } else {
    if (C->nrows != A->nrows || C->ncols != B->ncols) {
      mzd_free (BT);
      m4ri_die("mzd_mul_naiv: Provided return matrix has wrong dimensions.\n");
    }
  }
  _mzd_mul_naiv(C, A, BT, 1);
  mzd_free (BT);
  return C;
}

packedmatrix *mzd_addmul_naiv(packedmatrix *C, const packedmatrix *A, const packedmatrix *B) {
  packedmatrix *BT = mzd_transpose(NULL, B);

  if (C->nrows != A->nrows || C->ncols != B->ncols) {
    mzd_free (BT);
    m4ri_die("mzd_mul_naiv: Provided return matrix has wrong dimensions.\n");
  }
  _mzd_mul_naiv(C, A, BT, 0);
  mzd_free (BT);
  return C;
}

packedmatrix *_mzd_mul_naiv(packedmatrix *C, const packedmatrix *A, const packedmatrix *B, const int clear) {
  assert(A->offset == 0);
  assert(B->offset == 0);
  assert(C->offset == 0);
  size_t i, j, k, ii, eol;
  word *a, *b, *c;

  if (clear) {
    for (i=0; i<C->nrows; i++) {
      size_t truerow = C->rowswap[i];
      for (j=0; j<C->width-1; j++) {
  	C->values[truerow + j] = 0;
      }
      C->values[truerow + j] &= ~LEFT_BITMASK(C->ncols);
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
  const size_t wide = A->width;
  const size_t blocksize = MZD_MUL_BLOCKSIZE;
  size_t start;
  for (start = 0; start + blocksize <= C->nrows; start += blocksize) {
    for (i=start; i<start+blocksize; i++) {
      a = A->values + A->rowswap[i];
      c = C->values + C->rowswap[i];
      for (j=0; j<RADIX*eol; j+=RADIX) {
	for (k=0; k<RADIX; k++) {
          b = B->values + B->rowswap[j+k];
          parity[k] = a[0] & b[0];
          for (ii=wide-1; ii>=1; ii--)
	    parity[k] ^= a[ii] & b[ii];
        }
        c[j/RADIX] ^= parity64(parity);
      }
      
      if (eol != C->width) {
        for (k=0; k<(int)(C->ncols%RADIX); k++) {
          b = B->values + B->rowswap[RADIX*eol+k];
          parity[k] = a[0] & b[0];
          for (ii=1; ii<A->width; ii++)
            parity[k] ^= a[ii] & b[ii];
        }
        c[eol] ^= parity64(parity) & LEFT_BITMASK(C->ncols);
      }
    }
  }

  for (i=C->nrows - (C->nrows%blocksize); i<C->nrows; i++) {
    a = A->values + A->rowswap[i];
    c = C->values + C->rowswap[i];
    for (j=0; j<RADIX*eol; j+=RADIX) {
      for (k=0; k<RADIX; k++) {
        b = B->values + B->rowswap[j+k];
        parity[k] = a[0] & b[0];
        for (ii=wide-1; ii>=1; ii--)
          parity[k] ^= a[ii] & b[ii];
      }
      c[j/RADIX] ^= parity64(parity);
      }
    
    if (eol != C->width) {
      for (k=0; k<(int)(C->ncols%RADIX); k++) {
        b = B->values + B->rowswap[RADIX*eol+k];
        parity[k] = a[0] & b[0];
        for (ii=1; ii<A->width; ii++)
          parity[k] ^= a[ii] & b[ii];
      }
      c[eol] ^= parity64(parity) & LEFT_BITMASK(C->ncols);
    }
  }

  return C;
}

void mzd_randomize(packedmatrix *A) {
  assert(A->offset == 0);
  size_t i, j;
  for (i=0; i < A->nrows; i++) {
    for (j=0; j < A->ncols; j++) {
      mzd_write_bit(A, i, j, m4ri_coin_flip() );
    }
  }
}

void mzd_set_ui( packedmatrix *A, unsigned int value) {
  assert(A->offset == 0);

  size_t i,j;
  size_t stop = MIN(A->nrows, A->ncols);

  for (i=0; i< (A->nrows); i++) {
    for (j=0; j< (A->width); j++) {
      
      mzd_write_block(A, i, j*RADIX, 0);
    }
  }

  if(value%2 == 0)
    return;

  for (i=0; i<stop; i++) {
    mzd_write_bit(A, i, i, 1);
  }
}

BIT mzd_equal(const packedmatrix *A, const packedmatrix *B) {
  assert(A->offset == 0);
  assert(B->offset == 0);

  size_t i, j;
  word block1, block2;

  if (A->nrows != B->nrows) return FALSE;
  if (A->ncols != B->ncols) return FALSE;

  for (i=0; i< A->nrows; i++) {
    for (j=0; j< A->width; j++) {
      block1=mzd_read_block(A, i, j*RADIX);
      block2=mzd_read_block(B, i, j*RADIX);
      if (block1 != block2)
	return FALSE;
    }
  }
  return TRUE;
}

int mzd_cmp(const packedmatrix *A, const packedmatrix *B) {
  assert(A->offset == 0);
  assert(B->offset == 0);

  size_t i,j;

  if(A->nrows < B->nrows) return -1;
  if(B->nrows < A->nrows) return 1;
  if(A->ncols < B->ncols) return -1;
  if(B->ncols < A->ncols) return 1;

  for(i=0; i < A->nrows ; i++) {
    for(j=0 ; j< A->width ; j++) {
      if ( A->values[A->rowswap[i] + j] < B->values[B->rowswap[i] + j])
	return -1;
      else if( A->values[A->rowswap[i] + j] > B->values[B->rowswap[i] + j])
	return 1;
    }
  }
  return 0;
}

packedmatrix *mzd_copy(packedmatrix *n, const packedmatrix *p) {

  if (!p->offset){
    if (n == NULL) {
      n = mzd_init(p->nrows, p->ncols);
    } else {
      if (n == p) {
	return n;
      } else if (n->nrows < p->nrows || n->ncols < p->ncols) {
	m4ri_die("mzd_copy: Target matrix is too small.");
      }
    }
    size_t i, j, p_truerow, n_truerow;

    word mask = LEFT_BITMASK(p->ncols);
    for (i=0; i<p->nrows; i++) {
      p_truerow = p->rowswap[i];
      n_truerow = n->rowswap[i];
      for (j=0; j<p->width-1; j++) {
        n->values[n_truerow + j] = p->values[p_truerow + j];
      }
      n->values[n_truerow + j] = (n->values[n_truerow + j] & ~mask) | (p->values[p_truerow + j] & mask);
    }
  } else { // p->offset > 0
    if (n == NULL) {
      n = mzd_init(p->nrows, p->ncols+ p->offset);
      n->ncols -= p->offset;
    } else {
      if (n == p) {
	return n;
      } else if (n->nrows < p->nrows || n->ncols < p->ncols) {
	m4ri_die("mzd_copy: Target matrix is too small.");
      }
    }
    size_t i, j, p_truerow, n_truerow;
    /* TODO: This is wrong */
    int trailingdim =  RADIX - p->ncols - p->offset;

    if (trailingdim >= 0) {
      // All columns fit in one word
      word mask = ((ONE << p->ncols) - 1) << trailingdim;
      for (i=0; i<p->nrows; i++) {
	p_truerow = p->rowswap[i];
	n_truerow = n->rowswap[i];
	n->values[n_truerow] = (n->values[n_truerow] & ~mask) | (p->values[p_truerow] & mask);
      }
    } else {
      int r = (p->ncols + p->offset) % RADIX;
      word mask_begin = RIGHT_BITMASK(RADIX - p->offset); 
      word mask_end = LEFT_BITMASK(r);
      for (i=0; i<p->nrows; i++) {
	p_truerow = p->rowswap[i];
	n_truerow = n->rowswap[i];
	n->values[n_truerow] = (n->values[n_truerow] & ~mask_begin) | (p->values[p_truerow] & mask_begin);
	for (j=1; j<p->width-1; j++) {
	  n->values[n_truerow + j] = p->values[p_truerow + j];
	}
	n->values[n_truerow + j] = (n->values[n_truerow + j] & ~mask_end) | (p->values[p_truerow + j] & mask_end);
      }
    }
  }
  n->offset = p->offset;
  n->width=p->width;
  
  return n;
}

/* This is sometimes called augment */
packedmatrix *mzd_concat(packedmatrix *C, const packedmatrix *A, const packedmatrix *B) {
  assert(A->offset == 0);
  assert(B->offset == 0);
  size_t i, j, src_truerow, dst_truerow;
  
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
  assert(A->offset == 0);
  assert(B->offset == 0);
  size_t i, j, src_truerow, dst_truerow;

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
  assert(A->offset == 0);
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
  return _mzd_add(ret, left, right);
}

packedmatrix *_mzd_add(packedmatrix *C, const packedmatrix *A, const packedmatrix *B) {
  assert(C->offset == 0);
  assert(A->offset == 0);
  assert(B->offset == 0);
  size_t i;
  size_t nrows = MIN(MIN(A->nrows, B->nrows), C->nrows);
  const packedmatrix *tmp;

  if (C == B) { //swap
    tmp = A;
    A = B;
    B = tmp;
  }
  
  for(i=0; i<nrows; i++) {
    mzd_combine(C,i,0, A,i,0, B,i,0);
  }
  return C;
}

packedmatrix *mzd_submatrix(packedmatrix *S, const packedmatrix *M, const size_t startrow, const size_t startcol, const size_t endrow, const size_t endcol) {
  assert(M->offset == 0);
  size_t nrows, ncols, i, colword, x, y, block, spot, startword;
  size_t truerow;
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
        memcpy(S->values + S->rowswap[i], M->values + M->rowswap[x] + startword, 8*(ncols/RADIX));
      }
    }
    if (ncols%RADIX) {
      for(x = startrow, i=0; i<nrows; i++, x++) {
        /* process remaining bits */
	temp = M->values[M->rowswap[x] + startword + ncols/RADIX] & LEFT_BITMASK(ncols);
	S->values[S->rowswap[i] + ncols/RADIX] = temp;
      } 
    }
    /* startcol is not the beginning of a word */
  } else { 
    spot = startcol % RADIX;
    for(x = startrow, i=0; i<nrows; i++, x+=1) {
      truerow = M->rowswap[x];

      /* process full words first */
      for(y = startcol, colword=0; colword<(int)(ncols/RADIX); colword++, y+=RADIX) {
	block = truerow + colword + startword;
	temp = (M->values[block] << (spot)) | (M->values[block + 1] >> (RADIX-spot) ); 
	S->values[S->rowswap[i] + colword] = temp;
      }
      /* process remaining bits (lazy)*/
      colword = ncols/RADIX;
      for (y=0; y < (int)(ncols%RADIX); y++) {
	temp = mzd_read_bit(M, x, startcol + colword*RADIX + y);
	mzd_write_bit(S, i, colword*RADIX + y, (BIT)temp);
      }
    }
  }
  return S;
}

void mzd_combine( packedmatrix * C, const size_t c_row, const size_t c_startblock,
		  const packedmatrix * A, const size_t a_row, const size_t a_startblock, 
		  const packedmatrix * B, const size_t b_row, const size_t b_startblock) {
  assert(C->offset == 0);
  assert(A->offset == 0);
  assert(B->offset == 0);

  size_t i;
  size_t wide = A->width - a_startblock;

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
    for(i=0; i < wide; i++)
      a[i] ^= b[i];
    return;
    
  } else { /* C != A */
    word *c = C->values + c_startblock + C->rowswap[c_row];

    /* this is a corner case triggered by Strassen multiplication
       which assumes certain (virtual) matrix sizes */
    if (a_row >= A->nrows) {
      for(i = 0; i<wide; i++) {
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
    for(i = 0; i<wide; i++) {
      c[i] = a[i] ^ b[i];
    }
    return;
    }
  }
}


void mzd_col_swap(packedmatrix *M, const size_t cola, const size_t colb) {
  assert(M->offset == 0);
  if (cola == colb)
    return;

  const size_t dwa = cola/RADIX;
  const size_t dwb = colb/RADIX;
  const size_t dba = cola%RADIX;
  const size_t dbb = colb%RADIX;
  
  register word tmp;
  word *ptr_a, *ptr_b, *base;

  size_t i;
  
  for (i=0; i<M->nrows; i++) {
    base = M->values + M->rowswap[i];
    ptr_a = base + dwa;
    ptr_b = base + dwb;

    tmp = GET_BIT(*ptr_b, dbb);
    WRITE_BIT(*ptr_b, dbb, GET_BIT(*ptr_a, dba));
    WRITE_BIT(*ptr_a, dba, tmp);
  }
}

permutation *mzd_col_block_rotate(packedmatrix *M, size_t zs, size_t ze, size_t de, int zero_out, permutation *P) {
  assert(M->offset == 0);
  size_t i,j;

  const size_t ds = ze;
  const size_t ld_f = (de - ze)/RADIX;
  const size_t ld_r = (de - ds)%RADIX;

  const size_t lz_f = (ze - zs)/RADIX;
  const size_t lz_r = (ze - zs)%RADIX;

  word *tmp = (word*)m4ri_mm_calloc(DIV_CEIL(de-ze, RADIX), sizeof(word));

  for(i=0; i<M->nrows; i++) {
    /* copy out to tmp */
    for(j=0; j < ld_f; j++) {
      tmp[j] = mzd_read_bits(M, i, ds + j*RADIX, RADIX);
    }
    if (ld_r)
      tmp[ld_f] = mzd_read_bits(M, i, ds + ld_f*RADIX, ld_r);

    /* write to dst */
    for(j=0; j<ld_f; j++) {
      mzd_clear_bits(M, i, zs + j*RADIX, RADIX);
      mzd_write_zeroed_bits(M, i, zs + j*RADIX, RADIX, tmp[j]);
    }

    if(ld_r) {
      mzd_clear_bits(M, i, zs + ld_f*RADIX, ld_r);
      mzd_write_zeroed_bits(M, i, zs + ld_f*RADIX, ld_r, tmp[ld_f]);
    }
  }
  
  if (zero_out) {
    for(i=0; i<M->nrows; i++) {
      /* zero rest */
      for(j=0; j<lz_f; j++) {
        mzd_clear_bits(M, i, zs + (de - ds) + j*RADIX, RADIX);
      }
      if(lz_r)
        mzd_clear_bits(M, i, zs + (de - ds) + lz_f*RADIX, lz_r);
    }
  }

  if (P) {
    for(j=0; j<(de-ds); j++) {
      P->values[j] = P->values[de+j]; 
    }
  }
  m4ri_mm_free(tmp);
  return P;
}

void mzd_apply_p_left(packedmatrix *A, permutation *P) {
  assert(A->offset == 0);
  size_t i;
  for (i=0; i<P->length; i++) {
    if(P->values[i] != i) 
      mzd_row_swap(A, i, P->values[i]);
  }
}

void mzd_apply_p_left_trans(packedmatrix *A, permutation *P) {
  assert(A->offset == 0);
  size_t i;
  for (i=0; i<P->length; i++) {
    if(P->values[i] != i) 
      mzd_row_swap(A, i, P->values[i]);
  }
}

void mzd_apply_p_right_trans(packedmatrix *A, permutation *P) {
  assert(A->offset == 0);
  size_t i;
  for (i=0; i<P->length; i++) {
    if(P->values[i] != i) 
      mzd_col_swap(A, i, P->values[i]);
  }
}

void mzd_apply_p_right(packedmatrix *A, permutation *P) {
  assert(A->offset == 0);
  size_t i;
  for (i=0; 0<P->length; i++) {
    if(P->values[i] != i) 
      mzd_col_swap(A, i, P->values[i]);
  }
}
