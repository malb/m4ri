/******************************************************************************
*
*            M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org> 
*
*  Distributed under the terms of the GNU General Public License (GEL)
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
******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "packedmatrix.h"
#include "parity.h"

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

void mzd_print( const packedmatrix *M ) {
  size_t i, j;
  char temp[SAFECHAR];
  word *row;

  for (i=0; i< M->nrows; i++ ) {
    printf("[");
    row = M->values + M->rowswap[i];
    if(M->offset == 0) {
      for (j=0; j< M->width-1; j++) {
        m4ri_word_to_str(temp, row[j], 1);
        printf("%s ", temp);
      }
      row = row + M->width - 1;
      for (j=0; j< (size_t)(M->ncols%RADIX); j++) {
        if (GET_BIT(*row, j)) 
          printf("1");
        else
          printf(" ");
        if (((j % 4)==3) && (j!=RADIX-1))
          printf(":");
      }
    } else {
      for (j=0; j< M->ncols; j++) {
        if(mzd_read_bit(M, i, j))
          printf("1");
        else
          printf(" ");
        if (((j % 4)==3) && (j!=RADIX-1))
          printf(":");
      }
    }
    printf("]\n");
  }
}

void mzd_print_tight( const packedmatrix *M ) {
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

void mzd_row_add_offset(packedmatrix *M, size_t dstrow, size_t srcrow, size_t coloffset) {
  coloffset += M->offset;
  const size_t startblock= coloffset/RADIX;
  size_t wide = M->width - startblock;
  word *src = M->values + M->rowswap[srcrow] + startblock;
  word *dst = M->values + M->rowswap[dstrow] + startblock;

  word temp = *src++;
  if (coloffset%RADIX)
    temp = RIGHTMOST_BITS(temp, (RADIX-(coloffset%RADIX)-1));
  *dst++ ^= temp;
  wide--;

#ifdef HAVE_SSE2 
  if (ALIGNMENT(src,16)==8 && wide) {
    *dst++ ^= *src++;
    wide--;
  }
  __m128i *__src = (__m128i*)src;
  __m128i *__dst = (__m128i*)dst;
  const __m128i *eof = (__m128i*)((unsigned long)(src + wide) & ~0xF);
  __m128i xmm1;
  
  while(__src < eof) {
    xmm1 = _mm_xor_si128(*__dst, *__src++);
    *__dst++ = xmm1;
  }
  src  = (word*)__src;
  dst = (word*)__dst;
  wide = ((sizeof(word)*wide)%16)/sizeof(word);
#endif
  size_t i;
  for(i=0; i<wide; i++) {
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

int mzd_echelonize_naive(packedmatrix *m, int full) { 
  return mzd_gauss_delayed(m, 0, full); 
}

static inline packedmatrix *_mzd_transpose_direct(packedmatrix *DST, const packedmatrix *A) {
  size_t i,j,k, eol;
  word *temp;

  if(A->offset || DST->offset) {
    for(i=0; i<A->nrows; i++) {
      for(j=0; j<A->ncols; j++) {
        mzd_write_bit(DST, j, i, mzd_read_bit(A,i,j));
      }
    }
    return DST;
  }

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
  if (DST == NULL) {
    DST = mzd_init( A->ncols, A->nrows );
  } else {
    if (DST->nrows != A->ncols || DST->ncols != A->nrows) {
      m4ri_die("mzd_transpose: Wrong size for return matrix.\n");
    }
  }
  if(A->offset || DST->offset)
    return _mzd_transpose_direct(DST, A);
  else
    return _mzd_transpose(DST, A);
}

packedmatrix *mzd_mul_naive(packedmatrix *C, const packedmatrix *A, const packedmatrix *B) {
  packedmatrix *BT = mzd_transpose(NULL, B);

  if (C==NULL) {
    C=mzd_init(A->nrows, B->ncols);
  } else {
    if (C->nrows != A->nrows || C->ncols != B->ncols) {
      mzd_free (BT);
      m4ri_die("mzd_mul_naive: Provided return matrix has wrong dimensions.\n");
    }
  }
  _mzd_mul_naive(C, A, BT, 1);
  mzd_free (BT);
  return C;
}

packedmatrix *mzd_addmul_naive(packedmatrix *C, const packedmatrix *A, const packedmatrix *B) {
  packedmatrix *BT = mzd_transpose(NULL, B);

  if (C->nrows != A->nrows || C->ncols != B->ncols) {
    mzd_free (BT);
    m4ri_die("mzd_mul_naive: Provided return matrix has wrong dimensions.\n");
  }
  _mzd_mul_naive(C, A, BT, 0);
  mzd_free (BT);
  return C;
}

packedmatrix *_mzd_mul_naive(packedmatrix *C, const packedmatrix *A, const packedmatrix *B, const int clear) {
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

packedmatrix *_mzd_mul_va(packedmatrix *C, const packedmatrix *v, const packedmatrix *A, const int clear) {
  assert(C->offset == 0);
  assert(A->offset == 0);
  assert(v->offset == 0);

  if(clear)
    mzd_set_ui(C,0);

  size_t i,j;
  const size_t m=v->nrows;
  const size_t n=v->ncols;
  
  for(i=0; i<m; i++)
    
    for(j=0;j<n;j++)
      if (mzd_read_bit(v,i,j))
        mzd_combine(C,i,0,C,i,0,A,j,0);
  return C;
}

void mzd_randomize(packedmatrix *A) {
  size_t i, j;
  assert(A->offset == 0);

  for (i=0; i < A->nrows; i++) {
    for (j=0; j < A->ncols; j++) {
      mzd_write_bit(A, i, j, m4ri_coin_flip() );
    }
  }
}

void mzd_set_ui( packedmatrix *A, unsigned int value) {
  size_t i,j;
  size_t stop = MIN(A->nrows, A->ncols);

  word mask_begin = RIGHT_BITMASK(RADIX - A->offset);
  word mask_end = LEFT_BITMASK((A->offset + A->ncols)%RADIX);
  
  if(A->width==1) {
    for (i=0; i<A->nrows; i++) {
      for(j=0 ; j<A->ncols; j++)
        mzd_write_bit(A,i,j, 0);
    }
  } else {
    for (i=0; i<A->nrows; i++) {
      size_t truerow = A->rowswap[i];
      A->values[truerow] &= ~mask_begin;
      for(j=1 ; j<A->width-1; j++)
        A->values[truerow + j] = 0;
      A->values[truerow + A->width - 1] &= ~mask_end;
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

  if (A->nrows != B->nrows) return FALSE;
  if (A->ncols != B->ncols) return FALSE;

  for (i=0; i< A->nrows; i++) {
    for (j=0; j< A->width; j++) {
      if (A->values[A->rowswap[i] + j] != B->values[B->rowswap[i] + j])
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

packedmatrix *mzd_copy(packedmatrix *N, const packedmatrix *P) {
  if (N == P)
    return N;

  if (!P->offset){
    if (N == NULL) {
      N = mzd_init(P->nrows, P->ncols);
    } else {
      if (N->nrows < P->nrows || N->ncols < P->ncols)
	m4ri_die("mzd_copy: Target matrix is too small.");
    }
    size_t i, j, p_truerow, n_truerow;

    word mask = LEFT_BITMASK(P->ncols);
    for (i=0; i<P->nrows; i++) {
      p_truerow = P->rowswap[i];
      n_truerow = N->rowswap[i];
      for (j=0; j<P->width-1; j++) {
        N->values[n_truerow + j] = P->values[p_truerow + j];
      }
      N->values[n_truerow + j] = (N->values[n_truerow + j] & ~mask) | (P->values[p_truerow + j] & mask);
    }
  } else { // P->offset > 0
    if (N == NULL) {
      N = mzd_init(P->nrows, P->ncols+ P->offset);
      N->ncols -= P->offset;
      N->offset = P->offset;
      N->width=P->width;
    } else {
      if (N->nrows < P->nrows || N->ncols < P->ncols)
	m4ri_die("mzd_copy: Target matrix is too small.");
    }
    for(size_t i=0; i<P->nrows; i++) {
      mzd_copy_row(N, i, P, i);
    }
  }
/*     size_t i, j, p_truerow, n_truerow; */
/*     /\** */
/*      * \todo This is wrong  */
/*      *\/ */
/*     int trailingdim =  RADIX - P->ncols - P->offset; */

/*     if (trailingdim >= 0) { */
/*       // All columns fit in one word */
/*       word mask = ((ONE << P->ncols) - 1) << trailingdim; */
/*       for (i=0; i<P->nrows; i++) { */
/* 	p_truerow = P->rowswap[i]; */
/* 	n_truerow = N->rowswap[i]; */
/* 	N->values[n_truerow] = (N->values[n_truerow] & ~mask) | (P->values[p_truerow] & mask); */
/*       } */
/*     } else { */
/*       int r = (P->ncols + P->offset) % RADIX; */
/*       word mask_begin = RIGHT_BITMASK(RADIX - P->offset);  */
/*       word mask_end = LEFT_BITMASK(r); */
/*       for (i=0; i<P->nrows; i++) { */
/* 	p_truerow = P->rowswap[i]; */
/* 	n_truerow = N->rowswap[i]; */
/* 	N->values[n_truerow] = (N->values[n_truerow] & ~mask_begin) | (P->values[p_truerow] & mask_begin); */
/* 	for (j=1; j<P->width-1; j++) { */
/* 	  N->values[n_truerow + j] = P->values[p_truerow + j]; */
/* 	} */
/* 	N->values[n_truerow + j] = (N->values[n_truerow + j] & ~mask_end) | (P->values[p_truerow + j] & mask_end); */
/*       } */
/*     } */
  return N;
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

packedmatrix *mzd_invert_naive(packedmatrix *INV, packedmatrix *A, const packedmatrix *I) {
  assert(A->offset == 0);
  packedmatrix *H;
  int x;

  H = mzd_concat(NULL, A, I);

  x = mzd_echelonize_naive(H, TRUE);

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
    ret = mzd_init(left->nrows, left->ncols);
  } else if (ret != left) {
    if (ret->nrows != left->nrows || ret->ncols != left->ncols) {
      m4ri_die("mzd_add: rows and columns of returned matrix must match.\n");
    }
  }
  return _mzd_add(ret, left, right);
}

packedmatrix *_mzd_add(packedmatrix *C, const packedmatrix *A, const packedmatrix *B) {
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
  assert(M->offset == S->offset);

  startword = (M->offset + startcol) / RADIX;

  /* we start at the beginning of a word */
  if ((M->offset + startcol)%RADIX == 0) {
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
    spot = (M->offset + startcol) % RADIX;
    for(x = startrow, i=0; i<nrows; i++, x+=1) {
      truerow = M->rowswap[x];

      /* process full words first */
      for(colword=0; colword<(int)(ncols/RADIX); colword++) {
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

  size_t i;
  if(C->offset || A->offset || B->offset) {
    /**
     * \todo this code is slow if offset!=0 
     */
    for(i=0; i+RADIX<=A->ncols; i+=RADIX) {
      const word tmp = mzd_read_bits(A, a_row, i, RADIX) ^ mzd_read_bits(B, b_row, i, RADIX);
      for(size_t j=0; j<RADIX; j++) {
        mzd_write_bit(C, c_row, i*RADIX+j, GET_BIT(tmp, j));
      }
    }
    for( ; i<A->ncols; i++) {
      mzd_write_bit(C, c_row, i, mzd_read_bit(A, a_row, i) ^ mzd_read_bit(B, b_row, i));
    }
    return;
  }

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
  if (cola == colb)
    return;

  const size_t _cola = cola + M->offset;
  const size_t _colb = colb + M->offset;

  const size_t a_word = _cola/RADIX;
  const size_t b_word = _colb/RADIX;
  const size_t a_bit = _cola%RADIX;
  const size_t b_bit = _colb%RADIX;
  
  word a, b, *base;

  size_t i;
  
  if(a_word == b_word) {
    const word ai = RADIX - a_bit - 1;
    const word bi = RADIX - b_bit - 1;
    for (i=0; i<M->nrows; i++) {
      base = (M->values + M->rowswap[i] + a_word);
      register word b = *base;
      register word x = ((b >> ai) ^ (b >> bi)) & 1; // XOR temporary
      *base = b ^ ((x << ai) | (x << bi));
    }
    return;
  }

  const word a_bm = (ONE<<(RADIX - (a_bit) - 1));
  const word b_bm = (ONE<<(RADIX - (b_bit) - 1));

  if(a_bit > b_bit) {
    const size_t offset = a_bit - b_bit;
    for (i=0; i<M->nrows; i++) {
      base = M->values + M->rowswap[i];
      a = *(base + a_word);
      b = *(base + b_word);

      a ^= (b & b_bm) >> offset;
      b ^= (a & a_bm) << offset;
      a ^= (b & b_bm) >> offset;

      *(base + a_word) = a;
      *(base + b_word) = b;
    }
  } else {
    const size_t offset = b_bit - a_bit;
    for (i=0; i<M->nrows; i++) {
      base = M->values + M->rowswap[i];

      a = *(base + a_word);
      b = *(base + b_word);

      a ^= (b & b_bm) << offset;
      b ^= (a & a_bm) >> offset;
      a ^= (b & b_bm) << offset;
      *(base + a_word) = a;
      *(base + b_word) = b;
    }
  }

}

void mzd_col_swap_in_rows(packedmatrix *M, const size_t cola, const size_t colb, const size_t start_row, const size_t stop_row) {
  if (cola == colb)
    return;

  const size_t _cola = cola + M->offset;
  const size_t _colb = colb + M->offset;

  const size_t a_word = _cola/RADIX;
  const size_t b_word = _colb/RADIX;
  const size_t a_bit = _cola%RADIX;
  const size_t b_bit = _colb%RADIX;
  
  word a, b, *base;

  size_t i;
  
  if(a_word == b_word) {
    const word ai = RADIX - a_bit - 1;
    const word bi = RADIX - b_bit - 1;
    for (i=start_row; i<stop_row; i++) {
      base = (M->values + M->rowswap[i] + a_word);
      register word b = *base;
      register word x = ((b >> ai) ^ (b >> bi)) & 1; // XOR temporary
      *base = b ^ ((x << ai) | (x << bi));
    }
    return;
  }

  const word a_bm = (ONE<<(RADIX - (a_bit) - 1));
  const word b_bm = (ONE<<(RADIX - (b_bit) - 1));

  if(a_bit > b_bit) {
    const size_t offset = a_bit - b_bit;
    for (i=start_row; i<stop_row; i++) {
      base = M->values + M->rowswap[i];
      a = *(base + a_word);
      b = *(base + b_word);

      a ^= (b & b_bm) >> offset;
      b ^= (a & a_bm) << offset;
      a ^= (b & b_bm) >> offset;

      *(base + a_word) = a;
      *(base + b_word) = b;
    }
  } else {
    const size_t offset = b_bit - a_bit;
    for (i=start_row; i<stop_row; i++) {
      base = M->values + M->rowswap[i];
      a = *(base + a_word);
      b = *(base + b_word);

      a ^= (b & b_bm) << offset;
      b ^= (a & a_bm) >> offset;
      a ^= (b & b_bm) << offset;
      *(base + a_word) = a;
      *(base + b_word) = b;
    }
  }

}

int mzd_is_zero(packedmatrix *A) {
  /* Could be improved: stopping as the first non zero value is found (status!=0)*/
  size_t mb = A->nrows;
  size_t nb = A->ncols;
  size_t Aoffset = A->offset;
  size_t nbrest = (nb + Aoffset) % RADIX;
  int status=0;
  if (nb + Aoffset >= RADIX) {
          // Large A
    word mask_begin = RIGHT_BITMASK(RADIX-Aoffset);
    if (Aoffset == 0)
      mask_begin = ~mask_begin;
    word mask_end = LEFT_BITMASK(nbrest);
    size_t i;
    for (i=0; i<mb; ++i) {
        status |= A->values [A->rowswap [i]] & mask_begin;
        size_t j;
        for ( j = 1; j < A->width-1; ++j)
            status |= A->values [A->rowswap [i] + j];
        status |= A->values [A->rowswap [i] + A->width - 1] & mask_end;
    }
  } else {
          // Small A
    word mask = ((ONE << nb) - 1) ;
    mask <<= (RADIX-nb-Aoffset);

    size_t i;
    for (i=0; i < mb; ++i) {
        status |= A->values [A->rowswap [i]] & mask;
    }
  }
  
  return !status;
}

void mzd_copy_row(packedmatrix* B, size_t i, const packedmatrix* A, size_t j) {
  assert(B->offset == A->offset);
  assert(B->ncols >= A->ncols);
  size_t k;
  const size_t width= MIN(B->width, A->width) - 1;

  word* a=A->values + A->rowswap[j];
  word* b=B->values + B->rowswap[i];
 
  word mask_begin = RIGHT_BITMASK(RADIX - A->offset);
  word mask_end = LEFT_BITMASK( (A->offset + A->ncols)%RADIX );

  if (width != 0) {
    b[0] = (b[0] & ~mask_begin) | (a[0] & mask_begin);
    for(k = 1; k<width; k++)
      b[k] = a[k];
    b[width] = (b[width] & ~mask_end) | (a[width] & mask_end);
    
  } else {
    b[0] = (b[0] & ~mask_begin) | (a[0] & mask_begin & mask_end) | (b[0] & ~mask_end);
  }
}


void mzd_row_clear_offset(packedmatrix *M, size_t row, size_t coloffset) {
  coloffset += M->offset;
  size_t startblock= coloffset/RADIX;
  size_t i;
  word temp;

  /* make sure to start clearing at coloffset */
  if (coloffset%RADIX) {
    temp = M->values[M->rowswap[row] + startblock];
    temp &= RIGHT_BITMASK(RADIX - coloffset);
  } else {
    temp = 0;
  }
  M->values[M->rowswap[row] + startblock] = temp;
  temp=0;
  for ( i=startblock+1; i < M->width; i++ ) {
    M->values[M->rowswap[row] + i] = temp;
  }
}


int mzd_find_pivot(packedmatrix *A, size_t start_row, size_t start_col, size_t *r, size_t *c) { 
  assert(A->offset == 0);
  register size_t i = start_row;
  register size_t j = start_col;
  const size_t nrows = A->nrows;
  const size_t ncols = A->ncols;
  size_t row_candidate = 0;
  word data = 0;
  if(A->ncols - start_col < RADIX) {
    for(j=start_col; j<A->ncols; j+=RADIX) {
      const size_t length = MIN(RADIX, ncols-j);
      for(i=start_row; i<nrows; i++) {
        const word curr_data = (word)mzd_read_bits(A, i, j, length);
        if (curr_data > data && leftmost_bit(curr_data) > leftmost_bit(data)) {
          row_candidate = i;
          data = curr_data;
          if(GET_BIT(data,RADIX-length-1))
            break;
        }
      }
      if(data) {
        i = row_candidate;
        data <<=(RADIX-length);
        for(size_t l=0; l<length; l++) {
          if(GET_BIT(data, l)) {
            j+=l;
            break;
          }
        }
        *r = i, *c = j;
        return 1;
      }
    }
  } else {
    /* we definitely have more than one word */
    /* handle first word */
    const size_t bit_offset = (start_col % RADIX);
    const size_t word_offset = start_col / RADIX;
    const word mask_begin = RIGHT_BITMASK(RADIX-bit_offset);
    for(i=start_row; i<nrows; i++) {
      const word curr_data = A->values[A->rowswap[i] + word_offset] & mask_begin;
      if (curr_data > data && leftmost_bit(curr_data) > leftmost_bit(data)) {
        row_candidate = i;
        data = curr_data;
        if(GET_BIT(data,bit_offset)) {
          break;
        }
      }
    }
    if(data) {
      i = row_candidate;
      data <<=bit_offset;
      for(size_t l=0; l<(RADIX-bit_offset); l++) {
        if(GET_BIT(data, l)) {
          j+=l;
          break;
        }
      }
      *r = i, *c = j;
      return 1;
    }
    /* handle complete words */
    for(j=word_offset + 1; j<A->width - 1; j++) {
      for(i=start_row; i<nrows; i++) {
        const word curr_data = A->values[A->rowswap[i] + j];
        if (curr_data > data && leftmost_bit(curr_data) > leftmost_bit(data)) {
          row_candidate = i;
          data = curr_data;
          if(GET_BIT(data, 0))
            break;
        }
      }
      if(data) {
        i = row_candidate;
        for(size_t l=0; l<RADIX; l++) {
          if(GET_BIT(data, l)) {
            j=j*RADIX + l;
            break;
          }
        }
        *r = i, *c = j;
        return 1;
      }
    }
    /* handle last word */
    const size_t end_offset = A->ncols % RADIX ? (A->ncols%RADIX) : RADIX;
    const word mask_end = LEFT_BITMASK(end_offset);
    j = A->width-1;
    for(i=start_row; i<nrows; i++) {
      const word curr_data = A->values[A->rowswap[i] + j] & mask_end;
      if (curr_data > data && leftmost_bit(curr_data) > leftmost_bit(data)) {
        row_candidate = i;
        data = curr_data;
        if(GET_BIT(data,0))
          break;
      }
    }
    if(data) {
      i = row_candidate;
      for(size_t l=0; l<end_offset; l++) {
        if(GET_BIT(data, l)) {
          j=j*RADIX+l;
          break;
        }
      }
      *r = i, *c = j;
      return 1;
    }
  }
  return 0;
}


