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

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

#include <assert.h>

#include "brilliantrussian.h"
#include "grayflex.h"
#include "misc.h"

/**
 * Finds a pivot row between xstart and xstop. The column where this
 * pivot is search is y. Returns TRUE if such a pivot row was
 * found. Also, the appropriate row is swapped to the top (== xstart).
 *
 * @param m matrix to operate on
 * @param xstart start row
 * @param xstop stop row (exclusive)
 * @param y column to read
 *
 * @return True if a pivot row was found
 */

static int _mzd_force_non_zero(packedmatrix *m, int xstart, int xstop, int y);

/**
 * Performs Gaussian elimination on a submatrix of 3k x k starting at
 * point (homepoint, homepoint) of m.
 *
 * \param m matrix to operate on
 * \param homepoint row,col where to start
 * \param k
 *
 * \return rank of 3k x k submatrix.
 */

static int _mzd_prep(packedmatrix *m, const int ai, const int k);

/**
 * Get k bits starting a position (x,y) from the matrix m.
 *
 * \param m Source matrix.
 * \param x Starting row.
 * \param y Starting column.
 * \param k Number of bits.
 */ 

static inline int _mzd_get_bits(const packedmatrix *m, const int x, const int y, const int k);


/**-----------------------------------------------------------------------**/

static int _mzd_force_non_zero(packedmatrix *m, int xstart, int xstop, int y) {
  int i;

  for (i=xstart; i<xstop; i++) {
    if (mzd_read_bit(m, i, y)==1) {
      if (i!=xstart) mzd_row_swap(m, i, xstart);
      return TRUE;
    }
  }
  return FALSE;
}

static int _mzd_prep(packedmatrix *m, const int ai, const int k) {
  int pc; /* pivot column */
  int tr; /* target row */
  int good;

  int rank = 0;

  for (pc=ai; pc < MIN(ai+k,m->ncols); pc++) {
    /* Step one, find a pivot row in this column.*/
    good=_mzd_force_non_zero(m, pc, MIN( ai+k*3, m->nrows-1 ), pc);

    if (good == FALSE) 
      return rank;

    for (tr=ai; tr < MIN(ai+k*3, m->nrows); tr++) {
      /* Step two, add this pivot row to other rows as needed. */
      if (tr==pc) continue;
      
      if (mzd_read_bit(m, tr, pc)==0) continue;

      mzd_row_add_offset(m, pc, tr, ai);
    }
    rank++;
  }

  return rank;
}

static inline int _mzd_get_bits(const packedmatrix *m, const int x, const int y, const int k) {
  int truerow = m->rowswap[ x ];
  word temp;

  /* there are two possible situations. Either all bits are in one
   * word or they are spread across two words. */

  if ( (y%RADIX + k -1 ) < RADIX ) {
    /* everything happens in one word here */
    temp =  m->values[ y / RADIX + truerow ]; /* get the value */
    temp <<= y%RADIX; /* clear upper bits */
    temp >>= RADIX - k; /* clear lower bits and move to correct position.*/
    return (int)temp;

  } else { 
    /* two words are affected */
    const int block = y / RADIX + truerow; /* correct block */
    const int spot = (y + k ) % RADIX; /* correct offset */
    /* make room by shifting spot times to the right, and add stuff from the second word */
    temp = (m->values[block] << spot) | ( m->values[block + 1] >> (RADIX - spot) ); 
    return ((int)temp & ((1<<k)-1)); /* clear upper bits and return */
   }
}


void mzd_make_table( packedmatrix *m, int ai, int k,
		     packedmatrix *T, int *L, int full) {
  const int homeblock= full ? 0 : ai/RADIX;
  int i, rowneeded, id;
  int twokay= TWOPOW(k);

  L[0]=0;

  for (i=1; i<twokay; i++) {
    rowneeded=codebook[k]->inc[i-1]+ai;

    id=codebook[k]->ord[i];

    L[id]=i;

    mzd_combine( T,  i,         homeblock,
		 m,  rowneeded, homeblock, 
		 T,  i-1,       homeblock);
  }
}



void mzd_process_row(packedmatrix *m, const int row, const int homecol, const int k, const packedmatrix *T, const int *L) {
  const int blocknum = homecol/RADIX;

  const int value = _mzd_get_bits(m, row, homecol, k);

  const int tablerow = L[ value ];
  
  mzd_combine(m, row,      blocknum,
	      m, row,      blocknum,
	      T, tablerow, blocknum);
}

void mzd_process_rows(packedmatrix *m, int startrow, int stoprow, int startcol, int k, packedmatrix *T, int *L) {
  int i,j;
  const int blocknum=startcol/RADIX;
  int wide = m->width - blocknum;

  int value, tablerow;
  word *b1_ptr,*b2_ptr;

  /* for optimization reasons we distinguish several cases here. */

  switch(wide) {

  case 1:
    /* no loop needed as only one block is operated on. */
    for (i=startrow; i<=stoprow; i++) {
      value = _mzd_get_bits(m, i, startcol, k);
      tablerow = L[ value ];
      b1_ptr = m->values + blocknum + m->rowswap[ i ];
      b2_ptr = T->values + blocknum + T->rowswap[ tablerow ];
      *b1_ptr ^= *b2_ptr;
    }
    break;

  case 2:
    /* two blocks, no loop */
    for (i=startrow; i<=stoprow; i++) {
      value = _mzd_get_bits(m, i, startcol, k);
      tablerow = L[ value ];
      b1_ptr = m->values + blocknum + m->rowswap[ i ];
      b2_ptr = T->values + blocknum + T->rowswap[ tablerow ];
      *b1_ptr++ ^= *b2_ptr++;
      *b1_ptr ^= *b2_ptr;
    }
    break;

  default:
    /* the real deal more than two blocks. */
    for (i=startrow; i<=stoprow; i++) {
      const int tablerow = L[ _mzd_get_bits(m, i, startcol, k) ];
      word *m_ptr = m->values + blocknum + m->rowswap[i];
      word *T_ptr = T->values + blocknum + T->rowswap[tablerow];
#ifdef HAVE_SSE2
      /** check alignments **/
      if (wide > SSE2_CUTOFF) {
	if (ALIGNMENT(T_ptr,16) == ALIGNMENT(m_ptr, 16)) {
	  do {
	    *m_ptr++ ^= *T_ptr++;
	    wide--;
	  } while(ALIGNMENT(m_ptr,16) && wide);
	}
      
	if (ALIGNMENT(m_ptr,16)==0 && ALIGNMENT(T_ptr,16)==0) {
	  __m128i *__m_ptr = (__m128i*)m_ptr;
	  __m128i *__T_ptr = (__m128i*)T_ptr;
	  const __m128i *end_ptr = (__m128i*)((unsigned long)(m_ptr + wide) & ~0xF);
	  __m128i xmm1;

	  do {
	   xmm1 = _mm_load_si128(__m_ptr    );
	   const __m128i xmm2 = _mm_load_si128(__T_ptr    );
	   xmm1 = _mm_xor_si128(xmm1, xmm2);
	   _mm_store_si128(__m_ptr    , xmm1);
	   __m_ptr++;
	   __T_ptr++;
	  } while(__m_ptr < end_ptr);

	  m_ptr = (word*)__m_ptr;
	  T_ptr = (word*)__T_ptr;
	  wide = ((sizeof(word)*wide)%16)/sizeof(word);
	}
      }
#endif
      for(j=0; j<wide ; j++)
      	m_ptr[j] ^= T_ptr[j];
#ifdef HAVE_SSE2      
      wide = m->width - blocknum;
#endif
    } /* end row loop */
  } /* end switch case */
}

int mzd_step_m4ri(packedmatrix *m, int full, int k, int ai, 
		  packedmatrix *T, int *L) {
  int submatrixrank;
  
  /**
   * The algorithm works as follows:
   *
   * Step 1.Denote the first column to be processed in a given
   * iteration as \f$a_i\f$. Then, perform Gaussian elimination on the
   * first \f$3k\f$ rows after and including the \f$i\f$-th row to
   * produce an identity matrix in \f$a_{i,i} ... a_{i+k-1,i+k-1},\f$
   * and zeroes in \f$a_{i+k,i} ... a_{i+3k-1,i+k-1}\f$.
   */
  submatrixrank = _mzd_prep(m, ai, k);
  if (submatrixrank!=k) 
    return submatrixrank;

  /**
   * Step 2. Construct a table consisting of the \f$2^k\f$ binary strings of
   * length k in a Gray code.  Thus with only \f$2^k\f$ vector
   * additions, all possible linear combinations of these k rows
   * have been precomputed.
   */

  mzd_make_table(m, ai, k, T, L, 0);

  /**
   * Step 3. One can rapidly process the remaining rows from \f$i +
   * 3k\f$ until row \f$m\f$ (the last row) by using the table. For
   * example, suppose the \f$j\f$-th row has entries \f$a_{j,i}
   * ... a_{j,i+k-1}\f$ in the columns being processed. Selecting the
   * row of the table associated with this k-bit string, and adding it
   * to row j will force the k columns to zero, and adjust the
   * remaining columns from \f$ i + k\f$ to n in the appropriate way,
   * as if Gaussian elimination had been performed.
  */

  mzd_process_rows(m, ai+k*3, m->nrows-1, ai, k, T, L);

  /**
   * Step 4. While the above form of the algorithm will reduce a
   * system of boolean linear equations to unit upper triangular form,
   * and thus permit a system to be solved with back substitution, the
   * M4RI algorithm can also be used to invert a matrix, or put the
   * system into reduced row echelon form (RREF). Simply run Step 3
   * on rows \f$0 ... i-1\f$ as well as on rows \f$i + 3k
   * ... m\f$. This only affects the complexity slightly, changing the
   * 2.5 coeffcient to 3
   */

  if (full==TRUE)
    mzd_process_rows(m, 0, ai-1, ai, k, T, L);

  return submatrixrank;
}

int mzd_reduce_m4ri(packedmatrix *m, int full, int k, packedmatrix *T, int *L) {
  int i, submatrixrank;
  int stop = MIN(m->nrows, m->ncols); 

  int rank = 0;
  int simple = 0;

  if (k == 0) {
    k = m4ri_opt_k(m->nrows, m->ncols, 0);
  }

  if (T == NULL && L == NULL) {
    simple = 1;
    T = mzd_init( TWOPOW(k), m->ncols );
    L = (int *)m4ri_mm_calloc( TWOPOW(k), sizeof(int) );
  }
  
  /* main loop */
  for (i=0; i<stop; i+=k) {
    /* not enough room for M4RI left. */
    if ( ((i+k*3) > m->nrows) || ((i+k) > m->ncols) ) {
      rank += mzd_gauss_delayed(m, i, full);
      break;
    }
    
    submatrixrank=mzd_step_m4ri(m, full, k, i, T, L);

    if (submatrixrank!=k) {
      /* not full rank, use Gaussian elimination :-( */
      rank += mzd_gauss_delayed(m, i, full);
      break;
    }
    rank += submatrixrank;
  }

  if (simple) {
    m4ri_mm_free(L);
    mzd_free(T);
  }
  return rank; 
}

void mzd_top_reduce_m4ri(packedmatrix *m, int k, packedmatrix *T, int *L) {
  int i, submatrixrank;
  int stop = MIN(m->nrows, m->ncols); 
  int simple = 0;
  
  if (k == 0) {
    k = m4ri_opt_k(m->nrows, m->ncols, 0);
  }
  
  /* setup tables */
  if (T == NULL && L == NULL) {
    simple = 1;
    T = mzd_init( TWOPOW(k), m->ncols );
    L = (int *)m4ri_mm_calloc( TWOPOW(k), sizeof(int) );
  }
  
  /* main loop */
  for (i=0; i<stop; i+=k) {
    if ( (i+k > m->nrows) || (i+k > m->ncols) ) {
      mzd_gauss_delayed(m, i, 1);
      break;
    }
    
    submatrixrank = _mzd_prep(m, i, k);
    
    if (submatrixrank==k) {
      mzd_make_table(m, i, k, T, L, 0);
      mzd_process_rows(m, 0, i-1, i, k, T, L);
    } else {
      mzd_gauss_delayed(m, i, 1);
      break;
    }
  }
  
  /* clear tables */
  if (simple) {
    m4ri_mm_free(L);
    mzd_free(T);
  }
}

packedmatrix *mzd_invert_m4ri(packedmatrix *m, packedmatrix *I, int k) {
  packedmatrix *big = mzd_concat(NULL, m, I);
  int size=m->ncols;
  if (k == 0) {
    k = m4ri_opt_k(m->nrows, m->ncols, 0);
  }
  int twokay=TWOPOW(k);
  int i;
  packedmatrix *mytable=mzd_init(twokay, size*2);
  int *mylookups=(int *)m4ri_mm_calloc(twokay, sizeof(int));
  packedmatrix *answer;
  
  mzd_reduce_m4ri(big, TRUE, k, mytable, mylookups);
  
  for(i=0; i < size; i++) {
    if (!mzd_read_bit(big, i,i )) {
      answer = NULL;
      break;
    }
  }
  if (i == size)
    answer=mzd_submatrix(NULL, big, 0, size, size, size*2);
  
  m4ri_mm_free(mylookups);
  mzd_free(mytable);
  mzd_free(big);
  
  return answer;
}

packedmatrix *mzd_mul_m4rm_t(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k) {
  packedmatrix *AT, *BT, *CT;
  
  if(A->ncols != B->nrows) 
    m4ri_die("mzd_mul_m4rm_t: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);
  
  AT = mzd_transpose(NULL, A);
  BT = mzd_transpose(NULL, B);
  
  CT = mzd_init(B->ncols, A->nrows);
  CT = _mzd_mul_m4rm_impl(CT, BT, AT, k, NULL, NULL, 0);
  
  mzd_free(AT);
  mzd_free(BT);

  C = mzd_transpose(C, CT);
  mzd_free(CT);
  return C;
}

packedmatrix *mzd_mul_m4rm(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k, packedmatrix *T, int *L) {
  int a = A->nrows;
  int b = A->ncols;
  int c = B->ncols;

  if(A->ncols != B->nrows) 
    m4ri_die("mzd_mul_m4rm_t: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);
  if (C == NULL) {
    C = mzd_init(a, c);
  } else {
    if (C->nrows != a || C->ncols != c)
      m4ri_die("mzd_mul_m4rm: C (%d x %d) has wrong dimensions.\n", C->nrows, C->ncols);
  }
  return _mzd_mul_m4rm_impl(C, A, B, k, T, L, TRUE);

}

packedmatrix *mzd_addmul_m4rm(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k, packedmatrix *T, int *L) {
  int a = A->nrows;
  int b = A->ncols;
  int c = B->ncols;

  if(A->ncols != B->nrows) 
    m4ri_die("mzd_mul_m4rm A ncols (%d) need to match B nrows (%d) .\n", A->ncols, B->nrows);
  if (C == NULL) {
    C = mzd_init(a, c);
  } else {
    if (C->nrows != a || C->ncols != c)
      m4ri_die("mzd_mul_m4rm: C has wrong dimensions.\n");
  }
  return _mzd_mul_m4rm_impl(C, A, B, k, T, L, FALSE);
}

packedmatrix *_mzd_mul_m4rm_impl(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k, packedmatrix *T, int *L, int clear) {
  int i,j, a,b,c, simple;
  int truerow;
  unsigned int x;

  a = A->nrows;
  b = A->ncols;
  c = B->ncols;

  int wide = C->width;

  /** clear first **/
  if (clear) {
    for (i=0; i<C->nrows; i++) {
      truerow = C->rowswap[i];
      for (j=0; j<C->width; j++) {
  	C->values[truerow + j] = 0;
     }
    }
  }

  if (k == 0) {
    k = m4ri_opt_k(a,b,c);
  }

  if (T == NULL && L == NULL) {
    simple = 1;
    T = mzd_init(TWOPOW(k), c);
    L = (int *)m4ri_mm_calloc(TWOPOW(k), sizeof(int));
  }

  /**
   * The algorithm proceeds as follows:
   * 
   * Step 1. Make a Gray code table of all the \f$2^k\f$ linear combinations
   * of the \f$k\f$ rows of \f$B_i\f$.  Call the \f$x\f$-th row
   * \f$T_x\f$.
   *
   * Step 2. Read the entries 
   *    \f$a_{j,(i-1)k+1}, a_{j,(i-1)k+2} , ... , a_{j,(i-1)k+k}.\f$
   *
   * Let \f$x\f$ be the \f$k\f$ bit binary number formed by the
   * concatenation of \f$a_{j,(i-1)k+1}, ... , a_{j,ik}\f$.
   *
   * Step 3. for \f$h = 1,2, ... , c\f$ do
   *   calculate \f$C_{jh} = C_{jh} + T_{xh}\f$.
   */


  if (wide < SSE2_CUTOFF) {
    for(i=0; i < b/k; i++) {
      mzd_make_table( B, i*k, k, T, L, 1 );
      for(j = 0; j<a; j++) {
        x = L[ _mzd_get_bits(A, j, i*k, k) ];
        /* mzd_combine( C,j,0, C,j,0,  T,x,0); */
        word *C_ptr = C->values + C->rowswap[j];
        const word *T_ptr = T->values + T->rowswap[x];
        for(int ii=0; ii<wide ; ii++)
          C_ptr[ii] ^= T_ptr[ii];
      }
    }
  } else {
    for(i=0; i < b/k; i++) {
      mzd_make_table( B, i*k, k, T, L, 1 );
      for(j=0; j<a; j++) {
        x = L[ _mzd_get_bits(A, j, i*k, k) ];
        unsigned int togo = wide;
        word *C_ptr = C->values + C->rowswap[j];
        const word *T_ptr = T->values + T->rowswap[x];
#ifdef HAVE_SSE2
        /** check alignments **/
        if (ALIGNMENT(T_ptr,16) == ALIGNMENT(C_ptr,16)) {
          do {
            *C_ptr++ ^= *T_ptr++;
            togo--;
          } while(ALIGNMENT(C_ptr,16) && togo);
        }

        if (ALIGNMENT(C_ptr,16)==0 && ALIGNMENT(T_ptr,16)==0) {
          __m128i *dst_ptr = (__m128i*)C_ptr;
          __m128i *src_ptr = (__m128i*)T_ptr;
          const __m128i *end_ptr = (__m128i*)((unsigned long)(C_ptr + togo) & ~0xF);
          __m128i xmm1;
          
          do {
            xmm1 = _mm_load_si128(dst_ptr);
            const __m128i xmm2 = _mm_load_si128(src_ptr);
            xmm1 = _mm_xor_si128(xmm1, xmm2);
            _mm_store_si128(dst_ptr, xmm1);
            ++src_ptr;
            ++dst_ptr;
          } while(dst_ptr < end_ptr);
          
          C_ptr = (word*)dst_ptr;
          T_ptr = (word*)src_ptr;
          togo = ((sizeof(word)*togo)%16)/sizeof(word);
        }
#endif //HAVE_SSE2
        for(int ii=0; ii<togo; ii++)
          C_ptr[ii] ^= T_ptr[ii];
      }
    }
  }

  /* handle rest */
  if (b%k) {
    mzd_make_table( B, b/k * k , b%k, T, L, 1);
    
    for(j = 0; j<a; j++) {
      x = L[ _mzd_get_bits(A, j, i*k, b%k) ];
      mzd_combine(C,j,0, C,j,0, T,x,0);
    }
  }
  if (simple) {
    mzd_free(T);
    m4ri_mm_free(L);
  }
  return C;
}

