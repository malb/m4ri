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

static int _mzd_prep(packedmatrix *m, int ai, int k);

/**
 * Get k bits starting a position (x,y) from the matrix m.
 *
 * \param m Source matrix.
 * \param x Starting row.
 * \param y Starting column.
 * \param k Number of bits.
 */ 

static int _mzd_get_bits(packedmatrix *m, int x, int y, int k);



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

static int _mzd_prep(packedmatrix *m, int ai, int k) {
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

static int _mzd_get_bits(packedmatrix *m, int x, int y, int k) {
  int truerow = m->rowswap[ x ];
  int block;
  int spot;

  word temp;

  word *values = m->values;

  /* there are two possible situations. Either all bits are in one
   * word or they are spread across two words. */

  if ( (y%RADIX + k -1 ) < RADIX ) {
    /* everything happens in one word here */
    temp =  values[ y / RADIX + truerow ]; // get the value
    temp <<= y%RADIX; // clear upper bits
    temp >>= RADIX - k; // clear lower bits and move to correct position.
    return (int)temp;

  } else { 
    /* two words are affected */
    block = y / RADIX + truerow; // correct block
    spot = (y + k ) % RADIX; // correct offset
    // make room by shifting spot times to the right, and add stuff from the second word
    temp = (values[block] << spot) | ( values[block + 1] >> (RADIX - spot) ); 
    return ((int)temp & ((1<<k)-1)); // clear upper bits and return
   }
}

void mzd_combine( packedmatrix * dst, int row3, int startblock3,
		  packedmatrix * sc1, int row1, int startblock1, 
		  packedmatrix * sc2, int row2, int startblock2) {
  int i;
  int wide = sc1->width - startblock1;

  word *b1_ptr = sc1->values + startblock1 + sc1->rowswap[row1];
  word *b2_ptr = sc2->values + startblock2 + sc2->rowswap[row2];
  word *b3_ptr;
  
  /* this is a quite likely case, and we treat is specially to ensure
   * cache register friendlyness. (Keep in mind that the x86 has only
   * four general purpose registers)
   */

  if( dst == sc1 && row1 == row3 && startblock1 == startblock3) {
#ifdef HAVE_SSE2
    if(wide>100) {
      /** check alignments **/
      unsigned long alignment = (unsigned long)b1_ptr%16;
      if ((unsigned long)b2_ptr%16 == alignment) {
	do {
	  *b1_ptr++ ^= *b2_ptr++;
	  wide--;
	} while((unsigned long)b1_ptr%16 && wide);
      }

      if (((unsigned long)b1_ptr%16==0) && ((unsigned long)b2_ptr%16==0)) {
	__m128i *dst_ptr = (__m128i*)b1_ptr;
	__m128i *src_ptr = (__m128i*)b2_ptr;
	__m128i *end_ptr = (__m128i*)((unsigned long)(b1_ptr + wide) & ~0xF);
	__m128i xmm1;

	do {
	  xmm1 = _mm_load_si128(dst_ptr);
	  const __m128i xmm2 = _mm_load_si128(src_ptr);
	  xmm1 = _mm_xor_si128(xmm1, xmm2);
	  _mm_store_si128(dst_ptr, xmm1);
	  ++src_ptr;
	  ++dst_ptr;
	} while(dst_ptr < end_ptr);
	
	b1_ptr = (word*)dst_ptr;
	b2_ptr = (word*)src_ptr;
	wide = ((sizeof(word)*wide)%16)/sizeof(word);
      }
    }
#endif
    for(i=wide-1; i >= 0; i--)
      b1_ptr[i] ^= b2_ptr[i];
    return;
    
  } else { // dst != sc1
    b3_ptr = dst->values + startblock3 + dst->rowswap[row3];

    if (row1 >= sc1->nrows) {
	for(i = wide - 1 ; i >= 0 ; i--) {
	  b3_ptr[i] = b2_ptr[i];
	}
    } else {
      for(i = wide - 1 ; i >= 0 ; i--) {
	b3_ptr[i] = b1_ptr[i] ^ b2_ptr[i];
      }
    }
    return;
  }
}

void mzd_make_table( packedmatrix *m, int ai, int k,
		     packedmatrix *T, int *L, int full) {
  int homeblock= full ? 0 : ai/RADIX;
  int i, rowneeded, id;
  int twokay= TWOPOW(k);

  L[0]=0;

  for (i=1; i<twokay; i++) {
    rowneeded=codebook[k]->inc[i-1]+ai;

    id=codebook[k]->ord[i];

    L[id]=i;

    mzd_combine( T,  i,         homeblock,
		 m,            rowneeded, homeblock, 
		 T,  i-1,       homeblock);
  }
}



void mzd_process_row(packedmatrix *m, int row, int homecol, int k, packedmatrix *T, int *L) {
  int blocknum = homecol/RADIX;

  int value = _mzd_get_bits(m, row, homecol, k);

  int tablerow = L[ value ];

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

  // for optimization reasons we distinguish several cases here.

  switch(wide) {

  case 1:
    // no loop needed as only one block is operated on.
    for (i=startrow; i<=stoprow; i++) {
      value = _mzd_get_bits(m, i, startcol, k);
      tablerow = L[ value ];
      b1_ptr = m->values + blocknum + m->rowswap[ i ];
      b2_ptr = T->values + blocknum + T->rowswap[ tablerow ];
      *b1_ptr ^= *b2_ptr;
    }
    break;

  case 2:
    // two blocks, no loop
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
    // the real deal more than two blocks.
    for (i=startrow; i<=stoprow; i++) {
      const int tablerow = L[ _mzd_get_bits(m, i, startcol, k) ];
      word *m_ptr = m->values + blocknum + m->rowswap[i];
      word *T_ptr = T->values + blocknum + T->rowswap[tablerow];
#ifdef HAVE_SSE2
      /** check alignments **/
      wide = m->width - blocknum;
      if (wide>10) {
	unsigned long alignment = (unsigned long)m_ptr%16;
	if ((unsigned long)T_ptr%16 == alignment) {
	  do {
	    *m_ptr++ ^= *T_ptr++;
	    wide--;
	  } while((unsigned long)m_ptr%16 && wide);
	}
      
	if (((unsigned long)m_ptr%16==0) && ((unsigned long)T_ptr%16==0)) {
	  __m128i *__m_ptr = (__m128i*)m_ptr;
	  __m128i *__T_ptr = (__m128i*)T_ptr;
	  __m128i *end_ptr = (__m128i*)((unsigned long)(m_ptr + wide) & ~0xF);
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
      for(j=wide-1; j>=0 ; j--)
      	m_ptr[j] ^= T_ptr[j];
    } // end row loop
  } // end switch case
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
   * remaining columns from \$f i + k\$f to n in the appropriate way,
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
  
  // main loop
  for (i=0; i<stop; i+=k) {
    // not enough room for M4RI left.
    if ( ((i+k*3) > m->nrows) || ((i+k) > m->ncols) ) {
      rank += mzd_gauss_delayed(m, i, full);
      break;
    }
    
    submatrixrank=mzd_step_m4ri(m, full, k, i, T, L);

    if (submatrixrank!=k) {
      // not full rank, use Gaussian elimination :-(
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

packedmatrix *mzd_invert_m4ri(packedmatrix *m, 
			 packedmatrix *identity, int k) {
  packedmatrix *big = mzd_concat(m, identity);
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
    answer=mzd_submatrix(big, 0, size, size, size*2);
  
  m4ri_mm_free(mylookups);
  mzd_free(mytable);
  mzd_free(big);
  
  return answer;
}

packedmatrix *mzd_mul_m4rm_t(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k) {
  packedmatrix *AT, *BT, *CT;
  
  //if(A->ncols != B->nrows) m4ri_die("A cols need to match B rows");
  
  AT = mzd_transpose(NULL, A);
  BT = mzd_transpose(NULL, B);
  
  CT = mzd_mul_m4rm(NULL, BT,AT,k, NULL, NULL);
  
  mzd_free(AT);
  mzd_free(BT);

  C = mzd_transpose(C, CT);
  mzd_free(CT);
  return C;
}

packedmatrix *mzd_mul_m4rm(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k, packedmatrix *T, int *L) {
  int i,j,ii,  a,b,c, simple;
  unsigned int x;

  word *C_ptr, *T_ptr;
  
  //if(A->ncols != B->nrows) 
  //  m4ri_die("A cols need to match B rows");
  
  a = A->nrows;
  b = A->ncols;
  c = B->ncols;


  if (C == NULL) {
    C = mzd_init(a, c);
  } else {
    //if (C->nrows != a || C->ncols != c) {
    //  m4ri_die("C has wrong dimensions.\n");
    //}
    /** clear first **/
    for (i=0; i<C->nrows; i++) {
      for (j=0; j<C->width; j++) {
	C->values[ C->rowswap[i] + j ] = 0;
      }
    }
  }
  int wide = C->width;

  if (k == 0) {
    k = m4ri_opt_k(a,b,c);
    //printf("k: %d\n",k);
  }

  if (T == NULL && L == NULL) {
    simple = 1;
    T = mzd_init(TWOPOW(k), c);
    L = (int *)m4ri_mm_calloc(TWOPOW(k), sizeof(int));
  }

  for(i=0; i < b/k; i++) {

    /**
     * The algorithm proceeds as follows:
     * 
     * Step 1. Make a Gray code table of all the \f$2^k\f$ linear combinations
     * of the \f$k\f$ rows of \f$B_i\f$.  Call the \f$x\f$-th row
     * \f$T_x\f$.
     */
    mzd_make_table( B, i*k, k, T, L, 1 );
    
    for(j = 0; j<a; j++) {

      /**
       * Step 2. Read the entries 
       *    \f$a_{j,(i-1)k+1}, a_{j,(i-1)k+2} , ... , a_{j,(i-1)k+k}.\f$
       *
       * Let \f$x\f$ be the \f$k\f$ bit binary number formed by the
       * concatenation of \f$a_{j,(i-1)k+1}, ... , a_{j,ik}\f$.
       */

      x = L[ _mzd_get_bits(A, j, i*k, k) ];

      /**
       * Step 3. for \f$h = 1,2, ... , c\f$ do
       *   calculate \f$C_{jh} = C_{jh} + T_{xh}\f$.
       */
      //mzd_combine( C,j,0, C,j,0,  T,x,0);
      C_ptr = C->values + C->rowswap[j];
      T_ptr = T->values + T->rowswap[x];
      for(ii=wide-1; ii>=0 ; ii--)
	C_ptr[ii] ^= T_ptr[ii];
    }
  }

  //handle rest
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


static packedmatrix *_mzd_mul_strassen_impl(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  int a,b,c;
  int anr, anc, bnr, bnc;
  
  a = A->nrows;
  b = A->ncols;
  c = B->ncols;

  // adjust cutting numbers to work on words
  a += (a%RADIX) ? RADIX-(a%RADIX) : 0;
  b += (b%RADIX) ? RADIX-(b%RADIX) : 0;
  c += (c%RADIX) ? RADIX-(c%RADIX) : 0;

  anr = (a/RADIX >> 1) * RADIX;
  anc = (b/RADIX >> 1) * RADIX;
  bnr = anc;
  bnc = (c/RADIX >> 1) * RADIX;

  packedmatrix *A00 = mzd_init_window(A,   0,   0,   anr,   anc);
  packedmatrix *A01 = mzd_init_window(A,   0, anc,   anc, 2*anc);
  packedmatrix *A10 = mzd_init_window(A, anr,   0, 2*anr,   anc);
  packedmatrix *A11 = mzd_init_window(A, anr, anc, 2*anr, 2*anc);

  packedmatrix *B00 = mzd_init_window(B,   0,   0,   bnr,   bnc);
  packedmatrix *B01 = mzd_init_window(B,   0, bnc,   bnr, 2*bnc);
  packedmatrix *B10 = mzd_init_window(B, bnr,   0, 2*bnr,   bnc);
  packedmatrix *B11 = mzd_init_window(B, bnr, bnc, 2*bnr, 2*bnc);

  packedmatrix *U0 = mzd_init_window(C,   0,   0,   anr,   bnc);
  packedmatrix *U6 = mzd_init_window(C,   0, bnc,   anr, 2*bnc);
  packedmatrix *U3 = mzd_init_window(C, anr,   0, 2*anr,   bnc);
  packedmatrix *U4 = mzd_init_window(C, anr, bnc, 2*anr, 2*bnc);

  packedmatrix *tmp = mzd_init(7*anr + 4*anc, MAX(anc,bnc));

  int start_row = 0;
  packedmatrix *S0 = mzd_init_window(tmp, start_row, 0, start_row + anr, anc);

  start_row += anr;
  packedmatrix *S1 = mzd_init_window(tmp, start_row, 0, start_row + anr, anc);

  start_row += anr;
  packedmatrix *S2 = mzd_init_window(tmp, start_row, 0, start_row + anr, anc);
  
  start_row += anr;
  packedmatrix *S3 = mzd_init_window(tmp, start_row, 0, start_row + anr, anc);

  start_row += anr;
  packedmatrix *T0 = mzd_init_window(tmp, start_row, 0, start_row + anc, bnc);

  start_row += anc;
  packedmatrix *T1 = mzd_init_window(tmp, start_row, 0, start_row + anc, bnc);
  
  start_row += anc;
  packedmatrix *T2 = mzd_init_window(tmp, start_row, 0, start_row + anc, bnc);

  start_row += anc;
  packedmatrix *T3 = mzd_init_window(tmp, start_row, 0, start_row + anc, bnc);
  
  start_row += anc;
  packedmatrix *Q0 = mzd_init_window(tmp, start_row, 0, start_row + anr, bnc);

  start_row += anr;
  packedmatrix *Q1 = mzd_init_window(tmp, start_row, 0, start_row + anr, bnc);
  
  start_row += anr;
  packedmatrix *Q2 = mzd_init_window(tmp, start_row, 0, start_row + anr, bnc);

  _mzd_add_impl(S0, A10, A11);
  _mzd_add_impl(S1,  S0, A00);
  _mzd_add_impl(S2, A00, A10);
  _mzd_add_impl(S3, A01,  S1);

  _mzd_add_impl(T0, B01, B00);
  _mzd_add_impl(T1, B11,  T0);
  _mzd_add_impl(T2, B11, B01);
  _mzd_add_impl(T3, B10,  T1);

  if (anc <= cutoff || anc <= cutoff || bnc <= cutoff) {
    /* 
     * "This next chunk is arranged so that each output cell gets
     * written to exactly once. This is important because the output
     * blocks might be quite fragmented in memory, whereas our
     * temporary buffers (Q0, Q1, Q2) will be quite localised, so we
     * can afford to do a bit of arithmetic in them." (strassen.pyx)
     */
    
    mzd_mul_m4rm(Q0, A00, B00, 0, NULL, NULL); // now Q0 holds P0
    mzd_mul_m4rm(Q1, A01, B10, 0, NULL, NULL); // now Q1 holds P1
    
    _mzd_add_impl(U0, Q0, Q1); // now U0 is correct
    
    packedmatrix *S1T1 = mzd_mul_m4rm(NULL, S1, T1, 0, NULL, NULL);
    _mzd_add_impl(Q0, Q0, S1T1); // now Q0 holds U1
    mzd_free(S1T1);
    
    mzd_mul_m4rm(Q1, S2, T2, 0, NULL, NULL); // now Q1 holds P4
    
    _mzd_add_impl(Q1, Q1, Q0); // now Q1 holds U2
    mzd_mul_m4rm(Q2, A11, T3, 0, NULL, NULL); // now Q2 holds P6
    _mzd_add_impl(U3, Q1, Q2); // now U3 is correct
    
    mzd_mul_m4rm(Q2, S0, T0, 0, NULL, NULL); // now Q2 holds P2
    _mzd_add_impl(U4, Q2, Q1); // now U4 is correct
    
    _mzd_add_impl(Q0, Q0, Q2); // now Q0 holds U5
    mzd_mul_m4rm(Q2, S3, B11, 0, NULL, NULL); // now Q2 holds P5
    _mzd_add_impl(U6, Q0, Q2);// now U6 is correct

  } else{
    _mzd_mul_strassen_impl(Q0, A00, B00, cutoff); // now Q0 holds P0
    _mzd_mul_strassen_impl(Q1, A01, B10, cutoff); // now Q1 holds P1
    
    _mzd_add_impl(U0, Q0, Q1); // now U0 is correct
    
    packedmatrix *S1T1 = mzd_mul_strassen(NULL, S1, T1, cutoff);
    _mzd_add_impl(Q0, Q0, S1T1); // now Q0 holds U1
    mzd_free(S1T1);
    
    _mzd_mul_strassen_impl(Q1, S2, T2, cutoff); // now Q1 holds P4
    
    _mzd_add_impl(Q1, Q1, Q0); // now Q1 holds U2
    _mzd_mul_strassen_impl(Q2, A11, T3, cutoff); // now Q2 holds P6
    _mzd_add_impl(U3, Q1, Q2); // now U3 is correct
    
    _mzd_mul_strassen_impl(Q2, S0, T0, cutoff); // now Q2 holds P2
    _mzd_add_impl(U4, Q2, Q1); // now U4 is correct
    
    _mzd_add_impl(Q0, Q0, Q2); // now Q0 holds U5
    _mzd_mul_strassen_impl(Q2, S3, B11, cutoff); // now Q2 holds P5
    _mzd_add_impl(U6, Q0, Q2);// now U6 is correct
  }

  // deal with rest
  if (B->ncols > 2*bnc) {
    packedmatrix *B_last_col = mzd_init_window(B, 0, 2*bnc, A->ncols, B->ncols); 
    packedmatrix *C_last_col = mzd_init_window(C, 0, 2*bnc, A->nrows, C->ncols);
    mzd_mul_m4rm(C_last_col, A, B_last_col, 0, NULL, NULL);
    mzd_free_window(B_last_col);
    mzd_free_window(C_last_col);

  }
  if (A->nrows > 2*anr) {
    packedmatrix *A_last_row = mzd_init_window(A, 2*anr, 0, A->nrows, A->ncols);
    packedmatrix *C_last_row = mzd_init_window(C, 2*anr, 0, C->nrows, C->ncols);
    mzd_mul_m4rm(C_last_row, A_last_row, B, 0, NULL, NULL);
    mzd_free_window(A_last_row);
    mzd_free_window(C_last_row);

  }
  if (A->ncols > 2*anc) {
    packedmatrix *A_last_col = mzd_init_window(A,     0, 2*anc, 2*anr, A->ncols);
    packedmatrix *B_last_row = mzd_init_window(B, 2*bnr,     0, B->nrows, 2*bnc);
    packedmatrix *C_bulk = mzd_init_window(C, 0, 0, 2*anr, bnc*2);
    packedmatrix *AB = mzd_mul_m4rm(NULL, A_last_col, B_last_row, 0, NULL, NULL);
    _mzd_add_impl(C_bulk, C_bulk, AB);
    mzd_free(AB);
    mzd_free_window(A_last_col);
    mzd_free_window(B_last_row);
    mzd_free_window(C_bulk);
  }

  /** clean up **/
  mzd_free_window(A00); mzd_free_window(A01);
  mzd_free_window(A10); mzd_free_window(A11);

  mzd_free_window(B00); mzd_free_window(B01);
  mzd_free_window(B10); mzd_free_window(B11);

  mzd_free_window(U0); mzd_free_window(U6);
  mzd_free_window(U3); mzd_free_window(U4);
  
  mzd_free_window(S0); mzd_free_window(S1);
  mzd_free_window(S2); mzd_free_window(S3);

  mzd_free_window(T0); mzd_free_window(T1);
  mzd_free_window(T2); mzd_free_window(T3);

  mzd_free_window(Q0); mzd_free_window(Q1);
  mzd_free_window(Q2);

  mzd_free(tmp);

  return C;
}


packedmatrix *mzd_mul_strassen(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  int k;
  if(A->ncols != B->nrows) {
    m4ri_die("A ncols need to match B nrows.\n");
  }

  if (cutoff <= 0) {
    m4ri_die("cutoff must be > 0.\n");
  }
  cutoff = cutoff/RADIX * RADIX;
  if (cutoff == 0) {
    cutoff = RADIX;
  };

  if (C == NULL) {
    C = mzd_init(A->nrows, B->ncols);
  } else {
    if (C->nrows != A->nrows || C->ncols != B->ncols) {
      m4ri_die("C has wrong dimensions.\n");
    }
  }

  /** handle case first, where the input matrices are too small already */
  if (A->nrows <= cutoff || A->ncols <= cutoff || B->ncols <= cutoff) {
    k = m4ri_opt_k(A->nrows, A->ncols, B->ncols);
    C = mzd_mul_m4rm(C, A, B, k, NULL, NULL);
    return C;
  }

  return _mzd_mul_strassen_impl(C, A, B, cutoff);
}
