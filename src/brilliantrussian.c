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

#include "brilliantrussian.h"
#include <emmintrin.h>
#include "misc.h"

#define TWOPOW(i) (1<<(i))

int __do_print = 0;

int forceNonZero(packedmatrix *m, int xstart, int xstop, int y) {
  int i;

  for (i=xstart; i<=xstop; i++) {
    if (readCell(m, i, y)==1) {
      if (i!=xstart) rowSwap(m, i, xstart);
      return YES;
    }
  }
  return NO;
}

int prep(packedmatrix *m, int ai, int k) {
  int pc; /* pivot column */
  int tr; /* target row */
  int good;

  int rank = 0;

  for (pc=ai; pc<min(ai+k,m->ncols); pc++) {
    /* Step one, find a pivot row in this column.*/
    good=forceNonZero(m, pc, min( ai+k*3-1, m->nrows-1 ), pc);

    if (good==NO) return rank;

    for (tr=ai; tr<min(ai+k*3, m->nrows); tr++) {
      /* Step two, add this pivot row to other rows as needed. */
      if (tr==pc) continue;
      
      if (readCell(m, tr, pc)==0) continue;

      rowAddOffset(m, pc, tr, ai);
    }
    rank++;
  }

  return rank;
}

void m2t_combine( packedmatrix * sc1, int row1, int startblock1, 
		  packedmatrix * sc2, int row2, int startblock2,
		  packedmatrix * dst, int row3, int startblock3 ) {
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
    /** check alignments **/
    if (((unsigned long)b1_ptr%16==8) && ((unsigned long)b2_ptr%16==8)) {
      *b1_ptr++ ^= *b2_ptr++;
      wide--;
    }

    if (wide>32 && ((unsigned long)b1_ptr%16==0) && ((unsigned long)b2_ptr%16==0)) {
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

      // handle rest
      if (wide & 1)
	*b1_ptr ^= *b2_ptr;

    } else {
      int i;
      for(i = wide>>1 ; i > 0 ; i--) {
	*b1_ptr++ ^= *b2_ptr++;
	*b1_ptr++ ^= *b2_ptr++;
      }
      if (wide & 1)
	*b1_ptr ^= *b2_ptr;
    }

#else // no SSE2

    int i;
    for(i = wide>>1 ; i > 0 ; i--) {
      *b1_ptr++ ^= *b2_ptr++;
      *b1_ptr++ ^= *b2_ptr++;
    }
    if (wide & 1)
      *b1_ptr ^= *b2_ptr;

#endif
    
  } else { // dst != sc1
    int i;
    b3_ptr = dst->values + startblock3 + dst->rowswap[row3];

    if (row1 >= sc1->nrows) {
	for(i = 0 ; i < wide ; i++) {
	  *b3_ptr++ = *b2_ptr++;
	}
    } else {
      for(i = 0 ; i < wide ; i++) {
	*b3_ptr++ = *b1_ptr++ ^ *b2_ptr++;
      }
    }
    return;
  }
}

void makeTable( packedmatrix *m, int ai, int k,
		packedmatrix *tablepacked, int *lookuppacked, int full) {
  int homeblock= full ? 0 : ai/RADIX;
  int i, rowneeded, id;
  int twokay= TWOPOW(k);

  lookuppacked[0]=0;

  for (i=1; i<twokay; i++) {
    rowneeded=codebook[k]->inc[i-1]+ai;

    id=codebook[k]->ord[i];

    lookuppacked[id]=i;

    m2t_combine( m,            rowneeded, homeblock, 
		 tablepacked,  i-1,       homeblock, 
		 tablepacked,  i,         homeblock);
  }
  if (__do_print) {
    printf("Table Packed\n");
    printMatrix(tablepacked);
  }
}



static inline int getBits(packedmatrix *m, int x, int y, int k) {
  int truerow = m->rowswap[x];
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

void processRow(packedmatrix *m, int row, int homecol, int k, packedmatrix *tablepacked, int *lookuppacked) {
  int blocknum = homecol/RADIX;

  int value = getBits(m, row, homecol, k);

  int tablerow = lookuppacked[value];

  m2t_combine(m,                row, blocknum,
	      tablepacked, tablerow, blocknum,
	      m,                row, blocknum);
}

void process(packedmatrix *m, int startrow, int stoprow, int startcol, int k, packedmatrix *tablepacked, int *lookuppacked) {
  int i;
  int blocknum=startcol/RADIX;
  int value;
  int tablerow;
  word *b1_ptr,*b2_ptr;

  // for optimization reasons we distinguish several cases here.

  switch(m->width - startcol/RADIX) {

  case 1:
    // no loop needed as only one block is operated on.
    for (i=startrow; i<=stoprow; i++) {
      value = getBits(m, i, startcol, k);
      tablerow = lookuppacked[value];
      b1_ptr = m->values + blocknum + m->rowswap[i];
      b2_ptr = tablepacked->values + blocknum + tablepacked->rowswap[tablerow];
      *b1_ptr ^= *b2_ptr;
    }
    break;

  case 2:
    // two blocks, no loop
    for (i=startrow; i<=stoprow; i++) {
      value = getBits(m, i, startcol, k);
      tablerow = lookuppacked[value];
      b1_ptr = m->values + blocknum + m->rowswap[i];
      b2_ptr = tablepacked->values + blocknum + tablepacked->rowswap[tablerow];
      *b1_ptr++ ^= *b2_ptr++;
      *b1_ptr ^= *b2_ptr;
    }
    break;

  default:
    // the real deal more than two blocks.
    for (i=startrow; i<=stoprow; i++) {
      processRow(m, i, startcol, k, tablepacked, lookuppacked);
    }
    break;
  }
}

int stepM4RI(packedmatrix *m, int full, int k, int ai, 
	     packedmatrix *tablepacked, int *lookuppacked) {
  int submatrixrank;
  
  /*
   * Stage 1: Denote the first column to be processed in a given
   * iteration as ai. Then, perform Gaussian elimination on the
   * first 3k rows after and including the i-th row to produce an
   *
   * identity matrix in $a_{(i,i)} ... a_{(i+k-1),(i+k-1)},$ 
   *
   * and zeroes in $a_{(i+k),i} ... a_{(i+3k-1),(i+k-1)}$.
   */

  submatrixrank=prep(m, ai, k);


  if (submatrixrank!=k) return submatrixrank;

  /*
   * Stage 2: Construct a table consisting of the 2k binary strings of
   * length k in a Gray Code.  Thus with only 2k vector additions, all
   * possible linear combinations of these k rows have been
   * precomputed.
   */

  makeTable(m, ai, k, tablepacked, lookuppacked, 0);


  /*
   * Stage 3: One can rapidly process the remaining rows from i + 3k
   * until row m (the last row) by using the table. For example,
   * suppose the jth row has entries $a_{(j,i)} ... a_{(j,i+k-1)}$ in
   * the columns being processed. Selecting the row of the table
   * associated with this k-bit string, and adding it to row j will
   * force the k columns to zero, and adjust the remaining columns
   * from i + k to n in the appropriate way, as if Gaussian
   * Elimination had been performed.
  */

  process(m, ai+k*3, m->nrows-1, ai, k,
	  tablepacked, lookuppacked);

  /* While the above form of the algorithm will reduce a system of
   * boolean linear equations to unit upper triangular form, and thus
   * permit a system to be solved with back substitution, the M4RI
   * algorithm can also be used to invert a matrix, or put the system
   * into reduced row echelon form (RREF). Simply run Stage 3 on rows
   * 0 ... i - 1 as well as on rows i + 3k · · · m. This only affects
   * the complexity slightly, changing the 2.5 coeffcient to 3
   */

  if (full==YES) process(m, 0, ai-1, ai, k, 
			 tablepacked, lookuppacked);

  return submatrixrank;
}

int reduceM4RI(packedmatrix *m, int full, int k, packedmatrix *tablepacked, int *lookuppacked) {
  int i, submatrixrank;
  int stop = min(m->nrows, m->ncols); 

  int rank = 0;
  int simple = 0;

  if (k == 0) {
    k = m2t_opt_k(m->nrows, m->ncols, 0);
  }

  if (tablepacked == NULL && lookuppacked == NULL) {
    simple = 1;
    tablepacked = m2t_init( TWOPOW(k), m->ncols );
    lookuppacked = (int *)safeCalloc( TWOPOW(k), sizeof(int) );
  }
  
  // main loop
  for (i=0; i<stop; i+=k) {
    // not enough room for M4RI left.
    if ( ((i+k*3) > m->nrows) || ((i+k) > m->ncols) ) {
      rank += reduceGaussianDelayed(m, i, full);
      break;
    }
    
    submatrixrank=stepM4RI(m, full, k, i, tablepacked, lookuppacked);

    if (submatrixrank!=k) {
      // not full rank, use Gaussian elimination :-(
      rank += reduceGaussianDelayed(m, i, full);
      break;
    }
    rank += submatrixrank;
  }

  if (simple) {
    free(lookuppacked);
    m2t_free(tablepacked);
  }
  return rank; 
}

void topReduceM4RI(packedmatrix *m, int k, packedmatrix *tablepacked, int *lookuppacked) {
  int i, submatrixrank;
  int stop = min(m->nrows, m->ncols); 
  int simple = 0;
  
  if (k == 0) {
    k = m2t_opt_k(m->nrows, m->ncols, 0);
  }
  
  /* setup tables */
  if (tablepacked == NULL && lookuppacked == NULL) {
    simple = 1;
    tablepacked = m2t_init( TWOPOW(k), m->ncols );
    lookuppacked = (int *)safeCalloc( TWOPOW(k), sizeof(int) );
  }
  
  /* main loop */
  for (i=0; i<stop; i+=k) {
    if ( (i+k > m->nrows) || (i+k > m->ncols) ) {
      reduceGaussianDelayed(m, i, 1);
      break;
    }
    
    submatrixrank = prep(m, i, k);
    
    if (submatrixrank==k) {
      makeTable(m, i, k, tablepacked, lookuppacked, 0);
      process(m, 0, i-1, i, k, tablepacked, lookuppacked);
    } else {
      reduceGaussianDelayed(m, i, 1);
      break;
    }
  }
  
  /* clear tables */
  if (simple) {
    free(lookuppacked);
    m2t_free(tablepacked);
  }
}

packedmatrix *invertM4RI(packedmatrix *m, 
			 packedmatrix *identity, int k) {
  packedmatrix *big=concat(m, identity);
  int size=m->ncols;
  if (k == 0) {
    k = m2t_opt_k(m->nrows, m->ncols, 0);
  }
  int twokay=TWOPOW(k);
  int i;
  packedmatrix *mytable=m2t_init(twokay, size*2);
  int *mylookups=(int *)safeCalloc(twokay, sizeof(int));
  packedmatrix *answer;
  
  reduceM4RI(big, YES, k, mytable, mylookups);
  
  for(i=0; i < size; i++) {
    if (!readCell(big, i,i )) {
      answer=NULL;
      break;
    }
  }
  if (i == size)
    answer=m2t_submatrix(big, 0, size, size, size*2);
  
  free(mylookups);
  m2t_free(mytable);
  m2t_free(big);
  
  return answer;
}

packedmatrix *m2t_mul_m4rm_t(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k) {
  packedmatrix *AT, *BT, *CT;
  
  //if(A->ncols != B->nrows) die("A cols need to match B rows");
  
  AT = m2t_transpose(NULL, A);
  BT = m2t_transpose(NULL, B);
  
  CT = m2t_mul_m4rm(NULL, BT,AT,k, NULL, NULL);
  
  m2t_free(AT);
  m2t_free(BT);

  C = m2t_transpose(C, CT);
  m2t_free(CT);
  return C;
}

packedmatrix *m2t_mul_m4rm(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k, packedmatrix *tablepacked, int *lookuppacked) {
  int i,j,  a,b,c, simple;
  unsigned int x;
  
  //if(A->ncols != B->nrows) 
  //  die("A cols need to match B rows");
  
  a = A->nrows;
  b = A->ncols;
  c = B->ncols;


  if (C == NULL) {
    C = m2t_init(a, c);
  } else {
    //if (C->nrows != a || C->ncols != c) {
    //  die("C has wrong dimensions.\n");
    //}
    /** clear first **/
    for (i=0; i<C->nrows; i++) {
      for (j=0; j<C->width; j++) {
	C->values[C->rowswap[i] + j] = 0;
      }
    }
  }

  if (k == 0) {
    k = m2t_opt_k(a,b,c);
    //printf("k: %d\n",k);
  }

  if (tablepacked == NULL && lookuppacked == NULL) {
    simple = 1;
    tablepacked = m2t_init(TWOPOW(k), c);
    lookuppacked = (int *)safeCalloc(TWOPOW(k), sizeof(int));
  }

  for(i=0; i < b/k; i++) {

    /* Make a Gray Code table of all the $2^k$ linear combinations of
     * the $k$ rows of $B_i$.  Call the $x$-th row $T_x$.
     */
    makeTable( B, i*k, k, tablepacked, lookuppacked, 1 );
    
    for(j = 0; j<a; j++) {

      /* Read the entries 
       *
       *  $aj,(i-1)k+1, aj,(i-1)k+2 , ... , aj,(i-1)k+k.$
       *
       * Let $x$ be the $k$ bit binary number formed by the
       * concatenation of $aj,(i-1)k+1 , ... , aj,ik$.
       */

      x = lookuppacked[getBits(A, j, i*k, k)];

      /* for h = 1, 2, . . . , c do
       *   Calculate Cjh = Cjh + Txh.
       */
      m2t_combine( C,j,0,  tablepacked,x,0,  C,j,0 );
    }
  }

  //handle rest
  if (b%k) {
    makeTable( B, b/k * k , b%k, tablepacked, lookuppacked, 1);
    
    for(j = 0; j<a; j++) {
      x = lookuppacked[getBits(A, j, i*k, b%k)];
      m2t_combine(C,j,0, tablepacked,x,0,  C,j,0);
    }
  }
  if (simple) {
    m2t_free(tablepacked);
    free(lookuppacked);
  }
  return C;
}


static packedmatrix *_m2t_mul_strassen_impl(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
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

  packedmatrix *A00 = m2t_init_window(A,   0,   0,   anr,   anc);
  packedmatrix *A01 = m2t_init_window(A,   0, anc,   anc, 2*anc);
  packedmatrix *A10 = m2t_init_window(A, anr,   0, 2*anr,   anc);
  packedmatrix *A11 = m2t_init_window(A, anr, anc, 2*anr, 2*anc);

  packedmatrix *B00 = m2t_init_window(B,   0,   0,   bnr,   bnc);
  packedmatrix *B01 = m2t_init_window(B,   0, bnc,   bnr, 2*bnc);
  packedmatrix *B10 = m2t_init_window(B, bnr,   0, 2*bnr,   bnc);
  packedmatrix *B11 = m2t_init_window(B, bnr, bnc, 2*bnr, 2*bnc);

  packedmatrix *U0 = m2t_init_window(C,   0,   0,   anr,   bnc);
  packedmatrix *U6 = m2t_init_window(C,   0, bnc,   anr, 2*bnc);
  packedmatrix *U3 = m2t_init_window(C, anr,   0, 2*anr,   bnc);
  packedmatrix *U4 = m2t_init_window(C, anr, bnc, 2*anr, 2*bnc);

  packedmatrix *tmp = m2t_init(7*anr + 4*anc, max(anc,bnc));

  int start_row = 0;
  packedmatrix *S0 = m2t_init_window(tmp, start_row, 0, start_row + anr, anc);

  start_row += anr;
  packedmatrix *S1 = m2t_init_window(tmp, start_row, 0, start_row + anr, anc);

  start_row += anr;
  packedmatrix *S2 = m2t_init_window(tmp, start_row, 0, start_row + anr, anc);
  
  start_row += anr;
  packedmatrix *S3 = m2t_init_window(tmp, start_row, 0, start_row + anr, anc);

  start_row += anr;
  packedmatrix *T0 = m2t_init_window(tmp, start_row, 0, start_row + anc, bnc);

  start_row += anc;
  packedmatrix *T1 = m2t_init_window(tmp, start_row, 0, start_row + anc, bnc);
  
  start_row += anc;
  packedmatrix *T2 = m2t_init_window(tmp, start_row, 0, start_row + anc, bnc);

  start_row += anc;
  packedmatrix *T3 = m2t_init_window(tmp, start_row, 0, start_row + anc, bnc);
  
  start_row += anc;
  packedmatrix *Q0 = m2t_init_window(tmp, start_row, 0, start_row + anr, bnc);

  start_row += anr;
  packedmatrix *Q1 = m2t_init_window(tmp, start_row, 0, start_row + anr, bnc);
  
  start_row += anr;
  packedmatrix *Q2 = m2t_init_window(tmp, start_row, 0, start_row + anr, bnc);

  _m2t_add_impl(S0, A10, A11);
  _m2t_add_impl(S1,  S0, A00);
  _m2t_add_impl(S2, A00, A10);
  _m2t_add_impl(S3, A01,  S1);

  _m2t_add_impl(T0, B01, B00);
  _m2t_add_impl(T1, B11,  T0);
  _m2t_add_impl(T2, B11, B01);
  _m2t_add_impl(T3, B10,  T1);

  if (anc <= cutoff || anc <= cutoff || bnc <= cutoff) {
    /* 
     * "This next chunk is arranged so that each output cell gets
     * written to exactly once. This is important because the output
     * blocks might be quite fragmented in memory, whereas our
     * temporary buffers (Q0, Q1, Q2) will be quite localised, so we
     * can afford to do a bit of arithmetic in them." (strassen.pyx)
     */
    
    m2t_mul_m4rm(Q0, A00, B00, 0, NULL, NULL); // now Q0 holds P0
    m2t_mul_m4rm(Q1, A01, B10, 0, NULL, NULL); // now Q1 holds P1
    
    _m2t_add_impl(U0, Q0, Q1); // now U0 is correct
    
    packedmatrix *S1T1 = m2t_mul_m4rm(NULL, S1, T1, 0, NULL, NULL);
    _m2t_add_impl(Q0, Q0, S1T1); // now Q0 holds U1
    m2t_free(S1T1);
    
    m2t_mul_m4rm(Q1, S2, T2, 0, NULL, NULL); // now Q1 holds P4
    
    _m2t_add_impl(Q1, Q1, Q0); // now Q1 holds U2
    m2t_mul_m4rm(Q2, A11, T3, 0, NULL, NULL); // now Q2 holds P6
    _m2t_add_impl(U3, Q1, Q2); // now U3 is correct
    
    m2t_mul_m4rm(Q2, S0, T0, 0, NULL, NULL); // now Q2 holds P2
    _m2t_add_impl(U4, Q2, Q1); // now U4 is correct
    
    _m2t_add_impl(Q0, Q0, Q2); // now Q0 holds U5
    m2t_mul_m4rm(Q2, S3, B11, 0, NULL, NULL); // now Q2 holds P5
    _m2t_add_impl(U6, Q0, Q2);// now U6 is correct

  } else{
    _m2t_mul_strassen_impl(Q0, A00, B00, cutoff); // now Q0 holds P0
    _m2t_mul_strassen_impl(Q1, A01, B10, cutoff); // now Q1 holds P1
    
    _m2t_add_impl(U0, Q0, Q1); // now U0 is correct
    
    packedmatrix *S1T1 = m2t_mul_strassen(NULL, S1, T1, cutoff);
    _m2t_add_impl(Q0, Q0, S1T1); // now Q0 holds U1
    m2t_free(S1T1);
    
    _m2t_mul_strassen_impl(Q1, S2, T2, cutoff); // now Q1 holds P4
    
    _m2t_add_impl(Q1, Q1, Q0); // now Q1 holds U2
    _m2t_mul_strassen_impl(Q2, A11, T3, cutoff); // now Q2 holds P6
    _m2t_add_impl(U3, Q1, Q2); // now U3 is correct
    
    _m2t_mul_strassen_impl(Q2, S0, T0, cutoff); // now Q2 holds P2
    _m2t_add_impl(U4, Q2, Q1); // now U4 is correct
    
    _m2t_add_impl(Q0, Q0, Q2); // now Q0 holds U5
    _m2t_mul_strassen_impl(Q2, S3, B11, cutoff); // now Q2 holds P5
    _m2t_add_impl(U6, Q0, Q2);// now U6 is correct
  }

  // deal with rest
  if (B->ncols > 2*bnc) {
    packedmatrix *B_last_col = m2t_init_window(B, 0, 2*bnc, A->ncols, B->ncols); 
    packedmatrix *C_last_col = m2t_init_window(C, 0, 2*bnc, A->nrows, C->ncols);
    m2t_mul_m4rm(C_last_col, A, B_last_col, 0, NULL, NULL);
    m2t_free_window(B_last_col);
    m2t_free_window(C_last_col);

  }
  if (A->nrows > 2*anr) {
    packedmatrix *A_last_row = m2t_init_window(A, 2*anr, 0, A->nrows, A->ncols);
    packedmatrix *C_last_row = m2t_init_window(C, 2*anr, 0, C->nrows, C->ncols);
    m2t_mul_m4rm(C_last_row, A_last_row, B, 0, NULL, NULL);
    m2t_free_window(A_last_row);
    m2t_free_window(C_last_row);

  }
  if (A->ncols > 2*anc) {
    packedmatrix *A_last_col = m2t_init_window(A,     0, 2*anc, 2*anr, A->ncols);
    packedmatrix *B_last_row = m2t_init_window(B, 2*bnr,     0, B->nrows, 2*bnc);
    packedmatrix *C_bulk = m2t_init_window(C, 0, 0, 2*anr, bnc*2);
    packedmatrix *AB = m2t_mul_m4rm(NULL, A_last_col, B_last_row, 0, NULL, NULL);
    _m2t_add_impl(C_bulk, C_bulk, AB);
    m2t_free(AB);
    m2t_free_window(A_last_col);
    m2t_free_window(B_last_row);
    m2t_free_window(C_bulk);
  }

  /** clean up **/
  m2t_free_window(A00); m2t_free_window(A01);
  m2t_free_window(A10); m2t_free_window(A11);

  m2t_free_window(B00); m2t_free_window(B01);
  m2t_free_window(B10); m2t_free_window(B11);

  m2t_free_window(U0); m2t_free_window(U6);
  m2t_free_window(U3); m2t_free_window(U4);
  
  m2t_free_window(S0); m2t_free_window(S1);
  m2t_free_window(S2); m2t_free_window(S3);

  m2t_free_window(T0); m2t_free_window(T1);
  m2t_free_window(T2); m2t_free_window(T3);

  m2t_free_window(Q0); m2t_free_window(Q1);
  m2t_free_window(Q2);

  m2t_free(tmp);

  return C;
}


packedmatrix *m2t_mul_strassen(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff) {
  int k;
  if(A->ncols != B->nrows) {
    die("A ncols need to match B nrows.\n");
  }

  if (cutoff <= 0) {
    die("cutoff must be > 0.\n");
  }
  cutoff = cutoff/RADIX * RADIX;
  if (cutoff == 0) {
    cutoff = RADIX;
  };

  if (C == NULL) {
    C = m2t_init(A->nrows, B->ncols);
  } else {
    if (C->nrows != A->nrows || C->ncols != B->ncols) {
      die("C has wrong dimensions.\n");
    }
  }

  /** handle case first, where the input matrices are too small already */
  if (A->nrows <= cutoff || A->ncols <= cutoff || B->ncols <= cutoff) {
    k = m2t_opt_k(A->nrows, A->ncols, B->ncols);
    C = m2t_mul_m4rm(C, A, B, k, NULL, NULL);
    return C;
  }

  return _m2t_mul_strassen_impl(C, A, B, cutoff);
}
