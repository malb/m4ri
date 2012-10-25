#include "triangular_russian.h"
#include "graycode.h"
#include "brilliantrussian.h"
#include "ple_russian.h"
#include "xor.h"

void _mzd_trsm_upper_left_submatrix(mzd_t const *U, mzd_t *B, rci_t const start_row, int const k, word const mask_begin, word const mask_end) {
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < i; ++j) {
      if (mzd_read_bit(U, start_row+(k-i-1), start_row+(k-i)+j)) {
        word *a = B->rows[start_row+(k-i-1)];
        word *b = B->rows[start_row+(k-i)+j];

        *a++ ^= *b++ & mask_begin;
	wi_t ii;
        for(ii = 1; ii + 8 <= B->width - 1; ii += 8) {
          *a++ ^= *b++;
          *a++ ^= *b++;
          *a++ ^= *b++;
          *a++ ^= *b++;
          *a++ ^= *b++;
          *a++ ^= *b++;
          *a++ ^= *b++;
          *a++ ^= *b++;
        }
        switch(B->width - ii) {
        case 8:  *a++ ^= *b++;
        case 7:  *a++ ^= *b++;
        case 6:  *a++ ^= *b++;
        case 5:  *a++ ^= *b++;
        case 4:  *a++ ^= *b++;
        case 3:  *a++ ^= *b++;
        case 2:  *a++ ^= *b++;
        case 1:  *a++ ^= (*b++ & mask_end);
        }
      }
    }
  }

  __M4RI_DD_MZD(B);
}

void _mzd_trsm_upper_left_russian(mzd_t const *U, mzd_t *B, int k) {
  wi_t const wide = B->width;
  int const blocksize = __M4RI_MUL_BLOCKSIZE;

  word mask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix - B->offset);
  word mask_end = __M4RI_LEFT_BITMASK((B->ncols + B->offset) % m4ri_radix);

  if (B->width == 1)
    mask_begin = mask_begin & mask_end;

  if (k == 0) {
    k = m4ri_opt_k(blocksize, B->nrows, B->ncols);
    if (k > 3)
      k -= 2;
    /* reduce k further if that has a chance of hitting L1 */
    size_t const tsize = (int)(0.8 * (__M4RI_TWOPOW(k) * B->nrows));
    if(__M4RI_CPU_L1_CACHE < tsize && tsize <= 2 * __M4RI_CPU_L1_CACHE)
      k -= 1;
  }

  mzd_t *T0 = mzd_init(__M4RI_TWOPOW(k), B->ncols + B->offset);
  mzd_t *T1 = mzd_init(__M4RI_TWOPOW(k), B->ncols + B->offset);
  mzd_t *T2 = mzd_init(__M4RI_TWOPOW(k), B->ncols + B->offset);
  mzd_t *T3 = mzd_init(__M4RI_TWOPOW(k), B->ncols + B->offset);
  mzd_t *T4 = mzd_init(__M4RI_TWOPOW(k), B->ncols + B->offset);
  mzd_t *T5 = mzd_init(__M4RI_TWOPOW(k), B->ncols + B->offset);
  mzd_t *T6 = mzd_init(__M4RI_TWOPOW(k), B->ncols + B->offset);
  mzd_t *T7 = mzd_init(__M4RI_TWOPOW(k), B->ncols + B->offset);


  rci_t *L0 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L1 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L2 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L3 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L4 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L5 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L6 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *L7 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));

  int kk = 8 * k;

  rci_t i = 0;
  for (; i < B->nrows - kk; i += kk) {

    _mzd_trsm_upper_left_submatrix(U, B, B->nrows-i-kk, kk, mask_begin, mask_end);

    mzd_make_table(B, B->nrows - i - 8*k, 0, k, T7, L7);
    mzd_make_table(B, B->nrows - i - 7*k, 0, k, T6, L6);
    mzd_make_table(B, B->nrows - i - 6*k, 0, k, T5, L5);
    mzd_make_table(B, B->nrows - i - 5*k, 0, k, T4, L4);
    mzd_make_table(B, B->nrows - i - 4*k, 0, k, T3, L3);
    mzd_make_table(B, B->nrows - i - 3*k, 0, k, T2, L2);
    mzd_make_table(B, B->nrows - i - 2*k, 0, k, T1, L1);
    mzd_make_table(B, B->nrows - i - 1*k, 0, k, T0, L0);

    for(rci_t j = 0; j < B->nrows - i - kk; ++j) {
      rci_t const x7 = L7[ mzd_read_bits_int(U, j, B->nrows - i - 8*k, k) ];
      rci_t const x6 = L6[ mzd_read_bits_int(U, j, B->nrows - i - 7*k, k) ];
      rci_t const x5 = L5[ mzd_read_bits_int(U, j, B->nrows - i - 6*k, k) ];
      rci_t const x4 = L4[ mzd_read_bits_int(U, j, B->nrows - i - 5*k, k) ];
      rci_t const x3 = L3[ mzd_read_bits_int(U, j, B->nrows - i - 4*k, k) ];
      rci_t const x2 = L2[ mzd_read_bits_int(U, j, B->nrows - i - 3*k, k) ];
      rci_t const x1 = L1[ mzd_read_bits_int(U, j, B->nrows - i - 2*k, k) ];
      rci_t const x0 = L0[ mzd_read_bits_int(U, j, B->nrows - i - 1*k, k) ];


      word *b = B->rows[j];
      word *t7 = T7->rows[x7];
      word *t6 = T6->rows[x6];
      word *t5 = T5->rows[x5];
      word *t4 = T4->rows[x4];
      word *t3 = T3->rows[x3];
      word *t2 = T2->rows[x2];
      word *t1 = T1->rows[x1];
      word *t0 = T0->rows[x0];

      _mzd_combine8(b, t0, t1, t2, t3, t4, t5, t6, t7, wide);
    }
  }

  /* handle stuff that doesn't fit in multiples of kk */
  for ( ;i < B->nrows; i += k) {
    if (i > B->nrows - k)
      k = B->nrows - i;

    _mzd_trsm_upper_left_submatrix(U, B, B->nrows-i-k, k, mask_begin, mask_end);

    mzd_make_table(B, B->nrows - i - 1*k, 0, k, T0, L0);

    for(rci_t j = 0; j < B->nrows - i - k; ++j) {
      rci_t const x0 = L0[ mzd_read_bits_int(U, j, B->nrows - i - 1*k, k) ];

      word *b = B->rows[j];
      word *t0 = T0->rows[x0];

      for (wi_t ii = 0; ii < wide; ++ii)
        b[ii] ^= t0[ii];
    }
  }

  mzd_free(T0);
  mzd_free(T1);
  mzd_free(T2);
  mzd_free(T3);
  mzd_free(T4);
  mzd_free(T5);
  mzd_free(T6);
  mzd_free(T7);

  m4ri_mm_free(L0);
  m4ri_mm_free(L1);
  m4ri_mm_free(L2);
  m4ri_mm_free(L3);
  m4ri_mm_free(L4);
  m4ri_mm_free(L5);
  m4ri_mm_free(L6);
  m4ri_mm_free(L7);

  __M4RI_DD_MZD(B);
}

void mzd_make_table_trtri(mzd_t const *M, rci_t r, rci_t c, int k, mzd_t *T, rci_t *L) {
  assert(!(T->flags & mzd_flag_multiple_blocks));
  wi_t const blockoffset= c / m4ri_radix;
  int const twokay= __M4RI_TWOPOW(k);
  wi_t const wide = T->width - blockoffset;
  wi_t const count = (wide + 7) / 8;
  int const entry_point = wide % 8;
  wi_t const next_row_offset = blockoffset + T->rowstride - T->width;

  word *ti, *ti1, *m;

  ti1 = T->rows[0] + blockoffset;
  ti = ti1 + T->rowstride;

  L[0] = 0;
  for (int i = 1; i < twokay; ++i) {
    rci_t rowneeded = r + m4ri_codebook[k]->inc[i - 1];
    m = M->rows[rowneeded] + blockoffset;

    wi_t n = count;
    switch (entry_point) {
    case 0: do { *(ti++) = *(m++) ^ *(ti1++);
    case 7:      *(ti++) = *(m++) ^ *(ti1++);
    case 6:      *(ti++) = *(m++) ^ *(ti1++);
    case 5:      *(ti++) = *(m++) ^ *(ti1++);
    case 4:      *(ti++) = *(m++) ^ *(ti1++);
    case 3:      *(ti++) = *(m++) ^ *(ti1++);
    case 2:      *(ti++) = *(m++) ^ *(ti1++);
    case 1:      *(ti++) = *(m++) ^ *(ti1++);
      } while (--n > 0);
    }
    ti += next_row_offset;
    ti1 += next_row_offset;

    L[m4ri_codebook[k]->ord[i]] = i;
  }

  for(int i=1; i<twokay; ++i)
    mzd_xor_bits(T, i, c, k, (word)m4ri_codebook[k]->ord[i]);
}

#define __M4RI_TRTRI_NTABLES 4

static inline void _mzd_trtri_upper_submatrix(mzd_t *A, rci_t pivot_r, rci_t elim_r, const int k) {
  for(rci_t i=pivot_r; i<pivot_r+k; i++)
    for(rci_t j=elim_r; j<i; j++)
      if(mzd_read_bit(A,j,i) && (i+1)<A->ncols )
        mzd_row_add_offset(A, j, i, i+1);
}


mzd_t *mzd_trtri_upper_russian(mzd_t *A, int k) {
  assert(A->nrows == A->ncols && A->offset == 0);

  if (k == 0) {
    k = m4ri_opt_k(A->nrows, A->ncols, 0);
    if (k >= 7)
      k = 7;
    if (0.75 * __M4RI_TWOPOW(k) *A->ncols > __M4RI_CPU_L3_CACHE / 2.0)
      k -= 1;
  }

  const int kk = __M4RI_TRTRI_NTABLES*k;

  mzd_t *T[__M4RI_TRTRI_NTABLES];
  rci_t *L[__M4RI_TRTRI_NTABLES];
  mzd_t *U[__M4RI_TRTRI_NTABLES];
  for(int i=0; i<__M4RI_TRTRI_NTABLES; i++) {
    T[i] = mzd_init(__M4RI_TWOPOW(k), A->ncols);
    L[i] = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
    U[i] = mzd_init(k, A->ncols);
  }

  /** dummy offsets table for make_table_ple**/
  rci_t id[m4ri_radix];
  for(int i=0; i<m4ri_radix; i++) id[i] = i;

  rci_t r = 0;
  while(r+kk <= A->nrows) {

    /***
     * ----------------------------
     * [  ....................... ]
     * [  ... U00 U01 U02 U03 ... ]
     * [  ...     U10 U12 U13 ... ]
     * ---------------------------- r
     * [  ...         U22 U23 ... ]
     * [  ...             U33 ... ]
     * ----------------------------
     *
     * Assume [ U00 U01 ] was already inverted and multiplied with [ U02 U03 ... ]
     *        [     U10 ]                                          [ U12 U13 ... ]
     *
     * We then invert U22 and construct a table for [U22 U23 ... ], then we
     * invert [U33] and multiply it with [U23]. Then we construct a table for [U23 ... ]
     **/

    _mzd_trtri_upper_submatrix(A, r, r, k);
    _mzd_ple_to_e(U[0], A, r, r, k, id);
    mzd_make_table_trtri(U[0], 0, r,   k, T[0], L[0]);

    _mzd_trtri_upper_submatrix(A, r+k, r, k);
    _mzd_ple_to_e(U[1], A, r+k, r+k, k, id);
    mzd_make_table_trtri(U[1], 0, r+k, k, T[1], L[1]);

    _mzd_trtri_upper_submatrix(A, r+2*k, r, k);
    _mzd_ple_to_e(U[2], A, r+2*k, r+2*k, k, id);
    mzd_make_table_trtri(U[2], 0, r+2*k, k, T[2], L[2]);

    _mzd_trtri_upper_submatrix(A, r+3*k, r, k);
    _mzd_ple_to_e(U[3], A, r+3*k, r+3*k, k, id);
    mzd_make_table_trtri(U[3], 0, r+3*k, k, T[3], L[3]);

    mzd_process_rows4_ple(A, 0, r, r,
                          k, T[0], L[0],
                          k, T[1], L[1],
                          k, T[2], L[2],
                          k, T[3], L[3]);
    r += kk;
  }

  /** deal with the rest **/
  while(r < A->nrows) {
    if (A->nrows - r < k)
      k = A->nrows - r;
    for(rci_t i=0; i<k; i++)
      for(rci_t j=0; j<i; j++)
        if(mzd_read_bit(A,r+j,r+i) && (r+i+1)<A->ncols )
          mzd_row_add_offset(A, r+j, r+i, r+i+1);

    _mzd_ple_to_e(U[0], A, r, r, k, id);
    mzd_make_table_trtri(U[0], 0, r, k, T[0], L[0]);

    mzd_process_rows(A, 0, r, r, k, T[0], L[0]);
    r += k;
  }

  for(int i=0; i<__M4RI_TRTRI_NTABLES; i++) {
    mzd_free(T[i]);
    m4ri_mm_free(L[i]);
    mzd_free(U[i]);
  }
  __M4RI_DD_MZD(A);
  return A;
}
