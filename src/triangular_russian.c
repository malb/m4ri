#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "triangular_russian.h"
#include "graycode.h"
#include "brilliantrussian.h"
#include "ple_russian.h"
#include "xor.h"

/** the number of tables used in TRSM decomposition **/
#define __M4RI_TRSM_NTABLES 8

void _mzd_trsm_upper_left_submatrix(mzd_t const *U, mzd_t *B, rci_t const start_row, int const k, word const mask_end) {
  for (rci_t i = 0; i < k; ++i) {
    for (rci_t j = 0; j < i; ++j) {
      if (mzd_read_bit(U, start_row+(k-i-1), start_row+(k-i)+j)) {
        word *a = B->rows[start_row+(k-i-1)];
        word *b = B->rows[start_row+(k-i)+j];

	wi_t ii;
        for(ii = 0; ii + 8 <= B->width - 1; ii += 8) {
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

  word mask_end = __M4RI_LEFT_BITMASK(B->ncols % m4ri_radix);

  if(k == 0) {
    /* __M4RI_CPU_L2_CACHE == __M4RI_TRSM_NTABLES * 2^k * B->width * 8 */
    k = (int)log2((__M4RI_CPU_L2_CACHE/8)/(double)B->width/(double)__M4RI_TRSM_NTABLES);

    rci_t const klog = round(0.75 * log2_floor(MIN(B->nrows, B->ncols)));

    if(klog < k) k = klog;
    if (k<2)     k = 2;
    else if(k>8) k = 8;
  }


  int kk = __M4RI_TRSM_NTABLES * k;
  assert(kk <= m4ri_radix);

  mzd_t *T[__M4RI_TRSM_NTABLES];
  rci_t *L[__M4RI_TRSM_NTABLES];

#ifdef __M4RI_HAVE_SSE2
  mzd_t *Talign[__M4RI_TRSM_NTABLES];
  int b_align = (__M4RI_ALIGNMENT(B->rows[0], 16) == 8);
#endif

  for(int i=0; i<__M4RI_TRSM_NTABLES; i++) {
#ifdef __M4RI_HAVE_SSE2
    /* we make sure that T are aligned as C */
    Talign[i] = mzd_init(__M4RI_TWOPOW(k), B->ncols + m4ri_radix);
    T[i] = mzd_init_window(Talign[i], 0, b_align*m4ri_radix, Talign[i]->nrows, B->ncols + b_align*m4ri_radix);
#else
    T[i] = mzd_init(__M4RI_TWOPOW(k), B->ncols);
#endif
    L[i] = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  }

  rci_t i = 0;
  for (; i < B->nrows - kk; i += kk) {

    _mzd_trsm_upper_left_submatrix(U, B, B->nrows-i-kk, kk, mask_end);

    switch(__M4RI_TRSM_NTABLES) {
    case 8:  mzd_make_table(B, B->nrows - i - 8*k, 0, k, T[7], L[7]);
    case 7:  mzd_make_table(B, B->nrows - i - 7*k, 0, k, T[6], L[6]);
    case 6:  mzd_make_table(B, B->nrows - i - 6*k, 0, k, T[5], L[5]);
    case 5:  mzd_make_table(B, B->nrows - i - 5*k, 0, k, T[4], L[4]);
    case 4:  mzd_make_table(B, B->nrows - i - 4*k, 0, k, T[3], L[3]);
    case 3:  mzd_make_table(B, B->nrows - i - 3*k, 0, k, T[2], L[2]);
    case 2:  mzd_make_table(B, B->nrows - i - 2*k, 0, k, T[1], L[1]);
    case 1:  mzd_make_table(B, B->nrows - i - 1*k, 0, k, T[0], L[0]);
      break;
    default:
      m4ri_die("__M4RI_TRSM_NTABLES must be <= 8 but got %d", __M4RI_TRSM_NTABLES);
    }


    for(rci_t j = 0; j < B->nrows - i - kk; ++j) {
      rci_t x;
      const word *t[__M4RI_TRSM_NTABLES];

      switch(__M4RI_TRSM_NTABLES) {
      case 8: x = L[7][ mzd_read_bits_int(U, j, B->nrows - i - 8*k, k) ]; t[7] = T[7]->rows[x];
      case 7: x = L[6][ mzd_read_bits_int(U, j, B->nrows - i - 7*k, k) ]; t[6] = T[6]->rows[x];
      case 6: x = L[5][ mzd_read_bits_int(U, j, B->nrows - i - 6*k, k) ]; t[5] = T[5]->rows[x];
      case 5: x = L[4][ mzd_read_bits_int(U, j, B->nrows - i - 5*k, k) ]; t[4] = T[4]->rows[x];
      case 4: x = L[3][ mzd_read_bits_int(U, j, B->nrows - i - 4*k, k) ]; t[3] = T[3]->rows[x];
      case 3: x = L[2][ mzd_read_bits_int(U, j, B->nrows - i - 3*k, k) ]; t[2] = T[2]->rows[x];
      case 2: x = L[1][ mzd_read_bits_int(U, j, B->nrows - i - 2*k, k) ]; t[1] = T[1]->rows[x];
      case 1: x = L[0][ mzd_read_bits_int(U, j, B->nrows - i - 1*k, k) ]; t[0] = T[0]->rows[x];
        break;
      default:
        m4ri_die("__M4RI_TRSM_NTABLES must be <= 8 but got %d", __M4RI_TRSM_NTABLES);
      }

      word *b = B->rows[j];
      switch(__M4RI_TRSM_NTABLES) {
      case 8: _mzd_combine_8(b, t, wide); break;
      case 7: _mzd_combine_7(b, t, wide); break;
      case 6: _mzd_combine_6(b, t, wide); break;
      case 5: _mzd_combine_5(b, t, wide); break;
      case 4: _mzd_combine_4(b, t, wide); break;
      case 3: _mzd_combine_3(b, t, wide); break;
      case 2: _mzd_combine_2(b, t, wide); break;
      case 1: _mzd_combine(b, t[0], wide);
        break;
      default:
        m4ri_die("__M4RI_TRSM_NTABLES must be <= 8 but got %d", __M4RI_TRSM_NTABLES);
      }
    }
  }

  /* handle stuff that doesn't fit in multiples of kk */
  for ( ;i < B->nrows; i += k) {
    if (i > B->nrows - k)
      k = B->nrows - i;

    _mzd_trsm_upper_left_submatrix(U, B, B->nrows-i-k, k, mask_end);

    mzd_make_table(B, B->nrows - i - 1*k, 0, k, T[0], L[0]);

    for(rci_t j = 0; j < B->nrows - i - k; ++j) {
      rci_t const x0 = L[0][ mzd_read_bits_int(U, j, B->nrows - i - 1*k, k) ];

      word *b = B->rows[j];
      word *t0 = T[0]->rows[x0];

      for (wi_t ii = 0; ii < wide; ++ii)
        b[ii] ^= t0[ii];
    }
  }
  for(int i=0; i<__M4RI_TRSM_NTABLES; i++) {
    mzd_free(T[i]);
#ifdef __M4RI_HAVE_SSE2
    mzd_free(Talign[i]);
#endif
    m4ri_mm_free(L[i]);
  }

  __M4RI_DD_MZD(B);
}

void _mzd_trsm_lower_left_submatrix(mzd_t const *L, mzd_t *B, rci_t const start_row, int const k, word const mask_end) {
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < i; ++j) {
      if (mzd_read_bit(L, start_row+i, start_row+j)) {
        word *a = B->rows[start_row+i];
        word *b = B->rows[start_row+j];

	wi_t ii;
        for(ii = 0; ii + 8 <= B->width - 1; ii += 8) {
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

void _mzd_trsm_lower_left_russian(mzd_t const *L, mzd_t *B, int k) {
  wi_t const wide = B->width;

  if(k == 0) {
    /* __M4RI_CPU_L2_CACHE == __M4RI_TRSM_NTABLES * 2^k * B->width * 8 */
    k = (int)log2((__M4RI_CPU_L2_CACHE/8)/(double)B->width/(double)__M4RI_TRSM_NTABLES);

    rci_t const klog = round(0.75 * log2_floor(MIN(B->nrows, B->ncols)));

    if(klog < k) k = klog;
    if (k<2)     k = 2;
    else if(k>8) k = 8;
  }
  int kk = __M4RI_TRSM_NTABLES * k;
  assert(kk <= m4ri_radix);

  mzd_t *T[__M4RI_TRSM_NTABLES];
  rci_t *J[__M4RI_TRSM_NTABLES];

#ifdef __M4RI_HAVE_SSE2
    /* we make sure that T are aligned as B, this is dirty, we need a function for this */
  mzd_t *Talign[__M4RI_TRSM_NTABLES];
  int b_align = (__M4RI_ALIGNMENT(B->rows[0], 16) == 8);
#endif

  for(int i=0; i<__M4RI_TRSM_NTABLES; i++) {
#ifdef __M4RI_HAVE_SSE2
    Talign[i] = mzd_init(__M4RI_TWOPOW(k), B->ncols + m4ri_radix);
    T[i] = mzd_init_window(Talign[i], 0, b_align*m4ri_radix, Talign[i]->nrows, B->ncols + b_align*m4ri_radix);
#else
    T[i] = mzd_init(__M4RI_TWOPOW(k), B->ncols);
#endif
    J[i] = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  }

  const word mask = __M4RI_LEFT_BITMASK(k);
  rci_t i = 0;
  for (; i < B->nrows - kk; i += kk) {

    _mzd_trsm_lower_left_submatrix(L, B, i, kk, B->high_bitmask);

    switch(__M4RI_TRSM_NTABLES) {
    case 8:  mzd_make_table(B, i + 7*k, 0, k, T[7], J[7]);
    case 7:  mzd_make_table(B, i + 6*k, 0, k, T[6], J[6]);
    case 6:  mzd_make_table(B, i + 5*k, 0, k, T[5], J[5]);
    case 5:  mzd_make_table(B, i + 4*k, 0, k, T[4], J[4]);
    case 4:  mzd_make_table(B, i + 3*k, 0, k, T[3], J[3]);
    case 3:  mzd_make_table(B, i + 2*k, 0, k, T[2], J[2]);
    case 2:  mzd_make_table(B, i + 1*k, 0, k, T[1], J[1]);
    case 1:  mzd_make_table(B, i + 0*k, 0, k, T[0], J[0]);
      break;
    default:
      m4ri_die("__M4RI_TRSM_NTABLES must be <= 8 but got %d", __M4RI_TRSM_NTABLES);
    }


    for(rci_t j = i+kk; j < B->nrows; ++j) {
      const word *t[__M4RI_TRSM_NTABLES];

      word tmp = mzd_read_bits(L, j, i, kk);

      switch(__M4RI_TRSM_NTABLES) {
      case 8: t[7] = T[7]->rows[ J[7][ (tmp >> (7*k)) & mask ] ];
      case 7: t[6] = T[6]->rows[ J[6][ (tmp >> (6*k)) & mask ] ];
      case 6: t[5] = T[5]->rows[ J[5][ (tmp >> (5*k)) & mask ] ];
      case 5: t[4] = T[4]->rows[ J[4][ (tmp >> (4*k)) & mask ] ];
      case 4: t[3] = T[3]->rows[ J[3][ (tmp >> (3*k)) & mask ] ];
      case 3: t[2] = T[2]->rows[ J[2][ (tmp >> (2*k)) & mask ] ];
      case 2: t[1] = T[1]->rows[ J[1][ (tmp >> (1*k)) & mask ] ];
      case 1: t[0] = T[0]->rows[ J[0][ (tmp >> (0*k)) & mask ] ];
        break;
      default:
        m4ri_die("__M4RI_TRSM_NTABLES must be <= 8 but got %d", __M4RI_TRSM_NTABLES);
      }

      word *b = B->rows[j];
      switch(__M4RI_TRSM_NTABLES) {
      case 8: _mzd_combine_8(b, t, wide); break;
      case 7: _mzd_combine_7(b, t, wide); break;
      case 6: _mzd_combine_6(b, t, wide); break;
      case 5: _mzd_combine_5(b, t, wide); break;
      case 4: _mzd_combine_4(b, t, wide); break;
      case 3: _mzd_combine_3(b, t, wide); break;
      case 2: _mzd_combine_2(b, t, wide); break;
      case 1: _mzd_combine(b, t[0], wide);
        break;
      default:
        m4ri_die("__M4RI_TRSM_NTABLES must be <= 8 but got %d", __M4RI_TRSM_NTABLES);
      }
    }
  }

  /* handle stuff that doesn't fit in multiples of kk */
  for ( ;i < B->nrows; i += k) {
    if (i > B->nrows - k)
      k = B->nrows - i;

    _mzd_trsm_lower_left_submatrix(L, B, i, k, B->high_bitmask);

    mzd_make_table(B, i + 0*k, 0, k, T[0], J[0]);

    for(rci_t j = i+k; j < L->nrows; ++j) {
      rci_t const x0 = J[0][ mzd_read_bits_int(L, j, i, k) ];

      word *b = B->rows[j];
      word *t0 = T[0]->rows[x0];

      for (wi_t ii = 0; ii < wide; ++ii)
        b[ii] ^= t0[ii];
    }
  }
  for(int i=0; i<__M4RI_TRSM_NTABLES; i++) {
    mzd_free(T[i]);
#ifdef __M4RI_HAVE_SSE2
    mzd_free(Talign[i]);
#endif
    m4ri_mm_free(J[i]);
  }

  __M4RI_DD_MZD(B);
}


void mzd_make_table_trtri(mzd_t const *M, rci_t r, rci_t c, int k, ple_table_t *Tb, rci_t startcol) {
  mzd_t *T = Tb->T;
  rci_t *L = Tb->E;

  assert(!(T->flags & mzd_flag_multiple_blocks));
  wi_t const blockoffset  = c / m4ri_radix;
  wi_t const blockoffset0 = startcol / m4ri_radix;

  assert(blockoffset - blockoffset0 <= 1);

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
    T->rows[i][blockoffset0] = 0; /* we make sure that we can safely add from blockoffset0 */
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
  Tb->B[0] = 0;
  for(int i=1; i<twokay; ++i) {
    mzd_xor_bits(T, i, c, k, (word)m4ri_codebook[k]->ord[i]);
    Tb->B[i] = mzd_read_bits(T, i, startcol, m4ri_radix);
  }
}

#define __M4RI_TRTRI_NTABLES 4

static inline void _mzd_trtri_upper_submatrix(mzd_t *A, rci_t pivot_r, rci_t elim_r, const int k) {
  for(rci_t i=pivot_r; i<pivot_r+k; i++)
    for(rci_t j=elim_r; j<i; j++)
      if(mzd_read_bit(A,j,i) && (i+1)<A->ncols )
        mzd_row_add_offset(A, j, i, i+1);
}


mzd_t *mzd_trtri_upper_russian(mzd_t *A, int k) {
  assert(A->nrows == A->ncols);

  if (k == 0) {
    k = m4ri_opt_k(A->nrows, A->ncols, 0);
    if (k >= 7)
      k = 7;
    if (0.75 * __M4RI_TWOPOW(k) *A->ncols > __M4RI_CPU_L3_CACHE / 2.0)
      k -= 1;
  }

  const int kk = __M4RI_TRTRI_NTABLES*k;

  int k_[__M4RI_TRTRI_NTABLES];
  for (int i=0; i<__M4RI_TRTRI_NTABLES; i++)
    k_[i] = k;

  ple_table_t *T[__M4RI_TRTRI_NTABLES];
  mzd_t *U[__M4RI_TRTRI_NTABLES];
  for(int i=0; i<__M4RI_TRTRI_NTABLES; i++) {
    T[i] = ple_table_init(k, A->ncols);
    U[i] = mzd_init(k, A->ncols);
  }

  /** dummy offsets table for _mzd_ple_to_e**/
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
    mzd_make_table_trtri(U[0], 0, r,   k, T[0], r);

    _mzd_trtri_upper_submatrix(A, r+k, r, k);
    _mzd_ple_to_e(U[1], A, r+k, r+k, k, id);
    mzd_make_table_trtri(U[1], 0, r+k, k, T[1], r);

    _mzd_trtri_upper_submatrix(A, r+2*k, r, k);
    _mzd_ple_to_e(U[2], A, r+2*k, r+2*k, k, id);
    mzd_make_table_trtri(U[2], 0, r+2*k, k, T[2], r);

    _mzd_trtri_upper_submatrix(A, r+3*k, r, k);
    _mzd_ple_to_e(U[3], A, r+3*k, r+3*k, k, id);
    mzd_make_table_trtri(U[3], 0, r+3*k, k, T[3], r);

    _mzd_process_rows_ple_4(A, 0, r, r, k_, (const ple_table_t** const)T);
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
    mzd_make_table_trtri(U[0], 0, r, k, T[0], r);

    mzd_process_rows(A, 0, r, r, k, T[0]->T, T[0]->E);
    r += k;
  }

  for(int i=0; i<__M4RI_TRTRI_NTABLES; i++) {
    ple_table_free(T[i]);
    mzd_free(U[i]);
  }
  __M4RI_DD_MZD(A);
  return A;
}
