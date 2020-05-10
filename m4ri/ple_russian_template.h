#include <m4ri/misc.h>

void __M4RI_TEMPLATE_NAME(_mzd_process_rows_ple)(mzd_t *M, rci_t startrow, rci_t stoprow,
                                                 rci_t startcol, int const k[N],
                                                 const ple_table_t *table[N]) {
  assert(1 <= N && N <= 8);

  const mzd_t *T[N];
  const rci_t *E[N];
  const word *B[N];
  word bm[N];
  int sh[N];
  int x[N];
  const word *t[N];

  switch (N) { /* we rely on the compiler to optimise this switch away, it reads nicer than #if */
  case 8:
    T[7]  = table[7]->T;
    E[7]  = table[7]->E;
    B[7]  = table[7]->B;
    bm[7] = __M4RI_LEFT_BITMASK(k[7]);
    sh[7] = k[0] + k[1] + k[2] + k[3] + k[4] + k[5] + k[6];
  case 7:
    T[6]  = table[6]->T;
    E[6]  = table[6]->E;
    B[6]  = table[6]->B;
    bm[6] = __M4RI_LEFT_BITMASK(k[6]);
    sh[6] = k[0] + k[1] + k[2] + k[3] + k[4] + k[5];
  case 6:
    T[5]  = table[5]->T;
    E[5]  = table[5]->E;
    B[5]  = table[5]->B;
    bm[5] = __M4RI_LEFT_BITMASK(k[5]);
    sh[5] = k[0] + k[1] + k[2] + k[3] + k[4];
  case 5:
    T[4]  = table[4]->T;
    E[4]  = table[4]->E;
    B[4]  = table[4]->B;
    bm[4] = __M4RI_LEFT_BITMASK(k[4]);
    sh[4] = k[0] + k[1] + k[2] + k[3];
  case 4:
    T[3]  = table[3]->T;
    E[3]  = table[3]->E;
    B[3]  = table[3]->B;
    bm[3] = __M4RI_LEFT_BITMASK(k[3]);
    sh[3] = k[0] + k[1] + k[2];
  case 3:
    T[2]  = table[2]->T;
    E[2]  = table[2]->E;
    B[2]  = table[2]->B;
    bm[2] = __M4RI_LEFT_BITMASK(k[2]);
    sh[2] = k[0] + k[1];
  case 2:
    T[1]  = table[1]->T;
    E[1]  = table[1]->E;
    B[1]  = table[1]->B;
    bm[1] = __M4RI_LEFT_BITMASK(k[1]);
    sh[1] = k[0];
  case 1:
    T[0]  = table[0]->T;
    E[0]  = table[0]->E;
    B[0]  = table[0]->B;
    bm[0] = __M4RI_LEFT_BITMASK(k[0]);
    sh[0] = 0;
  }

  wi_t const block = startcol / m4ri_radix;
  wi_t const wide  = M->width - block;

  for (rci_t r = startrow; r < stoprow; ++r) {
    word bits = mzd_read_bits(M, r, startcol, sh[N - 1] + k[N - 1]);
    word *m   = M->rows[r] + block;

    switch (N) { /* we rely on the compiler to optimise this switch away, it reads nicer than #if */
    case 8:
      x[N - 8] = E[N - 8][(bits >> sh[N - 8]) & bm[N - 8]];
      bits ^= B[N - 8][x[N - 8]];
      t[N - 8] = T[N - 8]->rows[x[N - 8]] + block;
    case 7:
      x[N - 7] = E[N - 7][(bits >> sh[N - 7]) & bm[N - 7]];
      bits ^= B[N - 7][x[N - 7]];
      t[N - 7] = T[N - 7]->rows[x[N - 7]] + block;
    case 6:
      x[N - 6] = E[N - 6][(bits >> sh[N - 6]) & bm[N - 6]];
      bits ^= B[N - 6][x[N - 6]];
      t[N - 6] = T[N - 6]->rows[x[N - 6]] + block;
    case 5:
      x[N - 5] = E[N - 5][(bits >> sh[N - 5]) & bm[N - 5]];
      bits ^= B[N - 5][x[N - 5]];
      t[N - 5] = T[N - 5]->rows[x[N - 5]] + block;
    case 4:
      x[N - 4] = E[N - 4][(bits >> sh[N - 4]) & bm[N - 4]];
      bits ^= B[N - 4][x[N - 4]];
      t[N - 4] = T[N - 4]->rows[x[N - 4]] + block;
    case 3:
      x[N - 3] = E[N - 3][(bits >> sh[N - 3]) & bm[N - 3]];
      bits ^= B[N - 3][x[N - 3]];
      t[N - 3] = T[N - 3]->rows[x[N - 3]] + block;
    case 2:
      x[N - 2] = E[N - 2][(bits >> sh[N - 2]) & bm[N - 2]];
      bits ^= B[N - 2][x[N - 2]];
      t[N - 2] = T[N - 2]->rows[x[N - 2]] + block;
    case 1:
      x[N - 1] = E[N - 1][(bits >> sh[N - 1]) & bm[N - 1]];
      bits ^= B[N - 1][x[N - 1]];
      t[N - 1] = T[N - 1]->rows[x[N - 1]] + block;
    }

    __M4RI_TEMPLATE_NAME(_mzd_combine)(m, t, wide);
  }

  __M4RI_DD_MZD(M);
}

void __M4RI_TEMPLATE_NAME(_mzd_ple_a11)(mzd_t *A, rci_t const start_row, rci_t const stop_row,
                                        rci_t const start_col, wi_t const block, int const k[N],
                                        ple_table_t const *table[N]) {

  wi_t const wide = A->width - block;

  if (wide <= 0) return;

  const mzd_t *T[N];
  const rci_t *M[N];
  word bm[N];
  int sh[N];
  int x[N];
  const word *t[N];

  switch (N) { /* we rely on the compiler to optimise this switch away, it reads nicer than #if */
  case 8:
    T[7]  = table[7]->T;
    M[7]  = table[7]->M;
    bm[7] = __M4RI_LEFT_BITMASK(k[7]);
    sh[7] = k[0] + k[1] + k[2] + k[3] + k[4] + k[5] + k[6];
  case 7:
    T[6]  = table[6]->T;
    M[6]  = table[6]->M;
    bm[6] = __M4RI_LEFT_BITMASK(k[6]);
    sh[6] = k[0] + k[1] + k[2] + k[3] + k[4] + k[5];
  case 6:
    T[5]  = table[5]->T;
    M[5]  = table[5]->M;
    bm[5] = __M4RI_LEFT_BITMASK(k[5]);
    sh[5] = k[0] + k[1] + k[2] + k[3] + k[4];
  case 5:
    T[4]  = table[4]->T;
    M[4]  = table[4]->M;
    bm[4] = __M4RI_LEFT_BITMASK(k[4]);
    sh[4] = k[0] + k[1] + k[2] + k[3];
  case 4:
    T[3]  = table[3]->T;
    M[3]  = table[3]->M;
    bm[3] = __M4RI_LEFT_BITMASK(k[3]);
    sh[3] = k[0] + k[1] + k[2];
  case 3:
    T[2]  = table[2]->T;
    M[2]  = table[2]->M;
    bm[2] = __M4RI_LEFT_BITMASK(k[2]);
    sh[2] = k[0] + k[1];
  case 2:
    T[1]  = table[1]->T;
    M[1]  = table[1]->M;
    bm[1] = __M4RI_LEFT_BITMASK(k[1]);
    sh[1] = k[0];
  case 1:
    T[0]  = table[0]->T;
    M[0]  = table[0]->M;
    bm[0] = __M4RI_LEFT_BITMASK(k[0]);
    sh[0] = 0;
  };

  const rci_t bits_to_read = sh[N - 1] + k[N - 1];

  for (rci_t i = start_row; i < stop_row; ++i) {
    const word bits = mzd_read_bits(A, i, start_col, bits_to_read);
    word *m         = A->rows[i] + block;

    switch (N) { /* we rely on the compiler to optimise this switch away, it reads nicer than #if */
    case 8:
      x[N - 8] = M[N - 8][(bits >> sh[N - 8]) & bm[N - 8]];
      t[N - 8] = T[N - 8]->rows[x[N - 8]] + block;
    case 7:
      x[N - 7] = M[N - 7][(bits >> sh[N - 7]) & bm[N - 7]];
      t[N - 7] = T[N - 7]->rows[x[N - 7]] + block;
    case 6:
      x[N - 6] = M[N - 6][(bits >> sh[N - 6]) & bm[N - 6]];
      t[N - 6] = T[N - 6]->rows[x[N - 6]] + block;
    case 5:
      x[N - 5] = M[N - 5][(bits >> sh[N - 5]) & bm[N - 5]];
      t[N - 5] = T[N - 5]->rows[x[N - 5]] + block;
    case 4:
      x[N - 4] = M[N - 4][(bits >> sh[N - 4]) & bm[N - 4]];
      t[N - 4] = T[N - 4]->rows[x[N - 4]] + block;
    case 3:
      x[N - 3] = M[N - 3][(bits >> sh[N - 3]) & bm[N - 3]];
      t[N - 3] = T[N - 3]->rows[x[N - 3]] + block;
    case 2:
      x[N - 2] = M[N - 2][(bits >> sh[N - 2]) & bm[N - 2]];
      t[N - 2] = T[N - 2]->rows[x[N - 2]] + block;
    case 1:
      x[N - 1] = M[N - 1][(bits >> sh[N - 1]) & bm[N - 1]];
      t[N - 1] = T[N - 1]->rows[x[N - 1]] + block;
    }
    __M4RI_TEMPLATE_NAME(_mzd_combine)(m, t, wide);
  }
  __M4RI_DD_MZD(A);
}
