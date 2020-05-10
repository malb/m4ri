/*
 * bench_packedmatrix.c
 *
 * Application to test functionality of packedmatrix.c.
 *
 * Copyright (C) 2011  Carlo Wood  <carlo@alinoe.com>
 * RSA-1024 0x624ACAD5 1997-01-26                    Sign & Encrypt
 * Fingerprint16 = 32 EC A7 B6 AC DB 65 A6  F6 F6 55 DD 1C DC FF 61
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS 1
#endif

#include <m4ri/config.h>

#ifdef HAVE_LIBPAPI
#define _GNU_SOURCE
#include <errno.h>
#include <papi.h>
#include <sys/types.h>  // papi.h needs caddr_t
#endif

#include <ctype.h>
#include <inttypes.h>
#include <stdlib.h>

#include "benchmarking.h"
#include "cpucycles.h"
#include <m4ri/m4ri.h>

struct test_params {
  rci_t m;
  rci_t n;
  rci_t k;
  rci_t l;
  rci_t row[3];
  int rows;
  rci_t col[3];
  int cols;
  wi_t wrd[3];
  int wrds;
  uint64_t count;
  int cutoff;
  int boolean;
  int integer;
  char const *funcname;
};

typedef int (*run_type)(void *, unsigned long long *, int *);

static unsigned long long loop_calibration[32];

#ifdef HAVE_LIBPAPI

#define BENCHMARK_PREFIX(mzd_func)                                                                 \
  int run_##mzd_func(void *_p, unsigned long long *data, int *data_len) {                          \
    *data_len             = MIN(papi_array_len + 1, *data_len);                                    \
    struct test_params *p = (struct test_params *)_p;                                              \
    int papi_res;                                                                                  \
    do
#define TIME_BEGIN(mzd_func_with_ARGS)                                                             \
  do {                                                                                             \
    int array_len = *data_len - 1;                                                                 \
    mzd_func_with_ARGS;                                                                            \
    unsigned long long t0 = PAPI_get_virt_usec();                                                  \
  papi_res = PAPI_start_counters((int *)papi_events, array_len)
#define TIME_END                                                                                   \
  PAPI_stop_counters((long long *)&data[1], array_len);                                            \
  t0      = PAPI_get_virt_usec() - t0;                                                             \
  data[0] = t0;                                                                                    \
  for (int nv = 0; nv <= array_len; ++nv) {                                                        \
    if (data[nv] < loop_calibration[nv]) loop_calibration[nv] = data[nv];                          \
    data[nv] -= loop_calibration[nv];                                                              \
  }                                                                                                \
  }                                                                                                \
  while (0)
#define BENCHMARK_POSTFIX                                                                          \
  while (0)                                                                                        \
    ;                                                                                              \
  return papi_res;                                                                                 \
  }

#else  // HAVE_LIBPAPI

#define BENCHMARK_PREFIX(mzd_func)                                                                 \
  int run_##mzd_func(void *_p, unsigned long long *data, int *data_len) {                          \
    *data_len             = 2;                                                                     \
    struct test_params *p = (struct test_params *)_p;                                              \
    do
#define TIME_BEGIN(mzd_func_with_ARGS)                                                             \
  do {                                                                                             \
    mzd_func_with_ARGS;                                                                            \
    data[0] = walltime(0);                                                                         \
  data[1] = cpucycles()
#define TIME_END                                                                                   \
  data[1] = cpucycles() - data[1];                                                                 \
  data[0] = walltime(data[0]);                                                                     \
  }                                                                                                \
  while (0)
#define BENCHMARK_POSTFIX                                                                          \
  while (0)                                                                                        \
    ;                                                                                              \
  return 0;                                                                                        \
  }

#endif  // HAVE_LIBPAPI

#define TIME(mzd_func_with_ARGS)                                                                   \
  TIME_BEGIN(mzd_func_with_ARGS);                                                                  \
  for (uint64_t i = 0; i < loop_count; ++i) { mzd_func_with_ARGS; }                                \
  TIME_END

mzd_t *volatile vA;
rci_t volatile vrowa;
rci_t volatile vcola;
rci_t volatile vrowb;
rci_t volatile vcolb;
wi_t volatile vstartblock;
int volatile vn;
int volatile vint;
word volatile vword;
BIT volatile vbit;

BENCHMARK_PREFIX(bench_nothing) {
  mzd_t *const A = mzd_init(64, 64);
  mzd_randomize(A);
  uint64_t volatile loop_count = p->count;

  TIME();

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(_mzd_row_swap) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);

  vA          = A;
  vrowa       = p->row[0];
  vrowb       = p->row[1];
  vstartblock = p->wrd[0];

  uint64_t const loop_count = p->count;

  TIME(_mzd_row_swap(vA, vrowa, vrowb, vstartblock));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_row_swap) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const rowa          = p->row[0];
  rci_t const rowb          = p->row[1];
  uint64_t const loop_count = p->count;

  TIME(mzd_row_swap(A, rowa, rowb));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_copy_row) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_t *const B = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const rowa          = p->row[0];
  rci_t const rowb          = p->row[1];
  uint64_t const loop_count = p->count;

  TIME(mzd_copy_row(B, rowb, A, rowa));

  mzd_free(A);
  mzd_free(B);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_col_swap) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const cola          = p->col[0];
  rci_t const colb          = p->col[1];
  uint64_t const loop_count = p->count;

  TIME(mzd_col_swap(A, cola, colb));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_col_swap_in_rows) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);

  vA    = A;
  vcola = p->col[0];
  vcolb = p->col[1];
  vrowa = p->row[0];
  vrowb = p->row[1];

  uint64_t const loop_count = p->count;

  TIME(mzd_col_swap_in_rows(vA, vcola, vcolb, vrowa, vrowb));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_read_bit) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);

  vA    = A;
  vrowa = p->row[0];
  vcola = p->col[0];

  uint64_t const loop_count = p->count;

  TIME(vbit = mzd_read_bit(vA, vrowa, vcola));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_write_bit) {
  mzd_t *const A = mzd_init(p->m, p->n);

  vA    = A;
  vrowa = p->row[0];
  vcola = p->col[0];
  vbit  = 0;

  uint64_t const loop_count = p->count;

  TIME(mzd_write_bit(vA, vrowa, vcola, vbit); vbit = !vbit);

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_row_add_offset) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);

  vA    = A;
  vrowa = p->row[0];
  vrowb = p->row[1];
  vcola = p->col[0];

  uint64_t const loop_count = p->count;

  TIME(mzd_row_add_offset(vA, vrowa, vrowb, vcola));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_row_add) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const rowa          = p->row[0];
  rci_t const rowb          = p->row[1];
  uint64_t const loop_count = p->count;

  TIME(mzd_row_add(A, rowa, rowb));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_transpose) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_t *const B = mzd_init(p->n, p->m);
  mzd_randomize(A);
  uint64_t const loop_count = p->count;

  TIME(mzd_transpose(B, A));

  mzd_free(A);
  mzd_free(B);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_mul_naive) {
  mzd_t *const A = mzd_init(p->m, p->l);
  mzd_t *const B = mzd_init(p->l, p->n);
  mzd_t *const C = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_randomize(B);
  uint64_t const loop_count = p->count;

  TIME(mzd_mul_naive(C, A, B));

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_addmul_naive) {
  mzd_t *const A = mzd_init(p->m, p->l);
  mzd_t *const B = mzd_init(p->l, p->n);
  mzd_t *const C = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_randomize(B);
  uint64_t const loop_count = p->count;

  TIME(mzd_addmul_naive(C, A, B));

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(_mzd_mul_naive) {
  mzd_t *const A = mzd_init(p->m, p->l);
  mzd_t *const B = mzd_init(p->n, p->l);
  mzd_t *const C = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_randomize(B);
  int const clear           = p->boolean;
  uint64_t const loop_count = p->count;

  TIME(_mzd_mul_naive(C, A, B, clear));

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(_mzd_mul_va) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_t *const V = mzd_init(1, p->m);
  mzd_t *const C = mzd_init(1, p->n);
  mzd_randomize(A);
  mzd_randomize(V);
  int const clear           = p->boolean;
  uint64_t const loop_count = p->count;

  TIME(_mzd_mul_va(C, V, A, clear));

  mzd_free(A);
  mzd_free(V);
  mzd_free(C);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_gauss_delayed) {
  mzd_t **A                 = malloc(sizeof(mzd_t) * (p->count + 1));
  rci_t const cola          = p->col[0];
  int const full            = p->boolean;
  uint64_t const loop_count = p->count;
  rci_t result;

  for (int i = loop_count; i >= 0; --i) {
    A[i] = mzd_init(p->m, p->n);
    mzd_randomize(A[i]);
  }

  TIME_BEGIN(result = mzd_gauss_delayed(A[0], cola, full));
  for (int i = loop_count; i > 0; --i) { result = mzd_gauss_delayed(A[i], cola, full); }
  TIME_END;

  for (int i = 0; i <= loop_count; ++i) mzd_free(A[i]);
  free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_echelonize_naive) {
  mzd_t **A                 = malloc(sizeof(mzd_t) * (p->count + 1));
  int const full            = p->boolean;
  uint64_t const loop_count = p->count;
  rci_t result;

  for (int i = loop_count; i >= 0; --i) {
    A[i] = mzd_init(p->m, p->n);
    mzd_randomize(A[i]);
  }

  TIME_BEGIN(result = mzd_echelonize_naive(A[0], full));
  for (int i = loop_count; i > 0; --i) { result = mzd_echelonize_naive(A[i], full); }
  TIME_END;

  for (int i = 0; i <= loop_count; ++i) mzd_free(A[i]);
  free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_equal) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_t *const B            = mzd_copy(NULL, A);
  uint64_t const loop_count = p->count;
  int volatile result;

  TIME(result = mzd_equal(A, B));

  mzd_free(A);
  mzd_free(B);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_cmp) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_t *const B            = mzd_copy(NULL, A);
  uint64_t const loop_count = p->count;
  int volatile result;

  TIME(result = mzd_cmp(A, B));

  mzd_free(A);
  mzd_free(B);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_copy) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_t *const B = mzd_init(p->m, p->n);
  mzd_randomize(A);
  uint64_t const loop_count = p->count;

  TIME(mzd_copy(B, A));

  mzd_free(A);
  mzd_free(B);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_concat) {
  mzd_t *const A = mzd_init(p->m, p->k);
  mzd_t *const B = mzd_init(p->m, p->l);
  mzd_t *const C = mzd_init(p->m, p->k + p->l);
  mzd_randomize(A);
  mzd_randomize(B);
  uint64_t const loop_count = p->count;

  TIME(mzd_concat(C, A, B));

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_stack) {
  mzd_t *const A = mzd_init(p->k, p->n);
  mzd_t *const B = mzd_init(p->l, p->n);
  mzd_t *const C = mzd_init(p->k + p->l, p->n);
  mzd_randomize(A);
  mzd_randomize(B);
  uint64_t const loop_count = p->count;

  TIME(mzd_stack(C, A, B));

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_submatrix) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const rowa          = p->row[0];
  rci_t const cola          = p->col[0];
  rci_t const rowb          = p->row[1];
  rci_t const colb          = p->col[1];
  mzd_t *const S            = mzd_init(rowb - rowa, colb - cola);
  uint64_t const loop_count = p->count;

  TIME(mzd_submatrix(S, A, rowa, cola, rowb, colb));

  mzd_free(A);
  mzd_free(S);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_invert_naive) {
  mzd_t *const A = mzd_init(p->m, p->m);
  mzd_t *const I = mzd_init(p->m, p->m);
  mzd_t *const C = mzd_init(p->m, p->m);
  mzd_randomize(A);
  mzd_set_ui(I, 1);
  uint64_t const loop_count = p->count;

  TIME(mzd_invert_naive(C, A, I));

  mzd_free(A);
  mzd_free(I);
  mzd_free(C);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_add) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_t *const B = mzd_init(p->m, p->n);
  mzd_t *const C = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_randomize(B);
  uint64_t const loop_count = p->count;

  TIME(mzd_add(C, A, B));

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(_mzd_add) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_t *const B = mzd_init(p->m, p->n);
  mzd_t *const C = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_randomize(B);
  uint64_t const loop_count = p->count;

  TIME(_mzd_add(C, A, B));

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_combine) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_t *const B = mzd_init(p->m, p->n);
  mzd_t *const C = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_randomize(B);
  rci_t row1                = p->row[0];
  rci_t row2                = p->row[1];
  rci_t row3                = p->row[2];
  wi_t startblock           = p->wrd[0];
  uint64_t const loop_count = p->count;

  TIME(mzd_combine(C, row3, startblock, A, row1, startblock, B, row2, startblock));

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_read_bits) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);

  vA    = A;
  vrowa = p->row[0];
  vcola = p->col[0];
  vn    = p->integer;

  uint64_t const loop_count = p->count;

  TIME(vword = mzd_read_bits(vA, vrowa, vcola, vn));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_read_bits_int) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);

  vA    = A;
  vrowa = p->row[0];
  vcola = p->col[0];
  vn    = p->integer;

  uint64_t const loop_count = p->count;

  TIME(vint = mzd_read_bits_int(vA, vrowa, vcola, vn));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_xor_bits) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);

  vA    = A;
  vrowa = p->row[0];
  vcola = p->col[0];
  vn    = p->integer;
  vword = 0xffffffffffffffffULL;

  uint64_t const loop_count = p->count;

  TIME(mzd_xor_bits(vA, vrowa, vcola, vn, vword));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_and_bits) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);

  vA    = A;
  vrowa = p->row[0];
  vcola = p->col[0];
  vn    = p->integer;
  vword = 0xffffffffffffffffULL;

  uint64_t const loop_count = p->count;

  TIME(mzd_and_bits(vA, vrowa, vcola, vn, vword));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_clear_bits) {
  mzd_t *volatile A = mzd_init(p->m, p->n);
  mzd_randomize(A);

  vA    = A;
  vrowa = p->row[0];
  vcola = p->col[0];
  vn    = p->integer;

  uint64_t const loop_count = p->count;

  TIME(mzd_clear_bits(vA, vrowa, vcola, vn));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_is_zero) {
  mzd_t *const A            = mzd_init(p->m, p->n);
  uint64_t const loop_count = p->count;
  int volatile result;

  TIME(result = mzd_is_zero(A));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_find_pivot) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t row                 = p->row[0];
  rci_t col                 = p->col[0];
  uint64_t const loop_count = p->count;
  int volatile result;
  rci_t row_out;
  rci_t col_out;

  TIME(result = mzd_find_pivot(A, row, col, &row_out, &col_out));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_density) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  wi_t res                  = p->wrd[0];
  uint64_t const loop_count = p->count;
  double volatile result;

  TIME(result = mzd_density(A, res));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(_mzd_density) {
  mzd_t *const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t row                 = p->row[0];
  rci_t col                 = p->col[0];
  wi_t res                  = p->wrd[0];
  uint64_t const loop_count = p->count;
  double volatile result;

  TIME(result = _mzd_density(A, res, row, col));

  mzd_free(A);
}
BENCHMARK_POSTFIX

BENCHMARK_PREFIX(mzd_first_zero_row) {
  mzd_t *const A = mzd_init(p->m, p->m);
  mzd_set_ui(A, 1);
  uint64_t const loop_count = p->count;
  rci_t volatile result;

  TIME(result = mzd_first_zero_row(A));

  mzd_free(A);
}
BENCHMARK_POSTFIX

// Returns a number proportional with the ideal number of
// mathematical operations for the given code.
double complexity1(struct test_params *p, char code) {
  switch (code) {
  case 'k': return p->k;  // Linear with size 'k' of a matrix.
  case 'l': return p->l;  // Linear with size 'l' of a matrix.
  case 'm': return p->m;  // Linear with the number of rows of the matrix.
  case 'n': return p->n;  // Linear with the number of columns of the matrix.
  case 'W':
    assert(p->n > m4ri_radix * p->wrd[0]);  // Linear with the number of processed columns.
    return p->n - m4ri_radix * p->wrd[0];
  case 'D':
    assert(p->row[0] < p->row[1]);
    return p->row[1] - p->row[0];  // Linear with the number of rows between start_row and stop_row.
  case 'E':
    assert(p->col[0] < p->col[1]);
    return p->col[1] - p->col[0];  // Linear with the number of cols between start_col and stop_col.
  case 'C':
    assert(p->col[0] < p->n);
    return p->n - p->col[0];  // Linear with the number of columns of column col and beyond.
  }
  return 0.0;
}

char const *complexity1_human(struct test_params *p, char code) {
  switch (code) {
  case 'k': return "k";
  case 'l': return "l";
  case 'm': return "m";
  case 'n': return "n";
  case 'W': return "cols";
  case 'D': return "rows";
  case 'E': return "cols";
  case 'C': return "cols";
  }
  return "UNKNOWN";
}

double complexity(struct test_params *p, char const *cp) {
  double c = 1;
  while (*cp) {
    c *= complexity1(p, *cp);
    ++cp;
  }
  return c;
}

void print_complexity_human(struct test_params *p, char const *cp) {
  int first    = 1;
  char last_cp = 0;
  int power    = 0;
  while (*cp) {
    if (*cp != last_cp) {
      if (power > 1) printf("^%d", power);
      if (!first && isupper(*cp)) printf("*");
      printf("%s", complexity1_human(p, *cp));
      power   = 0;
      last_cp = *cp;
    }
    ++power;
    ++cp;
  }
  if (power > 1) printf("^%d", power);
}

struct function_st {
  char const *funcname;
  run_type run_func;
  char const *input_codes;
  char const *complexity_code;
  uint64_t count;
};

typedef struct function_st function_st;

static function_st function_mapper[] = {
    {"_mzd_row_swap", run__mzd_row_swap, "Rmn,ri,ri,wi", "W", 1000000000},
    {"mzd_row_swap", run_mzd_row_swap, "Rmn,ri,ri", "n", 1000000000},
    {"mzd_copy_row", run_mzd_copy_row, "Omn,ri,R,ri", "n", 1000000000},
    {"mzd_col_swap", run_mzd_col_swap, "Rmn,ci,ci", "m", 10000000},
    {"mzd_col_swap_in_rows", run_mzd_col_swap_in_rows, "Rmn,ci,ci,ri,ri", "D", 10000000},
    {"mzd_read_bit", run_mzd_read_bit, "Rmn,ri,ci", "", 100000000},
    {"mzd_write_bit", run_mzd_write_bit, "Omn,ri,ci", "", 100000000},
    {"mzd_row_add_offset", run_mzd_row_add_offset, "Rmn,ri,ri,ci", "C", 100000000},
    {"mzd_row_add", run_mzd_row_add, "Rmn,ri,ri", "n", 100000000},
    {"mzd_transpose", run_mzd_transpose, "Onm,Rmn", "mn", 10000000},
    {"mzd_mul_naive", run_mzd_mul_naive, "Omn,Rml,Rln", "mnl", 10000000},
    {"mzd_addmul_naive", run_mzd_addmul_naive, "Omn,Rml,Rln", "mnl", 10000000},
    {"_mzd_mul_naive", run__mzd_mul_naive, "Omn,Rml,Rnl,b", "mnl", 10000000},
    {"_mzd_mul_va", run__mzd_mul_va, "O1n,V1m,Amn,b", "mn", 1000000000},
    {"mzd_gauss_delayed", run_mzd_gauss_delayed, "Rmn,ci,b", "mC", 10000000},
    {"mzd_echelonize_naive", run_mzd_echelonize_naive, "Rmn,b", "mn", 10000000},
    {"mzd_equal", run_mzd_equal, "Rmn,Rmn", "mn", 1000000000},
    {"mzd_cmp", run_mzd_cmp, "Rmn,Rmn", "mn", 1000000000},
    {"mzd_copy", run_mzd_copy, "Omn,Rmn", "mn", 1000000000},
    {"mzd_concat", run_mzd_concat, "Omn,Rmk,Rml", "mn", 10000000},
    {"mzd_stack", run_mzd_stack, "Omn,Rkn,Rln", "mn", 1000000000},
    {"mzd_submatrix", run_mzd_submatrix, "O,Rmn,ri,ci,ri,ci", "DE", 10000000},
    {"mzd_invert_naive", run_mzd_invert_naive, "Omm,Rmm,Imm", "mmm", 10000000},
    {"mzd_add", run_mzd_add, "Omn,Rmn,Rmn", "mn", 10000000},
    {"_mzd_add", run__mzd_add, "Omn,Rmn,Rmn", "mn", 10000000},
    {"mzd_combine", run_mzd_combine, "Omn,ri,wi,R,ri,R,ri", "W", 10000000},
    {"mzd_read_bits", run_mzd_read_bits, "Rmn,ri,ci,n", "", 10000000},
    {"mzd_read_bits_int", run_mzd_read_bits_int, "Rmn,ri,ci,n", "", 10000000},
    {"mzd_xor_bits", run_mzd_xor_bits, "Rmn,ri,ci,n,w", "", 10000000},
    {"mzd_and_bits", run_mzd_and_bits, "Rmn,ri,ci,n,w", "", 10000000},
    {"mzd_clear_bits", run_mzd_clear_bits, "Rmn,ri,ci,n", "", 10000000},
    {"mzd_is_zero", run_mzd_is_zero, "Rmn", "mn", 10000000},
    {"mzd_find_pivot", run_mzd_find_pivot, "Rmn,ri,ci", "", 1000000},
    {"mzd_density", run_mzd_density, "Rmn,wi", "", 10000000},
    {"_mzd_density", run__mzd_density, "Rmn,wi,ri,ci", "", 10000000},
    {"mzd_first_zero_row", run_mzd_first_zero_row, "Rmm", "m", 10000000000},
    {"nothing", run_bench_nothing, "", "", 1}};

int decode_size(char var, struct test_params *params, int *argcp, char ***argvp) {
  if (*argcp < 1) {
    fprintf(stderr, "%s: Not enough arguments. Expected matrix size: %c\n", progname, var);
    return 1;
  }
  --(*argcp);
  switch (var) {
  case 'k': params->k = atoi((*argvp)[0]); break;
  case 'l': params->l = atoi((*argvp)[0]); break;
  case 'm': params->m = atoi((*argvp)[0]); break;
  case 'n': params->n = atoi((*argvp)[0]); break;
  }
  ++(*argvp);
  return 0;
}

int decode_index(char idx, struct test_params *params, int *argcp, char ***argvp) {
  if (*argcp < 1) {
    int count = 0;
    switch (idx) {
    case 'r': count = params->rows; break;
    case 'c': count = params->cols; break;
    case 'w': count = params->wrds; break;
    }
    fprintf(stderr, "%s: Not enough arguments. Expected ", progname);
    switch (idx) {
    case 'r': fprintf(stderr, "row"); break;
    case 'c': fprintf(stderr, "column"); break;
    case 'w': fprintf(stderr, "word"); break;
    }
    fprintf(stderr, " index : %c%d\n", idx, count + 1);
    return 1;
  }
  --(*argcp);
  switch (idx) {
  case 'r': params->row[params->rows++] = atoi((*argvp)[0]); break;
  case 'c': params->col[params->cols++] = atoi((*argvp)[0]); break;
  case 'w': params->wrd[params->wrds++] = atoi((*argvp)[0]); break;
  }
  ++(*argvp);
  return 0;
}

int decode_code(char idx, struct test_params *params, int *argcp, char ***argvp) {
  if (*argcp < 1) {
    fprintf(stderr, "%s: Not enough arguments. Expected ", progname);
    switch (idx) {
    case 'b': printf("boolean"); break;
    case 'n': printf("integer"); break;
    default: printf("%c", idx);
    }
    fprintf(stderr, ".\n");
    return 1;
  }
  --(*argcp);
  switch (idx) {
  case 'b':
    params->boolean = atoi((*argvp)[0]);
    if (params->boolean != 0 && params->boolean != 1) {
      fprintf(stderr, "%s: Expected boolean: %s\n", progname, (*argvp)[0]);
      return 1;
    }
    break;
  case 'n': params->integer = atoi((*argvp)[0]); break;
  }
  ++(*argvp);
  return 0;
}

int main(int argc, char **argv) {
  int opts = global_options(&argc, &argv);

  struct test_params params;
  unsigned long long data[8];
  int data_len;

#ifdef HAVE_LIBPAPI
  int papi_counters = PAPI_num_counters();
  if (papi_counters < papi_array_len) {
    fprintf(stderr, "%s: Warning: there are only %d hardware counters available!\n", progname,
            papi_counters);
    papi_array_len = papi_counters;
  }
  int res = PAPI_start_counters((int *)papi_events, papi_array_len);
  switch (res) {
  case 0: {
    long long *tmp = (long long *)malloc(papi_array_len * sizeof(long long));
    PAPI_stop_counters(tmp, papi_array_len);
    free(tmp);
    break;
  }
  case PAPI_ECNFLCT: {
    fprintf(stderr,
            "%s: %s: Conflicting event: The underlying counter hardware cannot count the specified "
            "events simultaneously.\n",
            progname, papi_event_name(papi_events[papi_array_len - 1]));
    fprintf(stderr, "Run `papi_event_chooser PRESET");
    for (int nv = 0; nv < papi_array_len - 1; ++nv)
      fprintf(stderr, " %s", papi_event_name(papi_events[nv]));
    fprintf(stderr, "` to get a list of possible events that can be added.\n");
    break;
  }
  case PAPI_ENOEVNT: {
    for (int nv = 0; nv < papi_array_len; ++nv)
      if ((res = PAPI_query_event(papi_events[nv])) != PAPI_OK) {
        fprintf(stderr, "%s: PAPI_start_counters: %s: %s.\n", progname,
                papi_event_name(papi_events[nv]), PAPI_strerror(res));
        break;
      }
    break;
  }
  case PAPI_ESYS:
    fprintf(stderr, "%s: PAPI_start_counters: %s\n", progname, strerror(errno));
    break;
  default: fprintf(stderr, "%s: PAPI_start_counters: %s.\n", progname, PAPI_strerror(res)); break;
  }
  if (res) return 1;
  for (int nv = 0; nv <= papi_array_len; ++nv) loop_calibration[nv] = 100000000;
  params.count = 1;
  data_len     = papi_array_len + 1;
  for (int i = 0; i < 100; ++i) run_bench_nothing((void *)&params, data, &data_len);
#endif

  int f;
  int found = 0;

  params.rows   = 0;
  params.cols   = 0;
  params.wrds   = 0;
  params.cutoff = -1;

  if (argc >= 2) {
    params.funcname = argv[1];

    for (f = 0; f < sizeof(function_mapper) / sizeof(function_mapper[0]); ++f) {
      if (strcmp(params.funcname, function_mapper[f].funcname) == 0) {
        found = 1;
        break;
      }
    }
  }
  if (!found) {
    if (argc >= 2)
      fprintf(stderr, "%s: function name \"%s\" not found.\n", progname, params.funcname);
    else {
      fprintf(stderr, "Usage: %s [OPTIONS] <funcname> [ARGS]\n", progname);
      bench_print_global_options(stderr);
    }
    fprintf(stderr, "Possible values for <funcname>:\n");
    for (f = 0; f < sizeof(function_mapper) / sizeof(function_mapper[0]); ++f) {
      if (f != 0 && f % 4 == 0) fprintf(stderr, "\n");
      fprintf(stderr, "%-22s", function_mapper[f].funcname);
    }
    fprintf(stderr, "\n");
    return 1;
  }

  argc -= 2;  // argc >= 1 if more arguments.
  argv += 2;  // Next argument in argv[0]
  char *input_codes = strdup(function_mapper[f].input_codes);
  char *input_code[10];
  char *p   = input_codes;
  int codes = 0;
  while (*p) {
    input_code[codes++] = p++;
    while (*p && *p != ',') ++p;
    if (*p == ',') *p++ = '\0';
  }
  int saw_var[4];
  for (int var_index = 0; var_index < 4; ++var_index) saw_var[var_index] = 0;
  int saw_vars = 0;
  char usage[64];
  char *usage_ptr = usage;
  int error       = 0;
  for (int c = 0;; ++c) {
    if (c < codes) {
      p = input_code[c];
      if (isupper(*p)) {
        while (*++p) {
          if (*p != '1') {
            int var_index = *p - 'k';
            assert(var_index >= 0 && var_index <= 3);  // 'k', 'l', 'm' or 'n'.
            saw_var[var_index] = 1;
            saw_vars           = 1;
          }
        }
        continue;
      }
    }
    if (saw_vars) {
      saw_vars = 0;
      for (int var_count = 2; var_count < 6; ++var_count) {
        int var_index = var_count % 4;
        if (saw_var[var_index] == 1) {
          *usage_ptr++       = ' ';
          *usage_ptr++       = 'k' + var_index;
          saw_var[var_index] = 2;
          if (!error && decode_size('k' + var_index, &params, &argc, &argv)) error = 1;
        }
      }
    }
    if (c == codes) break;
    if (p[1] == 'i') {
      *usage_ptr++ = ' ';
      *usage_ptr++ = *p;
      switch (*p) {
      case 'r':
        *usage_ptr++ = '1' + params.rows;
        if (error) ++params.rows;
        break;
      case 'c':
        *usage_ptr++ = '1' + params.cols;
        if (error) ++params.cols;
        break;
      case 'w':
        *usage_ptr++ = '1' + params.wrds;
        if (error) ++params.wrds;
        break;
      }
      if (!error && decode_index(*p, &params, &argc, &argv)) error = 1;
    } else {
      *usage_ptr++ = ' ';
      *usage_ptr++ = *p;
      if (!error && decode_code(*p, &params, &argc, &argv)) error = 1;
    }
  }
  *usage_ptr = '\0';
  if (argc != 0) error = 1;
  if (error) {
    if (argc != 0) fprintf(stderr, "%s %s: too many parameters.\n", progname, params.funcname);
    fprintf(stderr, "Usage: %s [OPTIONS] %s%s\n", progname, params.funcname, usage);
    if (opts <= 0) bench_print_global_options(stderr);
    return 1;
  }

  double cost  = complexity(&params, function_mapper[f].complexity_code);
  params.count = bench_count ? bench_count : function_mapper[f].count / cost;
  if (params.count < 1) params.count = 1;
  bench_count = params.count;

  srandom(17);

  data_len = run_bench(function_mapper[f].run_func, (void *)&params, data,
                       sizeof(data) / sizeof(unsigned long long));

  printf("function: %s, count: %" PRId64 ", ", params.funcname, params.count);
  if (saw_var[2]) printf("m: %d, ", params.m);
  if (saw_var[3]) printf("n: %d, ", params.n);
  if (saw_var[0]) printf("k: %d, ", params.k);
  if (saw_var[1]) printf("l: %d, ", params.l);
  for (int i = 0; i < 3; ++i) {
    if (i < params.rows) printf("row%c: %d, ", 'a' + i, params.row[i]);
    if (i < params.cols) printf("col%c: %d, ", 'a' + i, params.col[i]);
    if (i < params.wrds) printf("word%c: %d, ", 'a' + i, params.wrd[i]);
  }
  if (params.cutoff != -1) printf("cutoff: %d, ", params.cutoff);
  print_wall_time(data[0] / 1000000.0 / params.count);
  printf(", cpu cycles: %llu", (data[1] + params.count / 2) / params.count);
#ifndef HAVE_LIBPAPI
  printf(", cc/");
  print_complexity_human(&params, function_mapper[f].complexity_code);
  printf(": %f\n", data[1] / (params.count * cost));
#else
  printf("\n");
  for (int n = 1; n < data_len; ++n) {
    printf("%s (%f) per bit (divided by ", papi_event_name(papi_events[n - 1]),
           (double)data[n] / params.count);
    print_complexity_human(&params, function_mapper[f].complexity_code);
    printf("): %f\n", data[n] / (params.count * cost));
  }
#endif
}
