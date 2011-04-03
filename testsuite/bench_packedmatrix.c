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

#include <stdlib.h>
#include <ctype.h>
#include <inttypes.h>

#include "config.h"
#include "cpucycles.h"
#include "m4ri.h"
#include "benchmarketing.h"

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
  char const* funcname;
};

typedef int (*run_type)(void*, double*, unsigned long long*);

#define TIME(mzd_func, ARGS, count) do {\
    *wt = walltime(0.0);		\
    *cycles = cpucycles();		\
    for (uint64_t i = 0; i < count; ++i)	\
      mzd_func ARGS;			\
    *wt = walltime(*wt);		\
    *cycles = cpucycles() - *cycles;	\
  } while(0)

int run__mzd_row_swap(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const rowa = p->row[0];
  rci_t const rowb = p->row[1];
  wi_t const startblock = p->wrd[0];
  uint64_t const count = p->count;

  TIME(_mzd_row_swap, (A, rowa, rowb, startblock), count);

  mzd_free(A);
  return 0;
}

int run_mzd_row_swap(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const rowa = p->row[0];
  rci_t const rowb = p->row[1];
  uint64_t const count = p->count;

  TIME(mzd_row_swap, (A, rowa, rowb), count);

  mzd_free(A);
  return 0;
}

int run_mzd_copy_row(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_t* const B = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const rowa = p->row[0];
  rci_t const rowb = p->row[1];
  uint64_t const count = p->count;

  TIME(mzd_copy_row, (B, rowb, A, rowa), count);

  mzd_free(A);
  mzd_free(B);
  return 0;
}

int run_mzd_col_swap(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const cola = p->col[0];
  rci_t const colb = p->col[1];
  uint64_t const count = p->count;

  TIME(mzd_col_swap, (A, cola, colb), count);

  mzd_free(A);
  return 0;
}

int run_mzd_col_swap_in_rows(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const cola = p->col[0];
  rci_t const colb = p->col[1];
  rci_t const start_row = p->row[0];
  rci_t const stop_row = p->row[1];
  uint64_t const count = p->count;

  TIME(mzd_col_swap_in_rows, (A, cola, colb, start_row, stop_row), count);

  mzd_free(A);
  return 0;
}

int run_mzd_read_bit(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const rowa = p->row[0];
  rci_t const cola = p->col[0];
  uint64_t const count = p->count;
  BIT volatile result;

  TIME(result = mzd_read_bit, (A, rowa, cola), count);

  mzd_free(A);
  return 0;
}

int run_mzd_write_bit(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  rci_t const rowa = p->row[0];
  rci_t const cola = p->col[0];
  int bit = 0;
  uint64_t const count = p->count;

  TIME(mzd_write_bit, (A, rowa, cola, bit); bit = !bit, count);

  mzd_free(A);
  return 0;
}

int run_mzd_row_add_offset(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const rowa = p->row[0];
  rci_t const rowb = p->row[1];
  rci_t const cola = p->col[0];
  uint64_t const count = p->count;

  TIME(mzd_row_add_offset, (A, rowa, rowb, cola), count);

  mzd_free(A);
  return 0;
}

int run_mzd_row_add(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const rowa = p->row[0];
  rci_t const rowb = p->row[1];
  uint64_t const count = p->count;

  TIME(mzd_row_add, (A, rowa, rowb), count);

  mzd_free(A);
  return 0;
}

int run_mzd_transpose(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->m);
  mzd_t* const B = mzd_init(p->m, p->m);
  mzd_randomize(A);
  uint64_t const count = p->count;

  TIME(mzd_transpose, (B, A), count);

  mzd_free(A);
  mzd_free(B);
  return 0;
}

int run_mzd_mul_naive(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->l);
  mzd_t* const B = mzd_init(p->l, p->n);
  mzd_t* const C = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_randomize(B);
  uint64_t const count = p->count;

  TIME(mzd_mul_naive, (C, A, B), count);

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  return 0;
}

int run_mzd_addmul_naive(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->l);
  mzd_t* const B = mzd_init(p->l, p->n);
  mzd_t* const C = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_randomize(B);
  uint64_t const count = p->count;

  TIME(mzd_addmul_naive, (C, A, B), count);

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  return 0;
}

int run__mzd_mul_naive(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->l);
  mzd_t* const B = mzd_init(p->n, p->l);
  mzd_t* const C = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_randomize(B);
  int const clear = p->boolean;
  uint64_t const count = p->count;

  TIME(_mzd_mul_naive, (C, A, B, clear), count);

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  return 0;
}

int run__mzd_mul_va(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_t* const V = mzd_init(1, p->m);
  mzd_t* const C = mzd_init(1, p->n);
  mzd_randomize(A);
  mzd_randomize(V);
  int const clear = p->boolean;
  uint64_t const count = p->count;

  TIME(_mzd_mul_va, (C, V, A, clear), count);

  mzd_free(A);
  mzd_free(V);
  mzd_free(C);
  return 0;
}

int run_mzd_gauss_delayed(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const cola = p->col[0];
  int const full = p->boolean;
  uint64_t const count = p->count;
  rci_t result;

  TIME(result = mzd_gauss_delayed, (A, cola, full), count);

  mzd_free(A);
  return 0;
}

int run_mzd_echelonize_naive(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  int const full = p->boolean;
  uint64_t const count = p->count;
  rci_t result;

  TIME(result = mzd_echelonize_naive, (A, full), count);

  mzd_free(A);
  return 0;
}

int run_mzd_equal(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_t* const B = mzd_copy(NULL, A);
  uint64_t const count = p->count;
  int volatile result;

  TIME(result = mzd_equal, (A, B), count);

  mzd_free(A);
  mzd_free(B);
  return 0;
}

int run_mzd_cmp(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_t* const B = mzd_copy(NULL, A);
  uint64_t const count = p->count;
  int volatile result;

  TIME(result = mzd_cmp, (A, B), count);

  mzd_free(A);
  mzd_free(B);
  return 0;
}

int run_mzd_copy(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_t* const B = mzd_init(p->m, p->n);
  mzd_randomize(A);
  uint64_t const count = p->count;

  TIME(mzd_copy, (B, A), count);

  mzd_free(A);
  mzd_free(B);
  return 0;
}

int run_mzd_concat(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->k);
  mzd_t* const B = mzd_init(p->m, p->l);
  mzd_t* const C = mzd_init(p->m, p->k + p->l);
  mzd_randomize(A);
  mzd_randomize(B);
  uint64_t const count = p->count;

  TIME(mzd_concat, (C, A, B), count);

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  return 0;
}

int run_mzd_stack(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->k, p->n);
  mzd_t* const B = mzd_init(p->l, p->n);
  mzd_t* const C = mzd_init(p->k + p->l, p->n);
  mzd_randomize(A);
  mzd_randomize(B);
  uint64_t const count = p->count;

  TIME(mzd_stack, (C, A, B), count);

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  return 0;
}

int run_mzd_submatrix(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t const rowa = p->row[0];
  rci_t const cola = p->col[0];
  rci_t const rowb = p->row[1];
  rci_t const colb = p->col[1];
  mzd_t* const S = mzd_init(rowb - rowa, colb - cola);
  uint64_t const count = p->count;

  TIME(mzd_submatrix, (S, A, rowa, cola, rowb, colb), count);

  mzd_free(A);
  mzd_free(S);
  return 0;
}

int run_mzd_invert_naive(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->m);
  mzd_t* const I = mzd_init(p->m, p->m);
  mzd_t* const C = mzd_init(p->m, p->m);
  mzd_randomize(A);
  mzd_set_ui(I, 1);
  uint64_t const count = p->count;

  TIME(mzd_invert_naive, (C, A, I), count);

  mzd_free(A);
  mzd_free(I);
  mzd_free(C);
  return 0;
}

int run_mzd_add(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_t* const B = mzd_init(p->m, p->n);
  mzd_t* const C = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_randomize(B);
  uint64_t const count = p->count;

  TIME(mzd_add, (C, A, B), count);

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  return 0;
}

int run__mzd_add(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_t* const B = mzd_init(p->m, p->n);
  mzd_t* const C = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_randomize(B);
  uint64_t const count = p->count;

  TIME(_mzd_add, (C, A, B), count);

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  return 0;
}

int run_mzd_combine(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_t* const B = mzd_init(p->m, p->n);
  mzd_t* const C = mzd_init(p->m, p->n);
  mzd_randomize(A);
  mzd_randomize(B);
  rci_t row1 = p->row[0];
  rci_t row2 = p->row[1];
  rci_t row3 = p->row[2];
  wi_t startblock = p->wrd[0];
  uint64_t const count = p->count;

  TIME(mzd_combine, (C, row3, startblock, A, row1, startblock, B, row2, startblock), count);

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  return 0;
}

int run_mzd_read_bits(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t row = p->row[0];
  rci_t col = p->col[0];
  int n = p->integer;
  uint64_t const count = p->count;
  word volatile result;

  TIME(result = mzd_read_bits, (A, row, col, n), count);

  mzd_free(A);
  return 0;
}

int run_mzd_read_bits_int(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t row = p->row[0];
  rci_t col = p->col[0];
  int n = p->integer;
  uint64_t const count = p->count;
  int volatile result;

  TIME(result = mzd_read_bits_int, (A, row, col, n), count);

  mzd_free(A);
  return 0;
}

int run_mzd_xor_bits(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t row = p->row[0];
  rci_t col = p->col[0];
  int n = p->integer;
  word volatile const values = 0xffffffffffffffffULL;
  uint64_t const count = p->count;

  TIME(mzd_xor_bits, (A, row, col, n, values), count);

  mzd_free(A);
  return 0;
}

int run_mzd_and_bits(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t row = p->row[0];
  rci_t col = p->col[0];
  int n = p->integer;
  word volatile const values = 0xffffffffffffffffULL;
  uint64_t const count = p->count;

  TIME(mzd_and_bits, (A, row, col, n, values), count);

  mzd_free(A);
  return 0;
}

int run_mzd_clear_bits(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t row = p->row[0];
  rci_t col = p->col[0];
  int n = p->integer;
  uint64_t const count = p->count;

  TIME(mzd_clear_bits, (A, row, col, n), count);

  mzd_free(A);
  return 0;
}

int run_mzd_is_zero(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  uint64_t const count = p->count;
  int volatile result;

  TIME(result = mzd_is_zero, (A), count);

  mzd_free(A);
  return 0;
}

int run_mzd_row_clear_offset(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t row = p->row[0];
  rci_t col = p->col[0];
  uint64_t const count = p->count;

  TIME(mzd_row_clear_offset, (A, row, col), count);

  mzd_free(A);
  return 0;
}

int run_mzd_find_pivot(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t row = p->row[0];
  rci_t col = p->col[0];
  uint64_t const count = p->count;
  int volatile result;
  rci_t row_out;
  rci_t col_out;

  TIME(result = mzd_find_pivot, (A, row, col, &row_out, &col_out), count);

  mzd_free(A);
  return 0;
}

int run_mzd_density(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  wi_t res = p->wrd[0];
  uint64_t const count = p->count;
  double volatile result;

  TIME(result = mzd_density, (A, res), count);

  mzd_free(A);
  return 0;
}

int run__mzd_density(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_randomize(A);
  rci_t row = p->row[0];
  rci_t col = p->col[0];
  wi_t res = p->wrd[0];
  uint64_t const count = p->count;
  double volatile result;

  TIME(result = _mzd_density, (A, res, row, col), count);

  mzd_free(A);
  return 0;
}

int run_mzd_first_zero_row(void *_p, double *wt, unsigned long long *cycles)
{
  struct test_params *p = (struct test_params *)_p;

  mzd_t* const A = mzd_init(p->m, p->n);
  mzd_set_ui(A, 1);
  uint64_t const count = p->count;
  rci_t volatile result;

  TIME(mzd_first_zero_row, (A), count);

  mzd_free(A);
  return 0;
}

// Returns a number proportional with the ideal number of
// mathematical operations for the given code.
double complexity1(struct test_params *p, char code)
{
  switch(code)
  {
    case 'k':
      return p->k;			// Linear with size 'k' of a matrix.
    case 'l':
      return p->l;			// Linear with size 'l' of a matrix.
    case 'm':
      return p->m;			// Linear with the number of rows of the matrix.
    case 'n':
      return p->n;			// Linear with the number of columns of the matrix.
    case 'W':
      assert(p->n > RADIX * p->wrd[0]);	// Linear with the number of processed columns.
      return p->n - RADIX * p->wrd[0];
    case 'D':
      assert(p->row[0] < p->row[1]);
      return p->row[1] - p->row[0];	// Linear with the number of rows between start_row and stop_row.
    case 'E':
      assert(p->col[0] < p->col[1]);
      return p->col[1] - p->col[0];	// Linear with the number of cols between start_col and stop_col.
    case 'C':
      assert(p->col[0] < p->n);
      return p->n - p->col[0];		// Linear with the number of columns of column col and beyond.
  }
}

char const* complexity1_human(struct test_params *p, char code)
{
  switch(code)
  {
    case 'k':
      return "k";
    case 'l':
      return "l";
    case 'm':
      return "m";
    case 'n':
      return "n";
    case 'W':
      return "cols";
    case 'D':
      return "rows";
    case 'E':
      return "cols";
    case 'C':
      return "cols";
  }
}

double complexity(struct test_params *p, char const* cp)
{
  double c = 1;
  while (*cp)
  {
    c *= complexity1(p, *cp);
    ++cp;
  }
  return c;
}

void print_complexity_human(struct test_params *p, char const* cp)
{
  int first = 1;
  char last_cp = 0;
  int power = 0;
  while (*cp)
  {
    if (*cp != last_cp)
    {
      if (power > 1)
	printf("^%d", power);
      if (!first && isupper(*cp))
	printf("*");
      printf("%s", complexity1_human(p, *cp));
      power = 0;
      last_cp = *cp;
    }
    ++power;
    ++cp;
  }
  if (power > 1)
    printf("^%d", power);
}

struct function_st {
  char const* funcname;
  run_type run_func;
  char const* input_codes;
  char const* complexity_code;
  uint64_t count;
};

typedef struct function_st function_st;

static function_st function_mapper[] = {
  { "_mzd_row_swap",        run__mzd_row_swap,        "Rmn,ri,ri,wi",    "W",  1000000000 },
  { "mzd_row_swap",         run_mzd_row_swap,         "Rmn,ri,ri",       "n",  1000000000 },
  { "mzd_copy_row",         run_mzd_copy_row,         "Omn,ri,R,ri",     "n",  1000000000 },
  { "mzd_col_swap",         run_mzd_col_swap,         "Rmn,ci,ci",       "m",  10000000 },
  { "mzd_col_swap_in_rows", run_mzd_col_swap_in_rows, "Rmn,ci,ci,ri,ri", "D",  10000000 },
  { "mzd_read_bit",         run_mzd_read_bit,         "Rmn,ri,ci",       "",   100000000 },
  { "mzd_write_bit",        run_mzd_write_bit,        "Omn,ri,ci",       "",   100000000 },
  { "mzd_row_add_offset",   run_mzd_row_add_offset,   "Rmn,ri,ri,ci",    "C",  100000000 },
  { "mzd_row_add",          run_mzd_row_add,          "Rmn,ri,ri",       "n",  100000000 },
  { "mzd_transpose",        run_mzd_transpose,        "Omm,Rmm",         "mm", 10000000 },
  { "mzd_mul_naive",        run_mzd_mul_naive,        "Omn,Rml,Rln",     "mnl",10000000 },
  { "mzd_addmul_naive",     run_mzd_addmul_naive,     "Omn,Rml,Rln",     "mnl",10000000 },
  { "_mzd_mul_naive",       run__mzd_mul_naive,       "Omn,Rml,Rnl,b",   "mnl",10000000 },
  { "_mzd_mul_va",          run__mzd_mul_va,          "O1n,V1m,Amn,b",   "mn", 1000000000 },
  { "mzd_gauss_delayed",    run_mzd_gauss_delayed,    "Rmn,ci,b",        "mC", 10000000 },
  { "mzd_echelonize_naive", run_mzd_echelonize_naive, "Rmn,b",           "mn", 10000000 },
  { "mzd_equal",            run_mzd_equal,            "Rmn,Rmn",         "mn", 1000000000 },
  { "mzd_cmp",              run_mzd_cmp,              "Rmn,Rmn",         "mn", 1000000000 },
  { "mzd_copy",             run_mzd_copy,             "Omn,Rmn",         "mn", 10000000 },
  { "mzd_concat",           run_mzd_concat,           "Omn,Rmk,Rml",     "mn", 10000000 },
  { "mzd_stack",            run_mzd_stack,            "Omn,Rkn,Rln",     "mn", 10000000 },
  { "mzd_submatrix",        run_mzd_submatrix,      "O,Rmn,ri,ci,ri,ci", "DE", 10000000 },
  { "mzd_invert_naive",     run_mzd_invert_naive,     "Omm,Rmm,Imm",     "mmm",10000000 },
  { "mzd_add",              run_mzd_add,              "Omn,Rmn,Rmn",     "mn", 10000000 },
  { "_mzd_add",             run__mzd_add,             "Omn,Rmn,Rmn",     "mn", 10000000 },
  { "mzd_combine",          run_mzd_combine,      "Omn,ri,wi,R,ri,R,ri", "W",  10000000 },
  { "mzd_read_bits",        run_mzd_read_bits,        "Rmn,ri,ci,n",     "",   10000000 },
  { "mzd_read_bits",        run_mzd_read_bits_int,    "Rmn,ri,ci,n",     "",   10000000 },
  { "mzd_xor_bits",         run_mzd_xor_bits,         "Rmn,ri,ci,n,w",   "",   10000000 },
  { "mzd_and_bits",         run_mzd_and_bits,         "Rmn,ri,ci,n,w",   "",   10000000 },
  { "mzd_clear_bits",       run_mzd_clear_bits,       "Rmn,ri,ci,n",     "",   10000000 },
  { "mzd_is_zero",          run_mzd_is_zero,          "Rmn",             "mn", 10000000 },
  { "mzd_row_clear_offset", run_mzd_row_clear_offset, "Omn,ri,ci",       "C",  10000000 },
  { "mzd_find_pivot",       run_mzd_find_pivot,       "Rmn,ri,ci",       "",   10000000 },
  { "mzd_density",          run_mzd_density,          "Rmn,wi",          "",   10000000 },
  { "_mzd_density",         run__mzd_density,         "Rmn,wi,ri,ci",    "",   10000000 },
  { "mzd_first_zero_row",   run_mzd_first_zero_row,   "Rmn",             "m",  10000000 }
};

int decode_size(char var, struct test_params* params, int* argcp, char*** argvp)
{
  if (*argcp < 1)
  {
    fprintf(stderr, "%s: Not enough arguments. Expected matrix size: %c\n", progname, var);
    return 1;
  }
  --(*argcp);
  switch(var)
  {
    case 'k':
      params->k = atoi((*argvp)[0]);
      break;
    case 'l':
      params->l = atoi((*argvp)[0]);
      break;
    case 'm':
      params->m = atoi((*argvp)[0]);
      break;
    case 'n':
      params->n = atoi((*argvp)[0]);
      break;
  }
  ++(*argvp);
  return 0;
}

int decode_index(char idx, struct test_params* params, int* argcp, char*** argvp)
{
  if (*argcp < 1)
  {
    int count;
    switch(idx)
    {
      case 'r':
	count = params->rows;
	break;
      case 'c':
	count = params->cols;
	break;
      case 'w':
	count = params->wrds;
	break;
    }
    fprintf(stderr, "%s: Not enough arguments. Expected ", progname);
    switch(idx)
    {
      case 'r':
	fprintf(stderr, "row");
	break;
      case 'c':
	fprintf(stderr, "column");
	break;
      case 'w':
	fprintf(stderr, "word");
	break;
    }
    fprintf(stderr, " index : %c%d\n", idx, count + 1);
    return 1;
  }
  --(*argcp);
  switch(idx)
  {
    case 'r':
      params->row[params->rows++] = atoi((*argvp)[0]);
      break;
    case 'c':
      params->col[params->cols++] = atoi((*argvp)[0]);
      break;
    case 'w':
      params->wrd[params->wrds++] = atoi((*argvp)[0]);
      break;
  }
  ++(*argvp);
  return 0;
}

int decode_code(char idx, struct test_params* params, int* argcp, char*** argvp)
{
  if (*argcp < 1)
  {
    fprintf(stderr, "%s: Not enough arguments. Expected ", progname);
    switch(idx)
    {
      case 'b':
	printf("boolean");
	break;
      case 'n':
	printf("integer");
	break;
      default:
	printf("%c", idx);
    }
    fprintf(stderr, ".\n");
    return 1;
  }
  --(*argcp);
  switch(idx)
  {
    case 'b':
      params->boolean = atoi((*argvp)[0]);
      if (params->boolean != 0 && params->boolean != 1)
      {
	fprintf(stderr, "%s: Expected boolean: %s\n", progname, (*argvp)[0]);
	return 1;
      }
      break;
    case 'n':
      params->integer = atoi((*argvp)[0]);
      break;
  }
  ++(*argvp);
  return 0;
}

int main(int argc, char** argv)
{
  int opts = global_options(&argc, &argv);

  int f;
  struct test_params params;
  int found = 0;

  params.rows = 0;
  params.cols = 0;
  params.wrds = 0;
  params.cutoff = -1;

  if (argc >= 2)
  {
    params.funcname = argv[1];

    for (f = 0; f < sizeof(function_mapper) / sizeof(function_mapper[0]); ++f)
    {
      if (strcmp(params.funcname, function_mapper[f].funcname) == 0)
      {
	found = 1;
	break;
      }
    }
  }
  if (!found)
  {
    if (argc >= 2)
      fprintf(stderr, "%s: function name \"%s\" not found.\n", progname, params.funcname);
    else
    {
      fprintf(stderr, "Usage: %s [OPTIONS] <funcname> [ARGS]\n", progname);
      bench_print_global_options(stderr);
    }
    fprintf(stderr, "Possible values for <funcname>:\n");
    for (f = 0; f < sizeof(function_mapper) / sizeof(function_mapper[0]); ++f)
    {
      if (f != 0 && f % 4 == 0)
	fprintf(stderr, "\n");
      fprintf(stderr, "%-22s", function_mapper[f].funcname);
    }
    fprintf(stderr, "\n");
    return 1;
  }

  argc -= 2;	// argc >= 1 if more arguments.
  argv += 2;	// Next argument in argv[0]
  char* input_codes = strdup(function_mapper[f].input_codes);
  char* input_code[10];
  char* p = input_codes;
  int codes = 0;
  while(*p)
  {
    input_code[codes++] = p++;
    while(*p && *p != ',')
      ++p;
    if (*p == ',')
      *p++ = '\0';
  }
  int saw_var[4];
  for (int var_index = 0; var_index < 4; ++var_index)
    saw_var[var_index] = 0;
  int saw_vars = 0;
  char usage[64];
  char* usage_ptr = usage;
  int error = 0;
  for (int c = 0; ; ++c)
  {
    if (c < codes)
    {
      p = input_code[c];
      if (isupper(*p))
      {
	while(*++p)
	{
	  if (*p != '1')
	  {
	    int var_index = *p - 'k';
	    assert(var_index >= 0 && var_index <= 3); // 'k', 'l', 'm' or 'n'.
	    saw_var[var_index] = 1;
	    saw_vars = 1;
	  }
	}
	continue;
      }
    }
    if (saw_vars)
    {
      saw_vars = 0;
      for (int var_count = 2; var_count < 6; ++var_count)
      {
	int var_index = var_count % 4;
	if (saw_var[var_index] == 1)
	{
	  *usage_ptr++ = ' ';
	  *usage_ptr++ = 'k' + var_index;
	  saw_var[var_index] = 2;
	  if (!error && decode_size('k' + var_index, &params, &argc, &argv))
	    error = 1;
	}
      }
    }
    if (c == codes)
      break;
    if (p[1] == 'i')
    {
      *usage_ptr++ = ' ';
      *usage_ptr++ = *p;
      switch(*p)
      {
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
      if (!error && decode_index(*p, &params, &argc, &argv))
	error = 1;
    }
    else
    {
      *usage_ptr++ = ' ';
      *usage_ptr++ = *p;
      if (!error && decode_code(*p, &params, &argc, &argv))
	error = 1;
    }
  }
  *usage_ptr = '\0';
  if (argc != 0)
    error = 1;
  if (error)
  {
    if (argc != 0)
      fprintf(stderr, "%s %s: too many parameters.\n", progname, params.funcname);
    fprintf(stderr, "Usage: %s [OPTIONS] %s%s\n", progname, params.funcname, usage);
    if (opts <= 0)
      bench_print_global_options(stderr);
    return 1;
  }

  double cost = complexity(&params, function_mapper[f].complexity_code);
  params.count = bench_count ? bench_count : function_mapper[f].count / cost;
  if (params.count < 1)
    params.count = 1;

  srandom(17);
  unsigned long long t;
  double wt;

  run_bench(function_mapper[f].run_func, (void*)&params, &wt, &t);

  printf("function: %s, count: %" PRId64 ", ", params.funcname, params.count);
  if (saw_var[2])
    printf("m: %d, ", params.m);
  if (saw_var[3])
    printf("n: %d, ", params.n);
  if (saw_var[0])
    printf("k: %d, ", params.k);
  if (saw_var[1])
    printf("l: %d, ", params.l);
  for (int i = 0; i < 3; ++i)
  {
    if (i < params.rows)
      printf("row%c: %d, ", 'a' + i, params.row[i]);
    if (i < params.cols)
      printf("col%c: %d, ", 'a' + i, params.col[i]);
    if (i < params.wrds)
      printf("word%c: %d, ", 'a' + i, params.wrd[i]);
  }
  if (params.cutoff != -1)
    printf("cutoff: %d, ", params.cutoff);
  printf("cpu cycles per loop: %llu, wall time: %lf, cc/", t / params.count, wt);
  print_complexity_human(&params, function_mapper[f].complexity_code);
  printf(": %f\n", t / (params.count * cost));
}
