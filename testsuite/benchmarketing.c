/*
 * benchmarketing.c
 *
 * Benchmark engine.
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

/*
 * Example usage:
 *
 * ./bench_elimination -m 4 -c 90 -a 0.005 -d -t 30 -n 1000 1000 1000
 *
 * would run at most 30 seconds (-t) or 1000 times (-n), whichever comes
 * first, or stop after the real average falls with 90% certainty (-c) in
 * a range that is +/- 0.005 times the observed mean (-a: accuracry),
 * but no sooner than that at least 4 (-m: minimum) measurements have been
 * done. It would also print (-d: dump) each measurement.
 *
 * Example output.
 *
 * 0.002446
 * 0.002528
 * 0.002373
 * 0.002363
 * 0.002373
 * 0.002364
 * 0.002367
 * 0.002357
 * 0.002354
 * 0.002360
 * 0.002429
 * 0.002373
 * 0.002379
 * 0.002347
 * 0.002367
 * 0.002368
 * 0.002377
 * 0.002373
 * 0.002357
 * 0.002371
 * 0.002370
 * 0.002402
 * 0.002384
 * 0.002370
 * 0.002367
 * 0.002359
 * 0.002470
 * 0.002372
 * 0.002361
 * 0.002364
 * 0.002373
 * Total running time:  1.092 seconds.
 * Sample size: 31; mean: 0.002381; standard deviation: 0.000038
 * 90% confidence interval: +/- 0.000012 (0.5%): [0.002370..0.002393]
 *
 * The last three lines can be suppressed by passing the option -q (quiet).
 */

#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "benchmarketing.h"
#include "misc.h"

enum { C80, C90, C95, C98, C99 };

/*
 * Command line option decoding
 */

int bench_quiet = 0;			// Set if -q is used.
int bench_dump = 0;			// Set if -d is used.
int bench_minimum = 2;			// Minimum number of measurements. Set with -m <minimum>.
int bench_maximum = 1000;		// Maximum number of measurements. Set with -n <maximum>.
double bench_maxtime = 60.0;		// Maximum number of seconds to run. Set with -t <maxtime>
double bench_accuracy = 0.01;		// The +/- range (where 1.0 is 100%) within that we want the real population mean to be with the given confidence. Set with -a <bench_accuracy>
int bench_confidence_index = C99;	// The confidence that the real mean is within the given (or found) range.
char const* progname;			// Set to argv[0].

void global_options(int* argcp, char*** argvp)
{
  progname = (*argvp)[0];
  while((*argcp) > 1)
  {
    if ((*argvp)[1][0] != '-' || (*argvp)[1][1] == '\0' || (*argvp)[1][2] != '\0')
      return;
    switch((*argvp)[1][1])
    {
      case 'd':		// Dump
	bench_dump = 1;
	break;
      case 'q':		// Quiet
	bench_quiet = 1;
	break;
      case 'm':
	++*argvp;
	--*argcp;
	bench_minimum = atoi((*argvp)[1]);
	break;
      case 'n':
	++*argvp;
	--*argcp;
	bench_maximum = atoi((*argvp)[1]);
	break;
      case 't':
	++*argvp;
	--*argcp;
	bench_maxtime = strtod((*argvp)[1], NULL);
	break;
      case 'a':
	++*argvp;
	--*argcp;
	bench_accuracy = strtod((*argvp)[1], NULL);
	break;
      case 'c':
      {
	++*argvp;
	--*argcp;
	int confidence = atoi((*argvp)[1]);
	switch (confidence)
	{
	  case 80:
	    bench_confidence_index = C80;
	    break;
	  case 90:
	    bench_confidence_index = C90;
	    break;
	  case 95:
	    bench_confidence_index = C95;
	    break;
	  case 98:
	    bench_confidence_index = C98;
	    break;
	  case 99:
	    bench_confidence_index = C99;
	    break;
	  default:
	    m4ri_die("The only possible confidence percentages are 80, 90, 95, 98 and 99%\n");
	    break;
	}
	break;
      }
      default:
	return;
    }
    ++*argvp;
    --*argcp;
  }
}

/*
 * vector implementation
 *
 * vector_create:	Create vector of size s.
 * vector_destruct:	Destruct vector.
 * vector_resize:	Resize internal allocation.
 * vector_size:		Return number of elements.
 * vector_pushback:	Add one element at the end.
 * vector_get:		Get element at position index.
 */

struct vector_st {
  size_t alloc_size;
  size_t size;
  double* data;
};

typedef struct vector_st* vector;

vector vector_create(size_t s)
{
  vector v = (vector)malloc(sizeof(struct vector_st));
  v->alloc_size = s;
  v->data = s ? (double*)malloc(sizeof(double) * s) : NULL;
  v->size = 0;
  return v;
}

void vector_destruct(vector v)
{
  free(v->data);
}

void vector_resize(vector v, size_t s)
{
  v->data = (double*)realloc(v->data, sizeof(double) * s);
  v->alloc_size = s;
  if (v->size > v->alloc_size)
    v->size = v->alloc_size;
}

static inline size_t vector_size(vector v)
{
  return v->size;
}

void vector_pushback(vector v, double d)
{
  if (++(v->size) > v->alloc_size)
    vector_resize(v, v->alloc_size * 2);
  v->data[v->size - 1] = d;
}

static inline double vector_get(vector v, int index)
{
  return v->data[index];
}

/*
 * Normal distribution
 *
 * normal_calculate:	Calculate the mean and standard deviation of the data in vector v.
 *
 * Returns -1 on failure (not enough data points), 0 otherwise.
 */

struct normal_st {
  int size;
  double mean;
  double sigma;
};

typedef struct normal_st normal;

int normal_calculate(vector v, normal* dist)
{
  dist->size = vector_size(v);

  if (dist->size < 2)
    return -1;

  // Calculate the sum of all data.
  double sum = 0;
  for (int i = 0; i < dist->size; ++i)
    sum += vector_get(v, i);
  dist->mean = sum / dist->size;

  // Calculate the sum of the square of all differences with mean.
  sum = 0;
  for (int i = 0; i < dist->size; ++i)
  {
    double delta = vector_get(v, i) - dist->mean;
    sum += delta * delta;
  }
  dist->sigma = sqrt(sum / (dist->size - 1));

  return 0;
}

/*
 * T-Table
 */

static float student_t[5][34] = {
 { 3.078, 1.886, 1.638, 1.533,
   1.476, 1.440, 1.415, 1.397, 1.383,
   1.372, 1.363, 1.356, 1.350, 1.345,
   1.341, 1.337, 1.333, 1.330, 1.328,
   1.325, 1.323, 1.321, 1.319, 1.318,
   1.316, 1.315, 1.314, 1.313, 1.311,
   1.310, 1.303, 1.296, 1.289, 1.282 },
 { 6.314, 2.920, 2.353, 2.132,
   2.015, 1.943, 1.895, 1.860, 1.833,
   1.812, 1.796, 1.782, 1.771, 1.761,
   1.753, 1.746, 1.740, 1.734, 1.729,
   1.725, 1.721, 1.717, 1.714, 1.711,
   1.708, 1.706, 1.703, 1.701, 1.699,
   1.697, 1.684, 1.671, 1.658, 1.645 },
 { 12.706, 4.303, 3.182, 2.776,
   2.571, 2.447, 2.365, 2.306, 2.262,
   2.228, 2.201, 2.179, 2.160, 2.145,
   2.131, 2.120, 2.110, 2.101, 2.093,
   2.086, 2.080, 2.074, 2.069, 2.064,
   2.060, 2.056, 2.052, 2.048, 2.045,
   2.042, 2.021, 2.000, 1.980, 1.960 },
 { 31.821, 6.965, 4.541, 3.747,
   3.365, 3.143, 2.998, 2.896, 2.821,
   2.764, 2.718, 2.681, 2.650, 2.624,
   2.602, 2.583, 2.567, 2.552, 2.539,
   2.528, 2.518, 2.508, 2.500, 2.492,
   2.485, 2.479, 2.473, 2.467, 2.462,
   2.457, 2.423, 2.390, 2.358, 2.326 },
 { 63.657, 9.925, 5.841, 4.604,
   4.032, 3.707, 3.499, 3.355, 3.250,
   3.169, 3.106, 3.055, 3.012, 2.977,
   2.947, 2.921, 2.898, 2.878, 2.861,
   2.845, 2.831, 2.819, 2.807, 2.797,
   2.787, 2.779, 2.771, 2.763, 2.756,
   2.750, 2.704, 2.660, 2.617, 2.576 } };

static float student_t_certainty[5] = { 0.2, 0.1, 0.05, 0.02, 0.01 };        // Two-tails.

static float t_table(int confidence_index, int freedoms)
{
  if (freedoms <= 30)
    return student_t[confidence_index][freedoms - 1];
  double a, b, y1, y2, y3;
  long x1, x2;
  long x3 = 0;
  int i;
  if (freedoms <= 60)
  {
    i = 29;
    x1 = 30;
    x2 = 40;
    x3 = 60;
  }
  else if (freedoms <= 120)
  {
    i = 30;
    x1 = 40;
    x2 = 60;
    x3 = 120;
  }
  else
  {
    i = 31;
    x1 = 60;
    x2 = 120;
    /* x3 = infinity */
  }
  y1 = student_t[confidence_index][i];
  y2 = student_t[confidence_index][i + 1];
  y3 = student_t[confidence_index][i + 2];
  if (freedoms <= 120)
  {
    double c, d;
    d =   (x1 * x1 * (x3 - x2) + x2 * x2 * (x1 - x3) + x3 * x3 * (x2 - x1));
    a = - (x1      * (y3 - y2) + x2      * (y1 - y3) + x3      * (y2 - y1)) / d;
    b =   (x1 * x1 * (y3 - y2) + x2 * x2 * (y1 - y3) + x3 * x3 * (y2 - y1)) / d;
    c = y2 - a * x2 * x2 - b * x2;
    return (a * freedoms * freedoms + b * freedoms + c);
  }
  double ln1, ln2;
  ln1 = log(y2 - y3);
  ln2 = log(y1 - y3);
  a = - (     ln1 -      ln2) / (x1 - x2);
  b =   (x1 * ln1 - x2 * ln2) / (x1 - x2);
  return (y3 + exp(a * freedoms + b));
}

/*
 * walltime
 *
 * Author: Marin Albrecht.
 */

double walltime( double t0 )
{
  double mic, time;
  double mega = 0.000001;
  struct timeval tp;
  static long base_sec = 0;
  static long base_usec = 0;

  (void)gettimeofday(&tp, NULL);
  if (base_sec == 0)
  {
    base_sec = tp.tv_sec;
    base_usec = tp.tv_usec;
  }

  time = (double)(tp.tv_sec - base_sec);
  mic = (double)(tp.tv_usec - base_usec);
  time = (time + mic * mega) - t0;
  return time;
}

/*
 * run_bench
 *
 * Benchmark main loop.
 */

void run_bench(
    int (*f)(void* params, double* wtp, unsigned long long* ccp),
    void* params,
    double* wtp,
    unsigned long long* ccp)
{
  double const CONFIDENCE = 1.0 - student_t_certainty[bench_confidence_index];
  unsigned long long cc_sum = 0;
  double wt_diff = 0.0;
  vector wt_data = vector_create(128);
  normal stats;
  double start_walltime = walltime(0.0);
  for (int n = 1; n <= bench_maximum; ++n)
  {
    if (!bench_quiet && !bench_dump)
    {
      printf(".");
      fflush(stdout);
    }

    double wt = 0.0;
    unsigned long long cc = 0;
    int res = f(params, &wt, &cc);
    if (res)
      m4ri_die("benchmark function failed with exit code: %d\n", res);

    if (bench_dump)
    {
      printf("%f\n", wt);
      fflush(stdout);
    }

    cc_sum += cc;

    vector_pushback(wt_data, wt);

    if (n >= bench_minimum && normal_calculate(wt_data, &stats) == 0)
    {
      double standard_error = stats.sigma / sqrt(stats.size);
      double critical_value = t_table(bench_confidence_index, stats.size - 1);
      // Stop when the real mean lays with CONFIDENCE in the range [mean * (1 - bench_accuracy), mean * (1 + bench_accuracy)].
      // or when we're already running bench_maxtime seconds.
      if (standard_error * critical_value / stats.mean <= bench_accuracy ||
	  walltime(start_walltime) > bench_maxtime)
        break;
    }
  }

  *wtp = stats.mean;
  *ccp = cc_sum / stats.size;

  if (!bench_quiet)
  {
    if (!bench_quiet && !bench_dump)
      printf("\n");
    printf("Total running time: %6.3lf seconds.\n", walltime(start_walltime));
    printf("Sample size: %d; mean: %f; standard deviation: %f\n", stats.size, stats.mean, stats.sigma);
    double standard_error = stats.sigma / sqrt(stats.size);
    double critical_value = t_table(bench_confidence_index, stats.size - 1);
    double accuracy = standard_error * critical_value;
    printf("%2.0lf%% confidence interval: +/- %lf (%.1lf%%): [%lf..%lf]\n", (CONFIDENCE * 100), accuracy, accuracy / stats.mean * 100, stats.mean - accuracy, stats.mean + accuracy);
  }

  vector_destruct(wt_data);
}

/*
 * Random number generator
 */

static uint64_t bench_random_M;
static uint64_t bench_random_modulo;

uint64_t bench_random_init(uint64_t modulo)
{
  // Set bench_random_M to the largest multiple of modulo, minus one, that fits in an uint64_t.
  // A modulo of zero is interpreted as 2^64, and thus returns 0xffffffffffffffff.
  bench_random_M = modulo ? -modulo / modulo * modulo - 1 : -1;
  bench_random_M += modulo;
  bench_random_modulo = modulo;
}

// Returns a uniformly distributed random number in the range [0, bench_random_modulo>.
uint64_t bench_random()
{
  for(;;)
  {
    if (sizeof(long) == sizeof(uint64_t))
    {
      union { uint64_t R; long L; } x;
      x.L = random();
      if (x.R <= bench_random_M)
	return x.R % bench_random_modulo;
    }
    else
    {
      union { uint64_t R; long L1; long L2; } x;
      x.L1 = random();
      x.L2 = random();
      if (x.R <= bench_random_M)
	return x.R % bench_random_modulo;
    }
  }
}

// The same as m4ri_random_word. Duplicated here because it's
// not available in older revisions that we want to benchmark against.
word bench_random_word() {
  // random() only returns 31 bits, so we need three calls.
  word a0 = random();
  word a1 = random();
  word a2 = random();
  word v = a0 ^ (a1 << 24) ^ a2 << 48;
#ifdef BENCH_RANDOM_REVERSE
  v = ((v >>  1) & 0x5555555555555555ULL) | ((v & 0x5555555555555555ULL) << 1);
  v = ((v >>  2) & 0x3333333333333333ULL) | ((v & 0x3333333333333333ULL) << 2);
  v = ((v >>  4) & 0x0F0F0F0F0F0F0F0FULL) | ((v & 0x0F0F0F0F0F0F0F0FULL) << 4);
  v = ((v >>  8) & 0x00FF00FF00FF00FFULL) | ((v & 0x00FF00FF00FF00FFULL) << 8);
  v = ((v >> 16) & 0x0000FFFF0000FFFFULL) | ((v & 0x0000FFFF0000FFFFULL) << 16);
  v =  (v >> 32)                          |  (v                          << 32);
#endif
  return v;
}

// Needed for mzd_t.
#include "packedmatrix.h"

// The same as m4ri_randomize. Duplicated here because it's
// not available in older revisions that we want to benchmark against.
void bench_randomize(mzd_t *A) {
  wi_t const width = A->width - 1;
  int const offset = A->offset;
  if(offset) {
    if(width == 0) {
      word const mask = MIDDLE_BITMASK(A->ncols, offset);
      for(rci_t i = 0; i < A->nrows; ++i)
	A->rows[i][0] ^= (A->rows[i][0] ^ (bench_random_word() << offset)) & mask;
    } else {
      word const mask_begin = RIGHT_BITMASK(RADIX - offset);
      word const mask_end = LEFT_BITMASK((A->ncols + offset) % RADIX);
      int const need_last_bits = ((ONE << offset) & mask_end) != 0;
      for(rci_t i = 0; i < A->nrows; ++i) {
	word prev_random_word;
	word random_word = bench_random_word();
	A->rows[i][0] ^= (A->rows[i][0] ^ (random_word << offset)) & mask_begin;
	for(wi_t j = 1; j < width; ++j) {
	  prev_random_word = random_word;
	  random_word = bench_random_word();
	  A->rows[i][j] = (random_word << offset) | (prev_random_word >> (RADIX - offset));
	}
	prev_random_word = random_word;
	random_word = 0;
	if (need_last_bits)
	  random_word = bench_random_word();
	A->rows[i][width] ^= (A->rows[i][width] ^ ((random_word << offset) | (prev_random_word >> (RADIX - offset)))) & mask_end;
      }
    }
  } else {
    word const mask_end = LEFT_BITMASK(A->ncols % RADIX);
    for(rci_t i = 0; i < A->nrows; ++i) {
      for(wi_t j = 0; j < width; ++j)
	A->rows[i][j] = bench_random_word();
      A->rows[i][width] ^= (A->rows[i][width] ^ bench_random_word()) & mask_end;
    }
  }
}

