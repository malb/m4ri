/*
 * benchmarking.c
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
 * ./bench_elimination -s 0 -m 4 -c 90 -a 0.005 -d -t 30 -n 1000 1000 1000
 *
 * would run at most 30 seconds (-t) or 1000 times (-n), whichever comes
 * first, or stop after the real average of wall time (-s 0) falls with 90%
 * certainty (-c) in a range that is +/- 0.005 times the observed mean (-a: accuracry),
 * but no sooner than that at least 4 (-m: minimum) measurements have been
 * done. It would also print (-d: dump) each measurement (0:microseconds 1:cpuclocks).
 *
 * Example output.
 *
 * 2416 6441500
 * 2376 6335490
 * 2360 6294450
 * 2361 6295280
 * 2371 6321440
 * 2350 6266740
 * 2362 6298700
 * 2386 6362520
 * 2344 6249890
 * 2347 6260450
 * 2346 6254590
 * Total running time:  0.103 seconds.
 * Virtual time (s): Sample size: 11; mean: 0.002365; standard deviation: 0.000021
 * Virtual time (s): 90% confidence interval: +/- 0.000012 (0.5%): [0.002354..0.002377]
 *
 * The last three lines can be suppressed by passing the option -q (quiet).
 */

#include <m4ri/config.h>

#ifdef HAVE_LIBPAPI
#define _GNU_SOURCE
#include <errno.h>
#include <papi.h>
#include <sys/types.h>  // papi.h needs caddr_t
#endif

#include "benchmarking.h"
#include <ctype.h>
#include <m4ri/misc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

enum { C80, C90, C95, C98, C99 };

/*
 * Command line option decoding
 */

int bench_quiet                  = 0;     // Set if -q is used.
int bench_dump                   = 0;     // Set if -d is used.
int bench_minimum                = 2;     // Minimum number of measurements. Set with -m <minimum>.
int bench_maximum                = 1000;  // Maximum number of measurements. Set with -n <maximum>.
unsigned long long bench_maxtime = 60000000;  // Maximum number of microseconds to run. Set with -t
                                              // <maxtime>, in seconds (floating point).
double bench_accuracy =
    0.01;  // The +/- range (where 1.0 is 100%) within that we want the real population mean to be
           // with the given confidence. Set with -a <bench_accuracy>
int bench_confidence_index =
    C99;  // The confidence that the real mean is within the given (or found) range.
int bench_stats =
    1;  // The counter used for statistics (0 = realtime, 1 = cpuclocks). Set with -s <counter>.
int bench_dump_counter = -1;  // The counter to dump (see bench_stats). Set with -d <counter>. If
                              // not given all counters are dumped.
char const *progname;         // Set to argv[0].

/*
 * Command line option used by bench_packedmatrix.c
 */

uint64_t bench_count =
    0;  // Can be set by -x <count>, otherwise a reasonable default is being used.

#ifdef HAVE_LIBPAPI

int bench_disregard_L2_misses = 0;  // Set if -2 is used.

/*
 * PAPI events being counted.
 */

int papi_events[32] = {
    PAPI_TOT_CYC, /* Total cycles. This must always be the first entry. */
};

int papi_array_len = 1;
int bench_PAPI_L2_TCM_index;

char *papi_event_name(int event) {
  // PAPI needs to be initialized before calling PAPI_event_code_to_name.
  if (PAPI_is_initialized() == PAPI_NOT_INITED) {
    int res = PAPI_library_init(PAPI_VER_CURRENT);
    if (res != PAPI_OK && res != PAPI_VER_CURRENT) {
      fprintf(stderr, "%s: PAPI_library_init: error code %d %s\n", progname, res,
              PAPI_strerror(res));
      m4ri_die("PAPI failed to initialize.\n");
    }
  }

  static char buf[PAPI_MAX_STR_LEN];
  int res = PAPI_event_code_to_name(event, buf);
  if (res) snprintf(buf, PAPI_MAX_STR_LEN, "<unknown PAPI event code %d>", event);
  return buf;
}

int papi_add_event(char const *event_name) {
  // PAPI needs to be initialized before calling PAPI_event_name_to_code.
  if (PAPI_is_initialized() == PAPI_NOT_INITED) {
    int res = PAPI_library_init(PAPI_VER_CURRENT);
    if (res != PAPI_OK && res != PAPI_VER_CURRENT) {
      fprintf(stderr, "%s: PAPI_library_init: error code %d %s\n", progname, res,
              PAPI_strerror(res));
      m4ri_die("PAPI failed to initialize.\n");
    }
  }

  int event;
  int res = PAPI_event_name_to_code((char *)event_name, &event);

  if (res != PAPI_OK) {
    if (res == PAPI_ENOEVNT)
      fprintf(stderr, "%s: %s: No such event.\n", progname, event_name);
    else
      fprintf(stderr, "%s: PAPI_event_name_to_code(\"%s\"): %s\n", progname, event_name,
              PAPI_strerror(res));
    return res;
  }

  int found = 0;
  for (int nv = 0; nv < papi_array_len; ++nv) {
    if (papi_events[nv] == event) {
      found = 1;
      break;
    }
  }

  if (!found) papi_events[papi_array_len++] = event;

  return 0;
}

void papi_add_events(char *event_names) {
  char *tmpptr;
  char *name = strtok_r(event_names, ", ", &tmpptr);
  while (name) {
    papi_add_event(name);
    name = strtok_r(NULL, ", ", &tmpptr);
  }
}

#endif  // HAVE_LIBPAPI

int global_options(int *argcp, char ***argvp) {
  int result = 0;
  progname   = (*argvp)[0];
  while ((*argcp) > 1) {
    if ((*argvp)[1][0] != '-' || (*argvp)[1][1] == '\0' || (*argvp)[1][2] != '\0') return result;
    switch ((*argvp)[1][1]) {
    case 'd':
      bench_dump = 1;
      if (isdigit((*argvp)[2][0])) {
        ++*argvp;
        --*argcp;
        bench_dump_counter = atoi((*argvp)[1]);
      }
      break;
    case 'q': bench_quiet = 1; break;
#ifdef HAVE_LIBPAPI
    case '2': {
      bench_disregard_L2_misses = 1;
      if (papi_add_event("PAPI_L2_TCM")) {
        fprintf(stderr,
                "%s: Ignoring -2: Level 2 cache misses cannot be detected with the current set of "
                "PAPI events (-p).\n",
                progname);
        bench_disregard_L2_misses = 0;
      }
      for (int nv = 0; nv < papi_array_len; ++nv) {
        if (papi_events[nv] == PAPI_L2_TCM) {
          bench_PAPI_L2_TCM_index = nv + 1;  // +1 for in data[] inserted virtual time at index 0.
          break;
        }
      }
      break;
    }
    case 'p': {
      ++*argvp;
      --*argcp;
      papi_add_events((*argvp)[1]);
      break;
    }
#endif
    case 'm':
      ++*argvp;
      --*argcp;
      bench_minimum = atoi((*argvp)[1]);
      break;
    case 'n':
      ++*argvp;
      --*argcp;
      bench_maximum = atoi((*argvp)[1]);
      if (bench_maximum < bench_minimum) bench_minimum = bench_maximum;
      break;
    case 't':
      ++*argvp;
      --*argcp;
      bench_maxtime = 1000000 * strtod((*argvp)[1], NULL);
      break;
    case 'a':
      ++*argvp;
      --*argcp;
      bench_accuracy = strtod((*argvp)[1], NULL);
      break;
    case 'c': {
      ++*argvp;
      --*argcp;
      int confidence = atoi((*argvp)[1]);
      switch (confidence) {
      case 80: bench_confidence_index = C80; break;
      case 90: bench_confidence_index = C90; break;
      case 95: bench_confidence_index = C95; break;
      case 98: bench_confidence_index = C98; break;
      case 99: bench_confidence_index = C99; break;
      default:
        m4ri_die("The only possible confidence percentages are 80, 90, 95, 98 and 99%\n");
        break;
      }
      break;
    }
    case 'x':
      ++*argvp;
      --*argcp;
      bench_count = atoll((*argvp)[1]);
      break;
    case 's':
      ++*argvp;
      --*argcp;
      bench_stats = atoi((*argvp)[1]);
      break;
    default: return -1;
    }
    ++result;
    ++*argvp;
    --*argcp;
  }
  return result;
}

void bench_print_global_options(FILE *out) {
  fprintf(out, "OPTIONS\n");
  fprintf(out, "  -m <minimum>      Do at least <minimum> number of measurements. Default 2.\n");
  fprintf(out, "  -n <maximum>      Do at most <maximum> number of measurements. Default 1000.\n");
  fprintf(out, "  -t <max-time>     Stop after <max-time> seconds. Default 60.0 seconds.\n");
  fprintf(out,
          "  -a <accuracy>     Stop after <accuracy> has been reached. Default 0.01 (= 1%%).\n");
  fprintf(out, "  -c <confidence>   Stop when accuracy has been reached with this confidence. "
               "Default 99 (%%).\n");
  fprintf(out, "  -s <counter>      Counter to perform statistic over (0: realtime, 1: cpuclocks. "
               "Default: 1).\n");
  fprintf(out, "  -x <loop-count>   Call function <loop-count> times in the inner most loop (calls "
               "per measurement).\n");
  fprintf(out, "  -d [<counter>]    Dump measurements. Dump all or only <counter> when given.\n");
  fprintf(out, "  -q                Quiet. Suppress printing of statistics.\n");
#ifdef HAVE_LIBPAPI
  fprintf(out, "  -2                Disregard measurements with any level 2 cache misses.\n");
  fprintf(out, "  -p <PAPI-event>[,<PAPI-event>,...]\n");
  fprintf(out, "                    Count and report the given events. The list is comma or space "
               "separated,\n");
  fprintf(out, "                    for example -p \"PAPI_TOT_INS PAPI_L1_DCM\".\n");
  fprintf(out, "                    Run `papi_event_chooser PRESET PAPI_TOT_CYC [PAPI_*]` for more "
               "events.\n");
#endif
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
  double *data;
};

typedef struct vector_st *vector;

vector vector_create(size_t s) {
  vector v      = (vector)malloc(sizeof(struct vector_st));
  v->alloc_size = s;
  v->data       = s ? (double *)malloc(sizeof(double) * s) : NULL;
  v->size       = 0;
  return v;
}

void vector_destruct(vector v) {
  free(v->data);
  free(v);
}

void vector_resize(vector v, size_t s) {
  v->data       = (double *)realloc(v->data, sizeof(double) * s);
  v->alloc_size = s;
  if (v->size > v->alloc_size) v->size = v->alloc_size;
}

static inline size_t vector_size(vector v) { return v->size; }

void vector_pushback(vector v, double d) {
  if (++(v->size) > v->alloc_size) vector_resize(v, v->alloc_size * 2);
  v->data[v->size - 1] = d;
}

static inline double vector_get(vector v, int index) { return v->data[index]; }

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

int normal_calculate(vector v, normal *dist, double multiplier) {
  dist->size = vector_size(v);

  if (dist->size < 2) {
    dist->mean  = vector_get(v, 0) * multiplier;
    dist->sigma = 0.0;
    return 0;
  }
  // Calculate the sum of all data.
  double sum = 0;
  for (int i = 0; i < dist->size; ++i) sum += vector_get(v, i) * multiplier;
  dist->mean = sum / dist->size;

  // Calculate the sum of the square of all differences with mean.
  sum = 0;
  for (int i = 0; i < dist->size; ++i) {
    double delta = vector_get(v, i) * multiplier - dist->mean;
    sum += delta * delta;
  }
  dist->sigma = sqrt(sum / (dist->size - 1));

  return 0;
}

/*
 * T-Table
 */

static float student_t[5][34] = {
    {3.078, 1.886, 1.638, 1.533, 1.476, 1.440, 1.415, 1.397, 1.383, 1.372, 1.363, 1.356,
     1.350, 1.345, 1.341, 1.337, 1.333, 1.330, 1.328, 1.325, 1.323, 1.321, 1.319, 1.318,
     1.316, 1.315, 1.314, 1.313, 1.311, 1.310, 1.303, 1.296, 1.289, 1.282},
    {6.314, 2.920, 2.353, 2.132, 2.015, 1.943, 1.895, 1.860, 1.833, 1.812, 1.796, 1.782,
     1.771, 1.761, 1.753, 1.746, 1.740, 1.734, 1.729, 1.725, 1.721, 1.717, 1.714, 1.711,
     1.708, 1.706, 1.703, 1.701, 1.699, 1.697, 1.684, 1.671, 1.658, 1.645},
    {12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179,
     2.160,  2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064,
     2.060,  2.056, 2.052, 2.048, 2.045, 2.042, 2.021, 2.000, 1.980, 1.960},
    {31.821, 6.965, 4.541, 3.747, 3.365, 3.143, 2.998, 2.896, 2.821, 2.764, 2.718, 2.681,
     2.650,  2.624, 2.602, 2.583, 2.567, 2.552, 2.539, 2.528, 2.518, 2.508, 2.500, 2.492,
     2.485,  2.479, 2.473, 2.467, 2.462, 2.457, 2.423, 2.390, 2.358, 2.326},
    {63.657, 9.925, 5.841, 4.604, 4.032, 3.707, 3.499, 3.355, 3.250, 3.169, 3.106, 3.055,
     3.012,  2.977, 2.947, 2.921, 2.898, 2.878, 2.861, 2.845, 2.831, 2.819, 2.807, 2.797,
     2.787,  2.779, 2.771, 2.763, 2.756, 2.750, 2.704, 2.660, 2.617, 2.576}};

static float student_t_certainty[5] = {0.2, 0.1, 0.05, 0.02, 0.01};  // Two-tails.

static float t_table(int confidence_index, int freedoms) {
  if (freedoms <= 30) return student_t[confidence_index][freedoms - 1];
  double a, b, y1, y2, y3;
  long x1, x2;
  long x3 = 0;
  int i;
  if (freedoms <= 60) {
    i  = 29;
    x1 = 30;
    x2 = 40;
    x3 = 60;
  } else if (freedoms <= 120) {
    i  = 30;
    x1 = 40;
    x2 = 60;
    x3 = 120;
  } else {
    i  = 31;
    x1 = 60;
    x2 = 120;
    /* x3 = infinity */
  }
  y1 = student_t[confidence_index][i];
  y2 = student_t[confidence_index][i + 1];
  y3 = student_t[confidence_index][i + 2];
  if (freedoms <= 120) {
    double c, d;
    d = (x1 * x1 * (x3 - x2) + x2 * x2 * (x1 - x3) + x3 * x3 * (x2 - x1));
    a = -(x1 * (y3 - y2) + x2 * (y1 - y3) + x3 * (y2 - y1)) / d;
    b = (x1 * x1 * (y3 - y2) + x2 * x2 * (y1 - y3) + x3 * x3 * (y2 - y1)) / d;
    c = y2 - a * x2 * x2 - b * x2;
    return (a * freedoms * freedoms + b * freedoms + c);
  }
  double ln1, ln2;
  ln1 = log(y2 - y3);
  ln2 = log(y1 - y3);
  a   = -(ln1 - ln2) / (x1 - x2);
  b   = (x1 * ln1 - x2 * ln2) / (x1 - x2);
  return (y3 + exp(a * freedoms + b));
}

/*
 * walltime
 */

unsigned long long walltime(unsigned long long t0) {
  static time_t base_sec;
  struct timeval tp;
  gettimeofday(&tp, NULL);
  if (__M4RI_UNLIKELY(base_sec == 0)) base_sec = tp.tv_sec;
  return (tp.tv_sec - base_sec) * 1000000 + tp.tv_usec - t0;
}

/*
 * Printing doubles.
 */

int bench_precision(double sigma) {
  if (sigma < 1E-10) return 12;
  int log_sigma = log10(sigma);
  if (log_sigma >= 2) return 0;
  return 2 - log_sigma;
}

void print_double(double d, int precision) {
  switch (precision) {
  case 0: printf("%.0f", d); break;
  case 1: printf("%.1f", d); break;
  case 2: printf("%.2f", d); break;
  case 3: printf("%.3f", d); break;
  case 4: printf("%.4f", d); break;
  case 5: printf("%.5f", d); break;
  case 6: printf("%.6f", d); break;
  case 7: printf("%.7f", d); break;
  case 8: printf("%.8f", d); break;
  case 9: printf("%.9f", d); break;
  case 10: printf("%.10f", d); break;
  case 11: printf("%.11f", d); break;
  case 12: printf("%.12f", d); break;
  }
}

/*
 * run_bench
 *
 * Benchmark main loop.
 */

int run_bench(int (*f)(void *params, unsigned long long *data, int *data_len), void *params,
              unsigned long long *data, int data_len) {
  double const CONFIDENCE = 1.0 - student_t_certainty[bench_confidence_index];
  unsigned long long data_sum[32];
  memset(data_sum, 0, sizeof(data_sum));
  data_len          = MIN(data_len, sizeof(data_sum) / sizeof(unsigned long long));
  vector stats_data = vector_create(128);
  normal stats;
#ifdef HAVE_LIBPAPI
  int total_calls = 0;
#endif
  if (!bench_count) bench_count = 1;
  unsigned long long start_walltime = walltime(0);
  for (int n = 1; n <= bench_maximum; ++n) {
    if (!bench_quiet && !bench_dump) {
      printf(".");
      fflush(stdout);
    }

    do {
      int res = f(params, data, &data_len);
      if (res < 0) m4ri_die("benchmark function failed with exit code: %d\n", res);
#ifdef HAVE_LIBPAPI
      ++total_calls;
#endif
    }
#ifdef HAVE_LIBPAPI
    while (bench_disregard_L2_misses && data[bench_PAPI_L2_TCM_index]);
#else
    while (0);
#endif

    if (bench_dump) {
      if (bench_dump_counter >= 0 && bench_dump_counter < data_len)
        printf("%llu", data[bench_dump_counter]);
      else {
        printf("%llu", data[0]);
        for (int nv = 1; nv < data_len; ++nv) printf(" %llu", data[nv]);
      }
      printf("\n");
      fflush(stdout);
    }

    vector_pushback(stats_data, data[bench_stats]);
    for (int nv = 0; nv < data_len; ++nv) data_sum[nv] += data[nv];

    if (n >= bench_minimum &&
        normal_calculate(stats_data, &stats, (bench_stats == 0) ? 0.000001 : (1.0 / bench_count)) ==
            0) {
      double standard_error = stats.sigma / sqrt(stats.size);
      double critical_value = t_table(bench_confidence_index, stats.size - 1);
      // Stop when the real mean lays with CONFIDENCE in the range [mean * (1 - bench_accuracy),
      // mean * (1 + bench_accuracy)]. or when we're already running bench_maxtime seconds.
      if (standard_error * critical_value / stats.mean <= bench_accuracy ||
          walltime(start_walltime) > bench_maxtime)
        break;
    }
  }

  for (int nv = 0; nv < data_len; ++nv) data[nv] = (data_sum[nv] + stats.size / 2) / stats.size;

  if (!bench_quiet) {
    if (!bench_quiet && !bench_dump) printf("\n");
    printf("Total running time: %6.3f seconds.\n", walltime(start_walltime) / 1000000.0);
#ifdef HAVE_LIBPAPI
    if (bench_disregard_L2_misses)
      printf("Samples disregarded because of level 2 cache misses: %d\n", total_calls - stats.size);
#endif
    int precision = bench_precision(stats.sigma);
#ifdef HAVE_LIBPAPI
    if (bench_stats)
      printf("%s: ", papi_event_name(papi_events[bench_stats - 1]));
    else
      printf("Virtual time (s): ");
#endif
    printf("Sample size: %d; mean: ", stats.size);
    print_double(stats.mean, precision);
    printf("; standard deviation: ");
    print_double(stats.sigma, precision);
    printf("\n");
#ifdef HAVE_LIBPAPI
    if (bench_stats)
      printf("%s: ", papi_event_name(papi_events[bench_stats - 1]));
    else
      printf("Virtual time (s): ");
#endif
    double standard_error = stats.sigma / sqrt(stats.size);
    double critical_value = t_table(bench_confidence_index, stats.size - 1);
    double accuracy       = standard_error * critical_value;
    printf("%2.0f%% confidence interval: +/- ", CONFIDENCE * 100);
    print_double(accuracy, precision);
    printf(" (%.1f%%): [", accuracy / stats.mean * 100);
    print_double(stats.mean - accuracy, precision);
    printf("..");
    print_double(stats.mean + accuracy, precision);
    printf("]\n");
  }

  vector_destruct(stats_data);

  return data_len;
}

/*
 * Randomize
 */

// The same as m4ri_random_word. Duplicated here because it's
// not available in older revisions that we want to benchmark against.
word bench_random_word() {
  // random() only returns 31 bits, so we need three calls.
  word a0 = random();
  word a1 = random();
  word a2 = random();
  word v  = a0 ^ (a1 << 24) ^ a2 << 48;
#ifdef BENCH_RANDOM_REVERSE
  v = ((v >> 1) & 0x5555555555555555ULL) | ((v & 0x5555555555555555ULL) << 1);
  v = ((v >> 2) & 0x3333333333333333ULL) | ((v & 0x3333333333333333ULL) << 2);
  v = ((v >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((v & 0x0F0F0F0F0F0F0F0FULL) << 4);
  v = ((v >> 8) & 0x00FF00FF00FF00FFULL) | ((v & 0x00FF00FF00FF00FFULL) << 8);
  v = ((v >> 16) & 0x0000FFFF0000FFFFULL) | ((v & 0x0000FFFF0000FFFFULL) << 16);
  v = (v >> 32) | (v << 32);
#endif
  return v;
}

// Needed for mzd_t.
#include <m4ri/mzd.h>

// The same as m4ri_randomize. Duplicated here because it's
// not available in older revisions that we want to benchmark against.
void bench_randomize(mzd_t *A) {
  wi_t const width    = A->width - 1;
  int const offset    = 0;
  word const mask_end = __M4RI_LEFT_BITMASK(A->ncols % m4ri_radix);
  for (rci_t i = 0; i < A->nrows; ++i) {
    for (wi_t j = 0; j < width; ++j) A->rows[i][j] = bench_random_word();
    A->rows[i][width] ^= (A->rows[i][width] ^ bench_random_word()) & mask_end;
  }
}

/*
 * Random number generator
 */

static uint64_t bench_random_M;
static uint64_t bench_random_modulo;

void bench_random_init(uint64_t modulo) {
  // Set bench_random_M to the largest multiple of modulo, minus one, that fits in an uint64_t.
  // A modulo of zero is interpreted as 2^64, and thus returns 0xffffffffffffffff.
  bench_random_M = modulo ? -modulo / modulo * modulo - 1 : -1;
  bench_random_M += modulo;
  bench_random_modulo = modulo;
}

// Returns a uniformly distributed random number in the range [0, bench_random_modulo>.
uint64_t bench_random() {
  for (;;) {
    word R = bench_random_word();
    if (R <= bench_random_M) return R % bench_random_modulo;
  }
}

void print_wall_time(double seconds) {
  if (seconds >= 0.01)
    printf("wall time: %10.5f s", seconds);
  else if (seconds >= 0.00001)
    printf("wall time: %10.5f ms", 1000.0 * seconds);
  else
    printf("wall time: %10.5f us", 1000000.0 * seconds);
}

void print_cpu_time(double seconds) {
  if (seconds >= 0.01)
    printf("cpu time: %10.5f s", seconds);
  else if (seconds >= 0.00001)
    printf("cpu time: %10.5f ms", 1000.0 * seconds);
  else
    printf("cpu time: %10.5f us", 1000000.0 * seconds);
}

#ifdef HAVE_LIBPAPI
int papi_test(int *papi_events, int papi_array_len) {
  int res = PAPI_start_counters(papi_events, papi_array_len);
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
  return res;
}
#endif
