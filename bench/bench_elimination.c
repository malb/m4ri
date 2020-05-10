#include <stdlib.h>

#include "benchmarking.h"
#include "cpucycles.h"
#include <m4ri/config.h>
#include <m4ri/m4ri.h>

#ifdef HAVE_LIBPAPI
#define _GNU_SOURCE
#include <errno.h>
#include <papi.h>
#include <sys/types.h>  // papi.h needs caddr_t
#endif

struct elim_params {
  rci_t m;
  rci_t n;
  rci_t r;
  char const *algorithm;
};

static unsigned long long loop_calibration[32];

int run_nothing(void *_p, unsigned long long *data, int *data_len) {
  struct elim_params *p = (struct elim_params *)_p;

  mzd_t *A = mzd_init(p->m, p->n);

  if (p->r != 0) {
    mzd_t *L, *U;
    L = mzd_init(p->m, p->m);
    U = mzd_init(p->m, p->n);
    mzd_randomize(U);
    mzd_randomize(L);
    for (rci_t i = 0; i < p->m; ++i) {

      for (rci_t j = i + 1; j < p->m; j += m4ri_radix) {
        int const length = MIN(m4ri_radix, p->m - j);
        mzd_clear_bits(L, i, j, length);
      }
      mzd_write_bit(L, i, i, 1);

      for (rci_t j = 0; j < i; j += m4ri_radix) {
        int const length = MIN(m4ri_radix, i - j);
        mzd_clear_bits(U, i, j, length);
      }
      if (i < p->r) {
        mzd_write_bit(U, i, i, 1);
      } else {
        for (rci_t j = i; j < p->n; j += m4ri_radix) {
          int const length = MIN(m4ri_radix, p->n - i);
          mzd_clear_bits(U, i, j, length);
        }
      }
    }
    mzd_mul(A, L, U, 0);
    mzd_free(L);
    mzd_free(U);
  } else {
    mzd_randomize(A);
  }

#ifndef HAVE_LIBPAPI
  *data_len = 2;
#else
  *data_len             = MIN(papi_array_len + 1, *data_len);
#endif
  int papi_res;

#ifndef HAVE_LIBPAPI
  data[0] = walltime(0);
  data[1] = cpucycles();
#else
  int array_len         = *data_len - 1;
  unsigned long long t0 = PAPI_get_virt_usec();
  papi_res              = PAPI_start_counters((int *)papi_events, array_len);
  if (papi_res) m4ri_die("");
#endif

#ifndef HAVE_LIBPAPI
  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);
#else
  PAPI_stop_counters((long long *)&data[1], array_len);
  t0      = PAPI_get_virt_usec() - t0;
  data[0] = t0;
  for (int nv = 0; nv <= array_len; ++nv) {
    if (data[nv] < loop_calibration[nv]) loop_calibration[nv] = data[nv];
  }
#endif

  mzd_free(A);

  return (0);
}

int run(void *_p, unsigned long long *data, int *data_len) {
  struct elim_params *p = (struct elim_params *)_p;
#ifndef HAVE_LIBPAPI
  *data_len = 2;
#else
  *data_len             = MIN(papi_array_len + 1, *data_len);
#endif
  int papi_res;

  mzd_t *A = mzd_init(p->m, p->n);

  if (p->r != 0) {
    mzd_t *L, *U;
    L = mzd_init(p->m, p->m);
    U = mzd_init(p->m, p->n);
    mzd_randomize(U);
    mzd_randomize(L);
    for (rci_t i = 0; i < p->m; ++i) {

      for (rci_t j = i + 1; j < p->m; j += m4ri_radix) {
        int const length = MIN(m4ri_radix, p->m - j);
        mzd_clear_bits(L, i, j, length);
      }
      mzd_write_bit(L, i, i, 1);

      for (rci_t j = 0; j < i; j += m4ri_radix) {
        int const length = MIN(m4ri_radix, i - j);
        mzd_clear_bits(U, i, j, length);
      }
      if (i < p->r) {
        mzd_write_bit(U, i, i, 1);
      } else {
        for (rci_t j = i; j < p->n; j += m4ri_radix) {
          int const length = MIN(m4ri_radix, p->n - i);
          mzd_clear_bits(U, i, j, length);
        }
      }
    }
    mzd_mul(A, L, U, 0);
    mzd_free(L);
    mzd_free(U);
  } else {
    mzd_randomize(A);
  }

#ifndef HAVE_LIBPAPI
  data[0] = walltime(0);
  data[1] = cpucycles();
#else
  int array_len         = *data_len - 1;
  unsigned long long t0 = PAPI_get_virt_usec();
  papi_res              = PAPI_start_counters((int *)papi_events, array_len);
  if (papi_res) m4ri_die("");
#endif
  if (strcmp(p->algorithm, "m4ri") == 0)
    p->r = mzd_echelonize_m4ri(A, 1, 0);
  else if (strcmp(p->algorithm, "pluq") == 0)
    p->r = mzd_echelonize_pluq(A, 1);
  else if (strcmp(p->algorithm, "mmpf") == 0)
    p->r = _mzd_pluq_russian(A, mzp_init(A->nrows), mzp_init(A->ncols), 0);
  else if (strcmp(p->algorithm, "naive") == 0)
    p->r = mzd_echelonize_naive(A, 1);
#ifndef HAVE_LIBPAPI
  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);
#else
  PAPI_stop_counters((long long *)&data[1], array_len);
  t0      = PAPI_get_virt_usec() - t0;
  data[0] = t0;
  for (int nv = 0; nv <= array_len; ++nv) { data[nv] -= loop_calibration[nv]; }
#endif
  mzd_free(A);
  return 0;
}

void print_help_and_exit() {
  printf("Parameters m(, n, alg, r) expected.\n");
  printf(" m   -- integer > 0\n");
  printf(" n   -- integer > 0\n");
  printf(" alg -- 'm4ri', 'pluq', 'mmpf' or 'naive' (default: 'pluq')\n");
  printf(" r   -- target rank >= 0, if 0 then mzd_randomize() is called (default: MIN(m,n))\n");
  printf("\n");
  bench_print_global_options(stderr);
  m4ri_die("");
}

int main(int argc, char **argv) {
  int opts = global_options(&argc, &argv);
  int data_len;

#ifdef HAVE_LIBPAPI
  int papi_counters = PAPI_num_counters();
  if (papi_counters < papi_array_len) {
    fprintf(stderr, "%s: Warning: there are only %d hardware counters available!\n", progname,
            papi_counters);
    papi_array_len = papi_counters;
  }
  if (papi_test(papi_events, papi_array_len)) exit(1);

  for (int nv = 0; nv <= papi_array_len; ++nv) loop_calibration[nv] = 100000000;

  data_len = papi_array_len + 1;
#else
  data_len = 2;
#endif
  if (opts < 0 || argc < 2 || argc > 5) { print_help_and_exit(); }

  struct elim_params params;
  params.m = atoi(argv[1]);
  if (argc >= 3)
    params.n = atoi(argv[2]);
  else
    params.n = params.m;

  if (argc >= 4)
    params.algorithm = argv[3];
  else
    params.algorithm = "pluq";
  if (argc >= 5)
    params.r = atoi(argv[4]);
  else
    params.r = params.m;

  srandom(17);
  unsigned long long data[16];

  for (int i = 0; i < 4; ++i) run_nothing((void *)&params, data, &data_len);

  run_bench(run, (void *)&params, data, data_len);

  double cc_per_op =
      ((double)data[1]) / ((double)params.m * (double)params.n * powl((double)params.r, 0.807));

  printf("m: %5d, n: %5d, last r: %5d, cpu cycles: %12llu, cc/(mnr^0.807): %.5lf, ", params.m,
         params.n, params.r, data[1], cc_per_op);
  print_wall_time(data[0] / 1000000.0);
  printf(", ");
  print_cpu_time(data[1] / (double)cpucycles_persecond());
  printf("\n");
#ifdef HAVE_LIBPAPI
  for (int n = 1; n < data_len; ++n) {
    double tmp = ((double)data[n]) / powl((double)params.n, 2.807);
    printf("%20s (%20llu) per bit (divided by n^2.807): %15.5f\n",
           papi_event_name(papi_events[n - 1]), data[n], tmp);
  }
#endif
}
