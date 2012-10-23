#include <stdlib.h>

#include "config.h"
#include "cpucycles.h"
#include "m4ri.h"
#include "benchmarking.h"

#ifdef HAVE_LIBPAPI
#define _GNU_SOURCE
#include <sys/types.h>          // papi.h needs caddr_t
#include <papi.h>
#include <errno.h>
#endif

struct mul_params {
  rci_t n;
  int cutoff;
  uint64_t count;
};

static unsigned long long loop_calibration[32];

int run_nothing(void *_p, unsigned long long *data, int *data_len) {
  struct mul_params *p = (struct mul_params *)_p;

  mzd_t *A = mzd_init(p->n, p->n);
  mzd_t *B = mzd_init(p->n, p->n);

  mzd_randomize(A);
  mzd_randomize(B);

#ifndef HAVE_LIBPAPI
  *data_len = 2;
#else
  *data_len = MIN(papi_array_len + 1, *data_len);
#endif
  int papi_res;

#ifndef HAVE_LIBPAPI
  data[0] = walltime(0);
  data[1] = cpucycles();
#else
  int array_len = *data_len - 1;
  unsigned long long t0 = PAPI_get_virt_usec();
  papi_res = PAPI_start_counters((int*)papi_events, array_len);
  if(papi_res)
    m4ri_die("");
#endif

#ifndef HAVE_LIBPAPI
  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);
#else
  PAPI_stop_counters((long long*)&data[1], array_len);
  t0 = PAPI_get_virt_usec() - t0;
  data[0] = t0;
  for (int nv = 0; nv <= array_len; ++nv) {
    if (data[nv] < loop_calibration[nv])
      loop_calibration[nv] = data[nv];
  }
#endif

  mzd_free(A);
  mzd_free(B);

  return (0);
}


int run(void *_p, unsigned long long *data, int *data_len) {
  struct mul_params *p = (struct mul_params *)_p;

#ifndef HAVE_LIBPAPI
  *data_len = 2;
#else
  *data_len = MIN(papi_array_len + 1, *data_len);
#endif
  int papi_res;

  mzd_t *A = mzd_init(p->n, p->n);
  mzd_t *B = mzd_init(p->n, p->n);

  mzd_randomize(A);
  mzd_randomize(B);

#ifndef HAVE_LIBPAPI
  data[0] = walltime(0);
  data[1] = cpucycles();
#else
  int array_len = *data_len - 1;
  unsigned long long t0 = PAPI_get_virt_usec();
  papi_res = PAPI_start_counters((int*)papi_events, array_len);
  if (papi_res)
    m4ri_die("");
#endif
  mzd_t *C = mzd_mul(NULL, A, B, p->cutoff);
#ifndef HAVE_LIBPAPI
  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);
#else
  PAPI_stop_counters((long long*)&data[1], array_len);
  t0 = PAPI_get_virt_usec() - t0;
  data[0] = t0;
  for (int nv = 0; nv <= array_len; ++nv) {
    data[nv] -= loop_calibration[nv];
  }
#endif
  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  return (0);
}

void print_help_and_exit() {
  printf("Parameter(s) n(, cutoff) expected.\n");
  printf(" n      -- integer > 0\n");
  printf(" cutoff -- integer >= 0 (optional, default: 0).\n");
  printf("\n");
  bench_print_global_options(stderr);
  m4ri_die("");
}

int papi_test(int * papi_events, int papi_array_len) {
  int res = PAPI_start_counters(papi_events, papi_array_len);
  switch(res) {
  case 0: {
    long long* tmp = (long long*)malloc(papi_array_len * sizeof(long long));
    PAPI_stop_counters(tmp, papi_array_len);
    free(tmp);
    break;
  }
  case PAPI_ECNFLCT: {
    fprintf(stderr, "%s: %s: Conflicting event: The underlying counter hardware cannot count the specified events simultaneously.\n", progname, papi_event_name(papi_events[papi_array_len - 1]));
    fprintf(stderr, "Run `papi_event_chooser PRESET");
    for (int nv = 0; nv < papi_array_len - 1; ++nv)
      fprintf(stderr, " %s", papi_event_name(papi_events[nv]));
    fprintf(stderr, "` to get a list of possible events that can be added.\n");
    break;
  }
  case PAPI_ENOEVNT: {
    for (int nv = 0; nv < papi_array_len; ++nv)
      if ((res = PAPI_query_event(papi_events[nv])) != PAPI_OK) {
        fprintf(stderr, "%s: PAPI_start_counters: %s: %s.\n", progname, papi_event_name(papi_events[nv]), PAPI_strerror(res));
        break;
      }
    break;
  }
  case PAPI_ESYS:
    fprintf(stderr, "%s: PAPI_start_counters: %s\n", progname, strerror(errno));
    break;
  default:
    fprintf(stderr, "%s: PAPI_start_counters: %s.\n", progname, PAPI_strerror(res));
    break;
  }
  return res;
}

int main(int argc, char **argv) {
  int opts = global_options(&argc, &argv);
  int data_len;
  struct mul_params params;

#ifdef HAVE_LIBPAPI
  int papi_counters = PAPI_num_counters();
  if (papi_counters < papi_array_len) {
    fprintf(stderr, "%s: Warning: there are only %d hardware counters available!\n", progname, papi_counters);
    papi_array_len = papi_counters;
  }
  if (papi_test(papi_events, papi_array_len))
    exit(1);

  for (int nv = 0; nv <= papi_array_len; ++nv)
    loop_calibration[nv] = 100000000;

  data_len = papi_array_len + 1;
#else
  data_len = 2;
#endif

  if (opts < 0 || argc < 2 || argc > 3) {
    print_help_and_exit();
  }

  params.n = atoi(argv[1]);
  if (argc == 3)
    params.cutoff = atoi(argv[2]);
  else
    params.cutoff = 0;

  if (params.n <= 0) {
    m4ri_die("Parameter n must be > 0\n");
  }

  srandom(17);

  unsigned long long data[16];

  for (int i = 0; i < 100; ++i)
    run_nothing((void*)&params, data, &data_len);

  run_bench(run, (void*)&params, data, data_len);

  double cc_per_op = ((double)data[1])/ powl((double)params.n,2.807);
  printf("n: %5d, cutoff: %5d, cpu cycles: %12llu, cc/n^2.807: %.5lf, ", params.n, params.cutoff, data[1], cc_per_op);
  print_wall_time(data[0] / 1000000.0);
  printf("\n");

#ifdef HAVE_LIBPAPI
  for (int n = 1; n < data_len; ++n) {
    double tmp = ((double)data[n]) / powl((double)params.n,2.807);
    printf("%20s (%20llu) per bit (divided by n^2.807): %15.5f\n", papi_event_name(papi_events[n - 1]), data[n], tmp);
  }
#endif
}
