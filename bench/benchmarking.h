#ifndef BENCHMARKETING_H
#define BENCHMARKETING_H

#include <stdint.h>

/*
 * Command line options. See benchmarking.h for documentation.
 */
extern int bench_quiet;
extern int bench_dump;
extern int bench_minimum;
extern int bench_maximum;
extern unsigned long long bench_maxtime;
extern double bench_accuracy;
extern int bench_confidence_index;
extern char const *progname;
extern uint64_t bench_count;

unsigned long long walltime(unsigned long long t0);

int global_options(int *argcp, char ***argvp);
void bench_print_global_options(FILE *);

int run_bench(int (*f)(void *params, unsigned long long *data, int *data_len), void *params,
              unsigned long long *data, int data_len);

#ifdef HAVE_LIBPAPI
extern int papi_events[];
extern int papi_array_len;
char *papi_event_name(int event);
int papi_test(int *papi_events, int papi_array_len);
#endif

void print_wall_time(double seconds);
void print_cpu_time(double seconds);

#endif  // BENCHMARKETING_H
