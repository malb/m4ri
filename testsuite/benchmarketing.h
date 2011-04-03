#ifndef BENCHMARKETING_H
#define BENCHMARKETING_H

#include <stdint.h>

/*
 * Command line options. See benchmarketing.h for documentation.
 */
extern int bench_quiet;
extern int bench_dump;
extern int bench_minimum;
extern int bench_maximum;
extern double bench_maxtime;
extern double bench_accuracy;
extern int bench_confidence_index;
extern char const* progname;
extern uint64_t bench_count;

double walltime(double t0);

int global_options(int* argcp, char*** argvp);
void bench_print_global_options(FILE*);

void run_bench(
    int (*f)(void* params, double* wtp, unsigned long long* ccp),
    void* params,
    double* wtp,
    unsigned long long* ccp);

#endif //BENCHMARKETING_H
