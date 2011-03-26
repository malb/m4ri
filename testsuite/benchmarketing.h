#ifndef BENCHMARKETING_H
#define BENCHMARKETING_H

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

double walltime(double t0);

void global_options(int* argcp, char*** argvp);

void run_bench(
    int (*f)(void* params, double* wtp, unsigned long long* ccp),
    void* params,
    double* wtp,
    unsigned long long* ccp);

#endif //BENCHMARKETING_H
