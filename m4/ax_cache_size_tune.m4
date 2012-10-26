# SYNOPSIS
#
#   AX_CACHE_SIZE_TUNE
#
# DESCRIPTION
#
#   Find L1, L2, L3 caches size by running some timing experiments.
#   The results are available in the defines __M4RI_CPU_L1_CACHE,
#   __M4RI_CPU_L2_CACHE and __M4RI_CPU_L3_CACHE.
#
#   This macro depends on AC_PROG_SED, AC_PROG_CC.
#
# LAST MODIFICATION
#
#   2011-04-11
#
# COPYLEFT
#
#   Copyright (c) 2009,2010 Martin Albrecht <martinralbrecht@googlemail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_CACHE_SIZE_TUNE],
[ AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PROG_SED])

  AC_LANG_PUSH([C])
  AC_CACHE_CHECK(for cache sizes, ax_cv_cache_sizes,
  [AC_RUN_IFELSE([AC_LANG_PROGRAM([[
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>

double walltime(double t0) {
  double mic, time;
  double mega = 0.000001;
  struct timeval tp;
  static long base_sec = 0;
  static long base_usec = 0;

  (void) gettimeofday(&tp,NULL);
  if (base_sec == 0) {
    base_sec = tp.tv_sec;
    base_usec = tp.tv_usec;
  }

  time = (double) (tp.tv_sec - base_sec);
  mic = (double) (tp.tv_usec - base_usec);
  time = (time + mic * mega) - t0;
  return(time);
}

double run_experiment(size_t size, size_t trials) {
  size_t i,j;
  unsigned long *a = (unsigned long*)malloc(size/4);
  unsigned long *b = (unsigned long*)malloc(size/4);
  unsigned long *c = (unsigned long*)malloc(size/4);
  unsigned long *d = (unsigned long*)malloc(size/4);

  size_t n = size/4/(sizeof(unsigned long));

  /* we setup a lookup table with a random-ish pattern */
  a[0] = 1337;
  b[0] = 5345345;
  for(j=1; j<n ; j++) {
    a[j] = a[j-1] * 1073741827;
    a[j-1] %= n;
    b[j] = b[j-1] * 1073741827;
    b[j-1] %= n;
  }
  a[n-1] %= n;
  b[n-1] %= n;

  double wt = walltime(0.0);
  clock_t ct = clock();
  for(i=0; i<trials; i++) {
    for(j=0; j<n; j+=8) {
      d[b[j+0]] = c[a[j+0]];
      d[b[j+1]] = c[a[j+1]];
      d[b[j+2]] = c[a[j+2]];
      d[b[j+3]] = c[a[j+3]];
      d[b[j+4]] = c[a[j+4]];
      d[b[j+5]] = c[a[j+5]];
      d[b[j+6]] = c[a[j+6]];
      d[b[j+7]] = c[a[j+7]];
    }
  }
  ct = clock() - ct;
  wt = walltime(wt);
  free(a);
  free(b);
  free(c);
  free(d);
  return (double)wt;
}


#define NUMBER_OF_EXPERIMENTS 8

size_t cache_size(const size_t *candidates, const size_t n, size_t trials) {
  double times[NUMBER_OF_EXPERIMENTS][n];
  double dtimes[NUMBER_OF_EXPERIMENTS][n];
  size_t i,j;
  double wt, result;


  for(j=0; j<NUMBER_OF_EXPERIMENTS; j++) {
    size_t mult  = 1;
    size_t _trials = trials;
    run_experiment(candidates[0]*1024,_trials);
    wt = walltime(0.0);
    times[j][0] = run_experiment(candidates[0]*1024,_trials);
    wt = walltime(wt);
    dtimes[j][0] = 1.0;
    printf("s: %5zu, rx: %6.2f, x: %6.2f, wt: %6.2f, dx:    NaN\n",candidates[0],times[j][0],times[j][0],wt);
    fflush(NULL);

    for(i=1;i<n;i++) {
      run_experiment(candidates[i]*1024,_trials);
      wt = walltime(0.0);
      result = run_experiment(candidates[i]*1024,_trials);
      wt = walltime(wt);
      times[j][i] = mult*result;
      dtimes[j][i] = candidates[i-1]*times[j][i]/times[j][i-1]/candidates[i];

      printf("s: %5zu, rx: %6.2f, x: %6.2f, wt: %6.2f, dx: %6.2f\n",candidates[i],result,times[j][i],wt,dtimes[j][i]);
      fflush(NULL);

      while(wt > 0.25) {
        _trials = _trials/2;
        mult = 2*mult;
        wt /= 2.0;
        result /= 2.0;
      }
    }
    printf("\n");
  }
  for(i=0;i<n;i++) {
    double tmp = 0.0;
    for(j=0; j<NUMBER_OF_EXPERIMENTS; j++) {
      tmp += dtimes[j][i];
    }
    dtimes[0][i] = tmp/NUMBER_OF_EXPERIMENTS;
  }

  size_t max = 1;
  for(i=1;i<n;i++){
    if (dtimes[0][i] > dtimes[0][max] ) {
      max = i;
    }
  }
  return candidates[max-1];
}
      ]],
      [[
  const size_t c1[] = {   4,   8,  16,  32,  64, 128};
  const size_t c2[] = { 128, 256, 512};
  const size_t c3[] = {1024,1536,2048,3072,4096,6144,8192,16384,32768};

  FILE *f;
  printf("\n");
  size_t _l1 = cache_size(c1,  6, 1ULL<<15);
  size_t _l2 = cache_size(c2,  3, 1ULL<<12);
  size_t _l3 = cache_size(c3,  9, 1ULL<< 9);

  f = fopen("conftest_cache_sizes", "w"); if (!f) return 1;
  fprintf(f,"%lu:%lu:%lu\n",(unsigned long)(_l1*1024),(unsigned long)(_l2*1024),(unsigned long)(_l3*1024));
  fclose(f);
  return 0;
   ]])],
    [ax_cv_cache_sizes=`cat conftest_cache_sizes`; rm -f conftest_cache_sizes],
    [ax_cv_cache_sizes=unknown; rm -f conftest_cache_sizes],
    [ax_cv_cache_sizes=unknown])])
  AC_LANG_POP([C])
  AC_MSG_CHECKING(the L1 cache size)
  ax_l1_size=`echo $ax_cv_cache_sizes | cut -d ':' -f 1`
  AC_MSG_RESULT( $ax_l1_size Bytes)

  AC_MSG_CHECKING(the L2 cache size)
  ax_l2_size=`echo $ax_cv_cache_sizes | cut -d ':' -f 2`
  AC_MSG_RESULT( $ax_l2_size Bytes)

  AC_MSG_CHECKING(the L3 cache size)
  ax_l3_size=`echo $ax_cv_cache_sizes | cut -d ':' -f 3`
  AC_MSG_RESULT( $ax_l3_size Bytes)

  M4RI_CPU_L1_CACHE=${ax_l1_size}
  M4RI_CPU_L2_CACHE=${ax_l2_size}
  M4RI_CPU_L3_CACHE=${ax_l3_size}
  AC_SUBST(M4RI_CPU_L1_CACHE)
  AC_SUBST(M4RI_CPU_L2_CACHE)
  AC_SUBST(M4RI_CPU_L3_CACHE)
])
