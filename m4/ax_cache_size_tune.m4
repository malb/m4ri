# SYNOPSIS
#
#   AX_CACHE_SIZE_TUNE
#
# DESCRIPTION
#
#   Find L1 and L2 caches size by running some timing experiments.
#   The results are available in the defines CPU_L1_CACHE and
#   CPU_L2_CACHE.
#
#   This macro depends on AC_PROG_SED, AC_PROG_CC.
#
# LAST MODIFICATION
#
#   2009-11-01
#
# COPYLEFT
#
#   Copyright (c) 2009 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
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
      wt = walltime(wt);
      free(a);
      free(b);
      free(c);
      free(d);
      return wt;
    }
        
    size_t cache_size(const size_t *candidates, size_t trials, size_t n){
      size_t mult  = 1;
      double times[n];
      double dtimes[n];
      size_t i;
    
      run_experiment(candidates[0]*1024,trials);
      times[0] = run_experiment(candidates[0]*1024,trials);
      dtimes[0] = 1.0;
      printf("%5zu %9.3f       XXX\n",candidates[0],times[0]);
    
      for(i=1;i<n;i++) {
        run_experiment(candidates[i]*1024,trials);
        double result = run_experiment(candidates[i]*1024,trials);
        times[i] = mult*result;
        dtimes[i] = candidates[i-1]*times[i]/times[i-1]/candidates[i];
    
        printf("%5zu %9.3f %9.3f\n",candidates[i],times[i],dtimes[i]);
    
        if (result > 2) {
          trials = trials/2;
          mult = 2*mult;
        }
        if (result > 3) {
          trials = trials/3;
          mult = 3*mult;
        }
      }
      double avg = 0.0;
      for(i=0;i<n;i++)
        avg += dtimes[i];
      avg /= n;

      size_t min = n-1;
      for(i=0;i<n;i++)
       if(dtimes[i] < dtimes[min])
          min = i;
    
      size_t max = 1;
      for(i=1;i<n;i++){
        if (dtimes[i] > dtimes[max] && dtimes[i-1] <= 1.2*avg) {
          max = i;
        }
        if (dtimes[i] < 0.8*dtimes[max])
          break;
      }
      return candidates[max-1];
    }
      ]], 
      [[
      FILE *f;
      printf("\n");
      const size_t c1[] = {4,8,16,32,64,128,256,512};
      const size_t _l1 = cache_size(c1,50000,8);
      printf("\n");
      const size_t c2[] = {512,1024,1536,2048,3072,4096,6144,8192,16384,32768};
      const size_t _l2 = cache_size(c2,5000,10);
      f = fopen("conftest_cache_sizes", "w"); if (!f) return 1;
      fprintf(f,"%lu:%lu\n",(unsigned long)(_l1*1024),(unsigned long)(_l2*1024));
      fclose(f);
      return 0;
   ]])],
    [ax_cv_cache_sizes=`cat conftest_cache_sizes`; rm -f conftest_cache_sizes],
    [ax_cv_cache_sizes=unknown; rm -f conftest_cache_sizes],
    [ax_cv_cache_sizes=unknown])])
  AC_LANG_POP([C])
  AC_MSG_CHECKING(the L1 cache size)
  ax_l1_size=`echo $ax_cv_cache_sizes | $SED 's/\:.*//g'`
  AC_MSG_RESULT( $ax_l1_size Bytes)

  AC_MSG_CHECKING(the L2 cache size)
  ax_l2_size=`echo $ax_cv_cache_sizes | $SED 's/.*\://g'`
  AC_MSG_RESULT( $ax_l2_size Bytes)

  AC_DEFINE_UNQUOTED([CPU_L1_CACHE], ${ax_l1_size}, [L1 cache size (in Bytes)])
  AC_DEFINE_UNQUOTED([CPU_L2_CACHE], ${ax_l2_size}, [L2 cache size (in Bytes)])
])
