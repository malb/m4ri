#ifndef BENCHMARKETING_H
#define BENCHMARKETING_H

#include <stdio.h>
#include "walltime.h"

#define MAX_RUNTIME 60.0
#define CONFIDENCE 0.99

int run_bench(int (*f)(void *params, double *wt, unsigned long long *cc), 
              void *params,
              double *wt, 
              unsigned long long *cc) {
  
  double wt_diff = 0.0;
  double start_walltime = walltime(0.0);
  size_t i;

  for(i=0; i<1000; i++) {
    printf(".");
    fflush(stdout);
    double _wt = 0.0;
    unsigned long long _cc = 0;
    int r = f(params,&_wt,&_cc);
    if ( r ) {
      m4ri_die("benchmark function failed with exit code: %d\n",r);
    }

    if (i>=1) {
      wt_diff = ((*wt + _wt)/(i+2)) / ((*wt)/(i+1));
      if (wt_diff > 1.0) {
        wt_diff = wt_diff - 1.0;
      } else if (wt_diff < 1.0) {
        wt_diff = 1.0 - wt_diff;
      }

      if( wt_diff <  1.0 - CONFIDENCE  ) {
        *wt = *wt + _wt;
        *wt = *wt/(i+1);
        *cc = *cc + _cc;    
        *cc = *cc/(i+1);
        break;
      }
    }

    *wt = *wt + _wt;
    *cc = *cc + _cc;    

    if (walltime(start_walltime) > MAX_RUNTIME) {
      *wt = *wt / (i+1);
      *cc = *cc / (i+1);
      break;
    }
  };
  printf(" last delta: %3.2lf\%, total running time: %5.2lf\n",wt_diff*100,walltime(start_walltime));
}

#endif //BENCHMARKETING_H
