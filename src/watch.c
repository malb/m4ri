/******************************************************************************
*
*            M4RI: Method of the Four Russians Inversion
*
*       Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org> 
*
*  Distributed under the terms of the GNU General Public License (GPL)
*
*    This code is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    General Public License for more details.
*
*  The full text of the GPL is available at:
*
*                  http://www.gnu.org/licenses/
******************************************************************************/

#include "watch.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>

clock_t watchstart;
clock_t watchstop;
long watchrunning;

double logtotal, logtotalsquare, logcount;

void m4ri_watch_seed() {
  srand(clock()/10000);
}

double m4ri_sqrt(double n) {
  double guess=n/2;
  double miss, der;
  int i;

  if (n<0) {
    printf("\a Tried to take square root of %f.\n",n);
    exit(1);
  }

  for (i=0; i<40 ;i++) {
    miss=n-guess*guess;
    der=2.0*guess;
    guess+=miss/der;
  }

  return guess;
}

void m4ri_watch_start() {
  watchstart=clock();
  watchrunning=TRUE;
}

void m4ri_watch_stop() {
  watchstop=clock();
  watchrunning=FALSE;
}

clock_t m4ri_watch_get() {
  if (watchrunning==TRUE) 
    return clock()-watchstart;
  else 
    return watchstop-watchstart;
}

void m4ri_watch_clear_logs() {
  logtotal=0;
  logtotalsquare=0;
  logcount=0;
}

void m4ri_watch_store(clock_t watch) {
  double value=(double)watch;

  logtotal+=value;
  logtotalsquare+=value*value;
  logcount++;
}

double m4ri_watch_get_average() {
  return logtotal/(logcount*1.0);
}

double m4ri_watch_get_sigma() {
  double mean = m4ri_watch_get_average();
  double ex_sq = logtotalsquare/(logcount*1.0);
  double var = ex_sq-mean*mean;
  double sigma = m4ri_sqrt(var);
  return sigma;
}

long m4ri_watch_get_count() {
  return (long)logcount;
}




