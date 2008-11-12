/* wall clock */
/* see "Software Optimization for High Performance Computing" p. 135 */

#include <time.h>
#include <sys/time.h>

double walltime( double *t0 )
{
  double mic, time;
  double mega = 0.000001;
  struct timeval tp;
  static long base_sec = 0;
  static long base_usec = 0;

  (void) gettimeofday(&tp,NULL);
  if (base_sec == 0)
    {
      base_sec = tp.tv_sec;
      base_usec = tp.tv_usec;
    }

  time = (double) (tp.tv_sec - base_sec);
  mic = (double) (tp.tv_usec - base_usec);
  time = (time + mic * mega) - *t0;
  return(time);
}
