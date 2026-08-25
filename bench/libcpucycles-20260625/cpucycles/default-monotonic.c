// version 20251226
// public domain
// djb

// 20251226 djb: add ticks_close()
// 20251226 djb: doc reorg, same code
// 20230105 djb: adapted from supercop/cpucycles/monotonic.c

#include <time.h>
#include <sys/time.h>

long long ticks_setup(void)
{
  return 1000000000;
}

void ticks_close(void)
{
}

long long ticks(void)
{
  struct timespec t;
  long long result;
  clock_gettime(CLOCK_MONOTONIC,&t);
  result = t.tv_sec;
  result *= 1000000000;
  result += t.tv_nsec;
  return result;
}
