// version 20251226
// public domain
// djb

// 20251226 djb: add ticks_close()
// 20251226 djb: doc reorg, same code
// 20230105 djb: adapted from supercop/cpucycles/vct.c

#include "cpucycles_internal.h"

long long ticks(void)
{
  long long result;
  asm volatile("mrs %0, CNTVCT_EL0" : "=r" (result));
  return result;
}

void ticks_close(void)
{
}

long long ticks_setup(void)
{
  if (!cpucycles_works(ticks)) return cpucycles_SKIP;
  return cpucycles_FINDMULTIPLIER;
}
