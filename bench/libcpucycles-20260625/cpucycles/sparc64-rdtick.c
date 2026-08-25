// version 20251226
// public domain
// djb

// 20251226 djb: add ticks_close()
// 20230105 djb: adapted from supercop/cpucycles/sparccpuinfo.c

#include "cpucycles_internal.h"

long long ticks(void)
{
  long long result;
  asm volatile("rd %%tick,%0" : "=r" (result));
  return result;
}

void ticks_close(void)
{
}

long long ticks_setup(void)
{
  if (!cpucycles_works(ticks)) return cpucycles_SKIP;
  return cpucycles_CYCLECOUNTER;
}
