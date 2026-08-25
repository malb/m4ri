// version 20251226
// public domain
// djb

// 20251226 djb: add ticks_close()
// 20230105 djb: adapted from supercop/cpucycles/riscv.c, which has code from djb and Romain Dolbeau

#include "cpucycles_internal.h"

long long ticks(void)
{
  long long result;
  asm volatile("rdcycle %0" : "=r" (result));
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
