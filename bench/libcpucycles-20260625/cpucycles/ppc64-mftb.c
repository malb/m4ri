// version 20251226
// public domain
// djb

// 20251226 djb: add ticks_close()
// 20230105 djb: adapted from supercop/cpucycles/powerpccpuinfo.c

#include "cpucycles_internal.h"

long long ticks(void)
{
  unsigned int high, low, newhigh;
  unsigned long long result;

  do {
    asm volatile(
      "mftbu %0; mftb %1; mftbu %2"
      : "=r" (high), "=r" (low), "=r" (newhigh)
    );
  } while (newhigh != high);

  result = high;
  result <<= 32;
  result |= low;
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
