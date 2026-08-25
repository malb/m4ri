// version 20251226
// public domain
// djb

// 20251226 djb: add ticks_close()
// 20230105 djb: adapted from supercop/cpucycles/riscv.c, which has code from djb and Romain Dolbeau

#include "cpucycles_internal.h"

long long ticks(void)
{
  unsigned int low, high, newhigh;
  unsigned long long result;

  asm volatile( "start%=:\n"
                "rdcycleh %0\n"
                "rdcycle %1\n"
                "rdcycleh %2\n"
                "bne %0, %2, start%=\n"
                : "=r"(high), "=r"(low), "=r"(newhigh));
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
  return cpucycles_CYCLECOUNTER;
}
