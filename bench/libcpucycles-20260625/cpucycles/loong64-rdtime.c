// version 20251226
// public domain
// djb

// 20251226 djb: add ticks_close()
// 20250924 djb: first version

#include "cpucycles_internal.h"

long long ticks(void)
{
  long long result;
  asm volatile("rdtime.d %0, $zero\n" : "=r" (result));
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
