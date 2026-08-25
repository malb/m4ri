// version 20260105
// public domain
// djb

// 20260105 djb: rename from pmc to pmcff (for the use of fixed-function counters)
// 20251226 djb: add ticks_close()
// 20251226 djb: rename previous as amd64-perf-pmc.c; remove perf from this version
// 20230105 djb: adapted from supercop/cpucycles/amd64rdpmc.c

#include "cpucycles_internal.h"

long long ticks(void)
{
  long long result;
  asm volatile("rdpmc;shlq $32,%%rdx;orq %%rdx,%%rax"
    : "=a"(result) : "c"((1<<30)|1) : "%rdx");
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
