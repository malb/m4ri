// version 20251226
// public domain
// djb

// 20251226 djb: add ticks_close()
// 20251226 djb: use cpuid_tsc_invariant_frequency()
// 20230105 djb: adapted from supercop/cpucycles/amd64tscfreq.c

#include "cpucycles_internal.h"
#include "cpuid_intel.h"

long long ticks(void)
{
  unsigned long long result;
  asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
    : "=a"(result) :: "%rdx");
  return result;
}

void ticks_close(void)
{
}

long long ticks_setup(void)
{
  if (!cpucycles_works(ticks)) return cpucycles_SKIP;
  if (cpucycles_works(cpuid_tsc_invariant_frequency)) {
    long long freq = cpuid_tsc_invariant_frequency();
    if (freq) return freq;
  }
  return cpucycles_MAYBECYCLECOUNTER;
}
