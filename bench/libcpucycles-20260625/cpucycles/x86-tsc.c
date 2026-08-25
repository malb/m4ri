// version 20251226
// public domain
// djb

// 20251226 djb: add ticks_close()
// 20251226 djb: use cpuid_tsc_invariant_frequency()

#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

#include "cpucycles_internal.h"
#include "cpuid_intel.h"

long long ticks(void)
{
  return __rdtsc();
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
