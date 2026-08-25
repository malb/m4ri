// version 20251226
// public domain
// djb

// 20251226 djb: add ticks_close()
// 20230105 djb: first version

#include <mach/mach_time.h>
#include "cpucycles_internal.h"

long long ticks(void)
{
  return mach_absolute_time();
}

void ticks_close(void)
{
}

long long ticks_setup(void)
{
  if (!cpucycles_works(ticks)) return cpucycles_SKIP;
  return cpucycles_FINDMULTIPLIER;
}
