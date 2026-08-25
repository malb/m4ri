// version 20251226
// public domain
// djb

// 20251226 djb: add ticks_close()
// 20230105 djb: first version

#include "cpucycles_internal.h"

long long ticks_setup(void)
{
  return cpucycles_SKIP;
}

void ticks_close(void)
{
}

long long ticks(void)
{
  return 0;
}
