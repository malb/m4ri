#include <sys/types.h>
#include <sys/sysctl.h>
#include <mach/mach_time.h>

#define timebase mach_absolute_time

static int cpumib[2] = { CTL_HW, HW_CPU_FREQ } ;
static int tbmib[2] = { CTL_HW, HW_TB_FREQ } ;

static long myround(double u)
{
  long result = u;
  while (result + 0.5 < u) result += 1;
  while (result - 0.5 > u) result -= 1;
  return result;
}

static long tbcycles = 0;

static void init(void)
{
  int cpufrequency = 0; size_t cpufrequencylen = sizeof(int);
  int tbfrequency = 0; size_t tbfrequencylen = sizeof(int);
  sysctl(cpumib,2,&cpufrequency,&cpufrequencylen,0,0);
  sysctl(tbmib,2,&tbfrequency,&tbfrequencylen,0,0);
  if (tbfrequency > 0)
    tbcycles = myround((double) cpufrequency / (double) tbfrequency);
}

long long cpucycles_powerpcmacos(void)
{
  if (!tbcycles) init();
  return timebase() * tbcycles;
}

long long cpucycles_powerpcmacos_persecond(void)
{
  int result = 0; size_t resultlen = sizeof(int);
  sysctl(cpumib,2,&result,&resultlen,0,0);
  return result;
}
