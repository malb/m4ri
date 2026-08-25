// version 20251226
// public domain
// djb

// 20251226 djb: add ticks_close()
// 20251226 djb: wrap inside perf_event_open (otherwise current kernels disable rdcycle)
// 20251226 djb: fork from riscv64-rdcycle.c

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/syscall.h>
#include <linux/perf_event.h>
#include "cpucycles_internal.h"

long long ticks(void)
{
  long long result;
  asm volatile("rdcycle %0" : "=r" (result));
  return result;
}

#ifndef PERF_FLAG_FD_CLOEXEC
#define PERF_FLAG_FD_CLOEXEC 0
#endif

static int fdperf = -1;

void ticks_close(void)
{
  if (fdperf >= 0) close(fdperf);
  fdperf = -1;
}

long long ticks_setup(void)
{
  int xkernel,xhv;

  if (fdperf == -1) {
    // prefer to exclude if platform supports that
    // but fall back to non-exclude if necessary
    for (xkernel = 1;xkernel >= 0;--xkernel) {
      for (xhv = 1;xhv >= 0;--xhv) {
        static struct perf_event_attr attr;

        memset(&attr,0,sizeof attr);
        attr.type = PERF_TYPE_HARDWARE;
        attr.size = sizeof(struct perf_event_attr);
        attr.config = PERF_COUNT_HW_CPU_CYCLES;
        attr.exclude_kernel = xkernel;
        attr.exclude_hv = xhv;

        fdperf = syscall(__NR_perf_event_open,&attr,0,-1,-1,PERF_FLAG_FD_CLOEXEC);
        if (fdperf != -1) break;
      }
      if (fdperf != -1) break;
    }
    if (fdperf == -1) return cpucycles_SKIP;
  }

  if (!cpucycles_works(ticks)) return cpucycles_SKIP;
  return cpucycles_CYCLECOUNTER;
}
