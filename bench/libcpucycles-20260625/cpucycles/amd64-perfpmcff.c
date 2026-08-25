// version 20260625
// public domain
// djb

// 20260625 djb: add second fd for atom cores, sigh
// 20260105 djb: forked from amd64-perfpmc.c
// 20251226 djb: add ticks_close()
// 20251226 djb: include xkernel loop and xhv loop
// 20251226 djb: skip mmap
// 20251226 djb: always use rdpmc with PERF_FIXED_CTR1 (works around perf on e.g. 6.1.0 returning index=0 for e-cores)
// 20251225 djb: use PERF_FLAG_FD_CLOEXEC
// 20230105 djb: adapted from supercop/cpucycles/amd64rdpmc.c

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
  asm volatile("rdpmc;shlq $32,%%rdx;orq %%rdx,%%rax"
    : "=a"(result) : "c"((1<<30)|1) : "%rdx");
  return result;
}

#ifndef PERF_FLAG_FD_CLOEXEC
#define PERF_FLAG_FD_CLOEXEC 0
#endif

static int fdperf = -1;
static int fdperfatom = -1;

void ticks_close(void)
{
  if (fdperfatom >= 0) close(fdperfatom);
  fdperfatom = -1;
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
        if (fdperf != -1) {
#ifdef PERF_PMU_TYPE_SHIFT
          // XXX: maybe read /sys/bus/event_source/devices/cpu_atom/type
          // XXX: maybe allow more fds for more cpu_* directories
          // XXX: maybe ask linux to add simpler API
          attr.config |= ((long long) 10) << PERF_PMU_TYPE_SHIFT;
          fdperfatom = syscall(__NR_perf_event_open,&attr,0,-1,-1,PERF_FLAG_FD_CLOEXEC);
#endif
          break;
        }
      }
      if (fdperf != -1) break;
    }
    if (fdperf == -1) return cpucycles_SKIP;
  }

  if (!cpucycles_works(ticks)) return cpucycles_SKIP;
  return cpucycles_CYCLECOUNTER;
}
