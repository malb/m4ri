// version 20260105
// public domain
// djb

// 20260105 djb: disable GenuineIntel here given perf bug on e-cores (rather than trusting wrapper to detect perfpmcff being faster)
// 20260105 djb: go back to buf->index for AMD support
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
#include <sys/mman.h>
#include <sys/syscall.h>
#include <linux/perf_event.h>
#include "cpucycles_internal.h"
#include "cpuid_intel.h"

static struct perf_event_mmap_page *buf = 0;

long long ticks(void)
{
  long long result;
  unsigned int seq;
  long long index;
  long long offset;

  do {
    seq = buf->lock;
    asm volatile("" ::: "memory");
    index = buf->index;
    offset = buf->offset;
    asm volatile("rdpmc;shlq $32,%%rdx;orq %%rdx,%%rax"
      : "=a"(result) : "c"(index-1) : "%rdx");
    asm volatile("" ::: "memory");
  } while (buf->lock != seq);

  result += offset;
  result &= 0xffffffffffff;
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

static long long abort_on_intel(void)
{
  if (cpuid_says_genuineintel()) abort();
  return 0;
}

long long ticks_setup(void)
{
  int xkernel,xhv;

  if (!cpucycles_works(abort_on_intel)) return cpucycles_SKIP;

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
          buf = mmap(NULL,sysconf(_SC_PAGESIZE),PROT_READ,MAP_SHARED,fdperf,0);
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
