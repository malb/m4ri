#ifndef amd64_cpuid_h
#define amd64_cpuid_h

#ifdef __FILC__

#include <cpuid.h>
static void cpuid_0(unsigned int leaf,unsigned int *a,unsigned int *b,unsigned int *c,unsigned int *d)
{
  __get_cpuid(leaf,a,b,c,d);
}

#elif defined(_MSC_VER)

static void cpuid_0(unsigned int leaf,unsigned int *a,unsigned int *b,unsigned int *c,unsigned int *d)
{
  int abcd[4];
  __cpuidex(abcd,leaf,0);
  *a = abcd[0];
  *b = abcd[1];
  *c = abcd[2];
  *d = abcd[3];
}

#else

static void cpuid_0(unsigned int op,unsigned int *a,unsigned int *b,unsigned int *c,unsigned int *d)
{
  asm volatile("cpuid" : "=a"(*a),"=b"(*b),"=c"(*c),"=d"(*d) : "a"(op) : "memory");
}

#endif

__attribute__((unused))
static long long cpuid_says_genuineintel(void)
{
  unsigned int a,b,c,d;
  cpuid_0(0,&a,&b,&c,&d);
  return b == 0x756e6547 && d == 0x49656e69 && c == 0x6c65746e;
}

// 0 if not clear that tsc is invariant
__attribute__((unused))
static long long cpuid_tsc_invariant_frequency(void)
{
  unsigned int a,b,c,d,top,bot;
  cpuid_0(0x80000000,&a,&b,&c,&d);
  if (a < 0x80000007) return 0;
  cpuid_0(0x80000007,&a,&b,&c,&d);
  if (!(d & 0x100)) return 0;
  cpuid_0(0,&a,&b,&c,&d);
  if (a < 0x15) return 0;
  cpuid_0(0x15,&bot,&top,&c,&d);
  if (!c) {
    // see "Nominal Core Crystal Clock Frequency" exception table in intel manual
    cpuid_0(0x1,&a,&b,&c,&d);
    switch(a & 0x0fff0ff0) {
      case 0x000406e0: // skylake
      case 0x000506e0: // skylake
      case 0x000806e0: // kaby lake (or comet lake)
      case 0x000906e0: // kaby lake (or coffee lake)
      case 0x000a0650: // comet lake
      case 0x000a0660: // comet lake
        c = 24000000; break; // "6th and 7th generation Intel Core processors and Intel Xeon W Processor Family"
      case 0x00050650: // cascade lake
        c = 25000000; break; // "Intel Xeon Scalable Processor Family with CPUID signature 06_55H"
      case 0x000506c0: // goldmont
        c = 19200000; break; // "Next Generation Intel Atom processors based on Goldmont Microarchitecture with CPUID signature 06_5CH"
    }
  }
  return (top * (long long) c) / bot;
}

#endif
