#ifndef M4RI_M4RI_CONFIG_H
#define M4RI_M4RI_CONFIG_H

// Defines determined during configuration of m4ri.
#define __M4RI_HAVE_MM_MALLOC		0
#define __M4RI_HAVE_POSIX_MEMALIGN	0
#define __M4RI_HAVE_SSE2		0
#define __M4RI_HAVE_OPENMP		0
#define __M4RI_CPU_L1_CACHE		32768
#define __M4RI_CPU_L2_CACHE		262144
#define __M4RI_CPU_L3_CACHE		2147483648
#define __M4RI_DEBUG_DUMP		(0 || 0)
#define __M4RI_DEBUG_MZD		0

// Helper macros.
#define __M4RI_USE_MM_MALLOC		(__M4RI_HAVE_MM_MALLOC && __M4RI_HAVE_SSE2)
#define __M4RI_USE_POSIX_MEMALIGN	(__M4RI_HAVE_POSIX_MEMALIGN && __M4RI_HAVE_SSE2)
#define __M4RI_DD_QUIET			(0 && !0)

#endif // M4RI_M4RI_CONFIG_H
