AC_INIT(m4ri,20200125)

AC_CANONICAL_HOST

AC_CONFIG_SRCDIR(m4ri/brilliantrussian.c)

AM_INIT_AUTOMAKE

dnl Include maintainer mode targets.
AM_MAINTAINER_MODE

dnl Needed when reconfiguring with 'autoreconf -i -s'
AC_CONFIG_MACRO_DIR([m4])

dnl Compiling with per-target flags (test_elimination.c) requires AM_PROG_CC_C_O.
AM_PROG_CC_C_O

LT_INIT

AC_PROG_INSTALL

AC_CONFIG_HEADERS(m4ri/config.h)

AC_PROG_CC_C99()
if test "$ac_cv_prog_cc_c99" = "no"; then
  AC_MSG_ERROR([C99 support is required but not found.])
fi

# Find and set LIBM on the platform.
LT_LIB_M
AC_SUBST(LIBM)

# SSE2 support
AC_ARG_ENABLE([sse2],
        AS_HELP_STRING([--disable-sse2], [don't use SSE2 instruction set.]),
	, [if test "$m4ri_wrapword" = "yes"; then enable_sse2="no"; else enable_sse2="yes"; fi])

AS_IF([test "x$enable_sse2" != "xno"], [
   if test "$m4ri_wrapword" = "yes"; then
      AC_MSG_ERROR([SSE2 cannot be supported when wrapping word in a C++ class.])
   fi
   case $host_cpu in i[[3456]]86*|x86_64*)
          AX_EXT()
   esac
])
if test x"$ax_cv_have_sse2_ext" = x"yes"; then
  M4RI_HAVE_SSE2=1
else
  M4RI_HAVE_SSE2=0
fi
AC_SUBST(M4RI_HAVE_SSE2)


AC_ARG_WITH(papi,
    AS_HELP_STRING([--with-papi@<:@=PATH@:>@], [The PAPI install prefix, if configure can't find it.]),
    [m4ri_config_papi=$withval])

AC_ARG_WITH(cachesize,
        AS_HELP_STRING([--with-cachesize@<:@=VALUE@:>@], [L1,L2 and L3 cache sizes in bytes, separated by a colon. Overrides cache tuning.]),[m4ri_config_cachesize=$withval])

AC_CHECK_HEADER([mm_malloc.h],AC_DEFINE(HAVE_MM_MALLOC,,[Support aligned allocations]),)
if test "$ac_cv_header_mm_malloc_h" = "yes"; then
  M4RI_HAVE_MM_MALLOC=1
else
  M4RI_HAVE_MM_MALLOC=0
fi
AC_SUBST(M4RI_HAVE_MM_MALLOC)

# Correctly working posix_memalign
AX_FUNC_POSIX_MEMALIGN
if test "$ax_cv_func_posix_memalign_works" = "yes"; then
  M4RI_HAVE_POSIX_MEMALIGN=1
else
  M4RI_HAVE_POSIX_MEMALIGN=0
fi
AC_SUBST(M4RI_HAVE_POSIX_MEMALIGN)

# OpenMP support
AC_ARG_ENABLE([openmp],
        AS_HELP_STRING( [--enable-openmp],[add support for OpenMP multicore support.]))

AS_IF([test "x$enable_openmp" = "xyes"], [
   AX_OPENMP()
])
AC_SUBST(OPENMP_CFLAGS)
if test -n "$OPENMP_CFLAGS"; then
  M4RI_HAVE_OPENMP=1
else
  M4RI_HAVE_OPENMP=0
fi
AC_SUBST(M4RI_HAVE_OPENMP)


# Thread-Safety
M4RI_ENABLE_MZD_CACHE=1
M4RI_ENABLE_MMC=1

AC_ARG_ENABLE([thread-safe],
	AS_HELP_STRING([--enable-thread-safe], [make library thread safe.]))

if test "x$enable_thread_safe" = "xyes"; then
  M4RI_ENABLE_MMC=0
  M4RI_ENABLE_MZD_CACHE=0
fi

if test $M4RI_HAVE_OPENMP = 1; then
  M4RI_ENABLE_MZD_CACHE=0
fi

AC_SUBST(M4RI_ENABLE_MZD_CACHE)
AC_SUBST(M4RI_ENABLE_MMC)

# Debugging support
AC_ARG_ENABLE([debug],
	AS_HELP_STRING([--enable-debug], [enable assert() statements for debugging.]))

AC_ARG_ENABLE([debug-dump],
        AS_HELP_STRING([--enable-debug-dump], [dump output at exit of every function.]))

if test "x$enable_debug_dump" = "xyes"; then
  M4RI_DEBUG_DUMP=1
else
  M4RI_DEBUG_DUMP=0
fi
AC_SUBST(M4RI_DEBUG_DUMP)

AC_ARG_ENABLE([debug-mzd],
        AS_HELP_STRING([--enable-debug-mzd], [add consistency checks on matrix structures.]))

if test "x$enable_debug_mzd" = "xyes"; then
  M4RI_DEBUG_MZD=1
else
  M4RI_DEBUG_MZD=0
fi
AC_SUBST(M4RI_DEBUG_MZD)

if test "x$enable_debug" = x"yes"; then
  DEBUG_FLAGS="-g"
  AC_SUBST(DEBUG_FLAGS)
else
  if test "x$enable_debug_mzd" != "xyes"; then
    AC_DEFINE(NDEBUG,1,[define whether debugging is enabled])
  fi
fi

# For tests. Detect if PAPI is installed. See http://icl.cs.utk.edu/papi/ .

if test -z "$m4ri_config_papi"; then
   AC_CHECK_LIB(papi, PAPI_start_counters,
   [ 
    AX_GUESS_PATH_LIB(papi)
    AX_GUESS_PATH_HEADER(papi.h)
    if test -n "$LIBPAPI_PATH"; then
       PAPI_LDFLAGS="-Wl,-rpath,$LIBPAPI_PATH"
       PAPI_LIBS="-L$LIBPAPI_PATH -lpapi"
    else
        PAPI_LIBS="-lpapi"
        if ! test -e "/usr/lib/libpapi.so"; then
           AC_MSG_WARN([Could not find libpapi.so. Use --with-papi=<install_prefix> or set LD_LIBRARY_PATH correctly before running benchmark applications.])
        fi
    fi
    if test -n "$PAPI_H_PATH"; then
       PAPI_CFLAGS="-I$PAPI_H_PATH"
       AC_DEFINE_UNQUOTED([HAVE_LIBPAPI], 1, [Define when libpapi is available.])
    else
      AC_MSG_WARN([Could not find papi.h; Use --with-papi=<install_prefix> or add -I<install_prefix>/include to either CPPFLAGS 
                    or CFLAGS, or turn off papi all together by configuring with --without-papi.])
    fi
])
fi

if test x"$m4ri_config_papi" != x"no" && test -n "$m4ri_config_papi"; then
    LIBPAPI_PATH="`realpath -s $m4ri_config_papi/lib`"
    PAPI_H_PATH="`realpath -s $m4ri_config_papi/include`"
    PAPI_CFLAGS="-I$PAPI_H_PATH"
    PAPI_LDFLAGS="-Wl,-rpath,$LIBPAPI_PATH"
    PAPI_LIBS="-L$LIBPAPI_PATH -lpapi"
    AC_DEFINE_UNQUOTED([HAVE_LIBPAPI], 1, [Define when libpapi is available.])
fi

AC_SUBST(PAPI_LIBS)
AC_SUBST(PAPI_LDFLAGS)
AC_SUBST(PAPI_CFLAGS)

AC_ARG_ENABLE([cachetune], 
     AS_HELP_STRING([--enable-cachetune],[calculate cache size from timing information (deprecated).]))

# Cache Sizes
if test -z $m4ri_config_cachesize; then
   AX_CACHE_SIZE()
   AS_IF([test "x$enable_cachetune" = "xyes"], [AC_MSG_WARN(--enable-cachetune is deprecated since it usually does not provide optimal parameters.) AX_CACHE_SIZE_TUNE()])
else
  AS_IF([test "x$enable_cachetune" = "xyes"], [AC_MSG_WARN(Ignoring cache tuning since --with-cachesize was given.)])

  ax_l1_size=`echo $m4ri_config_cachesize | cut -d ":" -f 1`
  ax_l2_size=`echo $m4ri_config_cachesize | cut -d ":" -f 2`
  ax_l3_size=`echo $m4ri_config_cachesize | cut -d ":" -f 3`

  M4RI_CPU_L1_CACHE=${ax_l1_size}
  M4RI_CPU_L2_CACHE=${ax_l2_size}
  M4RI_CPU_L3_CACHE=${ax_l3_size}
  AC_SUBST(M4RI_CPU_L1_CACHE)
  AC_SUBST(M4RI_CPU_L2_CACHE)
  AC_SUBST(M4RI_CPU_L3_CACHE)

fi

# PNG

have_libpng="no"
AC_ARG_ENABLE([png],
   [AC_HELP_STRING([--disable-png], [disable PNG support @<:@default=enabled@:>@])],
   [
    if test "x${enableval}" = "xyes" ; then
       want_png="yes"
    else
       want_png="no"
    fi
   ],
   [want_png="yes"])

AC_MSG_CHECKING([whether to build with PNG support])
AC_MSG_RESULT([${want_png}])

if test "x${want_png}" = "xyes" ; then
   PKG_CHECK_MODULES([PNG], [libpng],
      [have_libpng="yes"; 
       LIBPNG_LIBADD=${PNG_LIBS};
       M4RI_USE_PNG_PC=libpng],
      [have_libpng="no"])
   if ! test "x${have_libpng}" = "xyes" ; then
      AC_CHECK_LIB([png],
         [png_create_write_struct],
         [have_libpng="yes"; LIBPNG_LIBADD="-lpng"],
         [AC_CHECK_LIB([png14],
            [png_create_write_struct],
            [have_libpng="yes"; LIBPNG_LIBADD="-lpng14"],
            [AC_CHECK_LIB([png12],
               [png_create_write_struct],
               [have_libpng="yes"; LIBPNG_LIBADD="-lpng12"],
               [AC_CHECK_LIB([png10],
                  [png_create_write_struct],
                  [have_libpng="yes"; LIBPNG_LIBADD="-lpng10"],
                  [have_libpng="no"])
               ])
            ])
        ])
       RAW_LIBPNG=${LIBPNG_LIBADD}
   fi
   if test "x${have_libpng}" = "xno" ; then
      AC_MSG_WARN([Can not find a usuable PNG library. Make sure that CPPFLAGS and LDFLAGS are correctly set.])
   fi
fi

if test "x${have_libpng}" = "xyes" ; then
   M4RI_HAVE_LIBPNG=1
   AC_SUBST(M4RI_HAVE_LIBPNG)
   AC_SUBST(LIBPNG_LIBADD)
   AC_SUBST(M4RI_USE_PNG_PC)
   AC_SUBST(RAW_LIBPNG)
else
   M4RI_HAVE_LIBPNG=0
   AC_SUBST(M4RI_HAVE_LIBPNG)
fi

M4RI_BUILDING_M4RI=1
AC_DEFINE([M4RI_BUILDING_M4RI], [1], [Define to indicate that m4ri is being built instead of being used])

RELEASE="AC_PACKAGE_VERSION"
AC_SUBST(RELEASE)

AC_PROG_MAKE_SET

AC_CONFIG_FILES([Makefile tests/Makefile bench/Makefile m4ri/m4ri_config.h m4ri.pc])
AC_OUTPUT

