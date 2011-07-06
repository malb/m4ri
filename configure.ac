AC_INIT(libm4ri,1.0.1)

AC_PREREQ([2.64])

AC_CANONICAL_HOST

AC_CONFIG_SRCDIR(src/brilliantrussian.c)

AM_INIT_AUTOMAKE

dnl Include maintainer mode targets.
AM_MAINTAINER_MODE

dnl Needed when reconfiguring with 'autoreconf -i -s'
AC_CONFIG_MACRO_DIR([m4])

dnl Compiling with per-target flags (test_elimination.c) requires AM_PROG_CC_C_O.
AM_PROG_CC_C_O

AC_PROG_LIBTOOL

AC_PROG_INSTALL

AC_CONFIG_HEADERS(src/config.h)

dnl Check if a C++ compiler was specified. If so, assume we want to wrap word in a C++ class.
AC_EGREP_CPP(YES, [
#ifdef __cplusplus
YES
#endif
], [m4ri_wrapword="yes"], [m4ri_wrapword="no"])

if test "$m4ri_wrapword" = "yes"; then
  AC_DEFINE([M4RI_WRAPWORD], [/**/],
      [Defined when compiling with a C++ compiler and word should be a C++ class.])
else
  dnl We can only run this test when we're using a C compiler.
  AC_PROG_CC_C99()
  if test "$ac_cv_prog_cc_c99" = "no"; then
    AC_MSG_ERROR([C99 support is required but not found.])
  fi
fi

# SSE2 support
AC_ARG_ENABLE([sse2],
        AS_HELP_STRING([--disable-sse2], [don't use SSE2 instruction set.]),
	, [if test "$m4ri_wrapword" = "yes"; then enable_sse2="no"; else enable_sse2="yes"; fi])

AS_IF([test "x$enable_sse2" != "xno"], [
   if test "$m4ri_wrapword" = "yes"; then
      AC_MSG_ERROR([SSE2 cannot be supported when wrapping word in a C++ class.])
   fi
   case $host_cpu in i[[3456]]86*|x86_64*)
          AX_CPU_VENDOR()
          if test "x$ax_cv_cpu_vendor" = "xIntel"; then
             AX_EXT() # SSE2 is slower on the Opteron
          fi
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
        AS_HELP_STRING([--with-cachesize@<:@=VALUE@:>@], [Manual cache sizes, separated by a colon. Overrides cache tuning.]),[m4ri_config_cachesize=$withval])

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

# Debugging support
AC_ARG_ENABLE([debug],
	AS_HELP_STRING([--enable-debug], [Enable assert() statements for debugging.]))

AC_ARG_ENABLE([debug-dump],
        AS_HELP_STRING([--enable-debug-dump], [Dump output at exit of every function.]))

if test "x$enable_debug_dump" = "xyes"; then
  M4RI_DEBUG_DUMP=1
else
  M4RI_DEBUG_DUMP=0
fi
AC_SUBST(M4RI_DEBUG_DUMP)

AC_ARG_ENABLE([debug-mzd],
        AS_HELP_STRING([--enable-debug-mzd], [Add consistency checks on matrix structures.]))

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
    AC_DEFINE(NDEBUG,1,[Define whether debugging is enabled])
  fi
fi

# For the testsuite. Detect if PAPI is installed. See http://icl.cs.utk.edu/papi/ .

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
     AS_HELP_STRING([--enable-cachetune],[calculate cache size from timing information.]))      

# L1 and L2 Cache Size
if test -z $m4ri_config_cachesize; then
   AX_CACHE_SIZE()
   AS_IF([test "x$enable_cachetune" = "xyes"], [AX_CACHE_SIZE_TUNE()]) 
else
  AS_IF([test "x$enable_cachetune" = "xyes"], [AC_MSG_WARN(Ignoring cache tuning since --with-cachesize was given.)]) 

  ax_l1_size=`echo $m4ri_config_cachesize | $SED 's/\:.*//g'`
  ax_l2_size=`echo $m4ri_config_cachesize | $SED 's/.*\://g'`

  M4RI_CPU_L1_CACHE=${ax_l1_size}
  M4RI_CPU_L2_CACHE=${ax_l2_size}
  AC_SUBST(M4RI_CPU_L1_CACHE)
  AC_SUBST(M4RI_CPU_L2_CACHE)

fi


AC_PROG_MAKE_SET

AC_CONFIG_FILES([Makefile testsuite/Makefile src/m4ri_config.h])
AC_OUTPUT

