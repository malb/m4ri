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

m4_include([m4/ax_gcc_x86_cpuid.m4])
m4_include([m4/ax_ext.m4])dnl
m4_include([m4/ax_cpu_vendor.m4])dnl
m4_include([m4/ax_cache_size.m4])dnl
m4_include([m4/ax_cache_size_tune.m4])dnl
m4_include([m4/ax_check_compiler_flags.m4])dnl
m4_include([m4/ax_openmp.m4])dnl

AC_PROG_CC_C99()

if test "$ac_cv_prog_cc_c99" = "no"; then
   AC_MSG_ERROR([C99 support is required but not found.])
fi

# SSE2 support
AC_ARG_ENABLE([sse2],
        AS_HELP_STRING([--disable-sse2], [don't use SSE2 instruction set.]))

AS_IF([test "x$enable_sse2" != "xno"], [
   case $host_cpu in i[[3456]]86*|x86_64*)
          AX_CPU_VENDOR()
          if test "x$ax_cv_cpu_vendor" = "xIntel"; then
             AX_EXT() # SSE2 is slower on the Opteron
          fi
   esac
])

AC_CHECK_HEADER([mm_malloc.h],AC_DEFINE(HAVE_MM_MALLOC,,[Support aligned allocations]),)

# OpenMP support
AC_ARG_ENABLE([openmp],
        AS_HELP_STRING( [--enable-openmp],[add support for OpenMP multicore support.]))

AS_IF([test "x$enable_openmp" = "xyes"], [
   AX_OPENMP()
])
AC_SUBST(OPENMP_CFLAGS)

# Debugging support
AC_ARG_ENABLE(debug, [  --enable-debug          Enable assert() statements for debugging.])

if test "x$enable_debug" = x"yes"; then
   DEBUG_FLAGS="-g"
   AC_SUBST(DEBUG_FLAGS)
else
   AC_DEFINE(NDEBUG,1,[Define whether debugging is enabled])
fi

# L1 and L2 Cache Size
AX_CACHE_SIZE()
AC_ARG_ENABLE([cachetune], 
        AS_HELP_STRING([--enable-cachetune],[calculate cache size from timing information.]))
AS_IF([test "x$enable_cachetune" = "xyes"], [
    AX_CACHE_SIZE_TUNE()
])

AC_PROG_MAKE_SET

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

