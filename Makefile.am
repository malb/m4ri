AUTOMAKE_OPTIONS = gnu
ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS=${SIMD_CFLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}

lib_LTLIBRARIES = libm4ri.la

libm4ri_la_SOURCES = \
	src/brilliantrussian.c \
	src/misc.c \
	src/mzd.c \
	src/graycode.c \
	src/strassen.c \
	src/mzp.c \
	src/triangular.c \
	src/triangular_russian.c \
	src/ple.c \
	src/ple_russian.c \
	src/solve.c \
	src/echelonform.c \
	src/mmc.c \
	src/debug_dump.c \
	src/io.c

BUILT_SOURCES = src/m4ri_config.h

pkgincludesubdir = $(includedir)/m4ri
pkgincludesub_HEADERS = src/m4ri.h \
	src/brilliantrussian.h \
	src/misc.h \
	src/mzd.h \
	src/graycode.h \
	src/strassen.h \
	src/parity.h \
	src/mzp.h \
	src/triangular.h \
	src/triangular_russian.h \
	src/ple.h \
	src/ple_russian.h \
	src/solve.h \
	src/echelonform.h \
	src/xor.h \
	src/mmc.h \
	src/debug_dump.h \
	src/io.h

nodist_pkgincludesub_HEADERS = src/m4ri_config.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = m4ri.pc

libm4ri_la_LDFLAGS = -release 0.0.20120615 -no-undefined
libm4ri_la_LIBADD = -lm $(LIBPNG_LIBADD)

check_PROGRAMS=test_multiplication test_elimination test_trsm test_ple test_solve test_kernel test_random test_smallops test_transpose test_colswap test_invert test_misc
test_multiplication_SOURCES=testsuite/test_multiplication.c
test_multiplication_LDFLAGS=-lm4ri -lm
test_multiplication_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

test_elimination_SOURCES=testsuite/test_elimination.c
test_elimination_LDFLAGS=-lm4ri -lm
test_elimination_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

test_trsm_SOURCES=testsuite/test_trsm.c
test_trsm_LDFLAGS=-lm4ri -lm
test_trsm_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

test_ple_SOURCES=testsuite/test_ple.c
test_ple_LDFLAGS=-lm4ri -lm
test_ple_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

test_solve_SOURCES=testsuite/test_solve.c
test_solve_LDFLAGS=-lm4ri -lm
test_solve_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

test_kernel_SOURCES=testsuite/test_kernel.c
test_kernel_LDFLAGS=-lm4ri -lm
test_kernel_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

test_random_SOURCES=testsuite/test_random.c
test_random_LDFLAGS=-lm4ri -lm
test_random_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

test_smallops_SOURCES=testsuite/test_smallops.c testsuite/testing.c testsuite/testing.h
test_smallops_LDFLAGS=-lm4ri -lm
test_smallops_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

test_transpose_SOURCES=testsuite/test_transpose.c
test_transpose_LDFLAGS=-lm4ri -lm
test_transpose_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

test_colswap_SOURCES=testsuite/test_colswap.c
test_colswap_LDFLAGS=-lm4ri -lm
test_colswap_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

test_invert_SOURCES=testsuite/test_invert.c
test_invert_LDFLAGS=-lm4ri -lm
test_invert_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

test_misc_SOURCES=testsuite/test_misc.c
test_misc_LDFLAGS=-lm4ri -lm
test_misc_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

TESTS = test_multiplication test_elimination test_trsm test_ple test_solve test_kernel test_random test_smallops test_transpose test_colswap test_invert test_misc
