AUTOMAKE_OPTIONS = foreign subdir-objects
ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS=${SIMD_CFLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}

lib_LTLIBRARIES = libm4ri.la

libm4ri_la_SOURCES = \
	m4ri/brilliantrussian.c \
	m4ri/misc.c \
	m4ri/mzd.c \
	m4ri/graycode.c \
	m4ri/strassen.c \
	m4ri/mzp.c \
	m4ri/triangular.c \
	m4ri/triangular_russian.c \
	m4ri/ple.c \
	m4ri/ple_russian.c \
	m4ri/solve.c \
	m4ri/echelonform.c \
	m4ri/mmc.c \
	m4ri/debug_dump.c \
	m4ri/io.c \
	m4ri/djb.c \
	m4ri/mp.c	

BUILT_SOURCES = m4ri/m4ri_config.h

pkgincludesubdir = $(includedir)/m4ri
pkgincludesub_HEADERS = m4ri/m4ri.h \
	m4ri/brilliantrussian.h \
	m4ri/misc.h \
	m4ri/mzd.h \
	m4ri/graycode.h \
	m4ri/strassen.h \
	m4ri/parity.h \
	m4ri/mzp.h \
	m4ri/triangular.h \
	m4ri/triangular_russian.h \
	m4ri/ple.h \
	m4ri/ple_russian.h \
	m4ri/ple_russian_template.h \
	m4ri/solve.h \
	m4ri/echelonform.h \
	m4ri/xor.h \
	m4ri/xor_template.h \
	m4ri/mmc.h \
	m4ri/debug_dump.h \
	m4ri/io.h \
	m4ri/djb.h \
	m4ri/mp.h

nodist_pkgincludesub_HEADERS = m4ri/m4ri_config.h

EXTRA_DIST=m4ri/Doxyfile

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = m4ri.pc

libm4ri_la_LDFLAGS = -release 0.0.$(RELEASE) -no-undefined
libm4ri_la_LIBADD = $(LIBPNG_LIBADD)

check_PROGRAMS=test_multiplication test_elimination test_trsm test_ple test_solve test_kernel test_random test_smallops test_transpose test_colswap test_invert test_misc
test_multiplication_SOURCES=testsuite/test_multiplication.c
test_multiplication_LDFLAGS=-lm4ri -lm
test_multiplication_CFLAGS=$(AM_CFLAGS)

test_elimination_SOURCES=testsuite/test_elimination.c
test_elimination_LDFLAGS=-lm4ri -lm
test_elimination_CFLAGS=$(AM_CFLAGS)

test_trsm_SOURCES=testsuite/test_trsm.c
test_trsm_LDFLAGS=-lm4ri -lm
test_trsm_CFLAGS=$(AM_CFLAGS)

test_ple_SOURCES=testsuite/test_ple.c
test_ple_LDFLAGS=-lm4ri -lm
test_ple_CFLAGS=$(AM_CFLAGS)

test_solve_SOURCES=testsuite/test_solve.c
test_solve_LDFLAGS=-lm4ri -lm
test_solve_CFLAGS=$(AM_CFLAGS)

test_kernel_SOURCES=testsuite/test_kernel.c
test_kernel_LDFLAGS=-lm4ri -lm
test_kernel_CFLAGS=$(AM_CFLAGS)

test_random_SOURCES=testsuite/test_random.c
test_random_LDFLAGS=-lm4ri -lm
test_random_CFLAGS=$(AM_CFLAGS)

test_smallops_SOURCES=testsuite/test_smallops.c testsuite/testing.c testsuite/testing.h
test_smallops_LDFLAGS=-lm4ri -lm
test_smallops_CFLAGS=$(AM_CFLAGS)

test_transpose_SOURCES=testsuite/test_transpose.c
test_transpose_LDFLAGS=-lm4ri -lm
test_transpose_CFLAGS=$(AM_CFLAGS)

test_colswap_SOURCES=testsuite/test_colswap.c
test_colswap_LDFLAGS=-lm4ri -lm
test_colswap_CFLAGS=$(AM_CFLAGS)

test_invert_SOURCES=testsuite/test_invert.c
test_invert_LDFLAGS=-lm4ri -lm
test_invert_CFLAGS=$(AM_CFLAGS)

test_misc_SOURCES=testsuite/test_misc.c
test_misc_LDFLAGS=-lm4ri -lm
test_misc_CFLAGS=$(AM_CFLAGS)

TESTS = test_multiplication test_elimination test_trsm test_ple test_solve test_kernel test_random test_smallops test_transpose test_colswap test_invert test_misc

