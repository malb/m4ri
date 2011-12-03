AUTOMAKE_OPTIONS = gnu
ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS=${SIMD_CFLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}

lib_LTLIBRARIES = libm4ri.la

libm4ri_la_SOURCES = \
	src/brilliantrussian.c \
	src/misc.c \
	src/packedmatrix.c \
	src/grayflex.c \
	src/strassen.c \
	src/permutation.c \
	src/trsm.c \
	src/pls.c \
	src/solve.c \
	src/pls_mmpf.c \
	src/echelonform.c \
	src/mmc.c \
	src/debug_dump.c \
	src/io.c

BUILT_SOURCES = src/m4ri_config.h

pkgincludesubdir = $(includedir)/m4ri
pkgincludesub_HEADERS = src/m4ri.h \
	src/brilliantrussian.h \
	src/misc.h \
	src/packedmatrix.h \
	src/grayflex.h \
	src/strassen.h \
	src/parity.h \
	src/permutation.h \
	src/trsm.h \
	src/pls.h \
	src/solve.h \
	src/pls_mmpf.h \
	src/echelonform.h \
	src/xor.h \
	src/mmc.h \
	src/debug_dump.h \
	src/io.h

nodist_pkgincludesub_HEADERS = src/m4ri_config.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = m4ri.pc

libm4ri_la_LDFLAGS = -release 0.0.20111203 -no-undefined
libm4ri_la_LIBADD = -lm $(LIBPNG_LIBADD)

check_PROGRAMS=test_multiplication test_elimination test_trsm test_pls test_solve test_kernel test_random test_smallops test_transpose test_colswap
test_multiplication_SOURCES=testsuite/test_multiplication.c
test_multiplication_LDFLAGS=-lm4ri -lm
test_multiplication_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

test_elimination_SOURCES=testsuite/test_elimination.c
test_elimination_LDFLAGS=-lm4ri -lm
test_elimination_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

test_trsm_SOURCES=testsuite/test_trsm.c
test_trsm_LDFLAGS=-lm4ri -lm
test_trsm_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

test_pls_SOURCES=testsuite/test_pluq.c
test_pls_LDFLAGS=-lm4ri -lm
test_pls_CFLAGS=-I$(srcdir)/src $(AM_CFLAGS)

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

TESTS = test_multiplication test_elimination test_trsm test_pls test_solve test_kernel test_random test_smallops test_transpose test_colswap
