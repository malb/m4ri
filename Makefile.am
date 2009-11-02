AUTOMAKE_OPTIONS = gnu

AM_CFLAGS=${C99_FLAG} ${SIMD_FLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}

lib_LTLIBRARIES = libm4ri.la

libm4ri_la_SOURCES = src/brilliantrussian.c src/misc.c src/packedmatrix.c src/grayflex.c src/strassen.c src/permutation.c src/trsm.c src/lqup.c src/solve.c src/pluq_mmpf.c

pkgincludesubdir = $(includedir)/m4ri
pkgincludesub_HEADERS = src/m4ri.h src/brilliantrussian.h src/misc.h src/packedmatrix.h src/grayflex.h src/strassen.h src/parity.h src/permutation.h src/config.h src/trsm.h src/lqup.h src/solve.h  src/pluq_mmpf.h

#libm4ri_la_LDFLAGS = -version-info 0:0:0
libm4ri_la_LDFLAGS = -release 0.0.20091101

check_PROGRAMS=test_multiplication test_elimination test_trsm test_lqup test_solve test_kernel
test_multiplication_SOURCES=testsuite/test_multiplication.c
test_multiplication_LDFLAGS=-lm4ri -lm

test_elimination_SOURCES=testsuite/test_elimination.c
test_elimination_LDFLAGS=-lm4ri -lm

test_trsm_SOURCES=testsuite/test_trsm.c
test_trsm_LDFLAGS=-lm4ri -lm

test_lqup_SOURCES=testsuite/test_lqup.c
test_lqup_LDFLAGS=-lm4ri -lm

test_solve_SOURCES=testsuite/test_solve.c
test_solve_LDFLAGS=-lm4ri -lm

test_kernel_SOURCES=testsuite/test_kernel.c
test_kernel_LDFLAGS=-lm4ri -lm

TESTS = test_multiplication test_elimination test_trsm test_lqup test_solve test_kernel
