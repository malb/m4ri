AUTOMAKE_OPTIONS = gnu

AM_CFLAGS=${C99_FLAG} ${SIMD_FLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}

lib_LTLIBRARIES = libm4ri.la

libm4ri_la_SOURCES = src/brilliantrussian.c src/misc.c src/packedmatrix.c src/grayflex.c src/strassen.c src/permutation.c src/trsm.c src/lqup.c

pkgincludesubdir = $(includedir)/m4ri
pkgincludesub_HEADERS = src/m4ri.h src/brilliantrussian.h src/misc.h src/packedmatrix.h src/grayflex.h src/strassen.h src/parity.h src/permutation.h src/config.h src/trsm.h src/lqup.h

#libm4ri_la_LDFLAGS = -version-info 0:0:0
libm4ri_la_LDFLAGS = -release 0.0.20081029

check_PROGRAMS=test_multiplication test_elimination test_trsm
test_multiplication_SOURCES=testsuite/test_multiplication.c
test_multiplication_LDFLAGS=-lm4ri -lm

test_elimination_SOURCES=testsuite/test_elimination.c
test_elimination_LDFLAGS=-lm4ri -lm

test_trsm_SOURCES=testsuite/test_trsm.c
test_trsm_LDFLAGS=-lm4ri -lm

TESTS = test_multiplication test_elimination test_trsm
