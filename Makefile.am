AUTOMAKE_OPTIONS = gnu

AM_CFLAGS=${C99_FLAG} ${SIMD_FLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}

lib_LTLIBRARIES = libm4ri.la

libm4ri_la_SOURCES = src/brilliantrussian.c src/misc.c src/packedmatrix.c src/grayflex.c src/strassen.c src/permutation.c src/trsm.c src/lqup.c

include_HEADERS = src/m4ri.h src/brilliantrussian.h src/misc.h src/packedmatrix.h src/grayflex.h src/strassen.h src/parity.h src/permutation.h src/config.h src/trsm.h src/lqup.h

#libm4ri_la_LDFLAGS = -version-info 0:0:0
libm4ri_la_LDFLAGS = -release 0.0.20080901
