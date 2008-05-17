AUTOMAKE_OPTIONS = gnu

AM_CFLAGS=${C99_FLAG} ${SIMD_FLAGS}

lib_LTLIBRARIES = libm4ri.la

libm4ri_la_SOURCES = src/brilliantrussian.c src/misc.c src/packedmatrix.c src/grayflex.c src/strassen.c src/watch.c

include_HEADERS = src/m4ri.h src/brilliantrussian.h src/misc.h src/packedmatrix.h src/grayflex.h src/watch.h src/strassen.h src/parity.h

libm4ri_la_LDFLAGS = -version-info 0:0:0
