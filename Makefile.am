AUTOMAKE_OPTIONS = gnu

AM_CFLAGS=-std=c99 ${SIMD_FLAGS}

lib_LTLIBRARIES = libm4ri.la

libm4ri_la_SOURCES = src/brilliantrussian.c src/misc.c src/packedmatrix.c src/grayflex.c src/strassen.c src/watch.c

include_HEADERS = src/m4ri.h src/brilliantrussian.h src/misc.h src/packedmatrix.h src/grayflex.h src/watch.h src/strassen.h

libm4ri_la_LDFLAGS = -version-info 0:0:0
