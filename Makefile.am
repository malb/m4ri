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
libm4ri_la_LIBADD = $(LIBPNG_LIBADD) $(LIBM)

SUBDIRS = . tests
DIST_SUBDIRS = $(SUBDIRS) bench


CLANGFORMAT       ?= clang-format
.PHONY: check-style
check-style:
	$(CLANGFORMAT) -i --style=file m4ri/*.{c,h} tests/*.{c,h} bench/*.{c,h}
