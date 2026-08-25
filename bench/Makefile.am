# TOPSRCDIR allows to compile this testsuite for a clone (and different revision) of this
# repository. If the environment variable topsrcdir is the top source directory of the clone,
# compile this testsuite as:
#
# make TOPSRCDIR="$topsrcdir"
#
# Finally, if you also leave out TOPSRCDIR and just run 'make' then the testsuite
# is compiled against the current source tree.

AUTOMAKE_OPTIONS = foreign
TOPSRCDIR = $(srcdir)/..
TOPBUILDDIR = $(builddir)/..
DEFINES =
# include TOPBUILDIR for m4ri_config.h;
# include $(srcdir)/$(CPUCYCLES_DIR)/cpucycles for cpucycles.h
AM_CFLAGS = -I$(TOPSRCDIR) -I$(TOPBUILDDIR) -I$(srcdir)/$(CPUCYCLES_DIR)/cpucycles -D_XOPEN_SOURCE=600 @CFLAGS@ $(DEFINES) @OPENMP_CFLAGS@ @PAPI_CFLAGS@ @LIBPNG_CFLAGS@

# Expanded lazily (with '='): $(TOPBUILDDIR)/.libs does not exist until the
# m4ri library has been built, and the library prerequisite below may create
# it during this very make run.
STAGEDIR = $(realpath -s $(TOPBUILDDIR)/.libs)
AM_LDFLAGS = -L$(STAGEDIR) -Wl,-rpath,$(STAGEDIR) libcpucycles.a -lm4ri -lm @M4RI_RT_LIBS@ @PAPI_LDFLAGS@ @PAPI_LIBS@
LDFLAGS = -no-install

# The benchmarks link against the m4ri shared library that libtool stages in
# $(TOPBUILDDIR)/.libs. Build it first if it is missing (fresh configure or
# after 'make clean'); without it the -L/-rpath above would be empty.
$(TOPBUILDDIR)/libm4ri.la:
	$(MAKE) -C $(TOPBUILDDIR) $(AM_MAKEFLAGS) libm4ri.la

$(BENCH): $(TOPBUILDDIR)/libm4ri.la

BENCH = bench_elimination \
	bench_multiplication \
        bench_m4rm \
	bench_ple \
	bench_trsm \
	bench_elimination_sparse \
	bench_mzd \
	bench_invert \
	bench_rank

# libcpucycles is a cycle counter library that the benchmarks link against.
# It is built with its own ./configure (a python3 script) and make. Only the
# static archive build/0/package/lib/libcpucycles.a is built: the default
# target additionally builds a shared library and command binaries that need
# -lrt and symlink support, which are unavailable on Windows. The public
# header is cpucycles/cpucycles.h.
CPUCYCLES_DIR = libcpucycles-20260625

PYTHON3 = @PYTHON3@

# Files needed to run libcpucycles' ./configure and make from the
# distribution tarball.
CPUCYCLES_EXTRA_DIST = $(srcdir)/$(CPUCYCLES_DIR)/Makefile \
                       $(srcdir)/$(CPUCYCLES_DIR)/configure \
                       $(srcdir)/$(CPUCYCLES_DIR)/version \
                       $(srcdir)/$(CPUCYCLES_DIR)/compilers/default \
                       $(srcdir)/$(CPUCYCLES_DIR)/cpucycles/options \
                       $(wildcard $(srcdir)/$(CPUCYCLES_DIR)/cpucycles/*.c) \
                       $(wildcard $(srcdir)/$(CPUCYCLES_DIR)/cpucycles/*.h) \
                       $(wildcard $(srcdir)/$(CPUCYCLES_DIR)/command/*) \
                       $(wildcard $(srcdir)/$(CPUCYCLES_DIR)/scripts-build/*) \
                       $(wildcard $(srcdir)/$(CPUCYCLES_DIR)/doc/man/*) \
                       $(wildcard $(srcdir)/$(CPUCYCLES_DIR)/doc/*.md)

EXTRA_DIST = $(CPUCYCLES_EXTRA_DIST)

bin_PROGRAMS = $(BENCH)

bench_elimination_SOURCES = bench_elimination.c benchmarking.c benchmarking.h 
bench_multiplication_SOURCES = bench_multiplication.c benchmarking.c benchmarking.h
bench_m4rm_SOURCES = bench_m4rm.c benchmarking.c benchmarking.h
bench_ple_SOURCES = bench_ple.c benchmarking.c benchmarking.h
bench_trsm_SOURCES = bench_trsm.c benchmarking.c benchmarking.h
bench_elimination_sparse_SOURCES = bench_elimination_sparse.c benchmarking.c benchmarking.h
bench_mzd_SOURCES = bench_mzd.c benchmarking.c benchmarking.h
bench_invert_SOURCES = bench_invert.c benchmarking.c benchmarking.h
bench_rank_SOURCES = bench_rank.c benchmarking.c benchmarking.h

BUILT_SOURCES = libcpucycles.a

libcpucycles.a:
	 (if [ $(srcdir) != $(builddir) ]; then \
		cp -r $(srcdir)/$(CPUCYCLES_DIR) $(builddir); \
	        chmod -R u+w $(CPUCYCLES_DIR)/; \
         fi; \
         cd $(CPUCYCLES_DIR); \
         if [ -z "$(PYTHON3)" ]; then \
            echo "error: python3 is required to build $(CPUCYCLES_DIR) but was not found at configure time" 1>&2; \
            exit 1; \
         fi; \
         rm -rf build; \
         $(PYTHON3) configure; \
         make -C build/0 package/lib/libcpucycles.a; \
         cp build/0/package/lib/libcpucycles.a ..; \
         cd ..; \
	 if [ $(srcdir) != $(builddir) ]; then \
		rm -rf $(CPUCYCLES_DIR); \
         fi)

distclean-local:
	-rm -f libcpucycles.a
	-if [ "$(srcdir)" != "$(builddir)" ]; then \
		rm -rf $(CPUCYCLES_DIR); \
	else \
		rm -rf $(CPUCYCLES_DIR)/build; \
	fi
