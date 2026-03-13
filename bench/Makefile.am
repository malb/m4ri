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
# include TOPBUILDIR for m4ri_config.h
AM_CFLAGS = -I$(TOPSRCDIR) -I$(TOPBUILDDIR) -D_XOPEN_SOURCE=600 @CFLAGS@ $(DEFINES) @OPENMP_CFLAGS@ @PAPI_CFLAGS@ @LIBPNG_CFLAGS@

STAGEDIR := $(realpath -s $(TOPBUILDDIR)/.libs)
AM_LDFLAGS = -L$(STAGEDIR) -Wl,-rpath,$(STAGEDIR) @PAPI_LDFLAGS@
LDADD = libcpucycles.a -lm4ri -lm -lrt @PAPI_LIBS@
LDFLAGS = -no-install

BENCH = bench_elimination \
	bench_multiplication \
        bench_m4rm \
	bench_ple \
	bench_trsm \
	bench_elimination_sparse \
	bench_mzd \
	bench_invert \
	bench_rank

CPUCYCLES_DIR = libcpucycles-20260105

CPUCYCLES_EXTRA_DIST = $(srcdir)/$(CPUCYCLES_DIR)

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

BUILT_SOURCES = cpucycles.h libcpucycles.a

cpucycles.h: libcpucycles.a

libcpucycles.a:
	 (if [ $(srcdir) != $(builddir) ]; then \
		cp -r $(srcdir)/$(CPUCYCLES_DIR) $(builddir); \
	        chmod -R +w $(CPUCYCLES_DIR)/; \
         fi; \
         cd $(CPUCYCLES_DIR); \
		 ./configure; \
	 $(MAKE); \
	 cp build/0/package/lib/libcpucycles.a ..; \
	 cp build/0/package/include/cpucycles.h ..; \
         cd ..; \
	 if [ $(srcdir) != $(builddir) ]; then \
		rm -r $(CPUCYCLES_DIR); \
         fi)

distclean-local:
	-rm -f cpucycles.h libcpucycles.a
	-rm -rf $(CPUCYCLES_DIR)/build
