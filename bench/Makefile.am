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
AM_CFLAGS = -I$(TOPSRCDIR) -I$(TOPBUILDDIR) -D_XOPEN_SOURCE=600 @CFLAGS@ $(DEFINES) @OPENMP_CFLAGS@ @PAPI_CFLAGS@

STAGEDIR := $(realpath -s $(TOPBUILDDIR)/.libs)
AM_LDFLAGS = -L$(STAGEDIR) -Wl,-rpath,$(STAGEDIR) cpucycles.o -lm4ri -lm @PAPI_LDFLAGS@ @PAPI_LIBS@
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

CPUCYCLES_DIR = cpucycles-20060326

CPUCYCLES_EXTRA_DIST = $(srcdir)/cpucycles-20060326/*.c \
                       $(srcdir)/cpucycles-20060326/*.h \
                       $(srcdir)/cpucycles-20060326/compile \
                       $(srcdir)/cpucycles-20060326/do 

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

BUILT_SOURCES = cpucycles.h

cpucycles.h: cpucycles.o

cpucycles.o:
	 (if [ $(srcdir) != $(builddir) ]; then \
		cp -r $(srcdir)/$(CPUCYCLES_DIR) $(builddir); \
	        chmod +w -R $(CPUCYCLES_DIR)/; \
         fi; \
         cd $(CPUCYCLES_DIR); \
         sh do; \
         cp cpucycles.o ..; \
         cp cpucycles.h ..; \
         cd ..; \
	 if [ $(srcdir) != $(builddir) ]; then \
		rm -r $(CPUCYCLES_DIR); \
         fi)

distclean-local:
	-rm -f cpucycles.h
