# TOPSRCDIR allows to compile this testsuite for a clone (and different revision) of this
# repository. If the environment variable topsrcdir is the top source directory of the clone,
# compile this testsuite as:
#
# make TOPSRCDIR="$topsrcdir"
#
# Finally, if you also leave out TOPSRCDIR and just run 'make' then the testsuite
# is compiled against the current source tree.

AUTOMAKE_OPTIONS = foreign
TOPSRCDIR = ..
DEFINES =
AM_CFLAGS = -I$(TOPSRCDIR) -D_XOPEN_SOURCE=600 @CFLAGS@ $(DEFINES) @OPENMP_CFLAGS@ @PAPI_CFLAGS@

STAGEDIR := $(realpath -s $(TOPSRCDIR)/.libs)
AM_LDFLAGS = -L$(STAGEDIR) -Wl,-rpath,$(STAGEDIR) cpucycles.o -lm4ri -lm @PAPI_LDFLAGS@ @PAPI_LIBS@
LDFLAGS = -no-install

BENCH = bench_elimination \
	bench_multiplication \
	bench_ple \
	bench_trsm \
	bench_elimination_sparse \
	bench_mzd \
	bench_invert \
	bench_rank

CPUCYCLES_DIR = ./cpucycles-20060326

bin_PROGRAMS = $(BENCH)

bench_elimination_SOURCES = bench_elimination.c benchmarking.c 
bench_multiplication_SOURCES = bench_multiplication.c benchmarking.c
bench_ple_SOURCES = bench_ple.c benchmarking.c
bench_trsm_SOURCES = bench_trsm.c benchmarking.c
bench_elimination_sparse_SOURCES = bench_elimination_sparse.c benchmarking.c
bench_mzd_SOURCES = bench_mzd.c benchmarking.c
bench_invert_SOURCES = bench_invert.c benchmarking.c
bench_rank_SOURCES = bench_rank.c benchmarking.c

BUILT_SOURCES = cpucycles.h

cpucycles.h: cpucycles.o

cpucycles.o:
	(cd $(CPUCYCLES_DIR); sh do; cp cpucycles.o ..; cp cpucycles.h ..)

