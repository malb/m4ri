CFLAGS=-I.. -std=c99
LDFLAGS=-L../.libs/ -lm4ri
DEBUG=-ggdb

TEST_PRGS=test_elimination test_multiplication test_trsm test_lqup test_solve
BENCH_PRGS=bench_elimination \
	bench_multiplication \
	bench_lqup \
	bench_trsm_lowerleft \
	bench_trsm_upperleft \
	bench_trsm_lowerright \
	bench_trsm_upperright

CPUCYCLES_DIR=./cpucycles-20060326

all: $(TEST_PRGS) $(BENCH_PRGS)

.c:	
	$(CC) $(DEBUG) $(CFLAGS) $@.c $(LDFLAGS) -o $@

bench_elimination: cpucycles.o bench_elimination.c
	$(CC) $(DEBUG) $(CFLAGS) bench_elimination.c cpucycles.o $(LDFLAGS) -o $@

bench_multiplication: cpucycles.o bench_multiplication.c
	$(CC) $(DEBUG) $(CFLAGS) bench_multiplication.c cpucycles.o $(LDFLAGS) -o $@

bench_lqup: cpucycles.o bench_lqup.c
	$(CC) $(DEBUG) $(CFLAGS) bench_lqup.c cpucycles.o $(LDFLAGS) -o $@

bench_trsm_lowerleft: cpucycles.o bench_trsm_lowerleft.c
	$(CC) $(DEBUG) $(CFLAGS) $@.c cpucycles.o $(LDFLAGS) -o $@

bench_trsm_lowerright: cpucycles.o bench_trsm_lowerright.c
	$(CC) $(DEBUG) $(CFLAGS) $@.c cpucycles.o $(LDFLAGS) -o $@

bench_trsm_upperleft: cpucycles.o bench_trsm_upperleft.c
	$(CC) $(DEBUG) $(CFLAGS) $@.c cpucycles.o $(LDFLAGS) -o $@

bench_trsm_upperright: cpucycles.o bench_trsm_upperright.c
	$(CC) $(DEBUG) $(CFLAGS) $@.c cpucycles.o $(LDFLAGS) -o $@


clean:
	rm -f $(TEST_PRGS)
	rm -f $(BENCH_PRGS)
	rm -f *.o cpucycles.h
	rm -f $(CPUCYCLES_DIR)/cpucycles.o
	rm -f $(CPUCYCLES_DIR)/cpucycles.h


cpucycles.h: cpucycles.o

cpucycles.o:
	(cd $(CPUCYCLES_DIR); sh do; cp cpucycles.o ..; cp cpucycles.h ..; cd ..)

