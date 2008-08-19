CFLAGS=-I../src 
LDFLAGS=-L../.libs/ -lm4ri
DEBUG=-ggdb



PRGS=test_elimination test_multiplication bench_elimination bench_multiplication bench_addition test_trsm

CPUCYCLES_DIR=./cpucycles-20060326

all: cpucycles.o $(PRGS)

# cpucycles.o should be a dependency here
.c:	
	$(CC) $(DEBUG) $(CFLAGS) $@.c cpucycles.o $(LDFLAGS) -o $@
clean:
	rm -f $(PRGS)
	rm -f cpucycles.o cpucycles.h
	rm -f $(CPUCYCLES_DIR)/cpucycles.o
	rm -f $(CPUCYCLES_DIR)/cpucycles.h

cpucycles.h: cpucycles.o

cpucycles.o:
	(cd $(CPUCYCLES_DIR); sh do; cp cpucycles.o ..; cp cpucycles.h ..; cd ..)

