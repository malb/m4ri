CFLAGS=-I../
LDFLAGS=-L../../.libs/ -lm4ri
DEBUG=-ggdb

PRGS=test_reduction test_multiplication bench_reduction bench_multiplication

all: $(PRGS)

.c:
	$(CC) $(DEBUG) $(CFLAGS) $@.c $(LDFLAGS) -o $@
clean:
	rm $(PRGS)

