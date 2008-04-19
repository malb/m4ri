INCLUDE=-I../
LIBS=-L../../.libs/ -lm4ri

all: test_reduction test_multiplication

test_reduction: test_reduction.c
	$(CC) -ggdb $(INCLUDE) $(LIBS) test_reduction.c -o test_reduction

test_multiplication: test_multiplication.c
	$(CC) -ggdb $(INCLUDE) $(LIBS) test_multiplication.c -o test_multiplication

