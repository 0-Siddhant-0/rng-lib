CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99 -I./include
LDFLAGS = -lm

all: librng.a test_rng

librng.a: src/rng.o
	ar rcs $@ $^

src/rng.o: src/rng.c include/rng.h
	$(CC) $(CFLAGS) -c $< -o $@

test_rng: src/test_rng.o librng.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

src/test_rng.o: src/test_rng.c include/rng.h
	$(CC) $(CFLAGS) -c $< -o $@

test: test_rng
	./test_rng

clean:
	rm -f src/*.o *.a test_rng
