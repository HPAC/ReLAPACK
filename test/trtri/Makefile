include ../../make.inc

.SECONDARY:

SRC = $(wildcard *.c)
TESTS = $(SRC:%.c=%.pass)
UTIL = ../util.o

test: $(TESTS)

%.pass: %.x
	@ $^

%.x: %.o ../../liblarpack.a
	@ $(CC) $(CFLAGS) $^ -o $@ $(UTIL) $(LINK_LAPACK) $(LINK_BLAS)

%.o: %.c
	@ $(CC) $(CFLAGS) -c $^ -o $@

%.o: %.f
	@ $(FC) $(FFLAGS) -c $^ -o $@

clean:
	rm -f *.pass *.o *.x
