include ../../make.inc

.SECONDARY:

SRC = $(wildcard *.c)
TESTS = $(SRC:%.c=%.pass)

test: $(TESTS)

%.pass: %.x
	@ $^
	@ touch  $@

%.x: %.o LAPACK_%.o ../../liblarpack.a
	@ $(CC) $(CFLAGS) $^ -o $@ -lgfortran ../../liblarpack.a $(LINK_LAPACK) $(LINK_BLAS)

%.o: %.c
	@ $(CC) $(CFLAGS) -c $^ -o $@

%.o: %.f
	@ $(FC) $(FFLAGS) -c $^ -o $@

clean:
	rm -f *.pass *.o *.x
