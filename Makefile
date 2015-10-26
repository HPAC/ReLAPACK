include make.inc

.PHONY: test

LIB = liblarpack.a

SRC = $(wildcard src/*.c) $(wildcard src/*.f)

TMP = $(SRC:%.c=%.o)
OBJS = $(TMP:%.f=%.o)

$(LIB): $(OBJS)
	$(AR) -r $@ $^

test: $(LIB) 
	cd test; make

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	rm -f $(LIB) $(OBJS)
	cd test; make clean
