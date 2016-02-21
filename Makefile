include make.inc

.PHONY: test

LIB = librelapack.a

CSRC = $(wildcard src/*.c)
FSRC = $(wildcard src/*.f)
OBJS = $(CSRC:%.c=%.o) $(FSRC:%.f=%.o)

OBJS = $(CSRC:%.c=%.o) $(FSRC:%.f=%.o)

$(LIB): $(OBJS)
	$(AR) -r $@ $^

%.o: %.c config.h
	$(CC) $(CFLAGS) -c $< -o $@

test: $(LIB)
	cd test; make

clean:
	rm -f $(LIB) $(OBJS)
	cd test; make clean
