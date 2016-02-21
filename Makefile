include make.inc

.PHONY: test

CSRC = $(wildcard src/*.c)
FSRC = $(wildcard src/*.f)
OBJS = $(CSRC:%.c=%.o) $(FSRC:%.f=%.o)

OBJS = $(CSRC:%.c=%.o) $(FSRC:%.f=%.o)

$(LIBNAME): $(OBJS)
	$(AR) -r $@ $^

%.o: %.c config.h
	$(CC) $(CFLAGS) -c $< -o $@

test: $(LIBNAME)
	cd test; make

clean:
	rm -f $(LIBNAME) $(OBJS)
	cd test; make clean
