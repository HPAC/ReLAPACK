include make.inc

.PHONY: test

SRC = $(wildcard src/*.c)
OBJS = $(SRC:%.c=%.o)

$(LIBNAME): $(OBJS)
	$(AR) -r $@ $^

%.o: %.c config.h
	$(CC) $(CFLAGS) -c $< -o $@

test: $(LIBNAME)
	cd test; make

clean:
	rm -f $(LIBNAME) $(OBJS)
	cd test; make clean
