include make.inc

.PHONY: test

LIB = liblarpack.a

SRC = $(wildcard src/*.c)
# LAUUM_SRC = $(wildcard src/lauum/*.c)
# TRTRI_SRC = $(wildcard src/trtri/*.c)
# POTRF_SRC = $(wildcard src/potrf/*.c)
# GETRF_SRC = $(wildcard src/getrf/*.c)
# TRSYL_SRC = $(wildcard src/trsyl/*.c)
# SYGST_SRC = $(wildcard src/sygst/*.c)

OBJS = $(SRC:%.c=%.o)

$(LIB): $(OBJS)
	$(AR) -r $@ $^

test: $(LIB) 
	cd test; make

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	rm -f $(LIB) $(OBJS)
	cd test; make clean
