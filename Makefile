include make.inc

.PHONY: test

LIB = liblarpack.a

LAUUM_SRC = $(wildcard src/lauum/*.c)
TRTRI_SRC = $(wildcard src/trtri/*.c)
POTRF_SRC = $(wildcard src/potrf/*.c)
GETRF_SRC = $(wildcard src/getrf/*.c)

SRC = $(LAUUM_SRC) $(TRTRI_SRC) $(POTRF_SRC) $(GETRF_SRC)

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
