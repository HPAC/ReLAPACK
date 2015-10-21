include make.inc

.PHONY: test

LIB = liblarpack.a

TRTRI_SRC = $(wildcard src/trtri/*.c)
POTRF_SRC = $(wildcard src/potrf/*.c)
GETRF_SRC = $(wildcard src/getrf/*.c)

SRC = $(TRTRI_SRC) $(POTRF_SRC) $(GETRF_SRC)

OBJS = $(SRC:%.c=%.o)

$(LIB): $(OBJS)
	$(AR) -r $@ $^

test: $(LIB) 
	cd test; make

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	rm $(LIB) src/*/*.o
	cd test; make clean
