include ../make.inc

ROUTINES = $(shell echo {s,d,c,z}{lauum,potrf,sytrf,getrf,trsyl,tgsyl,gemm_tr} {s,d}sygst {c,z}{hegst,hetrf})

TESTS = $(ROUTINES:%=%.pass)  # these are just dummies
EXES = $(ROTUINES:%=%.x)

.SECONDARY: $(EXES)

test: util.o $(TESTS)
	@echo "passed all tests"

%.pass: %.x
	@echo -n $*:
	@./$< > /dev/null && echo " pass" || (echo " FAIL" && ./$<)

s%.x: x%.c util.o ../$(LIBNAME) config.h test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=s $< util.o -o $@ $(LINK_TEST) ../librelapack.a

d%.x: x%.c util.o ../$(LIBNAME) config.h test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=d $< util.o -o $@ $(LINK_TEST) ../librelapack.a

c%.x: x%.c util.o ../$(LIBNAME) config.h test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=c $< util.o -o $@ $(LINK_TEST) ../librelapack.a

z%.x: x%.c util.o ../$(LIBNAME) config.h test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=z $< util.o -o $@ $(LINK_TEST) ../librelapack.a

clean:
	rm -f *.x
