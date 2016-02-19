include ../make.inc

# routine shorthand expansion
ROUTINES := $(ROUTINES:all=xlauum xsygst xtrtri xpotrf xsytrf xgetrf xtrsyl xtgsyl)
ROUTINES := $(ROUTINES:xlauum=slauum dlauum clauum zlauum)
ROUTINES := $(ROUTINES:xsygst=ssygst dsygst chegst zhegst)
ROUTINES := $(ROUTINES:xtrtri=strtri dtrtri ctrtri ztrtri)
ROUTINES := $(ROUTINES:xpotrf=spotrf dpotrf cpotrf zpotrf)
ROUTINES := $(ROUTINES:xsytrf=ssytrf dsytrf csytrf chetrf zsytrf zhetrf)
ROUTINES := $(ROUTINES:xgetrf=sgetrf dgetrf cgetrf zgetrf)
ROUTINES := $(ROUTINES:xtrsyl=strsyl dtrsyl ctrsyl ztrsyl)
ROUTINES := $(ROUTINES:xtgsyl=stgsyl dtgsyl ctgsyl ztgsyl)
# xsytrf[_rook] auxiliary routine xgemm_tr
ROUTINES += $(ROUTINES:%sytrf=%gemm_tr) $(ROUTINES:%hetrf=%gemm_tr)
ROUTINES += $(ROUTINES:%sytrf_rook=%gemm_tr) $(ROUTINES:%hetrf_rook=%gemm_tr)
# sort and remove duplicates
ROUTINES := $(sort $(ROUTINES))

TESTS = $(ROUTINES:%=%.pass)  # these are just dummies
EXES = $(ROTUINES:%=%.x)

.SECONDARY: $(EXES)

test: util.o $(TESTS)

%.pass: %.x
	@echo -n $*:
	@./$< > /dev/null && echo " pass" || (echo " FAIL" && ./$<)

s%.x: x%.c util.o ../librelapack.a config.h test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=s $^ -o $@ $(LINK_LAPACK) $(LINK_BLAS)

d%.x: x%.c util.o ../librelapack.a config.h test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=d $^ -o $@ $(LINK_LAPACK) $(LINK_BLAS)

c%.x: x%.c util.o ../librelapack.a config.h test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=c $^ -o $@ $(LINK_LAPACK) $(LINK_BLAS)

z%.x: x%.c util.o ../librelapack.a config.h test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=z $^ -o $@ $(LINK_LAPACK) $(LINK_BLAS)

#%.x: %.c util.o ../librelapack.a
#	$(CC) $(CFLAGS) $^ -o $@ $(LINK_LAPACK) $(LINK_BLAS)

clean:
	rm -f *.x
