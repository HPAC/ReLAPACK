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
# xsytrf auxiliary routine xgemm_tr
ROUTINES += $(ROUTINES:%sytrf=%gemm_tr) $(ROUTINES:%hetrf=%gemm_tr)
# sort and remove duplicates
ROUTINES := $(sort $(ROUTINES))

TESTS = $(ROUTINES:%=%.pass)  # these are just dummies
EXES = $(ROTUINES:%=%.x)

.SECONDARY: $(EXES)

test: util.o $(TESTS)

%.pass: %.x
	@ ./$< $(TEST_SIZE)

%.x: %.c util.o ../librelapack.a
	$(CC) $(CFLAGS) $^ -o $@ $(LINK_LAPACK) $(LINK_BLAS)

clean:
	rm -f *.x
