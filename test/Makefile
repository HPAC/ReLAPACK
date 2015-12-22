include ../make.inc

# routine shorthand expansion
ROUTINES := $(ROUTINES:all=xlauum xsygst xtrtri xpotrf xgetrf xtrsyl xsytrf)
ROUTINES := $(ROUTINES:xlauum=slauum dlauum clauum zlauum)
ROUTINES := $(ROUTINES:xsygst=ssygst dsygst chegst zhegst)
ROUTINES := $(ROUTINES:xtrtri=strtri dtrtri ctrtri ztrtri)
ROUTINES := $(ROUTINES:xpotrf=spotrf dpotrf cpotrf zpotrf)
ROUTINES := $(ROUTINES:xgetrf=sgetrf dgetrf cgetrf zgetrf)
ROUTINES := $(ROUTINES:xtrsyl=strsyl dtrsyl ctrsyl ztrsyl)
ROUTINES := $(ROUTINES:xsytrf=xsytrf xgemm_tr)
ROUTINES := $(ROUTINES:xsytrf=ssytrf)
ROUTINES := $(ROUTINES:xgemm_tr=sgemm_tr)
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
