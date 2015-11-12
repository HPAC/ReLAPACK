include ../make.inc

# routine shorthand expansion
ROUTINES := $(ROUTINES:all=xlauum xsygst xtrtri xpotrf xgetrf xtrsyl)
ROUTINES := $(ROUTINES:xlauum=slauum dlauum clauum zlauum)
ROUTINES := $(ROUTINES:xsygst=ssygst dsygst chegst zhegst)
ROUTINES := $(ROUTINES:xtrtri=strtri dtrtri ctrtri ztrtri)
ROUTINES := $(ROUTINES:xpotrf=spotrf dpotrf cpotrf zpotrf)
ROUTINES := $(ROUTINES:xgetrf=sgetrf dgetrf cgetrf zgetrf)
ROUTINES := $(ROUTINES:xtrsyl=strsyl dtrsyl ctrsyl ztrsyl)
# sort and remove duplicates
ROUTINES := $(sort $(ROUTINES))

TESTS = $(ROUTINES:%=%.pass)  # these are just dummies
EXES = $(ROTUINES:%=%.x)

.SECONDARY: $(EXES)

test: util.o $(TESTS)

%.pass: %.x
	@ ./$< $(TEST_SIZE)

%.x: %.c util.o
	$(CC) $(FLAGS) $< -o $@ util.o ../librelapack.a $(LINK_LAPACK) $(LINK_BLAS)

clean:
	rm -f *.x
