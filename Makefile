include make.inc

.PHONY: test

LIB = librelapack.a

# routine shorthand expansion
ROUTINES := $(ROUTINES:all=xlauum xsygst xtrtri xpotrf xsytrf xgetrf xtrsyl xtgsyl)
ROUTINES := $(ROUTINES:xlauum=slauum dlauum clauum zlauum)
ROUTINES := $(ROUTINES:xsygst=ssygst dsygst chegst zhegst)
ROUTINES := $(ROUTINES:xtrtri=strtri dtrtri ctrtri ztrtri)
ROUTINES := $(ROUTINES:xpotrf=spotrf dpotrf cpotrf zpotrf)
ROUTINES := $(ROUTINES:xsytrf=ssytrf dsytrf csytrf chetrf zsytrf zhetrf ssytrf_rook dsytrf_rook csytrf_rook chetrf_rook zsytrf_rook zhetrf_rook)
ROUTINES := $(ROUTINES:xgetrf=sgetrf dgetrf cgetrf zgetrf)
ROUTINES := $(ROUTINES:xtrsyl=strsyl dtrsyl ctrsyl ztrsyl)
ROUTINES := $(ROUTINES:xtgsyl=stgsyl dtgsyl ctgsyl ztgsyl)
# xtrsyl need unblocked xtrsy2
ROUTINES += $(ROUTINES:%trsyl=%trsy2)
# xsytrf[_rook] need xgemm_tr, and xsytrf[_rook]_rec
ROUTINES += $(ROUTINES:%sytrf=%sytrf_rec2) $(ROUTINES:%hetrf=%hetrf_rec2)
ROUTINES += $(ROUTINES:%sytrf_rook=%sytrf_rook_rec2) $(ROUTINES:%hetrf_rook=%hetrf_rook_rec2)
ROUTINES += $(ROUTINES:%sytrf=%gemm_tr) $(ROUTINES:%hetrf=%gemm_tr)
ROUTINES += $(ROUTINES:%sytrf_rook=%gemm_tr) $(ROUTINES:%hetrf_rook=%gemm_tr)
# sort and remove duplicates
ROUTINES := $(sort $(ROUTINES))

OBJS = $(ROUTINES:%=src/%.o)

$(LIB): $(OBJS)
	$(AR) -r $@ $^

%.o: %.c config.h
	$(CC) $(CFLAGS) -c $< -o $@

test: $(LIB) 
	cd test; make

clean:
	rm -f $(LIB) $(OBJS)
	cd test; make clean
