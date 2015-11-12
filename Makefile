include make.inc

.PHONY: test

LIB = librelapack.a

# routine shorthand expansion
ROUTINES := $(ROUTINES:all=xlauum xsygst xtrtri xpotrf xgetrf xtrsyl)
ROUTINES := $(ROUTINES:xlauum=slauum dlauum clauum zlauum)
ROUTINES := $(ROUTINES:xsygst=ssygst dsygst chegst zhegst)
ROUTINES := $(ROUTINES:xtrtri=strtri dtrtri ctrtri ztrtri)
ROUTINES := $(ROUTINES:xpotrf=spotrf dpotrf cpotrf zpotrf)
ROUTINES := $(ROUTINES:xgetrf=sgetrf dgetrf cgetrf zgetrf)
ROUTINES := $(ROUTINES:xtrsyl=strsyl dtrsyl ctrsyl ztrsyl)
# xtrsyl need unblocked xtrsy2
ROUTINES += $(ROUTINES:%trsyl=%trsy2)
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
