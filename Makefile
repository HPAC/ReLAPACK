include make.inc


SRC = $(wildcard src/*.c)
OBJS = $(SRC:%.c=%.o)

TEST_SUITS = \
	slauum dlauum clauum zlauum \
	spotrf dpotrf cpotrf zpotrf \
	ssygst dsygst chegst zhegst \
	ssytrf dsytrf csytrf chetrf zsytrf zhetrf \
	sgetrf dgetrf cgetrf zgetrf \
	strsyl dtrsyl ctrsyl ztrsyl \
	stgsyl dtgsyl ctgsyl ztgsyl \
	sgemm_tr dgemm_tr cgemm_tr zgemm_tr
TESTS = $(TEST_SUITS:%=test/%.pass)  # dummies
TEST_EXES = $(TEST_SUITS:%=test/%.x)

.SECONDARY: $(TEST_EXES)
.PHONY: test

# ReLAPACK compilation 

$(LIBNAME): $(OBJS)
	$(AR) -r $@ $^

%.o: %.c config.h
	$(CC) $(CFLAGS) -c $< -o $@


# ReLAPACK testing

test: $(TEST_EXES) $(TESTS)
	@echo "passed all tests"

test/%.pass: test/%.x
	@echo -n $*:
	@./$< > /dev/null && echo " pass" || (echo " FAIL" && ./$<)

test/s%.x: test/x%.c test/util.o $(LIBNAME) test/config.h test/test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=s $< test/util.o -o $@ $(LINK_TEST) $(LIBNAME) $(LINK_TEST)

test/d%.x: test/x%.c test/util.o $(LIBNAME) test/config.h test/test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=d $< test/util.o -o $@ $(LINK_TEST) $(LIBNAME) $(LINK_TEST)

test/c%.x: test/x%.c test/util.o $(LIBNAME) test/config.h test/test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=c $< test/util.o -o $@ $(LINK_TEST) $(LIBNAME) $(LINK_TEST)

test/z%.x: test/x%.c test/util.o $(LIBNAME) test/config.h test/test.h
	$(CC) $(CFLAGS) -DDT_PREFIX=z $< test/util.o -o $@ $(LINK_TEST) $(LIBNAME) $(LINK_TEST)


# cleaning up

clean:
	rm -f $(LIBNAME) $(OBJS) test/util.o test/*.x
