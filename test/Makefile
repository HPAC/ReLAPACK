
test: util.o $(TESTS)
	@ cd lauum; make
	@ cd trtri; make
	@ cd potrf; make
	@ cd getrf; make
	@ cd sygst; make
	@ cd trsyl; make

clean:
	rm *.o
	cd lauum; make clean
	cd trtri; make clean
	cd potrf; make clean
	cd getrf; make clean
	cd sygst; make clean
	cd trsyl; make clean
