
test: util.o $(TESTS)
	@ cd trtri; make
	@ cd potrf; make
	@ cd getrf; make

clean:
	rm *.o
	cd trtri; make clean
	cd potrf; make clean
	cd getrf; make clean
