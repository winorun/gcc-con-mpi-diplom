ifneq (,$(findstring WINDOWS,$(PATH)))
	LIB = -L"C:\Program Files\MPICH2\lib;." -lmpi -lreshotka
	INCL = -I"C:\Program Files\MPICH2\include"
	TYPE = .exe
	TYPELIB = lib
	CC=g++
else
	TYPELIB = a
	LIB = -L. -lreshotka
	CC=mpicxx
	#INCL = -I"C:\Program Files\MPICH2\include"
endif

FILENAME = main

$(FILENAME): libreshotka.$(TYPELIB) $(FILENAME).o
	$(CC) $(FILENAME).o  -o $(FILENAME)$(TYPE) $(LIB)

$(FILENAME).o : $(FILENAME).cpp reshotka.hpp
	$(CC) -c $(FILENAME).cpp -o $(FILENAME).o $(INCL)

libreshotka.a: reshotka.o
	ar rc libreshotka.$(TYPELIB) reshotka.o
	ranlib libreshotka.$(TYPELIB)

reshotka.o: reshotka.cpp reshotka.hpp
	$(CC) -c reshotka.cpp -o reshotka.o $(INCL)

clear:
	rm -rf *.o $(FILENAME).exe $(FILENAME)

#~ gcc -c hello_mpi.c -o hello_mpi.o -I"C:\Program Files\MPICH2\include" 
#~ gcc -o hello_mpi.exe hello_mpi.o -L"C:\Program Files\MPICH2\lib" -lmpi
