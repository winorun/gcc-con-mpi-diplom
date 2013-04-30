ifndef ComSpec
	LIB = -L"C:\Program Files\MPICH2\lib" -lmpi
	INCL = -I"C:\Program Files\MPICH2\include"
	TYPE = .exe
	CC=g++
else
	CC=mpicxx
	#INCL = -I"C:\Program Files\MPICH2\include"
endif

FILENAME = main

$(FILENAME): $(FILENAME).o reshotka.o
	$(CC) $(FILENAME).o reshotka.o -o $(FILENAME)$(TYPE) $(LIB)

$(FILENAME).o : $(FILENAME).cpp reshotka.hpp
	$(CC) -c $(FILENAME).cpp -o $(FILENAME).o $(INCL)

reshotka.o: reshotka.cpp reshotka.hpp
	$(CC) -c reshotka.cpp -o reshotka.o $(INCL)

clear:
	rm -rf *.o $(FILENAME).exe

#~ gcc -c hello_mpi.c -o hello_mpi.o -I"C:\Program Files\MPICH2\include" 
#~ gcc -o hello_mpi.exe hello_mpi.o -L"C:\Program Files\MPICH2\lib" -lmpi
