
# Это комментарий, который говорит, что переменная CC указывает компилятор, используемый для сборки
CC=mpicxx
#Это еще один комментарий. Он поясняет, что в переменной CFLAGS лежат флаги, которые передаются компилятору
CFLAGS=-c -Wall
FILENAME = main
	
$(FILENAME): $(FILENAME).o reshotka.o
	$(CC) $(FILENAME).o reshotka.o -o $(FILENAME)

$(FILENAME).o : $(FILENAME).cpp reshotka.hpp
	$(CC) -c $(FILENAME).cpp -o $(FILENAME).o

reshotka.o: reshotka.cpp reshotka.hpp
	$(CC) -c reshotka.cpp -o reshotka.o

clear:
	rm -rf *.o $(FILENAME)
