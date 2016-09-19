all: main.o solve.o
	g++ main.o solve.o -o solve -lm
main.o: main.cpp
	g++ -c main.cpp
solve.o: solve.cpp
	g++ -c solve.cpp 
clean:
	rm -rf *.o
