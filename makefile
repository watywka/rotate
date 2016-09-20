all: main.o solve.o
	gcc main.o solve.o -o solve -lm -lrt
debug: main.o solve.o
	gcc main.o solve.o -o solve -lm -lrt -g
main.o: main.c
	gcc -c main.c
solve.o: solve.c solve.h
	gcc -c solve.c
clean:
	rm -rf *.o
