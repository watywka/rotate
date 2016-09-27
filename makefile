all: main.o solve.o
	gcc main.o solve.o -o solve -lm -lrt 
debug: main.o solve.o
	gcc main.o solve.o -o solve -lm -lrt -g 
main.o: main.c solve.h
	gcc -c main.c -std=gnu99
solve.o: solve.c solve.h
	gcc -c solve.c -std=gnu99
clean:
	rm -rf *.o
