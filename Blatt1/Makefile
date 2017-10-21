all: pointer map
	
pointer: pointer.o
	gcc -o pointer pointer.o

map: map.o
	gcc -o map map.o

pointer.o: pointer.c
	gcc -c pointer.c

map.o: map.c
	gcc -c map.c

clean:
	rm -f all
	rm -f *.o
	rm -f *~
