all: pointer map
	
pointer: pointer.o
	gcc -o pointer pointer.o

map: map.o
	gcc -o map map.o

pointer.o: pointer.c
	gcc -c pointer.c

map.o: map.c
	gcc -c map.c

binaries=map pointer
clean:
	rm -f $(binaries)
	rm -f *.o
	rm -f *~
