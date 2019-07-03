## M A K E F I L E ##

## FLAGS
OP = -O3
LIB = `pkg-config --cflags --libs opencv`
FLG = -ansi -Wall $(OP) $(LIB)
CC = g++ $(FLG)
CO = g++ -c $(FLG)

## TARGETS

all: fourier.cc complex.o
	$(CC) fourier.cc complex.o -o fvis

complex.o: complex.h complex.c
	g++ complex.c -c $(OP) -o complex.o

clean:
	rm -f *.o
	rm -f ftest
