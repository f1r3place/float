binaries=float
CC=gcc

all: float

float: float.c
	$(CC) -o float float.c

clean:
	rm -f $(binaries) *.o
