default: dewey

CC=gcc-5
CFLAGS=-Wall -pedantic -std=c99 -g

%.o : %.c
	$(CC) -c $(CFLAGS) -o $@ $<

dewey: dewey.o placer.o blif.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -rf *.o

.PHONY: default clean
