default: dewey

CC = gcc-6
CFLAGS += -Wall -pedantic -std=c99 -g -I/usr/local/include -L/usr/local/lib -lyaml

%.o : %.c
	$(CC) -c $(CFLAGS) -o $@ $<

dewey: dewey.o placer.o blif.o cell_library.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -rf *.o

.PHONY: default clean
