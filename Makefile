default: dewey

CC = gcc-6
CFLAGS += -Wall -pedantic -std=c99 -g -I/usr/local/include -L/usr/local/lib -lyaml -lpng -pg

%.o : %.c
	$(CC) -c $(CFLAGS) -o $@ $<

dewey: dewey.o placer.o blif.o cell_library.o vis_png.o
	$(CC) $(CFLAGS) -o $@ $^

test: dewey
	./dewey counter.blif

clean:
	rm -rf *.o

.PHONY: default clean
