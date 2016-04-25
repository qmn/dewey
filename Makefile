default: dewey

CFLAGS=-Wall -pedantic

%.o : %.c
	gcc -c $(CFLAGS) -o $@ $<

dewey: dewey.o placer.o
	gcc $(CFLAGS) -o $@ $^

clean:
	rm -rf *.o

.PHONY: default clean
