CC = gcc-6
CFLAGS += -Wall -pedantic -std=c99 -g -I/usr/local/include -L/usr/local/lib -lyaml -lpng -pg -lgd -O2
BUILD_DIR = build
SRCS = $(wildcard *.c)
OBJS = $(foreach f,$(SRCS:.c=.o),$(BUILD_DIR)/$f)

default: dewey

$(BUILD_DIR) :
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o : %.c $(BUILD_DIR)
	$(CC) -c $(CFLAGS) -o $@ $<

dewey: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

test: dewey
	./dewey counter.blif

clean:
	rm -rf build dewey

.PHONY: default clean
