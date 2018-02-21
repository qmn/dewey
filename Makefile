CC = clang
CFLAGS += -Wall -g -pedantic -std=c99 -I/usr/local/include -Wmissing-field-initializers
EXEC_CFLAGS += -L/usr/local/lib -lyaml -lpng -pg -lgd
BUILD_DIR = build
SRCS = $(wildcard *.c)
OBJS = $(foreach f,$(SRCS:.c=.o),$(BUILD_DIR)/$f)

default: dewey

$(BUILD_DIR) :
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o : %.c $(BUILD_DIR)
	$(CC) -c $(CFLAGS) -o $@ $<

dewey: $(OBJS)
	$(CC) $(CFLAGS) $(EXEC_CFLAGS) -o $@ $^

test: dewey
	./dewey counter.blif

clean:
	rm -rf build dewey

tags: $(SRCS)
	ctags -R $(wildcard *.[ch])

.PHONY: default clean
