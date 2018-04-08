# the path to the PNG file containing Minecraft textures
# on a Mac, it will look like this
TEXTURES_FILE ?= "\"/Users/qmn/Library/Application Support/minecraft/textures_0.png\""

# compiler options
CC = gcc
CFLAGS += -Wall -g -pedantic -std=c99 -Wmissing-field-initializers -O3 -DTEXTURES_FILE=$(TEXTURES_FILE)
EXEC_CFLAGS += -lyaml -lpng -pg -lgd

BUILD_DIR = build

SRCS = $(wildcard *.c)
OBJS = $(foreach f,$(SRCS:.c=.o),$(BUILD_DIR)/$f)

# targets
default: dewey

$(BUILD_DIR) :
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o : %.c $(BUILD_DIR)
	$(CC) -c $(CFLAGS) -o $@ $<

dewey: $(OBJS)
	$(CC) $(CFLAGS) $(EXEC_CFLAGS) -o $@ $^

clean:
	rm -rf build dewey

tags: $(SRCS)
	ctags -R $(wildcard *.[ch])

.PHONY: default clean
