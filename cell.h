#ifndef __CELL_H__
#define __CELL_H__

typedef char block_t;
typedef char data_t;

struct coordinate {
	int y;
	int z;
	int x;
};

struct dimensions {
	unsigned int y;
	unsigned int z;
	unsigned int x;
};

enum pin_direction { INPUT, OUTPUT };
enum ordinal_direction { NORTH, EAST, SOUTH, WEST };

struct logic_cell_pin {
	enum pin_direction direction;
	enum ordinal_direction facing;
	struct coordinate coordinate;
	int level;
	int clock; /* for DFFs */
};

struct logic_cell {
	char *name;

	unsigned int n_pins;
	struct logic_cell_pin **pins;

	struct dimensions dimensions;
	block_t ***blocks;
	data_t ***data;

	int delay_combinational;
};

struct cell_library {
	char *name;

	unsigned int n_cells;
	struct logic_cell **cells;
};

struct cell_library *read_cell_library(FILE *);

void free_cell_library(struct cell_library *);
#endif /* __CELL_H__ */
