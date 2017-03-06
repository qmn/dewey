#ifndef __CELL_H__
#define __CELL_H__

typedef unsigned char block_t;
typedef unsigned char data_t;

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
	char *name;

	enum pin_direction direction;
	enum ordinal_direction facing;
	struct coordinate coordinate;
	int level;
	int clock; /* for DFFs */
};

typedef struct logic_cell_pin * struct_pin_p;
typedef block_t * block_t_p;
typedef data_t * data_t_p;

struct logic_cell {
	char *name;

	unsigned int n_pins;
	struct_pin_p pins[4];

	struct dimensions dimensions[4];
	block_t_p blocks[4];
	data_t_p data[4];

	int delay_combinational;
};

struct cell_library {
	char *name;

	unsigned int n_cells;
	struct logic_cell *cells;
};

struct cell_library *read_cell_library(FILE *);

void free_cell_library(struct cell_library *);
#endif /* __CELL_H__ */
