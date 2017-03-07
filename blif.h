#ifndef __BLIF_H__
#define __BLIF_H__

#include <stdio.h>

typedef unsigned int net_t;

struct pin {
	char *name;
	net_t net;
};


struct blif_cell {
	char *name;

	unsigned int n_pins;
	struct pin **pins;
};

struct blif {
	char *model;

	unsigned int n_nets;
	char **net_names;

	unsigned int n_inputs;
	net_t *inputs;

	unsigned int n_outputs;
	net_t *outputs;

	unsigned int n_cells;
	struct blif_cell **cells;
};

struct blif *read_blif(FILE *);

void free_blif(struct blif *);
#endif /* __BLIF_H__ */
