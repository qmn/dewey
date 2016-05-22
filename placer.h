#ifndef __PLACER_H__
#define __PLACER_H__

#include "blif.h"

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

struct placement {
	char *name;
	struct coordinate placement;
	unsigned long turns;
	unsigned long n_pins;
	struct pin **pins;
};

struct cell_placements {
	struct placement **placements;
	unsigned long length;
};

struct cell_placements *simulated_annealing_placement(struct cell_placements *,
		struct dimensions *,
		double,
		unsigned int, unsigned int);
#endif /* __PLACER_H__ */
