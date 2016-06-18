#ifndef __PLACER_H__
#define __PLACER_H__

#include "blif.h"
#include "cell.h"

struct placement {
	struct logic_cell *cell;
	struct coordinate placement;
	unsigned long turns;
};

struct cell_placements {
	struct placement **placements;
	unsigned long n_placements;
};

struct cell_placements *simulated_annealing_placement(struct cell_placements *,
		struct dimensions *,
		double,
		unsigned int, unsigned int);

struct cell_placements *placer_initial_place(struct blif *, struct cell_library *);

void print_cell_placements(struct cell_placements *);
#endif /* __PLACER_H__ */
