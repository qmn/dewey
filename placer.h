#ifndef __PLACER_H__
#define __PLACER_H__

#include "blif.h"
#include "cell.h"

struct segment {
	struct coordinate start;
	struct coordinate end;
};

struct segments {
	struct segment *segments;
	int n_segments;
};

struct placement {
	struct logic_cell *cell;
	struct coordinate placement;
	unsigned long turns;

	/* a mapping of net names to net_t and should be the same size as cell->n_pins */
	net_t *nets;
};

struct cell_placements {
	struct placement *placements;
	unsigned long n_placements;

	int n_nets;
};

struct cell_placements *simulated_annealing_placement(struct cell_placements *,
		struct dimensions *,
		double,
		unsigned int, unsigned int);

struct cell_placements *placer_initial_place(struct blif *, struct cell_library *);
struct dimensions compute_placement_dimensions(struct cell_placements *);

void print_cell_placements(struct cell_placements *);
#endif /* __PLACER_H__ */
