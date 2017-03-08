#ifndef __PLACER_H__
#define __PLACER_H__

#include "blif.h"
#include "cell.h"

struct pin_placements {
	int n_pins;
	struct placed_pin *pins;
};

struct extraction {
	struct dimensions dimensions;
	block_t *blocks;
	data_t *data;
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

struct net_pin_map {
	int n_nets;
	int *n_pins_for_net;
	struct placed_pin **pins;
};


struct cell_placements *simulated_annealing_placement(struct cell_placements *,
		struct dimensions *,
		double,
		unsigned int, unsigned int);

struct cell_placements *placer_initial_place(struct blif *, struct cell_library *);
struct dimensions compute_placement_dimensions(struct cell_placements *);

void print_cell_placements(struct cell_placements *);

struct pin_placements *placer_place_pins(struct cell_placements *);
void free_pin_placements(struct pin_placements *);

struct net_pin_map *placer_create_net_pin_map(struct pin_placements *);
void free_net_pin_map(struct net_pin_map *);

struct extraction *extract_placements(struct cell_placements *);
void free_extraction(struct extraction *);
#endif /* __PLACER_H__ */
