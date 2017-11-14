#ifndef __PLACER_H__
#define __PLACER_H__

#include "blif.h"
#include "cell.h"

#define CONSTR_NONE       0
#define CONSTR_NO_ROTATE  (1L << 0)
#define CONSTR_KEEP_LEFT  (1L << 1)
#define CONSTR_KEEP_RIGHT (1L << 2)

#define CONSTR_MASK_NO_INTERCHANGE (CONSTR_KEEP_LEFT | CONSTR_KEEP_RIGHT)

struct pin_placements {
	int n_pins;
	struct placed_pin *pins;
};

struct placement {
	struct logic_cell *cell;
	struct coordinate placement;
	unsigned long turns;

	/* a mapping of net names to net_t and should be the same size as cell->n_pins */
	net_t *nets;

	unsigned long constraints;
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

struct cell_placements *copy_placements(struct cell_placements *);
void placements_displace(struct cell_placements *, struct coordinate disp);
void placements_reconstrain(struct cell_placements *);
void free_cell_placements(struct cell_placements *);

struct cell_placements *placer_initial_place(struct blif *, struct cell_library *);
struct dimensions compute_placement_dimensions(struct cell_placements *);

void print_cell_placements(struct cell_placements *);

struct pin_placements *placer_place_pins(struct cell_placements *);
void free_pin_placements(struct pin_placements *);

struct net_pin_map *placer_create_net_pin_map(struct pin_placements *);
struct coordinate extend_pin(struct placed_pin *);
void free_net_pin_map(struct net_pin_map *);

struct extraction *extract_placements(struct cell_placements *);
void free_extraction(struct extraction *);
#endif /* __PLACER_H__ */
