#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <signal.h>

#include "blif.h"
#include "coord.h"
#include "extract.h"
#include "placer.h"
#include "segment.h"

#define MIN_MARGIN 2
#define EDGE_MARGIN MIN_MARGIN

#define MIN_WINDOW_WIDTH 8
#define MIN_WINDOW_HEIGHT 8

#define DISPLACE_INTERCHANGE_RATIO 5.0

enum placement_method {
	NONE,
	DISPLACE,
	REORIENT,
	INTERCHANGE
};

static int max(int a, int b)
{
	return a > b ? a : b;
}

static int min(int a, int b)
{
	return a < b ? a : b;
}

static struct dimensions compute_unconstrained_placement_dimensions(struct cell_placements *cp)
{
	int i;
	struct dimensions d = {0, 0, 0};

	for (i = 0; i < cp->n_placements; i++) {
		struct placement p = cp->placements[i];

		struct coordinate c = p.placement;
		struct dimensions pd = p.cell->dimensions[p.turns];

		int cell_x = c.x + pd.x + 1;
		int cell_y = c.y + pd.y + 1;
		int cell_z  = c.z + pd.z + 1;

		if (!(p.constraints & CONSTR_KEEP_RIGHT))
			d.x = max(cell_x, d.x);

		d.y = max(cell_y, d.y);
		d.z = max(cell_z, d.z);
	}

	return d;
}

void placements_reconstrain(struct cell_placements *cp)
{
	struct dimensions d = compute_unconstrained_placement_dimensions(cp);
	for (int i = 0; i < cp->n_placements; i++) {
		struct placement *p = &(cp->placements[i]);

		if (p->constraints & CONSTR_KEEP_LEFT)
			p->placement.x = 0;
		else if (p->constraints & CONSTR_KEEP_RIGHT)
			p->placement.x = d.x + EDGE_MARGIN;
	}
}

/*
 * Given an old placement, generate a new placement by either switching
 * the location of two cells or displacing a cell or rotating it. This
 * method modifies the original placement, i.e., no copy is made.
 * 
 * T is the current temperature; T_0 is the starting temperature. This
 * is used to scale the window for displacing a cell.
 * 
 * method can be "displace" or "reorient".
 * 
 * displace_interchange_ratio is the ratio of how often you displace
 * a cell and how often you interchange it with another cell.
 *
 * it returns the method used
 */
static enum placement_method generate(struct cell_placements *placements,
		struct dimensions dimensions,
		double t, double t_0,
		enum placement_method method)
{
	unsigned long cell_a_idx, cell_b_idx;
	struct placement *cell_a;
	double p, interchange_threshold;
	double scaling_factor;
	long window_height, window_width;

	/* compute the probabilty we change this cell */
	p = (double)random() / (double)(RAND_MAX);
	interchange_threshold = (1.0 / DISPLACE_INTERCHANGE_RATIO);

	/* select a random cell to interchange, displace, or reorient */
	cell_a_idx = (unsigned long)random() % placements->n_placements;
	cell_a = &(placements->placements[cell_a_idx]);

	// determine the window size
	scaling_factor = log(t) / log(t_0);

	/* figure the most this placement can move */
	window_height = max(lround(dimensions.z * scaling_factor), 20);
	window_width = max(lround(dimensions.x * scaling_factor), 20);

	window_height = max(window_height, MIN_WINDOW_HEIGHT);
	window_width = max(window_width, MIN_WINDOW_WIDTH);

	if (p > interchange_threshold) {
		/* select another cell_a if we can't interchange this one */
		while (cell_a->constraints & CONSTR_MASK_NO_INTERCHANGE) {
			cell_a_idx = (unsigned long)random() % placements->n_placements;
			cell_a = &(placements->placements[cell_a_idx]);
		}
		
		/* select another cell */
		struct placement *cell_b = NULL;
		do {
			cell_b_idx = (unsigned long)random() % placements->n_placements;
			cell_b = &(placements->placements[cell_b_idx]);
		} while (cell_b_idx == cell_a_idx || cell_b->constraints & CONSTR_MASK_NO_INTERCHANGE);

		// attempt an interchange only if a window would fit both cells
		struct coordinate center_a = cell_a->placement, center_b = cell_b->placement;
		if (abs(center_a.x - center_b.x) <= window_width && abs(center_a.z - center_b.z) <= window_height) {
			/* interchange the cells' placements */
			struct coordinate tmp = cell_a->placement;
			cell_a->placement = cell_b->placement;
			cell_b->placement = tmp;

			// printf("[placer] interchange %d (%d, %d, %d) with %d (%d, %d, %d)\n",
			//	cell_a_idx, cell_a->placement.y, cell_a->placement.z, cell_a->placement.x,
			//	cell_b_idx, cell_b->placement.y, cell_b->placement.z, cell_b->placement.x);

			return INTERCHANGE;
		}

		return NONE;
	} else if (method == DISPLACE) {

		/* displace */
		int dz = random() % (window_height * 2) - window_height;
		int dx = random() % (window_width * 2) - window_width;
		// printf("[placer] displace %d by dz = %d, dx = %d\n", cell_a_idx, dz, dx);
		cell_a->placement.z += dz;

		if (cell_a->constraints & CONSTR_KEEP_LEFT) {
			cell_a->placement.x = 0;
		} else if (cell_a->constraints & CONSTR_KEEP_RIGHT) {
			struct dimensions dd = compute_unconstrained_placement_dimensions(placements);
			cell_a->placement.x = dd.x + EDGE_MARGIN; // plus margin
		} else {
			cell_a->placement.x += dx;
		}

		return DISPLACE;
	} else {
		/* reorient */
		// printf("[placer] rotate\n");
		if (!(cell_a->constraints & CONSTR_NO_ROTATE)) {
			cell_a->turns = (cell_a->turns + 1) % 4;
			return REORIENT;
		}
		return NONE;
	}
}


/*
 * Performs a deep copy of cell_placements.
 */
struct cell_placements *copy_placements(struct cell_placements *old_placements)
{
	struct cell_placements *new_placements;
	int i;

	new_placements = malloc(sizeof(struct cell_placements));

	new_placements->n_placements = old_placements->n_placements;
	new_placements->n_nets = old_placements->n_nets;
	new_placements->placements = malloc(old_placements->n_placements * sizeof(struct placement));
	memcpy(new_placements->placements, old_placements->placements, sizeof(struct placement) * old_placements->n_placements);

	/* copy np->nets */
	for (i = 0; i < old_placements->n_placements; i++) {
		struct placement *p = &new_placements->placements[i];
		p->nets = malloc(sizeof(net_t) * p->cell->n_pins);
		memcpy(p->nets, old_placements->placements[i].nets, sizeof(net_t) * p->cell->n_pins);
	}

	return new_placements;
}

/* copies the cell placements and displaces them all by coordinate disp */
void placements_displace(struct cell_placements *cp, struct coordinate disp)
{
	for (int i = 0; i < cp->n_placements; i++) {
		cp->placements[i].placement = coordinate_add(cp->placements[i].placement, disp);
	}
}

void free_cell_placements(struct cell_placements *placements)
{
	free(placements->placements);
	free(placements);
}

/* based on the current placements, how large is the design? */
struct dimensions compute_placement_dimensions(struct cell_placements *cp)
{
	int i;
	struct dimensions d = {0, 0, 0};

	for (i = 0; i < cp->n_placements; i++) {
		struct placement p = cp->placements[i];

		struct coordinate c = p.placement;
		struct dimensions pd = p.cell->dimensions[p.turns];

		int cell_x = c.x + pd.x + 1;
		int cell_y = c.y + pd.y + 1;
		int cell_z  = c.z + pd.z + 1;

		d.x = max(cell_x, d.x);
		d.y = max(cell_y, d.y);
		d.z = max(cell_z, d.z);
	}

	return d;
}

static int overlap(int s1, int e1, int s2, int e2)
{
	assert(e1 >= s1 && e2 >= s2);
	int space = max(e1, e2) - min(s1, s2);
	int taken = (e1 - s1) + (e2 - s2);
	return taken > space ? (taken - space) : 0;
}

static struct dimensions placement_overlaps(struct placement p, struct placement q, int use_margin)
{
	struct coordinate pc = p.placement, qc = q.placement;
	struct dimensions pd = p.cell->dimensions[p.turns], qd = q.cell->dimensions[q.turns];
	int nt = 1; // nt = no touch -- cannot be within this area
	int m = use_margin ? min(p.margin, q.margin) : 0;

	// M    C    C+D   M
	// |<--X[__P__]X-->|
	//             |<--X[__Q__]X-->|
        //             M    C    C+D   M
	int yo = overlap(pc.y - nt - m, pc.y + pd.y + nt + m, qc.y - nt - m, qc.y + qd.y + nt + m);
	int zo = overlap(pc.z - nt - m, pc.z + pd.z + nt + m, qc.z - nt - m, qc.z + qd.z + nt + m);
	int xo = overlap(pc.x - nt - m, pc.x + pd.x + nt + m, qc.x - nt - m, qc.x + qd.x + nt + m);

	// if there is an overlap due to the margin, reduce overlap by the
	// size of the lesser margin (as far down to zero)
	yo = yo > 0 ? max(yo - m, 0) : 0;
	zo = zo > 0 ? max(zo - m, 0) : 0;
	xo = xo > 0 ? max(xo - m, 0) : 0;

	return (struct dimensions){yo, zo, xo};
}

// compute overlap penalty in a smarter way:
// compare placements cell-wise, to avoid large memory allocation
// if cells are more apart than the largest dimension of all of them,
// there is no overlap, and continue
// otherwise, only then do you create a overlap grid
struct overlap_penalty {
	int violations;
	int score;
};

static struct overlap_penalty compute_overlap_penalty_pairwise(struct cell_placements *cp)
{
	int i, j;
	struct placement p, q;
	int violations, score;

	score = 0;
	violations = 0;

	for (i = 0; i < cp->n_placements; i++) {
		p = cp->placements[i];

		for (j = i + 1; j < cp->n_placements; j++) {
			q = cp->placements[j];

			struct dimensions ov = placement_overlaps(p, q, 0);
			violations += ov.x * ov.y * ov.z;

			ov = placement_overlaps(p, q, 1);
			score += ov.x * ov.y * ov.z;
		}
	}

	return (struct overlap_penalty){violations, score};
}

struct coordinate extend_in_direction(enum ordinal_direction facing, struct coordinate c)
{
	switch (facing) {
		case EAST: c.x++; break;
		case WEST: c.x--; break;
		case NORTH: c.z--; break;
		case SOUTH: c.z++; break;
		default: break;
	}

	return c;
}

struct coordinate extend_pin(struct placed_pin *p)
{
	return extend_in_direction(p->cell_pin->facing, p->coordinate);
}

/* given a list of pin placements, generate a list collecting
 * pins based on the net they belong to */
struct net_pin_map *placer_create_net_pin_map(struct pin_placements *pp)
{
	struct net_pin_map *npm = malloc(sizeof(struct net_pin_map));

	/* count number of nets */
	int n_nets = 0;
	for (int i = 0; i < pp->n_pins; i++)
		n_nets = max(n_nets, pp->pins[i].net);

	assert(n_nets > 0);

	npm->n_nets = n_nets;

	/* create size of per-net maps */
	int *n_pins_for_net = calloc(n_nets + 1, sizeof(int));
	for (net_t i = 0; i < pp->n_pins; i++) {
		n_pins_for_net[pp->pins[i].net]++;
	}

	npm->n_pins_for_net = n_pins_for_net;

	/* create per-net pin lists */
	struct placed_pin **pins_for_net = calloc(n_nets + 1, sizeof(struct placed_pin *));
	for (net_t i = 1; i < n_nets + 1; i++) {
		pins_for_net[i] = calloc(n_pins_for_net[i], sizeof(struct placed_pin));
	}

	/* copy those placed pins into the map */
	int *this_pins_for_net = calloc(n_nets + 1, sizeof(int));
	for (int i = 0; i < pp->n_pins; i++) {
		net_t j = pp->pins[i].net;
		pins_for_net[j][this_pins_for_net[j]++] = pp->pins[i];
		assert(this_pins_for_net[j] <= n_pins_for_net[j]);
	}
	free(this_pins_for_net);

	npm->pins = pins_for_net;

	return npm;
}

void free_net_pin_map(struct net_pin_map *npm)
{
	for (int i = 1; i < npm->n_nets + 1; i++)
		free(npm->pins[i]);
	free(npm->pins);
	free(npm->n_pins_for_net);
	free(npm);
}

/* based on a set of cell placements, build the list of pins
 */
struct pin_placements *placer_place_pins(struct cell_placements *cp)
{
	int pin_count = 0;
	for (int i = 0; i < cp->n_placements; i++)
		pin_count += cp->placements[i].cell->n_pins;

	struct placed_pin *pins = calloc(pin_count, sizeof(struct placed_pin));

	int k = 0;
	for (int i = 0; i < cp->n_placements; i++) {
		struct placement pl = cp->placements[i];
		struct logic_cell *c = pl.cell;

		for (int j = 0; j < c->n_pins; j++) {
			struct coordinate b = pl.placement;
			struct coordinate p = c->pins[pl.turns][j].coordinate;

			struct placed_pin pin;
			pin.coordinate = coordinate_add(b, p);
			pin.cell = c;
			pin.cell_pin = &(c->pins[pl.turns][j]);
			pin.net = pl.nets[j];
			pin.parent = NULL;

			pins[k++] = pin;
		}
	}

	assert(k == pin_count);

	struct pin_placements *pp = malloc(sizeof(struct pin_placements));
	pp->n_pins = pin_count;
	pp->pins = pins;

	return pp;
}

void free_pin_placements(struct pin_placements *pp)
{
	free(pp->pins);
	free(pp);
}

/* determine the length of wire needed to connect all points, using
 * the minimal spanning tree that covers the wires. it's not a perfect metric,
 * but it is a good enough estimate
 */
#define PENALTY_MALLOC_START 4
static int compute_wire_length_penalty(struct cell_placements *cp)
{
	int penalty = 0;

	int (*distance_metric)(struct coordinate, struct coordinate) = distance_pythagorean;

	/* map pins to nets */
	struct pin_placements *pp = placer_place_pins(cp);
	struct net_pin_map *npm = placer_create_net_pin_map(pp);
	free_pin_placements(pp);

	/* for each net, compute the constituent pin coordinates */
	for (net_t i = 1; i < npm->n_nets; i++) {
		int n_pins = npm->n_pins_for_net[i];
		assert(n_pins >= 0);

		if (n_pins == 0 || n_pins == 1)
			continue;

		if (n_pins == 2) {
			penalty += distance_metric(npm->pins[i][0].coordinate, npm->pins[i][1].coordinate);
			continue;
		}

		struct coordinate *coords = malloc(n_pins * sizeof(struct coordinate));
		for (int j = 0; j < n_pins; j++)
			coords[j] = extend_pin(&npm->pins[i][j]);

		struct segments *mst = create_mst(coords, n_pins);
		for (int seg = 0; seg < mst->n_segments; seg++) {
			int d = distance_metric(mst->segments[seg].start, mst->segments[seg].end);
			penalty += d;
		}
		free_segments(mst);
		free(coords);
	}
	free_net_pin_map(npm);

	return penalty;
}

static int distance_outside_boundary(struct coordinate c, struct dimensions b)
{
	int dz = 0, dx = 0;

	if (c.z > b.z)
		dz = c.z - b.z;
	else if (c.z < 0)
		dz = -c.z;

	if (c.x > b.x)
		dx = c.x - b.x;
	else if (c.x < 0)
		dx = -c.x;

	/* compute pythagorean theoretic distance from boundary */
	return roundl(sqrt((double)(dx * dx)) + (double)(dz * dz));
}

static int compute_out_of_bounds_penalty(struct cell_placements *placements, struct dimensions boundary)
{
	int penalty = 0;

	/* add penalty based on pythagorean theoretic distance from boundary */
	for (int i = 0; i < placements->n_placements; i++) {
		struct coordinate c = placements->placements[i].placement;
		penalty += distance_outside_boundary(c, boundary);
	}

	return penalty;
}

static int compute_squareness_penalty(struct cell_placements *cp)
{
	struct dimensions d = compute_placement_dimensions(cp);;
	return max(d.x, d.z) * max(d.x, d.z) * d.y;
}

/* computes the area required to implement this design */
static int compute_design_size_penalty(struct cell_placements *placements)
{
	struct dimensions d = compute_placement_dimensions(placements);
	return d.x * d.z;
}

// #define PLACER_SCORE_DEBUG

/* requires placements be re-centered so that all numbers positive */
static int score(struct cell_placements *placements, struct dimensions boundary)
{
	struct overlap_penalty overlap = compute_overlap_penalty_pairwise(placements);
	int wire_length = compute_wire_length_penalty(placements);
	int bounds = compute_out_of_bounds_penalty(placements, boundary);
	int design_size = compute_design_size_penalty(placements);
	int squareness = compute_squareness_penalty(placements);
#ifdef PLACER_SCORE_DEBUG
	printf("[placer] score overlap: %d, wire_length: %d, out_of_bounds: %d, design_size: %d\n", overlap, wire_length, bounds, design_size);
#endif
	return (overlap.violations * overlap.violations * 10) + overlap.score * 2 + wire_length + bounds + design_size + squareness;
}

static int accept(int new_score, int old_score, double t)
{
	double ratio, acceptance_criterion;
	ratio = (double)(new_score - old_score) / t;

	acceptance_criterion = fmin(1.0, exp(-ratio));

	long r = random();
	long y = lround(acceptance_criterion * RAND_MAX);
	return r < y;
}

static double update(double t, double (*alpha)(double))
{
	return t * alpha(t);
}

static double fixed_alpha(double t)
{
	return 0.9;
}

static int interrupt_placement = 0;

static void placer_sigint_handler(int a)
{
	printf("Interrupt\n");
	interrupt_placement = 1;
}

// #define PLACER_GENERATION_DEBUG

/*
 * Performs simulated annealing (the Timberwolf algorithm)
 * to produce a placement.
 */
struct cell_placements *simulated_annealing_placement(struct cell_placements *initial_placements,
		struct dimensions *dimensions,
		double t_0,
		unsigned int iterations, unsigned int generations)
{
	double t;
	struct cell_placements *best_placements, *new_placements;
	unsigned int i, g;
	enum placement_method method, method_used;
	int new_score, old_score, taken_score;
	int violating_overlaps;

	taken_score = 0;

	// generations = 10 * initial_placements->n_placements;

	struct dimensions wanted;
	struct dimensions d = compute_placement_dimensions(initial_placements);

	wanted.x = 100;
	wanted.z = 100;
	wanted.y = 5;

	printf("[placer] beginning simulated annealing placement\n");

	t = t_0;
	best_placements = initial_placements;
	old_score = score(initial_placements, wanted);

	int match_iterations = 0;
	int match_score = old_score;
	int stop_iterations = 100;
	violating_overlaps = 0;

	interrupt_placement = 0;
	signal(SIGINT, placer_sigint_handler);

	i = 0;
	do {
#ifdef PLACER_GENERATION_DEBUG
		printf("[placer] iteration = %d\n", i);
#endif
		method = DISPLACE;

		for (g = 0; g < generations; ) {
#ifdef PLACER_GENERATION_DEBUG
			printf("[placer] generation = %d\n", g);
#endif
			/* make a copy of these placements */
			new_placements = copy_placements(best_placements);
			// print_cell_placements(new_placements);
#ifdef PLACER_GENERATION_DEBUG
			printf("[placer] made a copy\n");
#endif
			method_used = generate(new_placements, dimensions_piecewise_max(wanted, d), t, t_0, method);
			if (method_used != NONE)
				g++;
			recenter(new_placements, NULL, 0);
#ifdef PLACER_GENERATION_DEBUG
			printf("[placer] generated a new\n");
#endif
			new_score = score(new_placements, wanted);

#ifdef PLACER_GENERATION_DEBUG
			printf("[placer] old_score = %d, new_score = %d\n",
				old_score, new_score);
#endif

			if (accept(new_score, old_score, t)) {
#ifdef PLACER_GENERATION_DEBUG
				printf("[placer] placer accepts\n");
#endif
				/* 
				 * accept this new placement, free the old
				 * placements, and replace it with the new ones
				 */
				free_cell_placements(best_placements);
				best_placements = new_placements;
				taken_score = new_score;
				if (method_used == REORIENT)
					method = DISPLACE;
			} else {
#ifdef PLACER_GENERATION_DEBUG
				printf("[placer] placer rejects\n");
#endif
				/* reject the new placement */
				free_cell_placements(new_placements);
				taken_score = old_score;
				if (method_used == DISPLACE)
					method = REORIENT;
			}

		}
		// printf("[placer] T = %4.2f\n", t);

		if (taken_score == match_score) {
			match_iterations++;
		} else {
			match_score = taken_score;
			match_iterations = 0;
		}

		old_score = taken_score;

		d = compute_placement_dimensions(best_placements);
		violating_overlaps = compute_overlap_penalty_pairwise(best_placements).violations;

		printf("\rIteration: %4d, Score: %6u (violations: %6u, design size: %d x %d), Temperature: %6.0f", (i + 1), taken_score, violating_overlaps, d.z, d.x, t);
		fflush(stdout);
		// print_cell_placements(best_placements);

		t = update(t, fixed_alpha);
	} while ((i++ < iterations || match_iterations < stop_iterations || violating_overlaps > 0) && !interrupt_placement);

	signal(SIGINT, SIG_DFL);
	printf("\nPlacement complete\n");

	return best_placements;
}

/* tries to match a cell in the blif to a cell in the cell library,
 * returning NULL if this map fails */
static struct logic_cell *map_cell_to_library(struct blif_cell *blif_cell, struct cell_library *cl)
{
	int i;
	struct logic_cell *c;

	for (i = 0; i < cl->n_cells; i++) {
		c = &(cl->cells[i]);
		if (strcmp(c->name, blif_cell->name) == 0)
			return c;
	}

	return NULL;
}

static struct logic_cell *get_input_pin(struct cell_library *cl)
{
	struct blif_cell input_cell = {"input_pin", 0, NULL};
	struct logic_cell *lc = map_cell_to_library(&input_cell, cl);
	assert(lc);
	return lc;
}

static struct logic_cell *get_output_pin(struct cell_library *cl)
{
	struct blif_cell output_cell = {"output_pin", 0, NULL};
	struct logic_cell *lc = map_cell_to_library(&output_cell, cl);
	assert(lc);
	return lc;
}

/* for each cell in the blif, map it to a cell in the cell library,
 * and allocate struct cell_placements, but don't make any placements */
static struct cell_placements *map_blif_to_cell_library(struct blif *blif, struct cell_library *cl)
{
	struct cell_placements *placements;
	int i;

	placements = malloc(sizeof(struct cell_placements));
	placements->n_placements = blif->n_inputs + blif->n_cells + blif->n_outputs;
	placements->n_nets = blif->n_nets;

	placements->placements = malloc(sizeof(struct placement) * placements->n_placements);
	int cell_count = 0;

	struct logic_cell *input_pin = get_input_pin(cl);
	struct logic_cell *output_pin = get_output_pin(cl);
	/* place inputs */
	for (i = 0; i < blif->n_inputs; i++) {
		struct coordinate c = {0, 0, 0};
		net_t *nets = malloc(sizeof(net_t));
		nets[0] = blif->inputs[i];
		struct placement p = {input_pin, c, 0, nets, CONSTR_NO_ROTATE | CONSTR_KEEP_LEFT, EDGE_MARGIN};
		placements->placements[cell_count++] = p;
	}

	for (i = 0; i < blif->n_outputs; i++) {
		struct coordinate c = {0, 0, 0};
		net_t *nets = malloc(sizeof(net_t));
		nets[0] = blif->outputs[i];
		struct placement p = {output_pin, c, 0, nets, CONSTR_NO_ROTATE | CONSTR_KEEP_RIGHT, EDGE_MARGIN};
		placements->placements[cell_count++] = p;
	}

	/* place cells */
	for (i = 0; i < blif->n_cells; i++) {
		// printf("[placer] cell %d\n", i);
		struct blif_cell *c = blif->cells[i];
		struct placement p;
		p.nets = malloc(sizeof(net_t) * c->n_pins);


		p.cell = map_cell_to_library(c, cl);
		if (!p.cell)
			printf("[placer] could not map blif cell (%s) to cell library\n", c->name);
		else
			printf("[placer] map cell %d to %s\n", i, p.cell->name);

		/* copy blif pins net to p->nets */
		for (int j = 0; j < c->n_pins; j++) {
			// printf("[placer] placements[%d](%s)->nets[%d] = c->pins[%d]->net = %u\n", i, c->name, j, j, c->pins[j]->net);
			p.nets[j] = c->pins[j]->net;
		}

		p.placement.x = 0;
		p.placement.y = 0;
		p.placement.z = 0;

		p.turns = 0;

		p.constraints = CONSTR_NONE;

		p.margin = EDGE_MARGIN;

		placements->placements[cell_count++] = p;
	}

	return placements;
}

void print_cell_placements(struct cell_placements *cp)
{
	int i;
	for (i = 0; i < cp->n_placements; i++) {
		struct placement p = cp->placements[i];
		printf("[placer] placement: %s @ (y=%d, z=%d, x=%d), %lu turns, constraints: 0x%lx\n",
			p.cell->name, p.placement.y, p.placement.z, p.placement.x, p.turns, p.constraints);
	}
}

/* create an initial placement */
struct cell_placements *placer_initial_place(struct blif *blif, struct cell_library *cl)
{
	struct cell_placements *cp = map_blif_to_cell_library(blif, cl);

	int x_margin = 5;
	int z_margin = 5;
	int w = roundl(sqrt(cp->n_placements));

	for (int i = 0, x = 0; i < cp->n_placements; i++, x += x_margin) {
		cp->placements[i].placement.x = x % (x_margin * w);
		cp->placements[i].placement.z = x / (x_margin * w) * z_margin;
	}

	return cp;
}

