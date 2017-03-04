#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "placer.h"

#define MIN_WINDOW_WIDTH 2
#define MIN_WINDOW_HEIGHT 2

#define DISPLACE_INTERCHANGE_RATIO 5.0

enum placement_method {
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
		struct dimensions *dimensions,
		double t, double t_0,
		enum placement_method method)
{
	unsigned long cell_a_idx, cell_b_idx;
	struct placement *cell_a;
	double p, interchange_threshold;
	double scaling_factor;
	long window_height, window_width;

	/* select a random cell to interchange, displace, or reorient */
	cell_a_idx = (unsigned long)random() % placements->n_placements;
	cell_a = placements->placements[cell_a_idx];

	/* compute the probabilty we change this cell */
	p = (double)random() / (double)(RAND_MAX);
	interchange_threshold = (1.0 / DISPLACE_INTERCHANGE_RATIO);

	if (p > interchange_threshold) {
		/* select another cell */
		do {
			cell_b_idx = (unsigned long)random() % placements->n_placements;
		} while (cell_b_idx == cell_a_idx);

		/* interchange the cells' placements */
		struct coordinate tmp = cell_a->placement;
		cell_a->placement = placements->placements[cell_b_idx]->placement;
		placements->placements[cell_b_idx]->placement = tmp;

		// printf("[placer] interchange %p with %p\n", cell_a, placements[cell_b_idx]);

		return INTERCHANGE;
	} else if (method == DISPLACE) {
		scaling_factor = log(t) / log(t_0);

		/* figure the most this placement can move */
		window_height = lround(dimensions->z * scaling_factor);
		window_width = lround(dimensions->x * scaling_factor);

		window_height = max(window_height, MIN_WINDOW_HEIGHT);
		window_width = max(window_width, MIN_WINDOW_WIDTH);

		/* displace */
		int dz = random() % (window_height * 2) - window_height;
		int dx = random() % (window_width * 2) - window_width;
		// printf("[placer] displace cell#%d dz = %d, dx = %d\n", cell_a_idx, dz, dx);
		cell_a->placement.z += dz;
		cell_a->placement.x += dx;


		return DISPLACE;
	} else {
		/* reorient */
		// printf("[placer] rotate\n");
		cell_a->turns = (cell_a->turns + 1) % 4;
		return REORIENT;
	}
}

/*
 * Performs a deep copy of cell_placements.
 */
static struct cell_placements *copy_placements(struct cell_placements *old_placements)
{
	struct cell_placements *new_placements;
	struct placement **p;
	int i;

	new_placements = malloc(sizeof(struct cell_placements));

	new_placements->n_placements = old_placements->n_placements;
	new_placements->n_nets = old_placements->n_nets;
	p = malloc(old_placements->n_placements * sizeof(struct placement *));

	/* deep copy each placement */
	for (i = 0; i < old_placements->n_placements; i++) {
		struct placement *op = old_placements->placements[i];
		struct placement *np = malloc(sizeof(struct placement));
		np->cell = op->cell;
		np->placement = op->placement;
		np->turns = op->turns;
		np->nets = malloc(sizeof(net_t) * np->cell->n_pins);
		for (int j = 0; j < np->cell->n_pins; j++)
			np->nets[j] = op->nets[j];

		p[i] = np;
	}

	new_placements->placements = p;

	return new_placements;
}

void free_placements(struct cell_placements *placements)
{
	int i;

	for (i = 0; i < placements->n_placements; i++) {
		free(placements->placements[i]->nets);
		free(placements->placements[i]);
	}

	free(placements);
}

/* based on the current placements, how large is the design? */
struct dimensions compute_placement_dimensions(struct cell_placements *cp)
{
	int i;
	struct dimensions d = {0, 0, 0};

	for (i = 0; i < cp->n_placements; i++) {
		struct placement *p = cp->placements[i];
		struct coordinate c = p->placement;
		struct dimensions pd = p->cell->dimensions;

		int cell_x = c.x + pd.x;
		int cell_y = c.y + pd.y;
		int cell_z  = c.z + pd.z;

		d.x = max(cell_x, d.x);
		d.y = max(cell_y, d.y);
		d.z = max(cell_z, d.z);
	}

	return d;
}

/* move the entire design so that all coordinates are non-negative */
void recenter(struct cell_placements *cp)
{
	int i;
	struct coordinate d;

	for (i = 0; i < cp->n_placements; i++) {
		struct placement *p = cp->placements[i];
		struct coordinate c = p->placement;

		d.x = (i == 0) ? c.x : min(c.x, d.x);
		d.y = (i == 0) ? c.y : min(c.y, d.y);
		d.z = (i == 0) ? c.z : min(c.z, d.z);
	}

	for (i = 0; i < cp->n_placements; i++) {
		struct placement *p = cp->placements[i];
		p->placement.x -= d.x;
		p->placement.y -= d.y;
		p->placement.z -= d.z;
	}
}

static int distance_cityblock(struct coordinate a, struct coordinate b)
{
	return abs(a.x - b.x) + abs(a.z - b.z);
}

/*
 * Given a list of coordinates, construct the minimum spanning tree
 * with Kruskal's algorithm.
 */
static struct segments *create_mst(struct coordinate *locs, int n_locs)
{
	struct segments *mst = malloc(sizeof(struct segments));
	mst->n_segments = n_locs - 1;
	mst->segments = malloc(sizeof(struct segment *) * (mst->n_segments));
	for (int i = 0; i < mst->n_segments; i++)
		mst->segments[i] = malloc(sizeof(struct segment));

	/* create weight matrix */
	int *weights = malloc(sizeof(int) * n_locs * n_locs);
	for (int y = 0; y < n_locs; y++) {
		for (int x = 0; x < n_locs; x++) {
			weights[y * n_locs + x] = distance_cityblock(locs[x], locs[y]);
/*
			printf("[mst] weight (%d, %d) = %d\n",
				y, x, weights[y * n_locs + x]);
*/
		}
	}

	/* assign each location to its own set */
	int *set_idx = malloc(sizeof(int) * n_locs);
	for (int i = 0; i < n_locs; i++)
		set_idx[i] = i;

	/* repeat N-1 times */
	for (int seg = 0; seg < mst->n_segments; seg++) {
		/* find the minimum of the set that do not have the same set index; x > y */
		int seg_x, seg_y;
		int minimum = -1;

		for (int y = 0; y < n_locs; y++) {
			for (int x = y + 1; x < n_locs; x++) {
				int weight = weights[y * n_locs + x];
				if ((weight < minimum || minimum == -1) && set_idx[x] != set_idx[y]) {
					minimum = weight;
					seg_x = x;
					seg_y = y;
				}
			}
		}

		/* combine the two sets by making set_idx[x] == set_idx[y] */
/*
		printf("[mst] combined elt %d (set %d) with elt %d (set %d)\n",
			seg_x, set_idx[seg_x], seg_y, set_idx[seg_y]);
*/
		set_idx[seg_x] = set_idx[seg_y];

		/* add the segment */
		mst->segments[seg]->start = locs[seg_x];
		mst->segments[seg]->end = locs[seg_y];
	}

	free(weights);
	return mst;
}

static void free_struct_segments(struct segments *s)
{
	for (int i = 0; i < s->n_segments; i++)
		free(s->segments[i]);
	free(s->segments);
	free(s);
}

static char *overlap_tmp = NULL;
static int overlap_tmp_size = 0;

/* determine the number of overlaps in the grid */
static int compute_overlap_penalty(struct cell_placements *cp)
{
	int i, penalty;
	int size;
	struct placement *p;

	struct dimensions d = compute_placement_dimensions(cp);
#ifdef PLACER_SCORE_DEBUG
	printf("[compute_overlap_penalty] l=%d, w=%d, h=%d\n", d.x, d.z, d.y);
#endif

	penalty = 0;
	size = d.x * d.y * d.z;
	if (!overlap_tmp || size > overlap_tmp_size) {
		int new_size = size * 2;
		printf("[compute_overlap_penalty] resizing from %d to %d\n", overlap_tmp_size, new_size);
		free(overlap_tmp);
		overlap_tmp = calloc(new_size, sizeof(char));
		overlap_tmp_size = new_size;
	}
	for (i = 0; i < size; i++)
		overlap_tmp[i] = 0;

	for (i = 0; i < cp->n_placements; i++) {
		p = cp->placements[i];

		struct coordinate c = p->placement;
		struct dimensions pd = p->cell->dimensions;

		int cell_x = c.x + pd.x;
		int cell_y = c.y + pd.y;
		int cell_z  = c.z + pd.z;

		/* if another cell is there, increase the penalty */
		for (int y = c.y; y < cell_y; y++) {
			int by = y * d.z * d.x;
			for (int z = c.z; z < cell_z; z++) {
				int bz = z * d.x;
				for (int x = c.x; x < cell_x; x++) {
					if (overlap_tmp[by + bz + x]++)
						penalty += 1;
				}
			}
		}
	}

	return penalty;
}

/* determine the length of wire needed to connect all points, using
 * the minimal spanning tree that covers the wires. it's not a perfect metric,
 * but it is a good enough estimate
 */
#define PENALTY_MALLOC_START 4
static int compute_wire_length_penalty(struct cell_placements *cp)
{
	net_t i;
	int j;
	int penalty = 0;
	int n_coords = PENALTY_MALLOC_START;
	int found = 0;
	struct coordinate *coords = malloc(sizeof(struct coordinate) * n_coords);

	/* for each net, compute the constituent pin coordinates */
	for (i = 1; i < cp->n_nets; i++) {
		/* clear coords */
		for (j = 0; j < n_coords; j++) {
			coords[j].x = 0;
			coords[j].y = 0;
			coords[j].z = 0;
		}

		/* find placement pins on this net */
		for (j = 0; j < cp->n_placements; j++) {
			struct placement *pl = cp->placements[j];
			struct logic_cell *c = pl->cell;
			for (int k = 0; k < c->n_pins; k++) {
				if (pl->nets[k] == i) {
/*
					printf("[wire_length_penalty] placement %d, cell %s, pin %s\n",
						i, c->name, c->pins[k]->name);
*/
					/* offset pin and add to coords */
					struct coordinate b = pl->placement;
					struct coordinate p = c->pins[k]->coordinate;
					struct coordinate actual;

					actual.x = b.x + p.x;
					actual.y = b.y + p.y;
					actual.z = b.z + p.z;
					coords[found++] = actual;

					/* extend coords if too large */
					if (found >= n_coords) {
						n_coords *= 2;
						coords = realloc(coords, sizeof(struct coordinate) * n_coords);
					}
				}
			}
		}

		struct segments *mst = create_mst(coords, found);
		for (int seg = 0; seg < mst->n_segments; seg++) {
			int d = distance_cityblock(mst->segments[seg]->start, mst->segments[seg]->end);
			penalty += d;
		}
		free_struct_segments(mst);
	}

	free(coords);

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
		struct coordinate c = placements->placements[i]->placement;
		penalty += distance_outside_boundary(c, boundary);
	}

	return penalty;
}

/* computes the area required to implement this design */
static int compute_design_size_penalty(struct cell_placements *placements)
{
	struct dimensions d = compute_placement_dimensions(placements);
	return d.x * d.y * d.z;
}

// #define PLACER_SCORE_DEBUG

/* requires placements be re-centered so that all numbers positive */
static int score(struct cell_placements *placements, struct dimensions *dimensions, struct dimensions boundary)
{
	int overlap = compute_overlap_penalty(placements);
	int wire_length = compute_wire_length_penalty(placements);
	int bounds = compute_out_of_bounds_penalty(placements, boundary);
	int design_size = compute_design_size_penalty(placements);
#ifdef PLACER_SCORE_DEBUG
	printf("[placer] score overlap: %d, wire_length: %d, out_of_bounds: %d, design_size: %d\n", overlap, wire_length, bounds, design_size);
#endif
	return (overlap * overlap) + wire_length + bounds + design_size;
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

	struct dimensions wanted;
	wanted.x = 100;
	wanted.z = 100;
	wanted.y = 5;

	printf("[placer] beginning simulated annealing placement\n");

	t = t_0;
	best_placements = initial_placements;
	old_score = score(initial_placements, dimensions, wanted);

	for (i = 0; i < iterations; i++) {
		printf("[placer] iteration = %d\n", i);
		method = DISPLACE;

		for (g = 0; g < generations; g++) {
#ifdef PLACER_GENERATION_DEBUG
			printf("[placer] generation = %d\n", g);
#endif
			/* make a copy of these placements */
			new_placements = copy_placements(best_placements);
#ifdef PLACER_GENERATION_DEBUG
			printf("[placer] made a copy\n");
#endif
			method_used = generate(new_placements, &wanted, t, t_0, method);
			recenter(new_placements);
#ifdef PLACER_GENERATION_DEBUG
			printf("[placer] generated a new\n");
#endif
			new_score = score(new_placements, dimensions, wanted);

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
				free_placements(best_placements);
				best_placements = new_placements;
				taken_score = new_score;
				if (method_used == REORIENT)
					method = DISPLACE;
			} else {
#ifdef PLACER_GENERATION_DEBUG
				printf("[placer] placer rejects\n");
#endif
				/* reject the new placement */
				free_placements(new_placements);
				taken_score = old_score;
				if (method_used == DISPLACE)
					method = REORIENT;
			}
			old_score = taken_score;
		}
		t = update(t, fixed_alpha);
		printf("[placer] T = %4.2f\n", t);

		printf("Iteration: %d, Score: %d\n", i, taken_score);
		// print_cell_placements(best_placements);
	}

	/* free overlap_tmp structure */
	free(overlap_tmp);

	return best_placements;
}

/* tries to match a cell in the blif to a cell in the cell library,
 * returning NULL if this map fails */
static struct logic_cell *map_cell_to_library(struct blif_cell *blif_cell, struct cell_library *cl)
{
	int i;
	struct logic_cell *c;

	for (i = 0; i < cl->n_cells; i++) {
		c = cl->cells[i];
		if (strcmp(c->name, blif_cell->name) == 0)
			return c;
	}

	return NULL;
}

/* for each cell in the blif, map it to a cell in the cell library,
 * and allocate struct cell_placements, but don't make any placements */
static struct cell_placements *map_blif_to_cell_library(struct blif *blif, struct cell_library *cl)
{
	struct cell_placements *placements;
	int i;

	placements = malloc(sizeof(struct cell_placements));
	placements->n_placements = blif->n_cells;
	placements->n_nets = blif->n_nets;

	placements->placements = malloc(sizeof(struct placement) * placements->n_placements);

	/* for now, only place cells from the blif -- no inputs or outputs */
	for (i = 0; i < blif->n_cells; i++) {
		struct blif_cell *c = blif->cells[i];
		struct placement *p = malloc(sizeof(struct placement));
		p->nets = malloc(sizeof(net_t) * c->n_pins);

		p->cell = map_cell_to_library(c, cl);
		if (!p->cell)
			printf("[placer] could not map blif cell (%s) to cell library\n", c->name);

		/* copy blif pins net to p->nets */
		for (int j = 0; j < c->n_pins; j++)
			p->nets[j] = c->pins[j]->net;

		p->placement.x = 0;
		p->placement.y = 0;
		p->placement.z = 0;

		p->turns = 0;

		placements->placements[i] = p;
	}

	return placements;
}

void print_cell_placements(struct cell_placements *cp)
{
	int i;
	for (i = 0; i < cp->n_placements; i++) {
		struct placement *p = cp->placements[i];
		printf("[placer] placement: %s @ (y=%d, z=%d, x=%d), %lu turns\n",
			p->cell->name, p->placement.y, p->placement.z, p->placement.x, p->turns);
	}
}

/* create an initial placement */
struct cell_placements *placer_initial_place(struct blif *blif, struct cell_library *cl)
{
	struct cell_placements *cp = map_blif_to_cell_library(blif, cl);

	return cp;
}
