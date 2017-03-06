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
	cell_a = &(placements->placements[cell_a_idx]);

	/* compute the probabilty we change this cell */
	p = (double)random() / (double)(RAND_MAX);
	interchange_threshold = (1.0 / DISPLACE_INTERCHANGE_RATIO);

	if (p > interchange_threshold) {
		/* select another cell */
		do {
			cell_b_idx = (unsigned long)random() % placements->n_placements;
		} while (cell_b_idx == cell_a_idx);

		struct placement *cell_b = &(placements->placements[cell_b_idx]);

		/* interchange the cells' placements */
		struct coordinate tmp = cell_a->placement;
		cell_a->placement = cell_b->placement;
		cell_b->placement = tmp;

/*
		printf("[placer] interchange %d (%d, %d, %d) with %d (%d, %d, %d)\n",
			cell_a_idx, cell_a->placement.y, cell_a->placement.z, cell_a->placement.x,
			cell_b_idx, cell_b->placement.y, cell_b->placement.z, cell_b->placement.x);
*/

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
		// printf("[placer] displace %d by dz = %d, dx = %d\n", cell_a_idx, dz, dx);
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

void free_placements(struct cell_placements *placements)
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
		struct dimensions pd = p.cell->dimensions;

		int cell_x = c.x + pd.x + 1;
		int cell_y = c.y + pd.y + 1;
		int cell_z  = c.z + pd.z + 1;

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
		struct placement p = cp->placements[i];
		struct coordinate c = p.placement;

		d.x = (i == 0) ? c.x : min(c.x, d.x);
		d.y = (i == 0) ? c.y : min(c.y, d.y);
		d.z = (i == 0) ? c.z : min(c.z, d.z);
	}

	for (i = 0; i < cp->n_placements; i++) {
		struct placement *p = &cp->placements[i];
		p->placement.x -= d.x;
		p->placement.y -= d.y;
		p->placement.z -= d.z;
	}
}

static int distance_cityblock(struct coordinate a, struct coordinate b)
{
	return abs(a.x - b.x) + abs(a.z - b.z);
}

struct mst_heap {
	struct mst_heap_node *nodes;
	int n_nodes;
};

struct mst_heap_node {
	int weight;
	int x;
	int y;
};

struct mst_ubr_node {
	struct mst_ubr_node *parent;
	int rank;
};

static struct mst_ubr_node *mst_find(struct mst_ubr_node *n)
{
	if (n->parent != n)
		n->parent = mst_find(n->parent);

	return n->parent;
}

static void mst_union(struct mst_ubr_node *x, struct mst_ubr_node *y)
{
	struct mst_ubr_node *rx = mst_find(x);
	struct mst_ubr_node *ry = mst_find(y);

	if (rx == ry)
		return;

	rx->parent = ry;

	if (rx->rank == ry->rank)
		ry->rank++;
}

int mst_heapsort_cmp(const void *x, const void *y)
{
	return ((struct mst_heap_node *)x)->weight - ((struct mst_heap_node *)y)->weight;
}

static struct mst_heap *mst_heapsort(struct coordinate *locs, int n_locs)
{
	struct mst_heap *h = malloc(sizeof(struct mst_heap));
	h->n_nodes = n_locs * (n_locs - 1) / 2;
	h->nodes = calloc(h->n_nodes, sizeof(struct mst_heap_node));

	/* create weight matrix; x > y */
	int c = 0;
	for (int y = 0; y < n_locs; y++) {
		for (int x = y + 1; x < n_locs; x++) {
			struct mst_heap_node n = {distance_cityblock(locs[x], locs[y]), x, y};
			h->nodes[c++] = n;
		}
	}

	heapsort(h->nodes, h->n_nodes, sizeof(struct mst_heap_node), mst_heapsort_cmp);

	return h;
}

static struct segments *create_mst(struct coordinate *locs, int n_locs)
{
	struct segments *mst = malloc(sizeof(struct segments));
	mst->n_segments = n_locs - 1;
	mst->segments = calloc(sizeof(struct segment), (mst->n_segments));

	struct mst_heap *h = mst_heapsort(locs, n_locs);

	struct mst_ubr_node *s = calloc(sizeof(struct mst_ubr_node), n_locs);
	// make-set
	for (int i = 0; i < n_locs; i++) {
		s[i].parent = &s[i];
		s[i].rank = 0;
	}

	int n_segments = 0;
	for (int i = 0; i < h->n_nodes; i++) {
		int x = h->nodes[i].x;
		int y = h->nodes[i].y;
		if (mst_find(&s[x]) != mst_find(&s[y])) {
			struct segment new_seg = {locs[x], locs[y]};
			mst->segments[n_segments++] = new_seg;
			mst_union(&s[x], &s[y]);
			if (n_segments >= mst->n_segments)
				break;
		}
	}

	return mst;
}

static void free_struct_segments(struct segments *s)
{
	free(s->segments);
	free(s);
}

static char *overlap_tmp = NULL;
static int overlap_tmp_size = 0;
static int overlap_penalty = 0;

/* determine the number of overlaps in the grid */
static int compute_overlap_penalty(struct cell_placements *cp)
{
	int i, penalty;
	int size;
	struct placement p;

	struct dimensions d = compute_placement_dimensions(cp);
#ifdef PLACER_SCORE_DEBUG
	printf("[compute_overlap_penalty] l=%d, w=%d, h=%d\n", d.x, d.z, d.y);
#endif

	penalty = 0;
	size = d.x * d.y * d.z;
	if (!overlap_tmp || size > overlap_tmp_size) {
		int new_size = size * 2;
		printf("[compute_overlap_penalty] resizing from %d to %d (%d x %d x %d)\n", overlap_tmp_size, new_size, d.x, d.z, d.y);
		free(overlap_tmp);
		overlap_tmp = malloc(new_size * sizeof(char));
		overlap_tmp_size = new_size;
	}
	memset(overlap_tmp, 0, size);

	for (i = 0; i < cp->n_placements; i++) {
		p = cp->placements[i];

		struct coordinate c = p.placement;
		struct dimensions pd = p.cell->dimensions;

		int cell_x = c.x + pd.x;
		int cell_y = c.y + pd.y;
		int cell_z = c.z + pd.z;

/*
		printf("[compute_overlap_penalty] placing coordinate (y=%d, z=%d, x=%d) + (h=%d, w=%d, l=%d)\n",
			c.y, c.z, c.x, pd.y, pd.z, pd.x);
*/

		/* if another cell is there, increase the penalty */
		for (int y = c.y; y < cell_y; y++) {
			int by = y * d.z * d.x;
			for (int z = max(0, c.z - 1); z < min(d.z, cell_z + 1); z++) {
				int bz = z * d.x;
				for (int x = max(0, c.x - 1); x < min(d.x, cell_x + 1); x++) {
					// printf("[compute_overlap_penalty] y=%d z=%d x=%d\n", y, z, x);
					if (overlap_tmp[by + bz + x]++)
						penalty += 1;
				}
			}
		}
	}

	/* read by simulated_annealing_placement to see if we must continue */
	overlap_penalty = penalty;
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
	int penalty = 0;
	struct coordinate **coords = calloc(sizeof(struct coordinate *), cp->n_nets + 1);
	int *n_coords = calloc(sizeof(int), cp->n_nets + 1);
	int *found = calloc(sizeof(int), cp->n_nets + 1);

	for (i = 1; i < cp->n_nets + 1; i++) {
		found[i] = 0;
		n_coords[i] = PENALTY_MALLOC_START;
		coords[i] = calloc(sizeof(struct coordinate), n_coords[i]);
	}

	/* add each coordinate to its net */
	for (i = 0; i < cp->n_placements; i++) {
		struct placement pl = cp->placements[i];
		struct logic_cell *c = pl.cell;
		for (int j = 0; j < c->n_pins; j++) {
			net_t k = pl.nets[j];

			struct coordinate b = pl.placement;
			struct coordinate p = c->pins[j].coordinate;

			struct coordinate actual = {b.x + p.x, b.y + p.y, b.z + p.z};
			coords[k][found[k]++] = actual;

			/* resize coords if too large */
			if (found[k] >= n_coords[k]) {
				n_coords[k] *= 2;
				coords[k] = realloc(coords[k], sizeof(struct coordinate) * n_coords[k]);
			}
		}
	}

	/* for each net, compute the constituent pin coordinates */
	for (i = 1; i < cp->n_nets; i++) {
		if (found[i] == 0 || found[i] == 1)
			continue;

		if (found[i] == 2) {
			penalty += distance_cityblock(coords[i][0], coords[i][1]);
			continue;
		}

		struct segments *mst = create_mst(coords[i], found[i]);
		for (int seg = 0; seg < mst->n_segments; seg++) {
			int d = distance_cityblock(mst->segments[seg].start, mst->segments[seg].end);
			penalty += d;
		}
		free_struct_segments(mst);
	}


	free(found);
	free(n_coords);
	for (int i = 1; i < cp->n_nets + 1; i++)
		free(coords[i]);
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
		struct coordinate c = placements->placements[i].placement;
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

	int match_iterations = 0;
	int match_score = old_score;
	int stop_iterations = 3;

	for (i = 0; i < iterations || match_iterations < stop_iterations || overlap_penalty > 0; i++) {
#ifdef PLACER_GENERATION_DEBUG
		printf("[placer] iteration = %d\n", i);
#endif
		method = DISPLACE;

		for (g = 0; g < generations; g++) {
#ifdef PLACER_GENERATION_DEBUG
			printf("[placer] generation = %d\n", g);
#endif
			/* make a copy of these placements */
			new_placements = copy_placements(best_placements);
			// print_cell_placements(new_placements);
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

			if (taken_score == match_score) {
				match_iterations++;
			} else {
				match_score = taken_score;
				match_iterations = 0;
			}

			old_score = taken_score;
		}
		// printf("[placer] T = %4.2f\n", t);

		printf("Iteration: %d, Score: %d (overlap penalty: %d), Temperature: %4.2f\n", i, taken_score, overlap_penalty, t);
		// print_cell_placements(best_placements);

		t = update(t, fixed_alpha);
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
		c = &(cl->cells[i]);
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
		// printf("[placer] cell %d\n", i);
		struct blif_cell *c = blif->cells[i];
		struct placement p;
		p.nets = malloc(sizeof(net_t) * c->n_pins);


		p.cell = map_cell_to_library(c, cl);
		if (!p.cell)
			printf("[placer] could not map blif cell (%s) to cell library\n", c->name);

		/* copy blif pins net to p->nets */
		for (int j = 0; j < c->n_pins; j++) {
			// printf("[placer] placements[%d](%s)->nets[%d] = c->pins[%d]->net = %u\n", i, c->name, j, j, c->pins[j]->net);
			p.nets[j] = c->pins[j]->net;
		}

		p.placement.x = 0;
		p.placement.y = 0;
		p.placement.z = 0;

		p.turns = 0;

		placements->placements[i] = p;
	}

	return placements;
}

void print_cell_placements(struct cell_placements *cp)
{
	int i;
	for (i = 0; i < cp->n_placements; i++) {
		struct placement p = cp->placements[i];
		printf("[placer] placement: %s @ (y=%d, z=%d, x=%d), %lu turns\n",
			p.cell->name, p.placement.y, p.placement.z, p.placement.x, p.turns);
	}
}

/* create an initial placement */
struct cell_placements *placer_initial_place(struct blif *blif, struct cell_library *cl)
{
	struct cell_placements *cp = map_blif_to_cell_library(blif, cl);

	int margin = 6;
	int w = roundl(sqrt(cp->n_placements));

	for (int i = 0, x = 0; i < cp->n_placements; i++, x += margin) {
		cp->placements[i].placement.x = x % (margin * w);
		cp->placements[i].placement.z = x / (margin * w);
	}

	return cp;
}
