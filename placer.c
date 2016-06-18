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
	struct placement *cell_a, *tmp;
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
		cell_b_idx = cell_a_idx;
		while (cell_b_idx == cell_a_idx)
			cell_b_idx = (unsigned long)random() % placements->n_placements;
		
		/* interchange the cells */
		tmp = placements->placements[cell_b_idx];
		placements->placements[cell_b_idx] = cell_a;
		placements->placements[cell_a_idx] = tmp;

		return INTERCHANGE;
	} else if (method == DISPLACE) {
		scaling_factor = log(t) / log(t_0);

		/* figure the most this placement can move */
		window_height = lround(dimensions->z * scaling_factor);
		window_width = lround(dimensions->x * scaling_factor);

		window_height = max(window_height, MIN_WINDOW_HEIGHT);
		window_width = max(window_width, MIN_WINDOW_WIDTH);

		/* displace */
		cell_a->placement.z += random() % (window_height * 2) - window_height;
		cell_a->placement.x += random() % (window_width * 2) - window_width;

		return DISPLACE;
	} else {
		/* reorient */
		cell_a->turns = (cell_a->turns + 1 % 4);
		return REORIENT;
	}
}

/*
 * Performs a deep copy of cell_placements.
 */
static struct cell_placements *copy_placements(struct cell_placements *old_placements)
{
	struct cell_placements *new_placements;
	struct placement **p, **o;
	int i;

	new_placements = (struct cell_placements *)malloc(sizeof(struct cell_placements));

	new_placements->n_placements = old_placements->n_placements;
	p = (struct placement **)malloc(old_placements->n_placements * sizeof(struct placement *));
	o = old_placements->placements;

	/* deep copy each placement */
	for (i = 0; i < old_placements->n_placements; i++) {
		p[i]->cell = o[i]->cell;
		p[i]->placement = o[i]->placement;
		p[i]->turns = o[i]->turns;
	}

	new_placements->placements = p;

	return new_placements;
}

void free_placements(struct cell_placements *placements)
{
	int i;

	for (i = 0; i < placements->n_placements; i++) {
		free(placements->placements[i]);
	}

	free(placements);
}

static int score(struct cell_placements *placements, struct dimensions *dimensions)
{
	return 1;
}

static int accept(int new_score, int old_score, double t)
{
	double ratio, acceptance_criterion;
	ratio = (double)(new_score - old_score) / t;
	if (ratio > 1)
		return 1;

	acceptance_criterion = fmin(1.0, exp(ratio));
	return random() < (long)(acceptance_criterion * RAND_MAX);
}

static double update(double t, double (*alpha)(double))
{
	return t * alpha(t);
}

static double fixed_alpha(double t)
{
	return 0.9;
}

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

	t = t_0;
	best_placements = initial_placements;
	old_score = score(initial_placements, dimensions);

	for (i = 0; i < iterations; i++) {
		method = DISPLACE;

		for (g = 0; g < generations; g++) {
			/* make a copy of these placements */
			new_placements = copy_placements(best_placements);
			method_used = generate(new_placements, dimensions, t, t_0, method);
			new_score = score(new_placements, dimensions);

			if (accept(new_score, old_score, t)) {
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
				/* reject the new placement */
				free_placements(new_placements);
				taken_score = old_score;
				if (method_used == DISPLACE)
					method = REORIENT;
			}
		}
		t = update(t, fixed_alpha);

		printf("Iteration: %d, Score: %d\n", i, taken_score);
	}

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
	placements->n_placements = 0;

	/* for now, only place cells from the blif -- no inputs or outputs */
	for (i = 0; i < blif->n_cells; i++) {
		struct blif_cell *c = blif->cells[i];
		struct placement *p = malloc(sizeof(struct placement));

		p->cell = map_cell_to_library(c, cl);
		if (!p->cell)
			printf("[placer] could not map blif cell (%s) to cell library\n", c->name);

		p->placement.x = 0;
		p->placement.y = 0;
		p->placement.z = 0;

		p->turns = 0;

		/* extend the placements and insert the newly made placement at the end */
		if (placements->n_placements++) {
			placements->placements = realloc(placements->placements, sizeof(struct placement *) * placements->n_placements);
		} else {
			placements->placements = malloc(sizeof(struct placement *));
		}
		placements->placements[placements->n_placements - 1] = p;
	}

	return placements;
}

void print_cell_placements(struct cell_placements *cp)
{
	int i;
	for (i = 0; i < cp->n_placements; i++) {
		struct placement *p = cp->placements[i];
		printf("[placer] placement: %s @ (%d, %d, %d), %lu turns\n",
			p->cell->name, p->placement.y, p->placement.z, p->placement.x, p->turns);
	}
}

struct cell_placements *placer_initial_place(struct blif *blif, struct cell_library *cl)
{
	return map_blif_to_cell_library(blif, cl);
}
