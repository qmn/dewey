#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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
enum placement_method generate(struct cell_placements *placements,
		struct dimensions *dimensions,
		double t, double t_0,
		enum placement_method method)
{
	unsigned long cell_a_idx, cell_b_idx;
	struct placement *cell_a, *cell_b, *tmp;
	double p, interchange_threshold;
	enum placement_method method_used;
	double scaling_factor;
	long window_height, window_width;

	/* select a random cell to interchange, displace, or reorient */
	cell_a_idx = (unsigned long)random() % placements->length;
	cell_a = placements->placements[cell_a_idx];

	/* compute the probabilty we change this cell */
	p = (double)random() / (double)(RAND_MAX);
	interchange_threshold = (1.0 / DISPLACE_INTERCHANGE_RATIO);

	if (p > interchange_threshold) {
		/* select another cell */
		cell_b_idx = cell_a_idx;
		while (cell_b_idx == cell_a_idx)
			cell_b_idx = (unsigned long)random() % placements->length;
		
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
	int i, j;

	new_placements = (struct cell_placements *)malloc(sizeof(struct cell_placements));

	new_placements->length = old_placements->length;
	p = (struct placement **)malloc(old_placements->length * sizeof(struct placement *));
	o = old_placements->placements;

	/* deep copy each placement */
	for (i = 0; i < old_placements->length; i++) {
		p[i]->name = o[i]->name;
		p[i]->placement = o[i]->placement;
		p[i]->turns = o[i]->turns;

		/* pins, being sourced from the BLIF object, won't change */
		p[i]->n_pins = o[i]->n_pins;
		p[i]->pins = o[i]->pins;
	}

	new_placements->placements = p;

	return new_placements;
}

static void free_placements(struct cell_placements *placements)
{
	free(placements->placements);
	free(placements);
}

static int score(struct cell_placements *placements, struct dimensions *dimensions)
{
	return 1;
}

static int accept(int new_score, int old_score, double t)
{
	return 1;
}

static double update(double t)
{
	return t - 1.0;
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
				free(best_placements);
				best_placements = new_placements;
				taken_score = new_score;
				if (method_used == REORIENT)
					method = DISPLACE;
			} else {
				/* reject the new placement */
				free(new_placements);
				taken_score = old_score;
				if (method_used == DISPLACE)
					method = REORIENT;
			}
		}
		t = update(t);

		printf("Iteration: %d, Score: %d\n", i, taken_score);
	}

	return best_placements;
}
