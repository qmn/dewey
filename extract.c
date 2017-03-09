#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "extract.h"
#include "placer.h"
#include "router.h"

static struct coordinate placements_top_left_most_point(struct cell_placements *cp)
{
	/* select the first placement for a baseline point */
	struct coordinate d = cp->placements[0].placement;

	for (int i = 0; i < cp->n_placements; i++) {
		struct coordinate c = cp->placements[i].placement;

		d = coordinate_piecewise_min(c, d);
	}

	return d;
}

/* determine the top-left most point of routings */
static struct coordinate routings_top_left_most_point(struct routings *rt)
{
	/* select a baseline point */
	struct coordinate d = rt->routed_nets[1].routed_segments[0].coords[0];

	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++) {
			for (int k = 0; k < rt->routed_nets[i].routed_segments[j].n_coords; k++) {
				struct coordinate c = rt->routed_nets[i].routed_segments[j].coords[k];
				d = coordinate_piecewise_min(c, d);
			}
		}
	}

	return d;
}

/* move the entire design so that all coordinates are non-negative:
 * if routings are provided, consider any possible out-of-bounds routing as well,
 * and recenter the routings as well. */
void recenter(struct cell_placements *cp, struct routings *rt)
{
	struct coordinate disp = placements_top_left_most_point(cp);
	if (rt)
		disp = coordinate_piecewise_min(disp, routings_top_left_most_point(rt));

	placements_displace(cp, coordinate_neg(disp));

	if (rt)
		routings_displace(rt, coordinate_neg(disp));
}

struct extraction *extract(struct cell_placements *cp, struct routings *rt)
{
	struct cell_placements *ncp = copy_placements(cp);
	struct routings *nrt = rt ? copy_routings(rt) : NULL;
	recenter(ncp, nrt);

	struct dimensions cpd = compute_placement_dimensions(ncp);
	struct dimensions rtd = {0, 0, 0};
	if (nrt)
		rtd = compute_routings_dimensions(rt);
	struct dimensions d = dimensions_piecewise_max(cpd, rtd);
	assert(d.x > 0 && d.y > 0 && d.z > 0);
	printf("[extract] extraction dimensions: h=%d w=%d l=%d\n", d.y, d.z, d.x);

	struct extraction *e = malloc(sizeof(struct extraction));
	e->dimensions = d;

	int size = d.x * d.y * d.z;
	e->blocks = calloc(size, sizeof(block_t));
	e->data = calloc(size, sizeof(data_t));

	/* place blocks resulting from placement in image */
	for (int i = 0; i < cp->n_placements; i++) {
		struct placement p = cp->placements[i];

		struct coordinate c = p.placement;
		struct logic_cell *lc = p.cell;
		struct dimensions lcd = lc->dimensions[p.turns];

		for (int y = 0; y < lcd.y; y++) {
			for (int z = 0; z < lcd.z; z++) {
				for (int x = 0; x < lcd.x; x++) {
					int fd_off = (c.y + y) * d.z * d.x + (c.z + z) * d.x + (c.x + x);
					int b_off = y * lcd.z * lcd.x + z * lcd.x + x;

					assert(fd_off >= 0 && fd_off <= size);
					e->blocks[fd_off] = lc->blocks[p.turns][b_off];
					e->data[fd_off] = lc->data[p.turns][b_off];
				}
			}
		}
	}

	free_cell_placements(ncp);

	/* if that's all, return, otherwise proceed to add routings */
	if (!rt)
		return e;

	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++) {
			for (int k = 0; k < rt->routed_nets[i].routed_segments[j].n_coords; k++) {
				struct coordinate c = rt->routed_nets[i].routed_segments[j].coords[k];
				if (c.x > d.x || c.y > d.y || c.z > d.z || c.x < 0 || c.y < 0 || c.z < 0)
					continue;
				e->blocks[c.y * d.z * d.x + c.z * d.x + c.x] = 55;
			}
		}
	}

	free_routings(nrt);

	return e;
}

void free_extraction(struct extraction *e)
{
	free(e->blocks);
	free(e->data);
	free(e);
}
