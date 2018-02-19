#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "extract.h"
#include "placer.h"
#include "router.h"

struct coordinate placements_top_left_most_point(struct cell_placements *cp)
{
	/* select the first placement for a baseline point */
	int unfound = 1;
	struct coordinate d;

	for (int i = 0; i < cp->n_placements; i++) {
		struct coordinate c = cp->placements[i].placement;
		if (unfound) {
			d = c;
			unfound = 0;
		} else {
			d = coordinate_piecewise_min(c, d);
		}
		// printf("[ptlmp] c = (%d, %d, %d)\n", c.y, c.z, c.x);
	}

	return d;
}

/* determine the top-left most point of routings */
struct coordinate routings_top_left_most_point(struct routings *rt)
{
	int started = 0;
	struct coordinate d;

	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		for (struct routed_segment_head *rsh = rt->routed_nets[i].routed_segments; rsh; rsh = rsh->next) {
			struct coordinate c = rsh->rseg.seg.end;
			if (!started) {
				started = 1;
				d = c;
			} else {
				d = coordinate_piecewise_min(c, d);
			}

			for (int k = 0; k < rsh->rseg.n_backtraces; k++) {
				c = disp_backtrace(c, rsh->rseg.bt[k]);
				d = coordinate_piecewise_min(c, d);
			}

			d = coordinate_piecewise_min(rsh->rseg.seg.start, d);
		}
	}

	printf("[tlmp] (%d, %d, %d)\n", d.y, d.z, d.x);

	return d;
}

/* move the entire design so that all coordinates are non-negative:
 * if routings are provided, consider any possible out-of-bounds routing as well,
 * and recenter the routings as well. */
void recenter(struct cell_placements *cp, struct routings *rt, int xz_margin)
{
	struct coordinate disp = placements_top_left_most_point(cp);
	if (rt)
		disp = coordinate_piecewise_min(disp, routings_top_left_most_point(rt));

	struct coordinate xz_add = {0, -xz_margin, -xz_margin};
	disp = coordinate_add(disp, xz_add);

	placements_displace(cp, coordinate_neg(disp));

	if (rt)
		routings_displace(rt, coordinate_neg(disp));
}

void place_block(struct extraction *e, struct coordinate c)
{
	c.z += 2;
	c.x += 2;

	struct dimensions d = e->dimensions;

	if (c.x > d.x || c.y > d.y || c.z > d.z || c.x < 0 || c.y < 0 || c.z < 0)
		return;

	e->blocks[c.y * d.z * d.x + c.z * d.x + c.x] = 55;

	if (c.y == 3)
		e->blocks[(c.y - 1) * d.z * d.x + c.z * d.x + c.x] = 5;
}

struct extraction *extract(struct cell_placements *cp, struct routings *rt)
{
	struct coordinate disp = placements_top_left_most_point(cp);
	if (rt)
		disp = coordinate_piecewise_min(disp, routings_top_left_most_point(rt));

	struct dimensions cpd = compute_placement_dimensions(cp);
	struct dimensions rtd = {0, 0, 0};
	if (rt)
		rtd = compute_routings_dimensions(rt);
	struct dimensions d = dimensions_piecewise_max(cpd, rtd);

	// modify dimensions (we'll perform the displacement later)
	d.z += -disp.z + 4;
	d.x += -disp.x + 4;
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
					int fd_off = (c.y - disp.y + 2 + y) * d.z * d.x + (c.z - disp.z + 2 + z) * d.x + (c.x - disp.x + 2 + x);
					int b_off = y * lcd.z * lcd.x + z * lcd.x + x;

					assert(fd_off >= 0 && fd_off <= size);
					e->blocks[fd_off] = lc->blocks[p.turns][b_off];
					e->data[fd_off] = lc->data[p.turns][b_off];
				}
			}
		}
	}

	/* if that's all, return, otherwise proceed to add routings */
	if (!rt)
		return e;

	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		for (struct routed_segment_head *rsh = rt->routed_nets[i].routed_segments; rsh; rsh = rsh->next) {
			struct coordinate c = rsh->rseg.seg.end;
			for (int k = 0; k < rsh->rseg.n_backtraces; k++) {
				c = disp_backtrace(c, rsh->rseg.bt[k]);
				struct coordinate cc = coordinate_sub(c, disp);

				place_block(e, cc);
			}

			place_block(e, coordinate_sub(rsh->rseg.seg.start, disp));
			place_block(e, coordinate_sub(rsh->rseg.seg.end, disp));
		}
	}

	return e;
}

void free_extraction(struct extraction *e)
{
	free(e->blocks);
	free(e->data);
	free(e);
}
