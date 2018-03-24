#include "usage_matrix.h"
#include "util.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

int in_usage_bounds(struct usage_matrix *m, struct coordinate c)
{
	return c.x >= 0 && c.x < m->d.x &&
	       c.y >= 0 && c.y < m->d.y &&
	       c.z >= 0 && c.z < m->d.z;
}

static int is_toplevel_pin(struct placement *p)
{
	if (strncmp("input_pin", p->cell->name, 10) == 0 ||
	    strncmp("output_pin", p->cell->name, 11) == 0)
		return 1;
	return 0;
}

// marks a usage matrix so that no routing will occur on Y=0 to the margin
static void tunnel_to_margin(struct placement *p, struct usage_matrix *m)
{
	int sx = 0, ex = 0;
	int sz = p->placement.z - 1, ez = p->placement.z + 1;
	if (p->constraints & CONSTR_KEEP_LEFT) {
		sx = 0;
		ex = p->placement.x;
	} else if (p->constraints & CONSTR_KEEP_RIGHT) {
		sx = p->placement.x;
		ex = m->d.x;
	}

	for (int x = sx; x < ex; x++)
		for (int z = sz; z < ez; z++)
			usage_mark(m, (struct coordinate){0, z, x});
}

/* create a usage_matrix that marks where blocks from existing cell placements
   and routed nets occupy the grid. */
struct usage_matrix *create_usage_matrix(struct cell_placements *cp, struct routings *rt, int xz_margin)
{
	struct coordinate tlcp = placements_top_left_most_point(cp);
	struct coordinate tlrt = routings_top_left_most_point(rt);
	struct coordinate top_left_most = coordinate_piecewise_min(tlcp, tlrt);

	assert(top_left_most.x >= xz_margin && top_left_most.y >= 0 && top_left_most.z >= xz_margin);

	// size the usage matrix and allow for routing on y=0 and y=3
	struct dimensions d = dimensions_piecewise_max(compute_placement_dimensions(cp), compute_routings_dimensions(rt));
	assert(d.x > 0 && d.x < 1000 && d.z > 0 && d.z < 1000);
	d.y = max(d.y, 7);
	d.z += xz_margin;
	d.x += xz_margin;
	// printf("[usage_matrix] size is %dx%dx%d\n", d.y, d.z, d.x);

	struct usage_matrix *m = malloc(sizeof(struct usage_matrix));
	m->d = d;
	m->xz_margin = xz_margin;
	m->matrix = malloc(d.x * d.y * d.z * sizeof(unsigned char));
	memset(m->matrix, 0, d.x * d.y * d.z * sizeof(unsigned char));

	/* placements */
	for (int i = 0; i < cp->n_placements; i++) {
		struct placement p = cp->placements[i];
		struct coordinate c = p.placement;
		struct dimensions pd = p.cell->dimensions[p.turns];

		int cell_x = c.x + pd.x;
		int cell_y = c.y + pd.y;
		int cell_z = c.z + pd.z;

		int z1 = max(0, c.z), z2 = min(d.z, cell_z);
		int x1 = max(0, c.x), x2 = min(d.x, cell_x);
		int y2 = min(cell_y + 1, d.y);

		for (int y = c.y; y < y2; y++) {
			for (int z = z1; z < z2; z++) {
				for (int x = x1; x < x2; x++) {
					struct coordinate cc = {y, z, x};
					usage_mark(m, cc);
				}
			}
		}

		// tunnel to margin
		if (is_toplevel_pin(&p))
			tunnel_to_margin(&p, m);
	}


	/* routings */
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		for (struct routed_segment_head *rsh = rt->routed_nets[i].routed_segments; rsh; rsh = rsh->next) {
			struct coordinate c = rsh->rseg.seg.end;
			for (int k = 0; k < rsh->rseg.n_backtraces; k++) {
				// here
				c = disp_backtrace(c, rsh->rseg.bt[k]);
				usage_mark(m, c);

				// below
				struct coordinate c2 = c;
				c2.y--;
				if (in_usage_bounds(m, c2))
					usage_mark(m, c2);
			}
		}
	}

	return m;
}

int usage_matrix_violated(struct usage_matrix *m, struct coordinate c)
{
	// [yzx]2 DOES include that coordinate
	int y1 = max(c.y - 1, 0), y2 = min(c.y, m->d.y - 1);
	int z1 = max(c.z - 1, 0), z2 = min(c.z + 1, m->d.z - 1);
	int x1 = max(c.x - 1, 0), x2 = min(c.x + 1, m->d.x - 1);

	for (int y = y1; y <= y2; y++) {
		int dy = y * m->d.z * m->d.x;
		for (int z = z1; z <= z2; z++) {
			int dz = z * m->d.x;
			for (int x = x1; x <= x2; x++) {
				int idx = dy + dz + x;
				if (m->matrix[idx])
					return 1;
			}
		}
	}

	return 0;
}

