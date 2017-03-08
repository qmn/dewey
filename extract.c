#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "extract.h"
#include "placer.h"
#include "router.h"

/* produces a 3D array representing actual Minecraft block placements */
struct extraction *extract_placements(struct cell_placements *cp)
{
	struct extraction *e = malloc(sizeof(struct extraction));
	struct dimensions d = compute_placement_dimensions(cp);
	e->dimensions = d;

	int size = d.x * d.y * d.z;
	e->blocks = calloc(size, sizeof(block_t));
	e->data = calloc(size, sizeof(data_t));

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

	return e;
}

struct extraction *extract(struct cell_placements *cp, struct routings *rt)
{
	struct extraction *e = extract_placements(cp);

	if (!rt)
		return e;

	struct dimensions d = compute_placement_dimensions(cp);

	for (net_t i = 1; i < rt->n_routed_nets; i++) {
		for (int j = 0; j < rt->routed_nets[i].n_coords; j++) {
			struct coordinate c = rt->routed_nets[i].coords[j];
			if (c.x > d.x || c.y > d.y || c.z > d.z || c.x < 0 || c.y < 0 || c.z < 0)
				continue;
			e->blocks[c.y * d.z * d.x + c.z * d.x + c.x] = 55;
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
