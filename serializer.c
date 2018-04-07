#include <stdio.h>
#include <stdlib.h>

#include "placer.h"
#include "blif.h"
#include "serializer.h"
#include "base_router.h"

static void serialize_placement(FILE *f, struct placement *p, struct blif *blif)
{
	struct coordinate c = p->placement;
	int t = p->turns;
	// char *n = get_net_name(blif, p->net);

	fprintf(f, "  %s@%s:\n", p->cell->name, p->cell->lib->fn);
	fprintf(f, "    placement: (%d, %d, %d)\n", c.y, c.z, c.x);
	fprintf(f, "    turns: %d\n", t);
	fprintf(f, "    nets: [");
	for (int i = 0; i < p->cell->n_pins; i++)
		fprintf(f, "%s%s", get_net_name(blif, p->nets[i]), i != p->cell->n_pins - 1 ? ", " : "");
	fprintf(f, "]\n");
	fprintf(f, "    constraints: 0x%lx\n", p->constraints);
	fprintf(f, "    margin: %d\n", p->margin);
}

void serialize_placements(FILE *f, struct cell_placements *cp, struct blif *blif)
{
	fprintf(f, "placements:\n");
	for (int i = 0; i < cp->n_placements; i++) {
		struct placement *p = &cp->placements[i];
		serialize_placement(f, p, blif);
	}
}

static char serialize_backtrace(enum backtrace bt)
{
	switch (bt) {
	case BT_UP:
		return 'U';
	case BT_DOWN:
		return 'D';
	case BT_NORTH:
		return 'N';
	case BT_WEST:
		return 'W';
	case BT_EAST:
		return 'E';
	case BT_SOUTH:
		return 'S';
	default:
		return '_';
	}
}

static void serialize_routing(FILE *f, struct routed_net *rn)
{
	for (struct routed_segment_head *rsh = rn->routed_segments; rsh; rsh = rsh->next) {
		struct routed_segment rseg = rsh->rseg;
		struct coordinate s = rseg.seg.start, e = rseg.seg.end;
		fprintf(f, "  - (%d, %d, %d) -> (%d, %d, %d):\n", s.y, s.z, s.x, e.y, e.z, e.x);
		fprintf(f, "    backtrace: [");
		for (int i = 0; i < rseg.n_backtraces; i++)
			fprintf(f, "%c%s", serialize_backtrace(rseg.bt[i]), i < rseg.n_backtraces - 1 ? ", " : "");
		fprintf(f, "]\n");
	}
}

void serialize_routings(FILE *f, struct routings *rt, struct blif *blif)
{
	fprintf(f, "routings:\n");
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		struct routed_net *rn = &rt->routed_nets[i];
		char *n = get_net_name(blif, i);
		fprintf(f, "  %s:\n", n);
		serialize_routing(f, rn);
	}
}

void serialize_extraction(FILE *f, struct extraction *e)
{
	struct dimensions d = e->dimensions;

	fprintf(f, "extraction:\n");
	fprintf(f, "  dimensions: (%d, %d, %d)\n", d.y, d.z, d.x);
	fprintf(f, "  blocks: [");
	for (int i = 0; i < d.y * d.z * d.x; i++) {
		fprintf(f, "%d", e->blocks[i]);

		if (i + 1 < d.y * d.z * d.x)
			fprintf(f, ",");
		else
			fprintf(f, "]\n");
	}

	fprintf(f, "  data: [");
	for (int i = 0; i < d.y * d.z * d.x; i++) {
		fprintf(f, "%d", e->data[i]);
		if (i + 1 < d.y * d.z * d.x)
			fprintf(f, ",");
		else
			fprintf(f, "]\n");
	}
}
