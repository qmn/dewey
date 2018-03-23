#include <stdio.h>
#include <string.h>

#include "placer.h"
#include "blif.h"
#include "base_router.h"
#include "extract.h"

static void write_out_placement_array(FILE *f, struct placement p)
{
	struct coordinate c = p.placement;
	struct logic_cell *lc = p.cell;
	struct dimensions d = lc->dimensions[p.turns];

	for (int y = 0; y < d.y; y++) {
		int yy = y * d.z * d.x;
		for (int z = 0; z < d.z; z++) {
			int zz = z * d.x;
			for (int x = 0; x < d.x; x++) {
				int last = (y == d.y-1) && (z == d.z-1) && (x == d.x-1);
				int b_off = yy + zz + x; // block offset into logic cell
				block_t b = lc->blocks[p.turns][b_off];
				data_t da = lc->data[p.turns][b_off];

				fprintf(f, "[%d, %d, %d, %d, %d]", c.y, c.z, c.x, b, da);
				if (!last)
					fputc(',', f);
			}
		}
	}
}

static void vis_json_placements(FILE *f, struct cell_placements *cp)
{
	fprintf(f, "[\n");
	for (int i = 0; i < cp->n_placements; i++) {
		struct placement p = cp->placements[i];
		fprintf(f, "  {\"name\": \"%s\",\n", p.cell->name);
		fprintf(f, "   \"blocks\": [");
		write_out_placement_array(f, p);
		fprintf(f, "]}\n");
		if (i < cp->n_placements-1)
			fputc(',', f);
	}
	fprintf(f, "]\n");
}

static void vis_json_routings(FILE *f, struct blif *blif, struct routings *rt)
{
	struct coordinate disp = {0, 0, 0};
	fprintf(f, "[\n");
	for (net_t i = 1; i < rt->n_routed_nets; i++) {
		struct extracted_net *en = extract_net(&rt->routed_nets[i], disp);
		fprintf(f, "  {\"name\": \"");
		char *n = get_net_name(blif, i);
		for (int j = 0; j < strlen(n); j++)
			if (n[j] == '\\')
				fprintf(f, "\\\\");
			else
				fputc(n[j], f);
		fprintf(f, "\",\n   \"blocks\": [");
		for (int j = 0; j < en->n; j++) {
			fprintf(f, "[%d, %d, %d, %d, %d]", PRINT_COORD(en->c[j]), en->b[j], en->d[j]);
			if (j < en->n-1)
				fputc(',', f);
		}
		fprintf(f, "]}");
		if (i < rt->n_routed_nets - 1)
			fputc(',', f);
		fputc('\n', f);
	}
	fprintf(f, "]\n");
}

static struct dimensions vis_json_dimensions(struct cell_placements *cp, struct routings *rt)
{
	struct coordinate disp = placements_top_left_most_point(cp);
	if (rt)
		disp = coordinate_piecewise_min(disp, routings_top_left_most_point(rt));

	struct dimensions cpd = compute_placement_dimensions(cp);
	struct dimensions rtd = {0, 0, 0};
	if (rt)
		rtd = compute_routings_dimensions(rt);
	struct dimensions d = dimensions_piecewise_max(cpd, rtd);
	d.y += disp.y;
	d.x += disp.x;
	d.z += disp.z;

	return d;
}

void vis_json(FILE *f, struct blif *blif, struct cell_placements *cp, struct routings *rt)
{
	struct dimensions d = vis_json_dimensions(cp, rt);
	struct coordinate disp = placements_top_left_most_point(cp);
	if (rt)
		disp = coordinate_piecewise_min(disp, routings_top_left_most_point(rt));
	fprintf(f, "{\"dimensions\": [%d, %d, %d],\n \"disp\": [%d, %d, %d],\n", PRINT_COORD(d), PRINT_COORD(disp));
	fprintf(f, " \"placements\":\n");
	vis_json_placements(f, cp);
	fprintf(f, ",\n \"routings\":\n");
	vis_json_routings(f, blif, rt);
	fprintf(f, "}");
}
