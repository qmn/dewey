#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "extract.h"
#include "placer.h"
#include "router.h"
#include "base_router.h"

// mass movement routines
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

// within connection distace
static int adjacent(struct coordinate c1, struct coordinate c2)
{
	int dy = abs(c1.y - c2.y);
	int dz = abs(c1.z - c2.z);
	int dx = abs(c1.x - c2.x);

	return !coordinate_equal(c1, c2) && (dx <= 1 && dz <= 1) && (dy == 0 || dy == 3);
}

struct extracted_net {
	net_t net;

	int n;
	int sz;
	struct coordinate *c;
	block_t *b;
	data_t *d;
};

// extend the extracted net to fit new_size elements
static void extracted_net_append(struct extracted_net *en, struct coordinate c, block_t b, data_t d)
{
	if (en->n >= en->sz) {
		en->sz *= 2;
		en->c = realloc(en->c, sizeof(struct coordinate) * en->sz);
		en->b = realloc(en->b, sizeof(block_t) * en->sz);
		en->d = realloc(en->d, sizeof(data_t) * en->sz);
	}

	int i = en->n++;
	en->c[i] = c;
	en->b[i] = b;
	en->d[i] = d;
}

static void reverse_segment(struct routed_segment *rseg)
{
	rseg->seg = (struct segment){rseg->seg.end, rseg->seg.start};
	invert_backtrace_sequence(rseg->bt, rseg->n_backtraces);
}

enum movement {
	GO_NONE   = 0,
	GO_WEST   = 1 << 0,
	GO_SOUTH  = 1 << 1,
	GO_EAST   = 1 << 2,
	GO_NORTH  = 1 << 3,
	GO_DOWN   = 1 << 4,
	GO_UP     = 1 << 5,
	GO_REPEAT = 1 << 6,
	GO_FORBID_REPEAT = 1 << 7
};

static enum movement backtrace_to_movement(enum backtrace bt)
{
	switch (bt) {
	case BT_WEST:
		return GO_EAST;
	case BT_SOUTH:
		return GO_NORTH;
	case BT_EAST:
		return GO_WEST;
	case BT_NORTH:
		return GO_SOUTH;
	case BT_DOWN:
		return GO_UP;
	case BT_UP:
		return GO_DOWN;
	default:
		return GO_NONE;
	}
}

static int movement_cardinal(enum movement m)
{
	return (m & (BT_WEST | BT_SOUTH | BT_EAST | BT_NORTH));
}

static int movement_vertical(enum movement m)
{
	return (m & (BT_UP | BT_DOWN));
}

#define STONE 1
#define REDSTONE_TORCH 76
#define TORCH_UP 5
#define PISTON_DOWN 0
#define REDSTONE_DUST 55
#define PLANKS 5
#define DIRT 3
#define UNLIT_REDSTONE_REPEATER 94
#define STICKY_PISTON 29
#define REDSTONE_BLOCK 152
#define SLIME 165

static void place_movement(struct extracted_net *en, struct coordinate c, enum movement m)
{
	if (movement_cardinal(m)) {
		struct coordinate b = {c.y - 1, c.z, c.x};
		switch (c.y) {
		case 3:
			extracted_net_append(en, b, PLANKS, 0);
			break;
		case 6:
			extracted_net_append(en, b, DIRT, 0);
			break;
		default:
			break;
		}

		if (m & GO_REPEAT)
			extracted_net_append(en, c, UNLIT_REDSTONE_REPEATER, 0);
		else
			extracted_net_append(en, c, REDSTONE_DUST, 0);

	} else if (m & GO_UP) {
		extracted_net_append(en, c, STONE, 0); c.y++;
		extracted_net_append(en, c, REDSTONE_TORCH, TORCH_UP); c.y++;
		extracted_net_append(en, c, STONE, 0); c.y++;
		extracted_net_append(en, c, REDSTONE_TORCH, TORCH_UP);

	} else if (m & GO_DOWN) {
		extracted_net_append(en, c, STICKY_PISTON, PISTON_DOWN); c.y--;
		extracted_net_append(en, c, SLIME, 0); c.y--;
		extracted_net_append(en, c, REDSTONE_BLOCK, 0); c.y--;
		c.y--;
		extracted_net_append(en, c, STONE, 0);
	}
}

static struct coordinate movement_displace(struct coordinate c, enum movement m)
{
	switch (m) {
	case GO_WEST:
		c.x--;
		break;
	case GO_SOUTH:
		c.z++;
		break;
	case GO_EAST:
		c.x++;
		break;
	case GO_NORTH:
		c.z++;
		break;
	case GO_DOWN:
		c.y -= 3;
		break;
	case GO_UP:
		c.y += 3;
		break;
	default:
		break;
	}
	return c;
}

// determines whether this segment is repeatable:
// it must be a cardinal movement going twice in the same direction
static int is_repeatable(enum movement *m, int i, int n)
{
	if (n < 2 || i > n-1)
		return 0;

	return movement_cardinal(m[i]) && m[i] == m[i+1] && !(m[i] & GO_FORBID_REPEAT);
}

#define MAX_REDSTONE_STRENGTH 13
#define MIN_REDSTONE_STRENGTH 2
// to a net that has the correct direction (that is, the source is seg.start
// and the sink is seg.end), place repeaters -- we want to place a minimum
// of repeaters, to reduce delay, but enough to propagate the signal
static void add_repeaters(enum movement *m, int n, int strength)
{

	if (n == 0)
		return;

	int i = 0;
	while (i < n) {
		// vertical movements involve torches and redstone blocks,
		// which reset the strength to 15
		if (movement_vertical(m[i])) {
			add_repeaters(&m[i+1], n-(i+1), MAX_REDSTONE_STRENGTH);
			return;
		}

		if (strength < MIN_REDSTONE_STRENGTH) {
			// begin moving backwards until we have a repeatable part
			for (int j = i; j >= 0; j--) {
				if (is_repeatable(m, j, n)) {
					m[j] |= GO_REPEAT;
					// restart adding repeaters from this point at the maximum strength
					add_repeaters(&m[j+1], n-(j+1), MAX_REDSTONE_STRENGTH);
					return;
				}
			}

			// for loop terminated, which means it couldn't find a repeatable part
			printf("[extract] could not repeat this segment\n");
			exit(1);
		}

		strength--;
		i++;
	}
}

static int has_abutting_relative_segment(struct routed_segment *rseg, struct coordinate c)
{
	if (rseg->parent && (adjacent(rseg->parent->seg.start, c) || adjacent(rseg->parent->seg.end, c)))
		return 1;

	for (int i = 0; i < rseg->n_child_segments; i++) {
		if (adjacent(rseg->child_segments[i]->seg.start, c) || adjacent(rseg->child_segments[i]->seg.end, c))
			return 1;
	}

	return 0;
}

static void extract_segment(struct extracted_net *en, struct routed_segment *rseg, int strength)
{
	// create list of movements
	int n_bt = rseg->n_backtraces;
	enum movement *movts = malloc(sizeof(enum movement) * rseg->n_backtraces);
	struct coordinate c = rseg->seg.start;

	// iterate through the movements, calculating the coordinates
	// as we'll need to make sure to forbid repeating where another
	// segment abuts with this one
	for (int i = 0; i < n_bt; i++) {
		movts[i] = backtrace_to_movement(rseg->bt[n_bt-1-i]);
		c = movement_displace(c, movts[i]);
		if (has_abutting_relative_segment(rseg, c))
			movts[i] |= GO_FORBID_REPEAT;
	}

	add_repeaters(movts, n_bt, strength);

	rseg->extracted = 1;

	c = rseg->seg.start;
	place_movement(en, c, GO_NONE);
	strength--;
	for (int i = 0; i < n_bt; i++) {
		c = movement_displace(c, movts[i]);
		place_movement(en, c, movts[i]);

		if (movement_cardinal(movts[i]))
			strength--;
		else if (movement_vertical(movts[i]) || movts[i] & GO_REPEAT)
			strength = MAX_REDSTONE_STRENGTH;
		assert(strength > 0);

		// for parent and child segments, extract them if they abut here
		if (rseg->parent && !rseg->parent->extracted) {
			struct segment seg = rseg->parent->seg;
			if (adjacent(seg.start, c) || adjacent(seg.end, c))
				extract_segment(en, rseg->parent, strength);
		}

		for (int i = 0; i < rseg->n_child_segments; i++) {
			struct segment seg = rseg->child_segments[i]->seg;
			if (adjacent(seg.start, c) || adjacent(seg.end, c))
				extract_segment(en, rseg->parent, strength);
		}
	}

	free(movts);
}

struct extracted_net *extract_net(struct routed_net *rn)
{
	struct extracted_net *en = malloc(sizeof(struct extracted_net));
	en->net = rn->net;
	en->n = 0;
	en->sz = 4;
	en->c = malloc(sizeof(struct coordinate) * en->sz);
	en->b = malloc(sizeof(block_t) * en->sz);
	en->d = malloc(sizeof(data_t) * en->sz);

	for (struct routed_segment_head *rsh = rn->routed_segments; rsh; rsh = rsh->next)
		rsh->rseg.extracted = 0;

	// find the driving pin
	struct routed_segment *ds = NULL;
	struct placed_pin *dp = NULL;
	for (int i = 0; i < rn->n_pins; i++) {
		if (rn->pins[i].cell_pin->direction == OUTPUT) {
			// we can only have one driver
			assert(!ds && !dp && rn->pins[i].parent);
			ds = rn->pins[i].parent;
			dp = &rn->pins[i];
		}
	}

	// now, find the part of the segment where
	assert(adjacent(dp->coordinate, ds->seg.start) || adjacent(dp->coordinate, ds->seg.end));
	if (adjacent(dp->coordinate, ds->seg.end))
		reverse_segment(ds);

	extract_segment(en, ds, dp->cell_pin->level);

	return en;
}

// block placement routines
void place_block(struct extraction *e, struct coordinate c, int xz_margin)
{
	c.z += xz_margin;
	c.x += xz_margin;

	struct dimensions d = e->dimensions;

	if (c.x > d.x || c.y > d.y || c.z > d.z || c.x < 0 || c.y < 0 || c.z < 0)
		return;

	e->blocks[c.y * d.z * d.x + c.z * d.x + c.x] = 55;

	switch (c.y) {
	case 3:
		e->blocks[(c.y - 1) * d.z * d.x + c.z * d.x + c.x] = 5;
		break;
	case 6:
		e->blocks[(c.y - 1) * d.z * d.x + c.z * d.x + c.x] = 3;
		break;
	default:
		break;
	}
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

	int margin = 2;

	// modify dimensions (we'll perform the displacement later)
	d.z = d.z - disp.z + 2 * margin;
	d.x = d.x - disp.x + 2 * margin;
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
					// index into extraction matrix
					int fd_off = (c.y - disp.y + y) * d.z * d.x + (c.z - disp.z + margin + z) * d.x + (c.x - disp.x + margin + x);
					int b_off = y * lcd.z * lcd.x + z * lcd.x + x; // block offset into logic cell

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
				// if it's vertical, place a block at its current position too
				if (is_vertical(rsh->rseg.bt[k]))
					place_block(e, coordinate_sub(c, disp), margin);

				c = disp_backtrace(c, rsh->rseg.bt[k]);
				struct coordinate cc = coordinate_sub(c, disp);

				place_block(e, cc, margin);
			}

			place_block(e, coordinate_sub(rsh->rseg.seg.start, disp), margin);
			place_block(e, coordinate_sub(rsh->rseg.seg.end, disp), margin);
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
