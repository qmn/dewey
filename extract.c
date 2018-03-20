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
/*
static int adjacent(struct coordinate c1, struct coordinate c2)
{
	int dy = abs(c1.y - c2.y);
	int dz = abs(c1.z - c2.z);
	int dx = abs(c1.x - c2.x);

	return !coordinate_equal(c1, c2) && (dx <= 1 && dz <= 1) && (dy == 0 || dy == 3);
}
*/

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

#define MV_VERTICAL_MASK (GO_UP | GO_DOWN)
#define MV_CARDINAL_MASK (GO_WEST | GO_SOUTH | GO_EAST | GO_NORTH)

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

// identity conversion, essentially
static enum movement backtrace_IS_movement(enum backtrace bt)
{
	switch (bt) {
	case BT_WEST: return GO_WEST;
	case BT_SOUTH: return GO_SOUTH;
	case BT_EAST: return GO_EAST;
	case BT_NORTH: return GO_NORTH;
	case BT_DOWN: return GO_DOWN;
	case BT_UP: return GO_UP;
	default: return GO_NONE;
	}
}

static int movement_cardinal(enum movement m)
{
	return (m & MV_CARDINAL_MASK);
}

static int movement_vertical(enum movement m)
{
	return (m & MV_VERTICAL_MASK);
}

static data_t repeater_data(enum movement m)
{
	switch (m & MV_CARDINAL_MASK) {
	case GO_EAST:
		return 3;
	case GO_NORTH:
		return 2;
	case GO_WEST:
		return 1;
	default:
	case GO_SOUTH:
		return 0;
	}
}

#define AIR 0
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

static block_t base_block(struct coordinate c)
{
	if (c.y < 2)
		return STONE;
	else if (c.y < 5)
		return PLANKS;
	else if (c.y < 8)
		return DIRT;
	else
		return STONE;
}

// based on the movement done from here, place blocks at `c`
static void place_movement(struct extracted_net *en, struct coordinate c, enum movement m, struct coordinate disp)
{
	c = coordinate_add(c, disp);

	if (movement_cardinal(m) || m == GO_NONE) {
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
			extracted_net_append(en, c, UNLIT_REDSTONE_REPEATER, repeater_data(m));
		else
			extracted_net_append(en, c, REDSTONE_DUST, 0);

	} else if (m & GO_UP) {
		extracted_net_append(en, c, base_block(c), 0); c.y++; // 0
		extracted_net_append(en, c, REDSTONE_TORCH, TORCH_UP); c.y++; // 1
		extracted_net_append(en, c, base_block(c), 0); c.y++; // 2
		extracted_net_append(en, c, REDSTONE_TORCH, TORCH_UP);; // 3

	} else if (m & GO_DOWN) {
		extracted_net_append(en, c, STICKY_PISTON, PISTON_DOWN); c.y--; // 2
		extracted_net_append(en, c, REDSTONE_BLOCK, 0); c.y--; // 1
		extracted_net_append(en, c, AIR, 0); // 0
		extracted_net_append(en, c, base_block(c), 0); // -1
	}
}

static struct coordinate movement_displace(struct coordinate c, enum movement m)
{
	switch (m & (MV_CARDINAL_MASK | MV_VERTICAL_MASK)) {
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
		c.z--;
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
	if (n < 3 || i < 1)
		return 0;

	return movement_cardinal(m[i]) && m[i] == m[i-1] && !(m[i] & GO_FORBID_REPEAT);
}

#define MAX_REDSTONE_STRENGTH 15
#define MIN_REDSTONE_STRENGTH 2
// to a net that has the correct direction (that is, the source is seg.start
// and the sink is seg.end), place repeaters -- we want to place a minimum
// of repeaters, to reduce delay, but enough to propagate the signal
static void add_repeaters(enum movement *m, int* strengths, int n)
{
	if (n == 0)
		return;

	for (int i = 1; i < n; i++) {
		strengths[i] = strengths[i-1] - 1;

		// vertical movements involve torches and redstone blocks,
		// which reset the strength to 15
		if (movement_vertical(m[i])) {
			strengths[i] = 16;
			add_repeaters(&m[i], &strengths[i], n-i);
			return;
		}

		if (strengths[i] < MIN_REDSTONE_STRENGTH) {
			// begin moving backwards until we have a repeatable part
			for (int j = i; j >= 0; j--) {
				if (is_repeatable(m, j, n)) {
					m[j] |= GO_REPEAT;
					strengths[j] = 16;
					// restart adding repeaters from this point at the maximum strength
					add_repeaters(&m[j], &strengths[j], n-j);
					return;
				}
			}

			// for loop terminated, which means it couldn't find a repeatable part
			printf("[extract] could not repeat this segment\n");
			// exit(1);
		}
	}
}

struct neighbors {
	int n_neighbors;
	struct routed_segment **neighbors;
};

static struct neighbors find_neighbors(struct routed_net *rn, struct placed_pin *p, struct routed_segment *rseg)
{
	assert(!!p ^ !!rseg);

	struct neighbors n = {0, NULL};

	for (struct routed_segment_adjacency *rsa = rn->adjacencies; rsa; rsa = rsa->next) {
		// get ordinary segments as neighbors
		if (rseg && rsa->parent == rseg && rsa->child_type == SEGMENT) {
			n.neighbors = realloc(n.neighbors, sizeof(struct routed_segment *) * ++n.n_neighbors);
			n.neighbors[n.n_neighbors-1] = rsa->child.rseg;
		} else if (rseg && rsa->child_type == SEGMENT && rsa->child.rseg == rseg) {
			n.neighbors = realloc(n.neighbors, sizeof(struct routed_segment *) * ++n.n_neighbors);
			n.neighbors[n.n_neighbors-1] = rsa->parent;
		} else if (p && rsa->child_type == PIN && rsa->child.pin == p) {
			n.neighbors = realloc(n.neighbors, sizeof(struct routed_segment *) * ++n.n_neighbors);
			n.neighbors[n.n_neighbors-1] = rsa->parent;
		}
	}

	return n;
}

static int has_abutting_neighbor(struct neighbors neighbors, struct coordinate at)
{
	for (int i = 0; i < neighbors.n_neighbors; i++) {
		struct routed_segment *rseg = neighbors.neighbors[i];
		int n_bt = rseg->n_backtraces;

		struct coordinate c = rseg->seg.end;
		for (int j = -1; j <= n_bt; j++) {
			if (j >= 0 && j < n_bt)
				c = disp_backtrace(c, rseg->bt[j]);
			else if (j == n_bt)
				c = rseg->seg.start;

			if (coordinate_equal(c, at))
				return 1;
		}
	}

	return 0;
}

static enum movement ordinal_to_movement(enum ordinal_direction od)
{
	switch (od) {
		case NORTH: return GO_NORTH;
		case EAST: return GO_EAST;
		case SOUTH: return GO_SOUTH;
		case WEST: return GO_WEST;
		default: return GO_NONE;
	}
}

static void print_extracted_net(struct routed_net *rn)
{
	struct coordinate dbr = {0, 0, 0}, dtl = {0, 0, 0};
	struct routed_segment_head *rsh;

	for (rsh = rn->routed_segments; rsh; rsh = rsh->next) {
		struct routed_segment rseg = rsh->rseg;
		struct coordinate c = rseg.seg.end;
		for (int k = 0; k < rseg.n_backtraces; k++) {
			c = disp_backtrace(c, rseg.bt[k]);
			dbr = coordinate_piecewise_max(dbr, c);
			dtl = coordinate_piecewise_min(dtl, c);
		}

		dbr = coordinate_piecewise_max(dbr, rseg.seg.start);
		dbr = coordinate_piecewise_max(dbr, rseg.seg.end);
		dtl = coordinate_piecewise_min(dtl, rseg.seg.start);
		dtl = coordinate_piecewise_min(dtl, rseg.seg.end);
	}

	for (int i = 0; i < rn->n_pins; i++) {
		struct coordinate c = extend_pin(&rn->pins[i]);
		dbr = coordinate_piecewise_max(dbr, c);
		dtl = coordinate_piecewise_min(dtl, c);
	}

	struct coordinate disp = coordinate_neg(dtl);
	int w = dbr.x - dtl.x + 1;
	int h = dbr.z - dtl.z + 1;

	char *a = calloc(w * h, sizeof(char));

	char curr = 'A';

	for (rsh = rn->routed_segments; rsh; rsh = rsh->next) {
		struct routed_segment rseg = rsh->rseg;
		printf("segment %c: (%d, %d, %d) -> (%d, %d, %d)\n", curr, PRINT_COORD(rseg.seg.start), PRINT_COORD(rseg.seg.end));
		struct coordinate c = rseg.seg.end;
		for (int k = 0; k < rseg.n_backtraces; k++) {
			c = disp_backtrace(c, rseg.bt[k]);
			struct coordinate cc = coordinate_add(c, disp);
			int idx = cc.z * w + cc.x;
			if (a[idx])
				printf("  intersects with segment %c at (%d, %d, %d)\n", a[idx], PRINT_COORD(c));
			a[idx] = curr;
		}
		curr++;
	}

	curr = '1';
	for (int i = 0; i < rn->n_pins; i++) {
		struct coordinate c = extend_pin(&rn->pins[i]);
		int idx = c.z * w + c.x;
		if (a[idx])
			printf("  intersects with segment %c at (%d, %d, %d)\n", a[idx], PRINT_COORD(c));
		a[idx] = curr;
		curr++;
	}

	printf("\nextracted net follows:\n");
	for (int z = -3; z <= h; z++) {
		for (int x = -3; x <= w; x++) {
			if (z == -1 && x >= 0)
				putchar('-');
			else if (z == -3 && x >= 0)
				putchar(x / 10 + '0');
			else if (z == -2 && x >= 0)
				putchar(x % 10 + '0');
			else if (x == -3 && z >= 0)
				putchar(z / 10 + '0');
			else if (x == -2 && z >= 0)
				putchar(z % 10 + '0');
			else if (x == -1 && z >= 0)
				putchar('|');
			else if (x == -1 && z == -1)
				putchar('+');
			else if (z == h && x >= 0)
				putchar('-');
			else if (x == w && z >= 0)
				putchar('|');
			else if (a[z * w + x])
				putchar(a[z * w + x]);
			else
				putchar(' ');
		}
		putchar('\n');
	}
}

// forward declaration
static void extract_segment(struct extracted_net *, struct routed_net *rn, struct routed_segment *, struct coordinate, int, struct coordinate);

static void extract_pin(struct extracted_net *en, struct routed_net *rn, struct placed_pin *p, int strength, struct coordinate disp)
{
	struct neighbors neighbors = find_neighbors(rn, p, NULL);

	p->extracted = 1;

	struct coordinate c = extend_pin(p);
	place_movement(en, c, ordinal_to_movement(p->cell_pin->facing), disp);

	for (int i = 0; i < neighbors.n_neighbors; i++) {
		struct routed_segment *rseg = neighbors.neighbors[i];
		if (!rseg->extracted)
			extract_segment(en, rn, rseg, c, strength - 1, disp);
	}
}

static void propagate_extraction(struct extracted_net *en, struct routed_net *rn, struct routed_segment *rseg, struct neighbors neighbors, struct coordinate c, int strength, struct coordinate disp)
{
	for (int k = 0; k < neighbors.n_neighbors; k++) {
		struct routed_segment *nrseg = neighbors.neighbors[k];

		if (!nrseg->extracted)
			extract_segment(en, rn, nrseg, c, strength, disp);
	}

	// extract pins
	for (struct routed_segment_adjacency *rsa = rn->adjacencies; rsa; rsa = rsa->next) {
		if (rsa->parent == rseg && rsa->child_type == PIN) {
			struct placed_pin *p = rsa->child.pin;
			if (!p->extracted)
				extract_pin(en, rn, p, strength, disp);
		}
	}
}

// produce the next signal strength by making this movement
static int weaken(int strength, enum movement m)
{
	if (movement_cardinal(m))
		return strength - 1;
	else if (movement_vertical(m) || m & GO_REPEAT)
		return 16;
	else
		printf("[weaken] wat\n");

	return strength;
}

// extracts a segment `rseg`, outwards from coordinate `from`, displacing into the extraction with disp
static void extract_segment(struct extracted_net *en, struct routed_net *rn, struct routed_segment *rseg, struct coordinate from, int initial_strength, struct coordinate disp)
{
	int n_bt = rseg->n_backtraces;
	if (!n_bt)
		return;

	printf("i am segment (%d, %d, %d) -> (%d, %d, %d); from = (%d, %d, %d)\n", PRINT_COORD(rseg->seg.start), PRINT_COORD(rseg->seg.end), PRINT_COORD(from));

	// find the position of `from` in this segment
	struct coordinate c = rseg->seg.end;
	int i;
	for (i = 0; i <= n_bt; i++) {
		// printf("  checking (%d, %d, %d)\n", PRINT_COORD(c));
		if (coordinate_equal(c, from)) {
			printf("match!\n");
			break;
		}

		if (i < n_bt)
			c = disp_backtrace(c, rseg->bt[i]);
		else
			assert(coordinate_equal(c, rseg->seg.start));
	}

	// if it's strictly greater than n_bt, it's not found, return
	if (i > n_bt)
		return;

	struct neighbors neighbors = find_neighbors(rn, NULL, rseg);
	printf("i am segment (%d, %d, %d) -> (%d, %d, %d) and i have %d neighbors\n", PRINT_COORD(rseg->seg.start), PRINT_COORD(rseg->seg.end), neighbors.n_neighbors);
	rseg->extracted = 1;

	// catch anything else that's here
	propagate_extraction(en, rn, rseg, neighbors, from, initial_strength, disp);

	// do the extraction from `from` to rseg->seg.end
	int back_count = i;
	enum movement *back_movts = malloc(sizeof(enum movement) * back_count);
	int *back_strengths = malloc(sizeof(int) * back_count);
	back_strengths[0] = initial_strength;

	c = from;
	for (int j = 0; j < i; j++) {
		back_movts[j] = backtrace_to_movement(rseg->bt[i-j-1]);
		if (has_abutting_neighbor(neighbors, c))
			back_movts[j] |= GO_FORBID_REPEAT;
		c = movement_displace(c, back_movts[j]);
	}

	add_repeaters(back_movts, back_strengths, i);

	c = from;
	for (int j = 0; j < i; j++) {
		place_movement(en, c, back_movts[j], disp);
		c = movement_displace(c, back_movts[j]);
		propagate_extraction(en, rn, rseg, neighbors, c, back_strengths[j], disp);
	}

	assert(coordinate_equal(c, rseg->seg.end));
	place_movement(en, c, GO_NONE, disp);
	propagate_extraction(en, rn, rseg, neighbors, rseg->seg.end, back_strengths[i-1]-1, disp);

	// do the extraction from `from` to rseg->seg.start
	int fwd_count = n_bt - i;
	enum movement *fwd_movts = malloc(sizeof(enum movement) * fwd_count);
	int *fwd_strengths = malloc(sizeof(int) * fwd_count);
	int is_end = coordinate_equal(from, rseg->seg.end);
	fwd_strengths[0] = is_end ? initial_strength : weaken(initial_strength, backtrace_to_movement(rseg->bt[i]));

	c = from;
	for (int j = 0; j < fwd_count; j++) {
		fwd_movts[j] = backtrace_IS_movement(rseg->bt[j+i]);
		if (has_abutting_neighbor(neighbors, c))
			fwd_movts[j] |= GO_FORBID_REPEAT;
		c = movement_displace(c, fwd_movts[j]);
	}

	add_repeaters(fwd_movts, fwd_strengths, fwd_count);

	c = from;
	for (int j = 0; j < fwd_count; j++) {
		if (j != 0) // don't place `from` again
			place_movement(en, c, fwd_movts[j], disp);
		c = movement_displace(c, fwd_movts[j]);
		propagate_extraction(en, rn, rseg, neighbors, c, fwd_strengths[j], disp);
	}

	assert(coordinate_equal(c, rseg->seg.start));
	place_movement(en, rseg->seg.start, GO_NONE, disp);
	propagate_extraction(en, rn, rseg, neighbors, rseg->seg.start, fwd_strengths[n_bt-1]-1, disp);
}

/*
// extracts a segment `rseg`, outwards from coordinate `from`, displacing into the extraction with disp
static void extract_segment(struct extracted_net *en, struct routed_net *rn, struct routed_segment *rseg, struct coordinate from, int initial_strength, struct coordinate disp)
{
	int n_bt = rseg->n_backtraces;
	if (!n_bt)
		return;

	// find the index of `from` in this segment
	struct coordinate c = rseg->seg.end;
	int i;
	for (i = -1; i <= n_bt; i++) {
		if (i >= 0 && i < n_bt)
			c = disp_backtrace(c, rseg->bt[i]);
		else if (i == n_bt)
			c = rseg->seg.start;

		if (coordinate_equal(c, from))
			break;
	}

	// if i > n_bt (strictly greater than), we did not find `from`
	// in this segment
	if (i > n_bt) {
		return;
	} else if (i == -1) {
		from = rseg->seg.end;
		i = 0;
	} else if (i == n_bt) {
		from = rseg->seg.start;
	}

	struct neighbors neighbors = find_neighbors(rn, NULL, rseg);

	rseg->extracted = 1;

	//  n_bt         i       i+1         0
	// start < < < from < from_p1 < < < end
	// start < < < from > from_p1 > > > end
        //

	int strength = initial_strength;

	// do the backwards extraction (from `from`+1 to rseg->seg.end)
	if (i > -1) 
	enum movement *bmovts = malloc(sizeof(enum movement) * i);
	int *bstrengths = malloc(sizeof(int) * i);
	bstrengths[0] = initial_strength - 1;
	// set `from_p1` to the element after `from`, unless it is the start
	struct coordinate from_p1 = coordinate_equal(from, rseg->seg.end) ? rseg->seg.end : disp_backtrace(from, invert_backtrace(rseg->bt[i])); // this is from+1

	c = from_p1;
	for (int j = 0; j < i; j++) {
		bmovts[j] = backtrace_to_movement(invert_backtrace(rseg->bt[i-1-j]));
		c = movement_displace(c, bmovts[j]);
		if (has_abutting_neighbor(neighbors, c))
			bmovts[j] |= GO_FORBID_REPEAT;
	}

	add_repeaters(bmovts, bstrengths, i);

	// do backwards extract, and check for neighbors
	for (int j = -1; j < i; j++) {
		if (j < 0) {
			c = from_p1;
			strength = bstrengths[0];
		} else {
			c = movement_displace(c, bmovts[j]);
			place_movement(en, c, bmovts[j], disp);
			strength = bstrengths[i];
		}

		propagate_extraction(en, rn, rseg, neighbors, c, strength, disp);
	}


	// do the forward extraction (from `from` to rseg->seg.start)
	int fi = n_bt - i;
	enum movement *fmovts = malloc(sizeof(enum movement) * fi);
	int *fstrengths = malloc(sizeof(int) * fi);
	fstrengths[0] = initial_strength;
	c = from;
	for (int j = 0; j < fi; j++) {
		fmovts[j] = backtrace_to_movement(rseg->bt[j + i]);
		c = movement_displace(c, fmovts[j]);
		if (has_abutting_neighbor(neighbors, c))
			fmovts[j] |= GO_FORBID_REPEAT;
	}

	add_repeaters(fmovts, fstrengths, fi);

	c = from;
	for (int j = 0; j < fi; j++) {
		place_movement(en, c, fmovts[j], disp);
		strength = fstrengths[j];

		propagate_extraction(en, rn, rseg, neighbors, c, strength, disp);

		c = movement_displace(c, fmovts[j]);
	}

	// get start, end
	propagate_extraction(en, rn, rseg, neighbors, rseg->seg.start, fstrengths[fi-1]-1, disp);
	propagate_extraction(en, rn, rseg, neighbors, rseg->seg.end, bstrengths[i-1]-1, disp);

}
*/

struct extracted_net *extract_net(struct routed_net *rn, struct coordinate disp)
{
	struct extracted_net *en = malloc(sizeof(struct extracted_net));
	en->net = rn->net;
	en->n = 0;
	en->sz = 4;
	en->c = malloc(sizeof(struct coordinate) * en->sz);
	en->b = malloc(sizeof(block_t) * en->sz);
	en->d = malloc(sizeof(data_t) * en->sz);

	// find the driving pin
	struct routed_segment *ds = NULL;
	struct placed_pin *dp = NULL;
	for (struct routed_segment_adjacency *rsa = rn->adjacencies; rsa; rsa = rsa->next) {
		if (rsa->child_type == PIN && rsa->child.pin->cell_pin->direction == OUTPUT) {
			assert(!ds && !dp);
			dp = rsa->child.pin;
			ds = rsa->parent;
		}
	}
	assert(dp && ds);

	// reset all pins', segments' extraction
	for (struct routed_segment_head *rsh = rn->routed_segments; rsh; rsh = rsh->next)
		rsh->rseg.extracted = 0;

	for (int i = 0; i < rn->n_pins; i++)
		rn->pins[i].extracted = 0;

	extract_pin(en, rn, dp, dp->cell_pin->level - 1, disp);

	// ensure all segments extracted
	int extracted = 0, total = 0;
	for (struct routed_segment_head *rsh = rn->routed_segments; rsh; rsh = rsh->next, total++) {
		if (rsh->rseg.extracted)
			extracted++;

		// assert(rsh->rseg.extracted);
	}
	printf("[extraction] %d/%d extracted\n", extracted, total);

	return en;
}

// block placement routines
void place_block(struct extraction *e, struct coordinate c, block_t b, data_t da)
{
	struct dimensions d = e->dimensions;

	if (c.x > d.x || c.y > d.y || c.z > d.z || c.x < 0 || c.y < 0 || c.z < 0)
		return;

	int offset = c.y * d.z * d.x + c.z * d.x + c.x;

	e->blocks[offset] = b;
	e->data[offset] = da;
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

	struct coordinate rt_disp = {-disp.y, -disp.z + margin, -disp.x + margin};
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		struct extracted_net *en = extract_net(&rt->routed_nets[i], rt_disp);
		print_rsa(&rt->routed_nets[i]);
		print_extracted_net(&rt->routed_nets[i]);
		putchar('\n');
		for (int i = 0; i < en->n; i++)
			place_block(e, en->c[i], en->b[i], en->d[i]);
	}

	return e;
}

void free_extraction(struct extraction *e)
{
	free(e->blocks);
	free(e->data);
	free(e);
}
