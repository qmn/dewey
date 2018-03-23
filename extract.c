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

struct neighbor {
	struct coordinate at;
	enum {SEGMENT, PIN} tn;
	union {
		struct routed_segment *rseg;
		struct placed_pin *pin;
	} n;
};

struct neighbors {
	int n_neighbors;
	struct neighbor *neighbors;
};

static void add_neighbor_pin(struct neighbors *n, struct placed_pin *p, struct coordinate c)
{
	n->neighbors = realloc(n->neighbors, sizeof(struct neighbor) * ++n->n_neighbors);
	n->neighbors[n->n_neighbors-1].at = c;
	n->neighbors[n->n_neighbors-1].tn = PIN;
	n->neighbors[n->n_neighbors-1].n.pin = p;
}

static void add_neighbor_segment(struct neighbors *n, struct routed_segment *rseg, struct coordinate c)
{
	n->neighbors = realloc(n->neighbors, sizeof(struct neighbor) * ++n->n_neighbors);
	n->neighbors[n->n_neighbors-1].at = c;
	n->neighbors[n->n_neighbors-1].tn = SEGMENT;
	n->neighbors[n->n_neighbors-1].n.rseg = rseg;
}

static void check_adjacent(struct routed_net *rn, struct neighbors *n, struct coordinate c, void *skip)
{
	for (int i = 0; i < rn->n_pins; i++) {
		struct placed_pin *q = &rn->pins[i];
		if ((void *)q == skip)
			continue;

		if (coordinate_equal(c, extend_pin(q)))
			add_neighbor_pin(n, q, c);
	}

	for (struct routed_segment_head *rsh = rn->routed_segments; rsh; rsh = rsh->next) {
		struct routed_segment *rseg = &rsh->rseg;
		if ((void *)rseg == skip)
			continue;

		struct segment seg = rseg->seg;

		struct coordinate cc = seg.end;
		for (int i = 0; i < rseg->n_backtraces; i++) {
			if (coordinate_equal(c, cc))
				add_neighbor_segment(n, rseg, c);

			cc = disp_backtrace(cc, rseg->bt[i]);
		}

		assert(coordinate_equal(cc, seg.start));
		if (coordinate_equal(c, cc))
			add_neighbor_segment(n, rseg, c);
	}
}

static struct neighbors find_neighbors(struct routed_net *rn, struct placed_pin *p, struct routed_segment *rseg)
{
	assert(!!p ^ !!rseg);

	struct neighbors n = {0, NULL};

	if (p) {
		struct coordinate c = extend_pin(p);
		check_adjacent(rn, &n, c, (void *)p);

	} else if (rseg) {
		struct segment seg = rseg->seg;

		struct coordinate c = seg.end;
		for (int i = 0; i <= rseg->n_backtraces; i++) {
			check_adjacent(rn, &n, c, (void *)rseg);

			if (i < rseg->n_backtraces)
				c = disp_backtrace(c, rseg->bt[i]);
			else
				assert(coordinate_equal(c, seg.start));
		}
	}

	return n;
}

static int has_abutting_neighbor(struct neighbors neighbors, struct coordinate at)
{
	for (int i = 0; i < neighbors.n_neighbors; i++)
		if (coordinate_equal(neighbors.neighbors[i].at, at))
			return 1;

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

// forward declarations
static void propagate_extraction(struct extracted_net *, struct routed_net *, struct neighbors, struct coordinate, int, struct coordinate);
static void extract_segment(struct extracted_net *, struct routed_net *rn, struct routed_segment *, struct coordinate, int, struct coordinate);

static void extract_pin(struct extracted_net *en, struct routed_net *rn, struct placed_pin *p, int strength, struct coordinate disp)
{
	struct neighbors neighbors = find_neighbors(rn, p, NULL);

	p->extracted = 1;

	struct coordinate c = extend_pin(p);
	place_movement(en, c, ordinal_to_movement(p->cell_pin->facing), disp);

	propagate_extraction(en, rn, neighbors, c, strength, disp);
}

static void propagate_extraction(struct extracted_net *en, struct routed_net *rn, struct neighbors neighbors, struct coordinate c, int strength, struct coordinate disp)
{
	for (int k = 0; k < neighbors.n_neighbors; k++) {
		struct neighbor n = neighbors.neighbors[k];
		if (coordinate_equal(n.at, c)) {
			if (n.tn == PIN && !n.n.pin->extracted)
				extract_pin(en, rn, n.n.pin, strength, disp);
			else if (n.tn == SEGMENT && !n.n.rseg->extracted)
				extract_segment(en, rn, n.n.rseg, c, strength, disp);
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

	// printf("i am segment (%d, %d, %d) -> (%d, %d, %d); from = (%d, %d, %d)\n", PRINT_COORD(rseg->seg.start), PRINT_COORD(rseg->seg.end), PRINT_COORD(from));

	// find the position of `from` in this segment
	struct coordinate c = rseg->seg.end;
	int i;
	for (i = 0; i <= n_bt; i++) {
		// printf("  checking (%d, %d, %d)\n", PRINT_COORD(c));
		if (coordinate_equal(c, from)) {
			// printf("match!\n");
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
	// printf("i am segment (%d, %d, %d) -> (%d, %d, %d) and i have %d neighbors\n", PRINT_COORD(rseg->seg.start), PRINT_COORD(rseg->seg.end), neighbors.n_neighbors);
	rseg->extracted = 1;

	// catch anything else that's here
	propagate_extraction(en, rn, neighbors, from, initial_strength, disp);

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
		c = disp_movement(c, back_movts[j]);
	}

	add_repeaters(back_movts, back_strengths, i);

	c = from;
	for (int j = 0; j < i; j++) {
		if (j == 0 || !movement_vertical(back_movts[j-1]))
			place_movement(en, c, back_movts[j], disp);
		c = disp_movement(c, back_movts[j]);
		propagate_extraction(en, rn, neighbors, c, back_strengths[j], disp);
	}

	assert(coordinate_equal(c, rseg->seg.end));
	place_movement(en, c, GO_NONE, disp);
	propagate_extraction(en, rn, neighbors, rseg->seg.end, back_strengths[i-1]-1, disp);

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
		c = disp_movement(c, fwd_movts[j]);
	}

	add_repeaters(fwd_movts, fwd_strengths, fwd_count);

	c = from;
	for (int j = 0; j < fwd_count; j++) {
		if (j != 0 && !movement_vertical(fwd_movts[j-1])) // don't place `from` again
			place_movement(en, c, fwd_movts[j], disp);
		c = disp_movement(c, fwd_movts[j]);
		propagate_extraction(en, rn, neighbors, c, fwd_strengths[j], disp);
	}

	assert(coordinate_equal(c, rseg->seg.start));
	place_movement(en, rseg->seg.start, GO_NONE, disp);
	propagate_extraction(en, rn, neighbors, rseg->seg.start, fwd_strengths[n_bt-1]-1, disp);
}

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
		// print_rsa(&rt->routed_nets[i]);
		// print_extracted_net(&rt->routed_nets[i]);
		// putchar('\n');
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
