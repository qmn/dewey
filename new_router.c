#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "extract.h"
#include "segment.h"
#include "placer.h"
#include "router.h"
#include "heap.h"
#include "blif.h"

static int interrupt_routing = 0;

static void router_sigint_handler(int a)
{
	printf("Interrupt\n");
	interrupt_routing = 1;
}

void print_routed_segment(struct routed_segment *rseg)
{
	assert(rseg->n_coords >= 0);
	for (int i = 0; i < rseg->n_coords; i++)
		printf("(%d, %d, %d) ", rseg->coords[i].y, rseg->coords[i].z, rseg->coords[i].x);
	printf("\n");
}

void print_routed_net(struct routed_net *rn)
{
	int j = 0;
	for (struct routed_segment_head *rsh = rn->routed_segments; rsh; rsh = rsh->next, j++) {
		printf("[maze_route] net %d, segment %d: ", rn->net, j);
		print_routed_segment(&rsh->rseg);
	}
}

void print_routings(struct routings *rt)
{
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		print_routed_net(&rt->routed_nets[i]);
	}
}

void free_routings(struct routings *rt)
{
	free(rt->routed_nets);
	free(rt);
}

/* the pointer arithmetic to copy parent/child references works because
   they point exclusively into the routings structure and its sub-structures. */
struct routings *copy_routings(struct routings *old_rt)
{
	struct routings *new_rt = malloc(sizeof(struct routings));
	new_rt->n_routed_nets = old_rt->n_routed_nets;
	new_rt->routed_nets = malloc((new_rt->n_routed_nets + 1) * sizeof(struct routed_net));
	new_rt->npm = old_rt->npm;

	/* for each routed_net in routings */
	for (net_t i = 1; i < new_rt->n_routed_nets + 1; i++) {
		struct routed_net *rn = &(new_rt->routed_nets[i]);
		struct routed_net on = old_rt->routed_nets[i];
		rn->n_pins = on.n_pins;
		rn->pins = malloc(rn->n_pins * sizeof(struct placed_pin));
		memcpy(rn->pins, on.pins, sizeof(struct placed_pin) * rn->n_pins);
		rn->routed_segments = NULL;
/*
		for (int j = 0; j < rn->n_pins; j++)
			rn->pins[j].parent = on.pins[j].parent - on.routed_segments + rn->routed_segments;
*/

		abort(); // TODO: implement copying segments

/*
		// for each routed_segment in routed_net
		for (int j = 0; j < rn->n_routed_segments; j++) {
			struct routed_segment *rseg = &(rn->routed_segments[j]);
			struct routed_segment old_rseg = on.routed_segments[j];
			rseg->seg = old_rseg.seg;
			rseg->n_coords = old_rseg.n_coords;
			rseg->coords = malloc(sizeof(struct coordinate) * rseg->n_coords);
			memcpy(rseg->coords, old_rseg.coords, sizeof(struct coordinate) * rseg->n_coords);
			rseg->score = old_rseg.score;
			rseg->net = rn;
			rseg->n_child_segments = old_rseg.n_child_segments;
			rseg->n_child_pins = old_rseg.n_child_pins;

			// copy parent, child segments, and child pins
			rseg->parent = old_rseg.parent - on.routed_segments + rn->routed_segments;
			for (int k = 0; k < rseg->n_child_segments; k++)
				rseg->child_segments[k] = old_rseg.child_segments[k] - on.routed_segments + rn->routed_segments;
			for (int k = 0; k < rseg->n_child_pins; k++)
				rseg->child_pins[k] = old_rseg.child_pins[k] - on.pins + rn->pins;
			
		}
*/
	}

	return new_rt;
}

void routings_displace(struct routings *rt, struct coordinate disp)
{
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		struct routed_net *rn = &(rt->routed_nets[i]);
		for (struct routed_segment_head *rsh = rn->routed_segments; rsh; rsh = rsh->next) {
			struct routed_segment *rseg = &rsh->rseg;
			for (int k = 0; k < rseg->n_coords; k++) {
				rseg->coords[k] = coordinate_add(rseg->coords[k], disp);
			}
		}

		for (int j = 0; j < rn->n_pins; j++)
			rn->pins[j].coordinate = coordinate_add(rn->pins[j].coordinate, disp);

		/* displace pins separately, as segments refer to them possibly more than once */
		for (int j = 0; j < rt->npm->n_pins_for_net[i]; j++) {
			struct placed_pin *p = &(rt->npm->pins[i][j]);
			p->coordinate = coordinate_add(p->coordinate, disp);
		}
	}
}

struct dimensions compute_routings_dimensions(struct routings *rt)
{
	struct coordinate d = {0, 0, 0};
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		struct routed_net rn = rt->routed_nets[i];
		for (struct routed_segment_head *rsh = rn.routed_segments; rsh; rsh = rsh->next) {
			struct routed_segment rseg = rsh->rseg;
			for (int k = 0; k < rseg.n_coords; k++) {
				struct coordinate c = rseg.coords[k];
				d = coordinate_piecewise_max(d, c);
			}

			d = coordinate_piecewise_max(d, rseg.seg.start);
			d = coordinate_piecewise_max(d, rseg.seg.end);
		}
	}

	/* the dimension is the highest coordinate, plus 1 on each */
	struct dimensions dd = {d.y + 1, d.z + 1, d.x + 1};

	return dd;
}

static int max(int a, int b)
{
	return a > b ? a : b;
}

static int min(int a, int b)
{
	return a < b ? a : b;
}

/* coordinate manipulation routines */
enum backtrace {BT_NONE, BT_WEST, BT_SOUTH, BT_EAST, BT_NORTH, BT_DOWN, BT_UP, BT_START};
enum backtrace backtraces[] = {BT_WEST, BT_NORTH, BT_EAST, BT_SOUTH, BT_DOWN, BT_UP};
struct coordinate movement_offsets[] = {{0, 0, 1}, {0, 1, 0}, {0, 0, -1}, {0, -1, 0}, {3, 0, 0}, {-3, 0, 0}};
int n_movements = sizeof(movement_offsets) / sizeof(struct coordinate);
#define is_vertical(bt) (bt == BT_UP || bt == BT_DOWN)

static struct coordinate disp_backtrace(struct coordinate c, enum backtrace b)
{
	switch (b) {
	case BT_WEST:
		c.x--; break;
	case BT_EAST:
		c.x++; break;
	case BT_NORTH:
		c.z--; break;
	case BT_SOUTH:
		c.z++; break;
	case BT_UP:
		c.y += 3; break;
	case BT_DOWN:
		c.y -= 3; break;
	default:
		break;
	}

	return c;
}

/* usage matrix routines */
struct usage_matrix {
	struct dimensions d;
	int xz_margin;
	unsigned char *matrix;
};

static int in_usage_bounds(struct usage_matrix *m, struct coordinate c)
{
	return c.x >= 0 && c.x < m->d.x &&
	       c.y >= 0 && c.y < m->d.y &&
	       c.z >= 0 && c.z < m->d.z;
}

static int usage_idx(struct usage_matrix *m, struct coordinate c)
{
	return (c.y * m->d.z * m->d.x) + (c.z * m->d.x) + c.x;
}

static void usage_mark(struct usage_matrix *m, struct coordinate c)
{
	m->matrix[usage_idx(m, c)]++;
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
	d.y = max(d.y, 4);
	d.z += xz_margin;
	d.x += xz_margin;

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
	}

	/* routings */
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		for (struct routed_segment_head *rsh = rt->routed_nets[i].routed_segments; rsh; rsh = rsh->next) {
			for (int k = 0; k < rsh->rseg.n_coords; k++) {
				// here
				struct coordinate c = rsh->rseg.coords[k];
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

static struct coordinate check_offsets[] = {
	{0, 0, 0}, // here
	{-1, 0, 0}, // below
	{0, 1, 0}, // north
	{-1, 1, 0}, // below-north
	{0, 0, 1}, // east
	{-1, 0, 1}, // below-east
	{0, -1, 0}, // south
	{-1, -1, 0}, // below-south
	{0, 0, -1}, // west
	{-1, 0, -1} // below-west
};

static int usage_matrix_violated(struct usage_matrix *m, struct coordinate c)
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

/*
// ltl = lowest top left (smallest y, z, x); hbr = highest bottom right (largest y, z, x)
static int coordinate_within(struct coordinate llt, struct coordinate hbr, struct coordinate c)
{
	return c.y >= llt.y && c.y <= hbr.y &&
	       c.z >= llt.z && c.z <= hbr.z &&
	       c.x >= llt.x && c.x <= hbr.x;
}

// tests whether this coordinate joins with a point with a via
static int is_legal_net_join_point(struct routed_net *rn, struct coordinate c)
{
	// creates a 3 x 3 x 7 tall area centered on c
	int occupied[3 * 3 * 7];
	memset(occupied, 0, 3 * 3 * 7 * sizeof(int));

	struct coordinate disp = {c.y - 3, c.z - 1, c.x - 1};

	for (int i = 0; i < rn->n_routed_segments; i++) {
		struct routed_segment rseg = rn->routed_segments[i];
		for (int j = 0; j < rseg.n_coords; j++) {
			struct coordinate a = coordinate_sub(rseg.coords[j], disp);
			if (a.y >= 0 && a.y <= c.y + 3 && a.z >= 0 && a.z <= c.z + 1 && a.x >= 0 && <= c.x + 1) {
				int idx = a.y * (3 * 3) + a.z * 3 + a.x;
				assert(idx >= 0 && idx < 3 * 3 * 7);

				occupied[idx]++;
			}
		}
	}

	struct coordinate offsets[] = {{0, 0, 1}, {0, 0, -1}, {0, 1, 0}, {0, -1, 0}};
	int y_offsets[] = {3, 0, -3};

	for (int i = 0; i < sizeof(offsets) / sizeof(struct coordinate); i++) {
		struct coordinate cc = offsets[i];
		int intersections = 0;
		for (int j = 0; j < sizeof(y_offsets) / sizeof(int); j++) {
			cc.y = y_offsets[j];
		}
	}

	return 1;
}
*/

/* maze routing routines */
typedef unsigned int visitor_t;
visitor_t to_visitor_t(unsigned int i) {
	return i + 1;
}

unsigned int from_visitor_t(visitor_t i) {
	return i - 1;
}

/*
void unmark_group_in_visited(struct usage_matrix *m, visitor_t *visited, struct mst_ubr_node *groups, int g, int n_groups)
{
	int g_parent = mst_find(&groups[g])->me;
	for (int y = 0; y < m->d.y; y++) {
		for (int z = 0; z < m->d.z; z++) {
			for (int x = 0; x < m->d.x; x++) {
				struct coordinate c = {y, z, x};
				visitor_t v = visited[usage_idx(m, c)];
				if (!v)
					continue;

				unsigned int group_id = from_visitor_t(v);
				assert(group_id < n_groups);
				struct mst_ubr_node *v_g = &groups[group_id];
				if (mst_find(v_g)->me == g_parent)
					visited[usage_idx(m, c)] = 0;
			}
		}
	}
}
*/

// mark the usage matrix in a 3x3 zone centered on c to prevent subsequent routings
static void mark_via_violation_zone(struct usage_matrix *m, struct coordinate c)
{
	int y1 = max(c.y, 0), y2 = min(c.y, m->d.y - 1);
	int z1 = max(c.z - 1, 0), z2 = min(c.z + 1, m->d.z - 1);
	int x1 = max(c.x - 1, 0), x2 = min(c.x + 1, m->d.x - 1);

	for (int y = y1; y <= y2; y++) {
		for (int z = z1; z <= z2; z++) {
			for (int x = x1; x <= x2; x++) {
				struct coordinate cc = {y, z, x};

				usage_mark(m, cc);
			}
		}
	}
}

/* MAZE REROUTE */

/* a single routing group may consist of any number of pins or
   already-existing segments, and tracks the state of the
   wavefront in maze_reroute. extant pins/wires are marked
   BT_START in the backtrace structure, and cost 0 in the heap. */
struct routing_group {
	/* if parent points to this group, it's an independent routing group.
           when non-NULL, another routing group has subsumed this one. */
	struct routing_group *parent;

	/* heap of cost/coordinate pairs to visit */
	struct cost_coord_heap *heap;

	/* backtrace for this routing instance */
	enum backtrace *bt;

	/* cost matrix for this routing instance */
	unsigned int *cost;

	/* at most one of these can be non-null */
	struct routed_segment *origin_segment;
	struct placed_pin *origin_pin;
};

struct routing_group *alloc_routing_group(unsigned int usage_size)
{
	struct routing_group *rg = malloc(sizeof(struct routing_group));
	rg->parent = rg;
	rg->heap = NULL;
	rg->bt = NULL;
	rg->cost = NULL;

	rg->heap = create_cost_coord_heap();

	rg->bt = calloc(usage_size, sizeof(enum backtrace));
	for (int i = 0; i < usage_size; i++)
		rg->bt[i] = BT_NONE;

	rg->cost = malloc(usage_size * sizeof(unsigned int));
	memset(rg->cost, 0xff, usage_size * sizeof(unsigned int));

	rg->origin_pin = NULL;
	rg->origin_segment = NULL;

	return rg;
}

void init_routing_group_with_pin(struct routing_group *rg, struct usage_matrix *m, struct placed_pin *p, struct routing_group **visited)
{
	struct coordinate start = extend_pin(p);
	rg->bt[usage_idx(m, start)] = BT_START;
	rg->cost[usage_idx(m, start)] = 0;
	struct cost_coord start_cc = {0, start};
	cost_coord_heap_insert(rg->heap, start_cc);
	visited[usage_idx(m, extend_pin(p))] = rg;

	rg->origin_pin = p;
	rg->origin_segment = NULL;
}

void init_routing_group_with_segment(struct routing_group *rg, struct usage_matrix *m, struct routed_segment *rseg, struct routing_group **visited)
{
	for (int i = 0; i < rseg->n_coords; i++) {
		struct coordinate c = rseg->coords[i];
		assert(in_usage_bounds(m, c));

		struct cost_coord next = {0, c};
		cost_coord_heap_insert(rg->heap, next);
		rg->cost[usage_idx(m, c)] = 0;
		rg->bt[usage_idx(m, c)] = BT_START;
		visited[usage_idx(m, c)] = rg;
	}

	for (int i = 0; i < rseg->n_child_segments; i++) {
		assert(rseg != rseg->child_segments[i]);
		init_routing_group_with_segment(rg, m, rseg->child_segments[i], visited);
	}

	for (int i = 0; i < rseg->n_child_pins; i++)
		init_routing_group_with_pin(rg, m, rseg->child_pins[i], visited);

	rg->origin_segment = rseg;
	rg->origin_pin = NULL; // reset any origin_pin that may have been set
}

// each routing group has an originating child or segment
void routed_segment_add_child(struct routed_segment *rseg, struct routing_group *child)
{
	assert(!(child->origin_pin && child->origin_segment));
	if (child->origin_pin) {
		rseg->child_pins = realloc(rseg->child_pins, sizeof(struct placed_pin *) * ++rseg->n_child_pins);
		rseg->child_pins[rseg->n_child_pins - 1] = child->origin_pin;
		child->origin_pin->parent = rseg;

	} else if (child->origin_segment) {
		rseg->child_segments = realloc(rseg->child_segments, sizeof(struct routed_segment *) * ++rseg->n_child_segments);
		rseg->child_segments[rseg->n_child_segments - 1] = child->origin_segment;
		assert(rseg != child->origin_segment);
		child->origin_segment->parent = rseg;

	}
}

// union-by-rank's find() method adapted to routing groups
struct routing_group *routing_group_find(struct routing_group *rg)
{
	if (rg->parent != rg)
		rg->parent = routing_group_find(rg->parent);

	return rg->parent;
}

/* find the smallest active (i.e., parent = self) routing group that has a heap
   with elements */
struct routing_group *find_smallest_heap(struct routing_group **rgs, unsigned int n_groups)
{
	struct routing_group *smallest = NULL;
	int smallest_score = 0;

	for (int i = 0; i < n_groups; i++) {
		struct routing_group *rg = rgs[i];
		if (!rg || rg->parent != rg || rg->heap->n_elts == 0)
			continue;

		int heap_score = cost_coord_heap_peek(rg->heap).cost;
		if (!smallest || heap_score < smallest_score) {
			smallest_score = heap_score;
			smallest = rg;
		}
	}

	assert(smallest);

	return smallest;
}

#define USAGE_SIZE(m) (m->d.x * m->d.y * m->d.z)
#define ROUTER_PREFER_CONTINUE_IN_DIRECTION 0

static int extend_coords_from_backtrace(struct routed_segment *rseg, struct usage_matrix *m,
	struct coordinate c, enum backtrace *bt, int coords_size)
{
	assert(rseg);
	assert(rseg->coords);
	assert(rseg->n_coords < coords_size);

	struct coordinate curr = c;
	rseg->coords[rseg->n_coords++] = curr;
	if (rseg->n_coords >= coords_size) {
		coords_size *= 2;
		rseg->coords = realloc(rseg->coords, coords_size * sizeof(struct coordinate));
	}

#if ROUTER_PREFER_CONTINUE_IN_DIRECTION
	int prev_bt = BT_NONE;
#endif
	while (bt[usage_idx(m, curr)] != BT_START) {
		enum backtrace bt_ent = bt[usage_idx(m, curr)];
		assert(bt_ent != BT_NONE);

		struct coordinate next = disp_backtrace(curr, bt_ent);

#if ROUTER_PREFER_CONTINUE_IN_DIRECTION
		struct coordinate cont = disp_backtrace(curr, prev_bt);
		if (prev_bt != BT_NONE && in_usage_bounds(m, cont) && cost[usage_idx(m, cont)] <= usage_idx(m, next)) {
			curr = cont;
		else
			prev_bt = bt_ent;
#endif

		if (is_vertical(bt_ent))
			mark_via_violation_zone(m, curr);

		curr = next;
		assert(in_usage_bounds(m, curr));

		if (is_vertical(bt_ent))
			mark_via_violation_zone(m, curr);

		rseg->coords[rseg->n_coords++] = curr;
		if (rseg->n_coords >= coords_size) {
			coords_size *= 2;
			rseg->coords = realloc(rseg->coords, coords_size * sizeof(struct coordinate));
		}
	}

	return coords_size;
}

static void reverse_coords(struct coordinate *coords, int n_coords)
{
	struct coordinate tmp;
	for (int i = 0; i < n_coords / 2; i++) {
		tmp = coords[i];
		coords[i] = coords[n_coords - 1 - i];
		coords[n_coords - i - 1] = tmp;
	}
}

#define INITIAL_COORDS_SIZE 4
/* create a routed_segment based on two backtraces:
   from a, to a BT_START in a_bt, and from b, to a BT_START in b_bt.
   a and b should be adjacent. */
struct routed_segment make_segment_from_points(struct usage_matrix *m,
		struct coordinate a, struct coordinate b,
		enum backtrace *a_bt, enum backtrace *b_bt)
{
	int coords_size = INITIAL_COORDS_SIZE;
	struct routed_segment rseg = {{{0, 0, 0}, {0, 0, 0}}, 0, NULL, 0, NULL, NULL, 0, NULL, 0, NULL};
	rseg.coords = calloc(coords_size, sizeof(struct coordinate));

	// create points from A side, and reverse these
	coords_size = extend_coords_from_backtrace(&rseg, m, a, a_bt, coords_size);
	reverse_coords(rseg.coords, rseg.n_coords);
	rseg.seg.start = rseg.coords[0];

	// create points from B side, and do not reverse these
	extend_coords_from_backtrace(&rseg, m, b, b_bt, coords_size);
	rseg.seg.end = rseg.coords[rseg.n_coords - 1];

/*
	printf("[msfp] made segment: ");
	for (int i = 0; i < rseg.n_coords; i++) {
		struct coordinate c = rseg.coords[i];
		printf("(%d, %d, %d) ", c.y, c.z, c.x);
	}
*/

	return rseg;
}

void free_routing_group(struct routing_group *rg)
{
	free_cost_coord_heap(rg->heap);
	free(rg->bt);
	free(rg->cost);
}

int segment_routed(struct routed_segment *rseg)
{
	struct segment seg = rseg->seg;
	struct coordinate s = seg.start, e = seg.end;
	return (s.y | s.z | s.x | e.y | e.z | e.x);
}

/*
// find the first non-routed segment entry in *rn and add it there,
// returning the final place it returned it at, or extends the routed_segments
// array to accommodate it
struct routed_segment *add_segment(struct routed_net *rn, struct routed_segment rseg)
{
	rseg.net = rn;
	for (int i = 0; i < rn->sz_routed_segments; i++) {
		if (!segment_routed(&rn->routed_segments[i])) {
			rn->routed_segments[i] = rseg;
			rn->n_routed_segments++;
			return &(rn->routed_segments[i]);
		}
	}

	// a search for an empty entry failed, extend the array
	assert(rn->n_routed_segments == rn->sz_routed_segments);
	rn->routed_segments = realloc(rn->routed_segments, sizeof(struct routed_segment) * ++rn->sz_routed_segments);
	rn->routed_segments[rn->n_routed_segments++] = rseg;

	return &(rn->routed_segments[rn->n_routed_segments - 1]);
}
*/

int count_groups_to_merge(struct routed_net *rn)
{
	int count = 0;
	for (int i = 0; i < rn->n_pins; i++) {
		struct placed_pin *p = &rn->pins[i];
		if (!p->parent)
			count++;
	}

	for (struct routed_segment_head *rsh = rn->routed_segments; rsh; rsh = rsh->next) {
		struct routed_segment *rseg = &rsh->rseg;
		if (!rseg->parent && segment_routed(rseg))
			count++;
	}

	return count;
}

void routed_net_add_segment_node(struct routed_net *rn, struct routed_segment_head *rsh)
{
	assert(rsh->next == NULL);

	struct routed_segment_head *tail = rn->routed_segments;
	if (!tail) {
		rn->routed_segments = rsh;
		return;
	}

	while (tail->next)
		tail = tail->next;

	tail->next = rsh;
}

// see silk.md for a description of this algorithm
// accepts a routed_net object, with any combination of previously-routed
// segments and unrouted pins and uses Lee's algorithm to connect them
// assumes that all routed_segments are contiguously placed
void maze_reroute(struct cell_placements *cp, struct routings *rt, struct routed_net *rn, int xz_margin)
{
	if (rn->n_pins <= 1)
		return;

	// create usage matrix and ensure everything is in-bounds
	struct usage_matrix *m = create_usage_matrix(cp, rt, xz_margin);
	unsigned int usage_size = USAGE_SIZE(m);
	for (int i = 0; i < rn->n_pins; i++)
		assert(in_usage_bounds(m, rn->pins[i].coordinate));
	for (struct routed_segment_head *rsh = rn->routed_segments; rsh; rsh = rsh->next)
		for (int j = 0; j < rsh->rseg.n_coords; j++)
			assert(in_usage_bounds(m, rsh->rseg.coords[j]));

	// track visiting routing_groups; NULL if not-yet visited
	struct routing_group **visited = calloc(usage_size, sizeof(struct routing_group *));

	// at fewest we can have just one remaining group
	int n_groups = 0;
	int total_groups = count_groups_to_merge(rn) * 2 - 1;
	struct routing_group **rgs = calloc(total_groups, sizeof(struct routing_group *));

	// initialize parent-less pins
	for (int i = 0; i < rn->n_pins; i++) {
		struct placed_pin *p = &rn->pins[i];
		if (!p->parent) {
			struct routing_group *pin_rg = alloc_routing_group(usage_size);
			init_routing_group_with_pin(pin_rg, m, p, visited);
			rgs[n_groups++] = pin_rg;
			assert(!pin_rg->origin_segment);
		}
	}

	// intialize parent-less segments
	for (struct routed_segment_head *rsh = rn->routed_segments; rsh; rsh = rsh->next) {
		struct routed_segment *rseg = &rsh->rseg;
		if (!rseg->parent && segment_routed(rseg)) {
			struct routing_group *seg_rg = alloc_routing_group(usage_size);
			init_routing_group_with_segment(seg_rg, m, rseg, visited);
			rgs[n_groups++] = seg_rg;
			assert(!seg_rg->origin_pin);
		}
	}

	int remaining_groups = n_groups;

	// THERE CAN ONLY BE ONE-- i mean,
	// repeat until one group remains
	while (remaining_groups > 1) {
		// select the smallest non-empty heap that is also its own parent (rg->parent = rg)
		struct routing_group *rg = find_smallest_heap(rgs, total_groups);
		assert(rg == rg->parent);

		// expand this smallest heap
		assert(rg->heap->n_elts > 0);
		struct coordinate c = cost_coord_heap_delete_min(rg->heap).coord;
		assert(in_usage_bounds(m, c));

		// for each of the possible movements, expand in that direction
		for (int movt = 0; movt < n_movements; movt++) {
			struct coordinate cc = coordinate_add(c, movement_offsets[movt]);
			int movement_cost = is_vertical(backtraces[movt]) ? 10 : c.y == 3 ? 3 : 1;
			int violation_cost = 1000 + movement_cost;

			if (!in_usage_bounds(m, cc))
				continue;

			// skip this if it's been marked BT_START
			if (rg->bt[usage_idx(m, cc)] == BT_START)
				continue;

			int violation = usage_matrix_violated(m, cc);
			unsigned int cost_delta = violation ? violation_cost : movement_cost;
			unsigned int new_cost = rg->cost[usage_idx(m, c)] + cost_delta;

			// if the lowest min-heap element expands into another group that is "independent" (i.e., it is its own parent)
			struct routing_group *visited_rg = visited[usage_idx(m, cc)];
			if (visited_rg && visited_rg->parent == visited_rg && routing_group_find(visited_rg) != rg) {
				// we shouldn't expand vertically into another
				// net, but we won't overwrite the backtraces.
				if (is_vertical(backtraces[movt]))
					continue;

				// create a new segment arising from the merging of these two routing groups
				struct routed_segment_head *rsh = malloc(sizeof(struct routed_segment_head));
				rsh->next = NULL;
				rsh->rseg = make_segment_from_points(m, c, cc, rg->bt, visited_rg->bt);
				rsh->rseg.net = rn;
				routed_net_add_segment_node(rn, rsh);

				struct routed_segment *rseg = &rsh->rseg;
				assert(rseg);
				routed_segment_add_child(rseg, rg);
				routed_segment_add_child(rseg, visited_rg);

				// create a new routing group based on this segment
				struct routing_group *new_rg = alloc_routing_group(USAGE_SIZE(m));
				rg->parent = visited_rg->parent = new_rg;
				init_routing_group_with_segment(new_rg, m, rseg, visited);
				assert(!new_rg->origin_pin);

				rgs[n_groups++] = new_rg;

				assert(n_groups <= total_groups);

				remaining_groups--;
				break;
			}

			// if we haven't already visited this one, add it to the heap
			// (and by "we" i mean this exact routing group, not its children)
			if (visited_rg != rg) {
				struct cost_coord next = {new_cost, cc};
				cost_coord_heap_insert(rg->heap, next);
			}

			// if this location has a lower score, update the cost and backtrace
			visited[usage_idx(m, cc)] = rg;
			if (new_cost < rg->cost[usage_idx(m, cc)]) {
				rg->cost[usage_idx(m, cc)] = new_cost;
				rg->bt[usage_idx(m, cc)] = backtraces[movt];
			}

		}
	}
	// printf("[maze_reroute] n_routed_segments=%d, n_pins=%d\n", rn->n_routed_segments, rn->n_pins);
	// assert(rn->n_routed_segments >= rn->n_pins - 1);
	// printf("[maze_route] done\n");
}

/* net scoring routines */

static int max_net_score = -1;
static int min_net_score = -1;
static int total_nets = 0;

static int count_routings_violations(struct cell_placements *cp, struct routings *rt, FILE *log)
{
	total_nets = 0;
	max_net_score = min_net_score = -1;

	/* ensure no coordinate is < 0 */
	struct coordinate tlcp = placements_top_left_most_point(cp);
	struct coordinate tlrt = routings_top_left_most_point(rt);
	struct coordinate top_left_most = coordinate_piecewise_min(tlcp, tlrt);
/*
	printf("[count_routings_violations] tlcp x: %d, y: %d, z: %d\n", tlcp.x, tlcp.y, tlcp.z);
	printf("[count_routings_violations] tlrt x: %d, y: %d, z: %d\n", tlrt.x, tlrt.y, tlrt.z);
	printf("[count_routings_violations] top_left_most x: %d, y: %d, z: %d\n", top_left_most.x, top_left_most.y, top_left_most.z);
*/
	assert(top_left_most.x >= 0 && top_left_most.y >= 0 && top_left_most.z >= 0);

	struct dimensions d = dimensions_piecewise_max(compute_placement_dimensions(cp), compute_routings_dimensions(rt));

	int usage_size = d.y * d.z * d.x;
	unsigned char *matrix = malloc(usage_size * sizeof(unsigned char));
	memset(matrix, 0, usage_size * sizeof(unsigned char));

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

		for (int y = c.y; y < cell_y; y++) {
			for (int z = z1; z < z2; z++) {
				for (int x = x1; x < x2; x++) {
					int idx = y * d.z * d.x + z * d.x + x;

					matrix[idx]++;
				}
			}
		}
	}

	int total_violations = 0;

	/* segments */
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		struct routed_net *rnet = &(rt->routed_nets[i]);
		int score = 0;

		for (struct routed_segment_head *rsh = rnet->routed_segments; rsh; rsh = rsh->next, total_nets++) {
			int segment_violations = 0;

			struct routed_segment *rseg = &rsh->rseg;

			for (int k = 0; k < rseg->n_coords; k++) {
				struct coordinate c = rseg->coords[k];
				// printf("[crv] c = (%d, %d, %d)\n", c.y, c.z, c.x);

				int block_in_violation = 0;
				for (int m = 0; m < sizeof(check_offsets) / sizeof(struct coordinate); m++) {
					struct coordinate cc = coordinate_add(c, check_offsets[m]);
					// printf("[crv] cc = (%d, %d, %d)\n", cc.y, cc.z, cc.x);

					// ignore if checking out of bounds
					if (cc.y < 0 || cc.y >= d.y || cc.z < 0 || cc.z >= d.z || cc.x < 0 || cc.x >= d.x) {
						// printf("[crv] oob\n");
						continue;
					}

					// only ignore the start/end pins for first/last blocks on net
					if (k == 0 || k == rseg->n_coords - 1) {
						int skip = 0;
						for (int n = 0; n < rnet->n_pins; n++) {
							struct coordinate pin_cc = rnet->pins[n].coordinate;
							if (coordinate_equal(cc, pin_cc))
								skip++;
							pin_cc.y--;
							if (coordinate_equal(cc, pin_cc))
								skip++;
						}

						if (skip)
							continue;
					}

					int idx = (cc.y * d.z * d.x) + (cc.z * d.x) + cc.x;

					// do not mark or it will collide with itself
					if (matrix[idx]) {
						block_in_violation++;
						// printf("[crv] violation\n");
						fprintf(log, "[violation] by net %d, seg %p at (%d, %d, %d) with (%d, %d, %d)\n", i, (void *)rseg, c.y, c.z, c.x, cc.y, cc.z, cc.x);
					}
				}

				if (block_in_violation) {
					segment_violations++;
					total_violations++;
				}
			}

			// printf("[crv] segment_violations = %d\n", segment_violations);

			int segment_score = segment_violations * 1000 + rseg->n_coords;
			rseg->score = segment_score;
			score += segment_score;
			fprintf(log, "[crv] net %d seg %p score = %d\n", i, (void *)rseg, segment_score);
		}

		/* second loop actually marks segment in matrix */
		for (struct routed_segment_head *rsh = rnet->routed_segments; rsh; rsh = rsh->next) {
			struct routed_segment *rseg = &rsh->rseg;

			for (int k = 0; k < rseg->n_coords; k++) {
				struct coordinate c = rseg->coords[k];
				int idx = (c.y * d.z * d.x) + (c.z * d.x) + c.x;
				matrix[idx]++;

				if (c.y - 1 > 0) {
					idx = (c.y - 1) * d.z * d.x + c.z * d.x + c.x;
					matrix[idx]++;
				}
			}

		}

		if (max_net_score == -1) {
			max_net_score = min_net_score = score;
		} else {
			max_net_score = max(score, max_net_score);
			min_net_score = min(score, min_net_score);
		}
	}

	free(matrix);

	return total_violations;
}

void mark_routing_congestion(struct coordinate c, struct dimensions d, unsigned int *congestion, unsigned char *visited)
{
	int margin = 1;
	int z_start = max(c.z - margin, 0);
	int z_end = min(c.z + margin, d.z - 1);
	int x_start = max(c.x - margin, 0);
	int x_end = min(c.x + margin, d.x - 1);

	for (int z = z_start; z <= z_end; z++)
		for (int x = x_start; x <= x_end; x++)
			if (!visited[z * d.x + x]++)
				congestion[z * d.x + x]++;
}

/* create a table showing the congestion of a routing by showing where
   nets are currently being routed */
void print_routing_congestion(struct routings *rt)
{
	struct dimensions d = compute_routings_dimensions(rt);
	unsigned int *congestion = calloc(d.x * d.z, sizeof(unsigned int));

	// avoid marking a net over itself
	unsigned char *visited = calloc(d.x * d.z, sizeof(unsigned char));

	for (net_t i = 1; i < rt->n_routed_nets; i++)
		for (struct routed_segment_head *rsh = rt->routed_nets[i].routed_segments; rsh; rsh = rsh->next) {
			for (int k = 0; k < rsh->rseg.n_coords; k++)
				mark_routing_congestion(rsh->rseg.coords[k], d, congestion, visited);

			memset(visited, 0, sizeof(unsigned char) * d.x * d.z);
		}

	free(visited);

	printf("[routing_congestion] Routing congestion appears below:\n");
	printf("[routing_congestion] Z X ");
	for (int x = 0; x < d.x; x++)
		printf("%3d ", x);
	printf("\n");

	for (int z = 0; z < d.z; z++) {
		printf("[routing_congestion] %3d ", z);
		for (int x = 0; x < d.x; x++)
			printf("%3d ", congestion[z * d.x + x]);
		printf("\n");
	}
	printf("\n");

	free(congestion);
}

/* rip-up and natural selection routines */
struct rip_up_set {
	int n_ripped;
	struct routed_segment **rip_up;
};

struct routed_segment_head *remove_rsh(struct routed_segment *rseg)
{
	struct routed_net *rn = rseg->net;
	// find the previous element to delete this element
	struct routed_segment_head *node = NULL;

	if (!rn->routed_segments)
		return NULL;

	if (&(rn->routed_segments->rseg) == rseg) {
		node = rn->routed_segments;
		rn->routed_segments = node->next;
	} else {
		struct routed_segment_head *prev;
		for (prev = rn->routed_segments; prev; prev = prev->next) {
			if (&(prev->next->rseg) == rseg)
				break;
		}
		assert(prev);

		node = prev->next;
		prev->next = node->next;
	}

	node->next = NULL;
	return node;
}

// it's important to maintain the order of the routed segments in the routed
// net because the rip_up set struct relies on pointers
void rip_up_segment(struct routed_segment *rseg)
{
	// unlink segments' parent pointer here
	for (int i = 0; i < rseg->n_child_segments; i++) {
		struct routed_segment *child_rseg = rseg->child_segments[i];
		child_rseg->parent = NULL;
	}
	if (rseg->child_segments)
		free(rseg->child_segments);
	rseg->child_segments = NULL;
	rseg->n_child_segments = 0;

	// unlink pins' parent pointer here
	for (int i = 0; i < rseg->n_child_pins; i++) {
		struct placed_pin *child_pin = rseg->child_pins[i];
		child_pin->parent = NULL;
	}
	if (rseg->child_pins)
		free(rseg->child_pins);
	rseg->child_pins = NULL;
	rseg->n_child_pins = 0;

	// find and remove this segment from its parent
	struct routed_segment *parent_rseg = rseg->parent;
	if (parent_rseg) {
		for (int i = 0; i < parent_rseg->n_child_segments; i++) {
			if (parent_rseg->child_segments[i] == rseg) {
				if (i < parent_rseg->n_child_segments - 1)
					parent_rseg->child_segments[i] = parent_rseg->child_segments[parent_rseg->n_child_segments - 1];
				else
					parent_rseg->child_segments[i] = NULL;

				parent_rseg->n_child_segments--;
				break;
			}
			// this search should succeed
			assert(i < parent_rseg->n_child_segments - 1);
		}
	}

	assert(rseg->coords);
	free(rseg->coords);
	rseg->coords = NULL;
	rseg->n_coords = 0;
	struct segment zero = {{0, 0, 0}, {0, 0, 0}};
	rseg->seg = zero;
	rseg->score = 0;
}

// to sort in descending order, reverse the subtraction
int rseg_score_cmp(const void *a, const void *b)
{
	struct routed_segment *aa = *(struct routed_segment **)a;
	struct routed_segment *bb = *(struct routed_segment **)b;

	return bb->score - aa->score;
}

static struct rip_up_set natural_selection(struct routings *rt, FILE *log)
{
	int rip_up_count = 0;
	int rip_up_size = 4;
	struct routed_segment **rip_up = calloc(rip_up_size, sizeof(struct routed_segment *));

	int score_range = max_net_score - min_net_score;
	int bias = score_range / 8;
	int random_range = bias * 10;

	fprintf(log, "[natural_selection] adjusted_score = score - %d (min net score) + %d (bias)\n", min_net_score, bias);
	fprintf(log, "[natural_selection] net   seg                rip   rand(%5d)   adj. score\n", score_range);
	fprintf(log, "[natural_selection] ---   ----------------   ---   -----------   ----------\n");

	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		for (struct routed_segment_head *rsh = rt->routed_nets[i].routed_segments; rsh; rsh = rsh->next) {
			struct routed_segment *rseg = &rsh->rseg;
			if (!segment_routed(rseg))
				continue;

			int r = random() % random_range;
			int adjusted_score = rseg->score - min_net_score + bias;

			if (r < adjusted_score) {
				if (log)
					// fprintf(log, "[natural_selection] ripping up net %d (rand(%d) = %d < %d)\n", i, random_range, r, adjusted_score);
					fprintf(log, "[natural_selection] %3d   %p    X         %5d   %5d\n", i, (void *)rseg, r, adjusted_score);
#ifdef NATURAL_SELECTION_DEBUG
				printf("[natural_selection] ripping up net %2d, segment %p (rand(%d) = %d < %d)\n", i, rseg, random_range, r, adjusted_score);
#endif
				// print_routed_segment(&rt->routed_nets[i].routed_segments[j]);
				rip_up[rip_up_count++] = rseg;
				if (rip_up_count >= rip_up_size) {
					rip_up_size *= 2;
					rip_up = realloc(rip_up, rip_up_size * sizeof(struct routed_segment *));
				}
			} else {
#ifdef NATURAL_SELECTION_DEBUG
				printf("[natural_selection] leaving intact net %2d, segment %p (rand(%d) = %d >= %d)\n", i, rseg, random_range, r, adjusted_score);
#endif
				if (log)
					fprintf(log, "[natural_selection] %3d   %p               %5d   %5d\n", i, (void *)rseg, r, adjusted_score);
					// fprintf(log, "[natural_selection] leaving net %d intact (rand(%d) = %d >= %d)\n", i, random_range, r, adjusted_score);
			}
		}
	}

	struct rip_up_set rus = {rip_up_count, rip_up};

	return rus;
}

/* dumb_routing */
struct routed_segment cityblock_route(struct segment seg)
{
	struct coordinate a = seg.start;
	struct coordinate b = seg.end;

	int len = distance_cityblock(a, b) + 1; // to include block you start at
	int count = 0;
	struct coordinate *path = malloc(sizeof(struct coordinate) * len);
	path[count++] = a;

	assert(a.y == b.y);

	struct coordinate c = a;
	while (c.x != b.x || c.z != b.z) {
		if (c.x > b.x)
			c.x--;
		else if (c.x < b.x)
			c.x++;
		else if (c.z > b.z)
			c.z--;
		else if (c.z < b.z)
			c.z++;

		path[count++] = c;
	}

	assert(count == len);

	struct routed_segment rseg = {seg, len, path, 0, NULL, NULL, 0, NULL, 0, NULL};
	return rseg;
}

struct ubr_node {
	struct mst_node *x;
	struct mst_node *y;
	int score;
};

int ubr_node_cmp(const void *a, const void *b)
{
	struct ubr_node *aa = (struct ubr_node *)a, *bb = (struct ubr_node *)b;
	return aa->score - bb->score;
}

struct mst_node {
	struct mst_node *parent;
	struct placed_pin *pin;
	int rank;
};

struct mst_node *dumb_mst_find(struct mst_node *x)
{
	if (x->parent != x)
		x->parent = dumb_mst_find(x->parent);

	return x->parent;
}

void dumb_mst_union(struct mst_node *x, struct mst_node *y)
{
	struct mst_node *rx = dumb_mst_find(x);
	struct mst_node *ry = dumb_mst_find(y);

	if (rx == ry)
		return;

	rx->parent = ry;

	if (rx->rank == ry->rank)
		ry->rank++;
}

struct routed_segment *find_parent_rseg(struct placed_pin *p)
{
	struct routed_segment *rseg = p->parent;
	if (!rseg)
		return NULL;

	while (rseg->parent != NULL)
		rseg = rseg->parent;

	return rseg;
}

void add_child_pin(struct routed_segment *rseg, struct placed_pin *p)
{
	rseg->child_pins = realloc(rseg->child_pins, sizeof(struct placed_pin *) * ++rseg->n_child_pins);
	rseg->child_pins[rseg->n_child_pins - 1] = p;
	p->parent = rseg;
}

void add_child_segment(struct routed_segment *parent_rseg, struct routed_segment *rseg)
{
	parent_rseg->child_segments = realloc(parent_rseg->child_segments, sizeof(struct routed_segment *) * ++parent_rseg->n_child_segments);
	parent_rseg->child_segments[parent_rseg->n_child_segments - 1] = rseg;
	rseg->parent = parent_rseg;
}

// create routed segments for this net blindly (that is, without regard
// to other objects) using cityblock routing
void dumb_mst_route(struct routed_net *rn)
{
	// generate MST nodes
	struct mst_node *mst_nodes = calloc(rn->n_pins, sizeof(struct mst_node));
	for (int i = 0; i < rn->n_pins; i++)
		mst_nodes[i] = (struct mst_node){&mst_nodes[i], &rn->pins[i], 0};

	// generate union-by-rank nodes
	int n_ubr_nodes = rn->n_pins * (rn->n_pins - 1) / 2;
	struct ubr_node *ubr_nodes = calloc(n_ubr_nodes, sizeof(struct ubr_node));
	int k = 0;
	for (int i = 0; i < rn->n_pins; i++)
		for (int j = i + 1; j < rn->n_pins; j++)
			ubr_nodes[k++] = (struct ubr_node){&mst_nodes[i], &mst_nodes[j], distance_cityblock(extend_pin(mst_nodes[i].pin), extend_pin(mst_nodes[j].pin))};

	qsort(ubr_nodes, n_ubr_nodes, sizeof(struct ubr_node), ubr_node_cmp);

	// search through the list of union-by-rank nodes
	int count = 0;
	for (int i = 0; i < n_ubr_nodes && count < rn->n_pins - 1; i++) {
		struct ubr_node a = ubr_nodes[i];
		struct mst_node *x = a.x;
		struct mst_node *y = a.y;

		// if this node has a different parent than me, create
		// a new segment connecting these two parents.
		if (dumb_mst_find(x) != dumb_mst_find(y)) {
			struct segment seg = {extend_pin(x->pin), extend_pin(y->pin)};
			struct routed_segment_head *rsh = malloc(sizeof(struct routed_segment_head));
			rsh->next = NULL;
			rsh->rseg = cityblock_route(seg);
			rsh->rseg.net = rn;
			
			// if either object being joined has a parent segment already,
			// set that parent segment as a child of this newly-formed
			// segment

			struct routed_segment *xp_rseg = find_parent_rseg(x->pin);
			if (xp_rseg)
				add_child_segment(&rsh->rseg, xp_rseg);
			else
				add_child_pin(&rsh->rseg, x->pin);

			struct routed_segment *yp_rseg = find_parent_rseg(y->pin);
			if (yp_rseg)
				add_child_segment(&rsh->rseg, yp_rseg);
			else
				add_child_pin(&rsh->rseg, y->pin);

			routed_net_add_segment_node(rn, rsh);
			count++;
			dumb_mst_union(x, y);
		}
	}
}

/* generate the MST for this net to determine the order of connections,
 * then connect them all with a city*/
void dumb_route(struct routed_net *rn, struct blif *blif, struct net_pin_map *npm, net_t net)
{
	int n_pins = npm->n_pins_for_net[net];
	assert(n_pins > 0);

	rn->net = net;
	rn->routed_segments = NULL;

	rn->n_pins = n_pins;
	rn->pins = malloc(sizeof(struct placed_pin) * rn->n_pins);
	memcpy(rn->pins, npm->pins[net], sizeof(struct placed_pin) * rn->n_pins);

	if (n_pins == 1) {
		struct segment seg = {npm->pins[net][0].coordinate, npm->pins[net][0].coordinate};
		struct coordinate *coords = malloc(sizeof(struct coordinate));
		coords[0] = npm->pins[net][0].coordinate;

		struct placed_pin **child_pins = malloc(sizeof(struct placed_pin *));
		child_pins[0] = &npm->pins[net][0];

		struct routed_segment_head *rsh = malloc(sizeof(struct routed_segment_head));
		rsh->rseg = (struct routed_segment){seg, 1, coords, 0, rn, NULL, 0, NULL, 1, child_pins};
		rsh->next = NULL;
		child_pins[0]->parent = &rsh->rseg;
		
		rn->routed_segments = rsh;
	} else if (n_pins == 2) {
		struct segment seg = {extend_pin(&npm->pins[net][0]), extend_pin(&npm->pins[net][1])};

		struct placed_pin **child_pins = malloc(2 * sizeof(struct placed_pin *));
		child_pins[0] = &npm->pins[net][0];
		child_pins[1] = &npm->pins[net][1];

		struct routed_segment_head *rsh = malloc(sizeof(struct routed_segment_head));
		rsh->rseg = cityblock_route(seg);
		rsh->rseg.net = rn;
		rsh->rseg.n_child_pins = 2;
		rsh->rseg.child_pins = child_pins;
		child_pins[0]->parent = &rsh->rseg;
		child_pins[1]->parent = &rsh->rseg;
		rsh->next = NULL;

		rn->routed_segments = rsh;
	} else {
		dumb_mst_route(rn);
	}
}

static struct routings *initial_route(struct blif *blif, struct net_pin_map *npm)
{
	struct routings *rt = malloc(sizeof(struct routings));
	rt->n_routed_nets = npm->n_nets;
	rt->routed_nets = calloc(rt->n_routed_nets + 1, sizeof(struct routed_net));
	rt->npm = npm;

	for (net_t i = 1; i < npm->n_nets + 1; i++)
		dumb_route(&rt->routed_nets[i], blif, npm, i);

	return rt;
}

void assert_in_bounds(struct routed_net *rn)
{
	int arbitrary_max = 1000;
	for (struct routed_segment_head *rsh = rn->routed_segments; rsh; rsh = rsh->next) {
		struct routed_segment *rseg = &rsh->rseg;
		if (segment_routed(rseg)) {
			for (int j = 0; j < rseg->n_coords; j++) {
				struct coordinate c = rseg->coords[j];
				assert(c.y >= 0 && c.z >= 0 && c.x >= 0 && c.y < arbitrary_max && c.z < arbitrary_max && c.x < arbitrary_max);
			}
		}
	}
}

/* main route subroutine */
struct routings *route(struct blif *blif, struct cell_placements *cp)
{
	struct pin_placements *pp = placer_place_pins(cp);
	struct net_pin_map *npm = placer_create_net_pin_map(pp);

	struct routings *rt = initial_route(blif, npm);
	print_routings(rt);
	recenter(cp, rt, 2);

	int iterations = 0;
	int violations = 0;

	interrupt_routing = 0;
	signal(SIGINT, router_sigint_handler);
	FILE *log = fopen("new_router.log", "w");

	violations = count_routings_violations(cp, rt, log);

	printf("\n");
	while ((violations = count_routings_violations(cp, rt, log)) > 0 && !interrupt_routing) {
		struct rip_up_set rus = natural_selection(rt, log);
		// sort elements by highest score
		qsort(rus.rip_up, rus.n_ripped, sizeof(struct routed_segment *), rseg_score_cmp);

		struct routed_net **nets_ripped = calloc(rus.n_ripped, sizeof(struct routed_net *));

		printf("\r[router] Iterations: %4d, Violations: %d, Segments to re-route: %d", iterations + 1, violations, rus.n_ripped);
		fprintf(log, "\n[router] Iterations: %4d, Violations: %d, Segments to re-route: %d\n", iterations + 1, violations, rus.n_ripped);
		fflush(stdout);
		fflush(log);

		// rip up all segments in rip-up set
		for (int i = 0; i < rus.n_ripped; i++) {
			fprintf(log, "[router] Ripping up net %d, segment %p (score %d)\n", rus.rip_up[i]->net->net, (void *)rus.rip_up[i], rus.rip_up[i]->score);
			nets_ripped[i] = rus.rip_up[i]->net;
			rip_up_segment(rus.rip_up[i]);

			// remove segment from rt
			struct routed_segment_head *rsh = remove_rsh(rus.rip_up[i]);
			free(rsh);
		}

		// individually reroute all net instances that have had rip-ups occur
		for (int i = 0; i < rus.n_ripped; i++) {
			struct routed_net *net_to_reroute = nets_ripped[i];
			if (!net_to_reroute)
				continue;

			fprintf(log, "[router] Rerouting net %d\n", net_to_reroute->net);

			recenter(cp, rt, 2);
			maze_reroute(cp, rt, net_to_reroute, 2);

			// prevent subsequent reroutings of this net
			for (int j = i + 1; j < rus.n_ripped; j++)
				if (nets_ripped[j] == net_to_reroute)
					nets_ripped[j] = NULL;

			// printf("[maze_reroute] Rerouted net %d\n", net_to_reroute->net);
			// print_routed_net(net_to_reroute);
			assert_in_bounds(net_to_reroute);
		}
		free(rus.rip_up);
		rus.n_ripped = 0;

		recenter(cp, rt, 2);

		iterations++;
	}

	signal(SIGINT, SIG_DFL);

	printf("\n[router] Solution found! Optimizing...\n");
	fprintf(log, "\n[router] Solution found! Optimizing...\n");
	fclose(log);

	// rip up a net wholesale and reroute it
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		struct routed_segment_head *next = rt->routed_nets[i].routed_segments, *curr;
		do {
			curr = next;
			next = next->next;
			rip_up_segment(&curr->rseg);
			free(curr);
		} while (next);
		rt->routed_nets[i].routed_segments = NULL;

		recenter(cp, rt, 2);
		maze_reroute(cp, rt, &rt->routed_nets[i], 2);
	}

	printf("[router] Routing complete!\n");
	print_routings(rt);

	free_pin_placements(pp);
	// free_net_pin_map(npm); // screws with extract in vis_png

	return rt;
}
