#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "base_router.h"
#include "maze_router.h"
#include "heap.h"
#include "usage_matrix.h"
#include "util.h"

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

	// the (parentless) pin or segment that forms
	// the start from which a Lee's algo wavefront
	// begins
	enum {NONE, SEGMENT, PIN} origin_type;
	union {
		void *p;
		struct routed_segment *rseg;
		struct placed_pin *pin;
	} origin;
};

// also confusingly abbreviated MRI
struct maze_route_instance {
	struct routed_net *rn;

	struct usage_matrix *m;
	struct routing_group **visited;

	struct routing_group **rgs;

	int n_groups;         // groups we currently have
	int remaining_groups; // groups remaining to combine
};

static struct routing_group *alloc_routing_group(struct maze_route_instance *mri)
{
	unsigned int usage_size = USAGE_SIZE(mri->m);

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

	rg->origin_type = NONE;
	rg->origin.p = NULL;

	mri->rgs = realloc(mri->rgs, sizeof(struct routing_group *) * ++mri->n_groups);
	mri->rgs[mri->n_groups-1] = rg;

	return rg;
}

void free_routing_group(struct routing_group *rg)
{
	free_cost_coord_heap(rg->heap);
	free(rg->bt);
	free(rg->cost);
}

// union-by-rank's find() method adapted to routing groups
static struct routing_group *routing_group_find(struct routing_group *rg)
{
	if (rg->parent != rg)
		rg->parent = routing_group_find(rg->parent);

	return rg->parent;
}

// routines to initialize a routing group based on a pin or a segment
static void init_routing_group_with_pin(struct maze_route_instance *mri, struct routing_group *rg, struct placed_pin *p)
{
	struct coordinate start = extend_pin(p);
	rg->bt[usage_idx(mri->m, start)] = BT_START;
	rg->cost[usage_idx(mri->m, start)] = 0;
	struct cost_coord start_cc = {0, start};
	cost_coord_heap_insert(rg->heap, start_cc);
	mri->visited[usage_idx(mri->m, extend_pin(p))] = rg;

	rg->origin_type = PIN;
	rg->origin.pin = p;
}

// initializes a routing group with a segment
static void init_routing_group_with_segment(struct maze_route_instance *mri, struct routing_group *rg, struct routed_segment *rseg)
{
	struct coordinate c = rseg->seg.end;
	for (int i = 0; i < rseg->n_backtraces; i++) {
		c = disp_backtrace(c, rseg->bt[i]);
		assert(in_usage_bounds(mri->m, c));

		struct cost_coord next = {0, c};
		cost_coord_heap_insert(rg->heap, next);
		rg->cost[usage_idx(mri->m, c)] = 0;
		rg->bt[usage_idx(mri->m, c)] = BT_START;
		mri->visited[usage_idx(mri->m, c)] = rg;
	}

	rg->origin_type = SEGMENT;
	rg->origin.rseg = rseg;
}

// each routing group has an originating child or segment
static void routed_segment_add_child(struct routed_net *rn, struct routed_segment *rseg, struct routing_group *child, struct coordinate at)
{
	if (child->origin_type == PIN)
		add_adjacent_pin(rn, rseg, child->origin.pin);
	else if (child->origin_type == SEGMENT)
		add_adjacent_segment(rn, rseg, child->origin.rseg, at);
}


/* find the smallest active (i.e., parent = self) routing group that has a heap
   with elements */
static struct routing_group *find_smallest_heap(struct routing_group **rgs, unsigned int n_groups)
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

// TODO: implement this again
#define ROUTER_PREFER_CONTINUE_IN_DIRECTION 0

#define INITIAL_BT_SIZE 4

static int append_backtrace(enum backtrace ent, struct routed_segment *rseg, int bt_size)
{
	assert(rseg->n_backtraces < bt_size);
	rseg->bt[rseg->n_backtraces++] = ent;
	if (rseg->n_backtraces >= bt_size) {
		bt_size *= 2;
		rseg->bt = realloc(rseg->bt, bt_size * sizeof(enum backtrace));
	}
	return bt_size;
}

// starting from a coordinate, build a backtrace array tracing from `c` to the
// first instance of BT_START. the array is necessarily backwards
// create a routed_segment based on two backtraces:
// from a, to a BT_START in a_bt, and from b, to a BT_START in b_bt.
// a and b should be adjacent.
static struct routed_segment make_segment_from_backtrace(struct usage_matrix *m,
		struct coordinate a, struct coordinate b,
		enum backtrace *a_bt, enum backtrace *b_bt)
{
	int bt_size = INITIAL_BT_SIZE;
	struct routed_segment rseg = {{{0, 0, 0}, {0, 0, 0}}, 0, NULL, 0, NULL, 0};
	rseg.bt = calloc(bt_size, sizeof(enum backtrace));

	enum backtrace b_to_a = compute_backtrace(b, a);

	enum backtrace ent;

	// create backtrace from B side (the end)
	while (b_bt[usage_idx(m, b)] != BT_START) {
		assert(in_usage_bounds(m, b));
		ent = b_bt[usage_idx(m, b)];
		assert(ent != BT_NONE);

		if (is_vertical(ent))
			mark_via_violation_zone(m, b);

		b = disp_backtrace(b, ent);
		assert(in_usage_bounds(m, b));

		if (is_vertical(ent))
			mark_via_violation_zone(m, b);

		bt_size = append_backtrace(ent, &rseg, bt_size);
	}

	// now, at BT_START, we are at the end of the B side
	rseg.seg.end = b;

	// invert the B backtrace
	invert_backtrace_sequence(rseg.bt, rseg.n_backtraces);

	// add backtrace bridging (original) B and A
	bt_size = append_backtrace(b_to_a, &rseg, bt_size);

	// create backtrace to the A side (the start)
	while (a_bt[usage_idx(m, a)] != BT_START) {
		assert(in_usage_bounds(m, a));
		ent = a_bt[usage_idx(m, a)];
		assert(ent != BT_NONE);

		if (is_vertical(ent))
			mark_via_violation_zone(m, a);

		a = disp_backtrace(a, ent);
		assert(in_usage_bounds(m, a));

		if (is_vertical(ent))
			mark_via_violation_zone(m, a);

		bt_size = append_backtrace(ent, &rseg, bt_size);
	}

	// now, at BT_START again, we are at the start of the A side
	rseg.seg.start = a;

	// todo: realloc() to resize rseg->bt to size
	return rseg;
}

int segment_in_bounds(struct usage_matrix *m, struct routed_segment *rseg)
{
	struct coordinate c = rseg->seg.end;
	for (int i = 0; i < rseg->n_backtraces; i++)
		if (!in_usage_bounds(m, (c = disp_backtrace(c, rseg->bt[i]))))
			return 0;

	return 1;
}

// if the parent group is not represented by an origin, then
// create a new routing group for it; otherwise, add it to an existing group
// to search for a pin, set rseg to NULL; for a segment, set p to NULL
// this is the bottom-up version of populate_routing_group()
// returns 1 if it made a new group
void find_or_make_routing_group(struct maze_route_instance *mri, struct placed_pin *p, struct routed_segment *rseg)
{
	assert(!!p ^ !!rseg);

	struct routed_segment *parent = find_parent(mri->rn, p, rseg);

	// find parent's routing group
	struct routing_group *rg = NULL;
	for (int i = 0; i < mri->n_groups; i++)
		if ((!parent && mri->rgs[i]->origin_type == PIN && mri->rgs[i]->origin.pin == p) ||
		    (parent && mri->rgs[i]->origin_type == SEGMENT && mri->rgs[i]->origin.rseg == parent))
			rg = mri->rgs[i];

	// if a parent isn't found among mri->rgs, this is a new routing group
	if (!rg) {
		rg = alloc_routing_group(mri);
	}

	// expand this routing group with the segment or pin
	if (parent) {
		if (rseg)
			init_routing_group_with_segment(mri, rg, rseg);
		else if (p)
			init_routing_group_with_pin(mri, rg, p);
		else
			printf("[find_or_make_routing_group] wat\n");

		// if we have a parent rseg, make sure we set the
		// parent as the origin
		rg->origin_type = SEGMENT;
		rg->origin.rseg = parent;
	} else {
		// if there's no parent segment, the pin is the origin
		init_routing_group_with_pin(mri, rg, p);
	}
}

// for all of the routing groups other than the one specified,
// find all children routing_groups and add their origin pins or
// segments into this one; behavior for setting origin is undefined
void populate_routing_group(struct maze_route_instance *mri, struct routing_group *rg)
{
	for (int i = 0; i < mri->n_groups; i++) {
		struct routing_group *org = mri->rgs[i], *porg;
		if (org == rg)
			continue;

		porg = routing_group_find(org);

		if (porg == rg) {
			if (org->origin_type == PIN)
				init_routing_group_with_pin(mri, rg, org->origin.pin);
			else if (org->origin_type == SEGMENT)
				init_routing_group_with_segment(mri, rg, org->origin.rseg);
			else
				printf("[populate_routing_group] wat\n");
		}
	}
}

struct maze_route_instance create_maze_route_instance(struct cell_placements *cp, struct routings *rt, struct routed_net *rn, int xz_margin)
{
	struct maze_route_instance mri = {NULL, NULL, NULL, NULL, 0, 0};

	mri.rn = rn;

	// create usage matrix
	mri.m = create_usage_matrix(cp, rt, xz_margin);
	for (int i = 0; i < rn->n_pins; i++)
		assert(in_usage_bounds(mri.m, rn->pins[i].coordinate));
	for (struct routed_segment_head *rsh = rn->routed_segments; rsh; rsh = rsh->next)
		assert(segment_in_bounds(mri.m, &rsh->rseg));

	unsigned int usage_size = USAGE_SIZE(mri.m);

	// track visiting routing_groups; NULL if not-yet visited
	mri.visited = calloc(usage_size, sizeof(struct routing_group *));

	// at fewest we can have just one remaining group
	mri.n_groups = 0;
	mri.rgs = NULL;

	// initialize pins
	for (int i = 0; i < rn->n_pins; i++)
		find_or_make_routing_group(&mri, &rn->pins[i], NULL);

	// intialize segments
	for (struct routed_segment_head *rsh = rn->routed_segments; rsh; rsh = rsh->next)
		find_or_make_routing_group(&mri, NULL, &rsh->rseg);

	mri.remaining_groups = mri.n_groups;

	return mri;
}

void free_mri(struct maze_route_instance mri)
{
	// free things used in routing
	for (int i = 0; i < mri.n_groups; i++)
		free_routing_group(mri.rgs[i]);
	free(mri.rgs);
	free(mri.visited);
}

// visit coordinate cc from coordinate c (of routing group rg), by using
// backtrace bt, merging groups as needed; if it merged, return 1, otherwise
// return 0
static int mri_visit(struct maze_route_instance *mri, struct routing_group *rg, struct coordinate c, enum movement mv)
{
	struct usage_matrix *m = mri->m;

	enum backtrace bt = movement_to_backtrace(mv);
	struct coordinate cc = disp_movement(c, mv);

	if (!in_usage_bounds(m, cc))
		return 0;

	// skip this if it's been marked BT_START
	if (rg->bt[usage_idx(m, cc)] == BT_START)
		return 0;

	// do not allow up/down movements from a BT_START
	if (rg->bt[usage_idx(m, c)] == BT_START && is_vertical(bt))
		return 0;

	int violation = 0;
	// if this is a vertical movement, make sure its origin and the
	// origin's backtrace are the same (for proper signal pointing)
	enum backtrace my_bt = rg->bt[usage_idx(m, c)];
	enum backtrace b4_bt = rg->bt[usage_idx(m, disp_backtrace(c, my_bt))]; // ha ha, "before"
	if (is_vertical(bt) && is_cardinal(my_bt) && my_bt != b4_bt)
		violation++;

	// if the coordinate (c) that led to the exploration of this coordinate
	// (cc) was itself explored by a vertical movement, make sure that this
	// movement (for c->cc) is the same as the one for the vertical
	// movement to this one
	if (is_vertical(b4_bt) && is_cardinal(my_bt) && bt != my_bt)
		violation++;

	// if we are adjacent to a via we did not just come from, it is a violation
	if (is_cardinal(bt)) {
		struct coordinate via_checks[4] = {{0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
		for (int i = 0; i < 4; i++) {
			struct coordinate ccc = coordinate_add(cc, via_checks[i]);
			if (in_usage_bounds(m, ccc) && is_vertical(rg->bt[usage_idx(m, ccc)]) && !coordinate_equal(ccc, c))
				violation++;
		}
	}

	// dissuade turns
	int turn_cost = (is_cardinal(bt) && is_cardinal(my_bt) && bt != my_bt) ? 5 : 0;
	int via_cost = is_vertical(bt) ? 20 : 0;
	int y_cost = c.y / 2;

	// dissuade going too close to bounds
	int edge_margin = 2;
	int edge_cost = (cc.x < edge_margin || cc.x > m->d.x - edge_margin || cc.z < edge_margin || cc.z > m->d.z - edge_margin) ? 4 : 0;

	int preferred_direction_cost = 0;
	if ((cc.y == 0 || cc.y == 6) && (bt == BT_NORTH || bt == BT_SOUTH))
		preferred_direction_cost = 10;
	else if (cc.y == 3 && (bt == BT_WEST || bt == BT_EAST))
		preferred_direction_cost = 10;

	int movement_cost = 1 + turn_cost + via_cost + y_cost + edge_cost + preferred_direction_cost;

	if (usage_matrix_violated(m, cc))
		violation++;

	int violation_cost = 1000;

	unsigned int cost_delta = movement_cost + (violation ? violation_cost : 0);
	unsigned int new_cost = rg->cost[usage_idx(m, c)] + cost_delta;

	// if the lowest min-heap element expands into another group that is
	// "independent" (i.e., it is its own parent)
	struct routing_group *visited_rg = mri->visited[usage_idx(m, cc)];
	if (visited_rg && visited_rg->parent == visited_rg && routing_group_find(visited_rg) != rg) {
		// we shouldn't expand vertically into another
		// net, but we won't overwrite the backtraces.
		if (is_vertical(bt))
			return 0;

		// if these groups happened to be adjacent already, don't create a new segment;
		// instead, merge the adjacent group into this one and keep looking
		if (rg->bt[usage_idx(m, c)] == BT_START && visited_rg->bt[usage_idx(m, cc)] == BT_START) {
			visited_rg->parent = rg;
			visited_rg = rg;
			mri->remaining_groups--;
		} else {
			// create a new segment arising from the merging of these two routing groups
			struct routed_segment_head *rsh = malloc(sizeof(struct routed_segment_head));
			rsh->next = NULL;
			rsh->rseg = make_segment_from_backtrace(m, c, cc, rg->bt, visited_rg->bt);
			rsh->rseg.net = mri->rn;
			routed_net_add_segment_node(mri->rn, rsh);

			// add, as children, the two groups formed by this segment
			struct routed_segment *rseg = &rsh->rseg;
			assert(rseg);
			routed_segment_add_child(mri->rn, rseg, rg, rseg->seg.start);
			routed_segment_add_child(mri->rn, rseg, visited_rg, rseg->seg.end);

			// create a new routing group based on this segment
			struct routing_group *new_rg = alloc_routing_group(mri);
			rg->parent = visited_rg->parent = new_rg;

			init_routing_group_with_segment(mri, new_rg, rseg);
			populate_routing_group(mri, new_rg);
			new_rg->origin_type = SEGMENT;
			new_rg->origin.rseg = rseg;

			// add the group to the list
			mri->remaining_groups--;
			return 1;
		}
	}

	// if we haven't already visited this one, add it to the heap
	// (and by "we" i mean this exact routing group, not its children)
	if (visited_rg != rg) {
		struct cost_coord next = {new_cost, cc};
		cost_coord_heap_insert(rg->heap, next);
	}

	// if this location has a lower score, update the cost and backtrace
	mri->visited[usage_idx(m, cc)] = rg;
	if (new_cost < rg->cost[usage_idx(m, cc)]) {
		rg->cost[usage_idx(m, cc)] = new_cost;
		rg->bt[usage_idx(m, cc)] = bt;
	}

	return 0;
}

static enum movement movts[] = {GO_WEST, GO_NORTH, GO_EAST, GO_SOUTH, GO_DOWN, GO_UP};
static int n_movts = sizeof(movts) / sizeof(enum movement);

// see silk.md for a description of this algorithm
// accepts a routed_net object, with any combination of previously-routed
// segments and unrouted pins and uses Lee's algorithm to connect them
// assumes that all routed_segments are contiguously placed
void maze_reroute(struct cell_placements *cp, struct routings *rt, struct routed_net *rn, int xz_margin)
{
	if (rn->n_pins <= 1)
		return;

	struct maze_route_instance mri = create_maze_route_instance(cp, rt, rn, xz_margin);

	// THERE CAN ONLY BE ONE-- i mean,
	// repeat until one group remains
	while (mri.remaining_groups > 1) {
		// select the smallest non-empty heap that is also its own parent (rg->parent = rg)
		struct routing_group *rg = find_smallest_heap(mri.rgs, mri.n_groups);
		assert(rg == rg->parent);

		// expand this smallest heap
		assert(rg->heap->n_elts > 0);
		struct coordinate c = cost_coord_heap_delete_min(rg->heap).coord;
		assert(in_usage_bounds(mri.m, c));

		// for each of the possible movements, expand in that direction
		for (int i = 0; i < n_movts; i++) {
			int merge_occurred = mri_visit(&mri, rg, c, movts[i]);
			if (merge_occurred)
				break;
		}
	}

	free_mri(mri);
	// printf("[maze_reroute] n_routed_segments=%d, n_pins=%d\n", rn->n_routed_segments, rn->n_pins);
	// assert(rn->n_routed_segments >= rn->n_pins - 1);
	// printf("[maze_route] done\n");
}
