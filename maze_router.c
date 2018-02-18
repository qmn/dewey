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

	/* at most one of these can be non-null */
	struct routed_segment *origin_segment;
	struct placed_pin *origin_pin;
};

static struct routing_group *alloc_routing_group(unsigned int usage_size)
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

void free_routing_group(struct routing_group *rg)
{
	free_cost_coord_heap(rg->heap);
	free(rg->bt);
	free(rg->cost);
}

// routines to initialize a routing group based on a pin or a segment
static void init_routing_group_with_pin(struct routing_group *rg, struct usage_matrix *m, struct placed_pin *p, struct routing_group **visited)
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

static void init_routing_group_with_segment(struct routing_group *rg, struct usage_matrix *m, struct routed_segment *rseg, struct routing_group **visited)
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
static void routed_segment_add_child(struct routed_segment *rseg, struct routing_group *child)
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
static struct routing_group *routing_group_find(struct routing_group *rg)
{
	if (rg->parent != rg)
		rg->parent = routing_group_find(rg->parent);

	return rg->parent;
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

static int count_groups_to_merge(struct routed_net *rn)
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

struct coordinate disp_backtrace(struct coordinate c, enum backtrace b)
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

				// if these groups happened to be adjacent already, don't create a new segment;
				// instead, merge the adjacent group into this one and keep looking
				if (rg->bt[usage_idx(m, c)] == BT_START && visited_rg->bt[usage_idx(m, cc)] == BT_START) {
					visited_rg->parent = rg;
					visited_rg = rg;
					remaining_groups--;
				} else {
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

	// free things used in routing
	for (int i = 0; i < n_groups; i++)
		free_routing_group(rgs[i]);
	free(rgs);

	free(visited);

	// printf("[maze_reroute] n_routed_segments=%d, n_pins=%d\n", rn->n_routed_segments, rn->n_pins);
	// assert(rn->n_routed_segments >= rn->n_pins - 1);
	// printf("[maze_route] done\n");
}
