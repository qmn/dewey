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

void print_routings(struct routings *rt)
{
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++) {
			printf("[maze_route] net %d, segment %d: ", i, j);
			print_routed_segment(&rt->routed_nets[i].routed_segments[j]);
		}
	}
}

void free_routings(struct routings *rt)
{
	free(rt->routed_nets);
	free(rt);
}

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

		rn->n_routed_segments = on.n_routed_segments;
		rn->routed_segments = malloc(sizeof(struct routed_segment) * rn->n_routed_segments);
		
		/* for each routed_segment in routed_net */
		for (int j = 0; j < rn->n_routed_segments; j++) {
			struct routed_segment *rseg = &(rn->routed_segments[j]);
			struct routed_segment old_rseg = on.routed_segments[j];
			rseg->seg = old_rseg.seg;
			rseg->n_coords = old_rseg.n_coords;
			rseg->coords = malloc(sizeof(struct coordinate) * rseg->n_coords);
			memcpy(rseg->coords, old_rseg.coords, sizeof(struct coordinate) * rseg->n_coords);
			rseg->score = old_rseg.score;
			rseg->net = rn;
		}
	}

	return new_rt;
}

void routings_displace(struct routings *rt, struct coordinate disp)
{
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		struct routed_net *rn = &(rt->routed_nets[i]);
		for (int j = 0; j < rn->n_routed_segments; j++) {
			struct routed_segment *rseg = &(rn->routed_segments[j]);
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
		for (int j = 0; j < rn.n_routed_segments; j++) {
			struct routed_segment rseg = rn.routed_segments[j];
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

struct usage_matrix *create_usage_matrix(struct cell_placements *cp, struct routings *rt, int xz_margin)
{
	struct coordinate tlcp = placements_top_left_most_point(cp);
	struct coordinate tlrt = routings_top_left_most_point(rt);
	struct coordinate top_left_most = coordinate_piecewise_min(tlcp, tlrt);

	assert(top_left_most.x >= xz_margin && top_left_most.y >= 0 && top_left_most.z >= xz_margin);

	// size the usage matrix and allow for routing on y=0 and y=3
	struct dimensions d = dimensions_piecewise_max(compute_placement_dimensions(cp), compute_routings_dimensions(rt));
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
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++) {
			for (int k = 0; k < rt->routed_nets[i].routed_segments[j].n_coords; k++) {
				// here
				struct coordinate c = rt->routed_nets[i].routed_segments[j].coords[k];
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

struct routed_segment make_segment_from_backtrace(struct usage_matrix *m, enum backtrace *bt, struct coordinate start)
{
	int coords_size = 4;
	int n_coords = 0;
	struct coordinate *coords = calloc(coords_size, sizeof(struct coordinate));
	struct coordinate current = start;
	struct coordinate end = {0, 0, 0};
	coords[n_coords++] = current;

	// int prev_bt = BT_NONE;
	while (bt[usage_idx(m, current)] != BT_START) {

		enum backtrace bt_ent = bt[usage_idx(m, current)];
		assert(bt_ent != BT_NONE);

		// printf("[msfb] current = (%d, %d, %d), bt = %d\n", current.y, current.z, current.x, bt_ent);

		/* prefer moving in the same direction as opposed to the first thing the backtrace gives you */
		// struct coordinate cont = disp_backtrace(current, prev_bt);
		struct coordinate next = disp_backtrace(current, bt_ent);
/*
		if (prev_bt != BT_NONE && in_usage_bounds(m, cont) && cost[usage_idx(m, cont)] <= cost[usage_idx(m, next)]) {
			current = cont;
		} else {
*/
		if (is_vertical(bt_ent))
			mark_via_violation_zone(m, current);

		current = next;
		assert(in_usage_bounds(m, current));

		if (is_vertical(bt_ent))
			mark_via_violation_zone(m, current);
/*
			prev_bt = bt_ent;
		}
*/

		coords[n_coords++] = current;
		if (n_coords >= coords_size) {
			coords_size *= 2;
			coords = realloc(coords, coords_size * sizeof(struct coordinate));
		}
	}
	end = current;

	struct segment seg = {start, end};
	struct routed_segment rseg = {seg, n_coords, coords, 0, NULL};
/*
	printf("[msfb] adding segment (%d, %d, %d) -- (%d, %d, %d)\n  ", start.y, start.z, start.x, end.y, end.z, end.x);
	for (int i = 0; i < n_coords; i++)
		printf("(%d, %d, %d), ", coords[i].y, coords[i].z, coords[i].x);
	printf("\n");
*/

	return rseg;
}

/* MAZE REROUTE */

void maze_reroute(struct cell_placements *cp, struct routings *rt, struct routed_net *rn, int xz_margin)
{
	assert(rn->n_routed_segments == 0);
	assert(rn->routed_segments == NULL);
	assert(rn->n_pins > 1);

	struct usage_matrix *m = create_usage_matrix(cp, rt, xz_margin);
	for (int i = 0;  i < rn->n_pins; i++)
		assert(in_usage_bounds(m, rn->pins[i].coordinate));

	// there will be n_pins - 1 segments, each composed of 2 half-segments
	rn->n_routed_segments = 2 * (rn->n_pins - 1);
	rn->routed_segments = calloc(rn->n_routed_segments, sizeof(struct routed_segment));
	int segment_count = 0;

	unsigned int usage_size = m->d.x * m->d.y * m->d.z;

	// create group visited matrix, which is heap index + 1
	visitor_t *visited = calloc(usage_size, sizeof(visitor_t));

	// create heaps, backtrace matrices, and union-find structure
	unsigned int n_groups = rn->n_pins;
	struct cost_coord_heap **heaps = malloc(n_groups * sizeof(struct cost_coord_heap *));
	enum backtrace **bt = calloc(n_groups, sizeof(enum backtrace *));
	unsigned int **cost = calloc(n_groups, sizeof(unsigned int *));

	// initialize the heaps, costs, and the backtraces
	for (int i = 0; i < n_groups; i++) {
		heaps[i] = create_cost_coord_heap();
		bt[i] = calloc(usage_size, sizeof(enum backtrace));
		for (int j = 0; j < usage_size; j++)
			bt[i][j] = BT_NONE;
		cost[i] = malloc(usage_size * sizeof(unsigned int));
		memset(cost[i], 0xff, usage_size * sizeof(unsigned int));

		// extend the pin once and start there
		struct coordinate start = extend_pin(&rn->pins[i]);
		assert(in_usage_bounds(m, start));
		bt[i][usage_idx(m, start)] = BT_START;
		visited[usage_idx(m, start)] = to_visitor_t(i);
		cost[i][usage_idx(m, start)] = 0;
		struct cost_coord start_cc = {0, start};
		cost_coord_heap_insert(heaps[i], start_cc);
	}

	// make-set
	struct mst_ubr_node *groups = mst_make_set(n_groups);

	// repeat until one group remains
	for (int remaining_groups = n_groups; remaining_groups > 1; ) {

		// select the smallest heap with elements
		int smallest_heap_idx = -1;
		int smallest_heap_score = -1;
		int unset = 1;
		for (int i = 0; i < n_groups; i++) {
			int parent_heap = mst_find(&groups[i])->me;
			if (!heaps[parent_heap] || heaps[parent_heap]->n_elts == 0)
				continue;

			int heap_score = cost_coord_heap_peek(heaps[parent_heap]).cost;
			if (heap_score < smallest_heap_score || unset) {
				unset = 0;
				smallest_heap_score = heap_score;
				smallest_heap_idx = parent_heap;
			}
		}
		assert(!unset);

		// expand this smallest heap
		struct cost_coord_heap *smallest_heap = heaps[smallest_heap_idx];
		assert(smallest_heap->n_elts > 0);
		struct coordinate c = cost_coord_heap_delete_min(smallest_heap).coord;
		assert(in_usage_bounds(m, c));

		int did_merge = 0;

		for (int movt = 0; movt < n_movements; movt++) {
			struct coordinate cc = coordinate_add(c, movement_offsets[movt]);
			int movement_cost = is_vertical(backtraces[movt]) ? 10 : c.y == 3 ? 3 : 1;
			int violation_cost = 1000 + movement_cost;

			if (!in_usage_bounds(m, cc))
				continue;

			// skip this if it's been marked BT_START
			if (bt[smallest_heap_idx][usage_idx(m, cc)] == BT_START)
				continue;

			int violation = usage_matrix_violated(m, cc);
			unsigned int cost_delta = violation ? violation_cost : movement_cost;
			unsigned int new_cost = cost[smallest_heap_idx][usage_idx(m, c)] + cost_delta;

			// if the lowest min-heap element expands into a group
			// not currently in my union-find set...
			if (visited[usage_idx(m, cc)] != 0) {
				int their_heap_idx = from_visitor_t(visited[usage_idx(m, cc)]);
				int their_parent_group = mst_find(&groups[their_heap_idx])->me;

				if (mst_find(&groups[smallest_heap_idx]) != mst_find(&groups[their_heap_idx]) && !is_vertical(backtraces[movt])) {
					// create two new segments, meeting at this point
					// backtrace these points, filling in these segments
					struct coordinate meeting_point = cc;
					// printf("[maze_route] joining group %d and %d at (%d, %d, %d)\n",
					//	smallest_heap_idx, their_parent_group, cc.y, cc.z, cc.x);
					struct routed_segment backward_rseg = make_segment_from_backtrace(m, bt[smallest_heap_idx], c);
					struct routed_segment forward_rseg = make_segment_from_backtrace(m, bt[their_heap_idx], meeting_point);
					backward_rseg.net = rn;
					forward_rseg.net = rn;
					rn->routed_segments[segment_count++] = backward_rseg;
					rn->routed_segments[segment_count++] = forward_rseg;

					// union these groups
					mst_union(&groups[smallest_heap_idx], &groups[their_heap_idx]);

					// delete the min-heap corresponding to both sets
					assert(smallest_heap_idx != their_heap_idx);
					// printf("freeing cost coord heaps %d and %d\n", smallest_heap_idx, their_parent_group);
					free_cost_coord_heap(smallest_heap);
					smallest_heap = NULL;
					free_cost_coord_heap(heaps[their_parent_group]);
					heaps[smallest_heap_idx] = NULL;
					heaps[their_parent_group] = NULL;
					memset(cost[smallest_heap_idx], 0xff, usage_size * sizeof(unsigned int));
					memset(cost[their_parent_group], 0xff, usage_size * sizeof(unsigned int));
					for (int i = 0; i < usage_size; i++) {
						bt[smallest_heap_idx][i] = BT_NONE;
						bt[their_parent_group][i] = BT_NONE;
					}

					// create a new heap
					int my_new_group_idx = mst_find(&groups[smallest_heap_idx])->me;
					assert(my_new_group_idx == smallest_heap_idx || my_new_group_idx == their_parent_group);
					struct cost_coord_heap *new_heap = create_cost_coord_heap();
					heaps[my_new_group_idx] = new_heap;

					// zero out all points belonging to these groups in the visited matrix
					unmark_group_in_visited(m, visited, groups, my_new_group_idx, n_groups);

					// mark all points in these segments as visited
					// add these points to the new min-heap
					for (int i = 0; i < backward_rseg.n_coords; i++) {
						struct coordinate ccc = backward_rseg.coords[i];
						assert(in_usage_bounds(m, ccc));
						struct cost_coord next = {0, ccc};
						cost_coord_heap_insert(new_heap, next);
						bt[my_new_group_idx][usage_idx(m, ccc)] = BT_START;
						visited[usage_idx(m, ccc)] = to_visitor_t(my_new_group_idx);
					}

					for (int i = 0; i < forward_rseg.n_coords; i++) {
						struct coordinate ccc = forward_rseg.coords[i];
						assert(in_usage_bounds(m, ccc));
						struct cost_coord next = {0, ccc};
						cost_coord_heap_insert(new_heap, next);
						bt[my_new_group_idx][usage_idx(m, ccc)] = BT_START;
						visited[usage_idx(m, ccc)] = to_visitor_t(my_new_group_idx);
					}

					remaining_groups--;
					did_merge = 1;
					break;
				}
			} else {
				visited[usage_idx(m, cc)] = to_visitor_t(smallest_heap_idx);
				struct cost_coord next = {new_cost, cc};
				cost_coord_heap_insert(smallest_heap, next);
			}

			// also break out of this loop if we merged heaps
			if (did_merge)
				break;

			// if this location has a lower score, update the cost and backtrace
			if (new_cost < cost[smallest_heap_idx][usage_idx(m, cc)]) {
				cost[smallest_heap_idx][usage_idx(m, cc)] = new_cost;
				bt[smallest_heap_idx][usage_idx(m, cc)] = backtraces[movt];
			}

		}
	}
	assert(segment_count == rn->n_routed_segments);
	// printf("[maze_route] done\n");
}

/* net scoring routines */

static int max_net_score = -1;
static int min_net_score = -1;
static int total_nets = 0;

static int count_routings_violations(struct cell_placements *cp, struct routings *rt)
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
		total_nets += rnet->n_routed_segments;
		int score = 0;

		for (int j = 0; j < rnet->n_routed_segments; j++) {
			int segment_violations = 0;

			struct routed_segment *rseg = &(rnet->routed_segments[j]);

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
						// printf("[violation] by net %d, seg %d at (%d, %d, %d) with (%d, %d, %d)\n", i, j, c.y, c.z, c.x, cc.y, cc.z, cc.x);
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
			// printf("[crv] net %d seg %d score = %d\n", i, j, segment_score);
		}

		/* second loop actually marks segment in matrix */
		for (int j = 0; j < rnet->n_routed_segments; j++) {
			struct routed_segment *rseg = &(rnet->routed_segments[j]);

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

/* rip-up and natural selection routines */
struct rip_up_set {
	int n_ripped;
	struct routed_net **rip_up;
};

void rip_up(struct routed_net *rn)
{
	for (int i = 0; i < rn->n_routed_segments; i++) {
		struct routed_segment *rseg = &rn->routed_segments[i];
		rseg->n_coords = 0;
		free(rseg->coords);
		rseg->coords = NULL;
	}
	free(rn->routed_segments);
	rn->routed_segments = NULL;
	rn->n_routed_segments = 0;
}

int routed_net_score(struct routed_net *rn)
{
	int s = 0;
	for (int i = 0; i < rn->n_routed_segments; i++)
		s += rn->routed_segments[i].score;
	return s;
}

int routed_net_cmp(const void *a, const void *b)
{
	struct routed_net *aa = *(struct routed_net **)a;
	struct routed_net *bb = *(struct routed_net **)b;

	return routed_net_score(bb) - routed_net_score(aa);
}

static struct rip_up_set natural_selection(struct routings *rt)
{
	int rip_up_count = 0;
	int rip_up_size = 4;
	struct routed_net **rip_up = malloc(rip_up_size * sizeof(struct routed_net *));
	memset(rip_up, 0, rip_up_size * sizeof(struct routed_net *));

	int score_range = max_net_score - min_net_score;
	int bias = score_range / 8;
	int random_range = bias * 10;

	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		int score = 0;
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++)
			score += rt->routed_nets[i].routed_segments[j].score;

		int r = random() % random_range;
		int adjusted_score = score - min_net_score + bias;

		if (r < adjusted_score) {
#ifdef NATURAL_SELECTION_DEBUG
			printf("[natural_selection] ripping up net %d (rand(%d) = %d < %d)\n", i, random_range, r, adjusted_score);
#endif
			// print_routed_segment(&rt->routed_nets[i].routed_segments[j]);
			rip_up[rip_up_count++] = &rt->routed_nets[i];
			if (rip_up_count >= rip_up_size) {
				rip_up_size *= 2;
				rip_up = realloc(rip_up, rip_up_size * sizeof(struct routed_net *));
			}
		} else {
#ifdef NATURAL_SELECTION_DEBUG
			printf("[natural_selection] leaving net %d intact (rand(%d) = %d >= %d)\n", i, random_range, r, adjusted_score);
#endif
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

	struct routed_segment rseg = {seg, len, path, 0, NULL};
	return rseg;
}
/* generate the MST for this net to determine the order of connections,
 * then connect them all with a city*/
struct routed_net *dumb_route(struct blif *blif, struct net_pin_map *npm, net_t net)
{
	int n_pins = npm->n_pins_for_net[net];
	struct routed_net *rn = malloc(sizeof(struct routed_net));
	rn->net = net;

	rn->n_pins = n_pins;
	rn->pins = malloc(sizeof(struct placed_pin) * rn->n_pins);
	memcpy(rn->pins, npm->pins[net], sizeof(struct placed_pin) * rn->n_pins);

	rn->n_routed_segments = 0;

		
	if (n_pins == 1) {
		rn->n_routed_segments = 1;
		struct segment seg = {npm->pins[net][0].coordinate, npm->pins[net][0].coordinate};
		struct coordinate *coords = malloc(sizeof(struct coordinate));
		coords[0] = npm->pins[net][0].coordinate;

		struct routed_segment rseg = {seg, 1, coords, 0, rn};
		
		rn->routed_segments = malloc(sizeof(struct routed_segment));
		rn->routed_segments[0] = rseg;
	} else if (n_pins == 2) {
		struct segment seg = {extend_pin(&npm->pins[net][0]), extend_pin(&npm->pins[net][1])};
		rn->n_routed_segments = 1;
		rn->routed_segments = malloc(sizeof(struct routed_segment));
		rn->routed_segments[0] = cityblock_route(seg);
		rn->routed_segments[0].net = rn;

	} else {
		struct coordinate *coords = calloc(n_pins, sizeof(struct coordinate));
		for (int j = 0; j < n_pins; j++)
			coords[j] = extend_pin(&npm->pins[net][j]);

		struct segments *mst = create_mst(coords, n_pins);
		rn->n_routed_segments = n_pins - 1;
		rn->routed_segments = malloc(sizeof(struct routed_segment) * rn->n_routed_segments);
		for (int j = 0; j < mst->n_segments; j++) {
			struct routed_segment rseg = cityblock_route(mst->segments[j]);
			rseg.net = rn;
			rn->routed_segments[j] = rseg;
		}
		free_segments(mst);
		free(coords);
	}

	return rn;
}

static struct routings *initial_route(struct blif *blif, struct net_pin_map *npm)
{
	struct routings *rt = malloc(sizeof(struct routings));
	rt->n_routed_nets = npm->n_nets;
	rt->routed_nets = calloc(rt->n_routed_nets + 1, sizeof(struct routed_net));
	rt->npm = npm;

	for (net_t i = 1; i < npm->n_nets + 1; i++)
		rt->routed_nets[i] = *dumb_route(blif, npm, i);

	return rt;
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
	int violations = count_routings_violations(cp, rt);

	interrupt_routing = 0;
	signal(SIGINT, router_sigint_handler);

	printf("\n");
	while ((violations = count_routings_violations(cp, rt)) > 0 && !interrupt_routing) {
		struct rip_up_set rus = natural_selection(rt);
		// sort elements by highest score
		qsort(rus.rip_up, rus.n_ripped, sizeof(struct routed_net *), routed_net_cmp);

		printf("\r[router] Iterations: %4d, Violations: %d, Segments to re-route: %d", iterations + 1, violations, rus.n_ripped);
		fflush(stdout);


		for (int i = 0; i < rus.n_ripped; i++) {
			// printf("[router] rerouting net %d\n", rus.rip_up[i]->net);
			rip_up(rus.rip_up[i]);
		}

		for (int i = 0; i < rus.n_ripped; i++) {
			recenter(cp, rt, 2);
			// printf("[router] Rerouting a segment in net %d with score %d\n", rus.rip_up[i]->net->net, rus.rip_up[i]->score);
			maze_reroute(cp, rt, rus.rip_up[i], 2);
			// print_routed_segment(rus.rip_up[i]);
		}
		free(rus.rip_up);

		recenter(cp, rt, 2);

		iterations++;
	}

	signal(SIGINT, SIG_DFL);

	printf("\n[router] Solution found! Optimizing...\n");

	for (net_t i = 1; i < rt->n_routed_nets; i++) {
		rip_up(&rt->routed_nets[i]);
		recenter(cp, rt, 2);
		maze_reroute(cp, rt, &rt->routed_nets[i], 2);
	}

	printf("[router] Routing complete!\n");
	print_routings(rt);

	free_pin_placements(pp);
	// free_net_pin_map(npm); // screws with extract in vis_png

	return rt;
}
