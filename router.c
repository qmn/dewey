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

static int max(int a, int b)
{
	return a > b ? a : b;
}

static int min(int a, int b)
{
	return a < b ? a : b;
}

/*
static enum router_movements {
	EAST,
	NORTH,
	WEST,
	SOUTH,
	UP,
	DOWN
};
*/

struct coordinate extend_in_direction(enum ordinal_direction facing, struct coordinate c)
{
	switch (facing) {
		case EAST: c.x++; break;
		case WEST: c.x--; break;
		case NORTH: c.z--; break;
		case SOUTH: c.z++; break;
		default: break;
	}

	return c;
}

struct coordinate extend_pin(struct placed_pin *p)
{
	return extend_in_direction(p->cell_pin->facing, p->coordinate);
}

static int interrupt_routing = 0;

/* connect the two nets with a rectilinear path*/
struct routed_segment cityblock_route(struct segment seg)
{
	struct coordinate a = extend_pin(seg.start);
	struct coordinate b = extend_pin(seg.end);

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

void print_routed_segment(struct routed_segment *rseg)
{
	assert(rseg->n_coords >= 0);
	for (int i = 0; i < rseg->n_coords; i++)
		printf("(%d, %d, %d) ", rseg->coords[i].y, rseg->coords[i].z, rseg->coords[i].x);
	printf("\n");
}

void free_routed_net(struct routed_net *rn)
{
	free(rn->pins);
	for (int i = 0; i < rn->n_routed_segments; i++)
		free(rn->routed_segments[i].coords);
	free(rn);
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

	if (n_pins == 0) {
#ifdef DUMB_ROUTE_DEBUG
		printf("[dumb_route] net %d (%s) has 0 pins\n", net, get_net_name(blif, net));
#endif
	} else if (n_pins == 1) {
#ifdef DUMB_ROUTE_DEBUG
		printf("[dumb_route] net %d (%s) has 1 pin at (%d, %d, %d)\n", net, get_net_name(blif, net), c.y, c.z, c.x);
#endif
		struct segment seg = {&npm->pins[net][0], &npm->pins[net][0]};
		struct coordinate *coords = malloc(sizeof(struct coordinate));
		coords[0] = npm->pins[net][0].coordinate;

		struct routed_segment rseg = {seg, 1, coords, 0, rn};
		
		rn->routed_segments = malloc(sizeof(struct routed_segment));
		rn->routed_segments[0] = rseg;
	} else if (n_pins == 2) {
#ifdef DUMB_ROUTE_DEBUG
		printf("[dumb_route] net %d (%s):\n  ", net, get_net_name(blif, net));
#endif
		struct segment seg = {&npm->pins[net][0], &npm->pins[net][1]};
		rn->n_routed_segments = 1;
		rn->routed_segments = malloc(sizeof(struct routed_segment));
		rn->routed_segments[0] = cityblock_route(seg);
		rn->routed_segments[0].net = rn;
#ifdef DUMB_ROUTE_DEBUG
		print_routed_segment(rn);
#endif
	} else {
		struct segments *mst = create_mst(npm->pins[net], n_pins);
		rn->n_routed_segments = n_pins - 1;
		rn->routed_segments = malloc(sizeof(struct routed_segment) * rn->n_routed_segments);
#ifdef DUMB_ROUTE_DEBUG
		printf("[dumb_route] net %d (%s) has %d pins:\n", net, get_net_name(blif, net), n_pins);
#endif
		for (int j = 0; j < mst->n_segments; j++) {
			struct routed_segment rseg = cityblock_route(mst->segments[j]);
			rseg.net = rn;
#ifdef DUMB_ROUTE_DEBUG
			printf("  segment %d: ", j);
			print_routed_segment(rseg);
#endif
			rn->routed_segments[j] = rseg;
		}
		free_segments(mst);
	}

	return rn;
}

/* displace the entire design up, down, and right, by the amount given in displacement,
 * allocating new memory as needed, but displacement cannot be negative */
/*
struct extraction *extraction_displace(struct extraction *e, struct coordinate disp)
{
	assert(displacement.x >= 0 && displacement.y >= 0 && displacement.z >= 0)

	struct dimensions old_d = e->dimensions;
	struct dimensions new_d = {old_d.y + disp.y, old_d.z + disp.z, old_d.x + disp.x};

	block_t new_blocks = calloc(new_d.y * new_d.z * new_d.x, sizeof(block_t));
	data_t new_data = calloc(new_d.y * new_d.z * new_d.x, sizeof(data_t));

	for (int y = 0; y < old_d.y; y++) {
		for (int z = 0; z < old_d.z; z++) {
			int ny = y + disp.y;
			int nz = z + disy.z;
			int new_off = ny * new_d.z * new_d.x + nz * new_d.x + disp.x;
			int old_off = y * old_d.z * old_d.x + z * old_d.x;
			memcpy(&new_blocks[new_off], &e->blocks[old_off], old_d.x * sizeof(block_t));
			memcpy(&new_data[new_off], &e->data[old_off], old_d.x * sizeof(data_t));
		}
	}

	struct extraction *e2 = malloc(sizeof(struct extraction));

	e2->dimensions = new_d;
	e2->blocks = new_blocks;
	e2->data = new_data;

	return e2;
}
*/

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

			d = coordinate_piecewise_max(d, rseg.seg.start->coordinate);
			d = coordinate_piecewise_max(d, rseg.seg.end->coordinate);
		}
	}

	/* the dimension is the highest coordinate, plus 1 on each */
	struct dimensions dd = {d.y + 1, d.z + 1, d.x + 1};

	return dd;
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

		/* displace pins separately, as segments refer to them possibly more than once */
		for (int j = 0; j < rt->npm->n_pins_for_net[i]; j++) {
			struct placed_pin *p = &(rt->npm->pins[i][j]);
			p->coordinate = coordinate_add(p->coordinate, disp);
		}
	}
}

void free_routings(struct routings *rt)
{
	free(rt->routed_nets);
	free(rt);
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

/* creates a violation matrix;
 * use in_usage_bounds() and usage_idx() to check bounds and index
 * into the matrix, respectively.
 *
 * this requires that the top-left-most x and z points are at least
 * `xz_margin` out; use recenter() with a non-zero xz_margin to displace.
 */
static struct usage_matrix *create_usage_matrix(struct cell_placements *cp, struct routings *rt, int xz_margin)
{
	/* ensure no coordinate is < 0 */
	struct coordinate tlcp = placements_top_left_most_point(cp);
	struct coordinate tlrt = routings_top_left_most_point(rt);
	struct coordinate top_left_most = coordinate_piecewise_min(tlcp, tlrt);
/*
	printf("[create_usage_matrix] tlcp x: %d, y: %d, z: %d\n", tlcp.x, tlcp.y, tlcp.z);
	printf("[create_usage_matrix] tlrt x: %d, y: %d, z: %d\n", tlrt.x, tlrt.y, tlrt.z);
	printf("[create_usage_matrix] top_left_most x: %d, y: %d, z: %d\n", top_left_most.x, top_left_most.y, top_left_most.z);
*/
	assert(top_left_most.x >= xz_margin && top_left_most.y >= 0 && top_left_most.z >= xz_margin);

	struct dimensions d = dimensions_piecewise_max(compute_placement_dimensions(cp), compute_routings_dimensions(rt));
	// allow for routing on y=0 and y=3
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

		for (int y = c.y; y < cell_y; y++) {
			for (int z = z1; z < z2; z++) {
				for (int x = x1; x < x2; x++) {
					struct coordinate cc = {y, z, x};
					int idx = usage_idx(m, cc);

					// only update usage_matrix within the cell bounds
					// if (z >= c.z && z < cell_z && x >= c.x && x < cell_x)
					m->matrix[idx]++;

/*
					// only update violation_matrix within matrix bounds
					if (in_usage_bounds(m, cc))
						violation_matrix[idx]++;
*/
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
				m->matrix[usage_idx(m, c)]++;

				// below
				struct coordinate c2 = c;
				c2.y--;
				if (in_usage_bounds(m, c2))
					m->matrix[usage_idx(m, c2)]++;

/*
				// violation matrix calculation
				for (int m = 0; m < sizeof(check_offsets) / sizeof(struct coordinate); m++) {
					struct coordinate cc = coordinate_add(c, check_offsets[m]);
					if (in_usage_bounds(m, cc))
						violation_matrix[usage_idx(m, cc)]++;
				}
*/
			}
		}
	}

	return m;
}

static struct routed_segment *segment_is_actually_in_routing(struct routings *rt, struct routed_segment *rseg)
{
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++)
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++)
			if (rseg == &rt->routed_nets[i].routed_segments[j])
				return rseg;
	return NULL;
}

/*
static int cost_coord_cmp(const void *a, const void *b)
{
	// make heapsort sort from largest to smallest
	return -1 * (((struct cost_coord *)a)->cost - ((struct cost_coord *)b)->cost);
}
*/

static int usage_matrix_violated(struct usage_matrix *m, struct coordinate c, struct coordinate *ignore, int n_ignore)
{
	for (int j = 0; j < sizeof(check_offsets) / sizeof(struct coordinate); j++) {
		struct coordinate cc = coordinate_add(c, check_offsets[j]);
		if (!in_usage_bounds(m, cc))
			continue;

		// avoid the start or end pins or their locations immediately below it
		int skip = 0;
		for (int k = 0; k < n_ignore; k++) {
			if (coordinate_equal(ignore[k], cc))
				skip++;
		}

		if (skip)
			continue;

		if (m->matrix[usage_idx(m, cc)])
			return 1;
	}

	return 0;
}

static struct cost_coord_heap *heap = NULL;

// backtracing
enum backtrace {BT_NONE, BT_WEST, BT_SOUTH, BT_EAST, BT_NORTH, BT_DOWN, BT_UP, BT_START};
enum backtrace backtraces[] = {BT_WEST, BT_NORTH, BT_EAST, BT_SOUTH, BT_DOWN, BT_UP};
struct coordinate movement_offsets[] = {
	{0, 0, 1}, {0, 1, 0}, {0, 0, -1}, {0, -1, 0}, {3, 0, 0}, {-3, 0, 0}
};

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

static void maze_reroute_segment(struct cell_placements *cp, struct routings *rt, struct routed_segment *rseg)
{
	/* rip out this segment */
	rseg->n_coords = 0;
	if (rseg->coords)
		free(rseg->coords);
	rseg->coords = NULL;

	/* the additional xz_margin to each side for a possible route is 2 */
/*
	printf("[maze_route] attempting to maze-route from (%d, %d, %d) to (%d, %d, %d)\n", rseg->seg.start.y, rseg->seg.start.z, rseg->seg.start.x,
		rseg->seg.end.y, rseg->seg.end.z, rseg->seg.end.x);
*/
	assert(segment_is_actually_in_routing(rt, rseg));

	struct usage_matrix *m = create_usage_matrix(cp, rt, 2);
	assert(in_usage_bounds(m, rseg->seg.start->coordinate));
	assert(in_usage_bounds(m, rseg->seg.end->coordinate));

	unsigned int usage_size = m->d.x * m->d.y * m->d.z;

	// costs
	unsigned int *cost = malloc(usage_size * sizeof(unsigned int));
	memset(cost, 0xff, usage_size * sizeof(unsigned int)); // set all to highest value

	unsigned char *heap_touched = malloc(usage_size * sizeof(unsigned char));
	memset(heap_touched, 0, usage_size * sizeof(unsigned char));

	enum backtrace *bt = calloc(usage_size, sizeof(enum backtrace));
	for (int i = 0; i < usage_size; i++)
		bt[i] = BT_NONE;

	if (heap)
		clear_cost_coord_heap(heap);
	else
		heap = create_cost_coord_heap();

	// when considering violations, ignore the start pins
	/*
	struct coordinate s1, s2, e1, e2;
	s1 = s2 = rseg->seg.start->coordinate; s2.y--;
	e1 = e2 = rseg->seg.end->coordinate; e2.y--;
	struct coordinate ignore[] = {s1, s2, e1, e2};
	*/

	struct coordinate start = extend_pin(rseg->seg.start);
	struct coordinate end = extend_pin(rseg->seg.end);

	struct cost_coord first = {0, start};
	cost_coord_heap_insert(heap, first);
	cost[usage_idx(m, first.coord)] = 0;
	bt[usage_idx(m, first.coord)] = BT_START;
	heap_touched[usage_idx(m, first.coord)] = 1;

	// disallow vias here
	struct coordinate start_above = start;
	start_above.y += 3;
	struct coordinate end_above = end;
	end_above.y += 3;
	cost[usage_idx(m, start_above)] = ~0L;
	heap_touched[usage_idx(m, start_above)] = 1;
	cost[usage_idx(m, end_above)] = ~0L;
	heap_touched[usage_idx(m, end_above)] = 1;

#ifdef MAZE_ROUTE_DEBUG
	printf("[maze_route] starting...\n");
#endif
	while (heap->n_elts > 0 && !interrupt_routing) {
		struct coordinate c = cost_coord_heap_delete_min(heap).coord;
#ifdef MAZE_ROUTE_DEBUG
		sleep(1);
		printf("[maze_route] heap is\n");
		for (int i = 1; i < heap->n_elts + 1; i++) {
			printf("\t%4d:", i);
			for (int j = i; j > 0; j /= 2)
				putchar(' ');
			printf(" {coord: (%d, %d, %d), cost: %d}\n", heap->elts[i].coord.y, heap->elts[i].coord.z, heap->elts[i].coord.x, heap->elts[i].cost);
		}
		printf("\n");
		printf("[maze_route] popping (%d, %d, %d)\n", c.y, c.z, c.x);
#endif

		// from the current point, examine candidate neighbors
		for (int i = 0; i < sizeof(movement_offsets) / sizeof(struct coordinate); i++) {

			struct coordinate cc = coordinate_add(c, movement_offsets[i]);

			int movement_cost = backtraces[i] == BT_UP || backtraces[i] == BT_DOWN ? 10 :
			                    c.y == 3 ? 3 : 1;
			int violation_cost = 1000 + movement_cost;

			// out of bounds?
			if (!in_usage_bounds(m, cc))
				continue;

#ifdef MAZE_ROUTE_DEBUG
			printf("[maze_route] heap_size = %d, visiting (%d, %d, %d)\n", heap->n_elts, cc.y, cc.z, cc.x);
#endif

			// set new cost
			int violation = usage_matrix_violated(m, cc, NULL, 0); // ignore, sizeof(ignore) / sizeof(struct coordinate));
			unsigned int cost_delta = violation ? violation_cost : movement_cost;
			unsigned int new_cost = cost[usage_idx(m, c)] + cost_delta;
#ifdef MAZE_ROUTE_DEBUG
			printf("[maze_route] new_cost is %u + %d = %u (previously %u)\n", cost[usage_idx(m, c)], cost_delta, new_cost, cost[usage_idx(m, cc)]);
#endif

			if (cost[usage_idx(m, cc)] == ~0L || new_cost < cost[usage_idx(m, cc)]) {
#ifdef MAZE_ROUTE_DEBUG
				printf("[maze_route] cost[(%d, %d, %d)] <- %u (previously %u)\n", cc.y, cc.z, cc.x, new_cost, cost[usage_idx(m, cc)]);
#endif
				// printf("[maze_route] cost[usage_idx(m, c)] = %x, cost_delta = %d\n", cost[usage_idx(m, c)], cost_delta);
				// assert(new_cost != 0);
				cost[usage_idx(m, cc)] = new_cost;
				bt[usage_idx(m, cc)] = backtraces[i];

				if (!heap_touched[usage_idx(m, cc)]) {
					heap_touched[usage_idx(m, cc)] = 1;
					struct cost_coord next = {new_cost, cc};
					cost_coord_heap_insert(heap, next);
				}
			}
		}
	}

	free(heap_touched);

	if (interrupt_routing) {
		return;
	}

	struct coordinate current = end;
	assert(bt[usage_idx(m, current)] != BT_NONE);

	int n_coords = 0;
	int coords_size = 4;
	struct coordinate *coords = calloc(coords_size, sizeof(struct coordinate));
	coords[n_coords++] = current;

#ifdef MAZE_ROUTE_DEBUG
	printf("[maze_route] routed with cost %u\n", cost[usage_idx(m, rseg->seg.end->coordinate)]);
	printf("[maze_route] west=%d, east=%d, north=%d, south=%d, up=%d, down=%d, none=%d, start=%d\n",
		BT_WEST, BT_EAST, BT_NORTH, BT_SOUTH, BT_UP, BT_DOWN, BT_NONE, BT_START);
#endif

	int prev_bt = BT_NONE;
	while (!coordinate_equal(current, start)) {
		int bt_ent = bt[usage_idx(m, current)];
		assert(in_usage_bounds(m, current));
		assert(bt_ent != BT_NONE);

#ifdef MAZE_ROUTE_DEBUG
		printf("[maze_route] current=(%d, %d, %d) backtrace = %d\n", current.y, current.z, current.x, bt_ent);
#endif

#if 1
		/* prefer moving in the same direction as opposed to the first thing the backtrace gives you */
		struct coordinate cont = disp_backtrace(current, prev_bt);
		struct coordinate next = disp_backtrace(current, bt_ent);
		if (prev_bt != BT_NONE && in_usage_bounds(m, cont) && cost[usage_idx(m, cont)] <= cost[usage_idx(m, next)]) {
			current = cont;
		} else {
			current = next;
			prev_bt = bt_ent;
		}
#else
		current = disp_backtrace(current, bt_ent);
#endif

		coords[n_coords++] = current;
		if (n_coords >= coords_size) {
			coords_size *= 2;
			coords = realloc(coords, coords_size * sizeof(struct coordinate));
		}
	}

	coords = realloc(coords, n_coords * sizeof(struct coordinate));
	rseg->n_coords = n_coords;
	rseg->coords = coords;

	free(bt);
}

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
	net_t *net_matrix = malloc(usage_size * sizeof(net_t));
	memset(net_matrix, 0, usage_size * sizeof(net_t));

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
					net_matrix[idx] = -1;
				}
			}
		}
	}

	int total_violations = 0;

	/* segments */
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		struct routed_net *rnet = &(rt->routed_nets[i]);
		total_nets += rnet->n_routed_segments;

		for (int j = 0; j < rnet->n_routed_segments; j++) {
			int segment_violations = 0;

			struct routed_segment *rseg = &(rnet->routed_segments[j]);

			struct coordinate s1, s2, e1, e2;
			s1 = s2 = rseg->seg.start->coordinate; s2.y--;
			e1 = e2 = rseg->seg.end->coordinate; e2.y--;
			struct coordinate ignore[] = {s1, s2, e1, e2};

			for (int k = 0; k < rseg->n_coords; k++) {
				struct coordinate c = rseg->coords[k];
				// printf("[crv] c = (%d, %d, %d)\n", c.y, c.z, c.x);

				int block_in_violation = 0;
				for (int m = 0; m < sizeof(check_offsets) / sizeof(struct coordinate); m++) {
					struct coordinate cc = coordinate_add(c, check_offsets[m]);
					// printf("[crv] cc = (%d, %d, %d)\n", cc.y, cc.z, cc.x);

					if (cc.y < 0 || cc.y >= d.y || cc.z < 0 || cc.z >= d.z || cc.x < 0 || cc.x >= d.x) {
						// printf("[crv] oob\n");
						continue;
					}

					// only ignore the start/end pins for first/last blocks on net
					if (k == 0 || k == rseg->n_coords - 1) {
						int skip = 0;
						for (int n = 0; n < sizeof(ignore) / sizeof(struct coordinate); n++)
							if (coordinate_equal(cc, ignore[n]))
								skip++;
						if (skip) {
							// printf("[crv] skip\n");
							continue;
						}
					}

					int idx = (cc.y * d.z * d.x) + (cc.z * d.x) + cc.x;

					// do not mark or it will collide with itself
					if (matrix[idx] && net_matrix[idx] != i) {
						block_in_violation++;
						// printf("[crv] violation\n");
						// printf("[violation] by net %d, seg %d at (%d, %d, %d) with (%d, %d, %d) due to net %d\n", i, j, c.y, c.z, c.x, cc.y, cc.z, cc.x, net_matrix[idx]);
					}
				}

				if (block_in_violation) {
					segment_violations++;
					total_violations++;
				}
			}

			/* second loop actually marks segment in matrix */
			for (int k = 0; k < rseg->n_coords; k++) {
				struct coordinate c = rseg->coords[k];
				int idx = (c.y * d.z * d.x) + (c.z * d.x) + c.x;
				matrix[idx]++;
				net_matrix[idx] = i;

				if (c.y - 1 > 0) {
					idx = (c.y - 1) * d.z * d.x + c.z * d.x + c.x;
					matrix[idx]++;
					net_matrix[idx] = i;
				}
			}

			// printf("[crv] segment_violations = %d\n", segment_violations);

			int segment_score = segment_violations * 1000 + rseg->n_coords;
			if (max_net_score == -1) {
				max_net_score = min_net_score = segment_score;
			} else {
				max_net_score = max(segment_score, max_net_score);
				min_net_score = min(segment_score, min_net_score);
			}
			rseg->score = segment_score;
			// printf("[crv] net %d seg %d score = %d\n", i, j, segment_score);
		}
	}

	free(matrix);
	free(net_matrix);

	return total_violations;
}

/*
static int count_routed_segment_violations(struct routed_segment *rseg)
{
	// XXX redo
	int score = 0;
	assert(usage_matrix);

	// create list of pins to exclude
	struct coordinate s1, s2, e1, e2;
	s1 = s2 = extend_pin(rseg->seg.start);
	e1 = e2 = extend_pin(rseg->seg.end);
	s2.y--;
	e2.y--;

	for (int i = 0; i < rseg->n_coords; i++) {
		for (int j = 0; j < sizeof(check_offsets) / sizeof(struct coordinate); j++) {
			struct coordinate cc = coordinate_add(rseg->coords[i], check_offsets[j]);
			if (!in_usage_bounds(m, cc))
				continue;

			// avoid the start or end pins or their locations immediately below it
			if (coordinate_equal(cc, s1) ||
			    coordinate_equal(cc, s2) ||
			    coordinate_equal(cc, e1) ||
			    coordinate_equal(cc, e2))
				continue;

			if (usage_matrix[usage_idx(m, cc)] > 1)
				score++;
		}
	}

	return score;
}

static int score_routed_segment(struct routed_segment *rseg)
{
	// XXX redo
	assert(usage_matrix);
	int violations = count_routed_segment_violations(rseg);
	int wire_length = rseg->n_coords;

	return violations * 1000 + wire_length;
}
*/

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

/* this must have a violation_matrix computed beforehand to use */
/*
static int total_violations(struct routings *rt)
{
	int violations = 0;

	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++) {
			violations += count_routed_segment_violations(&rt->routed_nets[i].routed_segments[j]);
		}
	}

	return violations;
}
*/

static void score_nets(struct routings *rt)
{
	total_nets = 0;
	int unset = 1;

	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		total_nets += rt->routed_nets[i].n_routed_segments;
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++) {
			int score = 0; // score_routed_segment(&(rt->routed_nets[i].routed_segments[j]));
			rt->routed_nets[i].routed_segments[j].score = score;

			if (unset) {
				unset = 0;
				max_net_score = score;
				min_net_score = score;
			} else {
				max_net_score = max(score, max_net_score);
				min_net_score = min(score, min_net_score);
			}
			// printf("[router] net %d, segment %d has score %d\n", i, j, score);
			// print_routed_segment(&rt->routed_nets[i].routed_segments[j]);
		}
	}
}

/* pointers into a routings */
struct rip_up_set {
	int n_ripped;
	struct routed_segment **rip_up;
};

static struct rip_up_set natural_selection(struct routings *rt)
{
	int rip_up_count = 0;
	int rip_up_size = 4;
	struct routed_segment **rip_up = malloc(rip_up_size * sizeof(struct routed_segment *));
	memset(rip_up, 0, rip_up_size * sizeof(struct routed_segment *));

	int range = max_net_score - min_net_score;
	int lower = min_net_score - (range / 8);
	// int higher = range + (range / 8);

	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++) {
			int score = rt->routed_nets[i].routed_segments[j].score;

			int r = random() % (range / 8 * 10);

			if (r < score - lower) {
#ifdef NATURAL_SELECTION_DEBUG
				printf("[natural_selection] ripping up net %d, segment %d (rand(%d) = %d < %d)\n", i, j, (range / 8 * 10), r, score - lower);
#endif
				// print_routed_segment(&rt->routed_nets[i].routed_segments[j]);
				rip_up[rip_up_count++] = &rt->routed_nets[i].routed_segments[j];
				if (rip_up_count >= rip_up_size) {
					rip_up_size *= 2;
					rip_up = realloc(rip_up, rip_up_size * sizeof(struct routed_segment *));
				}
			} else {
#ifdef NATURAL_SELECTION_DEBUG
				printf("[natural_selection] leaving net %d, segment %d intact (rand(%d) = %d >= %d)\n", i, j, (range / 8 * 10), r, score - lower);
#endif
			}
		}
	}

	struct rip_up_set rus = {rip_up_count, rip_up};

	return rus;
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

static void router_sigint_handler(int a)
{
	printf("Interrupt\n");
	interrupt_routing = 1;
}

int routed_segment_cmp(const void *a, const void *b)
{
	return -1 * (((struct routed_segment *)a)->score - ((struct routed_segment *)b)->score);
}

struct routings *route(struct blif *blif, struct cell_placements *cp)
{
	struct pin_placements *pp = placer_place_pins(cp);
	struct net_pin_map *npm = placer_create_net_pin_map(pp);

	struct routings *rt = initial_route(blif, npm);
	// print_routings(rt);
	recenter(cp, rt, 2);

	int iterations = 0;
	int violations = count_routings_violations(cp, rt);
	printf("\n");
	interrupt_routing = 0;
	signal(SIGINT, router_sigint_handler);
	while ((violations = count_routings_violations(cp, rt)) > 0 && !interrupt_routing) {
		struct rip_up_set rus = natural_selection(rt);
		// sort elements by highest score
		qsort(rus.rip_up, rus.n_ripped, sizeof(struct routed_segment *), routed_segment_cmp);
		printf("\r[router] Iterations: %4d, Violations: %d, Segments to re-route: %d", iterations + 1, violations, rus.n_ripped);
		fflush(stdout);
		for (int i = 0; i < rus.n_ripped; i++) {
			rus.rip_up[i]->n_coords = 0;
			free(rus.rip_up[i]->coords);
			rus.rip_up[i]->coords = NULL;
		}

		for (int i = 0; i < rus.n_ripped; i++) {
			recenter(cp, rt, 2);
			// printf("[router] Rerouting a segment in net %d with score %d\n", rus.rip_up[i]->net->net, rus.rip_up[i]->score);
			maze_reroute_segment(cp, rt, rus.rip_up[i]);
			// print_routed_segment(rus.rip_up[i]);
		}
		free(rus.rip_up);

		iterations++;
	}
	signal(SIGINT, SIG_DFL);
	printf("\n[router] Routing complete!\n");
	print_routings(rt);

/*
	struct routed_segment *rseg = &rt->routed_nets[1].routed_segments[0];
	maze_reroute_segment(cp, rt, rseg);

	create_usage_matrix(cp, rt, 0);
	int score = score_routed_segment(rseg);
	print_routed_segment(rseg);
	printf("[router] new score: %d\n", score);
*/

	free_pin_placements(pp);
	// free_net_pin_map(npm);

	return rt;
}
