#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "extract.h"
#include "segment.h"
#include "placer.h"
#include "router.h"
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

static int interrupt_routing = 0;

/* connect the two nets with a rectilinear path*/
struct routed_segment cityblock_route(net_t net, struct segment s)
{
	struct coordinate a = s.start;
	struct coordinate b = s.end;

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

	struct routed_segment rseg = {s, len, path};
	return rseg;
}

void print_routed_segment(struct routed_segment *rseg)
{
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

struct coordinate extend_pin(enum ordinal_direction facing, struct coordinate c)
{
	switch (facing) {
		case EAST: c.x++; break;
		case WEST: c.x--; break;
		case NORTH: c.z++; break;
		case SOUTH: c.z--; break;
		default: break;
	}

	return c;
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
		struct coordinate c = npm->pins[net][0].coordinate;
#ifdef DUMB_ROUTE_DEBUG
		printf("[dumb_route] net %d (%s) has 1 pin at (%d, %d, %d)\n", net, get_net_name(blif, net), c.y, c.z, c.x);
#endif
		struct segment seg = {c, c};
		struct coordinate *coords = malloc(sizeof(struct coordinate));
		coords[0] = c;

		struct routed_segment rseg = {seg, 1, coords};
		
		rn->routed_segments = malloc(sizeof(struct routed_segment));
		rn->routed_segments[0] = rseg;
	} else if (n_pins == 2) {
		struct coordinate a = extend_pin(npm->pins[net][0].cell_pin->facing, npm->pins[net][0].coordinate);
		struct coordinate b = extend_pin(npm->pins[net][1].cell_pin->facing, npm->pins[net][1].coordinate);
#ifdef DUMB_ROUTE_DEBUG
		printf("[dumb_route] net %d (%s):\n  ", net, get_net_name(blif, net));
#endif
		struct segment s = {a, b};
		rn->n_routed_segments = 1;
		rn->routed_segments = malloc(sizeof(struct routed_segment));
		rn->routed_segments[0] = cityblock_route(net, s);
#ifdef DUMB_ROUTE_DEBUG
		print_routed_segment(rn);
#endif
	} else {
		struct coordinate *coords = calloc(n_pins, sizeof(struct coordinate));
		for (int j = 0; j < n_pins; j++)
			coords[j] = extend_pin(npm->pins[net][j].cell_pin->facing, npm->pins[net][j].coordinate);

		struct segments *mst = create_mst(coords, n_pins);
		rn->n_routed_segments = n_pins - 1;
		rn->routed_segments = malloc(sizeof(struct routed_segment) * rn->n_routed_segments);
#ifdef DUMB_ROUTE_DEBUG
		printf("[dumb_route] net %d (%s) has %d pins:\n", net, get_net_name(blif, net), n_pins);
#endif
		for (int j = 0; j < mst->n_segments; j++) {
			struct routed_segment rseg = cityblock_route(net, mst->segments[j]);
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
		}
	}

	return new_rt;
}

struct dimensions compute_routings_dimensions(struct routings *rt)
{
	struct dimensions d = {0, 0, 0};
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		struct routed_net rn = rt->routed_nets[i];
		for (int j = 0; j < rn.n_routed_segments; j++) {
			struct routed_segment rseg = rn.routed_segments[j];
			for (int k = 0; k < rseg.n_coords; k++) {
				struct coordinate c = rseg.coords[k];
				d.y = max(d.y, c.y + 1);
				d.z = max(d.z, c.z + 1);
				d.x = max(d.x, c.x + 1);
			}

			struct coordinate c = rseg.seg.start;
			d.y = max(d.y, c.y + 1);
			d.z = max(d.z, c.z + 1);
			d.x = max(d.x, c.x + 1);

			c = rseg.seg.end;
			d.y = max(d.y, c.y + 1);
			d.z = max(d.z, c.z + 1);
			d.x = max(d.x, c.x + 1);
		}
	}

	return d;
}

void routings_displace(struct routings *rt, struct coordinate disp)
{
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		struct routed_net *rn = &(rt->routed_nets[i]);
		for (int j = 0; j < rn->n_routed_segments; j++) {
			struct routed_segment *rseg = &(rn->routed_segments[j]);
			rseg->seg.start = coordinate_add(rseg->seg.start, disp);
			rseg->seg.end = coordinate_add(rseg->seg.end, disp);
			for (int k = 0; k < rseg->n_coords; k++) {
				rseg->coords[k] = coordinate_add(rseg->coords[k], disp);
			}
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

static unsigned char *usage_matrix = NULL;
static struct dimensions usage_dimensions;
static int usage_size = 0;
static int usage_xz_margin = 0;

static unsigned char *violation_matrix = NULL;

static int in_usage_bounds(struct coordinate c)
{
	return c.x >= 0 && c.x < usage_dimensions.x &&
	       c.y >= 0 && c.y < usage_dimensions.y &&
	       c.z >= 0 && c.z < usage_dimensions.z;
}

static int usage_idx(struct coordinate c)
{
	return (c.y * usage_dimensions.z * usage_dimensions.x) + (c.z * usage_dimensions.x) + c.x;
}

/* creates a violation matrix;
 * use in_usage_bounds() and usage_idx() to check bounds and index
 * into the matrix, respectively.
 *
 * this requires that the top-left-most x and z points are at least
 * `xz_margin` out; use recenter() with a non-zero xz_margin to displace.
 */
static void create_usage_matrix(struct cell_placements *cp, struct routings *rt, int xz_margin)
{
	/* ensure no coordinate is < 0 */
	struct coordinate tlcp = placements_top_left_most_point(cp);
	struct coordinate tlrt = routings_top_left_most_point(rt);
	struct coordinate top_left_most = coordinate_piecewise_min(tlcp, tlrt);
	assert(top_left_most.x >= xz_margin && top_left_most.y >= 0 && top_left_most.z >= xz_margin);

	struct dimensions d = dimensions_piecewise_max(compute_placement_dimensions(cp), compute_routings_dimensions(rt));
	d.x += xz_margin;
	d.z += xz_margin;
	d.y = max(d.y, 4); // allow for routing on y=0 and y=3

	usage_dimensions = d;
	usage_xz_margin = xz_margin;

	usage_size = d.x * d.y * d.z;

	usage_matrix = realloc(usage_matrix, usage_size * sizeof(unsigned char));
	violation_matrix = realloc(violation_matrix, usage_size * sizeof(unsigned char));
	memset(usage_matrix, 0, usage_size * sizeof(unsigned char));
	memset(violation_matrix, 0, usage_size * sizeof(unsigned char));

	/* placements */
	for (int i = 0; i < cp->n_placements; i++) {
		struct placement p = cp->placements[i];
		struct coordinate c = p.placement;
		struct dimensions pd = p.cell->dimensions[p.turns];

		int cell_x = xz_margin + c.x + pd.x;
		int cell_y = c.y + pd.y;
		int cell_z = xz_margin + c.z + pd.z;

		for (int y = c.y; y < cell_y; y++) {
			for (int z = max(0, c.z - 1) + xz_margin; z < min(d.z, cell_z); z++) {
				for (int x = max(0, c.x - 1) + xz_margin; x < min(d.x, cell_x); x++) {
					struct coordinate cc = {y, z, x};
					int idx = usage_idx(cc);

					// only update usage_matrix within the cell bounds
					if (z >= c.z && z < cell_z && x >= c.x && x < cell_x)
						usage_matrix[idx]++;

					// only update violation_matrix within matrix bounds
					if (in_usage_bounds(cc))
						violation_matrix[idx] = 1;
				}
			}
		}

/*
		for (int j = 0; j < p.cell->n_pins; j++) {
			struct coordinate pc = p.cell->pins[p.turns][j].coordinate;
			struct coordinate pcc = coordinate_add(c, pc);
			switch (p.cell->pins[p.turns][j].facing) {
				case EAST: pcc.x++; break;
				case WEST: pcc.x--; break;
				case NORTH: pcc.z++; break;
				case SOUTH: pcc.z--; break;
				default: break;
			}

			if (in_usage_bounds(pcc))
				violation_matrix[usage_idx(pcc)] = 0;
		}
*/
	}

	/* routings */
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++) {
			for (int k = 0; k < rt->routed_nets[i].routed_segments[j].n_coords; k++) {
				// here
				struct coordinate c = rt->routed_nets[i].routed_segments[j].coords[k];
				usage_matrix[usage_idx(c)]++;

				// below
				c.y--;
				if (in_usage_bounds(c))
					usage_matrix[usage_idx(c)]++;

				// violation matrix calculation
				for (int m = 0; m < sizeof(check_offsets) / sizeof(struct coordinate); m++) {
					struct coordinate cc = coordinate_add(c, check_offsets[m]);
					if (in_usage_bounds(cc))
						violation_matrix[usage_idx(cc)] = 1;
				}
			}
		}
	}
}

static struct routed_segment *segment_is_actually_in_routing(struct routings *rt, struct routed_segment *rseg)
{
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++)
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++)
			if (rseg == &rt->routed_nets[i].routed_segments[j])
				return rseg;
	return NULL;
}

struct cost_coord {
	unsigned int cost;
	struct coordinate coord;
};

static int cost_coord_cmp(const void *a, const void *b)
{
	// make heapsort sort from largest to smallest
	return -1 * (((struct cost_coord *)a)->cost - ((struct cost_coord *)b)->cost);
}

#define HEAP_INITIAL_SIZE 32

static void maze_reroute_segment(struct cell_placements *cp, struct routings *rt, struct routed_segment *rseg)
{
	static struct cost_coord *heap = NULL;
	static int heap_size = 0;

	/* rip out this segment */
	rseg->n_coords = 0;
	if (rseg->coords)
		free(rseg->coords);
	rseg->coords = NULL;

	/* the additional xz_margin to each side for a possible route is 2 */
	printf("[maze_route] attempting to maze-route from (%d, %d, %d) to (%d, %d, %d)\n", rseg->seg.start.y, rseg->seg.start.z, rseg->seg.start.x,
		rseg->seg.end.y, rseg->seg.end.z, rseg->seg.end.x);
	assert(segment_is_actually_in_routing(rt, rseg));

	create_usage_matrix(cp, rt, 2);
	assert(in_usage_bounds(rseg->seg.start));
	assert(in_usage_bounds(rseg->seg.end));

	// costs
	unsigned int *cost = malloc(usage_size * sizeof(unsigned int));
	memset(cost, 0xff, usage_size * sizeof(unsigned int)); // set all to highest value

	// backtracing
	enum backtrace {BT_NONE, BT_WEST, BT_SOUTH, BT_EAST, BT_NORTH, BT_DOWN, BT_UP, BT_START};
	enum backtrace backtraces[] = {BT_WEST, BT_SOUTH, BT_EAST, BT_NORTH, BT_DOWN, BT_UP};
	struct coordinate movement_offsets[] = {
		{0, 0, 1}, {0, 1, 0}, {0, 0, -1}, {0, -1, 0}, {3, 0, 0}, {-3, 0, 0}
	};

	enum backtrace *bt = calloc(usage_size, sizeof(enum backtrace));
	for (int i = 0; i < usage_size; i++)
		bt[i] = BT_NONE;

	// heap
	int heap_count = 0;

	if (!heap) {
		heap = calloc(HEAP_INITIAL_SIZE, sizeof(struct cost_coord));
		heap_size = HEAP_INITIAL_SIZE;
	} else {
		memset(heap, 0, heap_size * sizeof(struct cost_coord));
	}

	struct cost_coord first = {0, rseg->seg.start};
	heap[heap_count++] = first;
	cost[usage_idx(first.coord)] = 0;
	bt[usage_idx(first.coord)] = BT_START;

#ifdef MAZE_ROUTE_DEBUG
	printf("[maze_route] starting...\n");
#endif
	while (heap_count > 0 && !interrupt_routing) {
		struct coordinate c = heap[--heap_count].coord;
#ifdef MAZE_ROUTE_DEBUG
		printf("[maze_route] popping (%d, %d, %d)\n", c.y, c.z, c.x);
#endif
		int movement_cost = 1;
		int violation_cost = 1000 + movement_cost;

		// from the current point, examine candidate neighbors
		for (int i = 0; i < sizeof(movement_offsets) / sizeof(struct coordinate); i++) {
			struct coordinate cc = coordinate_add(c, movement_offsets[i]);

			// out of bounds?
			if (!in_usage_bounds(cc))
				continue;

#ifdef MAZE_ROUTE_DEBUG
			printf("[maze_route] heap_size = %d, visiting (%d, %d, %d) cost = %u\n", heap_count, cc.y, cc.z, cc.x, heap[heap_count].cost);
#endif

			// set new cost
			unsigned int cost_delta = violation_matrix[usage_idx(cc)] ? violation_cost : movement_cost;
			unsigned int new_cost = cost[usage_idx(c)] + cost_delta;
#ifdef MAZE_ROUTE_DEBUG
			printf("[maze_route] new_cost is %u + %d = %u (previously %u)\n", cost[usage_idx(c)], cost_delta, new_cost, cost[usage_idx(cc)]);
#endif

			if (cost[usage_idx(cc)] == ~0L || new_cost < cost[usage_idx(cc)]) {
#ifdef MAZE_ROUTE_DEBUG
				printf("[maze_route] cost[(%d, %d, %d)] <- %u (previously %u)\n", cc.y, cc.z, cc.x, new_cost, cost[usage_idx(cc)]);
#endif
				// printf("[maze_route] cost[usage_idx(c)] = %x, cost_delta = %d\n", cost[usage_idx(c)], cost_delta);
				// assert(new_cost != 0);
				cost[usage_idx(cc)] = new_cost;
				bt[usage_idx(cc)] = backtraces[i];

				// check presence in heap
				int in_heap = 0;
				for (int j = 0; j < heap_count; j++) {
					if (coordinate_equal(heap[j].coord, cc)) {
						in_heap = 1;
						break;
					}
				}

				if (!in_heap) {
					struct cost_coord next = {new_cost, cc};
					heap[heap_count++] = next;
					// printf("[maze_route] adding (%d, %d, %d)\n", cc.y, cc.z, cc.x);
					if (heap_count >= heap_size) {
						printf("[maze_route] doubling heap from %d\n", heap_size);
						heap_size *= 2;
						heap = realloc(heap, heap_size * sizeof(struct cost_coord));
					}
				}
			}
		}

		heapsort(heap, heap_count, sizeof(struct cost_coord), cost_coord_cmp);
	}

	if (interrupt_routing) {
		return;
	}

	struct coordinate current = rseg->seg.end;
	assert(bt[usage_idx(current)] != BT_NONE);

	int n_coords = 0;
	int coords_size = 4;
	struct coordinate *coords = calloc(coords_size, sizeof(struct coordinate));
	coords[n_coords++] = current;

#ifdef MAZE_ROUTE_DEBUG
	printf("[maze_route] routed with cost %u\n", cost[usage_idx(rseg->seg.end)]);
	printf("[maze_route] west=%d, east=%d, north=%d, south=%d, up=%d, down=%d, none=%d, start=%d\n",
		BT_WEST, BT_EAST, BT_NORTH, BT_SOUTH, BT_UP, BT_DOWN, BT_NONE, BT_START);
#endif

	while (!coordinate_equal(current, rseg->seg.start)) {
		int bt_ent = bt[usage_idx(current)];
		assert(in_usage_bounds(current));
		assert(bt_ent != BT_NONE);

#ifdef MAZE_ROUTE_DEBUG
		printf("[maze_route] current=(%d, %d, %d) backtrace = %d\n", current.y, current.z, current.x, bt_ent);
#endif

		switch (bt_ent) {
		case BT_WEST:
			current.x--; break;
		case BT_EAST:
			current.x++; break;
		case BT_NORTH:
			current.z++; break;
		case BT_SOUTH:
			current.z--; break;
		case BT_UP:
			current.y += 3; break;
		case BT_DOWN:
			current.y -= 3; break;
		default:
			break;
		}

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

static int count_routed_segment_violations(struct routed_segment *rseg)
{
	int score = 0;
	assert(usage_matrix);

	/* for each, but without start and end */
	for (int i = 1; i < rseg->n_coords - 1; i++) {
		if (violation_matrix[usage_idx(rseg->coords[i])])
			score++;
	}

	return score;
}

static int score_routed_segment(struct routed_segment *rseg)
{
	assert(usage_matrix);
	int violations = count_routed_segment_violations(rseg);
	int wire_length = rseg->n_coords;

	return violations * 1000 + wire_length;
}

static struct routings *initial_route(struct blif *blif, struct net_pin_map *npm)
{
	struct routings *rt = malloc(sizeof(struct routings));
	rt->n_routed_nets = npm->n_nets;
	rt->routed_nets = calloc(rt->n_routed_nets + 1, sizeof(struct routed_net));

	for (net_t i = 1; i < npm->n_nets + 1; i++)
		rt->routed_nets[i] = *dumb_route(blif, npm, i);

	return rt;
}

/* this must have a violation_matrix computed beforehand to use */
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

static int max_net_score = 0;
static int min_net_score = 0;
static int total_nets = 0;

static void score_nets(struct routings *rt)
{
	total_nets = 0;
	int unset = 1;

	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		total_nets += rt->routed_nets[i].n_routed_segments;
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++) {
			int score = score_routed_segment(&(rt->routed_nets[i].routed_segments[j]));
			rt->routed_nets[i].routed_segments[j].score = score;

			if (unset) {
				unset = 0;
				max_net_score = score;
				min_net_score = score;
			} else {
				max_net_score = max(score, max_net_score);
				min_net_score = min(score, min_net_score);
			}
			printf("[router] net %d, segment %d has score %d\n", i, j, score);
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
				printf("[natural_selection] ripping up net %d, segment %d (rand(%d) = %d < %d)\n", i, j, (range / 8 * 10), r, score - lower);
				// print_routed_segment(&rt->routed_nets[i].routed_segments[j]);
				rip_up[rip_up_count++] = &rt->routed_nets[i].routed_segments[j];
				if (rip_up_count >= rip_up_size) {
					rip_up_size *= 2;
					rip_up = realloc(rip_up, rip_up_size * sizeof(struct routed_segment *));
				}
			} else {
				printf("[natural_selection] leaving net %d, segment %d intact (rand(%d) = %d >= %d)\n", i, j, (range / 8 * 10), r, score - lower);
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

struct routings *route(struct blif *blif, struct cell_placements *cp)
{
	struct pin_placements *pp = placer_place_pins(cp);
	struct net_pin_map *npm = placer_create_net_pin_map(pp);

	struct routings *rt = initial_route(blif, npm);
	print_routings(rt);
	recenter(cp, rt, 2);

	create_usage_matrix(cp, rt, 2);
	score_nets(rt);

	int iterations = 0;
	int violations = 0;
	printf("\n");
	interrupt_routing = 0;
	signal(SIGINT, router_sigint_handler);
	while ((violations = total_violations(rt)) > 0 && !interrupt_routing) {
		printf("[router] Iterations: %4d, Violations: %d\n", iterations + 1, violations);
		struct rip_up_set rus = natural_selection(rt);
		for (int i = 0; i < rus.n_ripped; i++) {
			rus.rip_up[i]->n_coords = 0;
			free(rus.rip_up[i]->coords);
			rus.rip_up[i]->coords = NULL;
		}

		for (int i = 0; i < rus.n_ripped; i++) {
			recenter(cp, rt, 2);
			maze_reroute_segment(cp, rt, rus.rip_up[i]);
			print_routed_segment(rus.rip_up[i]);
		}
		free(rus.rip_up);

		create_usage_matrix(cp, rt, 0);
		score_nets(rt);
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
	free_net_pin_map(npm);

	return rt;
}
