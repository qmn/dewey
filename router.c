#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "extract.h"
#include "segment.h"
#include "placer.h"
#include "router.h"
#include "blif.h"

static struct coordinate add_coordinates(struct coordinate a, struct coordinate b)
{
	struct coordinate c = {a.y + b.y, a.z + b.z, a.x + b.x};
	return c;
}

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
		struct coordinate a = npm->pins[net][0].coordinate;
		struct coordinate b = npm->pins[net][1].coordinate;
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
			coords[j] = npm->pins[net][j].coordinate;

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
			for (int k = 0; k < rseg->n_coords; k++) {
				rseg->coords[k] = add_coordinates(rseg->coords[k], disp);
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

static int in_usage_bounds(struct coordinate c)
{
	return c.x >= 0 && c.x < usage_dimensions.x &&
	       c.y >= 0 && c.y < usage_dimensions.y &&
	       c.z >= 0 && c.z < usage_dimensions.z;
}

static int usage_idx(struct coordinate c)
{
	return c.y * usage_dimensions.z * usage_dimensions.x + (c.z + usage_xz_margin) * usage_dimensions.x + (c.x + usage_xz_margin);
}

static void create_usage_matrix(struct cell_placements *cp, struct routings *rt, int xz_margin)
{
	/* ensure no coordinate is < 0 */
	struct coordinate tlcp = placements_top_left_most_point(cp);
	struct coordinate tlrt = routings_top_left_most_point(rt);
	struct coordinate top_left_most = coordinate_piecewise_min(tlcp, tlrt);
	assert(top_left_most.x >= 0 && top_left_most.y >= 0 && top_left_most.z >= 0);

	struct dimensions d = dimensions_piecewise_max(compute_placement_dimensions(cp), compute_routings_dimensions(rt));
	d.x += xz_margin * 2;
	d.z += xz_margin * 2;

	usage_dimensions = d;
	usage_xz_margin = xz_margin;

	usage_size = d.x * d.y * d.z;

	usage_matrix = realloc(usage_matrix, usage_size * sizeof(unsigned char));
	memset(usage_matrix, 0, usage_size * sizeof(unsigned char));

	/* placements */
	for (int i = 0; i < cp->n_placements; i++) {
		struct placement p = cp->placements[i];
		struct coordinate c = p.placement;
		struct dimensions pd = p.cell->dimensions[p.turns];

		int cell_x = xz_margin + c.x + pd.x;
		int cell_y = c.y + pd.y;
		int cell_z = xz_margin + c.z + pd.z;

		for (int y = c.y; y < cell_y; y++) {
			int by = y * d.z * d.x;
			for (int z = max(0, c.z - 1); z < min(cell_z, d.z + 1); z++) {
				int bz = z * d.x;
				for (int x = max(0, c.x - 1); x < min(cell_x, d.x + 1); x++) {
					int idx = by + bz + x;
					assert(x >= 0 && x < usage_size);
					usage_matrix[idx]++;
				}
			}
		}
	}

	/* routings */
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++) {
			for (int k = 0; k < rt->routed_nets[i].routed_segments[j].n_coords; k++) {
				struct coordinate c = rt->routed_nets[i].routed_segments[j].coords[k];
				for (int l = 0; l < sizeof(check_offsets) / sizeof(struct coordinate); l++) {
					struct coordinate c_off = check_offsets[l];
					struct coordinate off = coordinate_add(c, c_off);
					if (!in_usage_bounds(off))
						continue;
					usage_matrix[usage_idx(off)]++;
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

static void maze_reroute_segment(struct cell_placements *cp, struct routings *rt, struct routed_segment *rseg)
{
	/* rip out this segment */
	rseg->n_coords = 0;
	free(rseg->coords);
	rseg->coords = NULL;

	/* the additional xz_margin to each side for a possible route is 2 */
	assert(segment_is_actually_in_routing(rt, rseg));
	create_usage_matrix(cp, rt, 2);

	// costs
	unsigned int *cost = malloc(usage_size * sizeof(unsigned int));
	memset(cost, 0xff, usage_size * sizeof(unsigned int)); // set all to highest value

	// backtracing
	enum backtrace {BT_NONE, BT_WEST, BT_SOUTH, BT_EAST, BT_NORTH, BT_DOWN, BT_UP};
	enum backtrace backtraces[] = {BT_WEST, BT_SOUTH, BT_EAST, BT_NORTH, BT_DOWN, BT_UP};
	enum movement {MV_NONE, MV_EAST, MV_NORTH, MV_WEST, MV_SOUTH, MV_UP, MV_DOWN};
	enum movement movements[] = {MV_EAST, MV_NORTH, MV_WEST, MV_SOUTH, MV_UP, MV_DOWN};
	struct coordinate movement_offsets[] = {
		{0, 0, 1}, {0, 1, 0}, {0, 0, -1}, {0, -1, 0}, {3, 0, 0}, {-3, 0, 0}
	};

	enum backtrace *bt = calloc(usage_size, sizeof(enum backtrace));
	for (int i = 0; i < usage_size; i++)
		bt[i] = BT_NONE;
	enum movement *mv = calloc(usage_size, sizeof(enum movement));

	// heap
	int heap_size = 32;
	int heap_count = 0;

	struct cost_coord *heap = calloc(heap_size, sizeof(struct cost_coord));
	struct cost_coord first = {0, rseg->seg.start};
	heap[heap_count++] = first;
	cost[usage_idx(first.coord)] = 0;

	printf("[maze_route] starting...\n");
	while (heap_count > 0) {
		struct coordinate c = heap[--heap_count].coord;
		int violation_cost = 1000;
		int movement_cost = 1;

		for (int i = 0; i < sizeof(movement_offsets) / sizeof(struct coordinate); i++) {
			struct coordinate cc = coordinate_add(c, movement_offsets[i]);
			printf("[maze_route] heap_size = %d, visiting (%d, %d, %d) cost = %u\n", heap_count, cc.y, cc.z, cc.x, heap[heap_count].cost);

			// out of bounds?
			if (!in_usage_bounds(cc))
				continue;

			// already visited?
			if (bt[usage_idx(cc)] != BT_NONE)
				continue;

			// set new cost
			unsigned int cost_delta = usage_matrix[usage_idx(cc)] ? violation_cost : movement_cost;
			unsigned int new_cost = cost[usage_idx(c)] + cost_delta;
			printf("[maze_route] new_cost is %u\n", new_cost);

			// see if this new path is more advantageous
			if (new_cost < cost[usage_idx(cc)]) {
				cost[usage_idx(cc)] = new_cost;
				bt[usage_idx(cc)] = backtraces[i];
			}

			struct cost_coord next = {new_cost, cc};
			heap[heap_count++] = next;
			if (heap_count >= heap_size) {
				printf("[maze_route] doubling heap from %d\n", heap_size);
				heap_size *= 2;
				heap = realloc(heap, heap_size * sizeof(struct cost_coord));
			}
		}

		heapsort(heap, heap_count, sizeof(struct cost_coord), cost_coord_cmp);
	}

	struct coordinate current = rseg->seg.end;
	if (bt[usage_idx(current)] != BT_NONE)
		printf("\nsuccess!\n");

	int n_coords = 0;
	int coords_size = 4;
	struct coordinate *coords = calloc(coords_size, sizeof(struct coordinate));
	coords[n_coords++] = current;

	while (!coordinate_equal(current, rseg->seg.start)) {
		int bt_ent = bt[usage_idx(current)];
		assert(bt_ent != BT_NONE);

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

	free(heap);
	free(bt);
	free(mv);
}

static int count_routed_segment_violations(struct routed_segment *rseg)
{
	int score = 0;

	/* temporarily remove self from usage matrix */
/*
	for (int i = 0; i < rseg->n_coords; i++) {
		for (int j = 0; j < sizeof(check_offsets) / sizeof(struct coordinate); j++)
		// here
		struct coordinate c = rseg->coords[i];
		usage_matrix[usage_idx(c)]--;

		// below
		struct coordinate cc = c;
		cc.y--;
		if (in_usage_bounds(cc))
			usage_matrix[usage_idx(cc)]--;
	}
*/

	assert(usage_matrix);
	/* for each */
	struct coordinate *prev = NULL;
	struct coordinate *next = NULL;
	for (int i = 0; i < rseg->n_coords; i++) {
		prev = i == 0 ? NULL : &rseg->coords[i-1];
		next = i == rseg->n_coords - 1 ? NULL : &rseg->coords[i+1];

		for (int j = 0; j < sizeof(check_offsets) / sizeof(struct coordinate); j++) {
			struct coordinate check_coord = coordinate_add(rseg->coords[i], check_offsets[j]);
			if (prev && coordinate_equal(check_coord, *prev))
				continue;

			if (next && coordinate_equal(check_coord, *next))
				continue;

			if (in_usage_bounds(check_coord) && usage_matrix[usage_idx(check_coord)] > 1)
				score++;
		}
	}

	/* re-add self to usage matrix */
/*
	for (int i = 0; i < rseg->n_coords; i++) {
		struct coordinate c = rseg->coords[i];
		usage_matrix[usage_idx(c)]++;

		struct coordinate cc = c;
		cc.y--;
		if (in_usage_bounds(cc))
			usage_matrix[usage_idx(cc)]++;
	}
*/

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

struct routings *route(struct blif *blif, struct cell_placements *cp)
{
	struct pin_placements *pp = placer_place_pins(cp);
	struct net_pin_map *npm = placer_create_net_pin_map(pp);

	struct routings *rt = initial_route(blif, npm);
/*
	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++) {
			print_routed_segment(&rt->routed_nets[i].routed_segments[j]);
		}
	}
*/
	recenter(cp, rt);

	create_usage_matrix(cp, rt, 0);

	for (net_t i = 1; i < rt->n_routed_nets + 1; i++) {
		for (int j = 0; j < rt->routed_nets[i].n_routed_segments; j++) {
			int score = score_routed_segment(&(rt->routed_nets[i].routed_segments[j]));
			printf("[router] net %d, segment %d has score %d\n", i, j, score);
			// print_routed_segment(&rt->routed_nets[i].routed_segments[j]);
		}
	}

/*
	struct routed_segment *rseg = &rt->routed_nets[1].routed_segments[0];
	maze_reroute_segment(cp, rt, rseg);

	create_usage_matrix(cp, rt, 0);
	int score = score_routed_segment(rseg);
	printf("[router] new score: %d\n", score);
*/

	free_pin_placements(pp);
	free_net_pin_map(npm);

	return rt;
}
