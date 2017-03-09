#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

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
	memcpy(new_rt->routed_nets, old_rt->routed_nets, (new_rt->n_routed_nets + 1) * sizeof(struct routed_net));

	/* for each routed_net in routings */
	for (net_t i = 1; i < new_rt->n_routed_nets + 1; i++) {
		struct routed_net *rn = &(new_rt->routed_nets[i]);
		struct routed_net on = old_rt->routed_nets[i];
		rn->n_pins = on.n_pins;
		memcpy(rn->pins, on.pins, sizeof(struct placed_pin) * rn->n_pins);

		rn->n_routed_segments = on.n_routed_segments;
		rn->routed_segments = malloc(sizeof(struct routed_segment) * rn->n_routed_segments);
		memcpy(rn->routed_segments, on.routed_segments, sizeof(struct routed_segment) * rn->n_routed_segments);
		
		/* for each routed_segment in routed_net */
		for (int j = 0; j < rn->n_routed_segments; j++) {
			struct routed_segment *rseg = &(rn->routed_segments[j]);
			struct routed_segment old_rseg = on.routed_segments[j];
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
				printf("[crd] y=%d z=%d x=%d\n", c.y, c.z, c.x);
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
				rseg->coords[j] = add_coordinates(rseg->coords[j], disp);
			}
		}
	}
}

void free_routings(struct routings *rt)
{
	free(rt->routed_nets);
	free(rt);
}

static struct routings *initial_route(struct blif *blif, struct net_pin_map *npm)
{
	struct routings *rt = malloc(sizeof(struct routings));
	rt->n_routed_nets = npm->n_nets;
	rt->routed_nets = calloc(rt->n_routed_nets + 1, sizeof(struct routed_net));

	for (net_t i = 1; i < npm->n_nets; i++)
		rt->routed_nets[i] = *dumb_route(blif, npm, i);

	return rt;
}

struct routings *route(struct blif *blif, struct cell_placements *cp)
{
	struct pin_placements *pp = placer_place_pins(cp);
	struct net_pin_map *npm = placer_create_net_pin_map(pp);

	struct routings *rt = initial_route(blif, npm);

	free_pin_placements(pp);
	free_net_pin_map(npm);

	return rt;
}
