#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

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

/* connect the two nets with a rectilinear path*/
struct routed_net *cityblock_route(net_t net, struct coordinate a, struct coordinate b)
{
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

	struct routed_net *rn = malloc(sizeof(struct routed_net));
	rn->net = net;
	rn->n_coords = len;
	rn->coords = path;

	return rn;
}

void print_routed_net(struct routed_net *rn)
{
	for (int i = 0; i < rn->n_coords; i++)
		printf("(%d, %d, %d) ", rn->coords[i].y, rn->coords[i].z, rn->coords[i].x);
	printf("\n");
}

void free_routed_net(struct routed_net *rn)
{
	free(rn->coords);
	free(rn);
}

/* generate the MST for this net to determine the order of connections,
 * then connect them all with a city*/
struct routed_net *dumb_route(struct blif *blif, struct net_pin_map *npm, net_t net)
{
	int n_pins = npm->n_pins_for_net[net];
	struct routed_net *rn = malloc(sizeof(struct routed_net));
	rn->net = net;
	rn->n_coords = n_pins;
	rn->coords = NULL;

	if (n_pins == 0) {
#ifdef DUMB_ROUTE_DEBUG
		printf("[dumb_route] net %d (%s) has 0 pins\n", net, get_net_name(blif, net));
#endif
	} else if (n_pins == 1) {
		struct coordinate c = npm->pins[net][0].coordinate;
#ifdef DUMB_ROUTE_DEBUG
		printf("[dumb_route] net %d (%s) has 1 pin at (%d, %d, %d)\n", net, get_net_name(blif, net), c.y, c.z, c.x);
#endif
		struct coordinate *coords = malloc(sizeof(struct coordinate));
		coords[0] = c;
		rn->coords = coords;
	} else if (n_pins == 2) {
		struct coordinate a = npm->pins[net][0].coordinate;
		struct coordinate b = npm->pins[net][1].coordinate;
#ifdef DUMB_ROUTE_DEBUG
		printf("[dumb_route] net %d (%s):\n  ", net, get_net_name(blif, net));
#endif
		free(rn);
		rn = cityblock_route(net, a, b);
#ifdef DUMB_ROUTE_DEBUG
		print_routed_net(rn);
#endif
	} else {
		struct coordinate *coords = calloc(n_pins, sizeof(struct coordinate));
		for (int j = 0; j < n_pins; j++)
			coords[j] = npm->pins[net][j].coordinate;

		struct segments *mst = create_mst(coords, n_pins);
#ifdef DUMB_ROUTE_DEBUG
		printf("[dumb_route] net %d (%s) has %d pins:\n", net, get_net_name(blif, net), n_pins);
#endif
		for (int j = 0; j < mst->n_segments; j++) {
			struct routed_net *sub_rn = cityblock_route(net, mst->segments[j].start, mst->segments[j].end);
#ifdef DUMB_ROUTE_DEBUG
			printf("  segment %d: ", j);
			print_routed_net(sub_rn);
#endif
			rn->coords = realloc(rn->coords, sizeof(struct coordinate) * (rn->n_coords + sub_rn->n_coords));
			memcpy(&rn->coords[rn->n_coords], sub_rn->coords, sizeof(struct coordinate) * sub_rn->n_coords);
			rn->n_coords += sub_rn->n_coords;

			free_routed_net(sub_rn);
		}
		free_segments(mst);
	}

	return rn;
}

struct routings *route(struct blif *blif, struct cell_placements *cp)
{
	struct pin_placements *pp = placer_place_pins(cp);
	struct net_pin_map *npm = placer_create_net_pin_map(pp);

	struct routings *rt = malloc(sizeof(struct routings));
	rt->n_routed_nets = npm->n_nets;
	rt->routed_nets = calloc(rt->n_routed_nets + 1, sizeof(struct routed_net));

	for (net_t i = 1; i < cp->n_nets; i++)
		rt->routed_nets[i] = *dumb_route(blif, npm, i);

	free_pin_placements(pp);
	free_net_pin_map(npm);

	return rt;
}
