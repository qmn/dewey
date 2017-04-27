#ifndef __ROUTER_H__
#define __ROUTER_H__

#include "coord.h"
#include "segment.h"
#include "blif.h"
#include "placer.h"

struct cost_coord {
	unsigned int cost;
	struct coordinate coord;
};

struct routings {
	int n_routed_nets;
	struct routed_net *routed_nets;

	struct net_pin_map *npm;
};

struct routed_segment {
	struct segment seg;

	int n_coords;
	struct coordinate *coords;

	int score;

	struct routed_net *net;
};

struct routed_net {
	net_t net;

	/* pins connected by this net;
	 * pins are references to a struct placed_pins
	 */
	int n_pins;
	struct placed_pin *pins;

	int n_routed_segments;
	struct routed_segment *routed_segments;
};

struct routings *route(struct blif *, struct cell_placements *);
struct routings *copy_routings(struct routings *);
struct dimensions compute_routings_dimensions(struct routings *);
void free_routings(struct routings *);

void routings_displace(struct routings *, struct coordinate);

#endif /* __ROUTER_H__ */
