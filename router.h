#ifndef __ROUTER_H__
#define __ROUTER_H__

#include "coord.h"
#include "blif.h"
#include "placer.h"

struct routings {
	int n_routed_nets;
	struct routed_net *routed_nets;
};

struct routed_net {
	net_t net;
	int n_coords;
	struct coordinate *coords;
};

struct routings *route(struct blif *, struct cell_placements *);
struct routings *copy_routings(struct routings *);
struct dimensions compute_routings_dimensions(struct routings *);
void free_routings(struct routings *);

void routings_displace(struct routings *, struct coordinate);

#endif /* __ROUTER_H__ */
