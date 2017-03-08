#ifndef __ROUTER_H__
#define __ROUTER_H__

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

#endif /* __ROUTER_H__ */
