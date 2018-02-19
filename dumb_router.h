#ifndef __DUMB_ROUTER_H__
#define __DUMB_ROUTER_H__

#include "blif.h"
#include "base_router.h"

void dumb_route(struct routed_net *, struct blif *, struct net_pin_map *, net_t);

#endif /* __DUMB_ROUTER_H__ */
