#ifndef __ROUTER_H__
#define __ROUTER_H__

#include "blif.h"
#include "placer.h"
#include "base_router.h"

struct routings *route(struct blif *, struct cell_placements *);
struct routings *copy_routings(struct routings *);
struct dimensions compute_routings_dimensions(struct routings *);
void free_routings(struct routings *);

int segment_routed(struct routed_segment *);
void routed_net_add_segment_node(struct routed_net *, struct routed_segment_head *);

void routings_displace(struct routings *, struct coordinate);

#endif /* __ROUTER_H__ */
