#include <stdio.h>
#include <stdlib.h>

#include "placer.h"
#include "router.h"

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

struct routing_sets *router_generate_routing_sets(struct cell_placements *cp)
{
	struct routing_sets *rs = malloc(sizeof(struct routing_sets));
	rs->n_sets = cp->n_nets;
	rs->sets = calloc(cp->n_nets, sizeof(struct routing_set));



	return rs;
}

void router_free_routing_sets(struct routing_sets *rs)
{
	free(rs->sets);
	free(rs);
}
