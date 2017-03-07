#ifndef __ROUTER_H__
#define __ROUTER_H__

/* one for each pin in a net */
struct routing_pin {
	int cell_index;
	char *name;
	struct coordinate coord;
	enum ordinal_direction facing;
	int is_output;
};

/* one for each net */
struct routing_set {
	char *name;

	int n_pins;
	struct routing_pin *pins;
};

/* one for a design */
struct routing_sets {
	int n_sets;
	struct routing_set *sets;
};

struct routed_net {
	net_t net;
	int n_coords;
	struct coordinate *coords;
};

#endif /* __ROUTER_H__ */
