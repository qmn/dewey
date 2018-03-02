#ifndef __BASE_ROUTER_H__
#define __BASE_ROUTER_H__

#include "blif.h"
#include "coord.h"
#include "segment.h"

// base_router.h contains things needed by all router implementations,
// but to avoid circular include chains (particularly of router.h)

// west: x-1, east: x+1, north: z-1, south: z+1
enum backtrace {BT_NONE, BT_WEST, BT_SOUTH, BT_EAST, BT_NORTH, BT_DOWN, BT_UP, BT_START};
#define is_vertical(bt) (bt == BT_UP || bt == BT_DOWN)
#define is_cardinal(bt) (bt == BT_WEST || bt == BT_SOUTH || bt == BT_EAST || bt == BT_NORTH)

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

	int n_backtraces;
	enum backtrace *bt;

	int score;

	struct routed_net *net;

	struct routed_segment *parent;
	int n_child_segments;
	struct routed_segment **child_segments;
	int n_child_pins;
	struct placed_pin **child_pins;
};

struct routed_segment_head {
	struct routed_segment_head *next;
	struct routed_segment rseg;
};

struct routed_net {
	net_t net;

	/* pins connected by this net;
	 * pins are references to a struct placed_pins
	 */
	int n_pins;
	struct placed_pin *pins;

	struct routed_segment_head *routed_segments;
};

struct coordinate disp_backtrace(struct coordinate, enum backtrace);
enum backtrace invert_backtrace(enum backtrace);
enum backtrace compute_backtrace(struct coordinate, struct coordinate);
#endif /* __BASE_ROUTER_H__ */
