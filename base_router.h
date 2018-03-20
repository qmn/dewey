#ifndef __BASE_ROUTER_H__
#define __BASE_ROUTER_H__

#include "blif.h"
#include "coord.h"
#include "segment.h"

// base_router.h contains things needed by all router implementations,
// but to avoid circular include chains (particularly of router.h)

// west: x-1, east: x+1, north: z-1, south: z+1
enum backtrace {
	BT_NONE  = 0,
	BT_WEST  = 1 << 0,
	BT_SOUTH = 1 << 1,
	BT_EAST  = 1 << 2,
	BT_NORTH = 1 << 3,
	BT_DOWN  = 1 << 4,
	BT_UP    = 1 << 5,
	BT_START = 1 << 6
};
#define BT_VERTICAL_MASK (BT_UP | BT_DOWN)
#define BT_CARDINAL_MASK (BT_WEST | BT_SOUTH | BT_EAST | BT_NORTH)
#define is_vertical(bt) (bt & BT_VERTICAL_MASK)
#define is_cardinal(bt) (bt & BT_CARDINAL_MASK)

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

	int extracted;
};

struct routed_segment_head {
	struct routed_segment_head *next;
	struct routed_segment rseg;
};

// a routed_segment_adjacency is a linked list of
// adjacencies between segments and other segments or pins
// b_type determines what's contained in the union
enum rsa_type { SEGMENT, PIN };
struct routed_segment_adjacency {
	struct routed_segment_adjacency *next;

	struct routed_segment *parent;
	enum rsa_type child_type;
	union {
		struct routed_segment *rseg;
		struct placed_pin *pin;
	} child;

	struct coordinate at;
};

struct routed_net {
	net_t net;

	/* pins connected by this net;
	 * pins are references to a struct placed_pins
	 */
	int n_pins;
	struct placed_pin *pins;

	struct routed_segment_head *routed_segments;

	// adjacency list expressing connections between
	// routed_segments and other segments or pins
	struct routed_segment_adjacency *adjacencies;
};

struct coordinate disp_backtrace(struct coordinate, enum backtrace);
enum backtrace invert_backtrace(enum backtrace);
void invert_backtrace_sequence(enum backtrace *, int);
enum backtrace compute_backtrace(struct coordinate, struct coordinate);

void add_adjacent_segment(struct routed_net *, struct routed_segment *, struct routed_segment *, struct coordinate);
void add_adjacent_pin(struct routed_net *, struct routed_segment *, struct placed_pin *);

struct routed_segment *find_parent(struct routed_net *rn, struct placed_pin *p, struct routed_segment *rseg);
#endif /* __BASE_ROUTER_H__ */
