#include "base_router.h"

#include <assert.h>
#include <stdlib.h>

#include "placer.h"

struct coordinate disp_backtrace(struct coordinate c, enum backtrace b)
{
	switch (b) {
	case BT_WEST:
		c.x--; break;
	case BT_EAST:
		c.x++; break;
	case BT_NORTH:
		c.z--; break;
	case BT_SOUTH:
		c.z++; break;
	case BT_UP:
		c.y += 3; break;
	case BT_DOWN:
		c.y -= 3; break;
	default:
		break;
	}

	return c;
}

enum backtrace invert_backtrace(enum backtrace b)
{
	switch (b) {
	case BT_WEST:
		return BT_EAST;
	case BT_EAST:
		return BT_WEST;
	case BT_NORTH:
		return BT_SOUTH;
	case BT_SOUTH:
		return BT_NORTH;
	case BT_UP:
		return BT_DOWN;
	case BT_DOWN:
		return BT_UP;
	default:
		return BT_NONE;
	}
}

// not only reverse the backtrace order but also invert the backtrace direction
// example: moving north, north, west, end:
// E<A
//   A
//   S
// becomes: moving east, south, south, end:
// S>V
//   V
//   E
void invert_backtrace_sequence(enum backtrace *bt, int n_bt)
{
	int i;
	enum backtrace tmp;

	for (i = 0; i < n_bt / 2; i++) {
		tmp = bt[i];
		bt[i] = bt[n_bt - 1 - i];
		bt[n_bt - i - 1] = tmp;
	}

	for (i = 0; i < n_bt; i++)
		bt[i] = invert_backtrace(bt[i]);
}

// figure out to get from here to there
enum backtrace compute_backtrace(struct coordinate here, struct coordinate there)
{
	int dx = there.x - here.x;
	int dz = there.z - here.z;
	int dy = there.y - here.y;

	// look away, linus
	if (dy == 0) {
		if (dx == 0) {
			if (dz == 1)
				return BT_SOUTH;
			else if (dz == -1)
				return BT_NORTH;
		} else if (dx == 1) {
			return BT_EAST;
		} else if (dx == -1) {
			return BT_WEST;
		}
	} else if (dy == 3) {
		return BT_UP;
	} else if (dy == -3) {
		return BT_DOWN;
	}

	return BT_NONE;
}

void add_adjacent_segment(struct routed_net *rn, struct routed_segment *sega, struct routed_segment *segb)
{
	assert(sega != segb);

	struct routed_segment_adjacency *rsa = malloc(sizeof(struct routed_segment_adjacency));
	rsa->parent = sega;
	rsa->child_type = SEGMENT;
	rsa->child.rseg = segb;

	rsa->next = rn->adjacencies;
	rn->adjacencies = rsa;
}

void add_adjacent_pin(struct routed_net *rn, struct routed_segment *seg, struct placed_pin *pin)
{
	struct routed_segment_adjacency *rsa = malloc(sizeof(struct routed_segment_adjacency));
	rsa->parent = seg;
	rsa->child_type = PIN;
	rsa->child.pin = pin;

	rsa->next = rn->adjacencies;
	rn->adjacencies = rsa;
}

// find the parent of the pin, or the routed segment
// for a pin, set rseg to NULL -- returns NULL if no parent, or the segment if there is
// for a segment, set p to NULL -- if there is no parent, it returns the rseg passed in
struct routed_segment *find_parent(struct routed_net *rn, struct placed_pin *p, struct routed_segment *rseg)
{
	for (struct routed_segment_adjacency *rsa = rn->adjacencies; rsa; rsa = rsa->next) {
		if ((rseg && rsa->child_type == SEGMENT && rsa->child.rseg == rseg) ||
		    (!rseg && rsa->child_type == PIN && rsa->child.pin == p))
			return find_parent(rn, p, rsa->parent);
	}

	return rseg;
}
