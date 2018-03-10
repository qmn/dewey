#include "base_router.h"

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
