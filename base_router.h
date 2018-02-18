#ifndef __BASE_ROUTER_H__
#define __BASE_ROUTER_H__

// base_router.h contains things needed by all router implementations,
// but to avoid circular include chains (particularly of router.h)

#include "coord.h"

enum backtrace {BT_NONE, BT_WEST, BT_SOUTH, BT_EAST, BT_NORTH, BT_DOWN, BT_UP, BT_START};
enum backtrace backtraces[] = {BT_WEST, BT_NORTH, BT_EAST, BT_SOUTH, BT_DOWN, BT_UP};
struct coordinate movement_offsets[] = {{0, 0, 1}, {0, 1, 0}, {0, 0, -1}, {0, -1, 0}, {3, 0, 0}, {-3, 0, 0}};
int n_movements = sizeof(movement_offsets) / sizeof(struct coordinate);
#define is_vertical(bt) (bt == BT_UP || bt == BT_DOWN)

#endif /* __BASE_ROUTER_H__ */
