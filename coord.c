#include "coord.h"

static int min(int a, int b)
{
	return a < b ? a : b;
}

static int max(int a, int b)
{
	return a > b ? a : b;
}

struct coordinate coordinate_neg(struct coordinate a)
{
	struct coordinate c = {-a.y, -a.z, -a.x};
	return c;
}

struct coordinate coordinate_add(struct coordinate a, struct coordinate b)
{
	struct coordinate c = {a.y + b.y, a.z + b.z, a.x + b.x};
	return c;
}

struct coordinate coordinate_sub(struct coordinate a, struct coordinate b)
{
	struct coordinate c = {a.y - b.y, a.z - b.z, a.x - b.x};
	return c;
}

struct coordinate coordinate_piecewise_min(struct coordinate a, struct coordinate b)
{
	struct coordinate c = {min(a.y, b.y), min(a.z, b.z), min(a.x, b.x)};
	return c;
}

struct coordinate coordinate_piecewise_max(struct coordinate a, struct coordinate b)
{
	struct coordinate c = {max(a.y, b.y), max(a.z, b.z), max(a.x, b.x)};
	return c;
}

struct dimensions dimensions_piecewise_max(struct dimensions a, struct dimensions b)
{
	struct dimensions d = {max(a.y, b.y), max(a.z, b.z), max(a.x, b.x)};
	return d;
}
