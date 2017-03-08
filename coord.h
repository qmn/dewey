#ifndef __COORD_H__
#define __COORD_H__

struct coordinate {
	int y;
	int z;
	int x;
};

struct dimensions {
	unsigned int y;
	unsigned int z;
	unsigned int x;
};

struct coordinate coordinate_neg(struct coordinate);
struct coordinate coordinate_add(struct coordinate, struct coordinate);
struct coordinate coordinate_sub(struct coordinate, struct coordinate);

struct coordinate coordinate_piecewise_min(struct coordinate, struct coordinate);
struct coordinate coordinate_piecewise_max(struct coordinate, struct coordinate);

struct dimensions dimensions_piecewise_max(struct dimensions, struct dimensions);

#endif /* __COORD_H__ */
