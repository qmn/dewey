#ifndef __SEGMENT_H__
#define __SEGMENT_H__

#include "cell.h"

struct segment {
	struct placed_pin *start;
	struct placed_pin *end;
};

struct segments {
	struct segment *segments;
	int n_segments;
};

struct mst_heap {
	struct mst_heap_node *nodes;
	int n_nodes;
};

struct mst_heap_node {
	int weight;
	int x; // indices into the array of locs
	int y;
};

struct mst_ubr_node {
	struct mst_ubr_node *parent;
	int rank;
};

int distance_pythagorean(struct coordinate, struct coordinate);
int distance_cityblock(struct coordinate, struct coordinate);

struct segments *create_mst(struct placed_pin *, int);
void free_segments(struct segments *);

#endif /* __SEGMENT_H__ */
