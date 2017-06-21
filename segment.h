#ifndef __SEGMENT_H__
#define __SEGMENT_H__

#include "cell.h"

struct segment {
	struct coordinate start;
	struct coordinate end;
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
	int me;
};

int distance_pythagorean(struct coordinate, struct coordinate);
int distance_cityblock(struct coordinate, struct coordinate);

struct mst_ubr_node *mst_make_set(int);
struct mst_ubr_node *mst_find(struct mst_ubr_node *);
void mst_union(struct mst_ubr_node *, struct mst_ubr_node *);

struct segments *create_mst(struct coordinate *, int);
void free_segments(struct segments *);

#endif /* __SEGMENT_H__ */
