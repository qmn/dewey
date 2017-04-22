#ifndef __HEAP_H__
#define __HEAP_H_

#include "coord.h"

struct cost_coord_heap {
	unsigned int size;
	unsigned int n_elts;
	struct cost_coord *elts;
};

struct cost_coord_heap *create_cost_coord_heap();
void free_cost_coord_heap(struct cost_coord_heap *);
void clear_cost_coord_heap(struct cost_coord_heap *);

void cost_coord_heap_insert(struct cost_coord_heap *, struct cost_coord);
struct cost_coord cost_coord_heap_delete_min(struct cost_coord_heap *);

int cost_coord_heap_contains_coordinate(struct cost_coord_heap *, struct coordinate);

#endif /* __HEAP_H__ */
