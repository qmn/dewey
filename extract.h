#ifndef __EXTRACT_H__
#define __EXTRACT_H__

#include "coord.h"
#include "cell.h"
#include "router.h"
#include "placer.h"

struct extraction {
	struct dimensions dimensions;
	block_t *blocks;
	data_t *data;
};

// void recenter_placements(struct cell_placements *);
void recenter(struct cell_placements *, struct routings *);

struct extraction *extract_placements(struct cell_placements *);
struct extraction *extract(struct cell_placements *, struct routings *);
void free_extraction(struct extraction *);

#endif /* __EXTRACT_H__ */
