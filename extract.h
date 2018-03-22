#ifndef __EXTRACT_H__
#define __EXTRACT_H__

#include "coord.h"
#include "cell.h"
#include "base_router.h"
#include "placer.h"

struct extraction {
	struct dimensions dimensions;
	block_t *blocks;
	data_t *data;
};

struct coordinate placements_top_left_most_point(struct cell_placements *);
struct coordinate routings_top_left_most_point(struct routings *);

void recenter(struct cell_placements *, struct routings *, int);

struct extraction *extract_placements(struct cell_placements *);
struct extraction *extract(struct cell_placements *, struct routings *);
void free_extraction(struct extraction *);

#endif /* __EXTRACT_H__ */
