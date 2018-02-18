#ifndef __USAGE_MATRIX_H__
#define __USAGE_MATRIX_H__

#include "coord.h"
#include "extract.h"
#include "placer.h"
#include "router.h"

#define USAGE_SIZE(m) (m->d.x * m->d.y * m->d.z)

struct usage_matrix {
	struct dimensions d;
	int xz_margin;
	unsigned char *matrix;
};

int in_usage_bounds(struct usage_matrix *, struct coordinate);

int usage_idx(struct usage_matrix *, struct coordinate);
void usage_mark(struct usage_matrix *, struct coordinate);

struct usage_matrix *create_usage_matrix(struct cell_placements *, struct routings *, int);

int usage_matrix_violated(struct usage_matrix *, struct coordinate);

#endif /* __USAGE_MATRIX_H__ */
