#ifndef __VIS_PNG_H__
#define __VIS_PNG_H__

#include "placer.h"
#include "router.h"

void vis_png_draw_placements(struct cell_placements *, struct routings *);

unsigned char *flatten(struct cell_placements *);

#endif /* __VIS_PNG_H__ */
