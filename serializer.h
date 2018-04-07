#ifndef __SERIALIZER_H__
#define __SERIALIZER_H__

#include "placer.h"
#include "blif.h"
#include "base_router.h"
#include "extract.h"

void serialize_placements(FILE *, struct cell_placements *, struct blif *);
void serialize_routings(FILE *, struct routings *, struct blif *);
void serialize_extraction(FILE *, struct extraction *);

#endif /* __SERIALIZER_H__ */
