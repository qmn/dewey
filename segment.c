#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "segment.h"

int distance_pythagorean(struct coordinate a, struct coordinate b)
{
	int dx = abs(a.x - b.x);
	int dz = abs(a.z - b.z);
	return roundl(sqrt(dx * dx + dz * dz));
}

int distance_cityblock(struct coordinate a, struct coordinate b)
{
	return abs(a.x - b.x) + abs(a.z - b.z);
}

int distance_cityblock_pins(struct placed_pin *start, struct placed_pin *end)
{
	return abs(start->coordinate.x - end->coordinate.x) + abs(start->coordinate.z - end->coordinate.z);
}

struct mst_ubr_node *mst_find(struct mst_ubr_node *n)
{
	if (n->parent != n)
		n->parent = mst_find(n->parent);

	return n->parent;
}

void mst_union(struct mst_ubr_node *x, struct mst_ubr_node *y)
{
	struct mst_ubr_node *rx = mst_find(x);
	struct mst_ubr_node *ry = mst_find(y);

	if (rx == ry)
		return;

	rx->parent = ry;

	if (rx->rank == ry->rank)
		ry->rank++;
}

struct mst_ubr_node *mst_make_set(int n)
{
	struct mst_ubr_node *ubr = calloc(sizeof(struct mst_ubr_node), n);

	for (int i = 0; i < n; i++) {
		ubr[i].parent = &ubr[i];
		ubr[i].rank = 0;
		ubr[i].me = i;
	}

	return ubr;
}

static int mst_heapsort_cmp(const void *x, const void *y)
{
	return ((struct mst_heap_node *)x)->weight - ((struct mst_heap_node *)y)->weight;
}

static struct mst_heap *mst_heapsort(struct coordinate *locs, int n_locs)
{
	struct mst_heap *h = malloc(sizeof(struct mst_heap));
	h->n_nodes = n_locs * (n_locs - 1) / 2;
	h->nodes = calloc(h->n_nodes, sizeof(struct mst_heap_node));

	/* create weight matrix; x > y */
	int c = 0;
	for (int y = 0; y < n_locs; y++) {
		for (int x = y + 1; x < n_locs; x++) {
			struct mst_heap_node n = {distance_cityblock(locs[x], locs[y]), x, y};
			h->nodes[c++] = n;
		}
	}

	/* surprise! use qsort instead */
	qsort(h->nodes, h->n_nodes, sizeof(struct mst_heap_node), mst_heapsort_cmp);

	return h;
}

struct segments *create_mst(struct coordinate *locs, int n_locs)
{
	struct segments *mst = malloc(sizeof(struct segments));
	mst->n_segments = n_locs - 1;
	mst->segments = calloc(sizeof(struct segment), (mst->n_segments));

	struct mst_heap *h = mst_heapsort(locs, n_locs);

	struct mst_ubr_node *s = calloc(sizeof(struct mst_ubr_node), n_locs);
	// make-set
	for (int i = 0; i < n_locs; i++) {
		s[i].parent = &s[i];
		s[i].rank = 0;
		s[i].me = i;
	}

	int n_segments = 0;
	for (int i = 0; i < h->n_nodes; i++) {
		int x = h->nodes[i].x;
		int y = h->nodes[i].y;
		if (mst_find(&s[x]) != mst_find(&s[y])) {
			struct segment new_seg = {locs[x], locs[y]};
			mst->segments[n_segments++] = new_seg;
			mst_union(&s[x], &s[y]);
			if (n_segments >= mst->n_segments)
				break;
		}
	}

	return mst;
}

void free_segments(struct segments *s)
{
	free(s->segments);
	free(s);
}
