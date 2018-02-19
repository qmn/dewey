#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "heap.h"
#include "base_router.h"

struct cost_coord_heap *create_cost_coord_heap()
{
	struct cost_coord_heap *h = malloc(sizeof(struct cost_coord_heap));
	h->elts = NULL;
	assert(h);

	h->size = 1;
	h->n_elts = 0;
	h->elts = malloc(h->size * sizeof(struct cost_coord));
	assert(h->elts);

	return h;
}

void free_cost_coord_heap(struct cost_coord_heap *h)
{
	free(h->elts);
	free(h);
}

static void swap(struct cost_coord_heap *h, unsigned int i, unsigned int j)
{
	struct cost_coord tmp = h->elts[i];
	h->elts[i] = h->elts[j];
	h->elts[j] = tmp;
}

static void min_heapify(struct cost_coord_heap *h, unsigned int i)
{
	unsigned int left = 2 * i;
	unsigned int right = 2 * i + 1;
	unsigned int smallest = i;

	if (left <= h->n_elts && h->elts[left].cost < h->elts[smallest].cost) {
		smallest = left;
	}

	if (right <= h->n_elts && h->elts[right].cost < h->elts[smallest].cost) {
		smallest = right;
	}

	if (smallest != i) {
		swap(h, i, smallest);
		min_heapify(h, smallest);
	}
}

void cost_coord_heap_insert(struct cost_coord_heap *h, struct cost_coord c)
{
	unsigned int i = h->n_elts + 1;
	if (i >= h->size) {
		h->size *= 2;
		h->elts = realloc(h->elts, h->size * sizeof(struct cost_coord));
#ifdef COST_COORD_HEAP_DEBUG
		printf("cost_coord_heap doubling in size to %d\n", h->size + 1);
#endif
		assert(h->elts);
	}

	h->elts[i] = c;
	h->n_elts++;

	int p = i / 2;
	while (h->elts[i].cost < h->elts[p].cost && p != 0) {
		swap(h, i, p);
		i = p;
		p = i / 2;
	}
}

struct cost_coord cost_coord_heap_delete_min(struct cost_coord_heap *h)
{
	assert(h->n_elts > 0);
	struct cost_coord elt = h->elts[1];

	/* replace the first entry with the last entry */
	h->elts[1] = h->elts[h->n_elts--];

	/* re-run min-heap on the remaining */
	min_heapify(h, 1);

	return elt;
}

struct cost_coord cost_coord_heap_peek(struct cost_coord_heap *h)
{
	assert(h->n_elts > 0);
	return h->elts[1];
}

/*
static int cost_coord_equals(struct cost_coord a, struct cost_coord b)
{
	return coordinate_equal(a.coord, b.coord) && a.cost == b.cost;
}
*/

int cost_coord_heap_contains_coordinate(struct cost_coord_heap *h, struct coordinate c)
{
	int i;

	for (i = 1; i < h->n_elts + 1; i++)
		if (coordinate_equal(h->elts[i].coord, c))
			return 1;

	return 0;
}

void clear_cost_coord_heap(struct cost_coord_heap *h)
{
	h->n_elts = 0;
}
