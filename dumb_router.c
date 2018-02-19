#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "router.h"
#include "base_router.h"

// route, based on a cityblock algorithm, without regard to Y or obstacles
static struct routed_segment cityblock_route(struct segment seg)
{
	struct coordinate a = seg.start;
	struct coordinate b = seg.end;

	int len = distance_cityblock(a, b); // does not include BT_START
	int count = 0;
	enum backtrace *path = malloc(sizeof(enum backtrace) * len);

	assert(a.y == b.y);

	struct coordinate c = b;
	while (c.x != a.x || c.z != a.z) {
		if (c.x > a.x) {
			c.x--;
			path[count++] = BT_WEST;
		} else if (c.x < a.x) {
			c.x++;
			path[count++] = BT_EAST;
		} else if (c.z > a.z) {
			c.z--;
			path[count++] = BT_NORTH;
		} else if (c.z < a.z) {
			c.z++;
			path[count++] = BT_SOUTH;
		}
	}

	assert(count == len);

	struct routed_segment rseg = {seg, len, path, 0, NULL, NULL, 0, NULL, 0, NULL};
	return rseg;
}


struct ubr_node {
	struct mst_node *x;
	struct mst_node *y;
	int score;
};

static int ubr_node_cmp(const void *a, const void *b)
{
	struct ubr_node *aa = (struct ubr_node *)a, *bb = (struct ubr_node *)b;
	return aa->score - bb->score;
}

struct mst_node {
	struct mst_node *parent;
	struct placed_pin *pin;
	int rank;
};

static struct mst_node *dumb_mst_find(struct mst_node *x)
{
	if (x->parent != x)
		x->parent = dumb_mst_find(x->parent);

	return x->parent;
}

static void dumb_mst_union(struct mst_node *x, struct mst_node *y)
{
	struct mst_node *rx = dumb_mst_find(x);
	struct mst_node *ry = dumb_mst_find(y);

	if (rx == ry)
		return;

	rx->parent = ry;

	if (rx->rank == ry->rank)
		ry->rank++;
}

struct routed_segment *find_parent_rseg(struct placed_pin *p)
{
	struct routed_segment *rseg = p->parent;
	if (!rseg)
		return NULL;

	while (rseg->parent != NULL)
		rseg = rseg->parent;

	return rseg;
}

void add_child_pin(struct routed_segment *rseg, struct placed_pin *p)
{
	rseg->child_pins = realloc(rseg->child_pins, sizeof(struct placed_pin *) * ++rseg->n_child_pins);
	rseg->child_pins[rseg->n_child_pins - 1] = p;
	p->parent = rseg;
}

void add_child_segment(struct routed_segment *parent_rseg, struct routed_segment *rseg)
{
	parent_rseg->child_segments = realloc(parent_rseg->child_segments, sizeof(struct routed_segment *) * ++parent_rseg->n_child_segments);
	parent_rseg->child_segments[parent_rseg->n_child_segments - 1] = rseg;
	rseg->parent = parent_rseg;
}


// create routed segments for this net blindly (that is, without regard
// to other objects) using cityblock routing
void dumb_mst_route(struct routed_net *rn)
{
	// generate MST nodes
	struct mst_node *mst_nodes = calloc(rn->n_pins, sizeof(struct mst_node));
	for (int i = 0; i < rn->n_pins; i++)
		mst_nodes[i] = (struct mst_node){&mst_nodes[i], &rn->pins[i], 0};

	// generate union-by-rank nodes
	int n_ubr_nodes = rn->n_pins * (rn->n_pins - 1) / 2;
	struct ubr_node *ubr_nodes = calloc(n_ubr_nodes, sizeof(struct ubr_node));
	int k = 0;
	for (int i = 0; i < rn->n_pins; i++)
		for (int j = i + 1; j < rn->n_pins; j++)
			ubr_nodes[k++] = (struct ubr_node){&mst_nodes[i], &mst_nodes[j], distance_cityblock(extend_pin(mst_nodes[i].pin), extend_pin(mst_nodes[j].pin))};

	qsort(ubr_nodes, n_ubr_nodes, sizeof(struct ubr_node), ubr_node_cmp);

	// search through the list of union-by-rank nodes
	int count = 0;
	for (int i = 0; i < n_ubr_nodes && count < rn->n_pins - 1; i++) {
		struct ubr_node a = ubr_nodes[i];
		struct mst_node *x = a.x;
		struct mst_node *y = a.y;

		// if this node has a different parent than me, create
		// a new segment connecting these two parents.
		if (dumb_mst_find(x) != dumb_mst_find(y)) {
			struct segment seg = {extend_pin(x->pin), extend_pin(y->pin)};
			struct routed_segment_head *rsh = malloc(sizeof(struct routed_segment_head));
			rsh->next = NULL;
			rsh->rseg = cityblock_route(seg);
			rsh->rseg.net = rn;
			
			// if either object being joined has a parent segment already,
			// set that parent segment as a child of this newly-formed
			// segment

			struct routed_segment *xp_rseg = find_parent_rseg(x->pin);
			if (xp_rseg)
				add_child_segment(&rsh->rseg, xp_rseg);
			else
				add_child_pin(&rsh->rseg, x->pin);

			struct routed_segment *yp_rseg = find_parent_rseg(y->pin);
			if (yp_rseg)
				add_child_segment(&rsh->rseg, yp_rseg);
			else
				add_child_pin(&rsh->rseg, y->pin);

			routed_net_add_segment_node(rn, rsh);
			count++;
			dumb_mst_union(x, y);
		}
	}
}

/* generate the MST for this net to determine the order of connections,
 * then connect them all with a city*/
void dumb_route(struct routed_net *rn, struct blif *blif, struct net_pin_map *npm, net_t net)
{
	int n_pins = npm->n_pins_for_net[net];
	assert(n_pins > 0);

	rn->net = net;
	rn->routed_segments = NULL;

	rn->n_pins = n_pins;
	rn->pins = malloc(sizeof(struct placed_pin) * rn->n_pins);
	memcpy(rn->pins, npm->pins[net], sizeof(struct placed_pin) * rn->n_pins);

	if (n_pins == 1) {
		struct segment seg = {npm->pins[net][0].coordinate, npm->pins[net][0].coordinate};
		enum backtrace *path = malloc(sizeof(enum backtrace));
		path[0] = BT_START;

		struct placed_pin **child_pins = malloc(sizeof(struct placed_pin *));
		child_pins[0] = &npm->pins[net][0];

		struct routed_segment_head *rsh = malloc(sizeof(struct routed_segment_head));
		rsh->rseg = (struct routed_segment){seg, 1, path, 0, rn, NULL, 0, NULL, 1, child_pins};
		rsh->next = NULL;
		child_pins[0]->parent = &rsh->rseg;
		
		rn->routed_segments = rsh;
	} else if (n_pins == 2) {
		struct segment seg = {extend_pin(&npm->pins[net][0]), extend_pin(&npm->pins[net][1])};

		struct placed_pin **child_pins = malloc(2 * sizeof(struct placed_pin *));
		child_pins[0] = &npm->pins[net][0];
		child_pins[1] = &npm->pins[net][1];

		struct routed_segment_head *rsh = malloc(sizeof(struct routed_segment_head));
		rsh->rseg = cityblock_route(seg);
		rsh->rseg.net = rn;
		rsh->rseg.n_child_pins = 2;
		rsh->rseg.child_pins = child_pins;
		child_pins[0]->parent = &rsh->rseg;
		child_pins[1]->parent = &rsh->rseg;
		rsh->next = NULL;

		rn->routed_segments = rsh;
	} else {
		dumb_mst_route(rn);
	}
}
