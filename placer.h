#ifndef __PLACER_H__
#define __PLACER_H__
struct coordinate {
	int y;
	int z;
	int x;
};

struct dimensions {
	unsigned int y;
	unsigned int z;
	unsigned int x;
};

struct pin {
	char *name;
	char *net;
};

struct placement {
	char *name;
	struct coordinate placement;
	unsigned long turns;
	unsigned long n_pins;
	struct pin **pins;
};

struct cell_placements {
	struct placement **placements;
	unsigned long length;
};
#endif /* __PLACER_H__ */
