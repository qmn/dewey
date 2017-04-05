#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "extract.h"
#include "blif.h"
#include "placer.h"
#include "router.h"
#include "cell.h"
#include "vis_png.h"

int main(int argc, char **argv)
{
	if (argc < 2) {
		/* print help text */
		printf("dewey -- a placer and router tool for Minecraft redstone circuits\n");
		printf("\n");
		printf("Usage: %s <input BLIF file>\n", argv[0]);
		return 1;
	}

	/* read input blif */
	FILE *blif_file;
	struct blif *blif;

	blif_file = fopen(argv[1], "r");

        if (!blif_file) {
                printf("[dewey] could not read %s: %s\n", argv[1], strerror(errno));
                return 2;
        }

	printf("[dewey] reading blif...\n");
	blif = read_blif(blif_file);
	fclose(blif_file);
	printf("[dewey] done reading blif\n");

	/* read cell library */
	FILE *cell_library_file;
	struct cell_library *cl;

	cell_library_file = fopen("quan.yaml", "rb");

	if (!cell_library_file) {
		printf("[dewey] could not read cell library file: %s\n", strerror(errno));
		return 3;
	}

	cl = read_cell_library(cell_library_file);
	fclose(cell_library_file);

	/* begin with initial placement */
	struct cell_placements *initial_placement = placer_initial_place(blif, cl);
	struct dimensions initial_dimensions = compute_placement_dimensions(initial_placement);
	print_cell_placements(initial_placement);

	printf("[dewey] intial dimensions: {x: %d, y: %d, z: %d}\n",
		initial_dimensions.x, initial_dimensions.y, initial_dimensions.z);

	struct cell_placements *new_placements = simulated_annealing_placement(initial_placement, &initial_dimensions, 100, 100, 100);
	// struct cell_placements *new_placements = initial_placement;
	// print_cell_placements(new_placements);

	struct dimensions placement_dimensions = compute_placement_dimensions(new_placements);
	printf("[dewey] placement dimensions: {x: %d, y: %d, z: %d}\n",
		placement_dimensions.x, placement_dimensions.y, placement_dimensions.z);

/*
	struct dimensions fd = compute_placement_dimensions(new_placements);
	block_t *flattened = extract_placements(new_placements);
	for (int y = 0; y < fd.y; y++) {
		for (int z = 0; z < fd.z; z++) {
			for (int x = 0; x < fd.x; x++) {
				printf("%3u ", flattened[y * fd.x * fd.z + z * fd.x + x]);
			}
			printf("\n");
		}
		printf("\n");
	}
	free(flattened);
*/

	struct routings *routings = route(blif, new_placements);

	recenter(new_placements, routings, 0);

	vis_png_draw_placements(blif, new_placements, routings);

        free_blif(blif);
	free_cell_library(cl);

	return 0;
}
