#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "blif.h"
#include "placer.h"
#include "cell.h"

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

	blif = read_blif(blif_file);
	fclose(blif_file);

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
	print_cell_placements(initial_placement);

        free_blif(blif);
	free_cell_library(cl);

	return 0;
}
