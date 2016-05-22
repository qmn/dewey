#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "blif.h"
#include "placer.h"

int main(int argc, char **argv)
{
	if (argc < 2) {
		/* print help text */
		printf("dewey -- a placer and router tool for Minecraft redstone circuits\n");
		printf("\n");
		printf("Usage: %s <input BLIF file>\n", argv[0]);
		return 1;
	}

	FILE *blif_file;
	struct blif *blif;

	/* read input blif */
	blif_file = fopen(argv[1], "r");

        if (!blif_file) {
                printf("[dewey] could not read %s: %s\n", argv[1], strerror(errno));
                return 2;
        }

	blif = read_blif(blif_file);
	fclose(blif_file);

        free_blif(blif);

	return 0;
}
