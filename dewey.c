#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/param.h>
#include <sys/stat.h>

#include "blif.h"
#include "cell.h"
#include "extract.h"
#include "placer.h"
#include "router.h"
#include "vis_png.h"
#include "vis_json.h"
#include "serializer.h"

void usage(char *argv0)
{
	printf("dewey -- a placer and router tool for Minecraft redstone circuits\n");
	printf("\n");
	printf("Usage: %s [options] <input BLIF file>\n", argv0);
	printf("Options:\n");
//	printf("  -l, --library=<yaml>       Cell library YAML file\n");
	printf("  -o, --output=<dir>         Directory to place output files\n");
	printf("  -s, --seed=<number>        Seed the random number generator\n");
}

int main(int argc, char **argv)
{
	// argv0 to remember binary name
	char *argv0 = argv[0];

	// input blif (required)
	char *input_blif = NULL;

	// output directory
	char output_dir[MAXPATHLEN];
	getcwd(output_dir, MAXPATHLEN);

	// random seed
	int seed = 0;

	// process long options
	static struct option longopts[] = {
		// {"library", optional_argument, NULL, 'l'},
		{"output" , optional_argument, NULL, 'o'},
		{"seed"   , optional_argument, NULL, 's'},
		{NULL,                      0, NULL,   0}
	};

	char c;
	int options_supplied = 0;
	while ((c = getopt_long(argc, argv, "i:o:s:", longopts, NULL)) != -1) {
		switch (c) {
		case 'o':
			realpath(optarg, output_dir);
			break;
		case 's':
			seed = atoi(optarg);
			break;
		default:
			usage(argv0);
			return 1;
		}
		argc -= optind;
		argv += optind;
		options_supplied++;
	}

	if (!options_supplied && argc == 2)
		input_blif = argv[1];
	else if (options_supplied && argc == 1)
		input_blif = argv[0];

	// process output dir
	strncat(output_dir, "/", MAXPATHLEN-1);
	printf("output dir is %s\n", output_dir);
	struct stat sb;
	if (stat(output_dir, &sb) != 0) {
		printf("[dewey] making output dir %s\n", output_dir);
		mkdir(output_dir, 0755);
	} else if (!S_ISDIR(sb.st_mode)) {
		printf("[dewey] output dir %s exists but isn't a directory, exiting\n", output_dir);
		return 1;
	}

	if (!input_blif) {
		usage(argv0);
		return 1;
	}

	/* read input blif */
	FILE *blif_file;
	struct blif *blif;

	blif_file = fopen(input_blif, "r");

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

	char *cl_fn = "quan.yaml";
	cell_library_file = fopen(cl_fn, "rb");

	if (!cell_library_file) {
		printf("[dewey] could not read cell library file: %s\n", strerror(errno));
		return 3;
	}

	cl = read_cell_library(cell_library_file, cl_fn);
	fclose(cell_library_file);

	printf("[dewey] seeding random generator to %d\n", seed);
	srandom(seed);

	/* begin with initial placement */
	struct cell_placements *initial_placement = placer_initial_place(blif, cl);
	struct dimensions initial_dimensions = compute_placement_dimensions(initial_placement);
	print_cell_placements(initial_placement);

	printf("[dewey] intial dimensions: {x: %d, y: %d, z: %d}\n",
		initial_dimensions.x, initial_dimensions.y, initial_dimensions.z);

	struct cell_placements *new_placements;
	new_placements = simulated_annealing_placement(initial_placement, &initial_dimensions, 100, 100, 100);
	// struct cell_placements *new_placements = initial_placement;
	// print_cell_placements(new_placements);

	vis_png_draw_placements(output_dir, blif, new_placements, NULL, 0);

	FILE *pf = fopen("placements.yaml", "w");
	serialize_placements(pf, new_placements, blif);
	fclose(pf);

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

	FILE *rf = fopen("routings.yaml", "w");
	serialize_routings(rf, routings, blif);
	fclose(rf);

	vis_png_draw_placements(output_dir, blif, new_placements, routings, 0);
	vis_png_draw_placements(output_dir, blif, new_placements, routings, 1);
	vis_png_draw_placements(output_dir, blif, new_placements, routings, 2);

	FILE *output_json = fopen("output.json", "w");
	vis_json(output_json, blif, new_placements, routings);
	vis_json(output_json, blif, new_placements, routings);
	vis_json(output_json, blif, new_placements, routings);

        free_blif(blif);
	free_cell_library(cl);

	return 0;
}
