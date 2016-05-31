#include <stdlib.h>
#include <stdio.h>
#include <yaml.h>

#include "cell.h"

// #define CELL_LIBRARY_DEBUG

void print_blocks(block_t *blocks, struct dimensions d)
{
	int y, z, x;
	int yoff, zoff;

	for (y = 0; y < d.y; y++) {
		yoff = y * d.z * d.x;
		for (z = 0; z < d.z; z++) {
			zoff = z * d.x;
			printf("    ");
			for (x = 0; x < d.x; x++) {
				printf("%3d ", blocks[yoff + zoff + x]);
			}
			printf("\n");
		}
		printf("\n");
	}
}

void print_data(data_t *data, struct dimensions d)
{
	int y, z, x;
	int yoff, zoff;

	for (y = 0; y < d.y; y++) {
		yoff = y * d.z * d.x;
		for (z = 0; z < d.z; z++) {
			zoff = z * d.x;
			printf("    ");
			for (x = 0; x < d.x; x++) {
				printf("%3d ", data[yoff + zoff + x]);
			}
			printf("\n");
		}
		printf("\n");
	}
}

void print_cell_information(struct logic_cell *lc)
{
	printf("logic cell: %s\n", lc->name);
	struct dimensions d = lc->dimensions;
	printf("dimensions: h: %d, w: %d, l: %d\n", d.y, d.z, d.x);
	printf("pins: %d\n", lc->n_pins);
	for (int i = 0; i < lc->n_pins; i++) {
		struct logic_cell_pin *pin = lc->pins[i];
		printf("  pin: %s\n", pin->name);
		printf("    direction: %d\n", pin->direction);
		printf("    facing: %d\n", pin->facing);
		struct coordinate c = pin->coordinate;
		printf("    coordinates: (%d, %d, %d)\n", c.y, c.z, c.x);
		printf("    level: %d\n", pin->level);
		printf("    clock: %d\n", pin->clock);
		printf("\n");
	}
	printf("blocks:\n");
	print_blocks(lc->blocks, lc->dimensions);

	printf("data:\n");
	print_data(lc->data, lc->dimensions);
}

/* Reads the next scalar from the yaml file, returning a copy
 * of the string from the yaml file. */
static char *read_mapped_string(yaml_parser_t *parser)
{
	yaml_event_t event;

	if (!yaml_parser_parse(parser, &event))
		return NULL;

	if (event.type != YAML_SCALAR_EVENT) {
		printf("[cell_library] read_mapped_string expected a scalar at line %lu\n", event.start_mark.line);
		return NULL;
	}

	return strdup((char *)event.data.scalar.value);
}

static int read_mapped_int(yaml_parser_t *parser)
{
	yaml_event_t event;
	if (!yaml_parser_parse(parser, &event))
		return 0;

	if (event.type != YAML_SCALAR_EVENT) {
		printf("[cell_library] read_mapped_int expected a scalar at line %lu\n", event.start_mark.line);
		return 0;
	}

	return atoi((char *)event.data.scalar.value);
}

static int calculate_layout_bytes(struct dimensions d)
{
	return d.y * d.z * d.x;
}

/* Reads a triple-nested sequence of scalars, which must be
 * numbers between 0 and 255, inclusive. Returns non-zero on success. */
static block_t *read_mapped_blocks(yaml_parser_t *parser, struct logic_cell *lc)
{
	struct dimensions dimensions;
	yaml_event_t event;

	enum parser_state { OUTSIDE, Y, Z, X } state = OUTSIDE;

	dimensions.y = 0;
	dimensions.z = 0;
	dimensions.x = 0;

	block_t *blocks = malloc(sizeof(block_t) * calculate_layout_bytes(dimensions));
	unsigned int copied = 0;

	do {
		if (!yaml_parser_parse(parser, &event))
			goto error;

		switch (event.type) {
		case YAML_SEQUENCE_START_EVENT:
			switch (state) {
			case OUTSIDE:
				state = Y;
				break;
			case Y:
				dimensions.y++;
				state = Z;
				break;
			case Z:
				if (dimensions.y == 1)
					dimensions.z++;
				state = X;
				break;
			default:
				printf("[cell_library] read_mapped_blocks got more sequence starts than expected\n");
				goto error;
			}
			break;

		case YAML_SEQUENCE_END_EVENT:
			switch (state) {
			case X:
				state = Z;
				break;
			case Z:
				state = Y;
				break;
			case Y:
				state = OUTSIDE;
				break;
			default:
				printf("[cell_library] read_mapped_blocks got more sequence ends than expected\n");
				break;
			}
			break;

		case YAML_SCALAR_EVENT:
			if (state == X) {
				/* only increase x when z is just 0 */
				if (dimensions.z == 1)
					dimensions.x++;

				/* reallocate the blocks array */
				int new_size = calculate_layout_bytes(dimensions);
				blocks = realloc(blocks, new_size * sizeof(block_t));

				/* set the value of the block */
				int b = atoi((char *)event.data.scalar.value);
				if (b < 0 || b > 255)
					printf("[cell_library] read_mapped_blocks block is < 0 or > 255: %d\n", b);
				blocks[copied++] = (block_t) b;
			}
			break;

		default:
			printf("[cell_library] read_mapped_blocks got unexpected event\n");
			break;
		}
	} while (state != OUTSIDE);

	lc->blocks = blocks;
	lc->dimensions = dimensions;

	return blocks;

error:
	free(blocks);
	return NULL;
}

static data_t *read_mapped_data(yaml_parser_t *parser, struct logic_cell *lc)
{
	struct dimensions dimensions;
	yaml_event_t event;

	enum parser_state { OUTSIDE, Y, Z, X } state = OUTSIDE;

	dimensions.y = 0;
	dimensions.z = 0;
	dimensions.x = 0;

	/* TODO: make data_t nibble-addressed, not byte-addressed */
	data_t *data = malloc(calculate_layout_bytes(dimensions));
	unsigned int copied = 0;

	do {
		if (!yaml_parser_parse(parser, &event)) {
			goto error;
		}

		switch (event.type) {
		case YAML_SEQUENCE_START_EVENT:
			switch (state) {
			case OUTSIDE:
				state = Y;
				break;
			case Y:
				dimensions.y++;
				state = Z;
				break;
			case Z:
				if (dimensions.y == 1)
					dimensions.z++;
				state = X;
				break;
			default:
				printf("[cell_library] read_mapped_data got more sequence starts than expected\n");
				goto error;
			}
			break;

		case YAML_SEQUENCE_END_EVENT:
			switch (state) {
			case X:
				state = Z;
				break;
			case Z:
				state = Y;
				break;
			case Y:
				state = OUTSIDE;
				break;
			default:
				printf("[cell_library] read_mapped_data got more sequence ends than expected\n");
				break;
			}
			break;

		case YAML_SCALAR_EVENT:
			if (state == X) {
				/* only increase x when z is just 0 */
				if (dimensions.z == 1)
					dimensions.x++;

				/* reallocate the blocks array */
				int new_size = calculate_layout_bytes(dimensions);
				data = realloc(data, new_size * sizeof(data_t));

				/* set the value of the block */
				int d = atoi((char *)event.data.scalar.value);
				if (d < 0 || d > 15)
					printf("[cell_library] read_mapped_data data is < 0 or > 15: %d\n", d);
				data[copied++] = (data_t) d;
			}
			break;

		default:
			printf("[cell_library] read_mapped_data got unexpected event\n");
			break;
		}
	} while (state != OUTSIDE);

	lc->data = data;
	lc->dimensions = dimensions;

	return data;

error:
	free(data);
	return NULL;
}

static enum pin_direction read_pin_direction(yaml_parser_t *parser)
{
	yaml_event_t event;
	if (!yaml_parser_parse(parser, &event)) {
		printf("[cell_library] read_pin_direction couldn't parse\n");
		return 0;
	}

	if (event.type != YAML_SCALAR_EVENT) {
		printf("[cell_library] read_pin_direction needs a scalar\n");
	}

	char *direction = (char *)event.data.scalar.value;
	if (strcmp(direction, "input") == 0) {
		return INPUT;
	} else if (strcmp(direction, "output") == 0) {
		return OUTPUT;
	} else {
		printf("[cell_library] read_pin_direction needs \"input\" or \"output\"; got %s\n", direction);
		return 0;
	}
}

static int read_logic_cell_delays(yaml_parser_t *parser, struct logic_cell *lc)
{
	yaml_event_t event;

	enum parser_state { OUTSIDE, MAP } state = OUTSIDE;
	
	do {
		if (!yaml_parser_parse(parser, &event))
			return 0;

		switch (state) {
		case OUTSIDE:
			if (event.type == YAML_MAPPING_START_EVENT)
				state = MAP;
			break;

		case MAP:
			if (event.type == YAML_MAPPING_END_EVENT) {
				state = OUTSIDE;
			} else if (event.type == YAML_SCALAR_EVENT) {
				if (strcmp((char *)event.data.scalar.value, "combinational") == 0) {
					lc->delay_combinational = read_mapped_int(parser);
				}
			}
			break;

		default:
			break;
		}
	} while (state != OUTSIDE);

	return 0;
}

static enum ordinal_direction read_pin_facing(yaml_parser_t *parser)
{
	yaml_event_t event;
	if (!yaml_parser_parse(parser, &event)) {
		printf("[cell_library] read_pin_facing couldn't parse\n");
		return 0;
	}

	if (event.type != YAML_SCALAR_EVENT) {
		printf("[cell_library] read_pin_facing needs a scalar\n");
	}

	char *facing = (char *)event.data.scalar.value;
	if (strcmp(facing, "north") == 0) {
		return NORTH;
	} else if (strcmp(facing, "east") == 0) {
		return EAST;
	} else if (strcmp(facing, "south") == 0) {
		return SOUTH;
	} else if (strcmp(facing, "west") == 0) {
		return WEST;
	} else {
		printf("[cell_library] read_pin_facing needs \"north\", \"east\", \"south\", or \"west\"; got %s\n", facing);
		return 0;
	}
}

static int read_pin_level(yaml_parser_t *parser)
{
	yaml_event_t event;
	if (!yaml_parser_parse(parser, &event)) {
		printf("[cell_library] read_pin_level couldn't parse\n");
		return 0;
	}

	if (event.type != YAML_SCALAR_EVENT) {
		printf("[cell_library] read_pin_level needs a scalar\n");
	}

	return atoi((char *)event.data.scalar.value);
}

static int read_pin_clock(yaml_parser_t *parser)
{
	yaml_event_t event;
	if (!yaml_parser_parse(parser, &event)) {
		printf("[cell_library] read_pin_clock couldn't parse\n");
		return 0;
	}

	if (event.type != YAML_SCALAR_EVENT) {
		printf("[cell_library] read_pin_clock needs a scalar\n");
	}

	char *clock = (char *)event.data.scalar.value;
	if (strcmp(clock, "true"))
		return 1;

	return 0;
}

static struct coordinate read_pin_coordinates(yaml_parser_t *parser)
{
	yaml_event_t event;
	enum parser_state { SEQ, Y, Z, X, DONE } state = SEQ;
	struct coordinate coordinate;


	do {
		if (!yaml_parser_parse(parser, &event)) {
			printf("[cell_library] read_pin_coordinates couldn't parse\n");
			goto error;
		}

		char *scalar_val = (char *)event.data.scalar.value;

		switch (state) {
		case SEQ:
			state = Y;
			break;
		case Y:
			coordinate.y = atoi(scalar_val);
			state = Z;
			break;
		case Z:
			coordinate.z = atoi(scalar_val);
			state = X;
			break;
		case X:
			if (event.type == YAML_SCALAR_EVENT)
				coordinate.x = atoi(scalar_val);
			else if (event.type == YAML_SEQUENCE_END_EVENT)
				state = DONE;
			break;
		default:
			break;
		}
	} while (state != DONE);

error:
	return coordinate;
}

static unsigned int read_mapped_logic_cell_pins(yaml_parser_t *parser, struct logic_cell *lc)
{
	unsigned int n_pins = 0;
	struct logic_cell_pin **pins = malloc(sizeof(struct logic_cell_pin *) * n_pins);

	yaml_event_t event;
	enum parser_state { OUTSIDE, MAP, PIN } state = OUTSIDE;

	struct logic_cell_pin *current = NULL;
	char *name;

	do {
		if (!yaml_parser_parse(parser, &event))
			goto error;

		switch (state) {
		case OUTSIDE:
			if (event.type == YAML_MAPPING_START_EVENT)
				state = MAP;
			break;

		case MAP:
			switch (event.type) {
			case YAML_SCALAR_EVENT:
				name = (char *)event.data.scalar.value;
				current = malloc(sizeof(struct logic_cell_pin));
				current->name = strdup(name);
				current->clock = 0;
				break;

			case YAML_MAPPING_START_EVENT:
				state = PIN;
				break;

			case YAML_MAPPING_END_EVENT:
				state = OUTSIDE;
				break;

			default:
				printf("[cell_library] read_mapped_logic_cell_pins got unexpected event\n");
				break;
			}
			break;

		case PIN:
			switch (event.type) {
			case YAML_SCALAR_EVENT:
				name = (char *)event.data.scalar.value;
				if (strcmp(name, "direction") == 0)
					current->direction = read_pin_direction(parser);
				else if (strcmp(name, "facing") == 0)
					current->facing = read_pin_facing(parser);
				else if (strcmp(name, "coordinates") == 0)
					current->coordinate = read_pin_coordinates(parser);
				else if (strcmp(name, "level") == 0)
					current->level = read_pin_level(parser);
				else if (strcmp(name, "clock") == 0)
					current->clock = read_pin_clock(parser);
				break;

			case YAML_MAPPING_END_EVENT:
				n_pins++;
				pins = realloc(pins, sizeof(struct logic_cell_pin *) * n_pins);
				pins[n_pins - 1] = current;
				state = MAP;
				break;

			default:
				printf("[cell_library] read_mapped_logic_cell_pins expected only scalars or mapping-ends in PIN state\n");
				goto error;
			}
			break;
		}
	} while (state != OUTSIDE);

	lc->n_pins = n_pins;
	lc->pins = pins;
	return n_pins;

error:
	if (current)
		free(current);

	return 0;
}

static unsigned int read_mapped_cells(yaml_parser_t *parser, struct cell_library *cl)
{
	unsigned int n_cells = 0;
	struct logic_cell **cells = malloc(sizeof(struct logic_cell *) * n_cells);

	yaml_event_t event;
	enum parser_state { OUTSIDE, MAP, CELL } state = OUTSIDE;

	struct logic_cell *current = NULL;
	char *name;

	do {
		if (!yaml_parser_parse(parser, &event))
			goto error;

		switch (state) {
		case OUTSIDE:
			switch (event.type) {
			case YAML_MAPPING_START_EVENT:
				state = MAP;
				break;

			default:
				printf("[cell_library] read_mapped_cells expected a mapping-start in OUTSIDE state\n");
				goto error;
			}
			break;

		case MAP:
			switch (event.type) {
			case YAML_SCALAR_EVENT:
				current = malloc(sizeof(struct logic_cell));
				current->name = strdup((char *)event.data.scalar.value);
				printf("[cell_library] reading in library cell \"%s\"\n", current->name);
				break;

			case YAML_MAPPING_START_EVENT:
				state = CELL;
				break;

			case YAML_MAPPING_END_EVENT:
				state = OUTSIDE;
				break;

			default:
				printf("[cell_library] read_mapped_cells expected a scalar, mapping-start, or mapping-end in MAP state\n");
				goto error;
			}
			break;

		case CELL:
			switch (event.type) {
			case YAML_MAPPING_START_EVENT:
				for (; event.type != YAML_MAPPING_END_EVENT; yaml_parser_parse(parser, &event));
				break;

			case YAML_MAPPING_END_EVENT:
				n_cells++;
				cells = realloc(cells, sizeof(struct logic_cell *) * n_cells);
				cells[n_cells - 1] = current;
#ifdef CELL_LIBRARY_DEBUG
				print_cell_information(current);
#endif
				state = MAP;
				break;

			case YAML_SCALAR_EVENT:
				name = (char *)event.data.scalar.value;
				if (strcmp(name, "blocks") == 0) {
					read_mapped_blocks(parser, current);
				} else if (strcmp(name, "data") == 0) {
					read_mapped_data(parser, current);
				} else if (strcmp(name, "pins") == 0) {
					read_mapped_logic_cell_pins(parser, current);
				} else if (strcmp(name, "delay") == 0) {
					read_logic_cell_delays(parser, current);
				} else {
					printf("[cell_library] read_mapped_cells got unexpected key: %s\n", name);
				}
				break;

			default:
				printf("[cell_library] read_mapped_cells expected scalar, mapping-start, or mapping-end in CELL state, got %d\n", event.type);
				goto error;
			}
			break;
		}
	} while (state != OUTSIDE);

	cl->n_cells = n_cells;
	cl->cells = cells;
	return n_cells;

error:
	if (current)
		free(current);

	return 0;
}

struct cell_library *read_cell_library(FILE *f)
{
	yaml_parser_t parser;
	yaml_event_t event;

	int done = 0;

	struct cell_library *cl = malloc(sizeof(struct cell_library));

	yaml_parser_initialize(&parser);
	yaml_parser_set_input_file(&parser, f);

	enum parser_state { OUTSIDE, INSIDE } state = OUTSIDE;
	char *value;

	while (!done) {
		if (!yaml_parser_parse(&parser, &event))
			goto error;

		switch (event.type) {
		case YAML_STREAM_START_EVENT:
#if CELL_LIBRARY_DEBUG
			printf("[cell_library] stream start\n");
#endif
			break;

		case YAML_DOCUMENT_START_EVENT:
#if CELL_LIBRARY_DEBUG
			printf("[cell_library] document start\n");
#endif
			break;

		case YAML_DOCUMENT_END_EVENT:
#if CELL_LIBRARY_DEBUG
			printf("[cell_library] document end\n");
#endif
			state = OUTSIDE;
			break;

		case YAML_ALIAS_EVENT:
#if CELL_LIBRARY_DEBUG
			printf("[cell_library] alias\n");
#endif
			break;

		case YAML_SCALAR_EVENT:
#if CELL_LIBRARY_DEBUG
			printf("[cell_library] scalar (anchor: %s, tag: %s, value: %s)\n", event.data.scalar.anchor, event.data.scalar.tag, event.data.scalar.value);
#endif
			value = (char *)event.data.scalar.value;
			if (state == INSIDE) {
				if (strcmp(value, "library_name") == 0) {
					if ((cl->name = read_mapped_string(&parser)))
						printf("[cell_library] library_name = %s\n", cl->name);
					else
						goto error;
				}

				if (strcmp(value, "cells") == 0) {
					if (!read_mapped_cells(&parser, cl))
						goto error;
				}
			}
			break;

		case YAML_MAPPING_START_EVENT:
			state = INSIDE;
			break;

		default:
			break;
		}

		done = (event.type == YAML_STREAM_END_EVENT);
		yaml_event_delete(&event);
	}

	return cl;

error:
	yaml_parser_delete(&parser);
	free(cl);

	return NULL;
}

static void free_logic_cell(struct logic_cell *logic_cell)
{
	int i;

	free(logic_cell->name);

	for (i = 0; i < logic_cell->n_pins; i++)
		free(logic_cell->pins[i]);

	free(logic_cell->pins);

	free(logic_cell->blocks);
	free(logic_cell->data);
}

void free_cell_library(struct cell_library *cl)
{
	int i;

	for (i = 0; i < cl->n_cells; i++)
		free_logic_cell(cl->cells[i]);

	free(cl->cells);
}
