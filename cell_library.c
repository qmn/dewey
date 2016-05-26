#include <stdlib.h>
#include <stdio.h>
#include <yaml.h>

#include "cell.h"

struct cell_library *read_cell_library(FILE *f)
{
	yaml_parser_t parser;
	yaml_event_t event;

	int done = 0;

	struct cell_library *cl = malloc(sizeof(struct cell_library));

	yaml_parser_initialize(&parser);
	yaml_parser_set_input_file(&parser, f);

	while (!done) {
		if (!yaml_parser_parse(&parser, &event))
			goto error;

		switch (event.type) {
		case YAML_STREAM_START_EVENT:
			printf("[cell_library] stream start\n");
			break;

		case YAML_DOCUMENT_START_EVENT:
			printf("[cell_library] document start\n");
			break;

		case YAML_DOCUMENT_END_EVENT:
			printf("[cell_library] document end\n");
			break;

		case YAML_ALIAS_EVENT:
			printf("[cell_library] alias\n");
			break;

		case YAML_SCALAR_EVENT:
			printf("[cell_library] scalar (anchor: %s, tag: %s, value: %s)\n", event.data.scalar.anchor, event.data.scalar.tag, event.data.scalar.value);
			break;

		case YAML_MAPPING_START_EVENT:
			printf("[cell_library] mapping-start (anchor: %s, tag: %s)\n", event.data.mapping_start.anchor, event.data.mapping_start.tag);
			break;

		case YAML_MAPPING_END_EVENT:
			printf("[cell_library] mapping-end\n");

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
