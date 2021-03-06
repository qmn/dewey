#include <png.h>
#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <gd.h>
#include <gdfonts.h>
#include <assert.h>
#include <sys/param.h>

#include "extract.h"
#include "router.h"
#include "vis_png.h"

// block id, col, row
static struct texture_0_coord {
	int id;
	int x;
	int y;
} texture_0_coords[] = {
	{1, 20, 9},    // stone
	{55, 18, 13},
	{76, 19, 2},
	{75, 19, 3},
	{94, 19, 6}, // lit repeater
	{93, 19, 5}, // unlit repeater
	{150, 2, 6},
	{149, 1, 6},
	{152, 18, 10},
	{5, 16, 4},
	{69, 10, 13},
	{123, 18, 15},
	{29, 11, 15}, // piston bottom
	{165, 20, 4}, // slime
	{3, 8, 6}       // dirt
};

static gdImagePtr redstone_mask(gdImagePtr textures_0)
{
	gdImagePtr rs = gdImageCreateTrueColor(16, 16);
	gdImageSaveAlpha(rs, 1);
	int transparent = gdImageColorAllocateAlpha(rs, 0xff, 0xff, 0xff, 0x7f);
	gdImageFill(rs, 0, 0, transparent);
	gdImageCopy(rs, textures_0, 0, 0, 18*16, 13*16, 16, 16);
	gdImageColor(rs, 0, -255, -255, 0);
	return rs;
}

#ifndef TEXTURES_FILE
# error "Define TEXTURES_FILE to use the vis_png routines"
#endif
static gdImagePtr load_textures_0()
{
	FILE *f = fopen(TEXTURES_FILE, "rb");
	if (!f) {
		printf("[vis_png] error opening textures_0\n");
		return NULL;
	}

	return gdImageCreateFromPng(f);
}

void free_textures_0(gdImagePtr textures_0)
{
	gdImageDestroy(textures_0);
}

/* draw block at block location x, y in the specified row_data */
void vis_png_draw_block(gdImagePtr im, gdImagePtr textures_0, block_t block, int x, int y, data_t data)
{
	if (block == 0)
		return;

	int im_start_x = x * 16;
	int im_start_y = y * 16;

	int texture_start_x = -1;
	int texture_start_y = -1;

	// find the block in textures_0
	for (int i = 0; i < (sizeof(texture_0_coords) / sizeof(int) / 3); i++) {
		if (texture_0_coords[i].id == block) {
			texture_start_x = texture_0_coords[i].x * 16;
			texture_start_y = texture_0_coords[i].y * 16;
		}
	}

	if (texture_start_x < 0 || texture_start_y < 0) { 
		printf("block id %u not found\n", block);
		return;
	}

	gdImagePtr rs = redstone_mask(textures_0);
	int repeater_angle[] = {0, 270, 180, 90}; // data & 0x3 -> angle
	int torch_angle[] = {0, 270, 90, 180, 0, 0};

	switch (block) {
	case 55: // redstone
		gdImageCopy(im, rs, im_start_x, im_start_y, 0, 0, 16, 16);
		gdImageDestroy(rs);
		break;
	case 75:
	case 76: // torch
		gdImageCopyRotated(im, textures_0, im_start_x + 8, im_start_y + 8, texture_start_x, texture_start_y, 16, 16, torch_angle[data]);
		break;
	case 93:
	case 94: // repeater
		gdImageCopyRotated(im, textures_0, im_start_x + 8, im_start_y + 8, texture_start_x, texture_start_y, 16, 16, repeater_angle[data & 0x3]);
		break;
	case 149:
	case 150: // comparator
		gdImageCopyRotated(im, textures_0, im_start_x + 8, im_start_y + 8, texture_start_x, texture_start_y, 16, 16, repeater_angle[data & 0x3]);
		break;
	default:
		gdImageCopy(im, textures_0, im_start_x, im_start_y, texture_start_x, texture_start_y, 16, 16);
		break;
	}
}

#define ARBITRARY_LIMIT 1000

#define min(x, y) (x <= y ? x : y)

void vis_png_draw_placements(char *output_dir, struct blif *blif, struct cell_placements *cp, struct routings *rt, int layer)
{
	struct extraction *e = extract(cp, rt);
	struct dimensions d = e->dimensions;

	switch (layer) {
	case 0: d.y = min(2, d.y); break;
	case 1: d.y = min(5, d.y); break;
	case 2: d.y = min(10, d.y); break;
	default: break;
	}

	int img_width = d.x * 16;
	int img_height = d.z * 16;

	printf("[vis_png_draw_placements] image dimensions to be %d x %d (%d by %d blocks)\n", img_width, img_height, d.x, d.z);
	assert(d.x > 0 && d.z > 0 && d.x < ARBITRARY_LIMIT && d.z < ARBITRARY_LIMIT);

	gdImagePtr im = gdImageCreateTrueColor(img_width, img_height);

	gdImageSaveAlpha(im, 1);
	int transparent = gdImageColorAllocateAlpha(im, 0xff, 0xff, 0xff, 0x7f);
	int black = gdImageColorAllocateAlpha(im, 0, 0, 0, 0);
	int white = gdImageColorAllocateAlpha(im, 0xff, 0xff, 0xff, 0);
	gdImageFill(im, 0, 0, transparent);

	gdImagePtr textures_0 = load_textures_0();
	if (!textures_0) {
		printf("[vis_png] punt\n");
		return;
	}

	/* draw blocks */
	for (int y = 0; y < d.y; y++) {
		for (int z = 0; z < d.z; z++) {
			for (int x = 0; x < d.x; x++) {
				int id = e->blocks[y * d.z * d.x + z * d.x + x];
				int data = e->data[y * d.z * d.x + z * d.x + x];
				if (y == 0 && id == 55)
					vis_png_draw_block(im, textures_0, 1, x, z, 0); // underlying stone
				vis_png_draw_block(im, textures_0, id, x, z, data);
			}
		}
	}
	free_extraction(e);

	/* draw input/output pin labels */
	for (int i = 0; i < cp->n_placements; i++) {
		if (strcmp("input_pin", cp->placements[i].cell->name) == 0 ||
		    strcmp("output_pin", cp->placements[i].cell->name) == 0) {
			struct coordinate c = cp->placements[i].placement;
			net_t n = cp->placements[i].nets[0];
			unsigned char *s = (unsigned char *)get_net_name(blif, n);
			gdFontPtr font = gdFontGetSmall();
			gdImageString(im, font, c.x * 16 + 8, c.z * 16 + 4, s, white);
			gdImageString(im, font, c.x * 16 + 7, c.z * 16 + 3, s, black);
		}
	}

	for (int x = 0; x < d.x; x++) {
		gdFontPtr font = gdFontGetSmall();
		char buf[3];
		snprintf(buf, 3, "%2d", x%100);
		gdImageString(im, font, x * 16 + 2, 0, (unsigned char *)buf, white);
	}

	for (int z = 1; z < d.z; z++) {
		gdFontPtr font = gdFontGetSmall();
		char buf[3];
		snprintf(buf, 3, "%2d", z%100);
		gdImageString(im, font, 2, 16 * z, (unsigned char *)buf, white);
	}

	free_textures_0(textures_0);

	char *bn;
	asprintf(&bn, "layer%d.png", layer);
	char fn[MAXPATHLEN];
	strncpy(fn, output_dir, MAXPATHLEN);
	strncat(fn, bn, MAXPATHLEN-strnlen(bn, MAXPATHLEN));
	FILE *f = fopen(fn, "wb");
	gdImagePng(im, f);
	fclose(f);

	printf("[vis_png] wrote to %s\n", fn);
}

