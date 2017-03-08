#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <gd.h>

#include "extract.h"
#include "router.h"
#include "vis_png.h"

static struct texture_0_coord {
	int id;
	int x;
	int y;
} texture_0_coords[] = {
	{1, 20, 9},
	{55, 18, 13},
	{76, 19, 2},
	{75, 19, 3},
	{94, 19, 6},
	{93, 1, 6},
	{150, 2, 6},
	{149, 1, 6},
	{152, 18, 10},
	{5, 16, 4},
	{69, 10, 13}
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

static gdImagePtr load_textures_0()
{
	FILE *f = fopen("/Users/qmn/Library/Application Support/minecraft/textures_0.png", "rb");
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
	int repeater_angle[] = {180, 90, 0, 270}; // data & 0x3 -> angle
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
	case 149:
	case 150: // comparator
		gdImageCopyRotated(im, textures_0, im_start_x + 8, im_start_y + 8, texture_start_x, texture_start_y, 16, 16, repeater_angle[data & 0x3]);
		break;
	default:
		gdImageCopy(im, textures_0, im_start_x, im_start_y, texture_start_x, texture_start_y, 16, 16);
		break;
	}
}

void vis_png_draw_placements(struct cell_placements *cp, struct routings *rt)
{
	struct dimensions d = compute_placement_dimensions(cp);

	int img_width = d.x * 16;
	int img_height = d.z * 16;

	gdImagePtr im = gdImageCreateTrueColor(img_width, img_height);

	gdImageSaveAlpha(im, 1);
	int transparent = gdImageColorAllocateAlpha(im, 0xff, 0xff, 0xff, 0x7f);
	gdImageFill(im, 0, 0, transparent);

	gdImagePtr textures_0 = load_textures_0();
	if (!textures_0) {
		printf("[vis_png] punt\n");
		return;
	}

	struct extraction *e = extract(cp, rt);
	for (int y = 0; y < d.y; y++) {
		for (int z = 0; z < d.z; z++) {
			for (int x = 0; x < d.x; x++) {
				int id = e->blocks[y * d.z * d.x + z * d.x + x];
				int data = e->data[y * d.z * d.x + z * d.x + x];
				vis_png_draw_block(im, textures_0, id, x, z, data);
			}
		}
	}
	free_extraction(e);

	free_textures_0(textures_0);

	FILE *f = fopen("placement.png", "wb");
	gdImagePng(im, f);
	fclose(f);

	printf("[vis_png] wrote to placement.png\n");
}

