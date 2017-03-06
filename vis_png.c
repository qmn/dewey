#include <png.h>
#include <stdio.h>
#include <stdlib.h>

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

png_bytepp load_textures_0()
{
	FILE *f = fopen("/Users/qmn/Library/Application Support/minecraft/textures_0.png", "rb");
	if (!f) {
		printf("[vis_png] error opening textures_0\n");
		return NULL;
	}

	char header[8];
	fread(&header, 1, 8, f);
	if (png_sig_cmp((png_const_bytep)&header, 0, 8)) {
		printf("[vis_png] not a PNG\n");
		goto err;
	}

	png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png) {
		printf("[vis_png] could not create read_struct\n");
		goto err;
	}

	png_infop png_info = png_create_info_struct(png);
	if (!png_info) {
		png_destroy_read_struct(&png, (png_infopp)NULL, (png_infopp)NULL);
		printf("[vis_png] could not create png_info\n");
		goto err;
	}

	png_infop png_end_info = png_create_info_struct(png);
	if (!png_end_info) {
		png_destroy_read_struct(&png, &png_info, (png_infopp)NULL);
		printf("[vis_png] could not create png_end_info\n");
		goto err;
	}

	if (setjmp(png_jmpbuf(png))) {
		printf("err\n");
		png_destroy_read_struct(&png, &png_info, &png_end_info);
		goto err;
	}

	png_init_io(png, f);
	png_set_sig_bytes(png, 8);

	png_set_keep_unknown_chunks(png, PNG_HANDLE_CHUNK_NEVER, NULL, 0);

	png_bytepp row_pointers = png_malloc(png, 512 * sizeof(png_bytep));
	for (int i = 0; i < 512; i++)
		row_pointers[i] = png_malloc(png, 512 * 4 * sizeof(png_byte));

	png_set_rows(png, png_info, row_pointers);

	png_read_png(png, png_info, PNG_TRANSFORM_IDENTITY, NULL);

	png_destroy_read_struct(&png, &png_info, &png_end_info);

	fclose(f);
	return row_pointers;

err:
	fclose(f);
	return NULL;
}

void free_textures_0(png_bytepp row_pointers)
{
	for (int i = 0; i < 512; i++)
		free(row_pointers[i]);
	free(row_pointers);
}

/* draw block at block location x, y in the specified row_data */
void vis_png_draw_block(png_bytepp row_data, png_bytepp textures_0, unsigned char block, int x, int y)
{
	if (block == 0)
		return;

	int row_data_start_x = x * 16 * 4;
	int row_data_start_y = y * 16;

	int block_start_x = -1;
	int block_start_y = -1;

	// find the block in textures_0
	for (int i = 0; i < (sizeof(texture_0_coords) / sizeof(int) / 3); i++) {
		if (texture_0_coords[i].id == block) {
			block_start_x = texture_0_coords[i].x * 16 * 4;
			block_start_y = texture_0_coords[i].y * 16;
		}
	}

	if (block_start_x < 0 || block_start_y < 0) { 
		printf("block id %u not found\n", block);
		return;
	}

	// copy image data from textures_0
	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++) {
			png_bytep loc = &row_data[row_data_start_y + i][row_data_start_x + j * 4];
			png_bytep tex = &textures_0[block_start_y + i][block_start_x + j * 4];
			if (block == 55) {
				if (tex[3]) { // alpha channel
					loc[0] = tex[0];
					loc[3] = 0xFF;
				}
			} else {
				if (tex[3]) {
					for (int k = 0; k < 4; k++) {
						loc[k] = tex[k];
					}
				}
			}
		}
	}
}

void vis_png_draw_placements(struct cell_placements *cp)
{
	FILE *f = fopen("placement.png", "wb");

	if (!f) {
		printf("[vis_png] error opening file\n");
		return;
	}

	png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png) {
		printf("[vis_png] error creating png_struct\n");
		return;
	}

	png_infop png_info = png_create_info_struct(png);
	if (!png_info) {
		printf("[vis_png] could not create png_info\n");
		png_destroy_write_struct(&png, (png_infopp)NULL);
		return;
	}

	if (setjmp(png_jmpbuf(png))) {
		png_destroy_write_struct(&png, &png_info);
		fclose(f);
		return;
	}

	png_init_io(png, f);

	struct dimensions d = compute_placement_dimensions(cp);

	int img_width = d.x * 16;
	int img_height = d.z * 16;
	png_set_IHDR(png, png_info, img_width, img_height, 8, PNG_COLOR_TYPE_RGB_ALPHA,
			PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

	png_write_info(png, png_info);

	png_bytepp textures_0 = load_textures_0();
	if (!textures_0) {
		printf("[vis_png] punt\n");
		return;
	}

	png_bytepp row_data = png_malloc(png, img_height * sizeof(png_bytep));
	for (int i = 0; i < img_height; i++)
		row_data[i] = 0;

	for (int i = 0; i < img_height; i++) {
		row_data[i] = png_malloc(png, img_width * 4 * sizeof(png_byte));
		for (int j = 0; j < img_width * 4 * sizeof(png_byte); j++)
			row_data[i][j] = 0;
	}


	unsigned char *flattened = flatten(cp);
	for (int y = 0; y < d.y; y++) {
		for (int z = 0; z < d.z; z++) {
			for (int x = 0; x < d.x; x++) {
				int id = flattened[y * d.z * d.x + z * d.x + x];
				vis_png_draw_block(row_data, textures_0, id, x, z);
			}
		}
	}

	png_write_rows(png, row_data, img_height);

	free_textures_0(textures_0);

	for (int i = 0; i < img_height; i++)
		png_free(png, row_data[i]);
	png_free(png, row_data);

	png_destroy_write_struct(&png, &png_info);
	fclose(f);

	printf("[vis_png] wrote to placement.png\n");
}

/* create a 3D array representing actual Minecraft block placements */
unsigned char *flatten(struct cell_placements *cp)
{
	struct dimensions d = compute_placement_dimensions(cp);
	int size = d.x * d.y * d.z;
	unsigned char *flat_data = calloc(size, sizeof(unsigned char));

	for (int i = 0; i < cp->n_placements; i++) {
		struct placement p = cp->placements[i];

		struct coordinate c = p.placement;
		struct logic_cell *lc = p.cell;
		struct dimensions lcd = lc->dimensions;

		for (int y = 0; y < lcd.y; y++) {
			for (int z = 0; z < lcd.z; z++) {
				for (int x = 0; x < lcd.x; x++) {
					int fd_off = (c.y + y) * d.z * d.x + (c.z + z) * d.x + (c.x + x);
					int b_off = y * lcd.z * lcd.x + z * lcd.x + x;
					if (fd_off > size)
						continue;
					flat_data[fd_off] = lc->blocks[b_off];
				}
			}
		}
	}

	return flat_data;
}
