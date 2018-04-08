#!/usr/bin/env python2.7

from __future__ import print_function

import sys
import yaml
import numpy as np

from nbt import nbt, world, region

# List of all block IDs in Minecraft, taken from
# http://minecraft.gamepedia.com/Data_values/Block_IDs
block_names = [ "air", "stone", "grass", "dirt", "cobblestone", "planks",
"sapling", "bedrock", "flowing_water", "water", "flowing_lava", "lava", "sand",
"gravel", "gold_ore", "iron_ore", "coal_ore", "log", "leaves", "sponge",
"glass", "lapis_ore", "lapis_block", "dispenser", "sandstone", "noteblock",
"bed", "golden_rail", "detector_rail", "sticky_piston", "web", "tallgrass",
"deadbush", "piston", "piston_head", "wool", "piston_extension",
"yellow_flower", "red_flower", "brown_mushroom", "red_mushroom", "gold_block",
"iron_block", "double_stone_slab", "stone_slab", "brick_block", "tnt",
"bookshelf", "mossy_cobblestone", "obsidian", "torch", "fire", "mob_spawner",
"oak_stairs", "chest", "redstone_wire", "diamond_ore", "diamond_block",
"crafting_table", "wheat", "farmland", "furnace", "lit_furnace",
"standing_sign", "wooden_door", "ladder", "rail", "stone_stairs", "wall_sign",
"lever", "stone_pressure_plate", "iron_door", "wooden_pressure_plate",
"redstone_ore", "lit_redstone_ore", "unlit_redstone_torch", "redstone_torch",
"stone_button", "snow_layer", "ice", "snow", "cactus", "clay", "reeds",
"jukebox", "fence", "pumpkin", "netherrack", "soul_sand", "glowstone",
"portal", "lit_pumpkin", "cake", "unpowered_repeater", "powered_repeater",
"stained_glass", "trapdoor", "monster_egg", "stonebrick",
"brown_mushroom_block", "red_mushroom_block", "iron_bars", "glass_pane",
"melon_block", "pumpkin_stem", "melon_stem", "vine", "fence_gate",
"brick_stairs", "stone_brick_stairs", "mycelium", "waterlily", "nethre_brick",
"nether_brick_fence", "nether_brick_stairs", "nether_wart", "enchanting_table",
"brewing_stand", "cauldron", "end_portal", "end_portal_frame", "end_stone",
"dragon_egg", "redstone_lamp", "lit_redstone_lamp", "double_wooden_slab",
"wooden_slab", "cocoa", "sandstone_stairs", "emerald_ore", "ender_chest",
"tripwire_hook", "tripwire", "emerald_block", "spruce_stairs", "birch_stairs",
"jungle_stairs", "command_block", "beacon", "cobblestone_wall", "flower_pot",
"carrots", "potatoes", "wooden_button", "skull", "anvil", "trapped_chest",
"light_weighted_pressure_plate", "heavy_weighted_pressure_plate",
"unpowered_comparator", "powered_comparator", "daylight_detector",
"redstone_block", "quartz_ore", "hopper", "quartz_block", "quartz_stairs",
"activator_rail", "dropper", "stained_hardened_clay", "stained_glass_pane",
"leaves2", "log2", "acacia_stairs", "dark_oak_stairs", "slime", "barrier",
"iron_trapdoor", "prismarine", "sea_lantern", "hay_block", "carpet",
"hardened_clay", "coal_block", "packed_ice", "double_plant", "standing_banner",
"wall_banner", "daylight_detector_inverted", "red_sandstone",
"red_sandstone_stairs", "double_stone_slab2", "stone_slab2",
"spruce_fence_gate", "birch_fence_gate", "jungle_fence_gate",
"dark_oak_fence_gate", "acacia_fence_gate", "spruce_fence", "birch_fence",
"jungle_fence", "dark_oak_fence", "acacia_fence", "spruce_door", "birch_door",
"jungle_door", "acacia_door", "dark_oak_door", "end_rod", "chorus_plant",
"chorus_flower", "purpur_block", "purpur_pillar", "purpur_stairs",
"purpur_double_slab", "purpur_slab", "end_bricks", "beetroots", "grass_path",
"end_gateway", "repeating_command_block", "chain_command_block", "frosted_ice"
]

class Region:
    """Convenience class for achieving some MPRT tasks."""
    def __init__(self, region):
        self.region = region
        self.chunks = {}

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        # Writing the affected chunks
        for k, v in self.chunks.iteritems():
            (chunk_x, chunk_z) = k
            self.region.write_chunk(chunk_x, chunk_z, v)

    def get_chunk(self, chunk_x, chunk_z):
        key = (chunk_x, chunk_z)
        if not key in self.chunks:
            try:
                self.chunks[key] = self.region.get_chunk(chunk_x, chunk_z)
            except region.InconceivedChunk:
                # create the chunk
                new_chunk = nbt.NBTFile()
                level_tag = nbt.TAG_Compound()
                level_tag.name = "Level"
                level_tag.tags.append(nbt.TAG_Int(name="xPos", value=chunk_x*32))
                level_tag.tags.append(nbt.TAG_Int(name="zPos", value=chunk_z*32))
                level_tag.tags.append(nbt.TAG_List(name="Sections", type=nbt.TAG_Compound))
                new_chunk.tags.append(level_tag)
                self.chunks[key] = new_chunk

        return self.chunks[key]

    def set_chunk(self, chunk_x, chunk_z, chunk):
        key = (chunk_x, chunk_z)
        # If we don't have the chunk in this wrapper, load it
        if key not in self.chunks:
            get_chunk(chunk_x, chunk_z)
        self.chunks[key] = chunk

    def create_empty_section(self, section_y):
        new_section = nbt.TAG_Compound()

        data = nbt.TAG_Byte_Array(u"Data")
        skylight = nbt.TAG_Byte_Array(u"SkyLight")
        blocklight = nbt.TAG_Byte_Array(u"BlockLight")
        y = nbt.TAG_Byte()
        blocks = nbt.TAG_Byte_Array(u"Blocks")

        # TAG_Byte_Array(u'Data'): [2048 byte(s)]
        # TAG_Byte_Array(u'SkyLight'): [2048 byte(s)]
        # TAG_Byte_Array(u'BlockLight'): [2048 byte(s)]
        # TAG_Byte(u'Y'): 0
        # TAG_Byte_Array(u'Blocks'): [4096 byte(s)]

        data.value = bytearray(2048)
        skylight.value = bytearray(2048)
        blocklight.value = bytearray(2048)
        y.name = u"Y"
        y.value = section_y
        blocks.value = bytearray(4096)

        new_section.tags.extend([data, skylight, blocklight, y, blocks])
        return new_section

    def create_section(self, section_y, chunk_x, chunk_z):
        """
        Creates a new section y in chunk (x, z).
        """
        new_section = self.create_empty_section(section_y)
        chunk = self.get_chunk(chunk_x, chunk_z)
        chunk["Level"]["Sections"].append(new_section)

        return new_section

    def get_section(self, section_y, chunk_x, chunk_z):
        """
        Gets the section y in chunk (x, z).
        """
        chunk = self.get_chunk(chunk_x, chunk_z)
        for section in chunk["Level"]["Sections"]:
            if section["Y"].value == section_y:
                return section
        return None

    def set_section(self, y, chunk_x, chunk_z, section):
        chunk = self.get_chunk(chunk_x, chunk_z)
        old_section = self.get_section(y, chunk_x, chunk_z)
        if old_section:
            chunk["Level"]["Sections"].remove(old_section)
        chunk["Level"]["Sections"].append(section)

    def set_block(self, x, y, z, id):
        """
        Sets the region-relative coordinate (x, y, z) to block id.
        """
        chunk_x, offset_x = divmod(x, 16)
        section_y, offset_y = divmod(y, 16)
        chunk_z, offset_z = divmod(z, 16)

        section = self.get_section(section_y, chunk_x, chunk_z)
        if not section:
            section = self.create_section(section_y, chunk_x, chunk_z)

        block_position = (offset_y * 16 * 16) + (offset_z * 16) + offset_x

        if (id > 0xFF):
            raise ValueError("Block id > 255")

        section["Blocks"][block_position] = id

    def get_block(self, x, y, z):
        chunk_x, offset_x = divmod(x, 16)
        section_y, offset_y = divmod(y, 16)
        chunk_z, offset_z = divmod(z, 16)

        section = self.get_section(section_y, chunk_x, chunk_z)
        if not section:
            section = self.create_section(section_y, chunk_x, chunk_z)

        block_position = (offset_y * 16 * 16) + (offset_z * 16) + offset_x

        return section["Blocks"][block_position]

    def set_section_blocks(self, section_y, chunk_x, chunk_z, blocks):
        section = self.get_section(section_y, chunk_x, chunk_z)
        if not section:
            section = self.create_section(section_y, chunk_x, chunk_z)

        section["Blocks"].value = bytearray(blocks)

    def set_redstone(self, x, y, z):
        """Sets a redstone dust piece above the block noted."""
        self.set_block(x, y+1, z, 55)

    def set_data(self, x, y, z, data):
        chunk_x, offset_x = divmod(x, 16)
        section_y, offset_y = divmod(y, 16)
        chunk_z, offset_z = divmod(z, 16)

        section = self.get_section(section_y, chunk_x, chunk_z)
        if not section:
            section = self.create_section(section_y, chunk_x, chunk_z)

        block_position = (offset_y * 16 * 16) + (offset_z * 16) + offset_x

        pos, off = divmod(block_position, 2)

        data_byte = section["Data"].value[pos]

        data_nibble = data & 0xF
        if off == 1:
            data_byte = data_byte & 0xF | data_nibble << 4
        else:
            data_byte = data_byte & 0xF0 | data_nibble

        section["Data"].value[pos] = data_byte

def place_block(world, y, z, x, i, d=0):
    region_x, offset_x = divmod(x, 32*16)
    region_z, offset_z = divmod(z, 32*16)

    with Region(world.get_region(region_x, region_z)) as region:
        region.set_block(offset_x, y, offset_z, i)
        region.set_data(offset_x, y, offset_z, d)

def insert_extracted_layout(world, extracted_layout, offset=(0, 0, 0)):
    """
    Places the extracted layout at offset (y, z, x) in the world.
    """

    blocks, data = extracted_layout

    len_z = blocks.shape[1]
    len_x = blocks.shape[2]
    start_y, start_z, start_x = offset

    height, width, length = blocks.shape
    count = 0

    # place base

    m = 32 * 16
    print("[inserter] placing bed of dirt...")
    for zz in xrange(width / m + 1):
        for xx in xrange(length / m + 1):
            with Region(world.get_region(xx, zz)) as region:
                for z in xrange(min(m, width - zz*m)):
                    for x in xrange(min(m, length - xx*m)):
                        region.set_block(x, start_y-1, z, block_names.index("dirt"))

    print("[inserter] placing actual blocks...")
    for zz in xrange(width / m + 1):
        for xx in xrange(length / m + 1):
            with Region(world.get_region(xx, zz)) as region:
                for y in xrange(height):
                    for z in xrange(min(m, width - zz*m)):
                        for x in xrange(min(m, length - xx*m)):
                            block = blocks[y, zz*m+z, xx*m+x]

                            if block == 0:
                                continue

                            datum = data[y, zz*m+z, xx*m+x]

                            region.set_block(start_x + xx*m + x, start_y + y, start_z + zz*m + z, block)
                            region.set_data(start_x + xx*m + x, start_y + y, start_z + zz*m + z, datum)

                            count += 1
                            msg = "Wrote {} blocks to Minecraft world".format(count)
                            sys.stdout.write("\b" * len(msg))
                            sys.stdout.write(msg)
                            sys.stdout.flush()

    sys.stdout.write(" ... done.\n")
    sys.stdout.flush()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: {} [extraction.yaml] [world folder]".format(sys.argv[0]))
        quit(1)

    extraction_yaml = sys.argv[1]
    world_folder = sys.argv[2]

    extraction = None
    try:
        print("[inserter] reading in extraction")
        with open(extraction_yaml) as f:
            extraction = yaml.load(f)["extraction"]
        print("[inserter] done.")
    except IOError as e:
        print("[inserter] could not read {}".format(extraction_yaml))
        quit(2)

    # process extraction data
    shape = tuple(extraction["dimensions"])
    blocks = np.reshape(extraction["blocks"], shape, order="C")
    data = np.reshape(extraction["data"], shape, order="C")

    print("[inserter] starting insertion...")
    world = world.WorldFolder(world_folder)
    offset = (4, 0, 0)
    insert_extracted_layout(world, (blocks, data), offset)

    print("[inserter] inserted extraction into {}".format(world_folder))

