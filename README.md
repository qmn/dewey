Dewey
=====
Dewey is intended to be a *fast* placer and router for Minecraft circuits. It
is the C port of [PERSHING](https://github.com/qmn/pershing). An external logic
synthesis program, like [yosys](http://www.clifford.at/yosys/), is needed to
synthesize Verilog into a Berkeley Logic Interchange File (BLIF). Dewey will
take such a file and produce a compacted layout of logic cells and redstone
wires needed to produce the corresponding circuit in Minecraft.

Requirements
------------
- `make`
- a C compiler (e.g., `gcc` or `clang`)
- [LibYAML](http://pyyaml.org/wiki/LibYAML)
- [LibGD](https://libgd.github.io/) and
  [libpng](http://libpng.org/pub/png/libpng.html), for visualizing designs
  without opening them in Minecraft
- [Yosys](http://www.clifford.at/yosys/) to produce BLIF files from source
  Verilog
- Python 2.7 and the PC version of Minecraft ("Java Edition"), to insert
  extracted designs and to interact with them in-game

Setup
-----
First, install the required libraries (LibYAML, LibGD, libpng). As part of
the process, Dewey produces PNG images to visualize the design without
having to open Minecraft. However, we cannot redistribute Minecraft
texture files -- so set the `TEXTURES_FILE` variable in the `Makefile`:

    $ vim Makefile # edit TEXTURES_FILE

Then, build Dewey from source:

    $ make

Usage
-----
To run Dewey, you must first create a BLIF file that corresponds to your
circuit. If you have an input Verilog file, use the `yosys.sh` convenience
script to generate a BLIF file. If you want to use your own synthesis
script, use the provided `quan.lib` file for standard cell mapping (`abc
-liberty quan.lib`, for instance).  If you do not have Yosys installed, we
have provided a 4-bit counter BLIF file for your convenience
(`counter.blif`). For this example, we will synthesize, place, and route a
four-bit counter (source provided in `counter.v`):

    $ scripts/yosys.sh examples/counter.v

We now have a file called `counter.blif`. Run `dewey`:

    $ dewey counter.blif

Dewey is split into, largely, three phases: placement, routing, and
optimization. Each phase can be interrupted by sending SIGINT (by pressing
Control-C). It's possible for a design to have no feasible routing --
Dewey cannot determine this, and may run forever. Do not leave Dewey
running unattended. (Be aware that Dewey is still very experimental. See
`Hacking` for details.)

Inserting your design into a Minecraft world
--------------------------------------------
At this point, you should have a visual representation of the circuit you
have placed-and-routed with Dewey as a PNG file. Also generated as part
of running Dewey is a file called `extraction.yaml`, which is a file
containing the grid of blocks to be placed in the Minecraft world. To
place these blocks into the Minecraft world, follow these instructions:

To read/write Minecraft worlds folders, we use the NBT package. Initialize
it with the `git submodule` command:

    $ git submodule update --init

Find the world folder you want to place your design in. _Make sure you
back-up this world folder! As part of the insertion process, Dewey WILL
overwrite your world data!_ On a Mac, you may find your world folders at
`~/Library/Application Support/minecraft/saves/...`. Then, run the
`inserter.py` Python script:

    $ ./inserter.py extraction.yaml path/to/your/world/folder
    [inserter] reading in extraction
    [inserter] done.
    [inserter] starting insertion...
    [inserter] placing bed of dirt...
    [inserter] placing actual blocks...
    Wrote 617 blocks to Minecraft world ... done.
    [inserter] inserted extraction into path/to/your/world/folder

Your circuit has been placed in the Minecraft world! Open it up and see
for yourself! By default, the design is placed at roughly (y, z, x) = (3,
0, 0). In Minecraft, use the `/tp` command to teleport yourself, taking
care to note that coordinates for the command are given in `<x> <y> <z>`:

    /tp 0 4 0

Hacking
-------
Dewey is experimental software. Here are some things that I have been
meaning to improve.

- Determine the best spacing between standard cells (see, in particular,
  the preprocessor macros `MIN_MARGIN`, `EDGE_MARGIN`, and definitions for
  minimum and maximum window height in `placer.c`). Also see the scoring
  functions.
- Appropriately route in the presence of vertical signal transmission. In
  particular, a signal must approach a via in a special manner, or else
  the proper connection will not occur. See `maze_router.c`, especially
  the routines describing violations occurring near vias.
- Implement better algorithms for detail routing, like the Mikami-Tabuchi
  algorithm.
- Parallelize placement (the `canneal` benchmark in the PARSEC suite is
  essentially this) and routing in a manner that compiles nicely across
  as many platforms as possible.

Why is it called Dewey?
-----------------------
The previous incarnation of this work, PERSHING, is named after a
[ballistic missile system](https://en.wikipedia.org/wiki/MGM-31_Pershing),
the successor of the [Redstone ballistic missile
system](https://en.wikipedia.org/wiki/PGM-11_Redstone). John Pershing held
the rank of General of the Armies, the highest rank in the U.S.
Army. George Dewey, who held the U.S. Navy equivalent rank of Admiral of
the Navy in roughly the same period, could be considered Pershing's peer.
As such, the name "Dewey" is especially fitting as Dewey can be seen as
the "sea" ("C") version of Pershing.

Other Notes
-----------
A publication resulting from the creation of PERSHING, which is the original
version of this project, appeared at the first annual conference on TBD
([SIGTBD'16](http://sigtbd.csail.mit.edu/)), a joke conference at MIT.

A publication describing the speed improvements gained through Dewey
appeared at the third annual conference on TBD (SIGTBD'18).
